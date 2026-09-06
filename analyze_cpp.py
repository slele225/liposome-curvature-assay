"""
Analyze native cme_detect results (detection_cpp.tsv) and export puncta
amplitude values as tab-separated text files.

This is the native-backend counterpart of analyze_matlab.py. It reads the
per-movie table written by native_detection/cme_detect under

    <condition>/cell*/<master>/Detection/detection_cpp.tsv

and applies exactly the same filter, through the same functions imported
from analyze_matlab.py:

  - master/lipid channel amplitude must exceed
    mean(c_master) + k_std * std(c_master)
  - hval_Ar handling is set by --hval-filter (default 'master': require
    hval == 1 in the lipid/master channel only; slave channels are NOT
    required to pass, so weak slave-channel signal is kept)

The two output files (raw_puncta_values.txt, filtered_puncta_A_values.txt)
have the same columns and layout as analyze_matlab.py's, so the downstream
scripts (dls_calibration.py, plot_curvature.py, plot_histograms.py, ...)
consume them unchanged.

Channel order: --channels must be the exact list that was passed to
cme_detect --channels (master/source channel first) and --lipid-channel
must equal cme_detect --master. The TSV header records the order
cme_detect used, so a mismatch is detected and rejected rather than
silently reinterpreted.
"""

import os
import sys
import argparse
import numpy as np

from analyze_matlab import (
    compute_filter_mask,
    describe_hval_policy,
    find_cell_dirs,
    parse_channel_list,
    write_filtered_table,
    write_raw_table,
)


DEFAULT_DETECTION_FILE = "detection_cpp.tsv"
HVAL_PREFIX = "hval_Ar_"


# ── TSV helpers ─────────────────────────────────────────────────────────

def read_detection_tsv(path: str):
    """Read a cme_detect per-movie TSV. Returns (header, columns, n_rows)."""
    with open(path, "r", encoding="utf-8") as f:
        header_line = f.readline()
        if not header_line:
            raise ValueError(f"empty detection file: {path}")
        header = header_line.rstrip("\r\n").split("\t")
        rows = [line.rstrip("\r\n").split("\t") for line in f if line.strip()]

    for i, r in enumerate(rows):
        if len(r) != len(header):
            raise ValueError(
                f"{path}: row {i} has {len(r)} fields, header has {len(header)}"
            )
    cols = {name: [r[i] for r in rows] for i, name in enumerate(header)}
    return header, cols, len(rows)


def tsv_channel_order(header):
    """Channel names in the order cme_detect wrote them (hval_Ar_<ch> columns)."""
    return [h[len(HVAL_PREFIX):] for h in header if h.startswith(HVAL_PREFIX)]


def load_channel_arrays(cols, n: int, channel_names):
    """Build (n, n_channels) A / c / hval_Ar arrays in --channels order."""
    k = len(channel_names)
    A = np.empty((n, k), dtype=float)
    c = np.empty((n, k), dtype=float)
    h = np.empty((n, k), dtype=float)
    for j, ch in enumerate(channel_names):
        for arr, key in ((A, f"A_{ch}"), (c, f"c_{ch}"), (h, f"{HVAL_PREFIX}{ch}")):
            if key not in cols:
                raise KeyError(f"column '{key}' not found in detection file")
            arr[:, j] = np.array([float(v) for v in cols[key]], dtype=float)
    return A, c, h


# ── Main ────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Analyze native cme_detect results and export filtered "
                    "puncta A values (same filter and outputs as "
                    "analyze_matlab.py)."
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Condition folder containing cell*/<master>/Detection/"
             "detection_cpp.tsv (the folder given to cme_detect --input)",
    )
    parser.add_argument(
        "--channels",
        required=True,
        help="Comma-separated channel folder names in the EXACT order they "
             "were passed to cme_detect --channels. The first entry is the "
             'master/source channel. E.g. "ch2,ch1" if ch2 was the master.',
    )
    parser.add_argument(
        "--lipid-channel",
        required=True,
        help="Name of the lipid/master channel. Must equal the FIRST entry "
             "of --channels and the value given to cme_detect --master.",
    )
    parser.add_argument(
        "--k-std",
        type=float,
        default=2.0,
        help="Threshold multiplier: mean(c) + k·std(c)  (default: 2.0)",
    )
    parser.add_argument(
        "--hval-filter",
        choices=["master", "none", "all"],
        default="master",
        help="How the hval_Ar hypothesis test is used for filtering. "
             "'master' (default): require hval == 1 in the lipid/master "
             "channel only — slave channels are NOT required to pass, so "
             "weak slave-channel signal is kept. "
             "'none': hval is not used as a filtering criterion. "
             "'all': require hval == 1 in every requested channel.",
    )
    parser.add_argument(
        "--output-name",
        default="filtered_puncta_A_values.txt",
        help="Filtered output filename (default: filtered_puncta_A_values.txt)",
    )
    parser.add_argument(
        "--raw-output-name",
        default="raw_puncta_values.txt",
        help="Unfiltered diagnostic output filename "
             "(default: raw_puncta_values.txt). Always written, even when "
             "no puncta pass the filter.",
    )
    parser.add_argument(
        "--detection-file",
        default=DEFAULT_DETECTION_FILE,
        help="Per-cell detection table name under <cell>/<master>/Detection/ "
             f"(default: {DEFAULT_DETECTION_FILE}). Regression use: "
             "detection_matlab.tsv as exported by "
             "native_detection/reference_tools/run_matlab_reference.m.",
    )
    args = parser.parse_args()

    condition_folder = args.input
    channel_names = parse_channel_list(args.channels)
    master_name = args.lipid_channel
    k_std = args.k_std
    hval_policy = args.hval_filter
    filtered_txt_name = args.output_name
    raw_txt_name = args.raw_output_name
    det_name = args.detection_file

    if not os.path.isdir(condition_folder):
        print(f"Error: input folder does not exist: {condition_folder}")
        sys.exit(1)

    if master_name not in channel_names:
        print(f"Error: lipid channel '{master_name}' not in --channels {channel_names}")
        sys.exit(1)

    if channel_names[0] != master_name:
        print(
            f"Error: --lipid-channel '{master_name}' is not the FIRST entry in "
            f"--channels {channel_names}.\n"
            "--channels must be supplied in the exact order the channels were "
            "passed to cme_detect --channels, whose first entry is the "
            "master/source channel (cme_detect --master), so --lipid-channel "
            "must equal the first entry of --channels.\n"
            "For example, if cme_detect was run with --channels ch2,ch1 "
            "--master ch2, use:  --channels ch2,ch1 --lipid-channel ch2\n"
            "Refusing to continue: a mismatched order would attach detection "
            "columns to the wrong physical channel."
        )
        sys.exit(1)

    master_idx = channel_names.index(master_name)
    other_names = [ch for ch in channel_names if ch != master_name]

    cell_dirs = find_cell_dirs(condition_folder)
    if not cell_dirs:
        print("Error: no cell directories found.")
        sys.exit(1)

    # ── Process each cell ───────────────────────────────────────────────
    all_rows = []
    filtered_rows = []
    total_cells = 0
    total_seen = 0
    total_kept = 0

    print("=" * 70)
    print("NATIVE (cme_detect) DETECTION ANALYSIS")
    print("=" * 70)
    print(f"Condition folder: {condition_folder}")
    print(f"Channels:         {channel_names}")
    print(f"Lipid channel:    {master_name}")
    print(f"Other channels:   {other_names}")
    print("cme_detect channel mapping:")
    for i, ch_name in enumerate(channel_names):
        tag = "  [MASTER / LIPID]" if ch_name == master_name else ""
        print(f"  column {i} -> {ch_name}{tag}")
    print(f"k_std:            {k_std}")
    print(f"hval policy:      {hval_policy} ({describe_hval_policy(hval_policy, master_name, channel_names)})")
    print("-" * 70)

    for cell_dir in cell_dirs:
        cell_name = os.path.basename(cell_dir)
        det_path = os.path.join(cell_dir, master_name, "Detection", det_name)

        if not os.path.isfile(det_path):
            print(f"{cell_name}: missing detection file at {det_path}")
            if os.path.isdir(os.path.join(condition_folder, "cme_detect_output")):
                print(
                    f"{cell_name}: cme_detect_output/ exists but the per-cell "
                    f"table is absent — was cme_detect run with "
                    f"--no-matlab-layout, or with a different --master?"
                )
            continue

        total_cells += 1

        # Channel-order check against what cme_detect actually wrote. This
        # is a hard error (not a per-cell warning) because it means the
        # whole run was described inconsistently.
        header, cols, n = read_detection_tsv(det_path)
        tsv_order = tsv_channel_order(header)
        if tsv_order[: len(channel_names)] != channel_names:
            print(
                f"Error: channel order mismatch in {det_path}.\n"
                f"cme_detect wrote this table with --channels "
                f"{','.join(tsv_order)} (master {tsv_order[0] if tsv_order else '?'}), "
                f"but --channels {','.join(channel_names)} --lipid-channel "
                f"{master_name} was given here.\n"
                "Pass exactly the channel list used for cme_detect (master "
                "first; trailing slave channels may be omitted)."
            )
            sys.exit(1)

        try:
            A, c, h = load_channel_arrays(cols, n, channel_names)

            passes_A, passes_h, mask, c_thr = compute_filter_mask(
                A, c, h, master_idx, k_std, hval_policy, len(channel_names)
            )
            kept = int(np.sum(mask))
            total_seen += n
            total_kept += kept

            # Every punctum goes into the raw diagnostic rows — no
            # thresholding, no hval filtering, negative/zero A kept.
            for idx in range(n):
                row = [f"{cell_name}|row{idx}"]
                for ch_idx in range(len(channel_names)):
                    row.extend([A[idx, ch_idx], c[idx, ch_idx], h[idx, ch_idx]])
                row.extend(
                    [int(passes_A[idx]), int(passes_h[idx]), int(mask[idx])]
                )
                all_rows.append(row)

            if kept > 0:
                for idx in np.where(mask)[0]:
                    row = [f"{cell_name}|row{idx}"]
                    for ch_idx in range(len(channel_names)):
                        row.append(A[idx, ch_idx])
                    filtered_rows.append(row)

            print(f"{cell_name}: kept {kept}/{n} | c_thr = {c_thr:.3f}")

        except Exception as e:
            print(f"{cell_name}: ERROR: {e}")

    # ── Summary ─────────────────────────────────────────────────────────
    print("-" * 70)
    print("SUMMARY")
    print(f"Cells processed:   {total_cells}")
    print(f"Total puncta seen: {total_seen}")
    print(f"Total puncta kept: {total_kept}")
    if total_seen > 0:
        print(f"Fraction kept:     {total_kept / total_seen:.4f}")

    # ── Write raw diagnostic output (ALWAYS, before any filter exit) ────
    raw_txt_path = os.path.join(condition_folder, raw_txt_name)
    write_raw_table(
        raw_txt_path, all_rows, channel_names, master_name, k_std, hval_policy,
        loaded_from=det_name,
        source_line=(
            f"cell*/{master_name}/Detection/{det_name} "
            f"(cme_detect columns A_<ch>, c_<ch>, hval_Ar_<ch>)"
        ),
    )

    print(f"\nRaw diagnostic values ({len(all_rows)} rows) saved to {raw_txt_path}")

    # ── Write filtered output ───────────────────────────────────────────
    if not filtered_rows:
        print("No puncta passed the filters.")
        print("(The raw diagnostic file above was still written.)")
        sys.exit(1)

    filtered_txt_path = os.path.join(condition_folder, filtered_txt_name)
    write_filtered_table(
        filtered_txt_path, filtered_rows, channel_names, master_name,
        k_std, hval_policy,
    )

    print(f"DONE: {len(filtered_rows)} points saved to {filtered_txt_path}")
    sys.exit(0)


if __name__ == "__main__":
    main()
