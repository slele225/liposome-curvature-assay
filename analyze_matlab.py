"""
Analyze MATLAB detection results (detection_v2.mat) and export filtered
puncta amplitude values as a tab-separated text file.

Filtering logic:
  - Lipid/master channel amplitude must exceed
    mean(c_master) + k_std * std(c_master)
  - Hypothesis test value (hval_Ar) must equal 1 in the lipid channel
"""

import os
import re
import glob
import sys
import argparse
import h5py
import numpy as np


# ── MAT-file helpers ────────────────────────────────────────────────────

def _deref_and_concat(f: h5py.File, ref_arr: np.ndarray):
    """Dereference an array of HDF5 object references and concatenate."""
    pieces = []
    for r in ref_arr.ravel(order="F"):
        if isinstance(r, h5py.Reference) and r:
            a = np.squeeze(np.array(f[r]))
            if a.size == 0:
                continue
            pieces.append(a)

    if not pieces:
        return np.array([])

    norm = [p.reshape(-1, 1) if p.ndim == 1 else p for p in pieces]
    return np.concatenate(norm, axis=0)


def load_frameinfo_field(mat_path: str, field: str):
    """Load a field from frameInfo in a MATLAB v7.3 .mat file."""
    with h5py.File(mat_path, "r") as f:
        ds = f["frameInfo"][field]
        arr = np.array(ds)

        if arr.dtype == h5py.ref_dtype or arr.dtype == object:
            out = _deref_and_concat(f, arr)
        else:
            out = np.squeeze(arr)

        if out.ndim == 1:
            out = out.reshape(-1, 1)

        return out


# ── Utilities ───────────────────────────────────────────────────────────

def find_cell_dirs(condition_folder: str):
    """Return sorted list of cell* directories under a condition folder."""
    cells = [
        c
        for c in glob.glob(os.path.join(condition_folder, "cell*"))
        if os.path.isdir(c)
    ]
    return sorted(
        cells,
        key=lambda p: (
            int(re.search(r"cell(\d+)", p, re.I).group(1))
            if re.search(r"cell(\d+)", p, re.I)
            else 10**12
        ),
    )


def parse_channel_list(s: str):
    parts = [p.strip() for p in s.split(",") if p.strip()]
    if not parts:
        raise ValueError("Channel list is empty.")
    return parts


# ── Main ────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Analyze MATLAB detection results and export filtered puncta A values."
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Condition folder containing cell*/ch*/Detection/detection_v2.mat",
    )
    parser.add_argument(
        "--channels",
        required=True,
        help='Comma-separated channel names, e.g. "ch1,ch2"',
    )
    parser.add_argument(
        "--lipid-channel",
        required=True,
        help='Name of the lipid/master channel, e.g. "ch1"',
    )
    parser.add_argument(
        "--k-std",
        type=float,
        default=2.0,
        help="Threshold multiplier: mean(c) + k·std(c)  (default: 2.0)",
    )
    parser.add_argument(
        "--output-name",
        default="filtered_puncta_A_values.txt",
        help="Output filename (default: filtered_puncta_A_values.txt)",
    )
    args = parser.parse_args()

    condition_folder = args.input
    channel_names = parse_channel_list(args.channels)
    master_name = args.lipid_channel
    k_std = args.k_std
    raw_txt_name = args.output_name

    if not os.path.isdir(condition_folder):
        print(f"Error: input folder does not exist: {condition_folder}")
        sys.exit(1)

    if master_name not in channel_names:
        print(f"Error: lipid channel '{master_name}' not in --channels {channel_names}")
        sys.exit(1)

    master_idx = channel_names.index(master_name)
    other_names = [ch for ch in channel_names if ch != master_name]

    cell_dirs = find_cell_dirs(condition_folder)
    if not cell_dirs:
        print("Error: no cell directories found.")
        sys.exit(1)

    # ── Process each cell ───────────────────────────────────────────────
    raw_rows = []
    total_cells = 0
    total_seen = 0
    total_kept = 0

    print("=" * 70)
    print("MATLAB DETECTION ANALYSIS")
    print("=" * 70)
    print(f"Condition folder: {condition_folder}")
    print(f"Channels:         {channel_names}")
    print(f"Lipid channel:    {master_name}")
    print(f"Other channels:   {other_names}")
    print(f"k_std:            {k_std}")
    print("-" * 70)

    for cell_dir in cell_dirs:
        cell_name = os.path.basename(cell_dir)
        det_path = os.path.join(
            cell_dir, master_name, "Detection", "detection_v2.mat"
        )

        if not os.path.isfile(det_path):
            print(f"{cell_name}: missing detection file at {det_path}")
            continue

        total_cells += 1

        try:
            A = load_frameinfo_field(det_path, "A")
            c = load_frameinfo_field(det_path, "c")
            h = load_frameinfo_field(det_path, "hval_Ar")

            n = min(A.shape[0], c.shape[0], h.shape[0])
            A, c, h = A[:n, :], c[:n, :], h[:n, :]

            for arr_name, arr in [("A", A), ("c", c), ("h", h)]:
                if master_idx >= arr.shape[1]:
                    raise IndexError(
                        f"Master index {master_idx} out of bounds for "
                        f"{arr_name} with {arr.shape[1]} columns."
                    )

            c_master = c[:, master_idx]
            A_master = A[:, master_idx]
            c_thr = float(np.mean(c_master) + k_std * np.std(c_master))

            mask = (A_master > c_thr) & (h[:, master_idx] == 1)
            kept = int(np.sum(mask))
            total_seen += n
            total_kept += kept

            if kept > 0:
                for idx in np.where(mask)[0]:
                    row = [f"{cell_name}|row{idx}"]
                    for ch_name in channel_names:
                        ch_idx = channel_names.index(ch_name)
                        if ch_idx >= A.shape[1]:
                            raise IndexError(
                                f"Channel '{ch_name}' index {ch_idx} out of "
                                f"bounds for A with {A.shape[1]} columns."
                            )
                        row.append(A[idx, ch_idx])
                    raw_rows.append(row)

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

    # ── Write output ────────────────────────────────────────────────────
    if not raw_rows:
        print("No puncta passed the filters.")
        sys.exit(1)

    raw_txt_path = os.path.join(condition_folder, raw_txt_name)
    headers = ["source_image"] + [f"A_{ch}" for ch in channel_names]

    with open(raw_txt_path, "w", encoding="utf-8") as fout:
        fout.write("# Filtered puncta A values\n")
        fout.write(f"# Lipid/master channel: {master_name}\n")
        fout.write(f"# Channels in order: {', '.join(channel_names)}\n")
        fout.write(f"# Threshold: A_master > mean(c) + {k_std} * std(c)\n")
        fout.write(f"# Validation: hval == 1 in {master_name}\n")
        fout.write("\t".join(headers) + "\n")

        for row in raw_rows:
            source = row[0]
            vals = [f"{v:.10g}" for v in row[1:]]
            fout.write("\t".join([source] + vals) + "\n")

    print(f"\nDONE: {len(raw_rows)} points saved to {raw_txt_path}")
    sys.exit(0)


if __name__ == "__main__":
    main()
