"""
Split multi-frame TIFFs into per-channel folders for MATLAB detection.

Parses Olympus FV3000 ImageJ metadata to detect which lasers were active
(transmissivity > 0), maps channels to lasers via the laserDataId
entries, reads per-channel PMT voltage, and groups output by the
unique (wavelength, transmissivity, voltage) combination across all
channels in each TIFF.

Folder name format:
    {wave1}nm_{trans1}pct_{volt1}V_{wave2}nm_{trans2}pct_{volt2}V_...
sorted by wavelength ascending. Channel folders ch1, ch2, ... follow
the same wavelength-ascending order, so ch1 is always the lowest-
wavelength active laser.

Lambda-phase channels (laserDataId containing "phase_Lambda") are
ignored; only main imaging phases ("phase_1") are kept.
"""

import os
import re
import sys
import argparse
import tifffile
import numpy as np


# ── Metadata parsing ───────────────────────────────────────────────────

LASER_TRANS_RE = re.compile(
    r"Laser\s+(LD(\d+))\s+transmissivity\s*=\s*([0-9.]+)"
)
CHANNEL_LASER_RE = re.compile(
    r"channel\s+laserDataId\s+#0*(\d+)\s*=\s*(LD\d+)_([^\s]+)"
)
PMT_VOLTAGE_RE = re.compile(
    r"pmt\s+voltage\s+#0*(\d+)\s*=\s*([0-9.]+)"
)


def parse_metadata(tiff_path):
    """
    Return active imaging channels from FV3000 metadata.

    Each channel: {laser_id, wavelength, transmissivity, voltage}.
    Sorted by wavelength ascending. Returns [] if no active lasers
    were found or the metadata is missing.
    """
    with tifffile.TiffFile(tiff_path) as tif:
        ij_meta = tif.imagej_metadata or {}
        info = ij_meta.get("Info", "")

    if not info:
        return []

    # 1. Active lasers: transmissivity > 0.
    lasers = {}
    for m in LASER_TRANS_RE.finditer(info):
        laser_id = m.group(1)            # e.g. "LD488"
        wavelength = int(m.group(2))     # e.g. 488
        trans = float(m.group(3))
        if trans > 0:
            lasers[laser_id] = (wavelength, trans)

    if not lasers:
        return []

    # 2. Channel index -> laser_id, skipping lambda-phase channels.
    channel_lasers = {}
    for m in CHANNEL_LASER_RE.finditer(info):
        idx = int(m.group(1))
        laser_id = m.group(2)
        suffix = m.group(3)
        if "phase_Lambda" in suffix:
            continue
        if laser_id in lasers:
            channel_lasers[idx] = laser_id

    # 3. Channel index -> voltage (integer).
    voltages = {}
    for m in PMT_VOLTAGE_RE.finditer(info):
        idx = int(m.group(1))
        voltages[idx] = int(float(m.group(2)))

    # 4. Build channel records.
    channels = []
    for idx, laser_id in channel_lasers.items():
        wavelength, trans = lasers[laser_id]
        volt = voltages.get(idx)
        if volt is None:
            continue
        channels.append({
            "laser_id": laser_id,
            "wavelength": wavelength,
            "transmissivity": trans,
            "voltage": volt,
        })

    channels.sort(key=lambda c: c["wavelength"])
    return channels


def group_name(channels):
    """Folder name like 488nm_5.0pct_580V_561nm_3.2pct_500V."""
    parts = [
        f"{c['wavelength']}nm_{c['transmissivity']:.1f}pct_{c['voltage']}V"
        for c in channels
    ]
    return "_".join(parts)


def channels_summary(channels):
    """Summary like '[488 (5.0%, 580V), 561 (3.2%, 500V)]'."""
    parts = [
        f"{c['wavelength']} ({c['transmissivity']:.1f}%, {c['voltage']}V)"
        for c in channels
    ]
    return "[" + ", ".join(parts) + "]"


# ── Frames argument ────────────────────────────────────────────────────

def parse_int_list(s):
    if not s.strip():
        return []
    parts = [p.strip() for p in s.split(",") if p.strip()]
    out = []
    for p in parts:
        if p.startswith("-") or not p.isdigit():
            raise ValueError(f"Invalid index: '{p}'")
        out.append(int(p))
    return out


# ── Image processing ───────────────────────────────────────────────────

def center_crop_2d(img2d, divisor):
    if divisor <= 1:
        return img2d
    h, w = img2d.shape[:2]
    ch, cw = h // divisor, w // divisor
    y0, x0 = (h - ch) // 2, (w - cw) // 2
    return img2d[y0 : y0 + ch, x0 : x0 + cw]


def process_images(input_folder, output_base, frame_order, crop_divisor, dry_run):
    tiff_files = sorted(
        f for f in os.listdir(input_folder) if f.lower().endswith((".tif", ".tiff"))
    )
    if not tiff_files:
        print("Error: No .tif/.tiff files found.")
        sys.exit(1)

    group_stats = {}
    processed = 0
    skipped = 0
    errors = 0

    print(f"\nFound {len(tiff_files)} TIFF files to process...")
    if dry_run:
        print("(DRY RUN — no files will be written)")
    print("=" * 60)

    for tiff_file in tiff_files:
        input_path = os.path.join(input_folder, tiff_file)
        base_name = os.path.splitext(tiff_file)[0]

        try:
            channels = parse_metadata(input_path)
            if not channels:
                print(f"  ⚠ {tiff_file}: no active lasers detected, skipping")
                skipped += 1
                continue

            v_group = group_name(channels)
            group_stats[v_group] = group_stats.get(v_group, 0) + 1
            cell_idx = group_stats[v_group]

            # Frame order for this TIFF.
            if frame_order is None:
                fo = list(range(len(channels)))
            elif len(frame_order) != len(channels):
                raise ValueError(
                    f"--frames has {len(frame_order)} entries but TIFF has "
                    f"{len(channels)} active channels"
                )
            else:
                fo = frame_order

            print(
                f"  {tiff_file}: lasers={channels_summary(channels)} "
                f"-> {v_group}/cell{cell_idx}/"
            )

            if dry_run:
                processed += 1
                continue

            cell_dir = os.path.join(output_base, v_group, f"cell{cell_idx}")
            os.makedirs(cell_dir, exist_ok=True)

            with tifffile.TiffFile(input_path) as tif:
                n_pages = len(tif.pages)
                bad_frames = [fi for fi in fo if fi < 0 or fi >= n_pages]
                if bad_frames:
                    raise IndexError(
                        f"Frames out of range: {bad_frames} (n_pages={n_pages})"
                    )

                for save_idx, frame_idx in enumerate(fo):
                    channel_num = save_idx + 1
                    channel_dir = os.path.join(cell_dir, f"ch{channel_num}")
                    os.makedirs(channel_dir, exist_ok=True)

                    frame = tif.pages[frame_idx].asarray()
                    img_to_save = center_crop_2d(frame, crop_divisor)

                    file_out = f"{base_name}_ch{channel_num}.tif"
                    tifffile.imwrite(
                        os.path.join(channel_dir, file_out), img_to_save
                    )

            processed += 1

        except Exception as e:
            print(f"  ✗ Error processing {tiff_file}: {e}")
            errors += 1

    print("=" * 60)
    if dry_run:
        print("\nDRY RUN COMPLETE — no files were written.")
    else:
        print("\nPROCESSING COMPLETE!")
    print(f"  ✓ Processed: {processed} files")
    if skipped:
        print(f"  ⚠ Skipped (no active lasers): {skipped} files")
    print(f"  ✗ Errors: {errors} files")
    if not dry_run:
        print(f"  Output saved to: {output_base}")

    if group_stats:
        print("\nGroups found:")
        for group, count in sorted(group_stats.items()):
            print(f"  • {group}: {count} cells")

    sys.exit(1 if errors > 0 else 0)


# ── CLI ────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Split TIFF channels by laser/voltage metadata "
                    "and prepare MATLAB input."
    )
    parser.add_argument("--input", required=True, help="Input folder containing TIFFs")
    parser.add_argument("--output", required=True, help="Base output folder")
    parser.add_argument(
        "--frames",
        default=None,
        help="Comma-separated frame indices, one per active channel in "
             "wavelength-ascending order (e.g. '0,1,2' for three channels). "
             "Optional — if omitted, defaults to 0,1,...,N-1 for each TIFF "
             "based on its active channel count.",
    )
    parser.add_argument(
        "--crop", type=int, default=1, help="Center crop divisor (1 = no crop)"
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Parse metadata and print what would be created, without "
             "writing any files. Use to verify metadata parsing first.",
    )
    args = parser.parse_args()

    if not os.path.isdir(args.input):
        print(f"Error: Input folder does not exist: {args.input}")
        sys.exit(1)

    if args.crop < 1:
        print("Error: --crop must be >= 1")
        sys.exit(1)

    if args.frames is None:
        frame_order = None
    else:
        try:
            frame_order = parse_int_list(args.frames)
            if not frame_order:
                print("Error: --frames must contain at least one index.")
                sys.exit(1)
        except ValueError as e:
            print(f"Error parsing --frames: {e}")
            sys.exit(1)

    if not args.dry_run:
        os.makedirs(args.output, exist_ok=True)

    print("=" * 60)
    print("IMAGE PROCESSOR — Generic Multi-Channel Splitter")
    print("=" * 60)
    print(f"Input folder:  {args.input}")
    print(f"Output folder: {args.output}")
    print(f"Frame order:   {frame_order if frame_order is not None else '(auto: 0,1,...,N-1 per TIFF)'}")
    print(f"Crop divisor:  {args.crop}")
    if args.dry_run:
        print("Mode:          DRY RUN")

    process_images(args.input, args.output, frame_order, args.crop, args.dry_run)


if __name__ == "__main__":
    main()
