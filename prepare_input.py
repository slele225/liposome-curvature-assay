"""
Split multi-frame TIFFs into per-channel folders for MATLAB detection.

Groups images by Olympus FV3000 detector voltage metadata and organises
them into the cell/channel directory structure expected by CMEanalysis.
"""

import os
import sys
import argparse
import tifffile
import numpy as np


def parse_int_list(s: str):
    """Parse comma-separated list of nonnegative integers."""
    if not s.strip():
        return []
    parts = [p.strip() for p in s.split(",") if p.strip()]
    out = []
    for p in parts:
        if p.startswith("-") or not p.isdigit():
            raise ValueError(f"Invalid index: '{p}'")
        out.append(int(p))
    return out


def center_crop_2d(img2d: np.ndarray, divisor: int):
    """Center-crop a 2-D image. divisor=1 means no crop."""
    if divisor <= 1:
        return img2d
    h, w = img2d.shape[:2]
    ch, cw = h // divisor, w // divisor
    y0, x0 = (h - ch) // 2, (w - cw) // 2
    return img2d[y0 : y0 + ch, x0 : x0 + cw]


def get_voltages(tiff_path: str):
    """
    Parse Olympus FV3000 metadata for detector voltages.

    Detector #3 → 488 nm (typically Ch 0)
    Detector #5 → 561 nm (typically Ch 1)
    """
    with tifffile.TiffFile(tiff_path) as tif:
        ij_meta = tif.imagej_metadata
        if not ij_meta:
            return "Unknown_Voltage"

        info = ij_meta.get("Info", "")
        v_488, v_561 = "X", "X"

        for line in info.split("\n"):
            if "Detector voltage #3 =" in line:
                v_488 = line.split("=")[1].strip().split(".")[0]
            if "Detector voltage #5 =" in line:
                v_561 = line.split("=")[1].strip().split(".")[0]

        return f"488nm_{v_488}V_561nm_{v_561}V"


def process_images(input_folder, output_base, frame_order, crop_divisor):
    """Process TIFF files: split channels, reorder, crop, save."""
    tiff_files = sorted(
        f for f in os.listdir(input_folder) if f.lower().endswith((".tif", ".tiff"))
    )

    if not tiff_files:
        print("Error: No .tif/.tiff files found in the selected input folder.")
        sys.exit(1)

    group_stats = {}
    processed_count = 0
    error_count = 0

    print(f"\nFound {len(tiff_files)} TIFF files to process...")
    print("=" * 60)

    for tiff_file in tiff_files:
        input_path = os.path.join(input_folder, tiff_file)
        base_name = os.path.splitext(tiff_file)[0]

        try:
            v_group = get_voltages(input_path)
            condition_dir = os.path.join(output_base, v_group)
            os.makedirs(condition_dir, exist_ok=True)

            group_stats[v_group] = group_stats.get(v_group, 0) + 1
            cell_idx = group_stats[v_group]
            cell_dir = os.path.join(condition_dir, f"cell{cell_idx}")
            os.makedirs(cell_dir, exist_ok=True)

            with tifffile.TiffFile(input_path) as tif:
                n_pages = len(tif.pages)

                bad_frames = [fi for fi in frame_order if fi < 0 or fi >= n_pages]
                if bad_frames:
                    raise IndexError(
                        f"Frames out of range: {bad_frames} (n_pages={n_pages})"
                    )

                for save_idx, frame_idx in enumerate(frame_order):
                    channel_num = save_idx + 1
                    channel_dir = os.path.join(cell_dir, f"ch{channel_num}")
                    os.makedirs(channel_dir, exist_ok=True)

                    frame = tif.pages[frame_idx].asarray()
                    img_to_save = center_crop_2d(frame, crop_divisor)

                    file_out = f"{base_name}_ch{channel_num}.tif"
                    tifffile.imwrite(os.path.join(channel_dir, file_out), img_to_save)

            mapping = [f"ch{i+1}=frame{frame_order[i]}" for i in range(len(frame_order))]
            print(f"  ✓ [{v_group}] {tiff_file} -> cell{cell_idx} ({', '.join(mapping)})")
            processed_count += 1

        except Exception as e:
            print(f"  ✗ Error processing {tiff_file}: {e}")
            error_count += 1

    # --- summary ---
    print("=" * 60)
    print("\nPROCESSING COMPLETE!")
    print(f"  ✓ Successfully processed: {processed_count} files")
    print(f"  ✗ Errors: {error_count} files")
    print(f"  Output saved to: {output_base}")

    if group_stats:
        print("\nVoltage groups found:")
        for group, count in sorted(group_stats.items()):
            print(f"  • {group}: {count} cells")

    sys.exit(1 if error_count > 0 else 0)


def main():
    parser = argparse.ArgumentParser(
        description="Split TIFF channels, reorder frames, and prepare MATLAB input."
    )
    parser.add_argument("--input", required=True, help="Input folder containing TIFFs")
    parser.add_argument("--output", required=True, help="Base output folder")
    parser.add_argument(
        "--frames", required=True, help="Comma-separated frame order, e.g. 2,0,1"
    )
    parser.add_argument(
        "--crop", type=int, default=1, help="Center crop divisor (1 = no crop)"
    )
    args = parser.parse_args()

    if not os.path.isdir(args.input):
        print(f"Error: Input folder does not exist: {args.input}")
        sys.exit(1)

    os.makedirs(args.output, exist_ok=True)

    if args.crop < 1:
        print("Error: --crop must be >= 1")
        sys.exit(1)

    try:
        frame_order = parse_int_list(args.frames)
        if not frame_order:
            print("Error: No frames provided.")
            sys.exit(1)
    except ValueError as e:
        print(f"Error parsing --frames: {e}")
        sys.exit(1)

    # --- preview ---
    print("=" * 60)
    print("IMAGE PROCESSOR — Channel Splitter and Reorderer")
    print("=" * 60)
    print(f"Input folder:  {args.input}")
    print(f"Output folder: {args.output}")
    print(f"Frame order:   {frame_order}")
    print(f"Crop divisor:  {args.crop}")
    print("-" * 60)

    tiff_files = sorted(
        f for f in os.listdir(args.input) if f.lower().endswith((".tif", ".tiff"))
    )
    if not tiff_files:
        print("Error: No .tif/.tiff files found.")
        sys.exit(1)

    first_path = os.path.join(args.input, tiff_files[0])
    with tifffile.TiffFile(first_path) as tif:
        n_pages = len(tif.pages)
        probe = tif.pages[0].asarray()

    print(f"\nExample file: {tiff_files[0]}")
    print(f"  Frames available: {n_pages}")
    print(f"  First frame shape: {probe.shape}")

    bad_frames = [fi for fi in frame_order if fi < 0 or fi >= n_pages]
    if bad_frames:
        print(f"Error: Requested frames out of range: {bad_frames}")
        sys.exit(1)

    print("\nChannel mapping preview:")
    for i, frame_idx in enumerate(frame_order):
        print(f"  frame {frame_idx} → ch{i+1}/")

    print("\nStarting processing...\n")
    process_images(args.input, args.output, frame_order, args.crop)


if __name__ == "__main__":
    main()
