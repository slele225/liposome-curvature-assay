"""
DLS calibration — compute mean liposome diameter from a Malvern Zetasizer
Excel export.

The Zetasizer exports three distribution types stacked vertically in one
sheet: Intensity, Volume, and Number. Each has a header row followed by
diameter bins with percentage values for each measurement record.

This script finds the Number distribution section, computes the
weighted-average diameter for each record, then averages across records
to produce a single mean diameter for use in the curvature assay pipeline.

TODO (lognormal fitting):
  - Fit a lognormal to the DLS number distribution
  - Fit a lognormal to the sqrt(lipid_A) distribution from imaging
  - Use ratio of medians of the two fits as the conversion factor
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd


def find_section_rows(df, n_cols):
    """
    Find the starting rows for Intensity, Volume, and Number sections.

    The Zetasizer export has section headers like "X Intensity", "X Volume",
    "X Number" in column 0. Each section has 70 data rows after the header.
    """
    sections = {}
    for i in range(df.shape[0]):
        val = str(df.iloc[i, 0]).strip()
        if val.startswith("X "):
            section_name = val.replace("X ", "").strip().lower()
            sections[section_name] = i
    return sections


def extract_section(df, start_row, next_start=None):
    """
    Extract diameter bins and measurement data from one section.

    Returns (diameters, data_array) where data_array has shape
    (n_bins, n_records).
    """
    header = df.iloc[start_row]
    record_names = [str(header.iloc[c]) for c in range(1, df.shape[1])]

    if next_start is not None:
        end_row = next_start
    else:
        end_row = df.shape[0]

    data = df.iloc[start_row + 1 : end_row].copy()
    data = data.dropna(how="all")
    data = data.apply(pd.to_numeric, errors="coerce")
    data = data.dropna(how="all")

    diameters = data.iloc[:, 0].values.astype(float)
    values = data.iloc[:, 1:].values.astype(float)

    return diameters, values, record_names


def weighted_mean_diameter(diameters, weights):
    """Compute weighted mean diameter from a single distribution."""
    mask = (weights > 0) & np.isfinite(weights) & np.isfinite(diameters)
    if not np.any(mask):
        return np.nan
    d = diameters[mask]
    w = weights[mask]
    return float(np.sum(d * w) / np.sum(w))


def main():
    parser = argparse.ArgumentParser(
        description="Compute mean liposome diameter from Malvern Zetasizer "
                    "DLS export (.xlsx)."
    )
    parser.add_argument(
        "--input", required=True, help="Path to DLS .xlsx file"
    )
    parser.add_argument(
        "--sheet",
        default=0,
        help="Sheet name or index (default: first sheet)",
    )
    parser.add_argument(
        "--distribution",
        choices=["number", "volume", "intensity"],
        default="number",
        help="Which distribution to use (default: number)",
    )
    args = parser.parse_args()

    if not os.path.isfile(args.input):
        print(f"Error: file not found: {args.input}")
        sys.exit(1)

    # Handle sheet argument — try as int first
    sheet = args.sheet
    try:
        sheet = int(sheet)
    except ValueError:
        pass

    df = pd.read_excel(args.input, sheet_name=sheet, header=None)

    print("=" * 50)
    print("DLS CALIBRATION")
    print("=" * 50)
    print(f"File:         {args.input}")
    print(f"Sheet:        {sheet}")
    print(f"Distribution: {args.distribution}")
    print(f"Data shape:   {df.shape}")

    # Find sections
    sections = find_section_rows(df, df.shape[1])
    print(f"Sections found: {list(sections.keys())}")

    if args.distribution not in sections:
        print(
            f"\nError: '{args.distribution}' section not found. "
            f"Available: {list(sections.keys())}"
        )
        sys.exit(1)

    # Determine section boundaries
    sorted_sections = sorted(sections.items(), key=lambda x: x[1])
    target_idx = [i for i, (name, _) in enumerate(sorted_sections)
                  if name == args.distribution][0]
    start_row = sorted_sections[target_idx][1]
    next_start = (sorted_sections[target_idx + 1][1]
                  if target_idx + 1 < len(sorted_sections) else None)

    diameters, values, record_names = extract_section(df, start_row, next_start)

    n_records = values.shape[1]
    print(f"\nRecords: {n_records}")
    print(f"Diameter bins: {len(diameters)}")
    print(f"Diameter range: {diameters.min():.1f} – {diameters.max():.1f} nm")

    # Compute weighted mean for each record
    print(f"\nPer-record weighted mean diameters:")
    record_means = []
    for j in range(n_records):
        wm = weighted_mean_diameter(diameters, values[:, j])
        record_means.append(wm)
        name = record_names[j] if j < len(record_names) else f"Record {j+1}"
        if np.isfinite(wm):
            print(f"  {name}: {wm:.2f} nm")
        else:
            print(f"  {name}: no signal")

    # Average across records (excluding NaN)
    valid_means = [m for m in record_means if np.isfinite(m)]
    if not valid_means:
        print("\nError: no valid measurements found.")
        sys.exit(1)

    overall_mean = float(np.mean(valid_means))
    overall_std = float(np.std(valid_means)) if len(valid_means) > 1 else 0.0

    print(f"\n{'=' * 50}")
    print(f"Mean diameter ({args.distribution}, averaged across "
          f"{len(valid_means)} records): {overall_mean:.2f} nm")
    if len(valid_means) > 1:
        print(f"Std across records: {overall_std:.2f} nm")
    print(f"\n→ Use this in plot_curvature.py:")
    print(f"  --dls-mean-diameter {overall_mean:.2f}")


if __name__ == "__main__":
    main()
