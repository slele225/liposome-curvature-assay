"""
Plot histograms of puncta amplitude values and derived size distributions.

Reads the filtered puncta A-value file produced by analyze_matlab.py and
generates:
  1. Raw amplitude histograms for each channel
  2. Estimated diameter distribution (sqrt(lipid_A) scaled by DLS calibration)
     optionally overlaid with the DLS distribution for visual comparison
"""

import os
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt


EPS = 1e-12


# ── Data loading ────────────────────────────────────────────────────────

def load_puncta_file(path: str, columns: list[str]):
    """
    Read specified columns from the TSV produced by analyze_matlab.py.

    Returns a dict of column_name -> 1-D numpy array.
    """
    headers = None
    rows = []

    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if headers is None:
                headers = line.split("\t")
                continue
            rows.append(line.split("\t"))

    if headers is None:
        raise ValueError(f"No header row found in {path}")

    col_indices = {}
    for col in columns:
        if col not in headers:
            raise ValueError(f"Column '{col}' not found. Available: {headers}")
        col_indices[col] = headers.index(col)

    result = {}
    for col, idx in col_indices.items():
        result[col] = np.array([float(r[idx]) for r in rows])

    return result


# ── Plotting ────────────────────────────────────────────────────────────

def apply_transform(values: np.ndarray, transform: str):
    """Apply amplitude transform and return (transformed_values, label)."""
    if transform == "sqrt":
        out = np.sqrt(np.clip(values, 0, None))
        label = "√A"
    elif transform == "log_sqrt":
        clipped = np.clip(values, EPS, None)
        out = np.log(np.sqrt(clipped))
        label = "log(√A)"
    else:
        out = values
        label = "A"
    return out, label


def plot_raw_amplitudes(data: dict, save_path: str, bins: int, transform: str = "raw"):
    """Histogram of A values for each channel, with optional transform."""
    n_cols = len(data)
    fig, axes = plt.subplots(1, n_cols, figsize=(6 * n_cols, 5))
    if n_cols == 1:
        axes = [axes]

    for ax, (col_name, values) in zip(axes, data.items()):
        transformed, t_label = apply_transform(values, transform)

        ax.hist(transformed, bins=bins, alpha=0.7, edgecolor="black", linewidth=0.3)
        ax.set_xlabel(f"{col_name}  ({t_label})")
        ax.set_ylabel("Count")
        ax.set_title(f"Distribution of {col_name}\n({t_label}, n = {len(values)})")

        mean_val = np.mean(transformed)
        median_val = np.median(transformed)
        ax.axvline(mean_val, color="red", linestyle="--", linewidth=1,
                    label=f"Mean = {mean_val:.4g}")
        ax.axvline(median_val, color="orange", linestyle="--", linewidth=1,
                    label=f"Median = {median_val:.4g}")
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.2)

    transform_names = {"raw": "Raw", "sqrt": "√A", "log_sqrt": "log(√A)"}
    fig.suptitle(f"Amplitude Distributions ({transform_names[transform]})",
                 fontsize=14, y=1.02)
    fig.tight_layout()
    fig.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved amplitude histograms ({transform}) → {save_path}")


def plot_diameter_distribution(lipid_A: np.ndarray, dls_mean_diameter: float,
                                save_path: str, bins: int):
    """
    Histogram of estimated diameters from lipid amplitudes, calibrated
    against DLS mean diameter.
    """
    lipid_A = np.clip(lipid_A, 0, None)
    sqrt_lipid = np.sqrt(lipid_A)
    mean_sqrt = float(np.mean(sqrt_lipid))

    if mean_sqrt <= 0:
        raise ValueError("mean(sqrt(lipid_A)) is non-positive — cannot calibrate.")

    scale = dls_mean_diameter / mean_sqrt
    diameters = sqrt_lipid * scale

    fig, ax = plt.subplots(figsize=(8, 5))

    ax.hist(diameters, bins=bins, alpha=0.7, edgecolor="black", linewidth=0.3,
            density=True, label="Estimated from imaging")

    mean_d = np.mean(diameters)
    median_d = np.median(diameters)
    std_d = np.std(diameters)

    ax.axvline(mean_d, color="red", linestyle="--", linewidth=1,
               label=f"Mean = {mean_d:.1f} nm")
    ax.axvline(median_d, color="orange", linestyle="--", linewidth=1,
               label=f"Median = {median_d:.1f} nm")
    ax.axvline(dls_mean_diameter, color="green", linestyle="-", linewidth=1.5,
               label=f"DLS mean = {dls_mean_diameter:.1f} nm")

    ax.set_xlabel("Estimated liposome diameter (nm)")
    ax.set_ylabel("Density")
    ax.set_title(
        f"Estimated Diameter Distribution (n = {len(diameters)})\n"
        f"Mean = {mean_d:.1f} nm, Std = {std_d:.1f} nm"
    )
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.2)

    fig.tight_layout()
    fig.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved diameter distribution → {save_path}")


def plot_protein_per_liposome(lipid_A: np.ndarray, protein_A: np.ndarray,
                               save_path: str, bins: int):
    """Histogram of protein amplitude per punctum."""
    fig, ax = plt.subplots(figsize=(8, 5))

    ax.hist(protein_A, bins=bins, alpha=0.7, edgecolor="black", linewidth=0.3)

    mean_val = np.mean(protein_A)
    median_val = np.median(protein_A)
    ax.axvline(mean_val, color="red", linestyle="--", linewidth=1,
               label=f"Mean = {mean_val:.2f}")
    ax.axvline(median_val, color="orange", linestyle="--", linewidth=1,
               label=f"Median = {median_val:.2f}")

    ax.set_xlabel("Protein channel amplitude")
    ax.set_ylabel("Count")
    ax.set_title(f"Protein Amplitude Distribution (n = {len(protein_A)})")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.2)

    fig.tight_layout()
    fig.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved protein amplitude histogram → {save_path}")


# ── CLI ─────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Plot histograms of puncta amplitudes and estimated size distributions."
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Path to filtered_puncta_A_values.txt from analyze_matlab.py",
    )
    parser.add_argument(
        "--lipid-col",
        default="A_ch1",
        help='Column name for lipid amplitude (default: "A_ch1")',
    )
    parser.add_argument(
        "--protein-col",
        default=None,
        help='Column name for protein amplitude (default: None). '
             'Omit for lipid-only analysis.',
    )
    parser.add_argument(
        "--dls-mean-diameter",
        type=float,
        default=None,
        help="Mean liposome diameter (nm) from DLS. If provided, generates "
             "a calibrated diameter distribution plot.",
    )
    parser.add_argument(
        "--bins",
        type=int,
        default=80,
        help="Number of histogram bins (default: 80)",
    )
    parser.add_argument(
        "--transform",
        choices=["raw", "sqrt", "log_sqrt"],
        default="raw",
        help="Transform for amplitude histograms: "
             "'raw' = plot A directly, "
             "'sqrt' = plot sqrt(A), "
             "'log_sqrt' = plot log(sqrt(A)). "
             "(default: raw)",
    )
    parser.add_argument(
        "--save-dir",
        default=None,
        help="Output directory for figures (default: same folder as input)",
    )
    args = parser.parse_args()

    if not os.path.isfile(args.input):
        print(f"Error: file not found: {args.input}")
        sys.exit(1)

    columns = [args.lipid_col]
    if args.protein_col:
        columns.append(args.protein_col)

    data = load_puncta_file(args.input, columns)

    lipid_A = data[args.lipid_col]
    protein_A = data.get(args.protein_col) if args.protein_col else None

    print(f"File: {args.input}")
    print(f"Points: {len(lipid_A)}")
    if protein_A is None:
        print("Mode: lipid-only (no protein channel)")

    # Output paths
    parent = os.path.basename(os.path.dirname(args.input))
    save_dir = args.save_dir or os.path.dirname(args.input)
    os.makedirs(save_dir, exist_ok=True)

    # 1. Amplitude histograms (with optional transform)
    suffix = {"raw": "amplitude", "sqrt": "sqrt_amplitude", "log_sqrt": "log_sqrt_amplitude"}
    raw_path = os.path.join(save_dir, f"{parent}__{suffix[args.transform]}_histograms.png")
    plot_raw_amplitudes(data, raw_path, args.bins, args.transform)

    # 2. Protein amplitude histogram (only if protein column provided)
    if protein_A is not None:
        prot_path = os.path.join(save_dir, f"{parent}__protein_amplitude_histogram.png")
        plot_protein_per_liposome(lipid_A, protein_A, prot_path, args.bins)

    # 3. Diameter distribution (only if DLS diameter provided)
    if args.dls_mean_diameter is not None:
        diam_path = os.path.join(save_dir, f"{parent}__diameter_distribution.png")
        plot_diameter_distribution(lipid_A, args.dls_mean_diameter, diam_path, args.bins)
    else:
        print("\nSkipping diameter distribution (no --dls-mean-diameter provided).")
        print("Add --dls-mean-diameter to see estimated liposome sizes.")


if __name__ == "__main__":
    main()
