"""
Plot protein surface density vs. liposome diameter (curvature sorting).

Reads the filtered puncta A-value file produced by analyze_matlab.py,
converts lipid amplitude → physical diameter using a DLS-derived mean
diameter, then plots protein_A / (πD²) vs D with binned averages.
"""

import os
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt


EPS = 1e-12  # avoid division by zero


# ── Data loading ────────────────────────────────────────────────────────

def load_puncta_file(path: str, lipid_col: str, protein_col: str):
    """
    Read the TSV produced by analyze_matlab.py.

    Returns lipid_A and protein_A as 1-D arrays.
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

    if lipid_col not in headers:
        raise ValueError(
            f"Column '{lipid_col}' not found. Available: {headers}"
        )
    if protein_col not in headers:
        raise ValueError(
            f"Column '{protein_col}' not found. Available: {headers}"
        )

    li = headers.index(lipid_col)
    pi = headers.index(protein_col)

    lipid_A = np.array([float(r[li]) for r in rows])
    protein_A = np.array([float(r[pi]) for r in rows])
    return lipid_A, protein_A


# ── Core computation ────────────────────────────────────────────────────

def amplitude_to_diameter(lipid_A: np.ndarray, dls_mean_diameter_nm: float = None,
                         conversion_factor: float = None):
    """
    Convert lipid amplitudes to physical diameters (nm).

    Two modes:
      1. --dls-mean-diameter: ratio-of-means method
         scale = dls_mean_diameter / mean(sqrt(lipid_A))
         diameter = sqrt(lipid_A) * scale

      2. --conversion-factor: direct mapping from sqrt(A) to diameter
         diameter = sqrt(lipid_A) * conversion_factor

    Returns (diameter_nm, scale_factor).
    """
    lipid_A = np.clip(lipid_A, 0, None)
    sqrt_lipid = np.sqrt(lipid_A)

    if conversion_factor is not None:
        scale = conversion_factor
    elif dls_mean_diameter_nm is not None:
        mean_sqrt = float(np.mean(sqrt_lipid))
        if mean_sqrt <= 0:
            raise ValueError("mean(sqrt(lipid_A)) is non-positive.")
        scale = dls_mean_diameter_nm / mean_sqrt
    else:
        raise ValueError("Must provide either dls_mean_diameter or conversion_factor.")

    diameter_nm = sqrt_lipid * scale
    return diameter_nm, scale


def compute_protein_density(protein_A: np.ndarray, diameter_nm: np.ndarray):
    """Protein surface density = protein_A / (π D²)."""
    surface_area = np.pi * (diameter_nm ** 2)
    density = protein_A / (surface_area + EPS)

    valid = np.isfinite(density) & (surface_area > 0)
    return density, valid


def bin_by_diameter(x: np.ndarray, y: np.ndarray, bin_width: float):
    """Equal-width bins along x, returning bin centres and mean y."""
    edges = np.arange(np.min(x), np.max(x) + bin_width, bin_width)
    centres, means = [], []

    for i in range(len(edges) - 1):
        mask = (x >= edges[i]) & (x < edges[i + 1])
        if np.any(mask):
            centres.append(0.5 * (edges[i] + edges[i + 1]))
            means.append(np.mean(y[mask]))

    return np.array(centres), np.array(means)


# ── Plotting ────────────────────────────────────────────────────────────

def make_plot(x, y, bx, by, title, save_path, bin_width, y_pad=0.3):
    fig, ax = plt.subplots(figsize=(10, 7))

    ax.scatter(x, y, s=6, alpha=0.08, label="Individual puncta")
    ax.scatter(bx, by, s=30, zorder=5, label=f"{bin_width} nm bin mean")

    ax.set_xlabel("Liposome diameter (nm)")
    ax.set_ylabel("Protein surface density:  A / (πD²)  [A per nm²]")
    ax.set_title(title)
    ax.grid(True, alpha=0.2)
    ax.legend()

    # Fit y-axis to bin means so the trend isn't squashed by scatter outliers
    if len(by) > 0:
        pad = (by.max() - by.min()) * y_pad
        if pad == 0:
            pad = by.max() * y_pad
        ax.set_ylim(by.min() - pad, by.max() + pad)

    fig.tight_layout()
    fig.savefig(save_path, dpi=300)
    plt.close(fig)
    print(f"Saved figure → {save_path}")


# ── CLI ─────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Plot protein curvature sorting from filtered puncta A values."
    )
    parser.add_argument(
        "--input",
        required=True,
        nargs="+",
        help="Path(s) to filtered_puncta_A_values.txt files",
    )
    parser.add_argument(
        "--dls-mean-diameter",
        type=float,
        default=None,
        help="Mean liposome diameter (nm) from DLS. Uses ratio-of-means "
             "conversion. Provide this OR --conversion-factor, not both.",
    )
    parser.add_argument(
        "--conversion-factor",
        type=float,
        default=None,
        help="Direct conversion factor: diameter_nm = sqrt(A) * factor. "
             "From dls_calibration.py. "
             "Provide this OR --dls-mean-diameter, not both.",
    )
    parser.add_argument(
        "--lipid-col",
        default="A_ch1",
        help='Column name for lipid amplitude (default: "A_ch1")',
    )
    parser.add_argument(
        "--protein-col",
        default="A_ch2",
        help='Column name for protein amplitude (default: "A_ch2")',
    )
    parser.add_argument(
        "--bin-width",
        type=float,
        default=0.5,
        help="Diameter bin width in nm for averaged curve (default: 0.5)",
    )
    parser.add_argument(
        "--y-pad",
        type=float,
        default=0.3,
        help="Y-axis padding factor above/below the bin means range "
             "(default: 0.3 = 30%%).",
    )
    parser.add_argument(
        "--diameter-cutoff",
        type=float,
        default=None,
        help="Exclude puncta with diameter above this value (nm). "
             "Useful for filtering sparse large-diameter tail.",
    )
    parser.add_argument(
        "--save-dir",
        default=None,
        help="Output directory for figures (default: same folder as input)",
    )
    args = parser.parse_args()

    if args.dls_mean_diameter is None and args.conversion_factor is None:
        parser.error("Must provide either --dls-mean-diameter or --conversion-factor.")
    if args.dls_mean_diameter is not None and args.conversion_factor is not None:
        parser.error("Provide only one of --dls-mean-diameter or --conversion-factor.")

    for path in args.input:
        if not os.path.isfile(path):
            print(f"[SKIP] File not found: {path}")
            continue

        try:
            lipid_A, protein_A = load_puncta_file(
                path, args.lipid_col, args.protein_col
            )

            diameter, scale = amplitude_to_diameter(
                lipid_A,
                dls_mean_diameter_nm=args.dls_mean_diameter,
                conversion_factor=args.conversion_factor,
            )
            density, valid = compute_protein_density(protein_A, diameter)

            x = diameter[valid]
            y = density[valid]

            if args.diameter_cutoff is not None:
                cutoff_mask = x <= args.diameter_cutoff
                n_before = len(x)
                x = x[cutoff_mask]
                y = y[cutoff_mask]
                print(f"  Diameter cutoff at {args.diameter_cutoff} nm: "
                      f"{n_before} -> {len(x)} points")

            print(f"\nFile: {path}")
            print(f"  Scale factor: {scale:.6f} nm / sqrt(A)")
            print(f"  Points plotted: {len(x)}")

            bx, by = bin_by_diameter(x, y, args.bin_width)

            # Build output path
            parent = os.path.basename(os.path.dirname(path))
            save_dir = args.save_dir or os.path.dirname(path)
            os.makedirs(save_dir, exist_ok=True)

            save_name = f"{parent}__protein_density_vs_diameter.png"
            save_path = os.path.join(save_dir, save_name)

            title = f"Protein surface density vs diameter\n{parent}"
            make_plot(x, y, bx, by, title, save_path, args.bin_width,
                      y_pad=args.y_pad)

        except Exception as e:
            print(f"[ERROR] {path}: {e}")


if __name__ == "__main__":
    main()
