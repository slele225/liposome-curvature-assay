"""
Overlay normalized curvature-sorting curves from multiple experiments.

Reads multiple filtered puncta A-value files (each with its own conversion
factor from dls_calibration.py), computes protein surface density vs. diameter
for each, bins the data, normalizes each curve so the rightmost bin
(largest diameter, lowest curvature) = 1, and overlays just the binned
averages on a single plot.

The y-axis shows fold-enrichment at high curvature relative to flat
membranes.
"""

import os
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt


EPS = 1e-12


# ── Data loading ────────────────────────────────────────────────────────

def load_puncta_file(path, lipid_col, protein_col):
    """Read lipid and protein A columns from the TSV."""
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

    for col in [lipid_col, protein_col]:
        if col not in headers:
            raise ValueError(f"Column '{col}' not found. Available: {headers}")

    li = headers.index(lipid_col)
    pi = headers.index(protein_col)

    lipid_A = np.array([float(r[li]) for r in rows])
    protein_A = np.array([float(r[pi]) for r in rows])
    return lipid_A, protein_A


# ── Core computation ────────────────────────────────────────────────────

def amplitude_to_diameter(lipid_A, conversion_factor):
    """
    Convert lipid amplitudes to physical diameters (nm).

    diameter = sqrt(lipid_A) * conversion_factor
    """
    lipid_A = np.clip(lipid_A, 0, None)
    sqrt_lipid = np.sqrt(lipid_A)
    return sqrt_lipid * conversion_factor


def compute_density(protein_A, diameter_nm):
    """Protein surface density = protein_A / (π D²)."""
    surface_area = np.pi * (diameter_nm ** 2)
    density = protein_A / (surface_area + EPS)

    valid = np.isfinite(density) & (surface_area > 0)
    return density, valid


def bin_by_diameter(x, y, bin_width):
    """Equal-width bins, returns centres and mean y."""
    edges = np.arange(np.min(x), np.max(x) + bin_width, bin_width)
    centres, means = [], []

    for i in range(len(edges) - 1):
        mask = (x >= edges[i]) & (x < edges[i + 1])
        if np.any(mask):
            centres.append(0.5 * (edges[i] + edges[i + 1]))
            means.append(np.mean(y[mask]))

    return np.array(centres), np.array(means)


def normalize_to_rightmost(centres, means):
    """Divide all bin means by the rightmost bin value so it equals 1."""
    rightmost_idx = np.argmax(centres)
    rightmost_val = means[rightmost_idx]

    if rightmost_val <= 0 or not np.isfinite(rightmost_val):
        raise ValueError(
            f"Rightmost bin value is {rightmost_val} -- cannot normalize."
        )

    return means / rightmost_val


# ── Input parsing ───────────────────────────────────────────────────────

def parse_input_pairs(input_args):
    """
    Parse input arguments as file:conversion_factor pairs.

    Accepts either:
      file1.txt:1.234 file2.txt:1.567
    or:
      file1.txt 1.234 file2.txt 1.567
    """
    pairs = []

    # Try colon-separated first
    if any(":" in a for a in input_args):
        for arg in input_args:
            if ":" not in arg:
                raise ValueError(
                    f"Expected file:conversion_factor format, got '{arg}'. "
                    f"Use path/to/file.txt:1.234"
                )
            path, factor_str = arg.rsplit(":", 1)
            pairs.append((path, float(factor_str)))
    else:
        # Alternating: file factor file factor
        if len(input_args) % 2 != 0:
            raise ValueError(
                "Expected alternating file/conversion_factor pairs. "
                "Use: file1.txt 1.234 file2.txt 1.567  OR  "
                "file1.txt:1.234 file2.txt:1.567"
            )
        for i in range(0, len(input_args), 2):
            pairs.append((input_args[i], float(input_args[i + 1])))

    return pairs


# ── CLI ─────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Overlay normalized curvature-sorting curves from "
                    "multiple experiments.",
        epilog=(
            "Example:\n"
            "  python plot_overlay.py \\\n"
            "    --input data/cond1/filtered.txt:1.234 "
            "data/cond2/filtered.txt:1.567 \\\n"
            "    --labels \"WT protein\" \"Mutant K58A\" \\\n"
            "    --save-dir figures/"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--input",
        required=True,
        nargs="+",
        help="Input files with conversion factors. Use file.txt:factor "
             "format (e.g., data/filtered.txt:1.234) or alternating "
             "file factor pairs. The conversion factor maps "
             "sqrt(A) to diameter in nm.",
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
        help="Bin width in nm (default: 0.5)",
    )
    parser.add_argument(
        "--labels",
        nargs="+",
        default=None,
        help="Custom labels for each curve (default: parent folder name)",
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
        default="figures/",
        help="Output directory for figure (default: figures/)",
    )
    parser.add_argument(
        "--output-name",
        default="normalized_curvature_overlay.png",
        help="Output filename (default: normalized_curvature_overlay.png)",
    )
    args = parser.parse_args()

    # Parse file:factor pairs
    try:
        pairs = parse_input_pairs(args.input)
    except ValueError as e:
        print(f"Error: {e}")
        sys.exit(1)

    n_curves = len(pairs)

    if args.labels and len(args.labels) != n_curves:
        print(
            f"Error: --labels has {len(args.labels)} entries but "
            f"there are {n_curves} input files."
        )
        sys.exit(1)

    # ── Process each file ───────────────────────────────────────────
    all_curves = []

    print("=" * 60)
    print("NORMALIZED CURVATURE OVERLAY")
    print("=" * 60)

    for i, (path, conv_factor) in enumerate(pairs):
        if not os.path.isfile(path):
            print(f"[SKIP] File not found: {path}")
            continue

        try:
            lipid_A, protein_A = load_puncta_file(
                path, args.lipid_col, args.protein_col
            )

            diameter = amplitude_to_diameter(lipid_A, conv_factor)
            density, valid = compute_density(protein_A, diameter)

            x = diameter[valid]
            y = density[valid]

            if args.diameter_cutoff is not None:
                cutoff_mask = x <= args.diameter_cutoff
                n_before = len(x)
                x = x[cutoff_mask]
                y = y[cutoff_mask]
                print(f"  Diameter cutoff at {args.diameter_cutoff} nm: "
                      f"{n_before} -> {len(x)} points")

            bx, by = bin_by_diameter(x, y, args.bin_width)
            by_norm = normalize_to_rightmost(bx, by)

            if args.labels:
                label = args.labels[i]
            else:
                label = os.path.basename(os.path.dirname(path))

            all_curves.append((bx, by_norm, label, conv_factor, len(x)))

            print(
                f"  {label}: {len(x)} points, conv={conv_factor:.4f}, "
                f"{len(bx)} bins, rightmost density = {by[-1]:.4g}"
            )

        except Exception as e:
            print(f"[ERROR] {path}: {e}")

    if not all_curves:
        print("No curves to plot.")
        sys.exit(1)

    # ── Plot ────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(10, 7))

    for bx, by_norm, label, conv_factor, n_pts in all_curves:
        ax.plot(
            bx, by_norm, "o-",
            markersize=5, linewidth=1.5,
            label=f"{label} (n={n_pts})",
        )

    ax.axhline(1.0, color="gray", linestyle="--", linewidth=0.8, alpha=0.5)
    ax.set_xlabel("Liposome diameter (nm)")
    ax.set_ylabel("Normalized protein density (fold enrichment)")
    ax.set_title("Curvature Sorting \u2014 Normalized Overlay")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.2)

    fig.tight_layout()

    os.makedirs(args.save_dir, exist_ok=True)
    save_path = os.path.join(args.save_dir, args.output_name)
    fig.savefig(save_path, dpi=300)
    plt.close(fig)
    print(f"\nSaved -> {save_path}")


if __name__ == "__main__":
    main()
