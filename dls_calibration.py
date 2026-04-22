"""
DLS calibration via distribution overlay.

Finds a single scalar k that converts sqrt(fluorescence amplitude) to
physical diameter in nm by minimizing the chi-squared difference between
the rebinned fluorescence histogram and the DLS number-weighted size
distribution.

    diameter = sqrt(A) * k

This is the calibration procedure introduced by Kunding et al. (2008),
simplified by Hatzakis et al. (2009) as d = Ccal * sqrt(I), and used in
its current form by Zeno et al. (2018) and Johnson et al. (2025).

Also reports the ratio-of-means conversion as a quick sanity check.

References:
  Kunding et al. (2008) Biophys J 95:1176
  Hatzakis et al. (2009) Nat Chem Biol 5:835
  Bhatia et al. (2009) EMBO J 28:3303
  Zeno et al. (2018) Nat Commun 9:4152
  Johnson et al. (2025) Commun Biol 8:1179
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar


# ── Data loading ────────────────────────────────────────────────────────

def load_dls_number(xlsx_path, sheet=0):
    """
    Load and average the number distribution from a Zetasizer export.

    Finds the "X Number" section header, reads diameter bins and
    number percentages, averages across all measurement records.
    Returns (diameters, avg_number_pct).
    """
    df = pd.read_excel(xlsx_path, sheet_name=sheet, header=None)

    start_row = None
    for i in range(df.shape[0]):
        val = str(df.iloc[i, 0]).strip()
        if val.lower() == "x number":
            start_row = i
            break

    if start_row is None:
        raise ValueError("Could not find 'X Number' section in spreadsheet.")

    end_row = df.shape[0]
    for i in range(start_row + 1, df.shape[0]):
        val = str(df.iloc[i, 0]).strip()
        if val.startswith("X ") and val != "X Number":
            end_row = i
            break

    data = df.iloc[start_row + 1 : end_row].copy()
    data = data.apply(pd.to_numeric, errors="coerce").dropna(how="all")

    diameters = data.iloc[:, 0].values.astype(float)
    records = data.iloc[:, 1:].values.astype(float)
    avg_number = np.nanmean(records, axis=1)

    valid = np.isfinite(diameters) & np.isfinite(avg_number)
    return diameters[valid], avg_number[valid]


def load_sqrt_A(puncta_path, lipid_col="A_ch1"):
    """
    Load filtered puncta file, return raw sqrt(lipid_A) values.
    """
    headers = None
    rows = []

    with open(puncta_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if headers is None:
                headers = line.split("\t")
                continue
            rows.append(line.split("\t"))

    if headers is None:
        raise ValueError(f"No header row found in {puncta_path}")
    if lipid_col not in headers:
        raise ValueError(f"Column '{lipid_col}' not found. Available: {headers}")

    idx = headers.index(lipid_col)
    lipid_A = np.array([float(r[idx]) for r in rows])
    sqrt_A = np.sqrt(np.clip(lipid_A, 0, None))
    sqrt_A = sqrt_A[sqrt_A > 0]

    return sqrt_A


# ── DLS bin edges ───────────────────────────────────────────────────────

def dls_bin_edges(dls_diameters):
    """
    Infer bin edges from DLS bin centers (which are logarithmically spaced).

    For N bin centers, returns N+1 edges.  Uses geometric midpoints between
    consecutive centers, with first/last edges extrapolated.
    """
    d = dls_diameters
    n = len(d)

    edges = np.empty(n + 1)
    for i in range(n - 1):
        edges[i + 1] = np.sqrt(d[i] * d[i + 1])
    if n >= 2:
        ratio = d[1] / d[0]
        edges[0] = d[0] / np.sqrt(ratio)
        ratio = d[-1] / d[-2]
        edges[-1] = d[-1] * np.sqrt(ratio)
    else:
        edges[0] = d[0] * 0.5
        edges[-1] = d[0] * 1.5

    return edges


# ── Distribution overlay cost ──────────────────────────────────────────

def overlay_cost(k, sqrt_A, dls_edges, dls_density_norm, dls_widths):
    """
    Chi-squared cost for a candidate conversion factor k.

    Converts sqrt(A) to trial diameters via D = sqrt(A) * k,
    histograms them on the DLS bin grid, peak-normalizes both
    density distributions, and returns sum of squared residuals.
    """
    D_trial = sqrt_A * k

    fluor_counts, _ = np.histogram(D_trial, bins=dls_edges)
    fluor_density = fluor_counts.astype(float) / dls_widths

    fmax = fluor_density.max()
    if fmax == 0:
        return np.inf

    fluor_density_norm = fluor_density / fmax

    return float(np.sum((fluor_density_norm - dls_density_norm) ** 2))


# ── Ratio-of-means ─────────────────────────────────────────────────────

def ratio_of_means_conversion(dls_diameters, dls_weights, sqrt_A):
    """
    Quick conversion via mean(DLS number diameter) / mean(sqrt(A)).
    """
    mask = (dls_weights > 0) & np.isfinite(dls_weights) & np.isfinite(dls_diameters)
    d = dls_diameters[mask]
    w = dls_weights[mask]

    if np.sum(w) <= 0:
        return np.nan, np.nan, float(np.mean(sqrt_A))

    dls_mean = float(np.sum(d * w) / np.sum(w))
    fluor_mean = float(np.mean(sqrt_A))

    return dls_mean / fluor_mean, dls_mean, fluor_mean


# ── CLI ─────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="DLS calibration via distribution overlay."
    )
    parser.add_argument(
        "--dls-input", required=True,
        help="Path to Zetasizer .xlsx file",
    )
    parser.add_argument(
        "--fluor-input", required=True,
        help="Path to filtered_puncta_A_values.txt",
    )
    parser.add_argument(
        "--lipid-col", default="A_ch1",
        help='Column name for lipid amplitude (default: "A_ch1")',
    )
    parser.add_argument(
        "--save-dir", default=None,
        help="Directory to save overlay plot",
    )
    parser.add_argument(
        "--bootstrap", type=int, default=0,
        help="Number of bootstrap resamples for variance estimate "
             "(default: 0 = off).  Try 200-500.",
    )
    parser.add_argument(
        "--bootstrap-k", type=int, default=None,
        help="Number of puncta to sample per bootstrap iteration "
             "(default: same as total puncta count).",
    )
    args = parser.parse_args()

    # ── Load data ───────────────────────────────────────────────────
    print("=" * 60)
    print("DLS CALIBRATION — DISTRIBUTION OVERLAY")
    print("=" * 60)

    num_x, num_y = load_dls_number(args.dls_input)
    sqrt_A = load_sqrt_A(args.fluor_input, args.lipid_col)

    print(f"  DLS file:         {args.dls_input}")
    print(f"  Fluor file:       {args.fluor_input}")
    print(f"  DLS number bins:  {len(num_x)}")
    print(f"  Fluor puncta:     {len(sqrt_A)}")

    # ── Prepare DLS grid ──────────────────────────────────────────
    edges = dls_bin_edges(num_x)
    widths = np.diff(edges)

    # Convert DLS counts to density and peak-normalize
    dls_density = num_y / widths
    dls_density_norm = dls_density / dls_density.max()

    # ── Ratio-of-means (quick sanity check) ───────────────────────
    rom_factor, dls_mean, fluor_mean = ratio_of_means_conversion(
        num_x, num_y, sqrt_A
    )

    print(f"\n  Ratio-of-means sanity check:")
    print(f"    DLS number-weighted mean:  {dls_mean:.2f} nm")
    print(f"    Mean sqrt(A):              {fluor_mean:.4f}")
    print(f"    k_rom = {rom_factor:.6f}")

    # ── Distribution overlay optimization ─────────────────────────
    # Bracket k around the ratio-of-means estimate
    k_low = rom_factor * 0.3
    k_high = rom_factor * 3.0

    result = minimize_scalar(
        overlay_cost,
        bounds=(k_low, k_high),
        args=(sqrt_A, edges, dls_density_norm, widths),
        method="bounded",
        options={"xatol": 1e-10},
    )

    k_best = result.x

    # k_best IS the conversion factor: diameter = sqrt(A) * k_best
    mean_diameter_overlay = float(np.mean(sqrt_A)) * k_best

    print(f"\n{'=' * 60}")
    print("DISTRIBUTION OVERLAY RESULT")
    print(f"{'=' * 60}")
    print(f"  k (nm per sqrt(A)):          {k_best:.6f}")
    print(f"  Implied mean diameter:       {mean_diameter_overlay:.2f} nm")
    print(f"  Chi-squared residual:        {result.fun:.6f}")

    pct_diff = abs(k_best - rom_factor) / rom_factor * 100
    print(f"\n  Ratio-of-means factor:       {rom_factor:.6f}  "
          f"({pct_diff:.1f}% difference from overlay)")

    print(f"\n{'=' * 60}")
    print("USE IN PIPELINE")
    print(f"{'=' * 60}")
    print(f"\n  python plot_curvature.py --conversion-factor {k_best:.6f} ...")
    print(f"\n  Equivalently, using mean diameter:")
    print(f"  python plot_curvature.py --dls-mean-diameter {mean_diameter_overlay:.2f} ...")

    # ── Overlay plot ──────────────────────────────────────────────
    if args.save_dir:
        os.makedirs(args.save_dir, exist_ok=True)

        # Fluorescence histogram on DLS grid at k_best
        D_best = sqrt_A * k_best
        fluor_counts, _ = np.histogram(D_best, bins=edges)
        fluor_density = fluor_counts.astype(float) / widths
        fluor_density_norm = fluor_density / (fluor_density.max() or 1)

        fig, axes = plt.subplots(1, 2, figsize=(14, 5))

        # Left: overlay on DLS grid
        ax = axes[0]
        ax.step(num_x, dls_density_norm, where="mid", color="blue",
                linewidth=2, label="DLS number distribution")
        ax.step(num_x, fluor_density_norm, where="mid", color="red",
                linewidth=2, linestyle="--",
                label=f"Fluorescence (k={k_best:.4f})")
        ax.axvline(dls_mean, color="blue", linestyle=":", alpha=0.5,
                   label=f"DLS mean = {dls_mean:.1f} nm")
        ax.axvline(mean_diameter_overlay, color="red", linestyle=":",
                   alpha=0.5,
                   label=f"Fluor mean = {mean_diameter_overlay:.1f} nm")
        ax.set_xlabel("Diameter (nm)")
        ax.set_ylabel("Normalized density")
        ax.set_title("Distribution Overlay")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.2)

        # Right: raw sqrt(A) histogram with nm axis
        ax = axes[1]
        hist_bins = 80
        counts_raw, edges_raw = np.histogram(sqrt_A, bins=hist_bins)
        centers_raw = 0.5 * (edges_raw[:-1] + edges_raw[1:])
        ax.bar(centers_raw, counts_raw,
               width=centers_raw[1] - centers_raw[0],
               alpha=0.6, color="gray", label="sqrt(A) histogram")
        ax.set_xlabel("sqrt(Lipid Amplitude)")
        ax.set_ylabel("Count")
        ax.set_title("Fluorescence sqrt(A) Distribution")

        # Second x-axis for nm
        ax2 = ax.twiny()
        ax2.set_xlim(ax.get_xlim()[0] * k_best, ax.get_xlim()[1] * k_best)
        ax2.set_xlabel("Diameter (nm)")

        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.2)

        fig.tight_layout()
        save_path = os.path.join(args.save_dir, "dls_overlay_calibration.png")
        fig.savefig(save_path, dpi=300)
        plt.close(fig)
        print(f"\nPlot saved to {save_path}")

    # ── Bootstrap ─────────────────────────────────────────────────
    if args.bootstrap > 0:
        n_boot = args.bootstrap
        k_per = args.bootstrap_k or len(sqrt_A)
        rng = np.random.default_rng(42)

        print(f"\n{'=' * 60}")
        print(f"BOOTSTRAP VARIANCE (n={n_boot}, k={k_per})")
        print(f"{'=' * 60}")

        overlay_factors = []
        rom_factors_boot = []

        for b in range(n_boot):
            sample = rng.choice(sqrt_A, size=k_per, replace=True)

            # Ratio-of-means
            rom_b = dls_mean / np.mean(sample)
            rom_factors_boot.append(rom_b)

            # Distribution overlay
            k_lo_b = rom_b * 0.3
            k_hi_b = rom_b * 3.0
            res_b = minimize_scalar(
                overlay_cost,
                bounds=(k_lo_b, k_hi_b),
                args=(sample, edges, dls_density_norm, widths),
                method="bounded",
                options={"xatol": 1e-10},
            )
            overlay_factors.append(res_b.x)

            if (b + 1) % 50 == 0 or b == 0:
                print(f"  Bootstrap {b+1}/{n_boot}...")

        overlay_factors = np.array(overlay_factors)
        rom_factors_boot = np.array(rom_factors_boot)

        ov_cv = np.std(overlay_factors) / np.mean(overlay_factors) * 100
        rom_cv = np.std(rom_factors_boot) / np.mean(rom_factors_boot) * 100

        print(f"\n  Overlay:         mean={np.mean(overlay_factors):.4f}, "
              f"std={np.std(overlay_factors):.4f}, CV={ov_cv:.2f}%")
        print(f"  Ratio-of-means:  mean={np.mean(rom_factors_boot):.4f}, "
              f"std={np.std(rom_factors_boot):.4f}, CV={rom_cv:.2f}%")

        if ov_cv < rom_cv:
            print(f"\n  -> Overlay has lower CV ({ov_cv:.2f}% vs {rom_cv:.2f}%)")
        else:
            print(f"\n  -> Ratio-of-means has lower CV "
                  f"({rom_cv:.2f}% vs {ov_cv:.2f}%)")

        # Bootstrap plot
        if args.save_dir:
            fig, ax = plt.subplots(figsize=(10, 5))

            ax.hist(overlay_factors, bins=40, alpha=0.5, color="red",
                    label=f"Overlay (CV={ov_cv:.1f}%)", density=True)
            ax.hist(rom_factors_boot, bins=40, alpha=0.5, color="blue",
                    label=f"Ratio-of-means (CV={rom_cv:.1f}%)", density=True)

            ax.axvline(np.mean(overlay_factors), color="red",
                       linestyle="--", linewidth=1.5)
            ax.axvline(np.mean(rom_factors_boot), color="blue",
                       linestyle="--", linewidth=1.5)

            ax.set_xlabel("Conversion factor (nm / sqrt(A))")
            ax.set_ylabel("Density")
            ax.set_title(f"Bootstrap Comparison (n={n_boot})")
            ax.legend(fontsize=9)
            ax.grid(True, alpha=0.2)
            fig.tight_layout()

            boot_path = os.path.join(args.save_dir, "bootstrap_variance.png")
            fig.savefig(boot_path, dpi=300)
            plt.close(fig)
            print(f"\n  Bootstrap plot saved to {boot_path}")


if __name__ == "__main__":
    main()
