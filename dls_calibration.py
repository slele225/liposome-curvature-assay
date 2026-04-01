"""
DLS calibration with log-normal fitting.

Simultaneously fits log-normal distributions to:
  1. The DLS intensity-weighted size distribution (from Zetasizer xlsx)
  2. The sqrt(lipid amplitude) distribution (from fluorescence data)
  3. A consistency term linking the two fits via a conversion factor

The conversion factor maps sqrt(A) units to nm:
  conversion = exp(mu_dls - mu_fluor)

This also computes the simple ratio-of-means conversion for comparison.

Loss = L1 + L2 + L3 where:
  L1 = sum( (lognormal_dls(d) - empirical_dls(d))^2 )
  L2 = sum( (lognormal_fluor(x) - empirical_fluor(x))^2 )
  L3 = sum( (lognormal_dls(d) - lognormal_fluor(d / conversion))^2 )
       evaluated on a shared nm axis
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import minimize


# ── Log-normal model ────────────────────────────────────────────────────

def lognormal(x, A, mu, sigma):
    """
    Log-normal PDF (unnormalized, for curve fitting).

    f(x) = A * exp(-0.5 * ((ln(x) - mu) / sigma)^2)

    Parameters: A (amplitude), mu (log-mean), sigma (log-std).
    Mode = exp(mu - sigma^2).
    """
    with np.errstate(divide="ignore", invalid="ignore"):
        logx = np.log(np.clip(x, 1e-12, None))
        return A * np.exp(-0.5 * ((logx - mu) / sigma) ** 2)


def lognormal_mode(mu, sigma):
    """Mode of a log-normal distribution: exp(mu - sigma^2)."""
    return np.exp(mu - sigma ** 2)


# ── Data loading ────────────────────────────────────────────────────────

def load_dls_intensity(xlsx_path, sheet=0):
    """
    Load and average the intensity distribution from a Zetasizer export.

    Finds the "X Intensity" section header, reads diameter bins and
    intensity percentages, averages across all measurement records.
    """
    df = pd.read_excel(xlsx_path, sheet_name=sheet, header=None)

    # Find intensity section
    start_row = None
    for i in range(df.shape[0]):
        val = str(df.iloc[i, 0]).strip()
        if val.lower() == "x intensity":
            start_row = i
            break

    if start_row is None:
        raise ValueError("Could not find 'X Intensity' section in spreadsheet.")

    # Find next section header to determine end
    end_row = df.shape[0]
    for i in range(start_row + 1, df.shape[0]):
        val = str(df.iloc[i, 0]).strip()
        if val.startswith("X ") and val != "X Intensity":
            end_row = i
            break

    data = df.iloc[start_row + 1 : end_row].copy()
    data = data.apply(pd.to_numeric, errors="coerce").dropna(how="all")

    diameters = data.iloc[:, 0].values.astype(float)
    records = data.iloc[:, 1:].values.astype(float)

    # Average across records
    avg_intensity = np.nanmean(records, axis=1)

    # Keep only bins with finite values
    valid = np.isfinite(diameters) & np.isfinite(avg_intensity)
    return diameters[valid], avg_intensity[valid]


def load_dls_number(xlsx_path, sheet=0):
    """
    Load and average the number distribution from a Zetasizer export.

    Finds the "X Number" section header, reads diameter bins and
    number percentages, averages across all measurement records.
    """
    df = pd.read_excel(xlsx_path, sheet_name=sheet, header=None)

    # Find number section
    start_row = None
    for i in range(df.shape[0]):
        val = str(df.iloc[i, 0]).strip()
        if val.lower() == "x number":
            start_row = i
            break

    if start_row is None:
        raise ValueError("Could not find 'X Number' section in spreadsheet.")

    # Find next section header to determine end
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

    # Average across records
    avg_number = np.nanmean(records, axis=1)

    # Keep only bins with finite values
    valid = np.isfinite(diameters) & np.isfinite(avg_number)
    return diameters[valid], avg_number[valid]


def load_sqrt_A(puncta_path, lipid_col="A_ch1", n_bins=80):
    """
    Load filtered puncta file, compute sqrt(lipid_A), and histogram it.

    Returns (bin_centers, counts) for the sqrt(A) distribution.
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

    # Remove zeros and outliers
    sqrt_A = sqrt_A[sqrt_A > 0]

    # Histogram
    counts, edges = np.histogram(sqrt_A, bins=n_bins)
    centers = 0.5 * (edges[:-1] + edges[1:])

    return centers, counts.astype(float), sqrt_A


# ── Loss functions ──────────────────────────────────────────────────────

def total_loss(params, dls_x, dls_y, fluor_x, fluor_y):
    """
    Combined loss from all three terms.

    params = [A_dls, mu_dls, sigma_dls, A_fluor, mu_fluor, sigma_fluor]
    """
    A_dls, mu_dls, sigma_dls, A_fluor, mu_fluor, sigma_fluor = params

    # Enforce positive sigma
    if sigma_dls <= 0 or sigma_fluor <= 0:
        return 1e12
    if A_dls <= 0 or A_fluor <= 0:
        return 1e12

    # L1: DLS lognormal fit to empirical DLS
    pred_dls = lognormal(dls_x, A_dls, mu_dls, sigma_dls)
    L1 = np.sum((pred_dls - dls_y) ** 2)

    # L2: Fluorescence lognormal fit to empirical sqrt(A)
    pred_fluor = lognormal(fluor_x, A_fluor, mu_fluor, sigma_fluor)
    L2 = np.sum((pred_fluor - fluor_y) ** 2)

    # L3: Consistency between the two lognormals on a shared nm axis
    # conversion maps fluor x-axis to nm: d_nm = x_fluor * conversion
    # conversion = exp(mu_dls) / exp(mu_fluor) = exp(mu_dls - mu_fluor)
    conversion = np.exp(mu_dls - mu_fluor)

    # Evaluate both lognormals on the DLS x-axis (nm)
    dls_on_shared = lognormal(dls_x, A_dls, mu_dls, sigma_dls)

    # Convert DLS x-axis to fluor units, evaluate fluor lognormal, then
    # scale amplitude to match DLS space
    fluor_x_converted = dls_x / conversion
    fluor_on_shared = lognormal(fluor_x_converted, A_fluor, mu_fluor, sigma_fluor)

    # Normalize both to unit area for shape comparison
    dls_norm = dls_on_shared / (np.sum(dls_on_shared) + 1e-12)
    fluor_norm = fluor_on_shared / (np.sum(fluor_on_shared) + 1e-12)

    L3 = np.sum((dls_norm - fluor_norm) ** 2)

    return L1 + L2 + L3


# ── Simple ratio-of-means (from old dls_calibration.py) ─────────────────

def ratio_of_means_conversion(num_diameters, num_weights, sqrt_A_raw):
    """
    Simple conversion factor using number-weighted mean diameter / mean(sqrt(A)).

    Uses the DLS number distribution (not intensity) for the weighted mean.
    """
    mask = (num_weights > 0) & np.isfinite(num_weights) & np.isfinite(num_diameters)
    d = num_diameters[mask]
    w = num_weights[mask]

    if np.sum(w) <= 0:
        return np.nan, np.nan, float(np.mean(sqrt_A_raw))

    dls_mean = float(np.sum(d * w) / np.sum(w))
    fluor_mean = float(np.mean(sqrt_A_raw))

    return dls_mean / fluor_mean, dls_mean, fluor_mean


# ── CLI ─────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="DLS calibration with log-normal fitting."
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
        "--z-avg", type=float, required=True,
        help="Z-average diameter (nm) from Zetasizer cumulants analysis",
    )
    parser.add_argument(
        "--pdi", type=float, required=True,
        help="PDI from Zetasizer cumulants analysis",
    )
    parser.add_argument(
        "--lipid-col", default="A_ch1",
        help='Column name for lipid amplitude (default: "A_ch1")',
    )
    parser.add_argument(
        "--n-bins", type=int, default=80,
        help="Number of bins for sqrt(A) histogram (default: 80)",
    )
    parser.add_argument(
        "--save-dir", default=None,
        help="Directory to save overlay plots (default: no plots)",
    )
    args = parser.parse_args()

    # ── Load data ───────────────────────────────────────────────────
    print("=" * 60)
    print("DLS CALIBRATION — LOG-NORMAL FITTING")
    print("=" * 60)

    dls_x, dls_y = load_dls_intensity(args.dls_input)
    num_x, num_y = load_dls_number(args.dls_input)
    fluor_x, fluor_y, sqrt_A_raw = load_sqrt_A(
        args.fluor_input, args.lipid_col, args.n_bins
    )

    print(f"DLS file:       {args.dls_input}")
    print(f"Fluor file:     {args.fluor_input}")
    print(f"Z-average:      {args.z_avg} nm")
    print(f"PDI:            {args.pdi}")
    print(f"DLS intensity bins: {len(dls_x)}")
    print(f"DLS number bins:    {len(num_x)}")
    print(f"Fluor bins:     {len(fluor_x)}")
    print(f"Fluor points:   {len(sqrt_A_raw)}")

    # ── Initial guesses ─────────────────────────────────────────────
    mu_dls_init = np.log(args.z_avg)
    sigma_dls_init = np.sqrt(args.pdi) if args.pdi > 0 else 0.3
    A_dls_init = float(np.max(dls_y))

    mu_fluor_init = float(np.log(np.mean(sqrt_A_raw)))
    sigma_fluor_init = float(np.std(np.log(sqrt_A_raw)))
    A_fluor_init = float(np.max(fluor_y))

    p0 = [
        A_dls_init, mu_dls_init, sigma_dls_init,
        A_fluor_init, mu_fluor_init, sigma_fluor_init,
    ]

    print(f"\nInitial guesses:")
    print(f"  DLS:   A={A_dls_init:.2f}, mu={mu_dls_init:.4f}, "
          f"sigma={sigma_dls_init:.4f}")
    print(f"  Fluor: A={A_fluor_init:.2f}, mu={mu_fluor_init:.4f}, "
          f"sigma={sigma_fluor_init:.4f}")

    # ── Optimize ────────────────────────────────────────────────────
    result = minimize(
        total_loss,
        p0,
        args=(dls_x, dls_y, fluor_x, fluor_y),
        method="Nelder-Mead",
        options={"maxiter": 50000, "xatol": 1e-10, "fatol": 1e-10},
    )

    if not result.success:
        print(f"\nWarning: optimizer did not fully converge: {result.message}")

    A_dls, mu_dls, sigma_dls, A_fluor, mu_fluor, sigma_fluor = result.x

    # ── Conversion factor ───────────────────────────────────────────
    conversion_lognormal = np.exp(mu_dls - mu_fluor)
    mode_dls = lognormal_mode(mu_dls, sigma_dls)
    mode_fluor = lognormal_mode(mu_fluor, sigma_fluor)
    conversion_mode_ratio = mode_dls / mode_fluor

    # Simple ratio-of-means for comparison (uses number distribution)
    conversion_simple, dls_num_mean, fluor_mean = ratio_of_means_conversion(
        num_x, num_y, sqrt_A_raw
    )

    # ── Report ──────────────────────────────────────────────────────
    print(f"\n{'=' * 60}")
    print("FIT RESULTS")
    print(f"{'=' * 60}")
    print(f"\nDLS log-normal fit:")
    print(f"  A     = {A_dls:.4f}")
    print(f"  mu    = {mu_dls:.6f}")
    print(f"  sigma = {sigma_dls:.6f}")
    print(f"  mode  = {mode_dls:.2f} nm")
    print(f"  median = {np.exp(mu_dls):.2f} nm")

    print(f"\nFluorescence log-normal fit:")
    print(f"  A     = {A_fluor:.4f}")
    print(f"  mu    = {mu_fluor:.6f}")
    print(f"  sigma = {sigma_fluor:.6f}")
    print(f"  mode  = {mode_fluor:.4f} sqrt(A) units")
    print(f"  median = {np.exp(mu_fluor):.4f} sqrt(A) units")

    print(f"\n{'=' * 60}")
    print("CONVERSION FACTORS")
    print(f"{'=' * 60}")
    print(f"\n  Log-normal (exp(mu_dls - mu_fluor)):  {conversion_lognormal:.6f} nm/sqrt(A)")
    print(f"  Mode ratio (mode_dls / mode_fluor):   {conversion_mode_ratio:.6f} nm/sqrt(A)")
    print(f"  Simple ratio of means (old method):   {conversion_simple:.6f} nm/sqrt(A)")

    print(f"\n  DLS number-weighted mean:       {dls_num_mean:.2f} nm")
    print(f"  Mean sqrt(A):                  {fluor_mean:.4f}")

    print(f"\n  Final loss: {result.fun:.6f}")

    print(f"\n{'=' * 60}")
    print("USE IN PIPELINE")
    print(f"{'=' * 60}")

    print(f"\n  Using log-normal conversion (recommended):")
    print(f"    python plot_curvature.py --conversion-factor {conversion_lognormal:.6f} ...")
    print(f"\n  Using mode ratio conversion:")
    print(f"    python plot_curvature.py --conversion-factor {conversion_mode_ratio:.6f} ...")
    print(f"\n  Using simple ratio-of-means (number distribution):")
    print(f"    python plot_curvature.py --dls-mean-diameter {dls_num_mean:.2f} ...")

    # ── Plots ───────────────────────────────────────────────────────
    if args.save_dir:
        os.makedirs(args.save_dir, exist_ok=True)

        # Plot 1: DLS fit
        fig, axes = plt.subplots(1, 3, figsize=(18, 5))

        ax = axes[0]
        ax.plot(dls_x, dls_y, "o", markersize=3, alpha=0.6, label="DLS data")
        d_fine = np.linspace(dls_x.min(), dls_x.max(), 500)
        ax.plot(d_fine, lognormal(d_fine, A_dls, mu_dls, sigma_dls),
                "r-", linewidth=2, label="Log-normal fit")
        ax.axvline(mode_dls, color="green", linestyle="--", linewidth=1,
                   label=f"Mode = {mode_dls:.1f} nm")
        ax.set_xlabel("Diameter (nm)")
        ax.set_ylabel("Intensity (%)")
        ax.set_title("DLS Intensity Distribution")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.2)

        # Plot 2: Fluorescence fit
        ax = axes[1]
        ax.bar(fluor_x, fluor_y, width=fluor_x[1] - fluor_x[0],
               alpha=0.5, label="sqrt(A) histogram")
        f_fine = np.linspace(fluor_x.min(), fluor_x.max(), 500)
        ax.plot(f_fine, lognormal(f_fine, A_fluor, mu_fluor, sigma_fluor),
                "r-", linewidth=2, label="Log-normal fit")
        ax.axvline(mode_fluor, color="green", linestyle="--", linewidth=1,
                   label=f"Mode = {mode_fluor:.2f}")
        ax.set_xlabel("sqrt(Lipid Amplitude)")
        ax.set_ylabel("Count")
        ax.set_title("Fluorescence sqrt(A) Distribution")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.2)

        # Plot 3: Overlay on shared nm axis
        ax = axes[2]
        dls_fit = lognormal(d_fine, A_dls, mu_dls, sigma_dls)
        dls_fit_norm = dls_fit / (np.max(dls_fit) + 1e-12)

        fluor_nm = f_fine * conversion_lognormal
        fluor_fit = lognormal(f_fine, A_fluor, mu_fluor, sigma_fluor)
        fluor_fit_norm = fluor_fit / (np.max(fluor_fit) + 1e-12)

        ax.plot(d_fine, dls_fit_norm, "b-", linewidth=2, label="DLS fit")
        ax.plot(fluor_nm, fluor_fit_norm, "r-", linewidth=2,
                label="Fluor fit (converted to nm)")
        ax.set_xlabel("Diameter (nm)")
        ax.set_ylabel("Normalized amplitude")
        ax.set_title("Overlay — Shared nm Axis")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.2)

        fig.tight_layout()
        save_path = os.path.join(args.save_dir, "dls_lognormal_calibration.png")
        fig.savefig(save_path, dpi=300)
        plt.close(fig)
        print(f"\nPlots saved to {save_path}")


if __name__ == "__main__":
    main()
