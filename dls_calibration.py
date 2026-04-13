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

def make_z_lognormal(dls_x, z_avg, pdi):
    """
    Generate a lognormal distribution from Z-average and PDI.

    The Z-average is the intensity-weighted harmonic mean from the
    cumulant fit. We use mu = ln(Z_avg) and sigma = sqrt(PDI) as the
    lognormal parameters, then evaluate on the DLS diameter bins.
    """
    mu = np.log(z_avg)
    sigma = np.sqrt(pdi) if pdi > 0 else 0.3
    with np.errstate(divide="ignore", invalid="ignore"):
        logx = np.log(np.clip(dls_x, 1e-12, None))
        vals = np.exp(-0.5 * ((logx - mu) / sigma) ** 2)
    # Normalize to same area as a percentage distribution
    vals = vals / (np.sum(vals) + 1e-12) * 100
    return vals


def total_loss(params, dls_x, num_y, z_lognormal_y, fluor_x, fluor_y):
    """
    Combined loss from all three terms, each normalized to relative error.

    params = [alpha, A_dls, mu_dls, sigma_dls, A_fluor, mu_fluor, sigma_fluor]

    alpha: weight for linear combination
        pseudo_dls = alpha * number_distribution + (1 - alpha) * z_lognormal
    """
    alpha, A_dls, mu_dls, sigma_dls, A_fluor, mu_fluor, sigma_fluor = params

    # Enforce constraints
    if sigma_dls <= 0 or sigma_fluor <= 0:
        return 1e12
    if A_dls <= 0 or A_fluor <= 0:
        return 1e12
    if alpha < 0 or alpha > 1:
        return 1e12

    # Construct pseudo-DLS as linear combination
    pseudo_dls = alpha * num_y + (1 - alpha) * z_lognormal_y

    # L1: Lognormal fit to pseudo-DLS (normalized by data magnitude)
    pred_dls = lognormal(dls_x, A_dls, mu_dls, sigma_dls)
    L1 = np.sum((pred_dls - pseudo_dls) ** 2) / (np.sum(pseudo_dls ** 2) + 1e-12)

    # L2: Fluorescence lognormal fit to empirical sqrt(A) (normalized)
    pred_fluor = lognormal(fluor_x, A_fluor, mu_fluor, sigma_fluor)
    L2 = np.sum((pred_fluor - fluor_y) ** 2) / (np.sum(fluor_y ** 2) + 1e-12)

    # L3: Consistency between the two lognormals on a shared nm axis
    conversion = np.exp(mu_dls - mu_fluor)

    dls_on_shared = lognormal(dls_x, A_dls, mu_dls, sigma_dls)
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
        "--z-avg", type=float, default=None,
        help="Z-average diameter (nm) from Zetasizer cumulants analysis. "
             "Required for log-normal fitting. Omit to use ratio-of-means only.",
    )
    parser.add_argument(
        "--pdi", type=float, default=None,
        help="PDI from Zetasizer cumulants analysis. "
             "Required for log-normal fitting. Omit to use ratio-of-means only.",
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
    parser.add_argument(
        "--bootstrap", type=int, default=0,
        help="Number of bootstrap resamples to compare variance of "
             "ratio-of-means vs lognormal conversion (default: 0 = off). "
             "Try 200-500. Requires --z-avg and --pdi.",
    )
    parser.add_argument(
        "--bootstrap-k", type=int, default=None,
        help="Number of puncta to sample per bootstrap iteration "
             "(default: same as total puncta count).",
    )
    args = parser.parse_args()

    # ── Load data ───────────────────────────────────────────────────
    do_lognormal = args.z_avg is not None and args.pdi is not None

    print("=" * 60)
    if do_lognormal:
        print("DLS CALIBRATION — LOG-NORMAL FITTING")
    else:
        print("DLS CALIBRATION — RATIO OF MEANS")
    print("=" * 60)

    num_x, num_y = load_dls_number(args.dls_input)
    fluor_x, fluor_y, sqrt_A_raw = load_sqrt_A(
        args.fluor_input, args.lipid_col, args.n_bins
    )

    print(f"DLS file:       {args.dls_input}")
    print(f"Fluor file:     {args.fluor_input}")
    print(f"DLS number bins:    {len(num_x)}")
    print(f"Fluor bins:     {len(fluor_x)}")
    print(f"Fluor points:   {len(sqrt_A_raw)}")

    # ── Ratio-of-means (always computed) ────────────────────────────
    conversion_simple, dls_num_mean, fluor_mean = ratio_of_means_conversion(
        num_x, num_y, sqrt_A_raw
    )

    print(f"\n{'=' * 60}")
    print("RATIO OF MEANS (number distribution)")
    print(f"{'=' * 60}")
    print(f"  DLS number-weighted mean:  {dls_num_mean:.2f} nm")
    print(f"  Mean sqrt(A):              {fluor_mean:.4f}")
    print(f"  Conversion factor:         {conversion_simple:.6f} nm/sqrt(A)")
    print(f"\n  python plot_curvature.py --dls-mean-diameter {dls_num_mean:.2f} ...")

    if not do_lognormal:
        if args.z_avg is None or args.pdi is None:
            print(f"\nSkipping log-normal fitting (--z-avg and --pdi not provided).")
            print("Add both to enable log-normal calibration.")
        return

    # ── Log-normal fitting (pseudo-DLS approach) ──────────────────
    dls_x, dls_y = num_x, num_y

    # Generate lognormal from Z-avg and PDI (cumulant representation)
    z_lognormal_y = make_z_lognormal(dls_x, args.z_avg, args.pdi)

    print(f"\n  Z-average:        {args.z_avg} nm")
    print(f"  PDI:              {args.pdi}")
    print(f"  DLS number bins:  {len(dls_x)}")
    print(f"  Pseudo-DLS = alpha * number_dist + (1-alpha) * Z-avg_lognormal")

    # ── Initial guesses (7 params: alpha + 3 DLS + 3 fluor) ────
    alpha_init = 0.5
    mu_dls_init = np.log(args.z_avg)
    sigma_dls_init = np.sqrt(args.pdi) if args.pdi > 0 else 0.3
    A_dls_init = float(np.max(dls_y))

    mu_fluor_init = float(np.log(np.mean(sqrt_A_raw)))
    sigma_fluor_init = float(np.std(np.log(sqrt_A_raw)))
    A_fluor_init = float(np.max(fluor_y))

    p0 = [
        alpha_init,
        A_dls_init, mu_dls_init, sigma_dls_init,
        A_fluor_init, mu_fluor_init, sigma_fluor_init,
    ]

    print(f"\nInitial guesses:")
    print(f"  alpha: {alpha_init}")
    print(f"  DLS:   A={A_dls_init:.2f}, mu={mu_dls_init:.4f}, "
          f"sigma={sigma_dls_init:.4f}")
    print(f"  Fluor: A={A_fluor_init:.2f}, mu={mu_fluor_init:.4f}, "
          f"sigma={sigma_fluor_init:.4f}")

    # ── Optimize ────────────────────────────────────────────────────
    result = minimize(
        total_loss,
        p0,
        args=(dls_x, dls_y, z_lognormal_y, fluor_x, fluor_y),
        method="Nelder-Mead",
        options={"maxiter": 100000, "xatol": 1e-10, "fatol": 1e-10},
    )

    if not result.success:
        print(f"\nWarning: optimizer did not fully converge: {result.message}")

    alpha, A_dls, mu_dls, sigma_dls, A_fluor, mu_fluor, sigma_fluor = result.x

    # Reconstruct the pseudo-DLS that was actually fit
    pseudo_dls_final = alpha * dls_y + (1 - alpha) * z_lognormal_y

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
    print("LOG-NORMAL FIT RESULTS")
    print(f"{'=' * 60}")

    print(f"\nPseudo-DLS mixing:")
    print(f"  alpha = {alpha:.4f}")
    print(f"  pseudo-DLS = {alpha:.2f} * number + {1-alpha:.2f} * Z-avg lognormal")

    print(f"\nDLS log-normal fit (to pseudo-DLS):")
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

    print(f"\n  Final loss: {result.fun:.6f}")

    print(f"\n{'=' * 60}")
    print("ALL CONVERSION FACTORS")
    print(f"{'=' * 60}")
    print(f"\n  Log-normal (exp(mu_dls - mu_fluor)):  {conversion_lognormal:.6f} nm/sqrt(A)")
    print(f"  Mode ratio (mode_dls / mode_fluor):   {conversion_mode_ratio:.6f} nm/sqrt(A)")
    print(f"  Ratio of means (number distribution): {conversion_simple:.6f} nm/sqrt(A)")

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

        # Plot 1: Pseudo-DLS and fit
        fig, axes = plt.subplots(1, 3, figsize=(18, 5))

        ax = axes[0]
        ax.plot(dls_x, pseudo_dls_final, "o", markersize=3, alpha=0.6,
                label="Pseudo-DLS")
        ax.plot(dls_x, dls_y, "s", markersize=2, alpha=0.3, color="gray",
                label="Number dist")
        ax.plot(dls_x, z_lognormal_y, "--", color="gray", alpha=0.4,
                linewidth=1, label="Z-avg lognormal")
        d_fine = np.linspace(dls_x.min(), dls_x.max(), 500)
        ax.plot(d_fine, lognormal(d_fine, A_dls, mu_dls, sigma_dls),
                "r-", linewidth=2, label="Log-normal fit")
        ax.axvline(mode_dls, color="green", linestyle="--", linewidth=1,
                   label=f"Mode = {mode_dls:.1f} nm")
        ax.set_xlabel("Diameter (nm)")
        ax.set_ylabel("Number (%)")
        ax.set_title(f"Pseudo-DLS (α={alpha:.2f})")
        ax.legend(fontsize=7)
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

    # ── Bootstrap variance comparison ──────────────────────────────
    if args.bootstrap > 0 and do_lognormal:
        n_boot = args.bootstrap
        k = args.bootstrap_k or len(sqrt_A_raw)
        rng = np.random.default_rng(42)

        print(f"\n{'=' * 60}")
        print(f"BOOTSTRAP VARIANCE COMPARISON (n={n_boot}, k={k})")
        print(f"{'=' * 60}")

        rom_factors = []
        ln_factors = []
        mode_factors = []

        for b in range(n_boot):
            # Resample fluorescence data
            sample = rng.choice(sqrt_A_raw, size=k, replace=True)

            # Ratio-of-means
            rom = dls_num_mean / np.mean(sample)
            rom_factors.append(rom)

            # Lognormal: re-histogram and re-optimize
            counts_b, edges_b = np.histogram(sample, bins=args.n_bins)
            centers_b = 0.5 * (edges_b[:-1] + edges_b[1:])
            counts_b = counts_b.astype(float)

            p0_b = [
                alpha,  # use previous alpha as init
                A_dls, mu_dls, sigma_dls,  # DLS side doesn't change
                float(np.max(counts_b)),
                float(np.log(np.mean(sample))),
                float(np.std(np.log(sample))),
            ]

            try:
                res_b = minimize(
                    total_loss,
                    p0_b,
                    args=(dls_x, dls_y, z_lognormal_y, centers_b, counts_b),
                    method="Nelder-Mead",
                    options={"maxiter": 50000, "xatol": 1e-8, "fatol": 1e-8},
                )
                _, _, mu_d_b, _, _, mu_f_b, sigma_f_b = res_b.x
                ln_factors.append(np.exp(mu_d_b - mu_f_b))
                mode_factors.append(
                    lognormal_mode(mu_d_b, res_b.x[3]) /
                    lognormal_mode(mu_f_b, sigma_f_b)
                )
            except Exception:
                pass  # skip failed optimizations

            if (b + 1) % 50 == 0 or b == 0:
                print(f"  Bootstrap {b+1}/{n_boot}...")

        rom_factors = np.array(rom_factors)
        ln_factors = np.array(ln_factors)
        mode_factors = np.array(mode_factors)

        print(f"\n  Ratio-of-means:  mean={np.mean(rom_factors):.4f}, "
              f"std={np.std(rom_factors):.4f}, "
              f"CV={np.std(rom_factors)/np.mean(rom_factors)*100:.2f}%")
        print(f"  Log-normal:      mean={np.mean(ln_factors):.4f}, "
              f"std={np.std(ln_factors):.4f}, "
              f"CV={np.std(ln_factors)/np.mean(ln_factors)*100:.2f}%")
        print(f"  Mode ratio:      mean={np.mean(mode_factors):.4f}, "
              f"std={np.std(mode_factors):.4f}, "
              f"CV={np.std(mode_factors)/np.mean(mode_factors)*100:.2f}%")

        # Plot bootstrap distributions
        if args.save_dir:
            fig, ax = plt.subplots(figsize=(10, 5))
            bins_hist = 40

            ax.hist(rom_factors, bins=bins_hist, alpha=0.5, color="blue",
                    label=f"Ratio-of-means (CV={np.std(rom_factors)/np.mean(rom_factors)*100:.1f}%)",
                    density=True)
            ax.hist(ln_factors, bins=bins_hist, alpha=0.5, color="red",
                    label=f"Log-normal (CV={np.std(ln_factors)/np.mean(ln_factors)*100:.1f}%)",
                    density=True)
            ax.hist(mode_factors, bins=bins_hist, alpha=0.5, color="green",
                    label=f"Mode ratio (CV={np.std(mode_factors)/np.mean(mode_factors)*100:.1f}%)",
                    density=True)

            ax.axvline(np.mean(rom_factors), color="blue", linestyle="--", linewidth=1.5)
            ax.axvline(np.mean(ln_factors), color="red", linestyle="--", linewidth=1.5)
            ax.axvline(np.mean(mode_factors), color="green", linestyle="--", linewidth=1.5)

            ax.set_xlabel("Conversion factor (nm / sqrt(A))")
            ax.set_ylabel("Density")
            ax.set_title(f"Bootstrap Comparison (n={n_boot}, k={k})")
            ax.legend(fontsize=9)
            ax.grid(True, alpha=0.2)
            fig.tight_layout()

            boot_path = os.path.join(args.save_dir, "bootstrap_variance.png")
            fig.savefig(boot_path, dpi=300)
            plt.close(fig)
            print(f"\n  Bootstrap plot saved to {boot_path}")


if __name__ == "__main__":
    main()