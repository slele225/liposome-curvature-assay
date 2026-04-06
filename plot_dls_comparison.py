"""
Side-by-side comparison of DLS distribution and fluorescence sqrt(A)
distributions for one or more channels.

Layout: one row with DLS on the left, then one panel per channel.
All panels zoomed to the central --zoom-pct of their data.

Example:
  python plot_dls_comparison.py \
      --dls-input data/dls.xlsx \
      --fluor-input data/.../filtered_puncta_A_values.txt \
      --channels 0 Lipid 1 EGFP \
      --zoom-pct 95 \
      --bins 200 \
      --save-dir figures/
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ── DLS loading ─────────────────────────────────────────────────────────

def load_dls_section(xlsx_path, section_name):
    """Load a DLS distribution section (number or intensity)."""
    df = pd.read_excel(xlsx_path, header=None)
    target = f"x {section_name}".lower()

    start_row = None
    for i in range(df.shape[0]):
        if str(df.iloc[i, 0]).strip().lower() == target:
            start_row = i
            break
    if start_row is None:
        raise ValueError(f"Could not find 'X {section_name}' section.")

    end_row = df.shape[0]
    for i in range(start_row + 1, df.shape[0]):
        val = str(df.iloc[i, 0]).strip()
        if val.startswith("X ") and val.lower() != target:
            end_row = i
            break

    data = df.iloc[start_row + 1 : end_row].apply(
        pd.to_numeric, errors="coerce"
    ).dropna(how="all")
    diameters = data.iloc[:, 0].values.astype(float)
    records = data.iloc[:, 1:].values.astype(float)
    avg = np.nanmean(records, axis=1)
    valid = np.isfinite(diameters) & np.isfinite(avg) & (avg > 0)
    return diameters[valid], avg[valid]


# ── Fluorescence loading ────────────────────────────────────────────────

def load_puncta_columns(puncta_path, col_names):
    """Load specified columns from filtered puncta file."""
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

    result = {}
    for col in col_names:
        if col not in headers:
            raise ValueError(f"Column '{col}' not found. Available: {headers}")
        idx = headers.index(col)
        vals = np.array([float(r[idx]) for r in rows])
        result[col] = np.sqrt(np.clip(vals, 0, None))

    return result


# ── Channel parsing ─────────────────────────────────────────────────────

def parse_channels(channel_args):
    """
    Parse channel arguments as index/label pairs.

    Input: ['0', 'Lipid', '1', 'EGFP']
    Output: [('A_ch1', 'Lipid'), ('A_ch2', 'EGFP')]
    """
    if len(channel_args) % 2 != 0:
        raise ValueError(
            "Channels must be index/label pairs: 0 Lipid 1 EGFP"
        )
    pairs = []
    for i in range(0, len(channel_args), 2):
        ch_idx = int(channel_args[i])
        label = channel_args[i + 1]
        col_name = f"A_ch{ch_idx + 1}"
        pairs.append((col_name, label))
    return pairs


# ── Zoom helper ─────────────────────────────────────────────────────────

def zoom_limits(data, weights, pct):
    """Return (lo, hi) x-limits showing central pct% of the distribution."""
    if pct >= 100:
        pad = (data[-1] - data[0]) * 0.05
        return data[0] - pad, data[-1] + pad

    tail = (100 - pct) / 2
    cumsum = np.cumsum(weights) / np.sum(weights)
    lo_idx = max(0, np.searchsorted(cumsum, tail / 100) - 1)
    hi_idx = min(len(data) - 1, np.searchsorted(cumsum, 1 - tail / 100) + 1)
    pad = (data[hi_idx] - data[lo_idx]) * 0.1
    return data[lo_idx] - pad, data[hi_idx] + pad


def zoom_limits_raw(values, pct):
    """Return (lo, hi) x-limits from percentiles of raw values."""
    if pct >= 100:
        pad = (np.max(values) - np.min(values)) * 0.05
        return np.min(values) - pad, np.max(values) + pad

    tail = (100 - pct) / 2
    lo = np.percentile(values, tail)
    hi = np.percentile(values, 100 - tail)
    pad = (hi - lo) * 0.1
    return lo - pad, hi + pad


# ── CLI ─────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Side-by-side DLS vs fluorescence distribution comparison.",
        epilog=(
            "Example:\n"
            "  python plot_dls_comparison.py \\\n"
            "      --dls-input data/dls.xlsx \\\n"
            "      --fluor-input data/.../filtered.txt \\\n"
            "      --channels 0 Lipid 1 EGFP \\\n"
            "      --zoom-pct 95 --bins 200 --save-dir figures/"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--dls-input", required=True, help="DLS .xlsx file")
    parser.add_argument("--fluor-input", required=True,
                        help="filtered_puncta_A_values.txt")
    parser.add_argument(
        "--channels", required=True, nargs="+",
        help="Channel index/label pairs: 0 Lipid 1 EGFP. "
             "Index 0 = A_ch1, index 1 = A_ch2, etc.",
    )
    parser.add_argument(
        "--dls-distribution", default="number",
        choices=["number", "intensity"],
        help="Which DLS distribution to plot (default: number)",
    )
    parser.add_argument("--bins", type=int, default=100,
                        help="Bins for sqrt(A) histograms (default: 100)")
    parser.add_argument("--zoom-pct", type=float, default=100.0,
                        help="Percent of data to show (default: 100 = no zoom)")
    parser.add_argument("--save-dir", default="figures/",
                        help="Output directory")
    args = parser.parse_args()

    # Parse channels
    try:
        channels = parse_channels(args.channels)
    except ValueError as e:
        print(f"Error: {e}")
        sys.exit(1)

    # Load DLS
    dls_d, dls_w = load_dls_section(args.dls_input, args.dls_distribution)

    # Load fluorescence columns
    col_names = [col for col, _ in channels]
    fluor_data = load_puncta_columns(args.fluor_input, col_names)

    n_panels = 1 + len(channels)  # DLS + each channel

    fig, axes = plt.subplots(1, n_panels, figsize=(6 * n_panels, 5))
    if n_panels == 1:
        axes = [axes]

    # ── DLS panel ───────────────────────────────────────────────────
    ax = axes[0]
    dls_w_norm = dls_w / (np.sum(dls_w) * np.mean(np.diff(dls_d)))
    widths = np.diff(np.append(dls_d, dls_d[-1] * 1.05))
    ax.bar(dls_d, dls_w_norm, width=widths, alpha=0.7,
           edgecolor="black", linewidth=0.3, align="edge", color="steelblue")

    dls_mean = np.sum(dls_d * dls_w) / np.sum(dls_w)
    ax.axvline(dls_mean, color="red", linestyle="--", linewidth=1.5,
               label=f"Mean = {dls_mean:.1f} nm")

    lo, hi = zoom_limits(dls_d, dls_w, args.zoom_pct)
    ax.set_xlim(lo, hi)
    ax.set_xlabel("Diameter (nm)")
    ax.set_ylabel("Density")
    ax.set_title(f"DLS {args.dls_distribution.capitalize()} Distribution")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.2)

    # ── Channel panels ──────────────────────────────────────────────
    for i, (col_name, label) in enumerate(channels):
        ax = axes[i + 1]
        sqrt_vals = fluor_data[col_name]
        sqrt_vals = sqrt_vals[sqrt_vals > 0]

        counts, edges = np.histogram(sqrt_vals, bins=args.bins, density=True)
        centers = 0.5 * (edges[:-1] + edges[1:])
        bin_widths = np.diff(edges)

        ax.bar(centers, counts, width=bin_widths, alpha=0.7,
               edgecolor="black", linewidth=0.3, color="coral")

        ch_mean = np.mean(sqrt_vals)
        ax.axvline(ch_mean, color="red", linestyle="--", linewidth=1.5,
                   label=f"Mean = {ch_mean:.2f}")

        lo, hi = zoom_limits_raw(sqrt_vals, args.zoom_pct)
        ax.set_xlim(lo, hi)
        ax.set_xlabel(f"sqrt(A) — {label}")
        ax.set_ylabel("Density")
        ax.set_title(f"Ch{int(col_name.split('ch')[1])} — {label}\n(n = {len(sqrt_vals)})")
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.2)

    fig.suptitle("DLS vs Fluorescence — Shape Comparison", fontsize=14)
    fig.tight_layout()

    os.makedirs(args.save_dir, exist_ok=True)
    out = os.path.join(args.save_dir, "dls_vs_fluorescence_comparison.png")
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved -> {out}")


if __name__ == "__main__":
    main()
