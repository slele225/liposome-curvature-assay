"""Plot DLS distribution - number or intensity, raw or log diameter."""

import os
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


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


parser = argparse.ArgumentParser(description="Plot DLS distribution.")
parser.add_argument("input", help="Path to DLS .xlsx file")
parser.add_argument("--distribution", default="number",
                    choices=["number", "intensity"],
                    help="Which distribution to plot (default: number)")
parser.add_argument("--log", action="store_true",
                    help="Plot log(diameter) instead of raw diameter")
parser.add_argument("--zoom-pct", type=float, default=95.0,
                    help="Percent of data to show (default: 95)")
args = parser.parse_args()

d, w = load_dls_section(args.input, args.distribution)

fig, ax = plt.subplots(figsize=(8, 5))

if args.log:
    x = np.log(d)
    xlabel = "log(Diameter nm)"
    suffix = "log"
else:
    x = d
    xlabel = "Diameter (nm)"
    suffix = "raw"

widths = np.diff(np.append(x, x[-1] + (x[-1] - x[-2])))
ax.bar(x, w, width=widths, alpha=0.7, edgecolor="black", linewidth=0.3,
       align="edge")

mean_x = np.sum(x * w) / np.sum(w)
ax.axvline(mean_x, color="red", linestyle="--", linewidth=1,
           label=f"Weighted mean = {mean_x:.3f}")

ax.set_xlabel(xlabel)
ax.set_ylabel(f"{args.distribution.capitalize()} %")
ax.set_title(f"DLS {args.distribution.capitalize()} Distribution"
             f"{' — log scale' if args.log else ''}")
ax.legend()
ax.grid(True, alpha=0.2)

# Zoom to central zoom_pct%
if args.zoom_pct < 100:
    cumsum = np.cumsum(w) / np.sum(w)
    tail = (100 - args.zoom_pct) / 2
    lo_idx = max(0, np.searchsorted(cumsum, tail / 100) - 1)
    hi_idx = min(len(x) - 1, np.searchsorted(cumsum, 1 - tail / 100) + 1)
    pad = (x[hi_idx] - x[lo_idx]) * 0.15
    ax.set_xlim(x[lo_idx] - pad, x[hi_idx] + pad)

fig.tight_layout()

os.makedirs("figures", exist_ok=True)
out = os.path.join("figures", f"dls_{args.distribution}_{suffix}_diameter.png")
fig.savefig(out, dpi=300)
print(f"Saved -> {out}")
plt.close()
