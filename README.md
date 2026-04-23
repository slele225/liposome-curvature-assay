# Single Liposome Curvature Assay Pipeline

A Python pipeline for analyzing protein curvature sensing on liposomes from
fluorescence microscopy images. Pairs with the
[CMEanalysis MATLAB detection code](https://github.com/DanuserLab/cmeAnalysis)
for sub-diffraction-limit spot detection.

> **New to this repo?** See the full step-by-step protocol document
> (`Single_Liposome_Curvature_Assay_Protocol.docx`) for beginner-friendly
> instructions with troubleshooting.

## For Non-Coders

If you're a wet-lab user running this pipeline for the first time and
you've never used Python or the command line, start with
**[PROTOCOL.md](PROTOCOL.md)** — a step-by-step walkthrough with
Windows/Mac commands side-by-side, placeholder paths, expected outputs,
and common-error troubleshooting. This README is the terse developer
reference; PROTOCOL.md is the friendly one.

## Getting Started

### 1. Install prerequisites

**Python 3.9+** — check with `python --version` (or `python3 --version` on
some Macs). If not installed, download from https://www.python.org/downloads/.
On Windows, check "Add python.exe to PATH" during install.

**uv** (package manager) — install by running:

macOS/Linux:
```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

Windows PowerShell:
```powershell
powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"
```

Close and reopen your terminal after installing uv.

### 2. Clone the repo

```bash
git clone https://github.com/slele225/liposome-curvature-assay.git
cd liposome-curvature-assay
```

If you don't have git, go to
https://github.com/slele225/liposome-curvature-assay, click the green
**Code** button, then **Download ZIP**. Unzip it and `cd` into the folder.

### 3. Set up the Python environment

```bash
uv venv
uv pip install -r requirements.txt
```

### 4. Activate the virtual environment

macOS/Linux:
```bash
source .venv/bin/activate
```

Windows PowerShell:
```powershell
.venv\Scripts\Activate
```

Your terminal prompt should now start with `(.venv)`. **You need to do this
every time you open a new terminal.**

### 5. Create your data folders and verify

```bash
mkdir data figures
python plot_curvature.py --help
```

### 6. Add your data

Copy your raw TIFFs and DLS `.xlsx` files into the `data/` folder:

```
liposome-curvature-assay/
└── data/
    ├── my_dls.xlsx
    └── my_experiment/
        ├── image001.tif
        ├── image002.tif
        └── ...
```

The `data/` and `figures/` folders are gitignored — your data stays on your
machine and will never be pushed to GitHub.

### Updating to the latest version

If the code gets updated, pull the changes:

```bash
cd liposome-curvature-assay
git pull
```

Your `data/` and `figures/` folders are unaffected.

## Overview

The assay works like this:

1. **Image liposomes** — Two-channel fluorescence microscopy: one channel is
   lipid dye, the other is bound protein.
2. **Detect spots** — MATLAB code (external) fits Gaussians to each
   sub-diffraction punctum and reports amplitudes.
3. **Calibrate sizes** — Dynamic light scattering (DLS) gives the true size
   distribution of the liposome stock. By comparing the mean DLS diameter to
   the mean sqrt(lipid amplitude), you get a conversion factor from
   fluorescence amplitude to physical diameter.
4. **Compute curvature sorting** — For each punctum, convert lipid amplitude
   to liposome diameter, then compute protein surface density =
   protein_A / (πD²).

## Pipeline

```
Raw TIFFs                ┌─────────────────────┐
  (microscope)     ───►  │ 1. prepare_input.py  │  Split channels, reorder,
                         └────────┬────────────┘  crop, organize by voltage
                                  │
                                  ▼
                         ┌─────────────────────┐
                         │ 2. MATLAB detection  │  External: CMEanalysis
                         │    (not in this repo)│  (see protocol doc)
                         └────────┬────────────┘
                                  │
                                  ▼
                         ┌─────────────────────┐
                         │ 3. analyze_matlab.py │  Read detection_v2.mat,
                         └────────┬────────────┘  filter puncta, export A vals
                                  │
                         ┌────────┴────────────────────┐
                         │                             │
                         ▼                             ▼
              ┌──────────────────────┐   ┌───────────────────────┐
              │ 4. plot_curvature.py │   │ 4b. plot_histograms.py│
              │    (needs protein)   │   │    (lipid-only OK)    │
              └──────────────────────┘   └───────────────────────┘

DLS .xlsx ──► dls_calibration.py ──► mean diameter ──► feeds into step 4
```

## Folder Structure

```
liposome-curvature-assay/
├── data/                       ← put your data here (gitignored)
│   ├── dls_data.xlsx
│   └── experiment_name/
│       └── 488nm_.../          ← created by prepare_input.py
│           ├── cell1/ch1/ ...
│           └── filtered_puncta_A_values.txt
├── figures/                    ← plots saved here (gitignored)
├── prepare_input.py
├── analyze_matlab.py
├── dls_calibration.py
├── plot_curvature.py
├── plot_histograms.py
├── plot_overlay.py
├── plot_dls.py
├── plot_dls_comparison.py
├── requirements.txt
├── .gitignore
└── README.md
```

## Usage

Every script uses `argparse` — run `python <script>.py --help` for full options.

### Step 1: Prepare images

Split multi-frame TIFFs into per-channel folders organized by detector voltage.

```bash
python prepare_input.py \
    --input  data/march_3_experiment \
    --output data/march_3_experiment_matlab \
    --frames 2,0 \
    --crop 1
```

| Argument    | Meaning |
|-------------|---------|
| `--input`   | Folder containing raw .tif files |
| `--output`  | Where to save the split channels (inside `data/`) |
| `--frames`  | Frame index order. `2,0` = frame 2 → ch1 (lipid), frame 0 → ch2 (protein) |
| `--crop`    | Center crop divisor. `1` = no crop, `2` = center quarter |

### Step 2: MATLAB detection (external)

See the protocol document for detailed CMEanalysis instructions. The detection
produces `detection_v2.mat` files inside the master channel's `Detection/`
subdirectory. Only the lipid/master channel (ch1) gets this folder.

### Step 3: Analyze MATLAB output

Filter puncta by intensity threshold and export amplitudes.

```bash
python analyze_matlab.py \
    --input  data/march_3_experiment_matlab/488nm_530V_561nm_500V \
    --channels ch1,ch2 \
    --lipid-channel ch1 \
    --k-std 2.0 \
    --output-name filtered_puncta_A_values.txt
```

| Argument          | Meaning |
|-------------------|---------|
| `--input`         | Condition folder (voltage group from Step 1) |
| `--channels`      | Channel names matching folder names. Use `ch1` for lipid-only |
| `--lipid-channel` | Which channel is lipid (used for thresholding) |
| `--k-std`         | Threshold: keep if A > mean(c) + k·std(c). Default `2.0` |
| `--output-name`   | Output filename, saved inside `--input` folder |

For **lipid-only** experiments (no protein channel):
```bash
python analyze_matlab.py \
    --input  data/.../488nm_530V_561nm_500V \
    --channels ch1 \
    --lipid-channel ch1
```

### Step 4a: Plot curvature sorting

Requires both lipid and protein channels. Accepts either a conversion factor
(from the log-normal DLS calibration) or a mean diameter (simple method).

```bash
# Using conversion factor (recommended, from dls_calibration.py)
python plot_curvature.py \
    --input data/.../filtered_puncta_A_values.txt \
    --conversion-factor 1.234567 \
    --save-dir figures/

# Using simple ratio-of-means method
python plot_curvature.py \
    --input data/.../filtered_puncta_A_values.txt \
    --dls-mean-diameter 80.12 \
    --save-dir figures/
```

| Argument              | Meaning |
|-----------------------|---------|
| `--input`             | Filtered puncta file from Step 3 (accepts multiple files) |
| `--conversion-factor` | Maps sqrt(A) to diameter in nm. From `dls_calibration.py` log-normal fit |
| `--dls-mean-diameter` | Alternative: mean diameter in nm (uses ratio-of-means internally) |
| `--lipid-col`         | Column for lipid amplitude (default: `A_ch1`) |
| `--protein-col`       | Column for protein amplitude (default: `A_ch2`) |
| `--bin-width`         | Diameter bin width in nm for averaged curve (default: `0.5`) |
| `--diameter-cutoff`   | Exclude puncta with diameter above this value in nm. Useful for filtering sparse large-diameter tail (optional) |
| `--y-pad`             | Y-axis padding factor around bin means. Default 0.3 (30%). Lower values (e.g., 0.1) zoom in tighter on the trend (optional) |
| `--save-dir`          | Output directory for figures |

Provide exactly one of `--conversion-factor` or `--dls-mean-diameter`.

### Step 4b: Plot histograms

Works with or without a protein channel. Useful for sanity-checking
distributions and for lipid-only experiments.

```bash
# Two-channel (lipid + protein)
python plot_histograms.py \
    --input data/.../filtered_puncta_A_values.txt \
    --lipid-col A_ch1 \
    --protein-col A_ch2 \
    --dls-mean-diameter 80.11 \
    --save-dir figures/

# Lipid-only (omit --protein-col)
python plot_histograms.py \
    --input data/.../filtered_puncta_A_values.txt \
    --lipid-col A_ch1 \
    --dls-mean-diameter 80.11 \
    --save-dir figures/
```

| Argument              | Meaning |
|-----------------------|---------|
| `--input`             | Filtered puncta file from Step 3 |
| `--lipid-col`         | Column for lipid amplitude (default: `A_ch1`) |
| `--protein-col`       | Column for protein amplitude. Omit for lipid-only |
| `--conversion-factor` | Conversion factor from `dls_calibration.py`. Omit to skip diameter plot |
| `--bins`              | Number of histogram bins (default: `80`) |
| `--transform`         | `raw`, `sqrt`, or `log_sqrt` (default: `raw`) |
| `--save-dir`          | Output directory for figures |

### Optional: DLS calibration

Computes a single scalar conversion factor that maps `sqrt(lipid amplitude)`
to physical diameter in nm. The script finds this factor by **overlaying**
the fluorescence sqrt(A) distribution onto the DLS number-weighted size
distribution — rebinning the fluorescence data onto the DLS bin grid and
minimizing the chi-squared difference. This formalizes the standard SLiC
calibration approach (Kunding 2008, Hatzakis 2009, Bhatia 2009, Zeno 2018,
Johnson 2025), which uses a single scalar to convert sqrt(intensity) to
diameter.

The ratio-of-means (DLS mean / mean(sqrt(A))) is also reported as a quick
sanity check. The two typically agree to within a few percent; the overlay
is more robust because it fits the full distribution shape rather than just
the first moment.

**DLS data preparation:** Export the size distribution from the Malvern
Zetasizer software and copy it into an Excel spreadsheet with the standard
Zetasizer format ("X Intensity", "X Volume", "X Number" section headers).
Only the number distribution is used.

```bash
# Standard usage
python dls_calibration.py \
    --dls-input data/dls_data.xlsx \
    --fluor-input data/.../filtered_puncta_A_values.txt \
    --save-dir figures/

# With bootstrap variance estimate
python dls_calibration.py \
    --dls-input data/dls_data.xlsx \
    --fluor-input data/.../filtered_puncta_A_values.txt \
    --bootstrap 200 \
    --save-dir figures/
```

| Argument          | Meaning |
|-------------------|---------|
| `--dls-input`     | Path to Zetasizer `.xlsx` export |
| `--fluor-input`   | Path to `filtered_puncta_A_values.txt` from Step 3 |
| `--bootstrap N`   | Bootstrap resamples for variance estimate (default: `0` = off) |
| `--bootstrap-k K` | Puncta per bootstrap iteration (default: same as total) |
| `--lipid-col`     | Column name for lipid amplitude (default: `A_ch1`) |
| `--save-dir`      | Save overlay plots (optional) |

**Using the output:** The script prints a conversion factor and an implied
mean diameter. Use `--conversion-factor` in the plotting scripts, or
equivalently `--dls-mean-diameter` with the implied mean.

**Bootstrap:** If you run with `--bootstrap`, the script resamples puncta
with replacement N times and compares the CV of the overlay vs ratio-of-means
conversion factors. A lower CV means the method is more robust to which
puncta were detected.

### Optional: Normalized overlay across experiments

Compare curvature sorting across conditions by normalizing each curve so
the largest-diameter bin = 1 (fold-enrichment at high curvature). Each input
file gets its own conversion factor (from `dls_calibration.py`).

```bash
python plot_overlay.py \
    --input data/cond1/filtered.txt:1.234 data/cond2/filtered.txt:1.567 \
    --labels "WT protein" "Mutant K58A" \
    --save-dir figures/
```

| Argument        | Meaning |
|-----------------|---------|
| `--input`       | Files with conversion factors, as `file.txt:factor` pairs |
| `--labels`      | Custom legend labels (default: parent folder name) |
| `--bin-width`   | Diameter bin width in nm (default: `0.5`) |
| `--normalize-to` | How to normalize curves: `rightmost` (default), `leftmost`, `minimum`, or `none`. Rightmost matches Bhatia 2009 / Zeno 2018 convention |
| `--diameter-cutoff` | Exclude puncta with diameter above this value in nm. Useful for filtering sparse large-diameter tail (optional) |
| `--y-pad`       | Y-axis padding factor around the plotted bin means range. Default 0.3 (30%). Lower values (e.g., 0.1) zoom in tighter on the trend (optional) |
| `--output-name` | Output filename (default: `normalized_curvature_overlay.png`) |
| `--save-dir`    | Output directory for figure |

The conversion factor for each file maps sqrt(A) to diameter in nm. Get it
from `dls_calibration.py`.

### Optional: DLS distribution plots

Quick visualization of the DLS distribution, with optional log x-axis.
Choose between number or intensity distribution.

```bash
python plot_dls.py data/dls_data.xlsx                          # number, raw diameter
python plot_dls.py data/dls_data.xlsx --log                     # number, log(diameter)
python plot_dls.py data/dls_data.xlsx --distribution intensity  # intensity distribution
python plot_dls.py data/dls_data.xlsx --zoom-pct 100            # show full range
```

| Argument           | Meaning |
|--------------------|---------|
| `input`            | Path to Zetasizer `.xlsx` file |
| `--distribution`   | `number` (default) or `intensity` |
| `--log`            | Plot log(diameter) instead of raw |
| `--zoom-pct`       | Percent of data to show (default: `95`) |

### Optional: DLS vs fluorescence comparison

Side-by-side comparison of DLS distribution and fluorescence sqrt(A)
distributions for one or more channels. One panel per channel, all
independently zoomed. Use this to check whether your detected liposome
population matches the DLS measurement.

```bash
python plot_dls_comparison.py \
    --dls-input data/dls_data.xlsx \
    --fluor-input data/.../filtered_puncta_A_values.txt \
    --channels 0 Lipid 1 EGFP \
    --zoom-pct 95 \
    --bins 200 \
    --save-dir figures/
```

| Argument              | Meaning |
|-----------------------|---------|
| `--dls-input`         | DLS `.xlsx` file |
| `--fluor-input`       | Filtered puncta file |
| `--channels`          | Index/label pairs: `0 Lipid 1 EGFP`. Index 0 = A_ch1, 1 = A_ch2, etc. |
| `--dls-distribution`  | `number` (default) or `intensity` for the DLS panel |
| `--bins`              | Bins for sqrt(A) histograms (default: `100`) |
| `--zoom-pct`          | Percent of data to show (default: `100` = no zoom) |
| `--save-dir`          | Output directory |

## File Descriptions

| File                      | Purpose |
|---------------------------|---------|
| `prepare_input.py`        | Split and reorder TIFF channels, organize for MATLAB |
| `analyze_matlab.py`       | Read MATLAB detection `.mat` files, filter puncta, export TSV |
| `dls_calibration.py`      | DLS-fluorescence distribution overlay to compute conversion factor |
| `plot_curvature.py`       | Convert amplitudes to diameters, plot protein density vs diameter |
| `plot_histograms.py`      | Plot amplitude histograms and estimated diameter distributions |
| `plot_overlay.py`         | Overlay normalized curvature-sorting curves across experiments |
| `plot_dls.py`             | Plot DLS distribution (number or intensity, raw or log) |
| `plot_dls_comparison.py`  | Side-by-side DLS vs fluorescence sqrt(A) per channel |
| `requirements.txt`        | Python dependencies: numpy, matplotlib, tifffile, h5py, openpyxl, pandas, scipy |

## Notes

- **MATLAB detection** is maintained separately
  ([DanuserLab/cmeAnalysis](https://github.com/DanuserLab/cmeAnalysis)).
  This pipeline reads its output: `detection_v2.mat` files containing a
  `frameInfo` struct with fields `A` (amplitude), `c` (background), and
  `hval_Ar` (hypothesis test).
- **Only the master/lipid channel** (ch1) gets a `Detection/` subfolder from
  CMEanalysis. The protein channel (ch2) will just contain the TIFF.
- **DLS data preparation:** Export the size distribution from the Malvern
  Zetasizer software. The spreadsheet must include the number distribution.
  Copy it into an Excel file with the standard Zetasizer section headers
  ("X Intensity", "X Volume", "X Number").
- **DLS calibration** uses a distribution overlay approach: the script
  rebins the fluorescence sqrt(A) data onto the DLS bin grid and finds the
  single scalar k that minimizes the chi-squared difference between the two
  distributions. This is the standard SLiC calibration procedure used from
  Kunding 2008 through Johnson/Zeno 2025. The ratio-of-means is also
  reported as a sanity check. See `DLS_Calibration_Notes.md` for background.
- **Lipid-only experiments** are fully supported. Run `analyze_matlab.py`
  with `--channels ch1 --lipid-channel ch1`, then use `plot_histograms.py`
  (skip `plot_curvature.py` since it requires protein data).
- The `data/` and `figures/` folders are gitignored. Your microscopy data
  stays local.
