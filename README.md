# Single Liposome Curvature Assay Pipeline

A Python pipeline for analyzing protein curvature sensing on liposomes from
fluorescence microscopy images. Pairs with the
[CMEanalysis MATLAB detection code](https://github.com/DanuserLab/cmeAnalysis)
for sub-diffraction-limit spot detection.

> **New to this repo?** See the full step-by-step protocol document
> (`Single_Liposome_Curvature_Assay_Protocol.docx`) for beginner-friendly
> instructions with troubleshooting.

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
   fluorescence amplitude to physical radius.
4. **Compute curvature sorting** — For each punctum, convert lipid amplitude
   to liposome radius, then compute protein surface density =
   protein_A / (4πR²).

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

Requires both lipid and protein channels.

```bash
python plot_curvature.py \
    --input data/.../filtered_puncta_A_values.txt \
    --dls-mean-diameter 80.11 \
    --lipid-col A_ch1 \
    --protein-col A_ch2 \
    --bin-width 0.5 \
    --save-dir figures/
```

| Argument              | Meaning |
|-----------------------|---------|
| `--input`             | Filtered puncta file from Step 3 (accepts multiple files) |
| `--dls-mean-diameter` | Mean liposome diameter in nm from DLS |
| `--lipid-col`         | Column for lipid amplitude (default: `A_ch1`) |
| `--protein-col`       | Column for protein amplitude (default: `A_ch2`) |
| `--bin-width`         | Radius bin width in nm for averaged curve (default: `0.5`) |
| `--save-dir`          | Output directory for figures |

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
| `--dls-mean-diameter` | Mean diameter from DLS. Omit to skip diameter plot |
| `--bins`              | Number of histogram bins (default: `80`) |
| `--transform`         | `raw`, `sqrt`, or `log_sqrt` (default: `raw`) |
| `--save-dir`          | Output directory for figures |

### Optional: DLS calibration

Compute mean diameter from a Malvern Zetasizer Excel export. The script
auto-detects the export format (Intensity, Volume, Number sections) and
uses the Number distribution by default, averaged across all measurement
records.

```bash
python dls_calibration.py --input data/dls_data.xlsx
```

| Argument           | Meaning |
|--------------------|---------|
| `--input`          | Path to Zetasizer `.xlsx` export |
| `--distribution`   | `number` (default), `volume`, or `intensity` |
| `--sheet`          | Sheet name or index (default: first sheet) |

Prints per-record weighted means and the overall average. Use that number
as `--dls-mean-diameter` in the plotting steps.

### Optional: Normalized overlay across experiments

Compare curvature sorting across conditions by normalizing each curve so
the largest-radius bin = 1 (fold-enrichment at high curvature). Each input
file gets its own DLS diameter.

```bash
python plot_overlay.py \
    --input data/cond1/filtered.txt:80.11 data/cond2/filtered.txt:95.0 \
    --labels "WT protein" "Mutant K58A" \
    --save-dir figures/
```

| Argument        | Meaning |
|-----------------|---------|
| `--input`       | Files with DLS diameters, as `file.txt:diameter` pairs |
| `--labels`      | Custom legend labels (default: parent folder name) |
| `--bin-width`   | Radius bin width in nm (default: `0.5`) |
| `--output-name` | Output filename (default: `normalized_curvature_overlay.png`) |
| `--save-dir`    | Output directory for figure |

## File Descriptions

| File                  | Purpose |
|-----------------------|---------|
| `prepare_input.py`    | Split and reorder TIFF channels, organize for MATLAB |
| `analyze_matlab.py`   | Read MATLAB detection `.mat` files, filter puncta, export TSV |
| `dls_calibration.py`  | Parse Malvern Zetasizer Excel export, compute mean diameter |
| `plot_curvature.py`   | Convert amplitudes to radii, plot protein density vs radius |
| `plot_histograms.py`  | Plot amplitude histograms and estimated diameter distributions |
| `plot_overlay.py`     | Overlay normalized curvature-sorting curves across experiments |
| `requirements.txt`    | Python dependencies: numpy, matplotlib, tifffile, h5py, openpyxl, pandas |

## Notes

- **MATLAB detection** is maintained separately
  ([DanuserLab/cmeAnalysis](https://github.com/DanuserLab/cmeAnalysis)).
  This pipeline reads its output: `detection_v2.mat` files containing a
  `frameInfo` struct with fields `A` (amplitude), `c` (background), and
  `hval_Ar` (hypothesis test).
- **Only the master/lipid channel** (ch1) gets a `Detection/` subfolder from
  CMEanalysis. The protein channel (ch2) will just contain the TIFF.
- **DLS calibration** uses the Number-weighted distribution by default,
  averaged across all measurement records. The ratio-of-means approach
  (DLS mean diameter / mean sqrt(lipid_A)) converts amplitudes to radii.
  A lognormal fit approach is planned — see the TODO in `dls_calibration.py`.
- **Lipid-only experiments** are fully supported. Run `analyze_matlab.py`
  with `--channels ch1 --lipid-channel ch1`, then use `plot_histograms.py`
  (skip `plot_curvature.py` since it requires protein data).
- The `data/` and `figures/` folders are gitignored. Your microscopy data
  stays local.
