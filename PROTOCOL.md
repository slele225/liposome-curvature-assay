# SLiC Analysis Protocol — Step-by-Step

A beginner-friendly walkthrough for running the Single Liposome Curvature
(SLiC) analysis pipeline end-to-end. Written for wet lab users with no
prior Python or command-line experience.

If you are comfortable with Python and the command line, use
[README.md](README.md) instead — this document covers the same commands
but with extra hand-holding.

---

## Introduction

This pipeline turns raw fluorescence microscope images of liposomes into
a plot showing how a protein of interest sorts onto membranes of
different curvature.

**What you put in:**
1. Two-channel TIFF images from the microscope — one channel is the lipid
   dye, the other is the protein of interest.
2. A dynamic light scattering (DLS) spreadsheet (`.xlsx`) from the Malvern
   Zetasizer, measuring the true size distribution of your liposome stock.

**What you get out:**
- A calibration factor (nm per √amplitude).
- A plot of **protein surface density vs. liposome diameter**, with each
  dot a single liposome. If your protein senses curvature, the curve
  rises toward small diameters (high curvature).
- Optional overlay plots comparing multiple experimental conditions
  (e.g., wild-type vs. mutant protein).

**What happens in the middle:**
- MATLAB code (external, see Step 2) fits a 2D Gaussian to each bright
  spot and reports its amplitude (the brightness of that spot).
- Python scripts filter those spots, compare to DLS, and produce plots.

---

## Setup

You only need to do this **once per computer**.

### 1. Install Python

1. Go to <https://www.python.org/downloads/> and download Python 3.9 or
   newer.
2. **Windows:** During install, check the box **"Add python.exe to PATH"**.
   This step is easy to miss and causes most "Python is not recognized"
   errors.
3. **Mac:** Just run the installer.

Verify by opening a terminal and running:

| Windows (PowerShell) | Mac (Terminal) |
|----------------------|----------------|
| `python --version`   | `python3 --version` |

You should see something like `Python 3.11.4`.

### 2. Open a terminal

- **Windows:** Press `Win + X`, then click **Terminal** or **PowerShell**.
- **Mac:** Press `Cmd + Space`, type `Terminal`, press Enter.

### 3. Navigate to the repo folder

Use `cd` ("change directory") to move into the folder where the code
lives. Replace `<path/to/liposome-curvature-assay>` with the actual
folder path.

| Windows (PowerShell) | Mac (Terminal) |
|----------------------|----------------|
| `cd C:\Users\yourname\liposome-curvature-assay` | `cd /Users/yourname/liposome-curvature-assay` |

**Tip:** If your path has spaces, wrap it in quotes:
`cd "C:\Users\your name\folder with spaces"`.

### 4. Install required packages

From inside the repo folder, run:

| Windows (PowerShell) | Mac (Terminal) |
|----------------------|----------------|
| `pip install -r requirements.txt` | `pip3 install -r requirements.txt` |

This installs numpy, pandas, matplotlib, scipy, tifffile, h5py, and
openpyxl. Takes about a minute on a fresh install.

### 5. Create data folders

| Windows (PowerShell) | Mac (Terminal) |
|----------------------|----------------|
| `mkdir data, figures` | `mkdir data figures` |

Copy your TIFF images into `data/<your_experiment_name>/` (for example,
`data/20240315_DOPC_EGFP/`) and your DLS spreadsheet into
`data/dls/`.

### 6. Confirm setup works

Run:

```
python plot_curvature.py --help
```

If you see a list of arguments, you are ready. If you see an error, see
the Troubleshooting section below.

---

## Running the pipeline

Each step produces an output that feeds into the next step. Run them
in order.

> **The paths shown below are examples — replace them with the paths to
> YOUR files.** All examples use a made-up experiment called
> `20240315_DOPC_EGFP` (DOPC liposomes with EGFP protein, imaged on
> 2024-03-15) and a DLS file `20240315_DOPC_LUV.xlsx`. Substitute your
> own filenames throughout. **Keep the double quotes if your paths
> contain spaces** — e.g. a Windows path with a space in it looks like
> `"data\march 15 experiment\filtered_puncta_A_values.txt"`.

### Step 1: Prepare TIFF images

**What this does:** Splits multi-frame microscope TIFFs into separate
per-channel folders, organized by detector voltage. This is the format
MATLAB's CMEanalysis expects.

**Command:**

| Windows (PowerShell) | Mac (Terminal) |
|----------------------|----------------|
| `python prepare_input.py --input data\20240315_DOPC_EGFP --output data\20240315_DOPC_EGFP_matlab --frames 2,0 --crop 1` | `python3 prepare_input.py --input data/20240315_DOPC_EGFP --output data/20240315_DOPC_EGFP_matlab --frames 2,0 --crop 1` |

**If your folder name has spaces,** wrap each path in double quotes:

| Windows | Mac |
|---------|-----|
| `python prepare_input.py --input "data\march 15 DOPC EGFP" --output "data\march 15 DOPC EGFP_matlab" --frames 2,0 --crop 1` | `python3 prepare_input.py --input "data/march 15 DOPC EGFP" --output "data/march 15 DOPC EGFP_matlab" --frames 2,0 --crop 1` |

**Arguments explained:**
- `--input` — the folder with your raw TIFFs (here,
  `data/20240315_DOPC_EGFP/` containing files like
  `image001.tif`, `image002.tif`, …).
- `--output` — a new folder (will be created).
- `--frames 2,0` — frame 2 is lipid, frame 0 is protein. Adjust if your
  acquisition is different.
- `--crop 1` — no cropping. Use `2` for center-quarter crop.

**Output to look for:** A new folder `data/20240315_DOPC_EGFP_matlab/`
containing subfolders like `488nm_580V_561nm_500V/cell1/ch1/` and
`.../ch2/`, each with one TIFF per cell.

**What to check:** Open one of the `ch1` TIFFs in Fiji/ImageJ — it
should look like a single lipid-channel frame, not a stack.

---

### Step 2: Run CMEanalysis in MATLAB (external)

**What this does:** Fits 2D Gaussians to each sub-diffraction bright
spot (punctum) and records the amplitude, background, and a quality
flag for each.

**This is not a Python script.** You run it inside MATLAB using the
CMEanalysis toolbox from the Danuser lab:
<https://github.com/DanuserLab/cmeAnalysis>

1. Install CMEanalysis in MATLAB per their instructions.
2. Point CMEanalysis at your
   `data/20240315_DOPC_EGFP_matlab/488nm_580V_561nm_500V/` folder from
   Step 1.
3. Set the **master channel** to `ch1` (lipid).
4. Run the detection step.

**Output to look for:** A file called `detection_v2.mat` appears inside
each cell's `ch1/Detection/` subfolder. Only the master channel (ch1)
gets a `Detection/` folder — this is expected.

**What to check:** Open one `detection_v2.mat` in MATLAB and confirm
the `frameInfo` struct contains fields named `A`, `c`, and `hval_Ar`.

---

### Step 3: Filter puncta

**What this does:** Reads all the MATLAB detection files, keeps only
puncta that pass an intensity threshold, and writes a single tab-
separated file with the kept amplitudes.

**Command (two-channel, lipid + protein):**

| Windows | Mac |
|---------|-----|
| `python analyze_matlab.py --input data\20240315_DOPC_EGFP_matlab\488nm_580V_561nm_500V --channels ch1,ch2 --lipid-channel ch1 --k-std 4 --output-name filtered_puncta_A_values.txt` | `python3 analyze_matlab.py --input data/20240315_DOPC_EGFP_matlab/488nm_580V_561nm_500V --channels ch1,ch2 --lipid-channel ch1 --k-std 4 --output-name filtered_puncta_A_values.txt` |

**Lipid-only command** (no protein channel): replace
`--channels ch1,ch2` with `--channels ch1` and keep
`--lipid-channel ch1`.

**Arguments explained:**
- `--k-std 4` — keep puncta with amplitude > mean(background) + 4·std.
  Raise for stricter filtering, lower for looser (default in the script
  is `2.0`).

**Output to look for:** A file called `filtered_puncta_A_values.txt`
inside the input folder. Columns: `A_ch1`, `A_ch2` (if two-channel).

**What to check:** Open the file in Excel. It should have thousands of
rows, all with positive amplitudes. Very low row counts (< 100) mean
the threshold is too strict or detection failed.

---

### Step 4: DLS calibration

**What this does:** Reads your Zetasizer spreadsheet and the filtered
puncta file, then finds a single number (the *conversion factor*) that
maps √amplitude to diameter in nm.

**Before running:** In Excel, confirm your DLS `.xlsx` file has the
standard Zetasizer sections — headers `X Intensity`, `X Volume`, and
`X Number` in column A. Only the number distribution is used.

**Command:**

| Windows | Mac |
|---------|-----|
| `python dls_calibration.py --dls-input data\dls\20240315_DOPC_LUV.xlsx --fluor-input data\20240315_DOPC_EGFP_matlab\488nm_580V_561nm_500V\filtered_puncta_A_values.txt --save-dir figures\` | `python3 dls_calibration.py --dls-input data/dls/20240315_DOPC_LUV.xlsx --fluor-input data/20240315_DOPC_EGFP_matlab/488nm_580V_561nm_500V/filtered_puncta_A_values.txt --save-dir figures/` |

**Output to look for:**
- Printed to the terminal: a line like
  `python plot_curvature.py --conversion-factor 4.149 ...`
  **Copy this number.** You need it for Step 5.
- `figures/dls_overlay_calibration.png` — a two-panel plot showing the
  DLS distribution (blue) and your fluorescence distribution (red,
  dashed) on the same axis.

**What to check:**
- The red (fluorescence) and blue (DLS) curves should overlap well.
- The "Implied mean diameter" printed should be close to the DLS
  number-weighted mean (usually within a few %).
- If they disagree badly, your detection threshold in Step 3 may be
  wrong — too many spurious detections bias the distribution.

**Optional extras:**
- `--bootstrap 200` — resamples the fluorescence data 200 times to
  estimate confidence intervals on the conversion factor. Use this when
  you want to report an uncertainty.
- `--bootstrap-k 500` — use only 500 puncta per bootstrap iteration
  instead of all of them. Useful for checking whether your result is
  stable under smaller samples.
- `--lipid-col A_ch1` — change the column name used for lipid amplitude
  if your filtered file does not use the default `A_ch1`.

---

### Step 5: Plot curvature sorting

**What this does:** Takes the conversion factor from Step 4 and makes
the final protein-density-vs-diameter plot.

**Command:**

| Windows | Mac |
|---------|-----|
| `python plot_curvature.py --input data\20240315_DOPC_EGFP_matlab\488nm_580V_561nm_500V\filtered_puncta_A_values.txt --conversion-factor 4.149 --save-dir figures\` | `python3 plot_curvature.py --input data/20240315_DOPC_EGFP_matlab/488nm_580V_561nm_500V/filtered_puncta_A_values.txt --conversion-factor 4.149 --save-dir figures/` |

The `4.149` is the exact number printed by `dls_calibration.py` in
Step 4 — replace it with whatever your own run printed.

**Output to look for:** A file in `figures/` ending with
`__protein_density_vs_diameter.png`. Light dots are individual puncta;
larger dots are bin-averaged means.

**What to check:**
- The x-axis range should roughly match the DLS size distribution
  (typically 50–200 nm). If you see a peak at 3 nm or 3000 nm, the
  conversion factor is wrong — re-run Step 4.

**Optional extras:**
- `--diameter-cutoff 120` — exclude puncta with diameter above 120 nm.
  Use this when a sparse noisy tail at large diameters is stretching
  the x-axis.
- `--y-pad 0.1` — shrink the y-axis padding from 30% to 10%. Use this
  when the trend is hard to see because a few y-axis outliers are
  squashing the plot.
- `--bin-width 5` — change the bin width for the averaged curve
  (default is 0.5 nm). Larger bins = smoother curve, fewer points.
- `--dls-mean-diameter 95.75` — alternative to `--conversion-factor`.
  Use the mean diameter printed by Step 4 instead of the conversion
  factor. Results are equivalent.
- `--lipid-col A_ch1` / `--protein-col A_ch2` — change the column
  names used for lipid/protein amplitudes if your filtered file does
  not use the defaults.

---

### Optional Step 6: Overlay multiple conditions

**What this does:** Plots normalized curvature-sorting curves from
several experiments on one figure (e.g., WT vs. mutant).

**Command:**

| Windows | Mac |
|---------|-----|
| `python plot_overlay.py --input data\20240315_DOPC_EGFP_matlab\488nm_580V_561nm_500V\filtered_puncta_A_values.txt:4.149 data\20240315_DOPC_K58A_matlab\488nm_580V_561nm_500V\filtered_puncta_A_values.txt:4.203 --labels "WT EGFP" "Mutant K58A" --save-dir figures\` | `python3 plot_overlay.py --input data/20240315_DOPC_EGFP_matlab/488nm_580V_561nm_500V/filtered_puncta_A_values.txt:4.149 data/20240315_DOPC_K58A_matlab/488nm_580V_561nm_500V/filtered_puncta_A_values.txt:4.203 --labels "WT EGFP" "Mutant K58A" --save-dir figures/` |

Each input is `path:conversion_factor` (note the colon, with no spaces
around it) — **run Step 4 separately for each condition** to get its
own conversion factor. Here `4.149` is for WT and `4.203` is for the
K58A mutant — use your own numbers.

**Output to look for:** `figures/normalized_curvature_overlay.png`.

**What to check:** All curves should end at `y = 1` on the right side
(the "flat" reference). Rising curves on the left mean curvature
sorting; flat curves mean no sorting.

**Optional extras:**
- `--normalize-to leftmost` — normalize to the smallest-diameter bin
  instead of the largest. Useful if your protein is a
  negative-curvature sensor (binds flat or convex membranes more than
  curved ones). Options: `rightmost` (default), `leftmost`, `minimum`,
  or `none` (raw density values).
- `--diameter-cutoff 120` — exclude puncta with diameter above 120 nm.
  Same purpose as in Step 5.
- `--y-pad 0.1` — shrink y-axis padding from 30% to 10%.
- `--bin-width 5` — change the bin width (default 0.5 nm).
- `--output-name WT_vs_K58A.png` — change the saved figure's filename
  (default is `normalized_curvature_overlay.png`).
- `--lipid-col A_ch1` / `--protein-col A_ch2` — change column names if
  your filtered files do not use the defaults.

---

## Troubleshooting

**"'python' is not recognized as an internal or external command"**
→ On Windows, Python was not added to PATH during install. Reinstall
Python and check the **"Add python.exe to PATH"** box. Or use the full
path to `python.exe`.

**"No module named pandas" (or numpy, scipy, etc.)**
→ You skipped Setup step 4. Run
`pip install -r requirements.txt` from inside the repo folder.

**"FileNotFoundError: [...] filtered_puncta_A_values.txt"**
→ Either the path is wrong, or Step 3 did not write the file. Paste
your path into File Explorer / Finder to confirm it exists. If the
path contains spaces, wrap it in double quotes.

**"'python3' is not recognized" on Windows**
→ On Windows, use `python` (not `python3`). The examples above use
`python3` only for Mac.

**`dls_calibration.py` runs but prints a tiny (< 10 nm) or huge
(> 1000 nm) implied mean diameter**
→ Something is wrong with your Step 3 output or DLS file. Check the
fluorescence histogram with `plot_histograms.py` — it should be
roughly unimodal and positive.

**Step 5 plot's x-axis range doesn't match the DLS mean**
→ You are applying the conversion factor the wrong way. Make sure you
copied the exact number printed by `dls_calibration.py` into
`--conversion-factor`, not a derived value.

**MATLAB detection produced no `Detection/` folder**
→ Check that ch1 is set as the master channel in CMEanalysis. Only the
master channel gets a `Detection/` folder.

**"UnicodeDecodeError" or weird text errors**
→ Your DLS `.xlsx` may have been opened in Google Sheets and saved as
`.csv`. Re-export from the Malvern software as a fresh `.xlsx`.

---

## Expected outputs

At the end of a run, `figures/` should contain at minimum:

| File | From which step |
|------|-----------------|
| `dls_overlay_calibration.png` | Step 4 |
| `488nm_580V_561nm_500V__protein_density_vs_diameter.png` | Step 5 |
| `normalized_curvature_overlay.png` (if you ran Step 6) | Step 6 |

And `data/20240315_DOPC_EGFP_matlab/488nm_580V_561nm_500V/` contains:

| File | From which step |
|------|-----------------|
| `cell1/ch1/.../detection_v2.mat` (one per cell) | Step 2 |
| `filtered_puncta_A_values.txt` | Step 3 |

---

## Further reading

- [README.md](README.md) — concise reference for every command, for
  users who prefer terse documentation.
- [DLS_Calibration_Notes.md](DLS_Calibration_Notes.md) — theoretical
  background on the distribution-overlay calibration: why it works,
  how it relates to published SLiC protocols (Kunding 2008, Hatzakis
  2009, Bhatia 2009, Zeno 2018, Johnson 2025), and what the reported
  statistics mean.
