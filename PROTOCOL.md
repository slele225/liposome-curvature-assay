# Single Liposome Curvature Assay — Image Analysis Protocol

**From raw microscopy images to curvature-sorting plots.**
A step-by-step protocol for users with no prior coding experience.

This document is the friendly, hand-holding companion to
[README.md](README.md). It walks through the entire pipeline end-to-end
with example paths, side-by-side macOS / Windows commands, and
troubleshooting. The README is the terse developer reference.

---

## What this pipeline does

Each liposome in your image is a sub-diffraction-limit fluorescent spot.
The lipid channel amplitude is proportional to surface area, so it
encodes the liposome's size. The protein channel amplitude tells you
how much protein is bound. By calibrating with DLS, you convert lipid
amplitude into a physical diameter, then compute protein density
= protein_A / (πD²) to see how protein binding depends on curvature.

**What goes in:**
1. Two-channel fluorescence TIFFs (lipid + protein) from the microscope.
2. A Malvern Zetasizer DLS export (`.xlsx`) of the liposome stock.

**What comes out:**
- A calibration factor (nm per √amplitude).
- A scatter + binned-mean plot of **protein surface density vs.
  liposome diameter**. Curves rising toward small diameters indicate
  positive-curvature sensing.
- Optional overlay plots comparing conditions (e.g. WT vs. mutant).

**Four pipeline stages:**
1. **Prepare images** — split microscope TIFFs into channels.
2. **MATLAB detection** — fit Gaussians to each spot (external CMEanalysis).
3. **Analyze MATLAB output** — filter puncta, export amplitudes.
4. **Calibrate and plot** — DLS overlay calibration, then plot.

---

## Quick Primer: Using the Terminal

Every step below involves typing commands into a terminal (also called
a "command line"). Open one like this:

- **macOS:** `Cmd + Space`, type `Terminal`, press Enter.
- **Windows:** Press the Windows key, type `PowerShell`, click
  **Windows PowerShell**.

A few things to know:

- `cd` means "change directory" — it's how you move into a folder.
  `cd Desktop` moves you into Desktop.
- Paste commands with `Cmd+V` (Mac) or `Ctrl+V` / right-click (Windows).
  You don't have to retype anything.
- Press Enter after each command to run it.
- If something goes wrong, read the error message — it usually tells
  you what's missing. The Troubleshooting section at the end covers
  the common ones.

> **Tip:** Many commands in this protocol span multiple lines, ending
> each line with `\` (Mac/Linux) or `` ` `` (Windows PowerShell). That
> is one command split across lines for readability — you can paste
> the whole block at once.

> **Tip:** You can type `cd ` (with a space), then drag-and-drop a
> folder from Finder / File Explorer into the terminal. The full path
> auto-fills.

> **Paths with spaces:** Wrap them in double quotes. Example:
> `cd "C:\Users\your name\my folder"`. This applies to every command
> in this protocol.

---

## Prerequisites

### Python

Check whether you already have Python:

| Windows (PowerShell) | Mac (Terminal) |
|----------------------|----------------|
| `python --version`   | `python --version` |

If you see `Python 3.9` or newer, you're fine. If not, **uv (next
section) will install Python for you** — you don't need to download
it separately. Skip to the uv step.

### uv (Python package manager)

uv is a fast, modern replacement for pip. It manages Python versions,
virtual environments, and packages in one tool.

**Mac / Linux (bash):**
```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

**Windows (PowerShell):**
```powershell
powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"
```

After installing, **close and reopen your terminal** so it picks up the
new command. Verify:

```
uv --version
```

You should see something like `uv 0.6.x`.

### MATLAB + CMEanalysis

Step 2 of the pipeline runs in MATLAB using the Danuser lab's
CMEanalysis toolbox. Skip ahead if MATLAB is already installed.

1. **Install MATLAB.** Penn students/faculty get free licenses —
   check your IT portal. During install, you'll be asked which
   toolboxes to add. Select these (they're required by CMEanalysis):
   - Curve Fitting Toolbox
   - Image Processing Toolbox
   - Optimization Toolbox
   - Statistics Toolbox
   - Parallel Computing Toolbox

2. **Install missing toolboxes** later via **Home → Add-Ons → Get Add-Ons**.
   Restart MATLAB after each.

3. **Download CMEanalysis** from
   <https://github.com/DanuserLab/cmeAnalysis>: click the green
   **Code** button → **Download ZIP**. Unzip somewhere stable on your
   computer (e.g. `~/Documents/cmeAnalysis-master/`).

### DLS data

Export the size distribution from the Malvern Zetasizer software into
an Excel spreadsheet that contains the standard `X Intensity`,
`X Volume`, and `X Number` section headers. Only the **number
distribution** is used for calibration, but the other sections are
fine to leave in. You do **not** need Z-average or PDI values — the
current pipeline uses a distribution-overlay method, not log-normal
fitting.

Place the file inside the repo at `data/dls/<your_dls_file>.xlsx`
(see Folder Structure below).

---

## Download the Code and Set Up

**You only need to do this once per computer.**

1. **Get the repo.** Either:
   - Clone with git:
     ```
     git clone https://github.com/slele225/liposome-curvature-assay.git
     ```
   - Or click the green **Code** button on
     <https://github.com/slele225/liposome-curvature-assay> →
     **Download ZIP**, then unzip.

2. **Move the folder somewhere stable** — Desktop or Documents is fine.
   These instructions assume Desktop.

3. **Open a terminal and navigate to the folder:**

   | Windows (PowerShell) | Mac (Terminal) |
   |----------------------|----------------|
   | `cd $HOME\Desktop\liposome-curvature-assay` | `cd ~/Desktop/liposome-curvature-assay` |

4. **Create a virtual environment.** This keeps the pipeline's Python
   packages isolated from anything else on your computer.

   ```
   uv venv
   ```

   You'll see `Creating virtual environment at: .venv`. A hidden
   `.venv/` folder now lives inside the repo.

5. **Activate the virtual environment.**

   **Windows (PowerShell):**
   ```powershell
   .\.venv\Scripts\Activate.ps1
   ```

   **Mac / Linux (bash):**
   ```bash
   source .venv/bin/activate
   ```

   Your prompt should now start with `(.venv)`. **You need to do this
   every time you open a new terminal** — the Quick Reference at the
   bottom is your every-time checklist.

6. **Install the Python packages.**

   ```
   uv pip install -r requirements.txt
   ```

   This installs numpy, matplotlib, tifffile, h5py, openpyxl, pandas,
   and scipy. Takes a few seconds.

7. **Create data folders.**

   | Windows (PowerShell) | Mac (Terminal) |
   |----------------------|----------------|
   | `mkdir data, figures` | `mkdir data figures` |

8. **Verify.**

   ```
   python plot_curvature.py --help
   ```

   You should see a help message listing all the options. If you see
   an error, check that `(.venv)` appears in your prompt — if not,
   redo step 5.

---

## Folder Structure: Where to Put Your Files

All your data and outputs live inside the repo folder. Throughout the
rest of this document, the examples use a fake experiment named
`20240315_DOPC_EGFP` (DOPC liposomes with EGFP, imaged on 2024-03-15)
and a DLS file `20240315_DOPC_LUV.xlsx`. Substitute your own names.

```
liposome-curvature-assay/
├── data/                                     ← YOUR DATA GOES HERE
│   ├── dls/
│   │   └── 20240315_DOPC_LUV.xlsx            ← DLS Excel export
│   ├── 20240315_DOPC_EGFP/                   ← raw TIFFs from microscope
│   │   ├── image001.tif
│   │   └── image002.tif
│   └── 20240315_DOPC_EGFP_matlab/            ← created by Step 1
│       └── 488nm_580V_561nm_500V/
│           ├── cell1/
│           │   ├── ch1/   ← lipid (master)
│           │   └── ch2/   ← protein
│           └── filtered_puncta_A_values.txt  ← created by Step 3
├── figures/                                  ← PLOTS SAVED HERE
├── prepare_input.py
├── analyze_matlab.py
├── dls_calibration.py
├── plot_curvature.py
├── plot_overlay.py
├── plot_histograms.py
├── plot_dls.py
├── plot_dls_comparison.py
├── requirements.txt
└── README.md
```

**`data/`** holds raw inputs and intermediate files. Each experiment
gets its own subfolder. **`figures/`** is where final plots are saved.
Both folders are gitignored — your data stays on your computer.

---

## Step 1: Prepare Images

**What this does:** Reads your raw multi-frame TIFFs from the
microscope, splits them into per-channel folders organized by detector
voltage (from the Olympus FV3000 metadata), and optionally
center-crops. The output structure is what MATLAB's CMEanalysis
expects.

**Before you start:**
- Put all the raw `.tif` files from one experiment into a subfolder of
  `data/` — for example, `data/20240315_DOPC_EGFP/`.
- Know which TIFF frame index is lipid and which is protein. Open one
  TIFF in Fiji/ImageJ if unsure.

**Commands** (with the venv activated and you in the repo folder):

**Mac / Linux (bash):**
```bash
python prepare_input.py \
    --input  data/20240315_DOPC_EGFP \
    --output data/20240315_DOPC_EGFP_matlab \
    --frames 2,0 \
    --crop 1
```

**Windows (PowerShell):**
```powershell
python prepare_input.py `
    --input  data\20240315_DOPC_EGFP `
    --output data\20240315_DOPC_EGFP_matlab `
    --frames 2,0 `
    --crop 1
```

**If your folder name has spaces,** wrap the paths in double quotes:

**Mac / Linux (bash):**
```bash
python prepare_input.py \
    --input  "data/march 15 DOPC EGFP" \
    --output "data/march 15 DOPC EGFP_matlab" \
    --frames 2,0 \
    --crop 1
```

**Windows (PowerShell):**
```powershell
python prepare_input.py `
    --input  "data\march 15 DOPC EGFP" `
    --output "data\march 15 DOPC EGFP_matlab" `
    --frames 2,0 `
    --crop 1
```

**Arguments:**

| Argument   | What it means |
|------------|---------------|
| `--input`  | Folder containing your raw `.tif` files. |
| `--output` | Folder to create with the split channels. Will be created. |
| `--frames` | Comma-separated frame order. `2,0` means frame 2 → ch1 (lipid), frame 0 → ch2 (protein). Frame indices start at 0. |
| `--crop`   | Center crop divisor. `1` = no crop. `2` = center quarter. Use if your image has dark borders. |

**Expected output:** A new folder
`data/20240315_DOPC_EGFP_matlab/` containing one subfolder per voltage
group (e.g. `488nm_580V_561nm_500V/`), each with `cell1/`, `cell2/`,
… and inside each cell `ch1/` (lipid) and `ch2/` (protein) folders
with one TIFF per cell. The script prints a summary at the end.

**What to check:** Open one of the `ch1` TIFFs in Fiji — it should be a
single lipid frame, not a stack.

---

## Step 2: Run MATLAB Detection

**What this does:** Fits a 2D Gaussian to every sub-diffraction
fluorescent spot in each cell's image and records its amplitude,
background, and a quality flag. This is **not** a Python script — it
runs inside MATLAB using the Danuser lab's CMEanalysis toolbox
(installed in Prerequisites).

**Steps in MATLAB:**

1. **Open MATLAB and navigate to the CMEanalysis software folder.** In
   the MATLAB command window, replace the path below with where you
   saved the unzipped CMEanalysis:
   ```matlab
   cd('/path/to/cmeAnalysis-master/software')
   addpath(genpath(pwd))
   ```

2. **Load the image folder.** Run:
   ```matlab
   data = loadConditionData();
   ```
   A folder picker opens. Select your voltage-group folder from Step 1
   — for example,
   `data/20240315_DOPC_EGFP_matlab/488nm_580V_561nm_500V/`.

3. **Select channels.** When prompted for the **first (master)
   channel**, pick the **lipid** folder (`ch1`). It will then ask for
   the remaining channels — pick `ch2` (protein).

4. **Enter fluorophore names.** For each channel, enter the marker:
   `texasred` for Texas Red, `egfp` for EGFP, etc. If a name is not
   recognized, MATLAB will ask for the emission wavelength in nm.

5. **Run the detection.** This estimates the Gaussian PSF sigma from
   the images themselves (the default mode):
   ```matlab
   runDetection(data);
   ```
   - If you've already run detection and want to redo it:
     ```matlab
     runDetection(data, 'Overwrite', true);
     ```
   - If you already know the PSF sigma from prior calibration:
     ```matlab
     runDetection(data, 'Overwrite', true, 'Sigma', sigma_value);
     ```
     Replace `sigma_value` with your number. Sigma is the standard
     deviation, not the variance.

6. **Wait for completion.** MATLAB shows progress and finishes in a
   few minutes at most. It reports detected sigma values per channel
   across all images and pops up some graphs — you can close them.

**Expected output:** Each `cell*/ch1/` folder gains a `Detection/`
subfolder containing `detection_v2.mat`. Only the master/lipid channel
gets this folder — `ch2/` keeps just the TIFF. Example:

```
data/20240315_DOPC_EGFP_matlab/488nm_580V_561nm_500V/
└── cell1/
    ├── ch1/
    │   └── Detection/
    │       └── detection_v2.mat   ← this file must exist
    └── ch2/
        └── (TIFF only, no Detection folder)
```

> **Note:** Each `.mat` contains a `frameInfo` struct with fields `A`
> (amplitude), `c` (background), and `hval_Ar` (hypothesis test). The
> next step reads these.

---

## Step 3: Analyze MATLAB Output

**What this does:** Reads every `detection_v2.mat` file in the voltage
folder, applies an intensity threshold to filter noise and poor fits,
and exports a tab-separated file containing the amplitude values for
every liposome that passed.

**Before you start:** Step 2 must have written `detection_v2.mat` into
each cell's `ch1/Detection/` folder.

**Commands** (venv activated, in the repo):

**Mac / Linux (bash):**
```bash
python analyze_matlab.py \
    --input data/20240315_DOPC_EGFP_matlab/488nm_580V_561nm_500V \
    --channels ch1,ch2 \
    --lipid-channel ch1 \
    --k-std 4 \
    --output-name filtered_puncta_A_values.txt
```

**Windows (PowerShell):**
```powershell
python analyze_matlab.py `
    --input data\20240315_DOPC_EGFP_matlab\488nm_580V_561nm_500V `
    --channels ch1,ch2 `
    --lipid-channel ch1 `
    --k-std 4 `
    --output-name filtered_puncta_A_values.txt
```

**Arguments:**

| Argument          | What it means |
|-------------------|---------------|
| `--input`         | Voltage-group folder from Step 1 (must contain `cellN/chX/Detection/`). |
| `--channels`      | Channel folder names in order, comma-separated. `ch1,ch2` for two-channel. Use `ch1` alone for lipid-only. |
| `--lipid-channel` | Which channel is the lipid (used as the threshold reference). |
| `--k-std`         | Filter strictness. Keep puncta with `A > mean(c) + k·std(c)`. `4` is fairly strict, `2.0` is the script default, `1.5` is loose. |
| `--output-name`   | Output filename. Saved inside the `--input` folder. |

**Lipid-only experiments** (no protein channel): replace
`--channels ch1,ch2` with `--channels ch1`, keep
`--lipid-channel ch1`.

**Expected output:** A tab-separated file at
`data/20240315_DOPC_EGFP_matlab/488nm_580V_561nm_500V/filtered_puncta_A_values.txt`
with columns `source_image`, `A_ch1` (lipid amplitude), `A_ch2`
(protein amplitude — if two-channel). Each row is one liposome. The
script prints a per-cell kept/total summary.

**What to check:** Thousands of rows with positive amplitudes. Fewer
than ~100 rows means the threshold is too strict or detection failed.

---

## Step 4: DLS Calibration

**What this does:** Computes the single scalar **conversion factor**
that maps √(lipid amplitude) to physical diameter in nm. The script
overlays the rebinned fluorescence √(A) histogram onto the DLS
number-weighted size distribution and finds the scale that minimizes
the chi-squared difference. This is the standard SLiC calibration
procedure (Kunding 2008, Hatzakis 2009, Bhatia 2009, Zeno 2018,
Johnson 2025). The simple ratio-of-means
(DLS-mean / mean(√A)) is also reported as a sanity check.

For background on the method, see
[DLS_Calibration_Notes.md](DLS_Calibration_Notes.md).

**Before you start:** Your DLS `.xlsx` must contain the standard
Zetasizer headers (`X Intensity`, `X Volume`, `X Number`). The number
distribution is the one used.

**Commands:**

**Mac / Linux (bash):**
```bash
python dls_calibration.py \
    --dls-input data/dls/20240315_DOPC_LUV.xlsx \
    --fluor-input data/20240315_DOPC_EGFP_matlab/488nm_580V_561nm_500V/filtered_puncta_A_values.txt \
    --save-dir figures/
```

**Windows (PowerShell):**
```powershell
python dls_calibration.py `
    --dls-input data\dls\20240315_DOPC_LUV.xlsx `
    --fluor-input data\20240315_DOPC_EGFP_matlab\488nm_580V_561nm_500V\filtered_puncta_A_values.txt `
    --save-dir figures\
```

**Arguments:**

| Argument          | What it means |
|-------------------|---------------|
| `--dls-input`     | Path to Zetasizer `.xlsx` export. |
| `--fluor-input`   | Path to `filtered_puncta_A_values.txt` from Step 3. |
| `--save-dir`      | Where to save the overlay plot (optional). |
| `--lipid-col`     | Column name for lipid amplitude (default: `A_ch1`). |
| `--bootstrap`     | Number of bootstrap resamples for variance estimate (default `0` = off). Try 200–500. |
| `--bootstrap-k`   | Puncta per bootstrap iteration (default: same as total). |

**Optional extras:**

- `--bootstrap 200` — resample 200 times to estimate confidence
  intervals on the conversion factor. Use this when you want to
  report an uncertainty.
- `--bootstrap-k 500` — use only 500 puncta per bootstrap iteration
  rather than all of them. Useful for checking stability under
  smaller sample sizes.
- `--lipid-col` — change if your filtered file uses a non-default
  column name.

**Expected output:**
- A line printed to the terminal like:
  ```
  python plot_curvature.py --conversion-factor 4.149 ...
  ```
  **Copy this number** — you need it for Step 5.
- `figures/dls_overlay_calibration.png` — a two-panel plot with the
  DLS number distribution (blue) and the rebinned fluorescence
  distribution (red) overlaid on a log diameter axis, plus a side
  panel showing the raw √(A) histogram with a nm scale.

**What to check:**
- Red and blue curves overlap well in the overlay panel.
- The "Implied mean diameter" printed matches the DLS number-weighted
  mean within a few %.
- A wildly different mean diameter (e.g. 3 nm or 3000 nm) means
  detection is picking up a spurious population — revisit Step 3
  threshold settings.

---

## Step 5: Plot Curvature Sorting

**What this does:** Converts every lipid amplitude into a physical
diameter using the conversion factor from Step 4, computes
`protein_density = protein_A / (πD²)`, and saves a scatter plot with
overlaid binned means.

**Before you start:** You need either the conversion factor (preferred)
or the implied mean diameter — both are printed by Step 4.

**Commands:**

**Mac / Linux (bash):**
```bash
python plot_curvature.py \
    --input data/20240315_DOPC_EGFP_matlab/488nm_580V_561nm_500V/filtered_puncta_A_values.txt \
    --conversion-factor 4.149 \
    --save-dir figures/
```

**Windows (PowerShell):**
```powershell
python plot_curvature.py `
    --input data\20240315_DOPC_EGFP_matlab\488nm_580V_561nm_500V\filtered_puncta_A_values.txt `
    --conversion-factor 4.149 `
    --save-dir figures\
```

The `4.149` is the exact number printed by `dls_calibration.py` —
substitute your own.

**Arguments:**

| Argument              | What it means |
|-----------------------|---------------|
| `--input`             | Filtered puncta file from Step 3 (accepts multiple files). |
| `--conversion-factor` | Maps √(A) to diameter in nm. From Step 4. |
| `--dls-mean-diameter` | Alternative to `--conversion-factor`. The "Implied mean diameter" from Step 4 (e.g. `95.75`). Equivalent result. |
| `--lipid-col`         | Column for lipid amplitude (default `A_ch1`). |
| `--protein-col`       | Column for protein amplitude (default `A_ch2`). |
| `--bin-width`         | Diameter bin width in nm for the binned-mean curve (default `0.5`). |
| `--diameter-cutoff`   | Exclude puncta with diameter above this value in nm. |
| `--y-pad`             | Y-axis padding factor around the binned means (default `0.3` = 30 %). |
| `--save-dir`          | Output directory for the PNG. |

Provide exactly one of `--conversion-factor` or `--dls-mean-diameter`.

**Optional extras:**

- `--diameter-cutoff 120` — drop puncta above 120 nm. Useful when a
  sparse noisy tail at large diameters is stretching the x-axis.
- `--y-pad 0.1` — shrink y-axis padding to 10 %. Useful when a few
  outliers squash the trend.
- `--bin-width 5` — wider bins → smoother but coarser averaged curve
  (default is `0.5` nm).
- `--dls-mean-diameter 95.75` — alternative to `--conversion-factor`.
- `--lipid-col` / `--protein-col` — change the column names if
  non-default.

**Expected output:** A PNG at
`figures/488nm_580V_561nm_500V__protein_density_vs_diameter.png`.
Faint dots are individual liposomes; larger dots are the binned means.

**Interpreting the plot:**

- **X-axis** = liposome diameter (nm). Smaller diameter = higher
  curvature.
- **Y-axis** = protein surface density (amplitude per nm²). Higher =
  more protein per unit area.
- **Binned mean rises toward small diameter** → your protein binds
  high-curvature membranes preferentially (positive curvature
  sensing).
- **Flat** → curvature-independent.
- **Falls toward small diameter** → protein avoids high curvature.

**Sanity check:** The x-axis range should roughly match the DLS size
distribution (typically 50–200 nm for LUVs). If the plot's x-axis
peaks at 3 nm or 3000 nm, the conversion factor is wrong — re-run
Step 4.

---

## Step 6: Optional — Overlay Multiple Conditions

**What this does:** Plots normalized curvature-sorting curves from
multiple experiments on a single figure (e.g. WT vs. mutant). By
default each curve is normalized so the largest-diameter bin = 1, so
the y-axis reads as fold-enrichment at high curvature relative to
flat membranes.

**Before you start:** Run Step 4 separately for each condition to get
its own conversion factor — different liposome stocks may have
different size distributions.

**Commands** (using two example conditions, WT EGFP and the K58A
mutant):

**Mac / Linux (bash):**
```bash
python plot_overlay.py \
    --input data/20240315_DOPC_EGFP_matlab/488nm_580V_561nm_500V/filtered_puncta_A_values.txt:4.149 \
            data/20240315_DOPC_K58A_matlab/488nm_580V_561nm_500V/filtered_puncta_A_values.txt:4.203 \
    --labels "WT EGFP" "Mutant K58A" \
    --save-dir figures/
```

**Windows (PowerShell):**
```powershell
python plot_overlay.py `
    --input data\20240315_DOPC_EGFP_matlab\488nm_580V_561nm_500V\filtered_puncta_A_values.txt:4.149 `
            data\20240315_DOPC_K58A_matlab\488nm_580V_561nm_500V\filtered_puncta_A_values.txt:4.203 `
    --labels "WT EGFP" "Mutant K58A" `
    --save-dir figures\
```

The format for each input is `path:conversion_factor` — note the
colon, with no spaces around it. `4.149` is the WT factor; `4.203` is
the K58A mutant. Use your own numbers.

**Arguments:**

| Argument          | What it means |
|-------------------|---------------|
| `--input`         | One or more `path:conversion_factor` pairs. |
| `--labels`        | Custom legend labels in matching order. Default: parent folder name of each file. |
| `--lipid-col`     | Column for lipid amplitude (default `A_ch1`). |
| `--protein-col`   | Column for protein amplitude (default `A_ch2`). |
| `--bin-width`     | Diameter bin width in nm (default `0.5`). |
| `--diameter-cutoff` | Exclude puncta above this diameter. |
| `--y-pad`         | Y-axis padding factor around the curves (default `0.3`). |
| `--normalize-to`  | `rightmost` (default), `leftmost`, `minimum`, or `none`. |
| `--output-name`   | Output filename (default `normalized_curvature_overlay.png`). |
| `--save-dir`      | Output directory. |

**Optional extras:**

- `--normalize-to leftmost` — normalize to the smallest-diameter bin
  instead of the largest. Useful for negative-curvature sensors.
- `--normalize-to none` — plot raw density instead of fold-enrichment.
- `--diameter-cutoff 120` — same purpose as in Step 5.
- `--y-pad 0.1` — shrink y-axis padding.
- `--bin-width 5` — wider bins → smoother curves.
- `--output-name WT_vs_K58A.png` — change the saved filename.
- `--lipid-col` / `--protein-col` — change column names if non-default.

**Expected output:** A PNG at
`figures/normalized_curvature_overlay.png` (or your `--output-name`).
A dashed line at `y = 1` marks the "flat" reference. Curves rising
above it on the small-diameter side indicate curvature sensing.

---

## Optional Utility: plot_histograms.py

A quick sanity-check plot of amplitude distributions and (if a
conversion factor is provided) estimated diameter distribution.
Useful before committing to Step 5, and for lipid-only experiments.

**Mac / Linux (bash):**
```bash
python plot_histograms.py \
    --input data/20240315_DOPC_EGFP_matlab/488nm_580V_561nm_500V/filtered_puncta_A_values.txt \
    --lipid-col A_ch1 \
    --protein-col A_ch2 \
    --conversion-factor 4.149 \
    --save-dir figures/
```

**Windows (PowerShell):**
```powershell
python plot_histograms.py `
    --input data\20240315_DOPC_EGFP_matlab\488nm_580V_561nm_500V\filtered_puncta_A_values.txt `
    --lipid-col A_ch1 `
    --protein-col A_ch2 `
    --conversion-factor 4.149 `
    --save-dir figures\
```

For lipid-only data, omit `--protein-col`. Add `--bins 100` to change
histogram resolution. The estimated diameter distribution panel is
especially useful — if it doesn't resemble your DLS, something is off
with calibration or filtering.

---

## Optional Utility: DLS visualization

Two scripts for sanity-checking the DLS data itself.

**Plot the DLS distribution alone:**

**Mac / Linux (bash):**
```bash
python plot_dls.py data/dls/20240315_DOPC_LUV.xlsx
```

**Windows (PowerShell):**
```powershell
python plot_dls.py data\dls\20240315_DOPC_LUV.xlsx
```

Add `--log` for log-diameter x-axis, `--distribution intensity` to
plot the intensity distribution instead of number, or `--zoom-pct 100`
to disable zooming.

**Compare DLS vs. fluorescence side-by-side:**

**Mac / Linux (bash):**
```bash
python plot_dls_comparison.py \
    --dls-input data/dls/20240315_DOPC_LUV.xlsx \
    --fluor-input data/20240315_DOPC_EGFP_matlab/488nm_580V_561nm_500V/filtered_puncta_A_values.txt \
    --channels 0 Lipid 1 EGFP \
    --zoom-pct 95 --bins 200 --save-dir figures/
```

**Windows (PowerShell):**
```powershell
python plot_dls_comparison.py `
    --dls-input data\dls\20240315_DOPC_LUV.xlsx `
    --fluor-input data\20240315_DOPC_EGFP_matlab\488nm_580V_561nm_500V\filtered_puncta_A_values.txt `
    --channels 0 Lipid 1 EGFP `
    --zoom-pct 95 --bins 200 --save-dir figures\
```

`--channels` takes pairs of `index label`. Index `0` corresponds to
`A_ch1`, `1` to `A_ch2`, etc. Each gets its own panel. Useful for
confirming your detection is picking up the same population the DLS
measured.

---

## Quick Reference: Every-Time Checklist

After the initial setup, every new dataset:

1. Open a terminal.
2. Navigate to the repo:
   `cd ~/Desktop/liposome-curvature-assay` (Mac) or
   `cd $HOME\Desktop\liposome-curvature-assay` (Windows).
3. Activate the venv:
   `source .venv/bin/activate` (Mac) or
   `.\.venv\Scripts\Activate.ps1` (Windows).
   Confirm `(.venv)` appears in the prompt.
4. Copy raw TIFFs into a new subfolder of `data/`.
5. Run `prepare_input.py` (Step 1).
6. Run MATLAB CMEanalysis (Step 2).
7. Run `analyze_matlab.py` (Step 3).
8. Run `dls_calibration.py` (Step 4) — copy the conversion factor it
   prints.
9. Run `plot_curvature.py` with that conversion factor (Step 5).
10. Find your plot in `figures/`.

---

## Troubleshooting

**"'python' is not recognized" / "command not found: python"**
→ The virtual environment isn't activated. Run the activation command
from Setup step 5 — your prompt should show `(.venv)` when active.

**"'uv' is not recognized" / "command not found: uv"**
→ Close and reopen your terminal after installing uv. If it's still
not found, reinstall using the commands in Prerequisites.

**Activation script blocked on Windows
("running scripts is disabled on this system")**
→ PowerShell's default execution policy blocks the script. Run:
```
Set-ExecutionPolicy -ExecutionPolicy RemoteSigned -Scope CurrentUser
```
then retry the activate command.

**"No module named tifffile" / pandas / numpy / etc.**
→ The venv isn't activated. Run the activation command. If that
doesn't fix it, redo `uv pip install -r requirements.txt`.

**"FileNotFoundError" on a path you typed**
→ Either the path is wrong, or the previous step didn't write the
file. Paste the path into Finder / File Explorer to confirm. If it
contains spaces, wrap it in double quotes.

**"No `.tif` files found"**
→ `--input` should point to the folder containing TIFFs, not a parent
or subfolder.

**"Frames out of range"**
→ The `--frames` indices don't match the TIFF. Open one in Fiji to
check how many frames it has — frame indices start at `0`.

**MATLAB detection produced no `Detection/` folder**
→ Check that you set `ch1` as the **master** channel in CMEanalysis.
Only the master channel gets a `Detection/` folder.

**"missing detection file" in Step 3**
→ MATLAB Step 2 hasn't been run on that cell/channel yet, or
overwrote in a different folder. Re-run Step 2.

**"No puncta passed the filters"**
→ Threshold too strict or signal is very low. Try `--k-std 2.0` or
`--k-std 1.5`.

**`dls_calibration.py` prints a tiny (< 10 nm) or huge (> 1000 nm)
implied mean diameter**
→ Step 3 output or DLS file is misformed. Run `plot_histograms.py`
to inspect the fluorescence distribution — it should be roughly
unimodal and positive.

**Step 5 plot's x-axis range doesn't match the DLS mean**
→ You used the wrong number for `--conversion-factor`. Copy the exact
value printed by `dls_calibration.py`, not a derived one.

**Very noisy plot in Step 5**
→ Increase `--bin-width` (try `1.0` or `5`). Or recheck calibration.

**"UnicodeDecodeError" / weird text errors reading the DLS file**
→ The `.xlsx` may have been opened and re-saved as `.csv` somewhere.
Re-export a fresh `.xlsx` from the Zetasizer software.

**"Permission denied" on macOS for `.venv`**
→ Delete the `.venv` folder and recreate with `uv venv`.

**"I closed my terminal and now nothing works"**
→ You need to `cd` back to the repo and re-activate the venv. See
the Quick Reference checklist above.

---

## Further reading

- [README.md](README.md) — concise reference for every command.
- [DLS_Calibration_Notes.md](DLS_Calibration_Notes.md) — theoretical
  background on the distribution-overlay calibration: why it works,
  how it relates to published SLiC protocols (Kunding 2008, Hatzakis
  2009, Bhatia 2009, Zeno 2018, Johnson 2025), and what the reported
  statistics mean.
