# Single Liposome Curvature Assay Pipeline

A Python pipeline for analyzing protein curvature sensing on liposomes from
fluorescence microscopy images. Sub-diffraction-limit spot detection uses
the CMEanalysis detection algorithm (Danuser Lab). The repository ships a
**native C++ CMEanalysis-compatible detector** (`native_detection/`,
recommended — no MATLAB required, see
[Native CMEanalysis-compatible detection](#native-cmeanalysis-compatible-detection-recommended))
and still supports the original
[CMEanalysis MATLAB code](https://github.com/DanuserLab/cmeAnalysis) as the
reference/fallback route.

> **New to this repo?** See [PROTOCOL.pdf](PROTOCOL.pdf) for a
> beginner-friendly step-by-step walkthrough with troubleshooting.

## For Non-Coders

If you're a wet-lab user running this pipeline for the first time and
you've never used Python or the command line, start with
**[PROTOCOL.pdf](PROTOCOL.pdf)** — a step-by-step walkthrough with
Windows/Mac commands side-by-side, placeholder paths, expected outputs,
and common-error troubleshooting. This README is the terse developer
reference; PROTOCOL.pdf (built from PROTOCOL.tex) is the friendly one.

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

### 6. Build the native detection backend (recommended)

This compiles `cme_detect`, the native spot detector, so that Step 2 of
the pipeline needs no MATLAB. On Windows install the free
[Build Tools for Visual Studio 2022](https://visualstudio.microsoft.com/downloads/)
with the "Desktop development with C++" workload (this includes CMake),
then:

```powershell
cd native_detection
cmake -S . -B build -G "Visual Studio 17 2022" -A x64
cmake --build build --config Release
cd ..
.\native_detection\build\Release\cme_detect.exe --help
```

Details, other platforms and troubleshooting: [native_detection/README.md](native_detection/README.md).
If you cannot build it, the MATLAB route (Step 2, fallback) still works.

### 7. Add your data

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
2. **Detect spots** — the CMEanalysis detection algorithm fits Gaussians to
   each sub-diffraction punctum and reports amplitudes. Run it with the
   native `cme_detect` backend in this repo (recommended) or with the
   original MATLAB CMEanalysis code (reference/fallback).
3. **Calibrate sizes** — Dynamic light scattering (DLS) gives the true size
   distribution of the liposome stock. By comparing the mean DLS diameter to
   the mean sqrt(lipid amplitude), you get a conversion factor from
   fluorescence amplitude to physical diameter.
4. **Compute curvature sorting** — For each punctum, convert lipid amplitude
   to liposome diameter, then compute protein surface density =
   protein_A / (πD²).

## Native CMEanalysis-compatible detection (recommended)

`native_detection/` contains `cme_detect`, a native C++ implementation of
the specific CMEanalysis detection path this project uses — the
dependency closure of the default `loadConditionData` + `rng(1)` +
`runDetection(data)` workflow (data-driven PSF sigma estimation,
`pointSourceDetection`, the `fitGaussian2D` MEX fitting, master/slave
channel fitting and the `frameInfo`/mask outputs). It is **not** a port of
all of CMEanalysis (no tracking, lifetime analysis, GUI or 3-D).

* It was **regression-validated against the original MATLAB CMEanalysis
  implementation**. On the tested data sets it reproduced the detection
  and filter decisions: identical detection counts and row order, identical
  `hval_Ar`/`hval_AD`/`isPSF` and mask decisions, master-channel positions
  and amplitudes within numerical tolerance, and the downstream SLiC filter
  kept exactly the same puncta (45,362 of 48,325 on a 10-cell reference
  set; 9,960 of 25,395 on a second, practical SLiC data set). The only
  differences are negligible floating-point deviations in a small number
  of numerically unstable fits. This is not a claim of bit-for-bit
  equivalence; see
  [native_detection/README.md](native_detection/README.md) and
  [PORTING_NOTES.md](native_detection/PORTING_NOTES.md) for the exact
  numbers.
* It is **substantially faster** and avoids the MATLAB / Parallel
  Computing Toolbox dependency. Runtime depends on image count, CPU and
  thread count; measured examples: the 10-cell reference benchmark took
  1903 s in MATLAB versus 366 s natively with 12 threads, and a practical
  SLiC condition ran in roughly 25 s natively versus several minutes with
  the prior MATLAB workflow.
* The original MATLAB workflow remains available as the reference
  implementation (Step 2, fallback, below).

Native pipeline:

```text
raw microscopy TIFFs
        ↓
prepare_input.py            (split channels, write channel_map.json)
        ↓
cme_detect.exe              (native_detection/, --channels <master first> --master <master>)
        ↓
analyze_cpp.py              (same filter as analyze_matlab.py)
        ↓
filtered_puncta_A_values.txt
        ↓
rest of SLiC analysis       (dls_calibration.py, plot_curvature.py, ...)
```

### Build (Windows, Visual Studio 2022)

```powershell
cd "<repo>\native_detection"

cmake -S . -B build `
  -G "Visual Studio 17 2022" `
  -A x64

cmake --build build --config Release
```

The executable is produced at `native_detection\build\Release\cme_detect.exe`
(with `gsl.dll` / `gslcblas.dll` copied next to it). The vendored
dependencies in `native_detection/third_party/` (GSL 2.8, libtiff 4.7.0,
GSL 1.16 Levenberg–Marquardt sources) are all that is needed on Windows;
other platforms use `-DCME_USE_SYSTEM_LIBS=ON`. `build/` is gitignored.

### Run

```powershell
& ".\native_detection\build\Release\cme_detect.exe" `
    --input "C:\path\to\prepared_condition" `
    --channels ch1,ch2 `
    --master ch1 `
    --seed 1 `
    --threads 12
```

```text
cme_detect --input <condition dir> --channels ch1,ch2[,...] --master ch1 [options]

  --input DIR          condition directory (contains the cell*/ movie folders)
  --channels a,b,...   channel folder names in the EXACT order loadConditionData received them
  --master NAME        must equal the first entry of --channels (master/source channel)
  --markers a,b,...    optional fluorophore names (metadata only; do not affect detection)
  --seed N             rng(N) seed used by the PSF sigma estimation (default 1)
  --sigma s1,s2,...    skip data-driven sigma estimation and use these values (runDetection 'Sigma')
  --output DIR         directory for detections_all.tsv / summary.tsv / sigma.tsv (default: <input>/cme_detect_output)
  --dump-dir DIR       write all intermediates for regression testing
  --no-masks           do not write Detection/dmasks.tif
  --no-matlab-layout   do not write Detection/detection_cpp.tsv under the master channel
  --movie-selector S   loadConditionData 'MovieSelector' (default 'cell')
  --threads N          OpenMP threads for frame-level parallelism (default 1)
```

* `--channels` is the exact analysis order; the **first entry is the
  master/source channel** and `--master` must equal it. The tool refuses
  any other arrangement rather than reordering silently. Folder names
  carry no master/slave meaning: `--channels ch2,ch1 --master ch2` is just
  as valid when `ch2` holds the lipid channel (check `channel_map.json`).
* `--seed 1` is the reproducible equivalent of the `rng(1)` used in the
  MATLAB workflow for the data-derived sigma estimation.
* `--threads` can be set to the number of cores available; results do not
  depend on it (output was byte-identical across thread counts).
* Output goes to `<input>\cme_detect_output` unless `--output` is given;
  in addition `<cell>\<master>\Detection\detection_cpp.tsv` is written in
  the MATLAB layout, which is what `analyze_cpp.py` reads.

Then filter and export, exactly as `analyze_matlab.py` would:

```powershell
python analyze_cpp.py `
    --input "C:\path\to\prepared_condition" `
    --channels ch1,ch2 `
    --lipid-channel ch1
```

`analyze_cpp.py` takes the same flags as `analyze_matlab.py` (`--k-std`,
`--hval-filter`, `--output-name`, `--raw-output-name`), imports the filter
and table writers from it, and writes the same
`filtered_puncta_A_values.txt` / `raw_puncta_values.txt` files, so every
downstream script works unchanged. It also checks that `--channels`
matches the order recorded in the `cme_detect` output.

### Complete native workflow

1. **Prepare microscope TIFFs** — `python prepare_input.py --input ... --output ... --crop 1`
   (Step 1 below). Read the `channel_map.json` it writes.
2. **Identify the master channel** from `channel_map.json` and the
   recorded wavelength/voltage/dye metadata. The lipid channel is the
   master; it is *not* necessarily `ch1`.
3. **Run native detection** — `cme_detect.exe --input <condition> --channels <master>,<slave> --master <master> --seed 1 --threads N`.
4. **Filter/export detections** — `python analyze_cpp.py --input <condition> --channels <master>,<slave> --lipid-channel <master>`.
5. **Continue the SLiC analysis** — `dls_calibration.py`, `plot_curvature.py`,
   `plot_histograms.py`, … as in Step 4 onwards below, on
   `filtered_puncta_A_values.txt`.

## Pipeline

```
Raw TIFFs                ┌─────────────────────┐
  (microscope)     ───►  │ 1. prepare_input.py  │  Split channels, reorder,
                         └────────┬────────────┘  crop, organize by voltage
                                  │
                                  ▼
             ┌────────────────────┴───────────────────────┐
             │ 2. spot detection (CMEanalysis algorithm)  │
             │  recommended: cme_detect (native_detection/)│
             │  fallback:    MATLAB CMEanalysis (external) │
             └────────────────────┬───────────────────────┘
                                  │
                                  ▼
             ┌────────────────────┴───────────────────────┐
             │ 3. analyze_cpp.py    (detection_cpp.tsv)   │  filter puncta,
             │    analyze_matlab.py (detection_v2.mat)    │  export A values
             └────────────────────┬───────────────────────┘
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
│           ├── cme_detect_output/  ← written by cme_detect
│           └── filtered_puncta_A_values.txt
├── figures/                    ← plots saved here (gitignored)
├── native_detection/           ← native C++ CMEanalysis-compatible detector
│   ├── CMakeLists.txt, src/, include/, tests/, third_party/
│   ├── README.md, PORTING_NOTES.md, THIRD_PARTY_NOTICES.md, LICENSE
│   └── build/                  ← created by cmake (gitignored)
├── tests/                      ← pytest suite
├── prepare_input.py
├── analyze_cpp.py              ← filter cme_detect output (native route)
├── analyze_matlab.py           ← filter MATLAB output (reference route)
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

Split multi-frame TIFFs into per-channel folders organized by laser
configuration. Folder names encode wavelength, transmissivity and PMT
voltage for each channel: e.g. `488nm_17.7pct_652V_561nm_15.3pct_876V/`.
Channel folders `ch1`, `ch2`, … follow wavelength-ascending order, so
`ch1` is always the lowest-wavelength channel.

#### How channels are resolved

The channel count comes from the **image data**, never from how many
lasers have nonzero transmissivity. `SizeC` is read from the Olympus
FV3000 metadata and cross-checked against the TIFF's actual page count;
a disagreement is a hard error naming both numbers.

Each channel is then resolved by **ID linkage**, not by position:

```
Channel CH<n> ID                 = <uuid>        ← the channel's UUID
Channel CH<n> linked laser index = <k>           ← 0-based into `laser name #<k+1>`
Detector linked channel ID #<m>  = <uuid>        ← find the m matching the channel
Detector voltage #<m>, Detector ID #<m>          ← that channel's PMT and detector
```

This matters because the detector block is **not** in the same order as
the channel block, and several `channel <field> #<i>` sub-blocks are
ragged (they describe different slices of the instrument configuration
and are not co-indexed with each other). Pairing any of these by index
silently mislabels channels. Detectors whose linked channel ID matches no
acquired channel are configured-but-unused and are ignored entirely.

Channels that are not fluorescence are dropped, each with a recorded
reason:

| Rule | Meaning |
|------|---------|
| Detector ID contains `LETD` | Transmitted-light/DIC channel |
| `channel deviceName` is `TD` | Transmitted-light channel (only applied when that sub-block has one entry per acquired channel, otherwise it cannot be trusted) |
| Laser data ID contains `Lambda` | Lambda-phase (spectral scan), not a main imaging phase |

Anything ambiguous or unparseable raises rather than guessing — a
wrong-but-silent channel assignment produces plausible-looking output and
quietly corrupts the analysis.

Dye names come from two sources. `dyeData excitationWavelength` matched
exactly to a channel's laser is unambiguous but incomplete (Texas Red
records 595 nm, so it never matches the 561 nm channel exciting it). The
`channel dyeName` block covers every channel but is ragged and
UUID-less, so it is used only when it has one entry per fluorescence
channel *and* every entry that can be cross-checked against an exact
excitation match agrees with it. One disagreement discards the whole
block. `dye_source` in the manifest records which route was used.

#### `channel_map.json`

Every output group folder gets a `channel_map.json` recording, for each
`ch<i>`: source page index, channel name, laser wavelength,
transmissivity, detector ID, PMT voltage, dye name and `dye_source`.
Excluded channels are listed with their reason. This makes the
`ch1`/`ch2` assignment auditable later without reopening the raw files
in Fiji.

The manifest is generated from the **same ordered slot list** used to
write the TIFFs, so `source_page` is always the page actually written to
that slot and the rest of the fields always belong to the channel that
page came from. It is never re-derived from wavelength.

`channel_order` states the real order and `frames_argument` records the
raw `--frames` value. Under an override that is not wavelength-ascending
it says so explicitly, and a warning is printed to stderr:

```json
"channel_order": "custom order via --frames 1,0 (NOT wavelength-ascending)",
"frames_argument": "1,0"
```

The **group folder name is built from the same slot list**, so it reads
in `ch1`, `ch2`, … order and `ch1`'s laser comes first. It therefore
tracks `--frames`:

```
(no --frames)   488nm_17.7pct_652V_561nm_15.3pct_876V/
--frames 1,0    561nm_15.3pct_876V_488nm_17.7pct_652V/
```

Excluded channels are not slots, so they never appear in the name.
`group_name_describes` in the manifest states this meaning. Note that a
remapped run lands in a *different* folder from a default run of the
same data, which is what keeps the two from silently merging.

`--frames` may also name **fewer** pages than there are resolved
channels, meaning "keep only these, in this order". The channels left
unnamed are dropped from the slots and the folder name, and are recorded
in `excluded_channels` with the reason `excluded by explicit --frames
selection` — so the manifest still accounts for every channel in the
file:

```
--frames 1      561nm_15.3pct_876V/     (CH1 recorded as deselected)
```

`--frames` may also **repeat** a page: `--frames 1,1` duplicates source
page 1 into both `ch1` and `ch2`, as a literal copy of the same image
pixels in both slots. The manifest stays truthful — both slots record
`source_page: 1` with the metadata of that one channel:

```
--frames 1,1    561nm_15.3pct_876V_561nm_15.3pct_876V/   (page 1 in ch1 AND ch2)
```

> **Warning.** Duplication is a compatibility/testing feature (e.g. to
> satisfy tooling that expects two channel folders). It must NOT be used
> when the second detector channel contains scientifically meaningful
> independent data — in that case keep the real slave-channel image and
> use the default `--hval-filter master` of `analyze_cpp.py` /
> `analyze_matlab.py` instead (see Step 3). A stderr warning is printed on
> every duplicated run.

If `--frames` forces a rule-excluded channel (e.g. the DIC page) into a
slot, that is allowed but recorded as `excluded_by_rule` on the slot and
warned about on stderr.

> **Note on channel roles.** The `ch1`/`ch2` ordering is by *wavelength*,
> which is not necessarily lipid-then-protein. In the bundled
> `data/olivia_data` acquisition `ch1` is 488/EGFP (protein) and `ch2` is
> 561/Texas Red (lipid) — the opposite of the `--lipid-col A_ch1` /
> `--protein-col A_ch2` defaults used downstream. Check
> `channel_map.json` and set those flags accordingly. Whichever folder
> holds the lipid channel is the one to put **first** in `cme_detect
> --channels` (and name in `--master`), or to select **first** (as master)
> in MATLAB `loadConditionData` — see Step 2.

```bash
python prepare_input.py \
    --input  data/march_3_experiment \
    --output data/march_3_experiment_matlab \
    --crop 1
```

Run with `--dry-run` first. It prints the resolved per-channel table and
every excluded channel with its reason, so the mapping is reviewable
before anything is written:

```
  dataHis eGFP SLiC post exess removal.tif: 2 channels -> 488nm_17.7pct_652V_561nm_15.3pct_876V/cell1/
      ch1  page 0  CH1  488nm  17.7%  652V  FV30-SD_D_1  EGFP
      ch2  page 1  CH2  561nm  15.3%  876V  FV30-SD_D_2
      excluded:
        CH3 (page 2, 561nm, FV31-LETD_D_1): transmitted-light/DIC channel: detector FV31-LETD_D_1 is an LETD detector
```

| Argument    | Meaning |
|-------------|---------|
| `--input`   | Folder containing raw .tif files |
| `--output`  | Where to save the split channels (inside `data/`) |
| `--frames`  | Comma-separated source page indices naming the pages to keep, in `ch1`, `ch2`, … order (e.g. `0,1`). May name **fewer** pages than there are resolved channels to keep a subset, and may **repeat** a page (`1,1` duplicates source page 1 into both `ch1` and `ch2` — compatibility/testing only, see warning above). Indices must be in range. Optional — resolution already yields the correct pages, so this is rarely needed |
| `--crop`    | Center crop divisor. `1` = no crop, `2` = center quarter |
| `--dry-run` | Parse metadata and print what would be created, without writing files |

### Step 2: Spot detection

Both routes run the same CMEanalysis detection algorithm; which physical
channel is the master/source channel is decided by you, from the metadata,
and is **not** necessarily `ch1`:

1. Check `channel_map.json` (and the recorded wavelength/voltage
   metadata) to see which physical acquisition channel each `ch1`/`ch2`
   folder holds.
2. Decide which physical channel is your lipid/**master** channel. Either
   folder may be the master.

#### Step 2, recommended: native `cme_detect`

Build once (Getting Started, step 6), then run with the master channel
first in `--channels` and named in `--master`:

```powershell
# ch1 holds the lipid (master) channel
& ".\native_detection\build\Release\cme_detect.exe" `
    --input "data\march_3_experiment_matlab\488nm_5.0pct_580V_561nm_3.2pct_500V" `
    --channels ch1,ch2 --master ch1 --seed 1 --threads 12

# ch2 holds the lipid (master) channel
& ".\native_detection\build\Release\cme_detect.exe" `
    --input "data\march_3_experiment_matlab\488nm_5.0pct_580V_561nm_3.2pct_500V" `
    --channels ch2,ch1 --master ch2 --seed 1 --threads 12
```

The run writes `<cell>/<master>/Detection/detection_cpp.tsv` for every
cell (only the master channel gets a `Detection/` folder) plus
`<condition>/cme_detect_output/`. Record the `--channels` order — Step 3
needs it verbatim. Full CLI, outputs and validation: see
[Native CMEanalysis-compatible detection](#native-cmeanalysis-compatible-detection-recommended)
and [native_detection/README.md](native_detection/README.md).

#### Step 2, reference/fallback: original MATLAB CMEanalysis workflow

See [PROTOCOL.pdf](PROTOCOL.pdf) for detailed CMEanalysis instructions. The
channel-selection order in MATLAB is load-bearing — it determines which
physical channel is the master/source channel and the column order of
everything CMEanalysis writes, and `analyze_matlab.py` must be given the
same order later:

3. In MATLAB `loadConditionData`, select the master channel **first**,
   then the slave channel(s).
4. For reproducibility, seed the RNG immediately before detection:

   ```matlab
   data = loadConditionData;
   rng(1);
   runDetection(data);
   ```

   The seed makes CMEanalysis's data-derived PSF sigma estimation
   reproducible; without it, repeated runs can yield slightly different
   sigma estimates, which propagate into the fitted amplitudes.
5. Record the selection order — Step 3 needs it verbatim.

The detection produces `detection_v2.mat` files inside the master
channel's `Detection/` subdirectory. Only the lipid/master channel (the
one selected **first**) gets this folder.

### Step 3: Analyze detection output

Filter puncta by intensity threshold and export amplitudes. Use the script
matching the Step 2 route — both apply the identical filter (the same
shared code) and write identical output files:

| Step 2 route | Script | Reads |
|---|---|---|
| native `cme_detect` (recommended) | `analyze_cpp.py` | `cell*/<master>/Detection/detection_cpp.tsv` |
| MATLAB CMEanalysis (fallback) | `analyze_matlab.py` | `cell*/<master>/Detection/detection_v2.mat` |

`--channels` means the **exact analysis order** — the list given to
`cme_detect --channels`, or the order the channels were selected in MATLAB
`loadConditionData` (the order stored in `frameInfo`) — *not*
numerical/sorted folder order. The first entry is the master/source
channel, so `--lipid-channel` must equal the first entry of `--channels`;
any other order is rejected with an error rather than silently
reinterpreted, because a mismatch would attach detection columns to the
wrong physical channel. `analyze_cpp.py` additionally compares
`--channels` with the order recorded in the TSV header and stops on a
mismatch.

```bash
# Native route: cme_detect was run with --channels ch1,ch2 --master ch1
python analyze_cpp.py \
    --input  data/march_3_experiment_matlab/488nm_5.0pct_580V_561nm_3.2pct_500V \
    --channels ch1,ch2 \
    --lipid-channel ch1

# Native route: cme_detect was run with --channels ch2,ch1 --master ch2
python analyze_cpp.py \
    --input  data/march_3_experiment_matlab/488nm_5.0pct_580V_561nm_3.2pct_500V \
    --channels ch2,ch1 \
    --lipid-channel ch2
```

```bash
# Example A: MATLAB selected ch1 first (master), ch2 second (slave)
python analyze_matlab.py \
    --input  data/march_3_experiment_matlab/488nm_5.0pct_580V_561nm_3.2pct_500V \
    --channels ch1,ch2 \
    --lipid-channel ch1 \
    --k-std 2.0 \
    --output-name filtered_puncta_A_values.txt
```

```bash
# Example B: MATLAB selected ch2 first (master), ch1 second (slave)
python analyze_matlab.py \
    --input  data/march_3_experiment_matlab/488nm_5.0pct_580V_561nm_3.2pct_500V \
    --channels ch2,ch1 \
    --lipid-channel ch2
```

Both scripts take the same arguments:

| Argument          | Meaning |
|-------------------|---------|
| `--input`         | Condition folder (voltage group from Step 1) |
| `--channels`      | Comma-separated channel folder names in the exact analysis order: as given to `cme_detect --channels`, or as selected in MATLAB `loadConditionData` / stored in `frameInfo`. The first entry is the master/source channel. Use a single name for lipid-only |
| `--lipid-channel` | The lipid/master channel (used for thresholding). Must equal the **first** entry of `--channels` |
| `--k-std`         | Threshold: keep if A > mean(c) + k·std(c). Default `2.0` |
| `--hval-filter`   | How the `hval_Ar` hypothesis test is used: `master` (default, require hval == 1 in the lipid/master channel only), `none` (hval not used for filtering), `all` (require hval == 1 in every requested channel) |
| `--output-name`   | Filtered output filename, saved inside `--input` folder |
| `--raw-output-name` | Unfiltered diagnostic filename (default: `raw_puncta_values.txt`) |
| `--detection-file` | `analyze_cpp.py` only: per-cell table name under `Detection/` (default `detection_cpp.tsv`; regression use: `detection_matlab.tsv`) |

#### hval filtering

The default `--hval-filter master` requires the hypothesis test to pass
only in the lipid/master channel — slave channels are deliberately **not**
required to pass. This matters whenever the master channel identifies a
real punctum whose signal in another channel is legitimately very weak:
requiring slave-channel `hval == 1` would censor exactly those low-signal
observations. Detection/localization comes from the master channel, and
the measured slave-channel amplitudes at those same locations are kept
regardless of their own hypothesis test:

```bash
# Default: master channel identifies puncta, weak slave signal retained
python analyze_matlab.py \
    --input  data/.../488nm_17.7pct_652V_561nm_15.3pct_876V \
    --channels ch1,ch2 \
    --lipid-channel ch1        # --hval-filter master is the default

# Ignore hval entirely / require it in every channel
python analyze_matlab.py --input data/... --channels ch1,ch2 \
    --lipid-channel ch1 --hval-filter none
python analyze_matlab.py --input data/... --channels ch1,ch2 \
    --lipid-channel ch1 --hval-filter all
```

The header comments of both output files state exactly which hval policy
was used. (`analyze_cpp.py` accepts `--hval-filter` identically.)

#### Raw diagnostic output

Every run (of either script) also writes `raw_puncta_values.txt` into the
condition folder — **every** punctum/candidate present in the detection
tables,
with `A`, `c` and `hval` for each channel and no Python-side filtering
whatsoever (no amplitude threshold, no hval filter; negative/zero A
values are retained). Columns for two channels:

```
source_image  A_ch1  c_ch1  hval_ch1  A_ch2  c_ch2  hval_ch2  passes_A_threshold_master  passes_hval_policy  passes_final_filter
```

The trailing `passes_*` columns are annotations only (would this row
pass the current filter settings?) — they never remove rows. The raw
file is written even when zero puncta pass the filter, so a too-strict
threshold can be diagnosed by inspecting it.

For **lipid-only** experiments (no protein channel):
```bash
python analyze_matlab.py \
    --input  data/.../488nm_5.0pct_580V_561nm_3.2pct_500V \
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
| `--lipid-col`   | Column name for lipid amplitude (default: `A_ch1`) |
| `--protein-col` | Column name for protein amplitude (default: `A_ch2`) |
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
| `prepare_input.py`        | Split and reorder TIFF channels, organize into the condition/cell/channel layout used by `cme_detect` and MATLAB |
| `native_detection/`       | Native C++ CMEanalysis-compatible detector `cme_detect` (recommended); own README, porting notes, tests, third-party notices |
| `analyze_cpp.py`          | Read `cme_detect` per-cell `detection_cpp.tsv` files, filter puncta (shared code with `analyze_matlab.py`), export TSV |
| `analyze_matlab.py`       | Read MATLAB detection `.mat` files, filter puncta, export TSV (reference/fallback route) |
| `dls_calibration.py`      | DLS-fluorescence distribution overlay to compute conversion factor |
| `plot_curvature.py`       | Convert amplitudes to diameters, plot protein density vs diameter |
| `plot_histograms.py`      | Plot amplitude histograms and estimated diameter distributions |
| `plot_overlay.py`         | Overlay normalized curvature-sorting curves across experiments |
| `plot_dls.py`             | Plot DLS distribution (number or intensity, raw or log) |
| `plot_dls_comparison.py`  | Side-by-side DLS vs fluorescence sqrt(A) per channel |
| `requirements.txt`        | Python dependencies: numpy, matplotlib, tifffile, h5py, openpyxl, pandas, scipy |

## Tests

Python (creates nothing outside `tmp`; the `cme_detect` end-to-end tests
are skipped automatically when the native backend is not built):

```powershell
pytest tests/test_analyze_matlab.py tests/test_analyze_cpp.py tests/test_prepare_channels.py
```

Native C++ unit/regression tests (266 checks against MATLAB/MEX reference
values), from a configured build:

```powershell
ctest --test-dir native_detection/build -C Release --output-on-failure
```

## Notes

- **Detection** uses the CMEanalysis algorithm
  ([DanuserLab/cmeAnalysis](https://github.com/DanuserLab/cmeAnalysis)),
  either through the native `cme_detect` port in `native_detection/`
  (recommended; writes `detection_cpp.tsv`) or through MATLAB
  (reference/fallback; writes `detection_v2.mat` with a `frameInfo` struct).
  In both cases the fields used downstream are `A` (amplitude), `c`
  (background), and `hval_Ar` (hypothesis test). The native port covers
  the default `loadConditionData` + `runDetection` path only, not the
  whole CMEanalysis package.
- **Only the master/lipid channel** — the channel passed *first* to
  `cme_detect --channels` (and named in `--master`), or selected *first*
  in MATLAB `loadConditionData`, which need not be `ch1` — gets a
  `Detection/` subfolder. Slave channel folders just contain the TIFF.
- **Licensing of `native_detection/`:** the native port re-implements
  GPL-3.0 algorithms from cmeAnalysis and links against GSL, so that
  subtree is GPL-3.0 (`native_detection/LICENSE`); attributions and the
  licences of the vendored GSL / libtiff builds are collected in
  `native_detection/THIRD_PARTY_NOTICES.md`. The rest of the repository
  carries no licence file at present.
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
- **Lipid-only experiments** are fully supported. Run `analyze_cpp.py` (or
  `analyze_matlab.py`) with `--channels ch1 --lipid-channel ch1`, then use `plot_histograms.py`
  (skip `plot_curvature.py` since it requires protein data).
- The `data/` and `figures/` folders are gitignored. Your microscopy data
  stays local.
