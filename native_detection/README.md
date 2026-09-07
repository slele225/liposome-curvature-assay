# native_detection — native C++ CMEanalysis-compatible detection (`cme_detect`)

`cme_detect` is a native C++ re-implementation of **exactly the cmeAnalysis
workflow this repository uses for spot detection**:

```matlab
data = loadConditionData;      % condition folder, channel folders in order, markers
rng(1);
runDetection(data);            % default options: data-driven PSF sigma, master/slave fitting
```

It lets the SLiC pipeline run without MATLAB, the Parallel Computing
Toolbox, or the cmeAnalysis MEX binaries. The original MATLAB workflow
remains available as the reference implementation (see the main
[README](../README.md) and [PROTOCOL.pdf](../PROTOCOL.pdf)).

## What was ported (and what was not)

The port covers the dependency closure of `loadConditionData` +
`runDetection(data)` **under default arguments**: condition/movie folder
discovery, TIFF loading, the data-driven PSF sigma estimation
(`getGaussianPSFsigmaFromData`, including the Gaussian-mixture fit and the
MATLAB `rng` stream it consumes), `pointSourceDetection` (LoG/prefilter,
candidate selection, `fitGaussian2D` MEX-equivalent Levenberg–Marquardt
fitting, amplitude and Anderson–Darling hypothesis tests), the master/slave
channel fitting in `runDetection`'s `main()`, and the `frameInfo` /
mask outputs.

**Not ported** (not on the default path): tracking, lifetime analysis, the
GUI, 3-D detection, mixture-model fitting, `runDetection` non-default
options other than `'Sigma'`. This is a port of the path SLiC uses, not of
the whole cmeAnalysis package. [PORTING_NOTES.md](PORTING_NOTES.md) records
the traced call graph, the classification of every dependency, and every
reconstruction decision (MEX internals, toolbox semantics, RNG).

## Build

### Windows (Visual Studio 2022, tested with VS 2022 Build Tools + CMake 3.30)

Pre-built x64/MSVC dependencies are vendored in `third_party/` (GSL 2.8 as
DLL + import library, libtiff 4.7.0 static; see
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md)). Nothing else needs to be
installed beyond the Visual Studio C++ tools and CMake (the CMake bundled
with Visual Studio is sufficient).

```powershell
cd "<repo>\native_detection"

cmake -S . -B build `
  -G "Visual Studio 17 2022" `
  -A x64

cmake --build build --config Release
```

The executables are produced in `build\Release\`:

```text
build\Release\cme_detect.exe      the detector
build\Release\cme_tests.exe       unit/regression tests vs MATLAB + MEX reference values
build\Release\fit_probe.exe       diagnostic: fitter vs recorded MEX results
build\Release\profile_sigma.exe   diagnostic: timing of the sigma-estimation pass
build\Release\gsl.dll, gslcblas.dll   copied next to the executables by the build
```

`build/` is gitignored; only sources, `CMakeLists.txt`, tests, reference
data and the vendored libraries are tracked.

### Other platforms / system libraries

Configure with `-DCME_USE_SYSTEM_LIBS=ON` and provide GSL and libtiff
through `find_package(GSL)` / `find_package(TIFF)`. The vendored binaries
are only used on Windows when `third_party/gsl/lib/gsl.lib` and
`third_party/tiff/lib/tiff.lib` exist. OpenMP is optional
(`-DCME_USE_OPENMP=OFF` to disable). Non-Windows builds have not been part
of the regression validation.

## Run

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
  --threads N          OpenMP threads (default 1); output does not depend on the thread count
  --par-level L        where the threads are used: 'candidate' (default: over candidate fits and
                       image columns within each frame) or 'movie' (over movies/frames/images)
  --profile            print a per-phase timing breakdown at the end (diagnostic)
```

(`cme_detect --help` prints the same text; that is the authoritative version.)

### Channel order and the master channel

* `--channels` is the **exact analysis order**, i.e. the order in which
  `loadConditionData` would have received the channel folders. The first
  entry is the master/source channel; every other entry is a slave channel
  fitted at the master's positions.
* `--master` must equal the first entry of `--channels`. The tool refuses
  inconsistent orderings instead of reordering silently.
* Channel folder names carry **no** master/slave meaning. `ch1`, `ch2`, …
  are wavelength-ascending slots written by `prepare_input.py`; which one
  holds the lipid (master) channel comes from `channel_map.json` and the
  recorded wavelength/voltage metadata. Both of these are valid:

  ```text
  --channels ch1,ch2 --master ch1      ch1 = 488 nm master, ch2 = 561 nm slave
  --channels ch2,ch1 --master ch2      ch2 = 561 nm master, ch1 = 488 nm slave
  ```
* `--seed N` corresponds to MATLAB `rng(N)` immediately before
  `runDetection`. It only influences the k-means++ initialisation of the
  Gaussian-mixture fit inside the PSF-sigma estimation; `--seed 1` is the
  reproducible equivalent of the project's `rng(1)` workflow.
* `--sigma s1,s2` bypasses the data-driven sigma estimation, like
  `runDetection(data, 'Sigma', ...)`.
* `--markers` is accepted for completeness only; markers do not affect the
  default detection path.

### Expected input layout (as for `loadConditionData`)

```text
condition/
    cell1/ch1/*.tif   cell1/ch2/*.tif
    cell2/ch1/*.tif   cell2/ch2/*.tif
```

Movie folders are those whose name contains `cell` (case-insensitive,
`--movie-selector`) and are sorted by the number following `cell`. A
channel folder may hold one multi-page TIFF or one TIFF per frame. This is
exactly what `prepare_input.py` produces.

## Outputs

* `<condition>/cme_detect_output/detections_all.tsv` — one row per detection
  with every `frameInfo` field of `detection_v2.mat` (per-channel columns
  suffixed by channel name, e.g. `A_ch1`, `c_ch1`, `hval_Ar_ch1`;
  coordinates are 1-based MATLAB pixel coordinates), plus `movie`, `frame`,
  `index`, `source_image`.
* `<condition>/cme_detect_output/summary.tsv` — `source_image, movie, frame,
  index, A_<ch>, c_<ch>, hval_<ch>, x_<ch>, y_<ch>`.
* `<condition>/cme_detect_output/sigma.tsv` — the PSF sigma per channel.
* `<cell>/<master>/Detection/detection_cpp.tsv` — the same columns as
  `detections_all.tsv`, one file per movie under the master channel, i.e.
  the MATLAB layout (`detection_v2.mat` lives in the same place). This is
  what [`analyze_cpp.py`](../analyze_cpp.py) reads. Disabled by
  `--no-matlab-layout`.
* `<cell>/<master>/Detection/dmasks.tif` (or `Detection/Masks/dmask_NNN.tif`
  for one-TIFF-per-frame movies) — detection masks, as MATLAB writes them.
  Disabled by `--no-masks`.

`detection_v2.mat` (HDF5) itself is not written; the TSV carries the same
columns. The per-channel column order in every table is the `--channels`
order, and `analyze_cpp.py` checks it against its own `--channels` argument.

## Downstream: `analyze_cpp.py`

[`analyze_cpp.py`](../analyze_cpp.py) (repository root) is the counterpart of
`analyze_matlab.py` for the TSV output. It imports the filter and table
writers from `analyze_matlab.py`, so both backends apply the identical
filter (`A_master > mean(c_master) + k·std(c_master)` and the
`--hval-filter` policy) and write identical `raw_puncta_values.txt` /
`filtered_puncta_A_values.txt` files:

```powershell
python analyze_cpp.py `
    --input "C:\path\to\prepared_condition" `
    --channels ch1,ch2 `
    --lipid-channel ch1
```

## Threading and determinism

* `--threads N` enables OpenMP parallelism. With the default
  `--par-level candidate` the threads work inside each frame: the
  separable convolutions, the local-maximum search and the per-pixel
  prefilter are split over image columns/pixels, and the Gaussian fits are
  distributed over candidates; the movie/frame loops run sequentially.
  `--par-level movie` instead distributes movies, frames and
  sigma-estimation images over the threads with sequential inner work
  (the original scheme). In both schemes every unit of work writes only
  its own output slots, no floating-point reduction is split across
  threads and the RNG-consuming GMM step is sequential, so **results do
  not depend on the thread count or the level**: in the regression runs
  `detections_all.tsv` was byte-identical for 1, 2, 4, 8 and 12 threads at
  both levels.
* The two levels are within a few percent of each other when there are
  more movies than threads; the candidate level also scales when there are
  only a few movies (or a single multi-frame movie) and needs less memory
  (one frame's intermediates at a time).
* With the same inputs, `--seed` and `--threads`, repeated runs produce
  identical output. The MATLAB RNG stream (`rng(N)`, Mersenne Twister with
  MATLAB's `rand`/`randi`/`randperm` semantics) is reproduced so that the
  GMM initialisation matches MATLAB's for the same seed.

## Regression validation against MATLAB cmeAnalysis

The backend was regression-validated against the untouched MATLAB
cmeAnalysis implementation (R2025a). Full details, per-stage tolerances and
the exact numbers are in [PORTING_NOTES.md §7](PORTING_NOTES.md).

* **Unit level** (`tests/unit_tests.cpp`, 266 checks): MATLAB RNG stream,
  `norminv`/`tcdf`/`normcdf`, Anderson–Darling test and the MEX's 940
  recorded hAD decisions, `conv2`, morphology, `fitGaussian2D` in every mode
  used (`xyAc`, `xyasc`, `Ac`, `xyac`, NaN-masked windows) against the
  original MEX, `pointSourceDetection`, `gmdistribution.fit`.
* **10-cell reference data set** (1024×1024, ch1 master / ch2 slave,
  `rng(1)`): identical detection counts in every cell (48,325 rows in the
  same order), identical `hval_Ar`, `hval_AD`, `isPSF`, `x_init`/`y_init`
  and mask decisions; master-channel x, y and A within numerical tolerance
  for every row; the downstream filter kept exactly the same 45,362 puncta
  as `analyze_matlab.py` on MATLAB's `detection_v2.mat`.
* **A second, practical SLiC data set** (488 nm master): MATLAB and
  `cme_detect` both produced 25,395 raw detections with the same cell/row
  ordering and the same `hval` decisions; the downstream filter retained the
  exact same 9,960 puncta from both, with zero inclusion/exclusion
  mismatches; master-channel A differed by at most ~3×10⁻⁵ absolute
  (relative errors around 10⁻⁷).
* Residual differences are limited to a small number of numerically
  unstable fits (typically slave-channel localized fits with huge position
  uncertainties, i.e. windows where the spot is absent in the slave
  channel) and to the last digits of near-zero background values. No
  scientifically meaningful downstream differences were observed. This is
  **not** a claim of bit-for-bit floating-point equivalence.

### Reproducing the comparison

1. `reference_tools/dump_reference.m` runs the MATLAB pipeline with the
   original cmeAnalysis functions and dumps all intermediates in the layout
   of `cme_detect --dump-dir`:
   ```matlab
   addpath(genpath('C:\path\to\cmeAnalysis-master\software'));
   addpath('native_detection/reference_tools');
   dump_reference('C:\data\condition', {'ch1','ch2'}, 'C:\ref_dump', 'NumImageMovies', 2);
   ```
2. `cme_detect ... --dump-dir C:\cpp_dump`
3. `python tests/compare_reference.py --cpp C:\cpp_dump --matlab C:\ref_dump`
4. End to end against the untouched `runDetection`:
   `reference_tools/run_matlab_reference.m` (also prints MATLAB timing,
   exports `detection_matlab.tsv`), then
   `python tests/compare_reference.py --e2e-cpp <Detection/detection_cpp.tsv> --e2e-matlab <Detection/detection_matlab.tsv>`,
   and `python analyze_cpp.py ... --detection-file detection_matlab.tsv`
   to run the SLiC filter on the MATLAB export.

The `.m` files in `reference_tools/` add the cmeAnalysis `software/` folder
to the path via a relative `addpath` that assumed the port lived inside a
cmeAnalysis checkout; from this repository, add cmeAnalysis to the MATLAB
path yourself first (as in step 1). `reference_tools/psd_dump.m` is a
modified copy of cmeAnalysis's `pointSourceDetection.m` (GPL-3.0, Francois
Aguet / Danuser Lab) with debug hooks; it exists only to produce golden
intermediates.

## Known numerical limitations

* MATLAB's `conv2` accumulation order could not be identified, so the
  prefilter / LoG / `A_est` / `c_est` images agree to ~10⁻¹⁰ relative
  (absolute differences ≤ 3×10⁻⁵ on values of order 10⁵). This floor
  propagates into the last digits of near-zero background (`c`) values.
* Levenberg–Marquardt fits whose iteration is numerically unstable (free-σ
  refits that do not converge within 500 iterations during sigma
  estimation; slave-channel fits with position uncertainties of tens of
  pixels) can land at slightly different points of the same solution
  basin. On the reference data this affected ~0.3 % of slave-channel rows
  and no master-channel A, x, y values beyond tolerance.
* The sigma estimate therefore agrees with MATLAB to ~10⁻⁹ relative rather
  than exactly; with identical `svect` the EM is bit-identical.

## Tests

Native unit/regression tests (from a configured build):

```powershell
ctest --test-dir build -C Release --output-on-failure
```

or directly (the argument locates `tests/reference_values.txt`, and
`tests/ad_reference*.txt` next to it):

```powershell
.\build\Release\cme_tests.exe .\tests\reference_values.txt
```

The Python side (`analyze_cpp.py`, and an end-to-end smoke test that runs
`cme_detect.exe` on synthetic images when the build exists) is covered by
the repository's pytest suite:

```powershell
pytest tests/test_analyze_cpp.py
```

`fit_bench` (built alongside) is an in-process micro-benchmark of
`fitGaussian2D` on recorded real candidate windows; with a results file it
doubles as a bit-identity check between two builds:

```powershell
.\build\Release\fit_bench.exe C:\path\to\mex_probe.txt 20 results.txt
```

## Performance

In project regression tests the native implementation reproduced
cmeAnalysis's detection/filter decisions while providing substantial
runtime reductions. Runtime depends on image count and size, CPU and
thread count. The numbers below were measured on one laptop (Intel Core
i7-1355U: 2 performance + 8 efficiency cores, 12 hardware threads, Windows
11); the "optimised" column is the current build, "initial port" the
first validated build (September 2026), both producing byte-identical
output.

| data set | run | initial port | optimised |
|---|---|---|---|
| 10-cell reference set (1024×1024, 2 channels; 80 sigma-estimation image instances of 20 distinct images) | MATLAB R2025a `loadConditionData` + `rng(1)` + `runDetection` (2-worker `parfor` pool) | 1903 s | — |
| | `cme_detect --threads 1` | 2451 s | 360.4 s |
| | `cme_detect --threads 12` | 366 s | 59.9 s |
| 21-cell SLiC condition (512×512, 2 channels; 42 instances of 21 distinct images) | `cme_detect --threads 1` | 702.7 s | 368.8 s |
| | `--threads 2` | 389.3 s | 171.2 s |
| | `--threads 4` | 190.8 s | 86.1 s |
| | `--threads 8` | 118.8 s | 56.0 s |
| | `--threads 12` | 102.1 s | 52.9 s |

Where the time goes (single thread, 21-cell set, initial port): Gaussian
fitting 95 % (Levenberg–Marquardt iterations 88 %, of which the pivoted QR
of the Jacobian is about 60 %), per-pixel prefilter t-test 2 %,
convolutions 0.6 %, everything else (TIFF I/O, local maxima, connected
components, output) below 0.5 % each. The optimisations that produced the
"optimised" column, all with bit-identical output, are documented in
[PORTING_NOTES.md §8](PORTING_NOTES.md); the largest single one is that the
sigma estimation no longer re-processes the identical image instances that
MATLAB's `round(linspace(1, L, nf))` frame sampling produces for short
movies. Scaling flattens beyond 8 threads on this CPU because only the two
performance cores are hyper-threaded.

## Layout

```text
native_detection/
  CMakeLists.txt, README.md, PORTING_NOTES.md, THIRD_PARTY_NOTICES.md, LICENSE
  include/cme/*.hpp, src/*.cpp        the port (PORTING_NOTES.md §1 maps them to the .m files)
  third_party/gsl, third_party/tiff   pre-built GSL 2.8 / libtiff 4.7.0 (MSVC x64) + licences
  third_party/gsl116                  GSL 1.16 Levenberg-Marquardt sources (verbatim, GPL-3) + g116_blas.h inline BLAS
  include/cme/profile.hpp, src/profile.cpp   --profile phase timers
  tests/unit_tests.cpp                regression tests vs MATLAB/MEX reference values (+ thread/screening identity checks)
  tests/reference_values.txt, tests/ad_reference*.txt   reference data produced by MATLAB
  tests/compare_reference.py          stage-by-stage / end-to-end comparison tool
  tests/fit_bench.cpp                 fitGaussian2D micro-benchmark / bit-identity check between builds
  tests/fit_probe.cpp, tests/profile_sigma.cpp          diagnostics
  reference_tools/*.m                 MATLAB scripts that generate the golden data
  ../analyze_cpp.py                   analyze_matlab.py equivalent for the TSV output (repo root)
```

## Relationship to upstream cmeAnalysis and licensing

`cme_detect` re-implements algorithms from
[cmeAnalysis](https://github.com/DanuserLab/cmeAnalysis) (Francois Aguet,
Danuser Lab, UT Southwestern), which is distributed under the GNU General
Public License v3 (`GPL-License.txt` in the upstream repository). No MATLAB
toolbox source is copied; toolbox functions were only read to determine the
semantics the port must reproduce. The port is therefore licensed under
**GPL-3.0** as well ([LICENSE](LICENSE)), and it links against GSL
(GPL-3.0) and libtiff (libtiff licence). Attribution and licence texts for
every vendored component are listed in
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
