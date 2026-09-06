# Porting notes: `loadConditionData` + `rng(1); runDetection(data)` to native C++

> **Location note.** This port was developed as `cpp_port/` inside a
> cmeAnalysis checkout and now lives in the SLiC repository as
> `native_detection/`. Path references below that say `cpp_port/` refer to
> this directory; `analyze_cpp.py` moved to the repository root.

This document records the traced dependency closure of the MATLAB workflow

```matlab
data = loadConditionData;
rng(1);
runDetection(data);
```

under **default arguments**, classifies every dependency, and documents every
place where the C++ port had to reconstruct behaviour that is not available as
MATLAB source (compiled MEX files, toolbox built-ins).

Reference MATLAB installation used for tracing and golden data: **R2025a Update 1**
(Statistics and Machine Learning Toolbox 25.1), Windows 11.

Licensing: cmeAnalysis is GPL-3.0 (`GPL-License.txt`). The port re-implements
GPL algorithms and is therefore GPL-3.0 as well (see `native_detection/README.md`).
No MATLAB toolbox source is copied; toolbox `.m` files were only read to
determine the exact semantics that the port must reproduce.

---

## 1. Traced call graph (default path only)

Legend of classes:

* **[M]** MATLAB `.m` implementation available in `software/`
* **[B]** MATLAB built-in / toolbox function (no source or toolbox source)
* **[X]** compiled MEX binary (`software/mex/*.mexw64`, **no C source in repo**)
* **[–]** not exercised on the default path (not ported)

### 1.1 `loadConditionData` (interactive, replaced by CLI arguments)

| Callee | Class | Ported? | Notes |
|---|---|---|---|
| `inputParser` | [B] | n/a | defaults: `Parameters=[1.49 108 6.5e-6]`, `MovieSelector='cell'`, `StrictSelector=false`, `IgnoreEmptyFolders=false`, `FrameRate=[]` |
| `uigetdir`, `input` | [B] | CLI | condition dir, channel names and markers come from CLI flags |
| `recursiveDir(condDir, 2)` | [M] | yes | lists `condDir`, its visible sub-dirs and sub-sub-dirs (depth ≤ 2), all with trailing separator; hidden dirs (leading `.`) skipped |
| `getDirFromPath` | [M] | yes | last directory name + parent |
| `regexpi(cellDirs,'cell','once')` | [B] | yes | any dir whose *name* contains `cell` (case-insensitive) is a movie; `StrictSelector=false` so position does not matter |
| `unique(cellPar)` + `sortStringsByToken(..., 'cell', 'post')` | [M]/[B] | yes | within each parent, movie dirs are sorted by the integer that directly follows the token `cell` (`(?<=cell)\d+`); dirs with no number are **dropped**; `sort` is stable |
| `getFluorPropStruct` | [M] | validation only | markers are only validated; they never influence the default (`SigmaSource='data'`) detection path |
| date / framerate regexes | [B] | yes (metadata only) | `\d{6}+` date, `_(\d+)?(\.)?\d+s`, `_\d+ms`, default 2 s. Not used by detection. |
| `dir`, `regexpi(tmp,'\.tif\|\.stk')` | [B] | yes | files (not dirs) whose name does not start with `.` and contains `.tif` or `.stk` anywhere (case-insensitive); MATLAB `dir` order = OS lexicographic (Windows: case-insensitive) |
| numeric re-sort `regexp(tmp,'\d+(?=\.)')` | [B] | yes | if the digit runs immediately before a `.` do not all have the same length, files are re-sorted numerically by that number |
| `imfinfo` | [B] | yes (libtiff) | `imagesize = [Height Width]`, `movieLength = numel(info)` for a single multi-page file, or number of files otherwise |
| `data.source = channels{1}` | | yes | **master = first channel given**. Physical folder names are irrelevant. |
| `maskPaths` | | yes | single file per channel → `Detection/dmasks.tif`; file list → `Detection/Masks/dmask_%0Nd.tif` |
| detection_v2.mat reload branch | [–] | no | only taken when frames are missing |

### 1.2 `runDetection(data)` defaults

`Sigma=[]`, `SigmaSource='data'`, `RemoveRedundant=true`, `Overwrite=false`,
`Master=[]` → `find(strcmpi(channels, source))` = 1, `Alpha=0.05`, `CellMask=[]`.

| Callee | Class | Ported? | Notes |
|---|---|---|---|
| `exist(.../Detection/detection_v2.mat)` skip logic | [B] | CLI `--overwrite` | port always recomputes unless told otherwise |
| `getGaussianPSFsigma` | [M] | **no** [–] | only for `SigmaSource='model'` |
| frame sampling for sigma | | yes | `nf = round(40/nd)`; `fidx = round(linspace(1, movieLength, nf))` per data set (**for single-frame movies the one frame is repeated `nf` times**); frames collected into an `nd × nf` cell and flattened **column-major** (`vertcat(frames(:))`): order is `f=1: cells 1..nd, f=2: cells 1..nd, ...` |
| `imread` / `readtiff` | [B]/[M] | yes (libtiff) | 16-bit grayscale TIFF → `double` |
| `getGaussianPSFsigmaFromData(frames, 'Display', true)` | [M] | yes | see §1.3; `Display` only plots |
| `sigma(sigma<1.1) = 1.1` | | yes | |
| `main(...)` per data set | [M] | yes | see §1.4 |
| `getShortPath` | [M] | cosmetic | printing only |
| `parfor` | [B] | OpenMP later | results are indexed by frame; order independent |

### 1.3 `getGaussianPSFsigmaFromData` (default path)

| Callee | Class | Ported? | Notes |
|---|---|---|---|
| `pointSourceDetection(img, 1.5, 'Mode', 'xyac')` | [M] | yes | fixed σ=1.5 first pass (all other options default) |
| `fitGaussians2D(img, x, y, A, 1.5, c, 'xyasc')` | [M] | yes | refit with **free σ** at the detected positions, no mask, `Alpha=AlphaT=0.05` |
| `isPSF = ~hval_AD & pval_Ar < 0.05` ; `svect = s(~isnan(s) & isPSF)` | | yes | |
| `gmdistribution.fit(svect', n, 'Options', statset('maxIter',200))` for n = 1,2,3 | [B] stats toolbox | **reimplemented** (§3) | R2025a defaults: `Start='plus'` (k-means++), `CovType='full'`, `SharedCov=false`, `Regularize=0`, `TolFun=1e-6`, `ProbabilityTolerance=1e-8`, `Replicates=1` |
| `min(BIC)` over the 3 fits | | yes | `BIC = 2*NlogL + (3k-1)*log(n)` for 1-D full covariance |
| `sort(mu)`; `svec = sqrt(Sigma)`; `amp = PComponents`; `argmax(amp./(sqrt(2*pi)*svec))` → `sigma = mu(idx)` | | yes | component with the highest **peak density** wins |
| `try/catch` → `sigma = mean(svect)` | | yes | any error inside (ill-conditioned covariance, zero variance, too few samples `n<=k`) falls back to the plain mean |
| `randsample`/`datasample` RNG use | [B] | reimplemented (§3.2) | this is the only place `rng(1)` matters |
| plotting (`setupFigure`, `histc`, `bar`, `pdf`, ...) | [–] | no | display only; `prctile` etc. never affect `sigma` |

### 1.4 `main()` inside `runDetection` (per data set, per frame)

| Callee | Class | Ported? | Notes |
|---|---|---|---|
| `pointSourceDetection(img, sigma(mCh), 'Alpha', 0.05, 'Mask', [], 'RemoveRedundant', true)` | [M] | yes | mode `'xyAc'`, `FitMixtures=false` |
| `pstruct.s = sigma` (all channels), `rmfield s_pstd` | | yes | |
| `dRange{c} = [min max]` of every channel image | | yes | |
| `bwconncomp(mask)` / `labelmatrix` | [B] IPT | yes | default **8-connectivity**; `maskN` = component pixel count, `maskA` = mean image intensity over the component containing `(y_init, x_init)` |
| slave fit 1: `fitGaussians2D(img, x_m, y_m, [], sigma_s, [], 'Ac')` | [M] | yes | fixed (x,y) from master, σ fixed, **A and c free**; `A=[]`,`c=[]` → `c_init` = mean of annulus `ceil(3σ)≤r≤ceil(4σ)`, `A_init = max(window)-c_init` |
| slave fit 2: `fitGaussians2D(img, x_m, y_m, A_1, sigma_s, c_1, 'xyAc')` | [M] | yes | localized refit initialised from fit 1 |
| choice: `dist(master, loc) < 3*sigma(mCh) & A_loc > A_fixed` → use fit 2 else fit 1 | | yes | fields copied: `x y A c x_pstd y_pstd A_pstd c_pstd sigma_r SE_sigma_r RSS pval_Ar hval_Ar hval_AD` |
| `nanIdx = isnan(pstructSlave.x)` → drop detection from **all** channels | | yes | slave-border / failed fixed fit removes the detection |
| `isPSF(ci,:) = ~hval_AD(ci,:)` | | yes | |
| `xCoord, yCoord, amp` | | yes | tracker fields: `[x x_pstd]`, `[y y_pstd]`, `[A A_pstd]` of master |
| `imwrite(uint8(255*mask), ..., 'lzw')` | [B] | yes (libtiff LZW) | |
| `save(... 'detection_v2.mat', '-v7.3')` | [B] | TSV instead | see §6 |

### 1.5 `pointSourceDetection` (defaults + `Mode`)

Defaults: `Alpha=0.05`, `Mask=[]`, `FitMixtures=false`, `RemoveRedundant=true`,
`RedundancyRadius=0.25`, `Prefilter=true`, `RefineMaskLoG=true`,
`RefineMaskValid=true`, `ConfRadius=[]`, `WindowSize=[]`.

| Step | Class | Notes |
|---|---|---|
| `w = ceil(4σ)`, `g = exp(-x²/(2σ²))`, `x=-w:w`, `u = ones` | | |
| `padarrayXT(img,[w w],'symmetric')` | [M] | mirror padding **without** edge repetition of the border pixel? No: MATLAB `symmetric` = `[1:M M:-1:1]` style. `padarrayXT` is a copy of IPT `padarray` **except** its `SymmetricPad` uses `dimNums = [1:M M-1:-1:2]` (period `2M-2`), i.e. the border pixel is **not** duplicated (whole-sample symmetric, "XT"). Ported literally. |
| `conv2(g',g,·,'valid')`, `conv2(u',u,·)`, `conv2(u',u,·.^2)` | [B] | separable full convolution then crop: output same size as `img` because pad = kernel half width |
| `imgLoG = (2 fg/σ² − (conv2(g,gx2,·)+conv2(gx2,g,·))/σ⁴) / (2πσ²)` with `gx2 = g.*x.^2` | [B] | `conv2(hcol,hrow,A)`: columns (y) with `hcol`, rows (x) with `hrow`; all kernels symmetric so orientation does not matter numerically except for summation order |
| `A_est`, `c_est` closed-form LSQ | | `A_est = (fg − gsum·fu/n)/(g2sum − gsum²/n)`, `c_est = (fu − A_est·gsum)/n`, `n=(2w+1)²` |
| Prefilter t-test | [B] `tcdf`, `norminv`, `inv` | `RSS = A²·g2sum − 2A(fg − c·gsum) + (fu2 − 2c·fu + n c²)`, clamp `<0 → 0`; `σ_e² = RSS/(n−3)`; `σ_A = sqrt(σ_e²·C(1,1))`, `C = inv(J'J)`, `J=[g(:) 1]`; `σ_res = sqrt(RSS/(n−1))`; `k = norminv(1−α/2)`; `SE_σc = σ_res/sqrt(2(n−1))·k`; Welch df `df2 = (n−1)(σ_A²+SE²)²/(σ_A⁴+SE⁴)`; `scomb = sqrt((σ_A²+SE²)/n)`; `T = (A − σ_res·k)/scomb`; `pval = tcdf(−T, df2)`; `mask = pval < 0.05` (hard-coded 0.05, not `alpha`) |
| `locmax2d(imgLoG, 2*ceil(σ)+1)` | [M] uses `ordfilt2` [B] | strict local maximum in odd square window (`fImg2==fImg` removes plateaus, i.e. max value must be unique in the window); border strip of half-window width set to 0 |
| `imgLM = allMax .* mask` | | |
| `RefineMaskLoG`: `logThreshold = min(imgLoG(imgLM~=0))`; `mask = mask \| (imgLoG >= logThreshold)`; re-select `imgLM = allMax .* mask` | | |
| `[lmy,lmx] = find(imgLM~=0)` | [B] | **column-major order** (x-major then y) — order of candidates is preserved into the output rows |
| `fitGaussians2D(img, lmx, lmy, A_est(lmIdx), σ, c_est(lmIdx), 'xyAc', 'mask', mask, 'alpha', α, 'ConfRadius', [], 'WindowSize', [])` | [M] | see §1.6 |
| remove `isnan(x)` | | |
| `idx = hval_Ar == 1` | | |
| `KDTreeBallQuery(pM, pM, 0.25)` | [X] | for every point, indices of all points (including itself) within 0.25 px. Replaced by an exact brute-force O(N²)/grid search — the result set is identical for `dist <= r` (boundary ties have probability 0 for continuous fitted coordinates) |
| keep only min-RSS point of each redundant group | | `idx(group(RSS ~= min(RSS))) = 0` applied per query point |
| `hval_Ar`, `hval_AD` → logical; `isPSF = ~hval_AD` | | |
| `RefineMaskValid`: keep only 8-connected components of `mask` containing at least one `(y_init,x_init)` | [B] `bwconncomp`,`labelmatrix` | |
| `fitGaussianMixtures2D`, `getMultiplicity` | [–] | mixture fitting disabled by default, **not ported** |

### 1.6 `fitGaussians2D` (wrapper around the MEX)

| Step | Class | Notes |
|---|---|---|
| `bwlabel(mask)` (8-conn) or zeros | [B] | only the *partition* matters (label numbers are never compared across calls) |
| `xi = round(x)`, `yi = round(y)` | [B] | MATLAB `round` = half away from zero (`std::round`) |
| `kLevel = norminv(1−Alpha/2)` | [B] | 1.959963984540054 |
| `iRange = [min max]` of the whole image | | |
| `estIdx = regexpi('xyAsc', ['[' mode ']'])` | [B] | indices of estimated params in canonical order x,y,A,s,c — **case-insensitive**, so `'xyac' ≡ 'xyAc'` |
| `sigma_max = max(sigma)`; `w2 = ceil(2σmax)`; `w4 = ceil(4σmax)` | | |
| annulus for `c=[]`: `ceil(3σ) ≤ r ≤ ceil(4σ)` on `meshgrid(-w4:w4)` | | |
| `g = exp(-(-w4:w4).²/(2σmax²))`, `g = g'*g` (used only for `mask_Ar`) | | |
| border test: `xi>w4 && xi<=nx−w4 && yi>w4 && yi<=ny−w4` (1-based) else all outputs stay NaN | | |
| `maskWindow`: labels in window, own label (value at window centre) zeroed; other components → `NaN` in `window` | | |
| `npx = #finite`; require `npx >= 10` | | |
| `fitGaussian2D(window, [x−xi, y−yi, A_init, σ, c_init], mode)` | [X] | §2 |
| accept iff `−w2 < dx < w2 && −w2 < dy < w2 && A < 2*diff(iRange)` | | |
| `x = xi+dx`, ...; `stdVect(estIdx) = prmStd` | | |
| `sigma_r = res.std`, `RSS = res.RSS`, `SE_sigma_r = res.std/sqrt(2(npx−1))`, `hval_AD = res.hAD` | | |
| `df2 = (npx−1)(σ_A²+SE²)²/(σ_A⁴+SE⁴)` with `SE = SE_sigma_r·k`; `scomb = sqrt((σ_A²+SE²)/npx)`; `T = (A − res.std·k)/scomb` | | for fixed-A modes `σ_A = 0` |
| `mask_Ar = sum(A·g > res.std·k)` | | |
| `pval_Ar = tcdf(−T, df2)` for **all** points (unfitted points have `T=0, df2=0` → `NaN`) ; `hval_Ar = pval_Ar < AlphaT` | [B] | |

### 1.7 Not on the path (not ported)

Tracking (`runTracking`, u-track), lifetime analysis, GUI, 3-D fitting
(`fitGaussian3D`), anisotropic fitting (`fitAnisoGaussian2D`), mixture fitting
(`fitGaussianMixture2D`, `fitGaussianMixtures2D`), `steerableDetector`,
`vectorialPSF`, `getGaussianPSFsigma` (model-based σ), `binterp`, `conv3fast`,
`createDistanceMatrix`, `mexLap`, `KDTreeRangeQuery`, plotting/`setupFigure`.

---

## 2. `fitGaussian2D` MEX reconstruction

**No C/C++ source for any MEX exists in this repository** (verified with a full
file search for `*.c`, `*.cpp`, `*.h`). Only the `.m` help stubs exist. The
reconstruction below is based on: the `.m` help stub (`software/mex/fitGaussian2D.m`),
the way `fitGaussians2D.m` consumes the outputs, `adtest1.m` (Aguet's MATLAB
version of the same Anderson–Darling test), and knowledge of the original
Aguet `fitGaussian2D.c` (GPL, distributed with cmeAnalysis/u-track source
releases) whose documented design is: **GSL Levenberg–Marquardt,
`gsl_multifit_fdfsolver_lmsder`, analytic Jacobian, `gsl_multifit_test_delta`
stopping rule, `gsl_multifit_covar` covariance**. The help stub itself says
"options: vector [maxIter eAbs eRel] ... See GSL documentation", which
confirms the GSL `lmsder`/`test_delta` design. GSL is therefore used in the port.

### 2.1 Model

Window `data` is an odd square `nx × nx` (`nx = 2*w4+1`), origin at the centre pixel
(`b = nx/2` integer division), **x runs along columns, y along rows** (MATLAB
`meshgrid` convention, column-major storage):

```
f(x, y) = A · exp(−((x − xp)² + (y − yp)²) / (2 s²)) + c
```

Parameter vector order: `[xp, yp, A, s, c]`. Residual `r_i = f(x_i,y_i) − data_i`
over the non-NaN pixels only (NaN pixels are masked). Mode string chooses the
free parameters (any subset of `x y a s c`, lower-cased before matching);
the others stay fixed at their initial values.

Jacobian (analytic):

```
∂r/∂xp = (x−xp)·A·g/s²      ∂r/∂yp = (y−yp)·A·g/s²
∂r/∂A  = g                  ∂r/∂s  = ((x−xp)²+(y−yp)²)·A·g/s³
∂r/∂c  = 1                  with g = exp(−((x−xp)²+(y−yp)²)/(2s²))
```

### 2.2 Optimizer

* `gsl_multifit_fdfsolver_lmsder` (scaled LM, MINPACK `lmder`) initialised at
  the free entries of `prmVect`.
* loop: `iterate`; stop when `gsl_multifit_test_delta(dx, x, eAbs, eRel)`
  returns `GSL_SUCCESS` (`|dx_i| < eAbs + eRel·|x_i|` for all i), or on an
  iterate error, or when `iter >= maxIter`.
* Defaults (as in the original MEX): `maxIter = 500`, `eAbs = 1e-8`, `eRel = 1e-8`.
* after the fit `s = |s|`.

### 2.3 Outputs

* `prmVect` – full 5-vector with the fitted values substituted.
* `prmStd` – `sqrt(RSS/(nValid − np − 1) · C_ii)` with `C = gsl_multifit_covar(J, 0)`
  (`J` = Jacobian at the final iterate, `np` = number of free parameters).
  **Verified against the original MEX** (`tests/reference_values.txt`, generated
  by `reference_tools/gen_unit_reference.m` on this machine): with n = 169 and
  np = 4 the divisor 164 reproduces `prmStd` to < 1e-9 relative, the divisor 165
  would be 0.3 % off. The full covariance matrix `C` also matches.
* Iteration stop: GSL's `lmsder` `iterate` returns `GSL_ETOLF` (29) once the
  predicted and actual reductions fall below its internal `ftol`; the loop then
  stops exactly like the MEX loop (`if (status) break;`). The MEX results for the
  synthetic windows ('xyAc', 'xyasc', 'Ac', masked window, x-shifted window) are
  reproduced to 1e-9 relative or better, which is only possible if the iteration
  path is identical.
* `res.data` – residual image (NaN where masked), `res.RSS = Σ r²`,
  `res.mean = mean(r)`, `res.std = sqrt(Σ(r−mean)²/(nValid−1))`.
* `res.hAD` – Anderson–Darling normality decision on the residuals.
  **Empirically established** (940 windows fitted with the original MEX, residuals
  and `hAD` recorded in `tests/ad_reference*.txt`, 540 of them deliberately
  within ±0.35 of the threshold): the MEX computes
  `A² = −n − (1/n)Σ(2i−1)(ln z_i + ln(1 − z_{n+1−i}))`, `z = normcdf(sorted r, mean, std)`
  with `std` = sample std (n−1), applies **no** small-sample correction, and
  rejects normality when `A² > 2.308`. That is the α = 0.05 critical value of the
  "μ known, σ estimated" row of the table in `adtest1.m` (`ctable(3,:)`), *not*
  the "both estimated" case (0.752 after correction) that `adtest1.m` would use.
  The admissible threshold interval from the data is (2.307741, 2.308986], so the
  identification is unambiguous; using the sample mean or μ = 0 is
  indistinguishable on this path because `c` is always a free parameter (mean
  residual ≈ 1e-15). All 940 MEX decisions are reproduced (unit test).
* `res.pval` – Kolmogorov–Smirnov p-value: **not used anywhere on the path,
  not ported** (`fitGaussians2D.m` only reads `std`, `RSS`, `hAD`).

### 2.4 GSL version and |σ|

The original `.mexw64` was built against GSL 1.x. Two details were pinned down
with `reference_tools/probe_mex_options.m` (357 real candidate windows of a
256×256 crop, free-σ mode `'xyasc'`, MEX results compared with the port by
`tests/fit_probe.cpp`):

* **Default options** are `[500 1e-8 1e-8]`: `maxIter=1000/200/100/50` change the
  MEX result for 2–5 windows, `maxIter=500` reproduces the default bit-for-bit;
  the tolerances are irrelevant in practice (`iterate` returns `GSL_ETOLF/ETOLX`
  before `test_delta` ever succeeds, for `eAbs,eRel ∈ {0, 1e-10 … 1e-6}` all
  357 results are identical).
* **`lmsder` of GSL 1.x vs 2.x**: GSL 2.x rescaled the LM column scaling
  (`dnrm2`) and fixed the 1.x `update_diag()` quirk (which sums only the first
  `p` rows of `J`). For ill-conditioned windows this changes the iteration path
  and the solution basin. The port therefore compiles the **unmodified GSL 1.16
  `lmder.c/lmiterate.c/lmpar.c/lmset.c/lmutil.c/qrsolv.c/covar.c/convergence.c`**
  (+ the five `linalg` routines they call, extracted verbatim from GSL 1.16)
  under the `g116_` prefix (`third_party/gsl116`, GPL-3), next to GSL 2.8 which
  is still used for vectors/BLAS/CDFs.
* **|σ| in model and Jacobian**: the MEX evaluates the Gaussian with
  `fabs(sigma)` *and* uses `|σ|²`, `|σ|³` in the Jacobian (so the derivative
  w.r.t. the signed parameter flips sign for σ<0). With the mathematically
  "correct" signed derivative 10/356 free-σ fits landed in a different basin
  than the MEX; with |σ| every fit lands in the MEX basin. Result on the probe
  set: 350/356 within 1e-6 relative in all five parameters, the remaining 6 are
  the same solution to 1e-5 relative (four are 500-iteration non-converging
  windows where the path is chaotic; the two others differ by ~2e-8 absolute in
  a near-zero coordinate). Fixed-σ fits (`'xyAc'`, `'Ac'`) agree to ≤ 1e-8 relative.

---

## 3. `gmdistribution.fit` and `rng(1)` reconstruction (R2025a semantics)

`gmdistribution.fit(X, k, 'Options', statset('maxIter', 200))` with 1-D data.
Read from the toolbox source (`@gmdistribution/fit.m`, `private/gmcluster.m`,
`private/wdensity.m`, `private/estep.m`, `datasample.m`, `randsample.m`):

### 3.1 Algorithm

1. Errors (caught by cmeAnalysis → `sigma = mean(svect)`): `n <= d`, `n <= k`,
   `var(X) < eps(var(X))` (zero variance), ill-conditioned covariance during EM.
2. **Initialisation `'plus'` (k-means++, R2025a default; older releases used `'randSample'`)**:
   * `PComponents = 1/k`, `Sigma_j = var(X)` (sample variance, n−1) for every j.
   * `C(1) = X(randi(n))` (one `rand` draw).
   * for `ii = 2..k`: `minDist = min(minDist, ((X − C(ii−1))² / var(X)))`;
     `p = minDist/Σ minDist`; `C(ii) = X(i)` where `i` is the first index with
     `cumsum(p) > u` (`edges = min([0 cumsum(p)],1); edges(end)=1`, `histcounts(u, edges)`),
     one `rand` draw `u` per component (`internal.stats.wswor(w,1)`, verified identical
     to inverse-CDF over 2000 random trials).
   * degenerate branch (`Σ minDist == 0` or `Inf`) → `datasample(...,'Replace',false)`
     is unreachable in practice (zero variance is rejected earlier) – implemented
     with `randperm(n,k)` via rejection sampling (verified for `randperm(100,2)`
     only; the small-`n` MATLAB branch differs and is unverified).
3. EM (`gmcluster_learn`), `MaxIter = 200`, `TolFun = 1e-6`, `probtol = 1e-8`:
   * E-step (`wdensity`+`estep`): `log_lh(:,j) = −0.5·(((x−mu_j)/L_j)² + 2 log L_j) + log p_j − 0.5·log(2π)`,
     `L_j = chol(Sigma_j) = sqrt(Sigma_j)`; error if `Sigma_j <= 0` or `L_j < eps(L_j)`.
     `maxll = max_j log_lh`, `post = exp(log_lh − maxll)`, `density = Σ_j post`,
     `ll = Σ_i (log density_i + maxll_i)`, `post /= density`; then (full covariance
     → `setSmallProbtoZero`) `post(post < 1e-8) = 0` and renormalise.
   * convergence test **before** the M-step: `llDiff = ll − ll_old`;
     converged iff `llDiff >= 0 && llDiff < 1e-6·|ll|`.
   * M-step: `N_j = Σ_i post_ij`; skip j with `N_j == 0`; `mu_j = Σ post_ij x_i / N_j`;
     `Sigma_j = Σ post_ij (x_i − mu_j)² / N_j` (computed from `sqrt(post).*Xcentered` products,
     restricted to `post > 0` rows when fewer than `floor(0.4 n)` are non-zero –
     numerically identical up to summation order); `p_j = N_j / Σ N_j`.
   * `Iters = iter`, `NlogL = −ll` (from the last E-step, i.e. the one that converged).
4. `BIC = 2·NlogL + (k + k + (k−1))·log(n)`.
5. `Converged=false` only produces a (silenced) warning; the fit is still used.

### 3.2 Random number generator

* `rng(1)` → Mersenne Twister `mt19937ar` seeded with `init_genrand(1)`
  (MATLAB treats seed 0 specially as 5489; seed 1 is plain). Doubles are produced
  with `genrand_res53` (`(a·2^26 + b)/2^53`, `a = next>>5`, `b = next>>6`).
  Verified: `rng(1); rand(1,3)` = `0.417022004702574, 0.7203244934421581, 0.00011437481734488664`.
* `randi(n)` ≡ `ceil(n·rand)` (verified on 100 000 draws, n = 5000).
* `randperm(n,k)` for `k ≪ n` ≡ rejection sampling with `randi` (verified for `(100,2)`).
  Only needed in an unreachable branch.
* RNG consumption per `getGaussianPSFsigmaFromData` call: `1 + 2 + 3 = 6` draws
  (k = 1, 2, 3). Channels are processed in order, so channel 2's draws follow
  channel 1's. Nothing else on the path consumes the global stream
  (`parfor` workers use independent streams and the detection code is deterministic).

---

## 4. MATLAB semantics that the port reproduces literally

* Images are stored **column-major** (`img(y, x)`, `sub2ind` = `y + (x−1)·ny`).
  The C++ `Image` type is column-major with the same linear indexing so that
  `find`, `sub2ind`, `bwconncomp` pixel lists and candidate ordering match.
* All coordinates in the output tables are **1-based MATLAB pixel coordinates**
  (`x` along columns, `y` along rows), exactly what `detection_v2.mat` contains.
* `round` = half away from zero; `ceil`, `linspace` as in MATLAB
  (`round(linspace(1,L,nf))`).
* `bwconncomp`/`bwlabel` 8-connectivity; `find` and `PixelIdxList` in
  column-major order.
* `ordfilt2` zero padding is irrelevant because `locmax2d` zeroes the border strip.
* `conv2(...,'valid')` after symmetric-XT padding → same-size output.
* `tcdf(x, ν)` with real (non-integer) ν and `ν = 0 → NaN`; `norminv` via `erfcinv`.
  The port uses GSL (`gsl_cdf_tdist_P`, `gsl_cdf_ugaussian_Pinv`, `gsl_sf_erfc`);
  agreement with MATLAB is checked by unit tests to ≤ 1e-14.
* `min`/`max` over images ignore nothing (no NaNs in raw data).
* `imwrite(uint8(255*mask), 'tif', 'compression', 'lzw')` → 8-bit LZW TIFF,
  one directory per frame.

---

## 5. Scope decisions / deviations (all intentional, none change numbers)

| Item | Decision |
|---|---|
| interactive dialogs | replaced by `--input`, `--channels`, `--markers` |
| `detection_v2.mat` (HDF5 v7.3) | replaced by TSV files (§6); an optional Python converter to `.mat` may be added later; not the first priority |
| `parfor` | sequential first; OpenMP over frames/cells once numerically validated |
| figure output of the σ histogram | omitted |
| `Overwrite=false` skip logic | `--overwrite` flag; default recompute |
| KS p-value in the MEX | not computed (unused) |

---

## 6. Outputs of the C++ tool

Per data set (`<cell>/<master>/Detection/`):

* `detection_cpp.tsv` – one row per detection with every `frameInfo` field,
  per-channel columns suffixed with the channel name as given on the command line.
* `dmasks.tif` / `Masks/dmask_NNN.tif` – identical layout to MATLAB.

Per condition (`--output`): `detections_all.tsv` (all cells, all rows) and
`summary.tsv` with `source_image, A_<ch>, c_<ch>, hval_<ch>` (+ x/y).

Debug dumps (`--dump-dir`) write every intermediate that the golden-reference
MATLAB script (`reference_tools/dump_reference.m`) also writes, so the Python
comparison script can diff stage by stage.

---

## 7. Verification status

### 7.1 Unit level (`tests/unit_tests.cpp`, 266 checks, all passing)

| Component | Reference | Agreement |
|---|---|---|
| MT19937 stream, `randi`, `randperm(100,2)` | MATLAB `rng(1)` | bit-exact |
| `norminv`, `normcdf`, `tcdf` (35 (x,ν) pairs incl. ν=0.5, 1, 1e8) | MATLAB | ≤ 2e-13 relative (`tcdf`), ≤ 1e-15 otherwise |
| `adtest1.m` case 3 | MATLAB | A² to 1e-12 |
| MEX `res.hAD` | 940 MEX residual sets | 940/940 decisions |
| `padarrayXT` | MATLAB | exact |
| `fitGaussian2D` modes `xyAc`, `xyasc`, `Ac`, `xyac`, NaN-masked window, x-shifted spot | original MEX | ≤ 7e-9 relative in prm, prmStd, C, res.std, res.RSS |
| `fitGaussians2D` (`xyAc`, `Ac`) and `pointSourceDetection` on a 40×50 synthetic image | MATLAB | ≤ 7e-9 relative; identical candidate set |
| `gmdistribution.fit` k=1,2,3 on 260 samples, `rng(1)` | MATLAB R2025a | μ, Σ, π, NlogL, BIC to 1e-12; identical iteration counts; RNG state after the three fits identical |
| free-σ fits on 356 real candidate windows | original MEX (`reference_tools/probe_mex_options.m`) | 350 within 1e-6 relative, 6 same-basin (≤ 1e-5), see §2.4 |

### 7.2 Pipeline level, cropped data set (3 cells × 256×256, channels ch1/ch2, `rng(1)`)

MATLAB golden dump: `reference_tools/dump_reference.m` (original functions +
`psd_dump.m`); end-to-end: untouched `runDetection` via `run_matlab_reference.m`.

* **Prefilter / LoG / A_est / c_est images**: agree to ~1e-10 relative to the
  magnitude of the convolved sums (absolute differences ≤ 3e-5 on values of
  order 1e5; a handful of pixels near zero exceed the 1e-6 relative tolerance).
  MATLAB's `conv2` accumulation order could not be identified (64 % of pixels
  bit-exact with the best of eight tried orders), so this floor remains.
* **Prefilter mask, combined mask, candidate lists (`lm.tsv`)**: identical.
* **Sigma-estimation candidates** (`pointSourceDetection(img,1.5,'xyac')` +
  `'xyasc'` refit, 78 image instances): identical candidate sets and row counts
  in every image; per image exactly one free-σ refit differs (a window whose
  LM fit does not converge within 500 iterations, e.g. σ→0.23, A_pstd ≈ 3e7 —
  the same solution basin, 1e-3…1e-2 relative), plus one `pval_Ar` NaN/value
  difference in ch2. All other refit values agree to ≤ 1e-6 relative.
* **Sigma**: ch1 2.3705947708858108 vs 2.3705947638152245 (3e-9 relative),
  ch2 2.4053736017628071 vs 2.4053736020168084 (1e-10 relative); the
  difference comes entirely from the unstable refits above (with identical
  `svect` the EM is bit-identical).
* **Detections (`runDetection` end-to-end, `detection_cpp.tsv` vs `detection_matlab.tsv`)**:

  | cell | rows MATLAB | rows C++ | columns beyond tol (atol 1e-8, rtol 1e-6) |
  |---|---|---|---|
  | Cell1 | 240 | 240 | 2 rows: slave-channel localized (`xyAc`) fits with `x_pstd ≈ 134` (degenerate, A ≈ 14.8 vs 14.80, 4e-4 relative); 1 row: master `c` −4.60393591 vs −4.60392951 |
  | Cell2 | 208 | 208 | 1 row: slave A at 1.4e-6 relative; 3 rows: master `c` at ≤ 3e-6 relative |
  | Cell3 | 198 | 198 | 1 row: master `c` at 3.3e-6 relative |

  Candidate sets, `hval_Ar`, `hval_AD`, `isPSF`, `x_init/y_init`, `maskN`,
  `maskA`, `mask_Ar` and the masks are identical; x, y, A agree to ≤ 1e-8
  relative for every row. The `c` deviations are absolute differences of
  ~5e-6 on a background level that happens to be near zero (image intensities
  are O(1e3)); they stem from the ~1e-10 relative `conv2` floor propagated
  through `A_est/c_est` initialisation and the LM stop criterion. The
  degenerate slave fits are windows whose LM iteration is numerically
  unstable: it was verified with `reference_tools/probe_candidates.m` +
  `tests/fit_probe.cpp` that the C++ fitter reproduces the MEX result of such
  a window to 1e-6 when given MATLAB's exact inputs, and that perturbing σ by
  3e-9 or `A_est` by 1e-8 relative moves the MEX/C++ solution by the amounts
  seen. (An earlier build of the port, with the GSL 2.8 solver and the signed
  σ-derivative, produced 3 extra rows in Cell2 and one differing candidate per
  sigma-estimation image — those were fixed by the GSL 1.16 solver and the |σ|
  convention of §2.4.)

### 7.3 Full data set (10 cells × 1024×1024, 16-bit, channels ch1/ch2, `rng(1)`)

MATLAB references: `dump_reference.m` (stage dumps, 2488 s) and the untouched
`data = loadConditionData(...); rng(1); runDetection(data)` via
`run_matlab_reference.m` (both give identical detections). C++: `cme_detect --seed 1`.

* **Sigma**: ch1 1.920628260736333 vs 1.9206282611477374 (2e-10 relative),
  ch2 1.9235122532708524 vs 1.9235122540473961 (4e-10 relative).
* **Sigma-estimation candidates**: identical candidate sets in all 80 image
  instances; per image 1–2 unstable free-σ refits differ (as in §7.2).
* **Detections, end-to-end against `detection_v2.mat` (exported by `run_matlab_reference.m`)**:

  | cell | rows (MATLAB = C++) | master-channel rows beyond tol | slave-channel rows beyond tol | discrete mismatches |
  |---|---|---|---|---|
  | Cell1 | 5127 | 0 | 15 | 0 |
  | Cell2 | 4488 | 0 | 15 | 0 |
  | Cell3 | 4456 | 6 | 15 | 0 |
  | Cell4 | 4593 | 1 | 14 | 0 |
  | Cell5 | 4753 | 5 | 16 | 0 |
  | Cell6 | 4831 | 2 | 14 | 0 |
  | Cell7 | 4950 | 1 | 14 | 0 |
  | Cell8 | 5005 | 1 | 12 | 0 |
  | Cell9 | 5108 | 2 | 17 | 0 |
  | Cell10 | 5014 | 1 | 17 | 0 |
  | **total** | **48325** | **19 (0.04 %)** | **149 (0.3 %)** | **0** |

  "Discrete" covers `hval_Ar`, `hval_AD`, `isPSF` (both channels), `x_init`,
  `y_init`, `maskN`, `mask_Ar`. Row order, source-image correspondence and the
  number of detections are identical in every cell.
  The 19 master-channel rows are all `c` (background) values with |c| < 0.3 on
  images of intensity O(1e3) — absolute differences ≤ 5e-7, i.e. the `conv2`
  floor; x, y, A, A_pstd, RSS, sigma_r, pval_Ar of the master channel are within
  tolerance for all 48325 rows (max relative difference of A: 4.5e-8).
  The 149 slave-channel rows are localized (`xyAc`) fits with huge position
  uncertainties (`x_pstd` 19…143 px, i.e. the spot is not there in ch2), whose LM
  path is unstable; the largest amplitude deviation is 0.75 % relative.
* **Downstream filter** (`analyze_matlab.py` logic: `A_master > mean(c)+2·std(c)`
  and `hval_Ar_master == 1`), run with `tests/analyze_cpp.py` on the C++ TSV,
  on the MATLAB TSV, and with the user's own `analyze_matlab.py` on the real
  `detection_v2.mat`: **45362 / 48325 puncta kept in all three, identical row
  keys**; `A_ch1` of the kept rows agrees to ≤ 4.5e-8 relative, `A_ch2` to
  ≤ 7.5e-3 relative (92 of 45362 rows above 1e-6, the unstable slave fits above).

### 7.4 Runtime (same 10-cell data set, same machine: 12 logical cores, Windows 11)

| run | load | σ estimation | detection + slave fits | output | total wall-clock |
|---|---|---|---|---|---|
| MATLAB R2025a `loadConditionData` + `rng(1)` + `runDetection` (default 2-worker `parfor` pool; a 4-thread C++ job was running concurrently for most of it) | 1.7 s | — | — | — | **1903 s** |
| `cme_detect --threads 1` | 0.02 s | 2164 s | 283 s | 3.5 s | **2451 s** |
| `cme_detect --threads 12` (frame-level parallelism only) | 0.03 s | 324 s | 284 s | 3.7 s | **611 s** |
| `cme_detect --threads 12` (+ movie-level parallelism, final build) | 0.02 s | 321 s | 42 s | 3.3 s | **366 s** (5.2× MATLAB) |

All three C++ runs produce byte-identical `detections_all.tsv`. Single-threaded,
the port is about as fast as MATLAB per image (the free-σ refit of ~6000
candidates per 1024² image, ~120k LM iterations, dominates: ~11 s; the σ=1.5
detection pass ~9 s); the speed-up comes from parallelism, which MATLAB only
gets through its parallel pool. Obvious remaining hot spots (not optimised, to
keep the numerics untouched): the per-candidate allocation of the LM state,
`exp` evaluation in the model, and the Jacobian QR in GSL 1.16 code.
