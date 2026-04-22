# DLS Calibration Notes

## The calibration problem

The SLiC (Single Liposome Curvature) assay requires converting
fluorescence amplitude into physical liposome diameter. Because each
vesicle is a thin spherical shell with a uniform mole fraction of
fluorophore, the integrated fluorescence intensity scales as the
surface area (D^2), so:

    diameter = C_cal × sqrt(intensity)

The calibration reduces to finding the single scalar C_cal. Dynamic
light scattering (DLS) provides the independent size reference.


## Standard procedure: distribution overlay

The procedure used throughout the SLiC literature is:

1. Collect the DLS number-weighted size distribution P_DLS(D) on the
   instrument's fixed bin grid.
2. Collect sqrt(fluorescence intensity) for each detected punctum.
3. For each candidate k, convert sqrt(F) → trial diameters via
   D = sqrt(F) / k, histogram them on the same DLS bin grid, and
   compute the chi-squared residual between the two distributions.
4. Minimize over k. This is a one-dimensional optimization.

The conversion factor is then 1/k (in units of nm per sqrt(A)).

This procedure was introduced by Kunding et al. (Biophys J 2008) and
simplified by Hatzakis et al. (Nat Chem Biol 2009) as d = C_cal × sqrt(I).
Zeno et al. (2018, Nat Commun 9:4152) describe the calibration as a single
linear scaling factor calibrated against the DLS mean diameter (Supplementary
Figs 5–7). Johnson et al. (Commun Biol 2025) describe the procedure as
"overlaying" the sqrt(fluorescence) and DLS distributions. None of these
papers specify the exact algorithm for choosing the scaling factor beyond
noting that it is determined from the DLS reference.

Our implementation formalizes the overlay as a chi-squared minimization:
rebin the fluorescence data onto the DLS bin grid for each trial k and
minimize the residual. This is slightly more robust than pure ratio-of-means
because it fits the full distribution shape rather than just the first
moment. In practice the two give nearly the same answer for reasonably
symmetric distributions.


## What this pipeline implements

`dls_calibration.py` performs the distribution overlay as the primary
method. Ratio-of-means is also reported as a quick sanity check. The
two typically agree to within a few percent.

Optional bootstrap resampling (`--bootstrap N`) estimates the variance
of both methods by resampling the fluorescence puncta with replacement.


## Previous exploration: log-normal fitting

An earlier version of this script explored fitting log-normal
distributions to both the DLS and fluorescence data simultaneously,
inspired by a conversation with Wade Zeno about Horiba DLS instruments
that output two complementary analysis modes ("standard" cumulant and
"broad" inversion). A pseudo-DLS was constructed as a linear combination
of the Malvern number distribution and a lognormal generated from the
Z-average and PDI.

This approach did not improve on ratio-of-means for Malvern data because:

- The Malvern number/volume/intensity distributions are mathematical
  transformations of the same underlying autocorrelation analysis, not
  independent modes — the optimizer always set alpha → 1.0, ignoring
  the Z-average component entirely.
- The fluorescence sqrt(A) distribution is not cleanly log-normal
  (hard left edge from detection threshold, right tail from
  oligolamellar puncta).
- Bootstrap variance comparison showed ratio-of-means had lower CV.

When we subsequently traced the actual published calibration methods
through Zeno's own supplementary materials, we found they all use a
single-scalar approach (distribution overlay or peak matching), not
log-normal fitting. The log-normal code has been removed from the
pipeline. See git history for the previous implementation.


## References

- **Kunding et al. (2008)** "A fluorescence-based technique to construct
  size distributions from single object measurements." Biophys J
  95:1176–1188. The original calibration procedure with cryoTEM
  validation and PSF correction for large vesicles.

- **Hatzakis et al. (2009)** "How curved membranes recruit amphipathic
  helices and protein anchoring motifs." Nat Chem Biol 5:835–841.
  Simplified the calibration to d = C_cal × sqrt(I).

- **Bhatia et al. (2009)** "Amphipathic motifs in BAR domains are
  essential for membrane curvature sensing." EMBO J 28:3303–3314.
  The foundational SLiC paper. Same calibration.

- **Jensen et al. (2011)** "Membrane curvature sensing by amphipathic
  helices: a single liposome study using α-synuclein and annexin B12."
  J Biol Chem 286:42603–42614. Same calibration approach.

- **Zeno et al. (2018)** "Synergy between intrinsically disordered
  domains and structured proteins amplifies membrane curvature sensing."
  Nat Commun 9:4152. Supplementary Figs 5–7 show the calibration: a
  single linear scaling factor (e.g. 5.69 nm/sqrt(intensity)) calibrated
  against the DLS mean diameter of the extruded population.

- **Johnson et al. (2025)** "Lipid packing defects are necessary and
  sufficient for membrane binding of α-synuclein." Commun Biol 8:1179.
  Describes calibration as "overlaying" the sqrt(fluorescence) and DLS
  distributions.
