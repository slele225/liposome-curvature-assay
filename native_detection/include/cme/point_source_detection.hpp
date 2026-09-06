// Port of pointSourceDetection.m (default options, no mixture fitting).
#pragma once
#include "cme/image.hpp"
#include "cme/fit_gaussians2d.hpp"
#include <string>
#include <vector>
#include <optional>

namespace cme {

struct PSDOptions {
    std::string mode = "xyAc";
    double alpha = 0.05;
    const ImageU8* exclusionMask = nullptr;   // 'Mask' option (cell mask); nullptr = all
    bool removeRedundant = true;
    double redundancyRadius = 0.25;
    bool prefilter = true;
    bool refineMaskLoG = true;
    bool refineMaskValid = true;
};

// Intermediates that the MATLAB reference dump also records.
struct PSDDebug {
    ImageD imgLoG;
    ImageD A_est, c_est;
    ImageD pvalPrefilter;          // t-test p-values (Prefilter)
    ImageU8 maskPrefilter;         // pval < 0.05
    ImageU8 maskCombined;          // after RefineMaskLoG (mask used for fitting)
    double logThreshold = 0.0;
    std::vector<double> lmx, lmy;  // candidate local maxima (1-based), MATLAB find() order
    PStruct fitAll;                // fitGaussians2D output before any filtering
    std::vector<char> keepNonNaN;  // ~isnan(x)
    std::vector<char> keepFinal;   // hval_Ar & non-redundant (on the non-NaN subset)
};

struct PSDResult {
    bool empty = true;             // pstruct = []
    PStruct pstruct;               // final detections (all have hval_Ar == 1)
    std::vector<char> isPSF;       // ~hval_AD
    ImageU8 mask;                  // returned mask (after RefineMaskValid)
};

PSDResult pointSourceDetection(const ImageD& img, double sigma, const PSDOptions& opt = PSDOptions(),
                               PSDDebug* dbg = nullptr);

} // namespace cme
