// Port of fitGaussians2D.m (wrapper that fits many candidates in an image).
#pragma once
#include "cme/image.hpp"
#include "cme/fit_gaussian2d.hpp"
#include <string>
#include <vector>
#include <optional>

namespace cme {

// pStruct fields of fitGaussians2D.m; all vectors have length np and use NaN /
// false / 0 for points that were not fitted (exactly like the MATLAB code).
struct PStruct {
    std::vector<double> x, y, A, s, c;
    std::vector<double> x_pstd, y_pstd, A_pstd, s_pstd, c_pstd;
    std::vector<double> x_init, y_init;          // rounded initial positions (1-based)
    std::vector<double> sigma_r, SE_sigma_r, RSS;
    std::vector<double> pval_Ar;
    std::vector<double> mask_Ar;
    std::vector<char> hval_Ar, hval_AD;
    // extra debug fields (not in MATLAB): LM iterations, AD statistic
    std::vector<int> iters;
    std::vector<double> A2;

    std::size_t size() const { return x.size(); }
    void resize(std::size_t np);
    // keep only the entries with keep[i] != 0 (MATLAB logical indexing)
    void filter(const std::vector<char>& keep);
};

struct FitGaussiansOptions {
    double alpha = 0.05;    // for kLevel
    double alphaT = 0.05;   // for hval_Ar
    const ImageU8* mask = nullptr;  // optional logical mask (bwlabel is applied)
    std::optional<int> confRadius;  // default ceil(2*sigma_max)
    std::optional<int> windowSize;  // default ceil(4*sigma_max)
    int threads = 1;                // > 1: fit candidates in parallel (each candidate is independent; results by index)
};

// x, y : initial (or fixed) positions, 1-based MATLAB coordinates
// A, c : initial (or fixed) values; empty vectors mean "estimate from window"
// sigma: per-point sigma (scalar broadcast handled by caller)
PStruct fitGaussians2D(const ImageD& img,
                       const std::vector<double>& x, const std::vector<double>& y,
                       const std::vector<double>& A, const std::vector<double>& sigma,
                       const std::vector<double>& c, const std::string& mode,
                       const FitGaussiansOptions& opt = FitGaussiansOptions());

} // namespace cme
