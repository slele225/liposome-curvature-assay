// Port of the fitGaussian2D MEX (Aguet, GPL) - see PORTING_NOTES.md §2.
//
// [prmVect prmStd C res J] = fitGaussian2D(data, prmVect, mode, options)
//   model  f(x,y) = A exp(-((x-xp)^2 + (y-yp)^2)/(2 s^2)) + c
//   origin at the centre of the odd square window; x along columns, y along rows
//   prmVect = [xp yp A s c]; mode = any subset of "xyasc" (case-insensitive)
//   options = [maxIter eAbs eRel] (default 500, 1e-8, 1e-8)
//   NaN pixels are masked.
#pragma once
#include "cme/image.hpp"
#include <array>
#include <string>
#include <vector>

namespace cme {

struct FitOptions {
    int maxIter = 500;
    double eAbs = 1e-8;
    double eRel = 1e-8;
};

struct FitResult {
    std::array<double, 5> prm{};      // [xp yp A s c]
    std::vector<double> prmStd;       // std of the estimated parameters (mode order x,y,A,s,c)
    std::vector<int> estIdx;          // indices (0..4) of estimated parameters
    double RSS = 0.0;
    double mean = 0.0;                // mean of residuals
    double std = 0.0;                 // sample std of residuals (n-1)
    bool hAD = false;                 // Anderson-Darling (case 3, alpha 0.05)
    double A2 = 0.0;                  // AD statistic (debug)
    int iterations = 0;               // LM iterations performed
    int status = 0;                   // final GSL status of test_delta / iterate
    int nValid = 0;
};

FitResult fitGaussian2D(const ImageD& window, const std::array<double, 5>& prm0, const std::string& mode,
                        const FitOptions& opt = FitOptions());

} // namespace cme
