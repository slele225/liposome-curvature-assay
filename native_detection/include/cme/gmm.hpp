// 1-D Gaussian mixture EM reproducing MATLAB R2025a
// gmdistribution.fit(X, k, 'Options', statset('maxIter', 200)) with default
// 'Start' = 'plus' (k-means++), 'CovType' = 'full', 'SharedCov' = false,
// 'Regularize' = 0, TolFun = 1e-6, ProbabilityTolerance = 1e-8.
#pragma once
#include "cme/mt19937.hpp"
#include <vector>
#include <stdexcept>
#include <string>

namespace cme {

struct GmmFitError : public std::runtime_error {
    explicit GmmFitError(const std::string& m) : std::runtime_error(m) {}
};

struct GmmResult {
    int k = 0;
    std::vector<double> mu;          // k
    std::vector<double> Sigma;       // k (variances)
    std::vector<double> PComponents; // k
    double NlogL = 0.0;
    double BIC = 0.0;
    int iters = 0;
    bool converged = false;
    std::vector<long> initIdx;       // 1-based indices chosen by k-means++ (debug)
    std::vector<double> initMu;      // initial means (debug)
};

// Throws GmmFitError for every condition in which MATLAB's fit would throw.
GmmResult gmdistribution_fit_1d(const std::vector<double>& X, int k, MatlabTwister& rng,
                                int maxIter = 200, double tolFun = 1e-6, double probTol = 1e-8);

} // namespace cme
