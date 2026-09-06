// Port of runDetection.m (default options) including the per-frame main().
#pragma once
#include "cme/condition_data.hpp"
#include "cme/image.hpp"
#include "cme/fit_gaussians2d.hpp"
#include "cme/psf_sigma.hpp"
#include "cme/point_source_detection.hpp"
#include <string>
#include <vector>

namespace cme {

struct RunOptions {
    std::vector<double> sigma;      // empty -> estimate from data
    bool removeRedundant = true;
    double alpha = 0.05;
    unsigned seed = 1;              // rng(seed)
    std::string dumpDir;            // optional intermediate dumps
    bool writeMasks = true;
    bool writeMatlabLayout = true;  // write Detection/ under the master channel like MATLAB
    int threads = 1;
};

// One frame's frameInfo entry (multi-channel fields are nCh x np, stored per channel).
struct FrameInfo {
    std::size_t frame = 0;
    std::size_t np = 0;
    std::size_t nCh = 0;
    // per channel vectors [c][p]
    std::vector<std::vector<double>> x, y, A, c, x_pstd, y_pstd, A_pstd, c_pstd, sigma_r, SE_sigma_r, RSS, pval_Ar;
    std::vector<std::vector<char>> hval_Ar, hval_AD, isPSF;
    std::vector<double> x_init, y_init, maskA, maskN, mask_Ar;
    std::vector<double> s;                    // sigma per channel
    std::vector<std::pair<double, double>> dRange;   // per channel [min max]
    // debug: slave fixed / localized fits and the chosen index
    std::vector<PStruct> slaveFixed, slaveLoc;   // indexed by channel (empty for master)
    std::vector<std::vector<char>> slaveUseLoc;
    std::vector<char> removedBySlave;            // over the master detections (before removal)
};

struct DetectionOutput {
    MovieData data;
    std::vector<FrameInfo> frames;
    std::vector<ImageU8> masks;
};

struct SigmaEstimate {
    std::vector<double> sigma;                 // per channel (after 1.1 clamp)
    std::vector<double> sigmaRaw;              // before clamp
    std::vector<PsfSigmaDebug> debug;          // per channel
};

SigmaEstimate estimate_sigma(const std::vector<MovieData>& data, const RunOptions& opt);

DetectionOutput run_detection_movie(const MovieData& d, const std::vector<double>& sigma, const RunOptions& opt,
                                    std::size_t movieIndex);

// runDetection(data): sigma estimation + all movies. Returns per movie outputs.
std::vector<DetectionOutput> runDetection(const std::vector<MovieData>& data, const RunOptions& opt, SigmaEstimate* sigmaOut = nullptr);

} // namespace cme
