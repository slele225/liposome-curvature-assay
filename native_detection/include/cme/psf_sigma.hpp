// Port of getGaussianPSFsigmaFromData.m (default path, Display ignored).
#pragma once
#include "cme/image.hpp"
#include "cme/gmm.hpp"
#include "cme/fit_gaussians2d.hpp"
#include "cme/mt19937.hpp"
#include "cme/parallel.hpp"
#include <vector>

namespace cme {

struct PsfSigmaDebug {
    std::vector<std::vector<double>> svectPerImage;   // per input image, in order
    std::vector<PStruct> refitPerImage;               // 'xyasc' refit tables per image (may be empty)
    std::vector<double> svect;                        // concatenated
    std::vector<GmmResult> gmm;                       // k = 1,2,3 (only those that ran)
    bool gmmFailed = false;                           // -> fallback mean(svect)
    int chosenK = 0;
    int chosenComponent = 0;                          // 1-based after sorting by mu
};

// Returns the estimated sigma; consumes the RNG exactly like MATLAB.
//
// `images` may contain the same pointer several times (MATLAB samples
// round(linspace(1, L, nf)) frames per movie, which repeats frames when the
// movie is short).  The detection + refit of an image is a deterministic
// function of its pixels alone (no RNG involved), so repeated pointers are
// processed once and the per-image results are replicated; the concatenated
// svect, and therefore the GMM input and RNG consumption, are unchanged.
double getGaussianPSFsigmaFromData(const std::vector<const ImageD*>& images, MatlabTwister& rng,
                                   PsfSigmaDebug* dbg = nullptr, int threads = 1,
                                   ParLevel level = ParLevel::Candidate);

} // namespace cme
