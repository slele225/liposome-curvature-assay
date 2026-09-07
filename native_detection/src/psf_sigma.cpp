#include "cme/psf_sigma.hpp"
#include "cme/point_source_detection.hpp"
#include "cme/profile.hpp"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <cstdio>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace cme {

double getGaussianPSFsigmaFromData(const std::vector<const ImageD*>& images, MatlabTwister& rng, PsfSigmaDebug* dbg,
                                   int threads, ParLevel level) {
    const std::size_t nd = images.size();
    std::vector<std::vector<double>> svectPer(nd);
    std::vector<PStruct> refits(nd);

    // Identical image instances (same pointer) yield identical results:
    // process each distinct image once, in order of first occurrence.
    std::vector<std::size_t> firstOf(nd), unique;
    for (std::size_t i = 0; i < nd; ++i) {
        std::size_t j = 0;
        while (j < i && images[j] != images[i]) ++j;
        firstOf[i] = j;
        if (j == i) unique.push_back(i);
    }
    const std::size_t nu = unique.size();

    // every image is processed independently; results are stored by index so
    // the concatenation order (and therefore the RNG-dependent GMM) is unchanged
    const int nthreads = threads > 0 ? threads : 1;
    const int outer = (level == ParLevel::Movie) ? nthreads : 1;
    const int inner = (level == ParLevel::Movie) ? 1 : nthreads;
    auto process_image = [&](std::size_t uu) {
        const std::size_t i = unique[uu];
#ifdef _OPENMP
        if (uu == 0) std::fprintf(stderr, "[sigma estimation: %d images (%d distinct), %d OpenMP threads over %s]\n",
                                  static_cast<int>(nd), static_cast<int>(nu), nthreads,
                                  level == ParLevel::Movie ? "images" : "candidates");
#endif
        const ImageD& img = *images[i];
        // First pass with fixed sigma
        PSDOptions o;
        o.mode = "xyac";
        o.threads = inner;
        PSDResult r;
        {
            prof::ContextGuard cg(prof::CtxSigmaDetect);
            r = pointSourceDetection(img, 1.5, o);
        }
        if (!r.empty) {
            const std::size_t np = r.pstruct.size();
            std::vector<double> sig(np, 1.5);
            prof::ContextGuard cg(prof::CtxSigmaRefit);
            FitGaussiansOptions fo;
            fo.threads = inner;
            PStruct P = fitGaussians2D(img, r.pstruct.x, r.pstruct.y, r.pstruct.A, sig, r.pstruct.c, "xyasc", fo);
            for (std::size_t q = 0; q < np; ++q) {
                const bool isPSF = (P.hval_AD[q] == 0) && (P.pval_Ar[q] < 0.05);
                if (!std::isnan(P.s[q]) && isPSF) svectPer[i].push_back(P.s[q]);
            }
            refits[i] = P;
        }
    };
    if (outer > 1) {
        #pragma omp parallel for schedule(dynamic) num_threads(outer)
        for (long uu = 0; uu < static_cast<long>(nu); ++uu) process_image(static_cast<std::size_t>(uu));
    } else {
        for (std::size_t uu = 0; uu < nu; ++uu) process_image(uu);
    }
    for (std::size_t i = 0; i < nd; ++i) {
        if (firstOf[i] != i) { svectPer[i] = svectPer[firstOf[i]]; refits[i] = refits[firstOf[i]]; }
    }
    std::vector<double> svect;
    for (auto& v : svectPer) svect.insert(svect.end(), v.begin(), v.end());
    if (dbg) { dbg->svectPerImage = svectPer; dbg->refitPerImage = refits; dbg->svect = svect; }

    // gmdistribution.fit for n = 1..3, min BIC, choose highest-peak component
    double sigma;
    prof::Scoped pGmm(prof::Gmm);
    try {
        std::vector<GmmResult> fits;
        for (int k = 1; k <= 3; ++k) {
            fits.push_back(gmdistribution_fit_1d(svect, k, rng, 200));
        }
        std::size_t best = 0;
        for (std::size_t j = 1; j < fits.size(); ++j) if (fits[j].BIC < fits[best].BIC) best = j;
        const GmmResult& g = fits[best];
        // [mu, idx] = sort(obj.mu); svec = sqrt(Sigma(idx)); amp = PComponents(idx)
        std::vector<std::size_t> order(g.mu.size());
        std::iota(order.begin(), order.end(), 0);
        std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return g.mu[a] < g.mu[b]; });
        double bestPeak = -std::numeric_limits<double>::infinity();
        std::size_t bestC = 0;
        for (std::size_t j = 0; j < order.size(); ++j) {
            const double svec = std::sqrt(g.Sigma[order[j]]);
            const double peak = g.PComponents[order[j]] / (std::sqrt(2.0 * 3.14159265358979323846) * svec);
            if (peak > bestPeak) { bestPeak = peak; bestC = j; }   // max() returns the first maximum
        }
        sigma = g.mu[order[bestC]];
        if (dbg) { dbg->gmm = fits; dbg->chosenK = g.k; dbg->chosenComponent = static_cast<int>(bestC) + 1; dbg->gmmFailed = false; }
    } catch (const GmmFitError&) {
        double m = 0.0;
        for (double v : svect) m += v;
        sigma = svect.empty() ? std::numeric_limits<double>::quiet_NaN() : m / static_cast<double>(svect.size());
        if (dbg) dbg->gmmFailed = true;
    }
    return sigma;
}

} // namespace cme
