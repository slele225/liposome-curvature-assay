#include "cme/point_source_detection.hpp"
#include "cme/conv.hpp"
#include "cme/morphology.hpp"
#include "cme/stats.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>

namespace cme {

namespace {
constexpr double PI = 3.14159265358979323846;
}

PSDResult pointSourceDetection(const ImageD& img, double sigma, const PSDOptions& opt, PSDDebug* dbg) {
    PSDResult R;
    const std::size_t ny = img.ny(), nx = img.nx();

    // Gaussian kernel
    const int w = static_cast<int>(std::ceil(4.0 * sigma));
    const int L = 2 * w + 1;
    std::vector<double> g(L), u(L, 1.0), gx2(L);
    for (int i = 0; i < L; ++i) {
        const double x = static_cast<double>(i - w);
        g[i] = std::exp(-x * x / (2.0 * sigma * sigma));
        gx2[i] = g[i] * x * x;
    }

    // convolutions
    ImageD imgXT = padarrayXT_symmetric(img, w);
    ImageD imgXT2(imgXT.ny(), imgXT.nx());
    for (std::size_t i = 0; i < imgXT.size(); ++i) imgXT2[i] = imgXT[i] * imgXT[i];
    ImageD fg = conv2_sep_valid(g, g, imgXT);
    ImageD fu = conv2_sep_valid(u, u, imgXT);
    ImageD fu2 = conv2_sep_valid(u, u, imgXT2);

    // Laplacian of Gaussian
    ImageD cxx = conv2_sep_valid(g, gx2, imgXT);   // conv2(g, gx2, ...): columns with g, rows with gx2
    ImageD cyy = conv2_sep_valid(gx2, g, imgXT);   // conv2(gx2, g, ...)
    ImageD imgLoG(ny, nx);
    const double s2 = sigma * sigma, s4 = s2 * s2;
    for (std::size_t i = 0; i < imgLoG.size(); ++i) {
        double v = 2.0 * fg[i] / s2 - (cxx[i] + cyy[i]) / s4;
        imgLoG[i] = v / (2.0 * PI * s2);
    }

    // 2-D kernel sums
    const double n = static_cast<double>(L) * static_cast<double>(L);
    double gsum = 0.0, g2sum = 0.0;
    for (int xx = 0; xx < L; ++xx) {
        for (int yy = 0; yy < L; ++yy) {
            const double gv = g[yy] * g[xx];
            gsum += gv;
            g2sum += gv * gv;
        }
    }

    // solution to linear system
    ImageD A_est(ny, nx), c_est(ny, nx);
    const double denom = g2sum - gsum * gsum / n;
    for (std::size_t i = 0; i < A_est.size(); ++i) {
        A_est[i] = (fg[i] - gsum * fu[i] / n) / denom;
        c_est[i] = (fu[i] - A_est[i] * gsum) / n;
    }

    ImageU8 mask(ny, nx, 1);
    ImageD pvalImg;
    if (opt.prefilter) {
        // J = [g(:) ones(n,1)]; C = inv(J'*J)  (2x2 inverse computed like MATLAB's inv via LU)
        // J'J = [g2sum gsum; gsum n]
        const double a = g2sum, b = gsum, d = n;
        // MATLAB inv() on a 2x2 uses LU with partial pivoting; for this SPD
        // matrix the result equals the closed form to within rounding.
        const double det = a * d - b * b;
        const double C11 = d / det;

        const double kLevel = norminv(1.0 - opt.alpha / 2.0);
        pvalImg = ImageD(ny, nx);
        for (std::size_t i = 0; i < A_est.size(); ++i) {
            const double A = A_est[i], c = c_est[i];
            const double f_c = fu2[i] - 2.0 * c * fu[i] + n * c * c;
            double RSS = A * A * g2sum - 2.0 * A * (fg[i] - c * gsum) + f_c;
            if (RSS < 0) RSS = 0;
            const double sigma_e2 = RSS / (n - 3.0);
            const double sigma_A = std::sqrt(sigma_e2 * C11);
            const double sigma_res = std::sqrt(RSS / (n - 1.0));
            const double SE_sigma_c = sigma_res / std::sqrt(2.0 * (n - 1.0)) * kLevel;
            const double sA2 = sigma_A * sigma_A, sc2 = SE_sigma_c * SE_sigma_c;
            const double df2 = (n - 1.0) * (sA2 + sc2) * (sA2 + sc2) / (sA2 * sA2 + sc2 * sc2);
            const double scomb = std::sqrt((sA2 + sc2) / n);
            const double T = (A - sigma_res * kLevel) / scomb;
            const double pval = tcdf(-T, df2);
            pvalImg[i] = pval;
            mask[i] = (pval < 0.05) ? 1 : 0;
        }
    }
    if (dbg) {
        dbg->imgLoG = imgLoG;
        dbg->A_est = A_est;
        dbg->c_est = c_est;
        dbg->pvalPrefilter = pvalImg;
        dbg->maskPrefilter = mask;
    }

    // all local max
    ImageD allMax = locmax2d(imgLoG, 2 * static_cast<int>(std::ceil(sigma)) + 1);

    // local maxima above threshold in image domain
    ImageD imgLM(ny, nx);
    double lmSum = 0.0;
    for (std::size_t i = 0; i < imgLM.size(); ++i) {
        imgLM[i] = allMax[i] * (mask[i] ? 1.0 : 0.0);
        lmSum += imgLM[i];
    }

    R.mask = mask;
    if (lmSum == 0.0) {   // sum(imgLM(:)) ~= 0 test
        if (dbg) dbg->maskCombined = mask;
        return R;
    }

    if (opt.refineMaskLoG) {
        double logThreshold = std::numeric_limits<double>::infinity();
        for (std::size_t i = 0; i < imgLM.size(); ++i) if (imgLM[i] != 0.0) logThreshold = std::min(logThreshold, imgLoG[i]);
        for (std::size_t i = 0; i < mask.size(); ++i) if (imgLoG[i] >= logThreshold) mask[i] = 1;
        if (dbg) dbg->logThreshold = logThreshold;
    }
    // re-select local maxima
    for (std::size_t i = 0; i < imgLM.size(); ++i) imgLM[i] = allMax[i] * (mask[i] ? 1.0 : 0.0);
    if (opt.exclusionMask) {
        for (std::size_t i = 0; i < imgLM.size(); ++i) if ((*opt.exclusionMask)[i] == 0) imgLM[i] = 0.0;
    }
    if (dbg) dbg->maskCombined = mask;

    // [lmy, lmx] = find(imgLM~=0)  -> column-major order
    std::vector<double> lmx, lmy;
    std::vector<std::size_t> lmIdx;
    for (std::size_t i = 0; i < imgLM.size(); ++i) {
        if (imgLM[i] != 0.0) {
            lmIdx.push_back(i);
            lmy.push_back(static_cast<double>(i % ny) + 1.0);
            lmx.push_back(static_cast<double>(i / ny) + 1.0);
        }
    }
    if (dbg) { dbg->lmx = lmx; dbg->lmy = lmy; }
    R.mask = mask;
    if (lmIdx.empty()) return R;

    std::vector<double> Ainit(lmIdx.size()), cinit(lmIdx.size()), sig(lmIdx.size(), sigma);
    for (std::size_t k = 0; k < lmIdx.size(); ++k) { Ainit[k] = A_est[lmIdx[k]]; cinit[k] = c_est[lmIdx[k]]; }

    FitGaussiansOptions fo;
    fo.alpha = opt.alpha;
    fo.mask = &mask;
    PStruct P = fitGaussians2D(img, lmx, lmy, Ainit, sig, cinit, opt.mode, fo);
    if (dbg) dbg->fitAll = P;

    // remove NaN values
    std::vector<char> keep(P.size());
    std::size_t nKeep = 0;
    for (std::size_t i = 0; i < P.size(); ++i) { keep[i] = !std::isnan(P.x[i]); nKeep += keep[i]; }
    if (dbg) dbg->keepNonNaN = keep;
    if (nKeep == 0) return R;
    P.filter(keep);

    // significant amplitudes
    std::vector<char> idx(P.size());
    for (std::size_t i = 0; i < P.size(); ++i) idx[i] = (P.hval_Ar[i] == 1);

    // eliminate duplicate positions (resulting from localization)
    if (opt.removeRedundant) {
        const double r = opt.redundancyRadius;
        const double r2 = r * r;
        const std::size_t N = P.size();
        // grid-accelerated ball query, equivalent to KDTreeBallQuery(pM, pM, r)
        std::map<std::pair<long, long>, std::vector<std::size_t>> grid;
        auto cell = [&](double v) { return static_cast<long>(std::floor(v / r)); };
        for (std::size_t i = 0; i < N; ++i) grid[{cell(P.x[i]), cell(P.y[i])}].push_back(i);
        for (std::size_t k = 0; k < N; ++k) {
            std::vector<std::size_t> nb;
            const long cx = cell(P.x[k]), cy = cell(P.y[k]);
            for (long dx = -1; dx <= 1; ++dx) {
                for (long dy = -1; dy <= 1; ++dy) {
                    auto it = grid.find({cx + dx, cy + dy});
                    if (it == grid.end()) continue;
                    for (std::size_t j : it->second) {
                        const double ddx = P.x[j] - P.x[k], ddy = P.y[j] - P.y[k];
                        if (ddx * ddx + ddy * ddy <= r2) nb.push_back(j);
                    }
                }
            }
            if (nb.size() > 1) {
                double minRSS = std::numeric_limits<double>::infinity();
                for (std::size_t j : nb) minRSS = std::min(minRSS, P.RSS[j]);
                for (std::size_t j : nb) if (P.RSS[j] != minRSS) idx[j] = 0;
            }
        }
    }
    if (dbg) dbg->keepFinal = idx;

    std::size_t nFinal = 0;
    for (char v : idx) nFinal += v;
    if (nFinal == 0) return R;
    P.filter(idx);
    R.pstruct = P;
    R.isPSF.resize(P.size());
    for (std::size_t i = 0; i < P.size(); ++i) R.isPSF[i] = P.hval_AD[i] ? 0 : 1;
    R.empty = false;

    if (opt.refineMaskValid) {
        ConnComp cc = bwconncomp8(mask);
        ImageI32 labels = labelmatrix(cc, ny, nx);
        std::vector<char> used(cc.numObjects + 1, 0);
        for (std::size_t i = 0; i < P.size(); ++i) {
            const std::size_t yy = static_cast<std::size_t>(P.y_init[i]) - 1;
            const std::size_t xx = static_cast<std::size_t>(P.x_init[i]) - 1;
            used[static_cast<std::size_t>(labels(yy, xx))] = 1;
        }
        ImageU8 m2(ny, nx, 0);
        for (std::size_t i = 0; i < m2.size(); ++i) {
            const int l = labels[i];
            m2[i] = (l != 0 && used[static_cast<std::size_t>(l)]) ? 1 : 0;
        }
        R.mask = m2;
    }
    return R;
}

} // namespace cme
