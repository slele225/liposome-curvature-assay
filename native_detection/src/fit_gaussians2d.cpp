#include "cme/fit_gaussians2d.hpp"
#include "cme/morphology.hpp"
#include "cme/stats.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace cme {

namespace {
const double NaN = std::numeric_limits<double>::quiet_NaN();

template <typename T>
void filter_vec(std::vector<T>& v, const std::vector<char>& keep) {
    std::vector<T> out;
    out.reserve(v.size());
    for (std::size_t i = 0; i < v.size(); ++i) if (keep[i]) out.push_back(v[i]);
    v.swap(out);
}
} // namespace

void PStruct::resize(std::size_t np) {
    x.assign(np, NaN); y.assign(np, NaN); A.assign(np, NaN); s.assign(np, NaN); c.assign(np, NaN);
    x_pstd.assign(np, NaN); y_pstd.assign(np, NaN); A_pstd.assign(np, NaN); s_pstd.assign(np, NaN); c_pstd.assign(np, NaN);
    x_init.assign(np, NaN); y_init.assign(np, NaN);
    sigma_r.assign(np, NaN); SE_sigma_r.assign(np, NaN); RSS.assign(np, NaN);
    pval_Ar.assign(np, NaN);
    mask_Ar.assign(np, 0.0);
    hval_Ar.assign(np, 0); hval_AD.assign(np, 0);
    iters.assign(np, 0); A2.assign(np, NaN);
}

void PStruct::filter(const std::vector<char>& keep) {
    filter_vec(x, keep); filter_vec(y, keep); filter_vec(A, keep); filter_vec(s, keep); filter_vec(c, keep);
    filter_vec(x_pstd, keep); filter_vec(y_pstd, keep); filter_vec(A_pstd, keep); filter_vec(s_pstd, keep); filter_vec(c_pstd, keep);
    filter_vec(x_init, keep); filter_vec(y_init, keep);
    filter_vec(sigma_r, keep); filter_vec(SE_sigma_r, keep); filter_vec(RSS, keep);
    filter_vec(pval_Ar, keep); filter_vec(mask_Ar, keep);
    filter_vec(hval_Ar, keep); filter_vec(hval_AD, keep);
    filter_vec(iters, keep); filter_vec(A2, keep);
}

PStruct fitGaussians2D(const ImageD& img,
                       const std::vector<double>& x, const std::vector<double>& y,
                       const std::vector<double>& Ain, const std::vector<double>& sigmaIn,
                       const std::vector<double>& cin, const std::string& mode,
                       const FitGaussiansOptions& opt) {
    const std::size_t np = x.size();
    if (y.size() != np) throw std::invalid_argument("fitGaussians2D: x/y size mismatch");
    std::vector<double> sigma = sigmaIn;
    if (sigma.size() == 1) sigma.assign(np, sigmaIn[0]);
    if (sigma.size() != np) throw std::invalid_argument("fitGaussians2D: sigma size mismatch");
    if (!Ain.empty() && Ain.size() != np) throw std::invalid_argument("fitGaussians2D: A size mismatch");
    if (!cin.empty() && cin.size() != np) throw std::invalid_argument("fitGaussians2D: c size mismatch");

    const long ny = static_cast<long>(img.ny());
    const long nx = static_cast<long>(img.nx());

    ImageI32 labels;
    if (opt.mask) labels = bwlabel8(*opt.mask);
    else labels = ImageI32(img.ny(), img.nx(), 0);

    PStruct P;
    P.resize(np);

    std::vector<long> xi(np), yi(np);
    for (std::size_t p = 0; p < np; ++p) {
        xi[p] = static_cast<long>(std::round(x[p]));
        yi[p] = static_cast<long>(std::round(y[p]));
        P.x_init[p] = static_cast<double>(xi[p]);
        P.y_init[p] = static_cast<double>(yi[p]);
    }

    const double kLevel = norminv(1.0 - opt.alpha / 2.0);
    const double iRange0 = img.minval(), iRange1 = img.maxval();
    const double diffRange = iRange1 - iRange0;

    // estIdx = regexpi('xyAsc', ['[' mode ']'])
    std::vector<int> estIdx;
    {
        std::string m = mode;
        for (char& ch : m) ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
        const char ref[] = "xyasc";
        for (int i = 0; i < 5; ++i) if (m.find(ref[i]) != std::string::npos) estIdx.push_back(i);
    }

    double sigma_max = 0.0;
    for (double sv : sigma) sigma_max = std::max(sigma_max, sv);
    const int w2 = opt.confRadius ? *opt.confRadius : static_cast<int>(std::ceil(2.0 * sigma_max));
    const int w4 = opt.windowSize ? *opt.windowSize : static_cast<int>(std::ceil(4.0 * sigma_max));
    const int W = 2 * w4 + 1;

    // annular mask for background estimation (only used when c is empty)
    ImageU8 annularMask;
    if (cin.empty()) {
        annularMask = ImageU8(W, W, 0);
        const double rin = std::ceil(3.0 * sigma_max), rout = std::ceil(4.0 * sigma_max);
        for (int xx = -w4; xx <= w4; ++xx) {
            for (int yy = -w4; yy <= w4; ++yy) {
                const double r = std::sqrt(static_cast<double>(xx * xx + yy * yy));
                if (r <= rout && r >= rin) annularMask(yy + w4, xx + w4) = 1;
            }
        }
    }

    // g = exp(-(-w4:w4).^2/(2*sigma_max^2)); g = g'*g;  (used for mask_Ar)
    std::vector<double> g1(W);
    for (int i = 0; i < W; ++i) {
        const double d = static_cast<double>(i - w4);
        g1[i] = std::exp(-d * d / (2.0 * sigma_max * sigma_max));
    }
    std::vector<double> g2(static_cast<std::size_t>(W) * W);
    for (int xx = 0; xx < W; ++xx) for (int yy = 0; yy < W; ++yy) g2[yy + xx * W] = g1[yy] * g1[xx];

    std::vector<double> T(np, 0.0), df2(np, 0.0);
    ImageD window(W, W);

    for (std::size_t p = 0; p < np; ++p) {
        // ignore points in border (1-based test in MATLAB)
        if (!(xi[p] > w4 && xi[p] <= nx - w4 && yi[p] > w4 && yi[p] <= ny - w4)) continue;

        const long x0 = xi[p] - 1 - w4;   // 0-based window origin
        const long y0 = yi[p] - 1 - w4;

        // label mask window; own label (centre) -> 0
        const int centreLabel = labels(static_cast<std::size_t>(yi[p] - 1), static_cast<std::size_t>(xi[p] - 1));
        std::vector<char> otherComp(static_cast<std::size_t>(W) * W, 0);
        for (int xx = 0; xx < W; ++xx) {
            for (int yy = 0; yy < W; ++yy) {
                int l = labels(static_cast<std::size_t>(y0 + yy), static_cast<std::size_t>(x0 + xx));
                if (l == centreLabel) l = 0;
                otherComp[yy + xx * W] = (l != 0);
                window(yy, xx) = img(static_cast<std::size_t>(y0 + yy), static_cast<std::size_t>(x0 + xx));
            }
        }

        double c_init;
        if (cin.empty()) {
            double s = 0.0; std::size_t n = 0;
            for (std::size_t i = 0; i < annularMask.size(); ++i) {
                if (annularMask[i] == 1 && !otherComp[i]) { s += window[i]; ++n; }
            }
            c_init = s / static_cast<double>(n);   // mean([]) = NaN in MATLAB (n==0 -> NaN here as well)
        } else {
            c_init = cin[p];
        }

        // set other components to NaN
        int npx = 0;
        for (std::size_t i = 0; i < window.size(); ++i) {
            if (otherComp[i]) window[i] = NaN;
            if (std::isfinite(window[i])) ++npx;
        }
        if (npx < 10) continue;

        double A_init;
        if (Ain.empty()) {
            double m = -std::numeric_limits<double>::infinity();
            for (std::size_t i = 0; i < window.size(); ++i) if (!std::isnan(window[i])) m = std::max(m, window[i]);
            A_init = m - c_init;
        } else {
            A_init = Ain[p];
        }

        std::array<double, 5> prm0 = {x[p] - static_cast<double>(xi[p]), y[p] - static_cast<double>(yi[p]), A_init, sigma[p], c_init};
        FitResult fr = fitGaussian2D(window, prm0, mode);

        const double dx = fr.prm[0], dy = fr.prm[1];
        if (dx > -w2 && dx < w2 && dy > -w2 && dy < w2 && fr.prm[2] < 2.0 * diffRange) {
            P.x[p] = static_cast<double>(xi[p]) + dx;
            P.y[p] = static_cast<double>(yi[p]) + dy;
            P.A[p] = fr.prm[2];
            P.s[p] = fr.prm[3];
            P.c[p] = fr.prm[4];

            double stdVect[5] = {0, 0, 0, 0, 0};
            for (std::size_t k = 0; k < estIdx.size(); ++k) stdVect[estIdx[k]] = fr.prmStd[k];
            P.x_pstd[p] = stdVect[0];
            P.y_pstd[p] = stdVect[1];
            P.A_pstd[p] = stdVect[2];
            P.s_pstd[p] = stdVect[3];
            P.c_pstd[p] = stdVect[4];

            P.sigma_r[p] = fr.std;
            P.RSS[p] = fr.RSS;
            P.SE_sigma_r[p] = fr.std / std::sqrt(2.0 * (npx - 1));
            const double SE_sigma_r = P.SE_sigma_r[p] * kLevel;
            P.hval_AD[p] = fr.hAD ? 1 : 0;
            P.A2[p] = fr.A2;
            P.iters[p] = fr.iterations;

            const double sigma_A = stdVect[2];
            const double A_est = fr.prm[2];
            const double sA2 = sigma_A * sigma_A, se2 = SE_sigma_r * SE_sigma_r;
            df2[p] = (npx - 1) * (sA2 + se2) * (sA2 + se2) / (sA2 * sA2 + se2 * se2);
            const double scomb = std::sqrt((sA2 + se2) / npx);
            T[p] = (A_est - fr.std * kLevel) / scomb;
            int cnt = 0;
            for (double gv : g2) if (A_est * gv > fr.std * kLevel) ++cnt;
            P.mask_Ar[p] = cnt;
        }
    }
    for (std::size_t p = 0; p < np; ++p) {
        P.pval_Ar[p] = tcdf(-T[p], df2[p]);
        P.hval_Ar[p] = (P.pval_Ar[p] < opt.alphaT) ? 1 : 0;
    }
    return P;
}

} // namespace cme
