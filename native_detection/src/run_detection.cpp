#include "cme/run_detection.hpp"
#include "cme/tiff_io.hpp"
#include "cme/morphology.hpp"
#include "cme/dump.hpp"
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <limits>
#include <stdexcept>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace fs = std::filesystem;

namespace cme {

namespace {
const double NaN = std::numeric_limits<double>::quiet_NaN();

// round(linspace(1, L, nf)) with MATLAB's linspace formula
std::vector<std::size_t> sample_frames(std::size_t L, std::size_t nf) {
    std::vector<std::size_t> out(nf);
    const double a = 1.0, b = static_cast<double>(L);
    const double n1 = static_cast<double>(nf) - 1.0;
    for (std::size_t i = 0; i < nf; ++i) {
        double v;
        if (nf == 1) v = b;   // MATLAB linspace(a,b,1) returns b
        else v = a + (b - a) / n1 * static_cast<double>(i);
        // MATLAB linspace computes y = d1 + (0:n1)*(d2-d1)/n1 and sets y(end)=d2
        if (i + 1 == nf) v = b;
        out[i] = static_cast<std::size_t>(std::round(v));
    }
    return out;
}

ImageD load_frame(const MovieData& d, std::size_t c, std::size_t f) {
    FrameRef r = frame_ref(d, c, f);
    return read_tiff_frame_double(r.path, r.directory);
}
} // namespace

SigmaEstimate estimate_sigma(const std::vector<MovieData>& data, const RunOptions& opt) {
    SigmaEstimate S;
    const std::size_t nd = data.size();
    std::size_t nCh = data.front().channels.size();
    for (const auto& d : data) if (d.channels.size() != nCh) throw std::runtime_error("mismatch in channel count between data sets");
    MatlabTwister rng(opt.seed);
    S.sigma.assign(nCh, 0.0);
    S.sigmaRaw.assign(nCh, 0.0);
    S.debug.resize(nCh);
    const std::size_t nf = static_cast<std::size_t>(std::round(40.0 / static_cast<double>(nd)));
    std::cout << "Determining Gaussian PSF parameters from data ... " << std::flush;
    for (std::size_t c = 0; c < nCh; ++c) {
        // frames = cell(nd, nf); frames(i,:) = ...; vertcat(frames(:)) -> column-major
        std::vector<std::vector<ImageD>> frames(nd);
        for (std::size_t i = 0; i < nd; ++i) {
            auto fidx = sample_frames(data[i].movieLength, nf);
            for (std::size_t f : fidx) frames[i].push_back(load_frame(data[i], c, f));
        }
        std::vector<const ImageD*> list;
        for (std::size_t f = 0; f < nf; ++f)
            for (std::size_t i = 0; i < nd; ++i) list.push_back(&frames[i][f]);
        S.sigmaRaw[c] = getGaussianPSFsigmaFromData(list, rng, &S.debug[c], opt.threads);
        S.sigma[c] = S.sigmaRaw[c];
    }
    std::cout << "done.\n";
    std::cout << "Gaussian PSF s.d. values:";
    for (double s : S.sigma) std::printf(" %.2f", s);
    std::cout << "\n";
    for (double& s : S.sigma) {
        if (s < 1.1) {
            std::cerr << "Sigma values < 1.1 were rounded to 1.1 to avoid poor localization performance.\n";
            s = 1.1;
        }
    }
    return S;
}

static void fill_frame_from_master(FrameInfo& F, const PStruct& P, const std::vector<char>& isPSF, std::size_t nCh, std::size_t mCh) {
    const std::size_t np = P.size();
    F.np = np;
    F.nCh = nCh;
    auto initD = [&](std::vector<std::vector<double>>& v) { v.assign(nCh, std::vector<double>(np, NaN)); };
    auto initL = [&](std::vector<std::vector<char>>& v) { v.assign(nCh, std::vector<char>(np, 0)); };
    initD(F.x); initD(F.y); initD(F.A); initD(F.c); initD(F.x_pstd); initD(F.y_pstd); initD(F.A_pstd); initD(F.c_pstd);
    initD(F.sigma_r); initD(F.SE_sigma_r); initD(F.RSS); initD(F.pval_Ar);
    initL(F.hval_Ar); initL(F.hval_AD); initL(F.isPSF);
    F.x[mCh] = P.x; F.y[mCh] = P.y; F.A[mCh] = P.A; F.c[mCh] = P.c;
    F.x_pstd[mCh] = P.x_pstd; F.y_pstd[mCh] = P.y_pstd; F.A_pstd[mCh] = P.A_pstd; F.c_pstd[mCh] = P.c_pstd;
    F.sigma_r[mCh] = P.sigma_r; F.SE_sigma_r[mCh] = P.SE_sigma_r; F.RSS[mCh] = P.RSS; F.pval_Ar[mCh] = P.pval_Ar;
    F.hval_Ar[mCh] = P.hval_Ar; F.hval_AD[mCh] = P.hval_AD; F.isPSF[mCh] = isPSF;
    F.x_init = P.x_init; F.y_init = P.y_init; F.mask_Ar = P.mask_Ar;
    F.maskA.assign(np, NaN); F.maskN.assign(np, NaN);
}

template <typename T>
static void remove_cols(std::vector<T>& v, const std::vector<char>& rm) {
    std::vector<T> out;
    for (std::size_t i = 0; i < v.size(); ++i) if (!rm[i]) out.push_back(v[i]);
    v.swap(out);
}

DetectionOutput run_detection_movie(const MovieData& d, const std::vector<double>& sigma, const RunOptions& opt, std::size_t movieIndex) {
    DetectionOutput out;
    out.data = d;
    const std::size_t nCh = d.channels.size();
    const std::size_t mCh = 0;   // master = first channel (data.source == channels{1})
    if (sigma.size() != nCh) throw std::runtime_error("sigma must have one entry per channel");
    const std::size_t L = d.movieLength;
    out.frames.resize(L);
    out.masks.resize(L);

    const bool dump = !opt.dumpDir.empty();
    std::string movieDump;
    if (dump) {
        movieDump = opt.dumpDir + "/" + dump_movie_name(d, movieIndex);
        fs::create_directories(movieDump);
    }

    #pragma omp parallel for schedule(dynamic) num_threads(opt.threads > 0 ? opt.threads : 1)
    for (long kk = 1; kk <= static_cast<long>(L); ++kk) {
        const std::size_t k = static_cast<std::size_t>(kk);
        FrameInfo& F = out.frames[k - 1];
        F.frame = k;
        F.nCh = nCh;
        F.s = sigma;
        F.dRange.assign(nCh, {NaN, NaN});
        F.slaveFixed.resize(nCh); F.slaveLoc.resize(nCh); F.slaveUseLoc.resize(nCh);

        ImageD img = load_frame(d, mCh, k);
        PSDOptions po;
        po.alpha = opt.alpha;
        po.removeRedundant = opt.removeRedundant;
        PSDDebug dbg;
        PSDResult R = pointSourceDetection(img, sigma[mCh], po, dump ? &dbg : nullptr);
        out.masks[k - 1] = R.mask;
        F.dRange[mCh] = {img.minval(), img.maxval()};
        if (dump) dump_psd_debug(movieDump, k, dbg, R);

        if (!R.empty) {
            PStruct& P = R.pstruct;
            fill_frame_from_master(F, P, R.isPSF, nCh, mCh);
            std::size_t np = P.size();

            // component size and intensity for each detection
            ConnComp CC = bwconncomp8(R.mask);
            ImageI32 labels = labelmatrix(CC, img.ny(), img.nx());
            std::vector<double> compSize(CC.numObjects), compInt(CC.numObjects);
            for (std::size_t j = 0; j < CC.numObjects; ++j) {
                compSize[j] = static_cast<double>(CC.pixelIdxList[j].size());
                double s = 0.0;
                for (std::size_t p : CC.pixelIdxList[j]) s += img[p];
                compInt[j] = s / static_cast<double>(CC.pixelIdxList[j].size());
            }
            for (std::size_t p = 0; p < np; ++p) {
                const int l = labels(static_cast<std::size_t>(F.y_init[p]) - 1, static_cast<std::size_t>(F.x_init[p]) - 1);
                if (l == 0) throw std::runtime_error("detection outside mask component (MATLAB would index with 0 and error)");
                F.maskN[p] = compSize[static_cast<std::size_t>(l - 1)];
                F.maskA[p] = compInt[static_cast<std::size_t>(l - 1)];
            }

            for (std::size_t ci = 0; ci < nCh; ++ci) {
                if (ci == mCh) continue;
                ImageD simg = load_frame(d, ci, k);
                F.dRange[ci] = {simg.minval(), simg.maxval()};
                std::vector<double> xm = F.x[mCh], ym = F.y[mCh];
                std::vector<double> sig(np, sigma[ci]);
                PStruct S1 = fitGaussians2D(simg, xm, ym, {}, sig, {}, "Ac");
                PStruct S2 = fitGaussians2D(simg, xm, ym, S1.A, sig, S1.c, "xyAc");
                std::vector<char> useLoc(np, 0);
                for (std::size_t p = 0; p < np; ++p) {
                    const double dx = xm[p] - S2.x[p], dy = ym[p] - S2.y[p];
                    const double dist = std::sqrt(dx * dx + dy * dy);
                    useLoc[p] = (dist < 3.0 * sigma[mCh] && S2.A[p] > S1.A[p]) ? 1 : 0;   // NaN comparisons are false
                }
                for (std::size_t p = 0; p < np; ++p) {
                    const PStruct& Q = useLoc[p] ? S2 : S1;
                    F.x[ci][p] = Q.x[p]; F.y[ci][p] = Q.y[p]; F.A[ci][p] = Q.A[p]; F.c[ci][p] = Q.c[p];
                    F.x_pstd[ci][p] = Q.x_pstd[p]; F.y_pstd[ci][p] = Q.y_pstd[p]; F.A_pstd[ci][p] = Q.A_pstd[p]; F.c_pstd[ci][p] = Q.c_pstd[p];
                    F.sigma_r[ci][p] = Q.sigma_r[p]; F.SE_sigma_r[ci][p] = Q.SE_sigma_r[p]; F.RSS[ci][p] = Q.RSS[p]; F.pval_Ar[ci][p] = Q.pval_Ar[p];
                    F.hval_Ar[ci][p] = Q.hval_Ar[p]; F.hval_AD[ci][p] = Q.hval_AD[p];
                }
                F.slaveFixed[ci] = S1; F.slaveLoc[ci] = S2; F.slaveUseLoc[ci] = useLoc;

                // points within slave channel border: remove from all channels
                std::vector<char> rm(np, 0);
                bool any = false;
                for (std::size_t p = 0; p < np; ++p) { rm[p] = std::isnan(S1.x[p]) ? 1 : 0; any = any || rm[p]; }
                F.removedBySlave = rm;
                if (any) {
                    for (std::size_t c = 0; c < nCh; ++c) {
                        remove_cols(F.x[c], rm); remove_cols(F.y[c], rm); remove_cols(F.A[c], rm); remove_cols(F.c[c], rm);
                        remove_cols(F.x_pstd[c], rm); remove_cols(F.y_pstd[c], rm); remove_cols(F.A_pstd[c], rm); remove_cols(F.c_pstd[c], rm);
                        remove_cols(F.sigma_r[c], rm); remove_cols(F.SE_sigma_r[c], rm); remove_cols(F.RSS[c], rm); remove_cols(F.pval_Ar[c], rm);
                        remove_cols(F.hval_Ar[c], rm); remove_cols(F.hval_AD[c], rm); remove_cols(F.isPSF[c], rm);
                    }
                    remove_cols(F.x_init, rm); remove_cols(F.y_init, rm); remove_cols(F.maskA, rm); remove_cols(F.maskN, rm); remove_cols(F.mask_Ar, rm);
                    np = F.x[mCh].size();
                    F.np = np;
                }
                for (std::size_t p = 0; p < np; ++p) F.isPSF[ci][p] = F.hval_AD[ci][p] ? 0 : 1;
            }
        } else {
            F.np = 0;
            for (std::size_t ci = 0; ci < nCh; ++ci) {
                if (ci == mCh) continue;
                ImageD simg = load_frame(d, ci, k);
                F.dRange[ci] = {simg.minval(), simg.maxval()};
            }
        }
        if (dump) dump_frame_info(movieDump, F, d);
    }
    return out;
}

std::vector<DetectionOutput> runDetection(const std::vector<MovieData>& data, const RunOptions& opt, SigmaEstimate* sigmaOut) {
    if (data.empty()) throw std::runtime_error("no data");
    std::vector<double> sigma = opt.sigma;
    SigmaEstimate S;
    if (sigma.empty()) {
        S = estimate_sigma(data, opt);
        sigma = S.sigma;
    } else {
        S.sigma = sigma;
        S.sigmaRaw = sigma;
        for (double& s : sigma) if (s < 1.1) s = 1.1;
        S.sigma = sigma;
    }
    if (sigmaOut) *sigmaOut = S;
    // Movies are independent: parallelise over movies (the per-frame loop gets
    // the threads only when there is a single movie).  Results are stored by
    // index, so output order and content do not depend on the schedule.
    const std::size_t nd = data.size();
    std::vector<DetectionOutput> byIndex(nd);
    const int outerThreads = (nd > 1 && opt.threads > 1) ? std::min<int>(opt.threads, static_cast<int>(nd)) : 1;
    RunOptions inner = opt;
    if (outerThreads > 1) inner.threads = 1;
    #pragma omp parallel for schedule(dynamic) num_threads(outerThreads)
    for (long ii = 0; ii < static_cast<long>(nd); ++ii) {
        const std::size_t i = static_cast<std::size_t>(ii);
        if (!data[i].hasFrames) continue;
        byIndex[i] = run_detection_movie(data[i], sigma, inner, i);
        #pragma omp critical
        std::cout << "Detection done for " << data[i].cellPath << " (" << byIndex[i].frames.size() << " frame(s))" << std::endl;
    }
    std::vector<DetectionOutput> outs;
    for (std::size_t i = 0; i < nd; ++i) {
        if (!data[i].hasFrames) { std::cerr << "Skipping (no frames): " << data[i].cellPath << std::endl; continue; }
        outs.push_back(std::move(byIndex[i]));
    }
    return outs;
}

} // namespace cme
