#include "cme/dump.hpp"
#include "cme/run_detection.hpp"
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <stdexcept>

namespace fs = std::filesystem;

namespace cme {

namespace {
void fprint_d(FILE* f, double v) { std::fprintf(f, "%.17g", v); }

struct FileCloser { void operator()(FILE* f) const { if (f) std::fclose(f); } };
using FilePtr = std::unique_ptr<FILE, FileCloser>;

FilePtr open_w(const std::string& path) {
    FILE* f = std::fopen(path.c_str(), "wb");
    if (!f) throw std::runtime_error("cannot write " + path);
    return FilePtr(f);
}
} // namespace

std::string dump_movie_name(const MovieData& d, std::size_t movieIndex) {
    std::string cp = d.cellPath;
    while (!cp.empty() && (cp.back() == '/' || cp.back() == '\\')) cp.pop_back();
    const auto pos = cp.find_last_of("/\\");
    std::string name = pos == std::string::npos ? cp : cp.substr(pos + 1);
    char buf[32];
    std::snprintf(buf, sizeof buf, "movie%03zu_", movieIndex + 1);
    return std::string(buf) + name;
}

void dump_image_bin(const std::string& path, const ImageD& img) {
    FilePtr f = open_w(path);
    std::fwrite(img.data(), sizeof(double), img.size(), f.get());
    FilePtr h = open_w(path + ".hdr");
    std::fprintf(h.get(), "float64 %zu %zu column-major\n", img.ny(), img.nx());
}

void dump_image_bin(const std::string& path, const ImageU8& img) {
    FilePtr f = open_w(path);
    std::fwrite(img.data(), 1, img.size(), f.get());
    FilePtr h = open_w(path + ".hdr");
    std::fprintf(h.get(), "uint8 %zu %zu column-major\n", img.ny(), img.nx());
}

void dump_pstruct_tsv(const std::string& path, const PStruct& P) {
    FilePtr f = open_w(path);
    std::fprintf(f.get(), "x\ty\tA\ts\tc\tx_pstd\ty_pstd\tA_pstd\ts_pstd\tc_pstd\tx_init\ty_init\tsigma_r\tSE_sigma_r\tRSS\tpval_Ar\tmask_Ar\thval_Ar\thval_AD\titers\tA2\n");
    for (std::size_t i = 0; i < P.size(); ++i) {
        const double* cols[] = {&P.x[i], &P.y[i], &P.A[i], &P.s[i], &P.c[i], &P.x_pstd[i], &P.y_pstd[i], &P.A_pstd[i], &P.s_pstd[i], &P.c_pstd[i],
                                &P.x_init[i], &P.y_init[i], &P.sigma_r[i], &P.SE_sigma_r[i], &P.RSS[i], &P.pval_Ar[i], &P.mask_Ar[i]};
        for (std::size_t k = 0; k < 17; ++k) { fprint_d(f.get(), *cols[k]); std::fputc('\t', f.get()); }
        std::fprintf(f.get(), "%d\t%d\t%d\t", static_cast<int>(P.hval_Ar[i]), static_cast<int>(P.hval_AD[i]), P.iters[i]);
        fprint_d(f.get(), P.A2[i]);
        std::fputc('\n', f.get());
    }
}

void dump_psd_debug(const std::string& dir, std::size_t frame, const PSDDebug& dbg, const PSDResult& res) {
    char pre[64];
    std::snprintf(pre, sizeof pre, "/frame%04zu_", frame);
    const std::string p = dir + pre;
    if (!dbg.imgLoG.empty()) {
        dump_image_bin(p + "imgLoG.bin", dbg.imgLoG);
        dump_image_bin(p + "A_est.bin", dbg.A_est);
        dump_image_bin(p + "c_est.bin", dbg.c_est);
        dump_image_bin(p + "pval_prefilter.bin", dbg.pvalPrefilter);
        dump_image_bin(p + "mask_prefilter.bin", dbg.maskPrefilter);
        dump_image_bin(p + "mask_combined.bin", dbg.maskCombined);
    }
    dump_image_bin(p + "mask_final.bin", res.mask);
    {
        FilePtr f = open_w(p + "lm.tsv");
        std::fprintf(f.get(), "lmx\tlmy\n");
        for (std::size_t i = 0; i < dbg.lmx.size(); ++i) std::fprintf(f.get(), "%.17g\t%.17g\n", dbg.lmx[i], dbg.lmy[i]);
    }
    {
        FilePtr f = open_w(p + "logThreshold.txt");
        std::fprintf(f.get(), "%.17g\n", dbg.logThreshold);
    }
    dump_pstruct_tsv(p + "master_fitall.tsv", dbg.fitAll);
    {
        FilePtr f = open_w(p + "master_keep.tsv");
        std::fprintf(f.get(), "keepNonNaN\n");
        for (char v : dbg.keepNonNaN) std::fprintf(f.get(), "%d\n", static_cast<int>(v));
        FilePtr g = open_w(p + "master_keepfinal.tsv");
        std::fprintf(g.get(), "keepFinal\n");
        for (char v : dbg.keepFinal) std::fprintf(g.get(), "%d\n", static_cast<int>(v));
    }
    dump_pstruct_tsv(p + "master_final.tsv", res.pstruct);
}

void dump_frame_info(const std::string& dir, const FrameInfo& F, const MovieData& d) {
    char pre[64];
    std::snprintf(pre, sizeof pre, "/frame%04zu_", F.frame);
    const std::string p = dir + pre;
    const std::size_t nCh = F.nCh;
    FilePtr f = open_w(p + "frameinfo.tsv");
    // header
    std::fprintf(f.get(), "idx");
    const char* dn[] = {"x", "y", "A", "c", "x_pstd", "y_pstd", "A_pstd", "c_pstd", "sigma_r", "SE_sigma_r", "RSS", "pval_Ar"};
    const char* ln[] = {"hval_Ar", "hval_AD", "isPSF"};
    for (std::size_t c = 0; c < nCh; ++c) {
        for (const char* n : dn) std::fprintf(f.get(), "\t%s_%zu", n, c + 1);
        for (const char* n : ln) std::fprintf(f.get(), "\t%s_%zu", n, c + 1);
    }
    std::fprintf(f.get(), "\tx_init\ty_init\tmaskA\tmaskN\tmask_Ar\n");
    for (std::size_t p = 0; p < F.np; ++p) {
        std::fprintf(f.get(), "%zu", p + 1);
        for (std::size_t c = 0; c < nCh; ++c) {
            const std::vector<std::vector<double>>* dv[] = {&F.x, &F.y, &F.A, &F.c, &F.x_pstd, &F.y_pstd, &F.A_pstd, &F.c_pstd, &F.sigma_r, &F.SE_sigma_r, &F.RSS, &F.pval_Ar};
            for (auto* v : dv) { std::fputc('\t', f.get()); fprint_d(f.get(), (*v)[c][p]); }
            std::fprintf(f.get(), "\t%d\t%d\t%d", static_cast<int>(F.hval_Ar[c][p]), static_cast<int>(F.hval_AD[c][p]), static_cast<int>(F.isPSF[c][p]));
        }
        std::fputc('\t', f.get()); fprint_d(f.get(), F.x_init[p]);
        std::fputc('\t', f.get()); fprint_d(f.get(), F.y_init[p]);
        std::fputc('\t', f.get()); fprint_d(f.get(), F.maskA[p]);
        std::fputc('\t', f.get()); fprint_d(f.get(), F.maskN[p]);
        std::fputc('\t', f.get()); fprint_d(f.get(), F.mask_Ar[p]);
        std::fputc('\n', f.get());
    }
    // slave fits
    for (std::size_t c = 0; c < nCh; ++c) {
        if (F.slaveFixed.size() <= c || F.slaveFixed[c].size() == 0) continue;
        char sb[64];
        std::snprintf(sb, sizeof sb, "slave_ch%zu_", c + 1);
        dump_pstruct_tsv(p + sb + "fixed.tsv", F.slaveFixed[c]);
        dump_pstruct_tsv(p + sb + "loc.tsv", F.slaveLoc[c]);
        FilePtr g = open_w(p + sb + "useloc.tsv");
        std::fprintf(g.get(), "useLoc\n");
        for (char v : F.slaveUseLoc[c]) std::fprintf(g.get(), "%d\n", static_cast<int>(v));
    }
    {
        FilePtr g = open_w(p + "dRange.tsv");
        std::fprintf(g.get(), "channel\tmin\tmax\n");
        for (std::size_t c = 0; c < nCh; ++c) std::fprintf(g.get(), "%zu\t%.17g\t%.17g\n", c + 1, F.dRange[c].first, F.dRange[c].second);
    }
    (void)d;
}

void dump_sigma_estimate(const std::string& dir, const SigmaEstimate& S) {
    fs::create_directories(dir);
    {
        FilePtr f = open_w(dir + "/sigma.tsv");
        std::fprintf(f.get(), "channel\tsigma_raw\tsigma\n");
        for (std::size_t c = 0; c < S.sigma.size(); ++c) std::fprintf(f.get(), "%zu\t%.17g\t%.17g\n", c + 1, S.sigmaRaw[c], S.sigma[c]);
    }
    for (std::size_t c = 0; c < S.debug.size(); ++c) {
        const PsfSigmaDebug& D = S.debug[c];
        char pre[64];
        std::snprintf(pre, sizeof pre, "/sigma_ch%zu_", c + 1);
        const std::string p = dir + pre;
        {
            FilePtr f = open_w(p + "svect.tsv");
            std::fprintf(f.get(), "s\n");
            for (double v : D.svect) std::fprintf(f.get(), "%.17g\n", v);
        }
        {
            FilePtr f = open_w(p + "svect_per_image.tsv");
            std::fprintf(f.get(), "image\tcount\n");
            for (std::size_t i = 0; i < D.svectPerImage.size(); ++i) std::fprintf(f.get(), "%zu\t%zu\n", i + 1, D.svectPerImage[i].size());
        }
        for (std::size_t i = 0; i < D.refitPerImage.size(); ++i) {
            char sb[64];
            std::snprintf(sb, sizeof sb, "refit_img%03zu.tsv", i + 1);
            dump_pstruct_tsv(p + sb, D.refitPerImage[i]);
        }
        {
            FilePtr f = open_w(p + "gmm.tsv");
            std::fprintf(f.get(), "k\tcomponent\tmu\tSigma\tPComponents\tNlogL\tBIC\titers\tconverged\tinitIdx\tinitMu\n");
            for (const auto& g : D.gmm) {
                for (int j = 0; j < g.k; ++j) {
                    std::fprintf(f.get(), "%d\t%d\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%d\t%d\t%ld\t%.17g\n", g.k, j + 1, g.mu[j], g.Sigma[j], g.PComponents[j],
                                 g.NlogL, g.BIC, g.iters, g.converged ? 1 : 0, g.initIdx[j], g.initMu[j]);
                }
            }
            std::fprintf(f.get(), "# gmmFailed=%d chosenK=%d chosenComponent=%d\n", D.gmmFailed ? 1 : 0, D.chosenK, D.chosenComponent);
        }
    }
}

} // namespace cme
