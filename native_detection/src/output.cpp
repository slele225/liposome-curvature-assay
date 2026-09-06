#include "cme/output.hpp"
#include "cme/tiff_io.hpp"
#include <cstdio>
#include <filesystem>
#include <memory>
#include <stdexcept>

namespace fs = std::filesystem;

namespace cme {

namespace {
struct FileCloser { void operator()(FILE* f) const { if (f) std::fclose(f); } };
using FilePtr = std::unique_ptr<FILE, FileCloser>;
FilePtr open_w(const std::string& path) {
    FILE* f = std::fopen(path.c_str(), "wb");
    if (!f) throw std::runtime_error("cannot write " + path);
    return FilePtr(f);
}

std::string movie_name(const MovieData& d) {
    std::string cp = d.cellPath;
    while (!cp.empty() && (cp.back() == '/' || cp.back() == '\\')) cp.pop_back();
    const auto pos = cp.find_last_of("/\\");
    return pos == std::string::npos ? cp : cp.substr(pos + 1);
}

std::string frame_source(const MovieData& d, std::size_t frame) {
    FrameRef r = frame_ref(d, 0, frame);
    return r.path;
}

void write_header(FILE* f, const std::vector<std::string>& ch) {
    std::fprintf(f, "movie\tframe\tindex\tsource_image");
    const char* dn[] = {"x", "y", "A", "c", "x_pstd", "y_pstd", "A_pstd", "c_pstd", "sigma_r", "SE_sigma_r", "RSS", "pval_Ar", "hval_Ar", "hval_AD", "isPSF", "s", "dRange_min", "dRange_max"};
    for (const auto& c : ch) for (const char* n : dn) std::fprintf(f, "\t%s_%s", n, c.c_str());
    std::fprintf(f, "\tx_init\ty_init\tmaskA\tmaskN\tmask_Ar\n");
}

void write_rows(FILE* f, const DetectionOutput& out, const std::vector<std::string>& ch) {
    const std::string mv = movie_name(out.data);
    for (const FrameInfo& F : out.frames) {
        const std::string src = frame_source(out.data, F.frame);
        for (std::size_t p = 0; p < F.np; ++p) {
            std::fprintf(f, "%s\t%zu\t%zu\t%s", mv.c_str(), F.frame, p + 1, src.c_str());
            for (std::size_t c = 0; c < ch.size(); ++c) {
                std::fprintf(f, "\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%d\t%d\t%d\t%.17g\t%.17g\t%.17g",
                             F.x[c][p], F.y[c][p], F.A[c][p], F.c[c][p], F.x_pstd[c][p], F.y_pstd[c][p], F.A_pstd[c][p], F.c_pstd[c][p],
                             F.sigma_r[c][p], F.SE_sigma_r[c][p], F.RSS[c][p], F.pval_Ar[c][p],
                             static_cast<int>(F.hval_Ar[c][p]), static_cast<int>(F.hval_AD[c][p]), static_cast<int>(F.isPSF[c][p]),
                             F.s[c], F.dRange[c].first, F.dRange[c].second);
            }
            std::fprintf(f, "\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\n", F.x_init[p], F.y_init[p], F.maskA[p], F.maskN[p], F.mask_Ar[p]);
        }
    }
}
} // namespace

void write_movie_tsv(const std::string& path, const DetectionOutput& out, const std::vector<std::string>& chNames) {
    FilePtr f = open_w(path);
    write_header(f.get(), chNames);
    write_rows(f.get(), out, chNames);
}

void write_condition_tables(const std::string& dir, const std::vector<DetectionOutput>& outs,
                            const std::vector<std::string>& chNames, const std::vector<double>& sigma) {
    fs::create_directories(dir);
    {
        FilePtr f = open_w(dir + "/detections_all.tsv");
        write_header(f.get(), chNames);
        for (const auto& o : outs) write_rows(f.get(), o, chNames);
    }
    {
        // Minimal table consumed by the downstream pipeline
        FilePtr f = open_w(dir + "/summary.tsv");
        std::fprintf(f.get(), "source_image\tmovie\tframe\tindex");
        for (const auto& c : chNames) std::fprintf(f.get(), "\tA_%s\tc_%s\thval_%s\tx_%s\ty_%s", c.c_str(), c.c_str(), c.c_str(), c.c_str(), c.c_str());
        std::fprintf(f.get(), "\n");
        for (const auto& o : outs) {
            const std::string mv = movie_name(o.data);
            for (const FrameInfo& F : o.frames) {
                const std::string src = frame_source(o.data, F.frame);
                for (std::size_t p = 0; p < F.np; ++p) {
                    std::fprintf(f.get(), "%s\t%s\t%zu\t%zu", src.c_str(), mv.c_str(), F.frame, p + 1);
                    for (std::size_t c = 0; c < chNames.size(); ++c)
                        std::fprintf(f.get(), "\t%.17g\t%.17g\t%d\t%.17g\t%.17g", F.A[c][p], F.c[c][p], static_cast<int>(F.hval_Ar[c][p]), F.x[c][p], F.y[c][p]);
                    std::fprintf(f.get(), "\n");
                }
            }
        }
    }
    {
        FilePtr f = open_w(dir + "/sigma.tsv");
        std::fprintf(f.get(), "channel\tsigma\n");
        for (std::size_t c = 0; c < chNames.size(); ++c) std::fprintf(f.get(), "%s\t%.17g\n", chNames[c].c_str(), sigma[c]);
    }
}

void write_masks(const DetectionOutput& out) {
    const MovieData& d = out.data;
    const std::string detDir = d.channels[0] + "Detection";
    fs::create_directories(detDir);
    if (d.singleFile) {
        const std::string path = detDir + "/dmasks.tif";
        for (std::size_t k = 0; k < out.masks.size(); ++k) {
            ImageU8 m = out.masks[k];
            for (std::size_t i = 0; i < m.size(); ++i) m[i] = m[i] ? 255 : 0;
            write_tiff_uint8(path, m, k > 0);
        }
    } else {
        fs::create_directories(detDir + "/Masks");
        const int width = static_cast<int>(std::ceil(std::log10(static_cast<double>(d.movieLength) + 1.0)));
        for (std::size_t k = 0; k < out.masks.size(); ++k) {
            char buf[64];
            std::snprintf(buf, sizeof buf, "/Masks/dmask_%0*zu.tif", width, k + 1);
            ImageU8 m = out.masks[k];
            for (std::size_t i = 0; i < m.size(); ++i) m[i] = m[i] ? 255 : 0;
            write_tiff_uint8(detDir + buf, m, false);
        }
    }
}

} // namespace cme
