#include "cme/output.hpp"
#include "cme/tiff_io.hpp"
#include "cme/profile.hpp"
#include <charconv>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <memory>
#include <stdexcept>
#include <string>

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

// Row buffer with a fast "%.17g" formatter.  std::to_chars(general, 17) is
// specified to produce the printf("%.17g") text for finite values; non-finite
// values (MSVC prints e.g. "nan", "-nan(ind)", "inf") go through snprintf so
// the output stays byte-identical to the original fprintf version.
struct RowBuf {
    std::string s;
    void reset() { s.clear(); }
    void tab() { s.push_back('\t'); }
    void str(const std::string& v) { s += v; }
    void uint(std::size_t v) {
        char b[32];
        auto r = std::to_chars(b, b + sizeof b, v);
        s.append(b, r.ptr);
    }
    void integer(int v) {
        char b[32];
        auto r = std::to_chars(b, b + sizeof b, v);
        s.append(b, r.ptr);
    }
    void g17(double v) {
        char b[64];
        if (std::isfinite(v)) {
            auto r = std::to_chars(b, b + sizeof b, v, std::chars_format::general, 17);
            s.append(b, r.ptr);
        } else {
            const int n = std::snprintf(b, sizeof b, "%.17g", v);
            s.append(b, static_cast<std::size_t>(n));
        }
    }
    void flush(FILE* f) {
        if (!s.empty() && std::fwrite(s.data(), 1, s.size(), f) != s.size()) throw std::runtime_error("write failed");
        s.clear();
    }
};

void write_rows(FILE* f, const DetectionOutput& out, const std::vector<std::string>& ch) {
    const std::string mv = movie_name(out.data);
    RowBuf rb;
    rb.s.reserve(1 << 20);
    for (const FrameInfo& F : out.frames) {
        const std::string src = frame_source(out.data, F.frame);
        for (std::size_t p = 0; p < F.np; ++p) {
            rb.str(mv); rb.tab(); rb.uint(F.frame); rb.tab(); rb.uint(p + 1); rb.tab(); rb.str(src);
            for (std::size_t c = 0; c < ch.size(); ++c) {
                rb.tab(); rb.g17(F.x[c][p]); rb.tab(); rb.g17(F.y[c][p]); rb.tab(); rb.g17(F.A[c][p]); rb.tab(); rb.g17(F.c[c][p]);
                rb.tab(); rb.g17(F.x_pstd[c][p]); rb.tab(); rb.g17(F.y_pstd[c][p]); rb.tab(); rb.g17(F.A_pstd[c][p]); rb.tab(); rb.g17(F.c_pstd[c][p]);
                rb.tab(); rb.g17(F.sigma_r[c][p]); rb.tab(); rb.g17(F.SE_sigma_r[c][p]); rb.tab(); rb.g17(F.RSS[c][p]); rb.tab(); rb.g17(F.pval_Ar[c][p]);
                rb.tab(); rb.integer(static_cast<int>(F.hval_Ar[c][p])); rb.tab(); rb.integer(static_cast<int>(F.hval_AD[c][p])); rb.tab(); rb.integer(static_cast<int>(F.isPSF[c][p]));
                rb.tab(); rb.g17(F.s[c]); rb.tab(); rb.g17(F.dRange[c].first); rb.tab(); rb.g17(F.dRange[c].second);
            }
            rb.tab(); rb.g17(F.x_init[p]); rb.tab(); rb.g17(F.y_init[p]); rb.tab(); rb.g17(F.maskA[p]); rb.tab(); rb.g17(F.maskN[p]); rb.tab(); rb.g17(F.mask_Ar[p]);
            rb.s.push_back('\n');
            if (rb.s.size() > (1 << 20)) rb.flush(f);
        }
    }
    rb.flush(f);
}
} // namespace

void write_movie_tsv(const std::string& path, const DetectionOutput& out, const std::vector<std::string>& chNames) {
    prof::Scoped p(prof::TsvWrite);
    FilePtr f = open_w(path);
    write_header(f.get(), chNames);
    write_rows(f.get(), out, chNames);
}

void write_condition_tables(const std::string& dir, const std::vector<DetectionOutput>& outs,
                            const std::vector<std::string>& chNames, const std::vector<double>& sigma) {
    prof::Scoped p(prof::TsvWrite);
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
        RowBuf rb;
        for (const auto& o : outs) {
            const std::string mv = movie_name(o.data);
            for (const FrameInfo& F : o.frames) {
                const std::string src = frame_source(o.data, F.frame);
                for (std::size_t p = 0; p < F.np; ++p) {
                    rb.str(src); rb.tab(); rb.str(mv); rb.tab(); rb.uint(F.frame); rb.tab(); rb.uint(p + 1);
                    for (std::size_t c = 0; c < chNames.size(); ++c) {
                        rb.tab(); rb.g17(F.A[c][p]); rb.tab(); rb.g17(F.c[c][p]); rb.tab(); rb.integer(static_cast<int>(F.hval_Ar[c][p]));
                        rb.tab(); rb.g17(F.x[c][p]); rb.tab(); rb.g17(F.y[c][p]);
                    }
                    rb.s.push_back('\n');
                    if (rb.s.size() > (1 << 20)) rb.flush(f.get());
                }
            }
        }
        rb.flush(f.get());
    }
    {
        FilePtr f = open_w(dir + "/sigma.tsv");
        std::fprintf(f.get(), "channel\tsigma\n");
        for (std::size_t c = 0; c < chNames.size(); ++c) std::fprintf(f.get(), "%s\t%.17g\n", chNames[c].c_str(), sigma[c]);
    }
}

void write_masks(const DetectionOutput& out) {
    prof::Scoped p(prof::MaskWrite);
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
