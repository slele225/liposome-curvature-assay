// Phase timing of the sigma-estimation work on one image:
//   pointSourceDetection(img, 1.5, 'xyac') and the 'xyasc' refit.
#include "cme/tiff_io.hpp"
#include "cme/point_source_detection.hpp"
#include "cme/fit_gaussians2d.hpp"
#include <chrono>
#include <cstdio>
#include <map>

int main(int argc, char** argv) {
    if (argc < 2) { std::printf("usage: profile_sigma <tiff> [sigma]\n"); return 2; }
    const double sigma = argc > 2 ? std::atof(argv[2]) : 1.5;
    using clock = std::chrono::steady_clock;
    auto ms = [](clock::time_point a, clock::time_point b) { return std::chrono::duration<double, std::milli>(b - a).count(); };
    auto t0 = clock::now();
    cme::ImageD img = cme::read_tiff_frame_double(argv[1], 1);
    auto t1 = clock::now();
    cme::PSDOptions o; o.mode = "xyac";
    cme::PSDDebug dbg;
    cme::PSDResult r = cme::pointSourceDetection(img, sigma, o, &dbg);
    auto t2 = clock::now();
    std::size_t np = r.empty ? 0 : r.pstruct.size();
    std::vector<double> sig(np, sigma);
    cme::PStruct P;
    if (np) P = cme::fitGaussians2D(img, r.pstruct.x, r.pstruct.y, r.pstruct.A, sig, r.pstruct.c, "xyasc");
    auto t3 = clock::now();
    std::map<int, int> hist;
    long itersTotal = 0;
    for (std::size_t i = 0; i < P.size(); ++i) { itersTotal += P.iters[i]; hist[P.iters[i] >= 500 ? 500 : (P.iters[i] / 50) * 50]++; }
    std::printf("read %.0f ms | psd(1.5) %.0f ms (candidates %zu, kept %zu) | refit xyasc %.0f ms (total LM iterations %ld)\n",
                ms(t0, t1), ms(t1, t2), dbg.lmx.size(), np, ms(t2, t3), itersTotal);
    for (auto& kv : hist) std::printf("  iters bucket %3d: %d fits\n", kv.first, kv.second);
    return 0;
}
