// In-process micro-benchmark of fitGaussian2D on recorded real candidate
// windows (the 'C' lines of reference_tools/probe_mex_options.m output).
//
//   fit_bench <mex_probe.txt> [repetitions] [results.txt] [mode ...]
//
// Every window is fitted `repetitions` times in each mode; the per-fit and
// per-LM-iteration times are printed, and (optionally) all fit results are
// written to results.txt with %.17g so that two builds can be compared for
// bit-identical behaviour (diff the two files).
#include "cme/fit_gaussian2d.hpp"
#include "cme/profile.hpp"
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

struct Window {
    cme::ImageD img;
    std::array<double, 5> prm0;
};

int main(int argc, char** argv) {
    if (argc < 2) { std::printf("usage: fit_bench <mex_probe.txt> [repetitions] [results.txt] [mode ...]\n"); return 2; }
    const int reps = argc > 2 ? std::atoi(argv[2]) : 20;
    const std::string resultsPath = argc > 3 ? argv[3] : "";
    std::vector<std::string> modes;
    for (int i = 4; i < argc; ++i) modes.push_back(argv[i]);
    if (modes.empty()) modes = {"xyasc", "xyAc", "xyac", "Ac"};

    std::ifstream in(argv[1]);
    std::string line;
    std::vector<Window> wins;
    while (std::getline(in, line)) {
        if (line.rfind("C ", 0) != 0) continue;
        std::istringstream ss(line);
        std::string tag, w;
        int p;
        double init[5], ref[5];
        ss >> tag >> p >> w;
        for (double& v : init) ss >> v;
        ss >> w;
        for (double& v : ref) ss >> v;
        std::getline(in, line);
        std::istringstream d(line);
        std::vector<double> vals;
        std::string tokv;
        while (d >> tokv) vals.push_back(std::strtod(tokv.c_str(), nullptr));
        const int nx = static_cast<int>(std::lround(std::sqrt(static_cast<double>(vals.size()))));
        Window W;
        W.img = cme::ImageD(nx, nx);
        for (std::size_t i = 0; i < vals.size(); ++i) W.img[i] = vals[i];
        W.prm0 = {init[0], init[1], init[2], init[3], init[4]};
        wins.push_back(std::move(W));
    }
    if (wins.empty()) { std::printf("no windows found\n"); return 1; }
    std::printf("%zu windows (%zux%zu), %d repetitions\n", wins.size(), wins[0].img.ny(), wins[0].img.nx(), reps);
    if (std::getenv("CME_PROFILE")) cme::prof::enabled = true;   // phase breakdown (adds timer overhead)

    std::FILE* res = resultsPath.empty() ? nullptr : std::fopen(resultsPath.c_str(), "w");
    using clock = std::chrono::steady_clock;
    for (const std::string& mode : modes) {
        long iters = 0;
        double checksum = 0.0;
        const auto t0 = clock::now();
        for (int r = 0; r < reps; ++r) {
            for (std::size_t i = 0; i < wins.size(); ++i) {
                cme::FitResult f = cme::fitGaussian2D(wins[i].img, wins[i].prm0, mode);
                iters += f.iterations;
                checksum += f.prm[2];
                if (res && r == 0) {
                    std::fprintf(res, "%s %zu iters=%d status=%d hAD=%d", mode.c_str(), i, f.iterations, f.status, f.hAD ? 1 : 0);
                    for (double v : f.prm) std::fprintf(res, " %.17g", v);
                    for (double v : f.prmStd) std::fprintf(res, " %.17g", v);
                    std::fprintf(res, " %.17g %.17g %.17g %.17g\n", f.RSS, f.mean, f.std, f.A2);
                }
            }
        }
        const double secs = std::chrono::duration<double>(clock::now() - t0).count();
        const double nfits = static_cast<double>(reps) * static_cast<double>(wins.size());
        std::printf("mode %-6s  %8.1f us/fit  %7.2f us/LM-iteration  (%.1f iterations/fit, checksum %.6g)\n",
                    mode.c_str(), 1e6 * secs / nfits, 1e6 * secs / static_cast<double>(iters),
                    static_cast<double>(iters) / nfits, checksum);
    }
    if (res) std::fclose(res);
    if (cme::prof::enabled) cme::prof::report(stdout, 0.0);
    return 0;
}
