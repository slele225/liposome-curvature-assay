// Compare fitGaussian2D('xyasc') against MEX results recorded by
// reference_tools/probe_mex_options.m (windows + init + MEX prm).
#include "cme/fit_gaussian2d.hpp"
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

int main(int argc, char** argv) {
    if (argc < 2) { std::printf("usage: fit_probe <mex_probe.txt> [mode]\n"); return 2; }
    std::string mode = argc > 2 ? argv[2] : "xyasc";
    std::ifstream in(argv[1]);
    std::string line;
    int n = 0, exact = 0, close = 0, bothNaN = 0, mism = 0;
    double worst = 0;
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
        cme::ImageD win(nx, nx);
        for (std::size_t i = 0; i < vals.size(); ++i) win[i] = vals[i];
        std::array<double, 5> prm0 = {init[0], init[1], init[2], init[3], init[4]};
        cme::FitResult f = cme::fitGaussian2D(win, prm0, mode);
        ++n;
        bool allExact = true, allClose = true;
        double maxRel = 0;
        for (int k = 0; k < 5; ++k) {
            if (std::isnan(ref[k]) && std::isnan(f.prm[k])) continue;
            if (f.prm[k] != ref[k]) allExact = false;
            const double rel = std::fabs(f.prm[k] - ref[k]) / std::max(std::fabs(ref[k]), 1e-12);
            maxRel = std::max(maxRel, rel);
            if (rel > 1e-6) allClose = false;
        }
        if (allExact) ++exact;
        else if (allClose) ++close;
        else {
            ++mism;
            if (mism <= 15)
                std::printf("cand %d: iters=%d status=%d  cpp [%.6g %.6g %.6g %.6g %.6g]  mex [%.6g %.6g %.6g %.6g %.6g]\n", p, f.iterations, f.status,
                            f.prm[0], f.prm[1], f.prm[2], f.prm[3], f.prm[4], ref[0], ref[1], ref[2], ref[3], ref[4]);
        }
        worst = std::max(worst, maxRel);
    }
    std::printf("%d candidates: %d bit-exact, %d within 1e-6 rel, %d mismatched (worst rel among close %.3g)\n", n, exact, close, mism, worst);
    return 0;
}
