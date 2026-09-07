// cme_detect: native replacement for
//   data = loadConditionData; rng(seed); runDetection(data);
#include "cme/condition_data.hpp"
#include "cme/run_detection.hpp"
#include "cme/output.hpp"
#include "cme/dump.hpp"
#include "cme/profile.hpp"
#include <chrono>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace fs = std::filesystem;

namespace {

void usage() {
    std::cout <<
        "cme_detect --input <condition dir> --channels ch1,ch2[,...] --master ch1 [options]\n"
        "\n"
        "  --input DIR          condition directory (contains the cell*/ movie folders)\n"
        "  --channels a,b,...   channel folder names in the EXACT order loadConditionData received them\n"
        "  --master NAME        must equal the first entry of --channels (master/source channel)\n"
        "  --markers a,b,...    optional fluorophore names (metadata only; do not affect detection)\n"
        "  --seed N             rng(N) seed used by the PSF sigma estimation (default 1)\n"
        "  --sigma s1,s2,...    skip data-driven sigma estimation and use these values (runDetection 'Sigma')\n"
        "  --output DIR         directory for detections_all.tsv / summary.tsv / sigma.tsv (default: <input>/cme_detect_output)\n"
        "  --dump-dir DIR       write all intermediates for regression testing\n"
        "  --no-masks           do not write Detection/dmasks.tif\n"
        "  --no-matlab-layout   do not write Detection/detection_cpp.tsv under the master channel\n"
        "  --movie-selector S   loadConditionData 'MovieSelector' (default 'cell')\n"
        "  --threads N          OpenMP threads (default 1); output does not depend on the thread count\n"
        "  --par-level L        where the threads are used: 'candidate' (default: over candidate fits and\n"
        "                       image columns within each frame) or 'movie' (over movies/frames/images)\n"
        "  --profile            print a per-phase timing breakdown at the end (diagnostic)\n";
}

std::vector<std::string> split(const std::string& s, char sep) {
    std::vector<std::string> out;
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, sep)) if (!item.empty()) out.push_back(item);
    return out;
}

} // namespace

int main(int argc, char** argv) {
    std::string input, master, markersArg, sigmaArg, output, dumpDir, selector = "cell";
    std::vector<std::string> channels;
    unsigned seed = 1;
    bool writeMasks = true, matlabLayout = true;
    int threads = 1;
    cme::ParLevel parLevel = cme::ParLevel::Candidate;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto need = [&](const char* name) -> std::string {
            if (i + 1 >= argc) { std::cerr << "missing value for " << name << "\n"; std::exit(2); }
            return argv[++i];
        };
        if (a == "--input") input = need("--input");
        else if (a == "--channels") channels = split(need("--channels"), ',');
        else if (a == "--master") master = need("--master");
        else if (a == "--markers") markersArg = need("--markers");
        else if (a == "--seed") seed = static_cast<unsigned>(std::stoul(need("--seed")));
        else if (a == "--sigma") sigmaArg = need("--sigma");
        else if (a == "--output") output = need("--output");
        else if (a == "--dump-dir") dumpDir = need("--dump-dir");
        else if (a == "--no-masks") writeMasks = false;
        else if (a == "--no-matlab-layout") matlabLayout = false;
        else if (a == "--movie-selector") selector = need("--movie-selector");
        else if (a == "--threads") threads = std::stoi(need("--threads"));
        else if (a == "--par-level") {
            const std::string v = need("--par-level");
            if (v == "candidate") parLevel = cme::ParLevel::Candidate;
            else if (v == "movie") parLevel = cme::ParLevel::Movie;
            else { std::cerr << "unknown --par-level: " << v << " (expected candidate|movie)\n"; return 2; }
        }
        else if (a == "--profile") cme::prof::enabled = true;
        else if (a == "-h" || a == "--help") { usage(); return 0; }
        else { std::cerr << "unknown argument: " << a << "\n"; usage(); return 2; }
    }
    if (input.empty() || channels.empty() || master.empty()) { usage(); return 2; }
#ifdef _OPENMP
    // Exactly one loop level is parallel (--par-level); the other level runs as
    // a plain loop, so no OpenMP region is ever nested inside another one and
    // the runtime reuses a single thread pool.  Keep the team size fixed.
    omp_set_dynamic(0);
#endif
    if (channels.front() != master) {
        std::cerr << "error: --master (" << master << ") must equal the first entry of --channels (" << channels.front()
                  << "). Channels are not reordered silently; pass them in the order the master should come first.\n";
        return 2;
    }
    for (std::size_t i = 0; i < channels.size(); ++i)
        for (std::size_t j = i + 1; j < channels.size(); ++j)
            if (channels[i] == channels[j]) { std::cerr << "error: duplicate channel name " << channels[i] << "\n"; return 2; }

    std::vector<std::string> markers = markersArg.empty() ? std::vector<std::string>(channels.size(), "unknown") : split(markersArg, ',');
    if (markers.size() != channels.size()) { std::cerr << "error: --markers must have one entry per channel\n"; return 2; }

    try {
        using clock = std::chrono::steady_clock;
        const auto t0 = clock::now();
        cme::ConditionOptions co;
        co.movieSelector = selector;
        std::vector<cme::MovieData> data = cme::loadConditionData(input, channels, markers, co);
        const auto t1 = clock::now();

        cme::RunOptions ro;
        ro.seed = seed;
        ro.dumpDir = dumpDir;
        ro.threads = threads;
        ro.parLevel = parLevel;
        ro.writeMasks = writeMasks;
        if (!sigmaArg.empty()) {
            for (const auto& s : split(sigmaArg, ',')) ro.sigma.push_back(std::stod(s));
            if (ro.sigma.size() != channels.size()) { std::cerr << "error: --sigma needs one value per channel\n"; return 2; }
        }
        cme::SigmaEstimate S;
        std::vector<cme::DetectionOutput> outs;
        // sigma estimation timing is inside runDetection; measure it separately
        const auto t2 = clock::now();
        if (ro.sigma.empty()) {
            S = cme::estimate_sigma(data, ro);
            ro.sigma = S.sigma;
        } else {
            S.sigma = ro.sigma; S.sigmaRaw = ro.sigma;
        }
        const auto t3 = clock::now();
        if (!dumpDir.empty()) cme::dump_sigma_estimate(dumpDir, S);
        outs = cme::runDetection(data, ro, nullptr);
        const auto t4 = clock::now();

        if (output.empty()) {
            std::string in = input;
            while (!in.empty() && (in.back() == '/' || in.back() == '\\')) in.pop_back();
            output = in + "/cme_detect_output";
        }
        cme::write_condition_tables(output, outs, channels, S.sigma);
        // per-movie files are independent of each other: write them in parallel
        std::string outputError;
        #pragma omp parallel for schedule(dynamic) num_threads(threads > 1 ? threads : 1) if(threads > 1)
        for (long oi = 0; oi < static_cast<long>(outs.size()); ++oi) {
            const auto& o = outs[static_cast<std::size_t>(oi)];
            try {
                if (matlabLayout) {
                    const std::string det = o.data.channels[0] + "Detection";
                    fs::create_directories(det);
                    cme::write_movie_tsv(det + "/detection_cpp.tsv", o, channels);
                }
                if (writeMasks) cme::write_masks(o);
            } catch (const std::exception& e) {
                #pragma omp critical
                { if (outputError.empty()) outputError = e.what(); }
            }
        }
        if (!outputError.empty()) throw std::runtime_error(outputError);
        const auto t5 = clock::now();
        auto ms = [](clock::time_point a, clock::time_point b) { return std::chrono::duration<double, std::milli>(b - a).count(); };
        std::printf("Timing (ms): loadConditionData %.1f | sigma estimation %.1f | detection %.1f | output %.1f | total %.1f\n",
                    ms(t0, t1), ms(t2, t3), ms(t3, t4), ms(t4, t5), ms(t0, t5));
        std::size_t total = 0;
        for (const auto& o : outs) for (const auto& F : o.frames) total += F.np;
        std::printf("Detections: %zu rows written to %s\n", total, output.c_str());
        if (cme::prof::enabled) cme::prof::report(stdout, ms(t0, t5) / 1000.0);
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
