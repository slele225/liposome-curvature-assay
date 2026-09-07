#include "cme/profile.hpp"
#include <mutex>
#include <vector>
#include <cstring>

namespace cme {
namespace prof {

bool enabled = false;
thread_local int context = CtxMaster;

// third_party/gsl116/g116_lm.c refers to these phases by number
static_assert(LM_qrpt == 19, "update CME_PROF_LM_QRPT in g116_lm.c");
static_assert(LM_qtvec == 20, "update CME_PROF_LM_QTVEC in g116_lm.c");

namespace {

struct Table {
    double phase[NumPhases] = {};
    long count[NumPhases] = {};
    double fitTime[NumContexts] = {};
    long fitCalls[NumContexts] = {};
    long fitIters[NumContexts] = {};
    long fitFevals[NumContexts] = {};
    long fitDfevals[NumContexts] = {};
};

std::mutex g_mutex;
std::vector<Table*> g_tables;

Table& local_table() {
    thread_local Table* t = nullptr;
    if (!t) {
        t = new Table();
        std::lock_guard<std::mutex> lock(g_mutex);
        g_tables.push_back(t);
    }
    return *t;
}

const char* const kPhaseNames[NumPhases] = {
    "TIFF load", "pad + separable convolutions", "LoG / A_est / c_est arithmetic",
    "prefilter t-test (tcdf per pixel)", "locmax2d", "connected components / labels",
    "fit: window extraction, mask_Ar", "fit: LM setup (alloc, set)", "fit: LM iterations",
    "fit: covariance, residual stats", "fit: Anderson-Darling test", "fit: pval_Ar tcdf",
    "redundancy (ball query)", "GMM (sigma estimation)", "mask writing", "TSV writing", "other",
    "  LM: model f evaluation", "  LM: Jacobian evaluation", "  LM: QRPT decomposition", "  LM: Q^T f",
};
constexpr int kFirstSubPhase = LM_f;   // sub-phases are excluded from the total
const char* const kContextNames[NumContexts] = {
    "master detection", "sigma: detection pass (xyac)", "sigma: free-sigma refit (xyasc)",
    "slave: fixed-position (Ac)", "slave: localized (xyAc)",
};

} // namespace

const char* phase_name(int p) { return (p >= 0 && p < NumPhases) ? kPhaseNames[p] : "?"; }
const char* context_name(int c) { return (c >= 0 && c < NumContexts) ? kContextNames[c] : "?"; }

void add(int phase, double seconds) {
    Table& t = local_table();
    t.phase[phase] += seconds;
    t.count[phase] += 1;
}

void add_count(int phase, long n) {
    local_table().count[phase] += n;
}

void add_fit(int ctx, double seconds, long iterations, long fevals, long dfevals) {
    Table& t = local_table();
    t.fitTime[ctx] += seconds;
    t.fitCalls[ctx] += 1;
    t.fitIters[ctx] += iterations;
    t.fitFevals[ctx] += fevals;
    t.fitDfevals[ctx] += dfevals;
}

void report(std::FILE* out, double wallTotal) {
    Table sum;
    {
        std::lock_guard<std::mutex> lock(g_mutex);
        for (const Table* t : g_tables) {
            for (int p = 0; p < NumPhases; ++p) { sum.phase[p] += t->phase[p]; sum.count[p] += t->count[p]; }
            for (int c = 0; c < NumContexts; ++c) {
                sum.fitTime[c] += t->fitTime[c]; sum.fitCalls[c] += t->fitCalls[c];
                sum.fitIters[c] += t->fitIters[c]; sum.fitFevals[c] += t->fitFevals[c]; sum.fitDfevals[c] += t->fitDfevals[c];
            }
        }
    }
    double total = 0.0;
    for (int p = 0; p < kFirstSubPhase; ++p) total += sum.phase[p];
    std::fprintf(out, "\n=== profile: thread-time by phase (sum over threads; %zu thread table(s)) ===\n", g_tables.size());
    std::fprintf(out, "%-40s %12s %8s %12s\n", "phase", "seconds", "%", "calls");
    for (int p = 0; p < NumPhases; ++p) {
        if (sum.count[p] == 0 && sum.phase[p] == 0.0) continue;
        if (p == kFirstSubPhase) std::fprintf(out, "%-40s\n", "(sub-phases of the LM fit, not added to the total)");
        std::fprintf(out, "%-40s %12.3f %7.1f%% %12ld\n", kPhaseNames[p], sum.phase[p], total > 0 ? 100.0 * sum.phase[p] / total : 0.0, sum.count[p]);
    }
    std::fprintf(out, "%-40s %12.3f %7.1f%%\n", "total instrumented thread-time", total, 100.0);
    std::fprintf(out, "%-40s %12.3f\n", "wall-clock (process)", wallTotal);
    std::fprintf(out, "\n=== profile: Gaussian fits by caller (LM setup+iterate+post+AD) ===\n");
    std::fprintf(out, "%-36s %10s %10s %12s %12s %12s %10s\n", "context", "seconds", "fits", "LM iters", "f evals", "df evals", "us/fit");
    for (int c = 0; c < NumContexts; ++c) {
        if (sum.fitCalls[c] == 0) continue;
        std::fprintf(out, "%-36s %10.3f %10ld %12ld %12ld %12ld %10.1f\n", kContextNames[c], sum.fitTime[c], sum.fitCalls[c],
                     sum.fitIters[c], sum.fitFevals[c], sum.fitDfevals[c], 1e6 * sum.fitTime[c] / static_cast<double>(sum.fitCalls[c]));
    }
    std::fflush(out);
}

} // namespace prof
} // namespace cme

extern "C" int cme_prof_enabled_c(void) { return cme::prof::enabled ? 1 : 0; }
extern "C" double cme_prof_now_c(void) {
    return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}
extern "C" void cme_prof_add_c(int phase, double t0) {
    cme::prof::add(phase, cme_prof_now_c() - t0);
}
