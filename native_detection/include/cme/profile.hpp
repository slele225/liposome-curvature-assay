// Lightweight phase profiler (enabled with `cme_detect --profile`).
//
// Every thread accumulates wall-clock time per phase in its own thread_local
// table (no synchronisation inside the hot loops); the tables are summed when
// the report is printed.  Overhead per scoped timer is one steady_clock read
// on entry and exit (~30 ns), so timers are placed around whole loops or
// around per-candidate work (>= tens of microseconds), never per pixel.
//
// When profiling is disabled the timers compile to a branch on a global bool.
#pragma once
#include <chrono>
#include <cstdint>
#include <cstdio>

namespace cme {
namespace prof {

enum Phase : int {
    TiffLoad = 0,
    PadConv,          // padarrayXT + the five separable convolutions
    LoGArith,         // LoG / A_est / c_est element-wise arithmetic
    Prefilter,        // per-pixel t-test (tcdf) -> prefilter mask
    LocMax,           // locmax2d
    ConnComp,         // bwconncomp / labelmatrix / bwlabel
    FitPrep,          // fitGaussians2D per-candidate window extraction, mask_Ar
    LMSetup,          // fitGaussian2D: DataStruct + solver allocation + set
    LMIterate,        // fitGaussian2D: LM iterations
    LMPost,           // fitGaussian2D: covariance, residual statistics
    ADTest,           // adtest_mex
    PvalT,            // fitGaussians2D final tcdf loop
    Redundancy,       // removeRedundant ball query
    Gmm,              // gmdistribution_fit_1d
    MaskWrite,
    TsvWrite,
    Other,
    // LM internals (subsets of LMIterate / LMSetup, reported separately)
    LM_f,             // model evaluation callback
    LM_df,            // Jacobian evaluation callback
    LM_qrpt,          // QRPT decomposition of J
    LM_qtvec,         // Q^T f
    NumPhases
};

// Which caller the LM fit is serving (accumulated separately).
enum Context : int {
    CtxMaster = 0,    // pointSourceDetection during runDetection (master channel)
    CtxSigmaDetect,   // pointSourceDetection(img, 1.5, 'xyac') in sigma estimation
    CtxSigmaRefit,    // 'xyasc' refit in sigma estimation
    CtxSlaveFixed,    // slave channel 'Ac'
    CtxSlaveLoc,      // slave channel 'xyAc'
    NumContexts
};

extern bool enabled;
extern thread_local int context;

const char* phase_name(int p);
const char* context_name(int c);

void add(int phase, double seconds);
void add_fit(int context, double seconds, long iterations, long fevals, long dfevals);
void add_count(int phase, long n);

// Sum of all threads' tables; printed as a table to `out`.
void report(std::FILE* out, double wallTotal);

class Scoped {
public:
    explicit Scoped(int phase) : phase_(phase), active_(enabled) {
        if (active_) t0_ = std::chrono::steady_clock::now();
    }
    ~Scoped() { stop(); }
    void stop() {
        if (!active_) return;
        active_ = false;
        const auto t1 = std::chrono::steady_clock::now();
        add(phase_, std::chrono::duration<double>(t1 - t0_).count());
    }
private:
    int phase_;
    bool active_;
    std::chrono::steady_clock::time_point t0_;
};

struct ContextGuard {
    int prev;
    explicit ContextGuard(int c) : prev(context) { context = c; }
    ~ContextGuard() { context = prev; }
};

} // namespace prof
} // namespace cme

// C-linkage hooks used by the C solver wrapper (third_party/gsl116/g116_lm.c).
extern "C" {
int cme_prof_enabled_c(void);
double cme_prof_now_c(void);
void cme_prof_add_c(int phase, double t0);
}
