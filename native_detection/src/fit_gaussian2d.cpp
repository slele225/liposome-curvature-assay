#include "cme/fit_gaussian2d.hpp"
#include "cme/stats.hpp"
#include "cme/profile.hpp"
#include "g116_lm.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>
#include <cctype>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <stdexcept>

// The Levenberg-Marquardt solver is the GSL 1.16 'lmsder' implementation
// (third_party/gsl116), i.e. the GSL generation the original MEX was linked
// against.  GSL 2.x changed the column scaling of lmsder, which changes the
// iteration path of ill-conditioned fits (see PORTING_NOTES.md §2.4).

namespace cme {

namespace {

// GSL's default error handler aborts the process.  The solver reports
// non-convergence through return codes, never through the handler, and the
// handler is a process-wide setting (setting/restoring it per fit from several
// threads would race), so it is switched off once for the whole process.
[[maybe_unused]] gsl_error_handler_t* const kGslHandlerOff = gsl_set_error_handler_off();

// The MEX evaluates both the model and the Jacobian with |sigma| (verified on
// real candidate windows: with the signed derivative 10/356 free-sigma fits
// diverge from the MEX, with |sigma| every fit lands in the MEX's basin).

struct DataStruct {
    int nx = 0;                 // window side length
    int np = 0;                 // number of estimated parameters
    const double* pixels = nullptr;
    double* gx = nullptr;       // 1-D Gaussian factors (nx each)
    double* gy = nullptr;
    const int* idx = nullptr;   // linear indices of valid (non-NaN) pixels
    const int* col = nullptr;   // idx / nx, precomputed once per fit
    const int* row = nullptr;   // idx % nx
    int nValid = 0;
    int estIdx[5] = {};         // which of [xp yp A s c] are estimated
    double prmVect[5];
    long nF = 0, nDF = 0;       // evaluation counters (profiling)
};

// Per-thread reusable workspace: the LM solver state for the current
// (n, p), the covariance matrix and the index/factor buffers.  Reusing the
// solver is exact: g116_solver_set() re-initialises every field the
// iteration reads (x, f, J, dx, r, tau, perm, diag, delta, par, fnorm,
// rptdx, w, f_trial); the remaining vectors are fully overwritten before
// they are read in every iteration.
struct Workspace {
    g116_solver* solver = nullptr;
    std::size_t n = 0, p = 0;
    gsl_matrix* covar = nullptr;
    std::size_t covp = 0;
    std::vector<int> idx, col, row;
    std::vector<double> gx, gy, res, xInit;
    ~Workspace() {
        if (solver) g116_solver_free(solver);
        if (covar) gsl_matrix_free(covar);
    }
    g116_solver* get_solver(std::size_t nn, std::size_t pp) {
        if (!solver || n != nn || p != pp) {
            if (solver) g116_solver_free(solver);
            solver = g116_solver_alloc(nn, pp);
            n = nn; p = pp;
        }
        return solver;
    }
    gsl_matrix* get_covar(std::size_t pp) {
        if (!covar || covp != pp) {
            if (covar) gsl_matrix_free(covar);
            covar = gsl_matrix_alloc(pp, pp);
            covp = pp;
        }
        return covar;
    }
};
thread_local Workspace tls;

// Evaluate 1-D factors for the current parameters.
inline void update_factors(DataStruct* d, const gsl_vector* x) {
    for (int i = 0; i < d->np; ++i) d->prmVect[d->estIdx[i]] = gsl_vector_get(x, i);
    const int nx = d->nx;
    const int b = nx / 2;
    const double xp = d->prmVect[0];
    const double yp = d->prmVect[1];
    const double sigma = std::fabs(d->prmVect[3]);
    const double dd = 2.0 * sigma * sigma;
    for (int i = 0; i < nx; ++i) {
        const double k = static_cast<double>(i - b);
        const double xi = k - xp;
        d->gx[i] = std::exp(-xi * xi / dd);
        const double yi = k - yp;
        d->gy[i] = std::exp(-yi * yi / dd);
    }
}

int gaussian_f(const gsl_vector* x, void* params, gsl_vector* f) {
    DataStruct* d = static_cast<DataStruct*>(params);
    ++d->nF;
    prof::Scoped pf(prof::LM_f);
    update_factors(d, x);
    const double A = d->prmVect[2];
    const double c = d->prmVect[4];
    const double* gx = d->gx;
    const double* gy = d->gy;
    for (int i = 0; i < d->nValid; ++i) {
        const int idx = d->idx[i];
        // column-major: column = x index, row = y index
        const int col = d->col[i];
        const int row = d->row[i];
        gsl_vector_set(f, i, A * gx[col] * gy[row] + c - d->pixels[idx]);
    }
    return GSL_SUCCESS;
}

int gaussian_df(const gsl_vector* x, void* params, gsl_matrix* J) {
    DataStruct* d = static_cast<DataStruct*>(params);
    ++d->nDF;
    prof::Scoped pdf(prof::LM_df);
    update_factors(d, x);
    const int nx = d->nx;
    const int b = nx / 2;
    const double xp = d->prmVect[0];
    const double yp = d->prmVect[1];
    const double A = d->prmVect[2];
    const double sigma = std::fabs(d->prmVect[3]);
    const double sigma2 = sigma * sigma;
    const double sigma3 = sigma2 * sigma;
    const double* gx = d->gx;
    const double* gy = d->gy;
    const int np = d->np;
    for (int i = 0; i < d->nValid; ++i) {
        const int col = d->col[i];
        const int row = d->row[i];
        const double xi = static_cast<double>(col - b) - xp;
        const double yi = static_cast<double>(row - b) - yp;
        const double g = gx[col] * gy[row];
        for (int k = 0; k < np; ++k) {
            double v;
            switch (d->estIdx[k]) {
                case 0: v = xi * A * g / sigma2; break;
                case 1: v = yi * A * g / sigma2; break;
                case 2: v = g; break;
                case 3: v = (xi * xi + yi * yi) * A * g / sigma3; break;
                default: v = 1.0; break;
            }
            gsl_matrix_set(J, i, k, v);
        }
    }
    return GSL_SUCCESS;
}

int gaussian_fdf(const gsl_vector* x, void* params, gsl_vector* f, gsl_matrix* J) {
    gaussian_f(x, params, f);
    gaussian_df(x, params, J);
    return GSL_SUCCESS;
}

} // namespace

FitResult fitGaussian2D(const ImageD& window, const std::array<double, 5>& prm0, const std::string& modeIn,
                        const FitOptions& opt) {
    if (window.nx() != window.ny()) throw std::invalid_argument("fitGaussian2D: window must be square");
    const int nx = static_cast<int>(window.nx());
    const int N = nx * nx;

    const auto tFit0 = prof::enabled ? std::chrono::steady_clock::now() : std::chrono::steady_clock::time_point();
    prof::Scoped pSetup(prof::LMSetup);
    Workspace& ws = tls;
    DataStruct d;
    d.nx = nx;
    d.pixels = window.data();
    ws.gx.assign(static_cast<std::size_t>(nx), 0.0);
    ws.gy.assign(static_cast<std::size_t>(nx), 0.0);
    d.gx = ws.gx.data();
    d.gy = ws.gy.data();
    ws.idx.clear(); ws.col.clear(); ws.row.clear();
    for (int i = 0; i < N; ++i) {
        if (!std::isnan(d.pixels[i])) { ws.idx.push_back(i); ws.col.push_back(i / nx); ws.row.push_back(i % nx); }
    }
    d.idx = ws.idx.data();
    d.col = ws.col.data();
    d.row = ws.row.data();
    d.nValid = static_cast<int>(ws.idx.size());
    for (int i = 0; i < 5; ++i) d.prmVect[i] = prm0[i];

    // mode -> estIdx, canonical order x y a s c
    std::string mode = modeIn;
    for (char& ch : mode) ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
    static const char refMode[] = "xyasc";
    d.np = 0;
    for (int i = 0; i < 5; ++i) {
        if (mode.find(refMode[i]) != std::string::npos) d.estIdx[d.np++] = i;
    }
    if (d.np == 0) throw std::invalid_argument("fitGaussian2D: no parameters to estimate");
    if (d.nValid <= d.np) throw std::invalid_argument("fitGaussian2D: not enough valid pixels");

    ws.xInit.resize(static_cast<std::size_t>(d.np));
    for (int i = 0; i < d.np; ++i) ws.xInit[static_cast<std::size_t>(i)] = d.prmVect[d.estIdx[i]];

    g116_solver* s = ws.get_solver(static_cast<std::size_t>(d.nValid), static_cast<std::size_t>(d.np));
    if (!s) throw std::runtime_error("fitGaussian2D: solver allocation failed");
    g116_function_fdf f;
    f.f = &gaussian_f;
    f.df = &gaussian_df;
    f.fdf = &gaussian_fdf;
    f.n = static_cast<size_t>(d.nValid);
    f.p = static_cast<size_t>(d.np);
    f.params = &d;
    gsl_vector_view xv = gsl_vector_view_array(ws.xInit.data(), static_cast<std::size_t>(d.np));
    g116_solver_set(s, &f, &xv.vector);
    pSetup.stop();

    prof::Scoped pIter(prof::LMIterate);
    int status = GSL_CONTINUE;
    int iter = 0;
    do {
        ++iter;
        status = g116_solver_iterate(s, &f);
        if (status) break;
        status = g116_test_delta(g116_solver_dx(s), g116_solver_x(s), opt.eAbs, opt.eRel);
    } while (status == GSL_CONTINUE && iter < opt.maxIter);
    pIter.stop();

    prof::Scoped pPost(prof::LMPost);
    FitResult r;
    r.iterations = iter;
    r.status = status;
    r.nValid = d.nValid;
    r.estIdx.assign(d.estIdx, d.estIdx + d.np);
    const gsl_vector* xs = g116_solver_x(s);
    for (int i = 0; i < d.np; ++i) d.prmVect[d.estIdx[i]] = gsl_vector_get(xs, i);
    d.prmVect[3] = std::fabs(d.prmVect[3]);
    for (int i = 0; i < 5; ++i) r.prm[i] = d.prmVect[i];

    // residuals (as stored in the solver at the final iterate)
    const gsl_vector* fs = g116_solver_f(s);
    std::vector<double>& res = ws.res;
    res.resize(static_cast<std::size_t>(d.nValid));
    double RSS = 0.0;
    for (int i = 0; i < d.nValid; ++i) {
        res[static_cast<std::size_t>(i)] = gsl_vector_get(fs, i);
        RSS += res[static_cast<std::size_t>(i)] * res[static_cast<std::size_t>(i)];
    }
    r.RSS = RSS;

    // covariance from the Jacobian at the final iterate
    gsl_matrix* covar = ws.get_covar(static_cast<std::size_t>(d.np));
    g116_covar(g116_solver_J(s), 0.0, covar);
    const double iRSS = RSS / static_cast<double>(d.nValid - d.np - 1);
    r.prmStd.resize(static_cast<std::size_t>(d.np));
    for (int i = 0; i < d.np; ++i) r.prmStd[static_cast<std::size_t>(i)] = std::sqrt(iRSS * gsl_matrix_get(covar, i, i));

    double mean = 0.0;
    for (double v : res) mean += v;
    mean /= d.nValid;
    double sd = 0.0;
    for (double v : res) sd += (v - mean) * (v - mean);
    sd = std::sqrt(sd / (d.nValid - 1));
    r.mean = mean;
    r.std = sd;
    pPost.stop();

    {
        prof::Scoped pAD(prof::ADTest);
        ADResult ad = adtest_mex(res, mean, sd, 0.05);
        r.hAD = ad.h;
        r.A2 = ad.A2;
    }

    if (prof::enabled) {
        const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - tFit0).count();
        prof::add_fit(prof::context, secs, iter, d.nF, d.nDF);
    }
    return r;
}

} // namespace cme
