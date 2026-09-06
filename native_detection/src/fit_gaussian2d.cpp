#include "cme/fit_gaussian2d.hpp"
#include "cme/stats.hpp"
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

// The MEX evaluates both the model and the Jacobian with |sigma| (verified on
// real candidate windows: with the signed derivative 10/356 free-sigma fits
// diverge from the MEX, with |sigma| every fit lands in the MEX's basin).

struct DataStruct {
    int nx = 0;                 // window side length
    int np = 0;                 // number of estimated parameters
    const double* pixels = nullptr;
    std::vector<double> gx, gy; // 1-D Gaussian factors
    std::vector<int> idx;       // linear indices of valid (non-NaN) pixels
    int nValid = 0;
    std::vector<int> estIdx;    // which of [xp yp A s c] are estimated
    double prmVect[5];
};

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
    update_factors(d, x);
    const int nx = d->nx;
    const double A = d->prmVect[2];
    const double c = d->prmVect[4];
    for (int i = 0; i < d->nValid; ++i) {
        const int idx = d->idx[i];
        // column-major: column = x index, row = y index
        const int col = idx / nx;
        const int row = idx % nx;
        gsl_vector_set(f, i, A * d->gx[col] * d->gy[row] + c - d->pixels[idx]);
    }
    return GSL_SUCCESS;
}

int gaussian_df(const gsl_vector* x, void* params, gsl_matrix* J) {
    DataStruct* d = static_cast<DataStruct*>(params);
    update_factors(d, x);
    const int nx = d->nx;
    const int b = nx / 2;
    const double xp = d->prmVect[0];
    const double yp = d->prmVect[1];
    const double A = d->prmVect[2];
    const double sigma = std::fabs(d->prmVect[3]);
    const double sigma2 = sigma * sigma;
    const double sigma3 = sigma2 * sigma;
    for (int i = 0; i < d->nValid; ++i) {
        const int idx = d->idx[i];
        const int col = idx / nx;
        const int row = idx % nx;
        const double xi = static_cast<double>(col - b) - xp;
        const double yi = static_cast<double>(row - b) - yp;
        const double g = d->gx[col] * d->gy[row];
        for (int k = 0; k < d->np; ++k) {
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

void silent_gsl_handler(const char*, const char*, int, int) {}

} // namespace

FitResult fitGaussian2D(const ImageD& window, const std::array<double, 5>& prm0, const std::string& modeIn,
                        const FitOptions& opt) {
    if (window.nx() != window.ny()) throw std::invalid_argument("fitGaussian2D: window must be square");
    const int nx = static_cast<int>(window.nx());
    const int N = nx * nx;

    DataStruct d;
    d.nx = nx;
    d.pixels = window.data();
    d.gx.assign(nx, 0.0);
    d.gy.assign(nx, 0.0);
    for (int i = 0; i < N; ++i) if (!std::isnan(d.pixels[i])) d.idx.push_back(i);
    d.nValid = static_cast<int>(d.idx.size());
    for (int i = 0; i < 5; ++i) d.prmVect[i] = prm0[i];

    // mode -> estIdx, canonical order x y a s c
    std::string mode = modeIn;
    for (char& ch : mode) ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
    static const char refMode[] = "xyasc";
    for (int i = 0; i < 5; ++i) {
        if (mode.find(refMode[i]) != std::string::npos) d.estIdx.push_back(i);
    }
    d.np = static_cast<int>(d.estIdx.size());
    if (d.np == 0) throw std::invalid_argument("fitGaussian2D: no parameters to estimate");
    if (d.nValid <= d.np) throw std::invalid_argument("fitGaussian2D: not enough valid pixels");

    std::vector<double> x_init(d.np);
    for (int i = 0; i < d.np; ++i) x_init[i] = d.prmVect[d.estIdx[i]];

    gsl_error_handler_t* old = gsl_set_error_handler(&silent_gsl_handler);

    g116_solver* s = g116_solver_alloc(d.nValid, d.np);
    if (!s) { gsl_set_error_handler(old); throw std::runtime_error("fitGaussian2D: solver allocation failed"); }
    g116_function_fdf f;
    f.f = &gaussian_f;
    f.df = &gaussian_df;
    f.fdf = &gaussian_fdf;
    f.n = static_cast<size_t>(d.nValid);
    f.p = static_cast<size_t>(d.np);
    f.params = &d;
    gsl_vector_view xv = gsl_vector_view_array(x_init.data(), d.np);
    g116_solver_set(s, &f, &xv.vector);

    int status = GSL_CONTINUE;
    int iter = 0;
    do {
        ++iter;
        status = g116_solver_iterate(s, &f);
        if (status) break;
        status = g116_test_delta(g116_solver_dx(s), g116_solver_x(s), opt.eAbs, opt.eRel);
    } while (status == GSL_CONTINUE && iter < opt.maxIter);

    FitResult r;
    r.iterations = iter;
    r.status = status;
    r.nValid = d.nValid;
    r.estIdx = d.estIdx;
    const gsl_vector* xs = g116_solver_x(s);
    for (int i = 0; i < d.np; ++i) d.prmVect[d.estIdx[i]] = gsl_vector_get(xs, i);
    d.prmVect[3] = std::fabs(d.prmVect[3]);
    for (int i = 0; i < 5; ++i) r.prm[i] = d.prmVect[i];

    // residuals (as stored in the solver at the final iterate)
    const gsl_vector* fs = g116_solver_f(s);
    std::vector<double> res(d.nValid);
    double RSS = 0.0;
    for (int i = 0; i < d.nValid; ++i) {
        res[i] = gsl_vector_get(fs, i);
        RSS += res[i] * res[i];
    }
    r.RSS = RSS;

    // covariance from the Jacobian at the final iterate
    gsl_matrix* covar = gsl_matrix_alloc(d.np, d.np);
    g116_covar(g116_solver_J(s), 0.0, covar);
    const double iRSS = RSS / static_cast<double>(d.nValid - d.np - 1);
    r.prmStd.resize(d.np);
    for (int i = 0; i < d.np; ++i) r.prmStd[i] = std::sqrt(iRSS * gsl_matrix_get(covar, i, i));
    gsl_matrix_free(covar);

    double mean = 0.0;
    for (double v : res) mean += v;
    mean /= d.nValid;
    double sd = 0.0;
    for (double v : res) sd += (v - mean) * (v - mean);
    sd = std::sqrt(sd / (d.nValid - 1));
    r.mean = mean;
    r.std = sd;

    ADResult ad = adtest_mex(res, mean, sd, 0.05);
    r.hAD = ad.h;
    r.A2 = ad.A2;

    g116_solver_free(s);
    gsl_set_error_handler(old);
    return r;
}

} // namespace cme
