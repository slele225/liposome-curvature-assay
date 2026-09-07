#include "cme/stats.hpp"
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_errno.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace cme {

double norminv(double p) {
    if (!(p >= 0.0 && p <= 1.0)) return std::numeric_limits<double>::quiet_NaN();
    return gsl_cdf_ugaussian_Pinv(p);
}

double normcdf(double x, double mu, double sigma) {
    // MATLAB: 0.5 * erfc(-z / sqrt(2))
    const double z = (x - mu) / sigma;
    return 0.5 * gsl_sf_erfc(-z / std::sqrt(2.0));
}

double tcdf(double x, double nu) {
    if (std::isnan(x) || !(nu > 0.0)) return std::numeric_limits<double>::quiet_NaN();
    if (x == 0.0) return 0.5;
    if (nu == 1.0) {
        // p = xpos + acot(-x)/pi
        const double acot = std::atan(1.0 / (-x));
        return (x > 0 ? 1.0 : 0.0) + acot / M_PI;
    }
    if (nu > 1e7) return normcdf(x);
    // MATLAB: for x<0 (or via symmetry) p = betainc(nu/(nu+x^2), nu/2, 1/2)/2
    // when nu < x^2, else betainc(x^2/(nu+x^2), 1/2, nu/2, 'upper')/2; then
    // p = 1-p for x > 0.  gsl_cdf_tdist_P implements the same regularised
    // incomplete-beta evaluation.
    return gsl_cdf_tdist_P(x, nu);
}

static double ad_statistic(const std::vector<double>& xin, double mu, double sigma) {
    const std::size_t n = xin.size();
    if (n < 5) throw std::invalid_argument("adtest requires at least 5 samples");
    // per-thread scratch (called once per Gaussian fit)
    thread_local std::vector<double> x, z;
    x.assign(xin.begin(), xin.end());
    std::sort(x.begin(), x.end());
    z.resize(n);
    for (std::size_t i = 0; i < n; ++i) z[i] = normcdf(x[i], mu, sigma);
    // A2 = -n - 1/n * sum_{i=1}^n (2i-1) * (log(z_i) + log(1 - z_{n+1-i}))
    double s = 0.0;
    for (std::size_t i = 1; i <= n; ++i) {
        s += (2.0 * static_cast<double>(i) - 1.0) * (std::log(z[i - 1]) + std::log(1.0 - z[n - i]));
    }
    const double nn = static_cast<double>(n);
    return -nn - s / nn;
}

ADResult adtest1_case3(const std::vector<double>& xin, double mu, double sigma, double alpha) {
    ADResult r;
    const double nn = static_cast<double>(xin.size());
    double A2 = ad_statistic(xin, mu, sigma);
    A2 = A2 * (1.0 + 0.75 / nn + 2.25 / (nn * nn));
    double cval;
    if (alpha == 0.5) cval = 0.341;
    else if (alpha == 0.25) cval = 0.470;
    else if (alpha == 0.15) cval = 0.561;
    else if (alpha == 0.10) cval = 0.631;
    else if (alpha == 0.05) cval = 0.752;
    else if (alpha == 0.025) cval = 0.873;
    else if (alpha == 0.01) cval = 1.035;
    else if (alpha == 0.005) cval = 1.159;
    else throw std::invalid_argument("unsupported alpha for adtest");
    r.A2 = A2;
    r.h = A2 > cval;
    return r;
}

ADResult adtest_mex(const std::vector<double>& xin, double mu, double sigma, double alpha) {
    ADResult r;
    const double A2 = ad_statistic(xin, mu, sigma);
    // adtest1.m ctable(3,:): alpha 0.25 0.15 0.10 0.05 0.025 0.01 0.005 0.0025
    double cval;
    if (alpha == 0.25) cval = 1.072;
    else if (alpha == 0.15) cval = 1.430;
    else if (alpha == 0.10) cval = 1.743;
    else if (alpha == 0.05) cval = 2.308;
    else if (alpha == 0.025) cval = 2.898;
    else if (alpha == 0.01) cval = 3.702;
    else if (alpha == 0.005) cval = 4.324;
    else if (alpha == 0.0025) cval = 4.954;
    else throw std::invalid_argument("unsupported alpha for adtest");
    r.A2 = A2;
    r.h = A2 > cval;
    return r;
}

} // namespace cme
