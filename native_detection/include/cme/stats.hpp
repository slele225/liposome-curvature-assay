// Statistical primitives with MATLAB parameterisation.
#pragma once
#include <vector>
#include <cstddef>

namespace cme {

// norminv(p, 0, 1)
double norminv(double p);
// normcdf(x, mu, sigma)
double normcdf(double x, double mu = 0.0, double sigma = 1.0);
// tcdf(x, nu) for real nu; nu <= 0 or NaN -> NaN (MATLAB semantics)
double tcdf(double x, double nu);

struct ADResult {
    bool h = false;     // true = normality rejected
    double A2 = 0.0;    // statistic that was compared with the critical value
};

// Anderson-Darling normality test as implemented in cmeAnalysis adtest1.m for
// the case "mu and sigma estimated" (case 3): modified statistic
// A2*(1 + 0.75/n + 2.25/n^2) against ctable(4,:).  Used only by unit tests.
ADResult adtest1_case3(const std::vector<double>& x, double mu, double sigma, double alpha = 0.05);

// Anderson-Darling decision exactly as produced by the original fitGaussian2D
// MEX (res.hAD): raw A2 (no small-sample correction) compared with the
// critical values of the "mu known, sigma estimated" case (adtest1.m
// ctable(3,:); 2.308 at alpha = 0.05).  Established empirically from ~940
// MEX decisions (tests/ad_reference*.txt), see PORTING_NOTES.md section 2.3.
ADResult adtest_mex(const std::vector<double>& x, double mu, double sigma, double alpha = 0.05);

} // namespace cme
