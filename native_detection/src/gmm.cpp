#include "cme/gmm.hpp"
#include <algorithm>
#include <cmath>
#include <limits>

namespace cme {

namespace {
constexpr double PI = 3.14159265358979323846;

double eps_of(double x) {
    // MATLAB eps(x): spacing of doubles at |x|
    x = std::fabs(x);
    if (x == 0.0) return std::numeric_limits<double>::denorm_min();
    int e;
    std::frexp(x, &e);
    return std::ldexp(1.0, e - 53);
}

double var_sample(const std::vector<double>& X) {
    const std::size_t n = X.size();
    double m = 0.0;
    for (double v : X) m += v;
    m /= static_cast<double>(n);
    double s = 0.0;
    for (double v : X) s += (v - m) * (v - m);
    return s / static_cast<double>(n - 1);
}
} // namespace

GmmResult gmdistribution_fit_1d(const std::vector<double>& X, int k, MatlabTwister& rng,
                                int maxIter, double tolFun, double probTol) {
    const std::size_t n = X.size();
    const int d = 1;
    if (static_cast<long>(n) <= d) throw GmmFitError("stats:gmdistribution:TooFewN");
    if (static_cast<long>(n) <= k) throw GmmFitError("stats:gmdistribution:TooManyClusters");
    for (double v : X) if (std::isnan(v)) throw GmmFitError("NaN in data (MATLAB would drop rows; not expected on this path)");

    const double varX = var_sample(X);
    if (varX < eps_of(varX)) throw GmmFitError("stats:gmdistribution:ZeroVariance");

    GmmResult R;
    R.k = k;

    // ---- plusInitParam (k-means++) ----
    std::vector<double> mu(k), Sigma(k, varX), p(k, 1.0 / k);
    std::vector<long> index(k, 0);
    {
        // [C(1,:), index(1)] = datasample(X,1)  -> randi(n)
        index[0] = rng.randi(static_cast<long>(n));
        mu[0] = X[static_cast<std::size_t>(index[0] - 1)];
        std::vector<double> minDist(n, std::numeric_limits<double>::infinity());
        for (int ii = 1; ii < k; ++ii) {
            double denominator = 0.0;
            for (std::size_t i = 0; i < n; ++i) {
                const double dd = (X[i] - mu[ii - 1]) / std::sqrt(varX);
                minDist[i] = std::min(minDist[i], dd * dd);
                denominator += minDist[i];
            }
            if (denominator == 0.0 || std::isinf(denominator)) {
                // C(ii:k,:) = datasample(X, k-ii+1, 1, 'Replace', false)  -> randperm(n, k-ii+1)
                std::vector<long> rp = rng.randperm_k(static_cast<long>(n), k - ii);
                for (int jj = ii; jj < k; ++jj) {
                    index[jj] = rp[static_cast<std::size_t>(jj - ii)];
                    mu[jj] = X[static_cast<std::size_t>(index[jj] - 1)];
                }
                break;
            }
            // sampleProbability = minDist/denominator; datasample(...,'Weights',sampleProbability)
            // -> internal.stats.wswor(w,1): inverse CDF with one rand draw
            // p = w/sum(w); edges = min([0 cumsum(p)],1); edges(end)=1; histcounts(u, edges)
            const double u = rng.rand();
            double sumw = 0.0;
            for (std::size_t i = 0; i < n; ++i) sumw += minDist[i] / denominator;
            std::vector<double> edges(n + 1);
            edges[0] = 0.0;
            double cs = 0.0;
            for (std::size_t i = 0; i < n; ++i) {
                cs += (minDist[i] / denominator) / sumw;
                edges[i + 1] = std::min(cs, 1.0);
            }
            edges[n] = 1.0;
            // histcounts: bin i (1-based) if edges(i) <= u < edges(i+1); last bin includes right edge
            long bin = 0;
            for (std::size_t i = 0; i < n; ++i) {
                const bool last = (i + 1 == n);
                if (u >= edges[i] && (u < edges[i + 1] || (last && u <= edges[i + 1]))) { bin = static_cast<long>(i) + 1; break; }
            }
            if (bin == 0) throw GmmFitError("wswor: u outside edges");
            index[ii] = bin;
            mu[ii] = X[static_cast<std::size_t>(bin - 1)];
        }
    }
    R.initIdx = index;
    R.initMu = mu;

    // ---- gmcluster_learn ----
    const std::size_t th = static_cast<std::size_t>(std::floor(static_cast<double>(n) * 0.4));
    (void)th;  // the subset optimisation is numerically equivalent (zeros add nothing)
    double ll_old = -std::numeric_limits<double>::infinity();
    double ll = 0.0;
    std::vector<double> post(n * static_cast<std::size_t>(k));
    std::vector<double> log_lh(n * static_cast<std::size_t>(k));
    int iter = 0;
    bool converged = false;
    for (iter = 1; iter <= maxIter; ++iter) {
        // wdensity
        for (int j = 0; j < k; ++j) {
            const double S = Sigma[static_cast<std::size_t>(j)];
            // [L,f] = chol(S): fails if S <= 0 (or NaN)
            if (!(S > 0.0)) throw GmmFitError("stats:gmdistribution:IllCondCovIter");
            const double Lc = std::sqrt(S);
            if (std::fabs(Lc) < eps_of(std::fabs(Lc)) * 1.0) throw GmmFitError("stats:gmdistribution:IllCondCovIter");
            const double logDetSigma = 2.0 * std::log(Lc);
            const double lp = std::log(p[static_cast<std::size_t>(j)]);
            for (std::size_t i = 0; i < n; ++i) {
                const double z = (X[i] - mu[static_cast<std::size_t>(j)]) / Lc;
                double v = z * z;
                v = -0.5 * (v + logDetSigma);
                log_lh[i + static_cast<std::size_t>(j) * n] = v + lp - d * std::log(2.0 * PI) / 2.0;
            }
        }
        // estep with probability tolerance
        ll = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            double maxll = -std::numeric_limits<double>::infinity();
            for (int j = 0; j < k; ++j) maxll = std::max(maxll, log_lh[i + static_cast<std::size_t>(j) * n]);
            double density = 0.0;
            for (int j = 0; j < k; ++j) {
                const double e = std::exp(log_lh[i + static_cast<std::size_t>(j) * n] - maxll);
                post[i + static_cast<std::size_t>(j) * n] = e;
                density += e;
            }
            ll += std::log(density) + maxll;
            for (int j = 0; j < k; ++j) post[i + static_cast<std::size_t>(j) * n] /= density;
            // post(post < probtol) = 0; renormalise
            double dens2 = 0.0;
            for (int j = 0; j < k; ++j) {
                double& pv = post[i + static_cast<std::size_t>(j) * n];
                if (pv < probTol) pv = 0.0;
                dens2 += pv;
            }
            for (int j = 0; j < k; ++j) post[i + static_cast<std::size_t>(j) * n] /= dens2;
        }
        const double llDiff = ll - ll_old;
        if (llDiff >= 0.0 && llDiff < tolFun * std::fabs(ll)) { converged = true; break; }
        ll_old = ll;

        // M-step
        std::vector<double> Nj(k, 0.0);
        for (int j = 0; j < k; ++j) {
            double s = 0.0;
            for (std::size_t i = 0; i < n; ++i) s += post[i + static_cast<std::size_t>(j) * n];
            Nj[static_cast<std::size_t>(j)] = s;
        }
        for (int j = 0; j < k; ++j) {
            const double N = Nj[static_cast<std::size_t>(j)];
            if (N == 0.0) continue;
            const double* pj = &post[static_cast<std::size_t>(j) * n];
            double m = 0.0;
            for (std::size_t i = 0; i < n; ++i) m += pj[i] * X[i];
            m /= N;
            mu[static_cast<std::size_t>(j)] = m;
            // Xcentered = sqrt(post_j) .* (X - mu); Sigma = Xcentered'*Xcentered / N
            double s = 0.0;
            for (std::size_t i = 0; i < n; ++i) {
                if (pj[i] > 0.0) {
                    const double xc = std::sqrt(pj[i]) * (X[i] - m);
                    s += xc * xc;
                }
            }
            Sigma[static_cast<std::size_t>(j)] = s / N;
        }
        double sumN = 0.0;
        for (double v : Nj) sumN += v;
        for (int j = 0; j < k; ++j) p[static_cast<std::size_t>(j)] = Nj[static_cast<std::size_t>(j)] / sumN;
    }
    if (iter > maxIter) iter = maxIter;   // MATLAB: optimInfo.Iters = iter (loop variable after completion)
    R.iters = iter;
    R.converged = converged;
    R.mu = mu;
    R.Sigma = Sigma;
    R.PComponents = p;
    R.NlogL = -ll;
    const double nParam = static_cast<double>(k) * d * (d + 1) / 2.0 + (k - 1) + k * d;
    R.BIC = 2.0 * R.NlogL + nParam * std::log(static_cast<double>(n));
    return R;
}

} // namespace cme
