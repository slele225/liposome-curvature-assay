#include "cme/conv.hpp"
#include <stdexcept>

// Parallelisation note: the loops below are split over *output columns*
// only.  Each output element is produced by the same sequence of
// multiply-adds in the same order as in the single-threaded version, so the
// results are bit-identical for any thread count.

namespace cme {

static std::vector<std::size_t> symmetric_xt_indices(std::size_t M, int p) {
    // dimNums = [1:M M-1:-1:2]; div = 2M-2; idx = dimNums(mod(-p:M+p-1, div)+1)
    std::vector<std::size_t> dimNums;
    std::size_t div;
    if (M > 1) {
        dimNums.reserve(2 * M - 2);
        for (std::size_t i = 0; i < M; ++i) dimNums.push_back(i);
        for (std::size_t i = M - 1; i >= 2; --i) dimNums.push_back(i - 1);
        div = 2 * M - 2;
    } else {
        dimNums = {0, 0};
        div = 2;
    }
    std::vector<std::size_t> idx;
    idx.reserve(M + 2 * static_cast<std::size_t>(p));
    for (long k = -p; k < static_cast<long>(M) + p; ++k) {
        long m = k % static_cast<long>(div);
        if (m < 0) m += static_cast<long>(div);
        idx.push_back(dimNums[static_cast<std::size_t>(m)]);
    }
    return idx;
}

ImageD padarrayXT_symmetric(const ImageD& a, int p, int threads) {
    if (p < 0) throw std::invalid_argument("negative pad");
    auto ry = symmetric_xt_indices(a.ny(), p);
    auto rx = symmetric_xt_indices(a.nx(), p);
    ImageD b(ry.size(), rx.size());
    const long NX = static_cast<long>(rx.size());
    const std::size_t NY = ry.size();
    #pragma omp parallel for schedule(static) num_threads(threads > 1 ? threads : 1) if(threads > 1)
    for (long xx = 0; xx < NX; ++xx) {
        const std::size_t x = static_cast<std::size_t>(xx);
        const std::size_t sx = rx[x];
        double* out = b.data() + x * NY;
        const double* src = a.data() + sx * a.ny();
        for (std::size_t y = 0; y < NY; ++y) out[y] = src[ry[y]];
    }
    return b;
}

ImageD conv2_sep_valid(const std::vector<double>& hcol, const std::vector<double>& hrow, const ImageD& a, int threads) {
    const std::size_t ny = a.ny(), nx = a.nx();
    const std::size_t kc = hcol.size(), kr = hrow.size();
    if (ny < kc || nx < kr) throw std::invalid_argument("conv2 valid: kernel larger than image");
    const std::size_t oy = ny - kc + 1;
    const std::size_t ox = nx - kr + 1;
    const int nth = threads > 1 ? threads : 1;

    // Column pass: tmp(y', x) = sum_j hcol(j) * a(y' + kc-1 - j, x)
    ImageD tmp(oy, nx);
    #pragma omp parallel for schedule(static) num_threads(nth) if(nth > 1)
    for (long xx = 0; xx < static_cast<long>(nx); ++xx) {
        const std::size_t x = static_cast<std::size_t>(xx);
        const double* col = a.data() + x * ny;
        double* out = tmp.data() + x * oy;
        for (std::size_t y = 0; y < oy; ++y) {
            double s = 0.0;
            for (std::size_t j = 0; j < kc; ++j) s += hcol[j] * col[y + kc - 1 - j];
            out[y] = s;
        }
    }
    // Row pass: r(y, x') = sum_j hrow(j) * tmp(y, x' + kr-1 - j)
    ImageD r(oy, ox);
    #pragma omp parallel for schedule(static) num_threads(nth) if(nth > 1)
    for (long xx = 0; xx < static_cast<long>(ox); ++xx) {
        const std::size_t x = static_cast<std::size_t>(xx);
        double* out = r.data() + x * oy;
        for (std::size_t y = 0; y < oy; ++y) out[y] = 0.0;
        for (std::size_t j = 0; j < kr; ++j) {
            const double h = hrow[j];
            const double* src = tmp.data() + (x + kr - 1 - j) * oy;
            for (std::size_t y = 0; y < oy; ++y) out[y] += h * src[y];
        }
    }
    return r;
}

} // namespace cme
