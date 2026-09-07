#include "cme/morphology.hpp"
#include <algorithm>
#include <numeric>
#include <stdexcept>

namespace cme {

ImageD locmax2d(const ImageD& img, int maskSize, int threads) {
    if (maskSize % 2 == 0) maskSize += 1;
    const int b = (maskSize - 1) / 2;
    const std::size_t ny = img.ny(), nx = img.nx();
    ImageD out(ny, nx, 0.0);
    if (static_cast<int>(ny) <= 2 * b || static_cast<int>(nx) <= 2 * b) return out;

    // Only pixels outside the border strip can be non-zero; for those the
    // window lies fully inside the image, so ordfilt2's zero padding never
    // enters the computation.
    const long xEnd = static_cast<long>(nx) - b;
    #pragma omp parallel for schedule(static) num_threads(threads > 1 ? threads : 1) if(threads > 1)
    for (long xx = b; xx < xEnd; ++xx) {
        const std::size_t x = static_cast<std::size_t>(xx);
        for (std::size_t y = static_cast<std::size_t>(b); y + static_cast<std::size_t>(b) < ny; ++y) {
            const double v = img(y, x);
            // max and second max over the window
            double m1 = -HUGE_VAL, m2 = -HUGE_VAL;
            for (std::size_t xx = x - b; xx <= x + b; ++xx) {
                const double* col = img.data() + xx * ny;
                for (std::size_t yy = y - b; yy <= y + b; ++yy) {
                    const double w = col[yy];
                    if (w > m1) { m2 = m1; m1 = w; }
                    else if (w > m2) { m2 = w; }
                }
            }
            // fImg = m1; fImg(m2 == m1) = 0; fImg(fImg ~= img) = 0
            if (m1 == m2) continue;
            if (m1 != v) continue;
            out(y, x) = m1;
        }
    }
    return out;
}

// bwconncomp for 2-D uses run-length based union-find in column-major scan
// order; components are numbered in the order their first pixel (lowest
// linear index) is encountered.  We reproduce that numbering: a simple
// column-major scan with a flood fill gives components ordered by their
// minimum linear index, which is exactly the order of first encounter.
ConnComp bwconncomp8(const ImageU8& bw) {
    const std::size_t ny = bw.ny(), nx = bw.nx();
    ConnComp cc;
    std::vector<int> label(ny * nx, 0);
    std::vector<std::size_t> stack;
    const long NY = static_cast<long>(ny), NX = static_cast<long>(nx);
    for (std::size_t i = 0; i < ny * nx; ++i) {
        if (!bw[i] || label[i]) continue;
        const int id = static_cast<int>(cc.pixelIdxList.size()) + 1;
        std::vector<std::size_t> pixels;
        stack.clear();
        stack.push_back(i);
        label[i] = id;
        while (!stack.empty()) {
            const std::size_t p = stack.back();
            stack.pop_back();
            pixels.push_back(p);
            const long x = static_cast<long>(p / ny), y = static_cast<long>(p % ny);
            for (long dx = -1; dx <= 1; ++dx) {
                const long xx = x + dx;
                if (xx < 0 || xx >= NX) continue;
                for (long dy = -1; dy <= 1; ++dy) {
                    const long yy = y + dy;
                    if (yy < 0 || yy >= NY) continue;
                    const std::size_t q = static_cast<std::size_t>(yy) + static_cast<std::size_t>(xx) * ny;
                    if (bw[q] && !label[q]) {
                        label[q] = id;
                        stack.push_back(q);
                    }
                }
            }
        }
        std::sort(pixels.begin(), pixels.end());
        cc.pixelIdxList.push_back(std::move(pixels));
    }
    cc.numObjects = cc.pixelIdxList.size();
    return cc;
}

ImageI32 labelmatrix(const ConnComp& cc, std::size_t ny, std::size_t nx) {
    ImageI32 L(ny, nx, 0);
    for (std::size_t k = 0; k < cc.pixelIdxList.size(); ++k) {
        for (std::size_t p : cc.pixelIdxList[k]) L[p] = static_cast<int>(k + 1);
    }
    return L;
}

ImageI32 bwlabel8(const ImageU8& bw) {
    return labelmatrix(bwconncomp8(bw), bw.ny(), bw.nx());
}

} // namespace cme
