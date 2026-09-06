// locmax2d (via ordfilt2 semantics) and bwconncomp/bwlabel (8-connectivity).
#pragma once
#include "cme/image.hpp"
#include <vector>
#include <cstddef>

namespace cme {

// locmax2d(img, maskSize) with keepFlat = 0:
//   fImg = ordfilt2(img, numEl, ones(m))          (window max)
//   fImg2 = ordfilt2(img, numEl-1, ones(m))       (second largest)
//   fImg(fImg2 == fImg) = 0;  fImg(fImg ~= img) = 0;  border strip (half window) = 0
// maskSize is made odd if even.  Returns an image equal to img at strict
// local maxima and 0 elsewhere (exactly like the MATLAB function).
ImageD locmax2d(const ImageD& img, int maskSize);

struct ConnComp {
    std::size_t numObjects = 0;
    // PixelIdxList in MATLAB order: components numbered as bwconncomp does
    // (see morphology.cpp); each list holds linear (column-major) indices in
    // ascending order.
    std::vector<std::vector<std::size_t>> pixelIdxList;
};

// bwconncomp(BW) with default 8-connectivity for a 2-D logical image.
ConnComp bwconncomp8(const ImageU8& bw);

// labelmatrix(CC): 0 for background, component number (1-based) otherwise.
ImageI32 labelmatrix(const ConnComp& cc, std::size_t ny, std::size_t nx);

// bwlabel(BW) (8-connectivity) == labelmatrix(bwconncomp(BW)) up to label
// numbering; only the partition is used by the port.
ImageI32 bwlabel8(const ImageU8& bw);

} // namespace cme
