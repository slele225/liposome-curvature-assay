// padarrayXT('symmetric') and conv2(hcol, hrow, A, 'valid') as used by
// pointSourceDetection.m.
#pragma once
#include "cme/image.hpp"
#include <vector>

namespace cme {

// padarrayXT(a, [p p], 'symmetric'): whole-sample symmetric extension
// (index sequence [1:M M-1:-1:2], period 2M-2), i.e. the border pixel is NOT
// repeated.  This differs from IPT padarray('symmetric'), which repeats it.
ImageD padarrayXT_symmetric(const ImageD& a, int p);

// conv2(hcol, hrow, A, 'valid'): convolve every column of A with hcol, then
// every row with hrow, return the fully-overlapping ('valid') part.
// hcol/hrow are applied as true convolutions (kernel flipped) - all kernels
// used on the path are symmetric so this only matters for consistency.
ImageD conv2_sep_valid(const std::vector<double>& hcol, const std::vector<double>& hrow, const ImageD& a);

} // namespace cme
