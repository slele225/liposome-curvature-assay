// Column-major image container mirroring MATLAB's memory layout.
//
// img(y, x) with 0-based y (row) and x (column); linear index = y + x*ny,
// identical to MATLAB's sub2ind(size(img), y+1, x+1) - 1.  Keeping the
// column-major layout means that MATLAB's find(), bwconncomp() pixel lists and
// candidate ordering are reproduced without any index translation.
#pragma once
#include <cstddef>
#include <vector>
#include <algorithm>
#include <stdexcept>

namespace cme {

template <typename T>
class Image {
public:
    Image() = default;
    Image(std::size_t ny, std::size_t nx, T fill = T()) : ny_(ny), nx_(nx), data_(ny * nx, fill) {}

    std::size_t ny() const { return ny_; }
    std::size_t nx() const { return nx_; }
    std::size_t size() const { return data_.size(); }
    bool empty() const { return data_.empty(); }

    T& operator()(std::size_t y, std::size_t x) { return data_[y + x * ny_]; }
    const T& operator()(std::size_t y, std::size_t x) const { return data_[y + x * ny_]; }
    T& operator[](std::size_t i) { return data_[i]; }
    const T& operator[](std::size_t i) const { return data_[i]; }

    T* data() { return data_.data(); }
    const T* data() const { return data_.data(); }
    std::vector<T>& vec() { return data_; }
    const std::vector<T>& vec() const { return data_; }

    // MATLAB sub2ind (0-based in, 0-based out)
    std::size_t idx(std::size_t y, std::size_t x) const { return y + x * ny_; }

    template <typename U>
    Image<U> cast() const {
        Image<U> out(ny_, nx_);
        for (std::size_t i = 0; i < data_.size(); ++i) out[i] = static_cast<U>(data_[i]);
        return out;
    }

    T minval() const { return *std::min_element(data_.begin(), data_.end()); }
    T maxval() const { return *std::max_element(data_.begin(), data_.end()); }

private:
    std::size_t ny_ = 0, nx_ = 0;
    std::vector<T> data_;
};

using ImageD = Image<double>;
using ImageU8 = Image<unsigned char>;
using ImageI32 = Image<int>;

} // namespace cme
