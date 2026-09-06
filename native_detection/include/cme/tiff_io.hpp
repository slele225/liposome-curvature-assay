// TIFF access through libtiff: replaces imfinfo/imread/readtiff/imwrite.
#pragma once
#include "cme/image.hpp"
#include <string>
#include <vector>

namespace cme {

struct TiffInfo {
    std::size_t width = 0;
    std::size_t height = 0;
    std::size_t numDirectories = 0;   // numel(imfinfo(file))
    int bitsPerSample = 0;
    int sampleFormat = 1;             // 1 uint, 2 int, 3 float
    int samplesPerPixel = 1;
};

TiffInfo tiff_info(const std::string& path);

// double(imread(path)) / double(readtiff(path, frame)); frame is 1-based.
ImageD read_tiff_frame_double(const std::string& path, std::size_t frame = 1);

// imwrite(uint8(255*mask), path, 'tif', 'compression', 'lzw' [, 'writemode', 'append'])
void write_tiff_uint8(const std::string& path, const ImageU8& img, bool append);

} // namespace cme
