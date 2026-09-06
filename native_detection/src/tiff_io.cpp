#include "cme/tiff_io.hpp"
#include <tiffio.h>
#ifdef _WIN32
#include <windows.h>
#endif
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <cstring>

namespace cme {

namespace {

struct TiffCloser {
    void operator()(TIFF* t) const { if (t) TIFFClose(t); }
};
using TiffPtr = std::unique_ptr<TIFF, TiffCloser>;

void quiet_handler(const char*, const char*, va_list) {}

TiffPtr open_tiff(const std::string& path, const char* mode) {
    TIFFSetWarningHandler(quiet_handler);   // MATLAB readtiff silences unknown-tag warnings too
#ifdef _WIN32
    // libtiff on Windows accepts UTF-8 paths with TIFFOpenW only; use the
    // wide variant when the path is not plain ASCII.
    bool ascii = true;
    for (unsigned char c : path) if (c >= 0x80) { ascii = false; break; }
    if (!ascii) {
        int n = MultiByteToWideChar(CP_UTF8, 0, path.c_str(), -1, nullptr, 0);
        std::wstring w(static_cast<std::size_t>(n), L'\0');
        MultiByteToWideChar(CP_UTF8, 0, path.c_str(), -1, &w[0], n);
        std::wstring wmode(mode, mode + std::strlen(mode));
        return TiffPtr(TIFFOpenW(w.c_str(), mode));
    }
#endif
    return TiffPtr(TIFFOpen(path.c_str(), mode));
}

} // namespace

TiffInfo tiff_info(const std::string& path) {
    TiffPtr t = open_tiff(path, "r");
    if (!t) throw std::runtime_error("Cannot open TIFF: " + path);
    TiffInfo info;
    uint32_t w = 0, h = 0;
    uint16_t bps = 0, spp = 1, fmt = SAMPLEFORMAT_UINT;
    TIFFGetField(t.get(), TIFFTAG_IMAGEWIDTH, &w);
    TIFFGetField(t.get(), TIFFTAG_IMAGELENGTH, &h);
    TIFFGetField(t.get(), TIFFTAG_BITSPERSAMPLE, &bps);
    TIFFGetFieldDefaulted(t.get(), TIFFTAG_SAMPLESPERPIXEL, &spp);
    TIFFGetFieldDefaulted(t.get(), TIFFTAG_SAMPLEFORMAT, &fmt);
    info.width = w;
    info.height = h;
    info.bitsPerSample = bps;
    info.samplesPerPixel = spp;
    info.sampleFormat = fmt;
    info.numDirectories = 0;
    do { ++info.numDirectories; } while (TIFFReadDirectory(t.get()));
    return info;
}

template <typename T>
static void unpack_rows(TIFF* t, ImageD& out, std::size_t ny, std::size_t nx, std::size_t spp) {
    std::vector<T> buf(static_cast<std::size_t>(TIFFScanlineSize(t)) / sizeof(T) + 1);
    for (std::size_t y = 0; y < ny; ++y) {
        if (TIFFReadScanline(t, buf.data(), static_cast<uint32_t>(y)) < 0)
            throw std::runtime_error("TIFF scanline read failed");
        for (std::size_t x = 0; x < nx; ++x) out(y, x) = static_cast<double>(buf[x * spp]);
    }
}

ImageD read_tiff_frame_double(const std::string& path, std::size_t frame) {
    TiffPtr t = open_tiff(path, "r");
    if (!t) throw std::runtime_error("Cannot open TIFF: " + path);
    if (frame < 1) throw std::invalid_argument("frame must be >= 1");
    if (!TIFFSetDirectory(t.get(), static_cast<uint32_t>(frame - 1)))
        throw std::runtime_error("TIFF directory not found: " + path);
    uint32_t w = 0, h = 0;
    uint16_t bps = 0, spp = 1, fmt = SAMPLEFORMAT_UINT, planar = PLANARCONFIG_CONTIG;
    TIFFGetField(t.get(), TIFFTAG_IMAGEWIDTH, &w);
    TIFFGetField(t.get(), TIFFTAG_IMAGELENGTH, &h);
    TIFFGetField(t.get(), TIFFTAG_BITSPERSAMPLE, &bps);
    TIFFGetFieldDefaulted(t.get(), TIFFTAG_SAMPLESPERPIXEL, &spp);
    TIFFGetFieldDefaulted(t.get(), TIFFTAG_SAMPLEFORMAT, &fmt);
    TIFFGetFieldDefaulted(t.get(), TIFFTAG_PLANARCONFIG, &planar);
    if (planar != PLANARCONFIG_CONTIG) throw std::runtime_error("Unsupported planar TIFF: " + path);
    if (spp != 1) throw std::runtime_error("Only grayscale TIFFs are supported: " + path);
    ImageD out(h, w);
    if (fmt == SAMPLEFORMAT_UINT || fmt == SAMPLEFORMAT_VOID) {
        if (bps == 8) unpack_rows<uint8_t>(t.get(), out, h, w, spp);
        else if (bps == 16) unpack_rows<uint16_t>(t.get(), out, h, w, spp);
        else if (bps == 32) unpack_rows<uint32_t>(t.get(), out, h, w, spp);
        else throw std::runtime_error("Unsupported bit depth");
    } else if (fmt == SAMPLEFORMAT_INT) {
        if (bps == 8) unpack_rows<int8_t>(t.get(), out, h, w, spp);
        else if (bps == 16) unpack_rows<int16_t>(t.get(), out, h, w, spp);
        else if (bps == 32) unpack_rows<int32_t>(t.get(), out, h, w, spp);
        else throw std::runtime_error("Unsupported bit depth");
    } else if (fmt == SAMPLEFORMAT_IEEEFP) {
        if (bps == 32) unpack_rows<float>(t.get(), out, h, w, spp);
        else if (bps == 64) unpack_rows<double>(t.get(), out, h, w, spp);
        else throw std::runtime_error("Unsupported float depth");
    } else {
        throw std::runtime_error("Unsupported sample format");
    }
    return out;
}

void write_tiff_uint8(const std::string& path, const ImageU8& img, bool append) {
    TiffPtr t = open_tiff(path, append ? "a" : "w");
    if (!t) throw std::runtime_error("Cannot open TIFF for writing: " + path);
    const uint32_t w = static_cast<uint32_t>(img.nx()), h = static_cast<uint32_t>(img.ny());
    TIFFSetField(t.get(), TIFFTAG_IMAGEWIDTH, w);
    TIFFSetField(t.get(), TIFFTAG_IMAGELENGTH, h);
    TIFFSetField(t.get(), TIFFTAG_BITSPERSAMPLE, 8);
    TIFFSetField(t.get(), TIFFTAG_SAMPLESPERPIXEL, 1);
    TIFFSetField(t.get(), TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField(t.get(), TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
    TIFFSetField(t.get(), TIFFTAG_COMPRESSION, COMPRESSION_LZW);
    TIFFSetField(t.get(), TIFFTAG_ROWSPERSTRIP, TIFFDefaultStripSize(t.get(), h));
    std::vector<uint8_t> row(w);
    for (uint32_t y = 0; y < h; ++y) {
        for (uint32_t x = 0; x < w; ++x) row[x] = img(y, x);
        if (TIFFWriteScanline(t.get(), row.data(), y) < 0) throw std::runtime_error("TIFF write failed");
    }
    TIFFWriteDirectory(t.get());
}

} // namespace cme
