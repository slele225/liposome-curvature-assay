// Text/binary dumps of intermediates for regression against MATLAB.
#pragma once
#include "cme/image.hpp"
#include "cme/point_source_detection.hpp"
#include "cme/condition_data.hpp"
#include "cme/psf_sigma.hpp"
#include <string>
#include <vector>

namespace cme {

struct FrameInfo;
struct SigmaEstimate;

std::string dump_movie_name(const MovieData& d, std::size_t movieIndex);

// raw little-endian float64 column-major with a small text header file
void dump_image_bin(const std::string& path, const ImageD& img);
void dump_image_bin(const std::string& path, const ImageU8& img);

void dump_psd_debug(const std::string& dir, std::size_t frame, const PSDDebug& dbg, const PSDResult& res);
void dump_pstruct_tsv(const std::string& path, const PStruct& P);
void dump_frame_info(const std::string& dir, const FrameInfo& F, const MovieData& d);
void dump_sigma_estimate(const std::string& dir, const SigmaEstimate& S);

} // namespace cme
