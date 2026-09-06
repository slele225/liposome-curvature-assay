// TSV outputs replacing detection_v2.mat
#pragma once
#include "cme/run_detection.hpp"
#include <string>
#include <vector>

namespace cme {

// Detailed per-movie table: one row per detection with every frameInfo field.
void write_movie_tsv(const std::string& path, const DetectionOutput& out, const std::vector<std::string>& chNames);

// Condition-level tables.
void write_condition_tables(const std::string& dir, const std::vector<DetectionOutput>& outs,
                            const std::vector<std::string>& chNames, const std::vector<double>& sigma);

// Masks (dmasks.tif or Masks/dmask_NNN.tif) in the MATLAB layout.
void write_masks(const DetectionOutput& out);

} // namespace cme
