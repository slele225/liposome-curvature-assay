// Port of the non-interactive parts of loadConditionData.m
#pragma once
#include <string>
#include <vector>
#include <cstddef>

namespace cme {

struct MovieData {
    std::string cellPath;                 // with trailing separator
    std::vector<std::string> channels;    // channel directories (trailing separator), given order
    std::string source;                   // channels[0]
    std::string date;
    double framerate = 2.0;
    std::size_t imageHeight = 0, imageWidth = 0;
    std::size_t movieLength = 0;
    // framePaths[c]: either one entry (single multi-page file, singleFile=true)
    // or one entry per frame
    std::vector<std::vector<std::string>> framePaths;
    bool singleFile = false;
    std::vector<std::string> markers;
    double NA = 1.49, M = 108, pixelSize = 6.5e-6;
    bool hasFrames = false;
};

struct ConditionOptions {
    std::string movieSelector = "cell";
    bool strictSelector = false;
    bool ignoreEmptyFolders = false;
};

// condDir: condition directory; chNames: channel folder names in the exact
// order the user would have selected them (first = master).
std::vector<MovieData> loadConditionData(const std::string& condDir, const std::vector<std::string>& chNames,
                                         const std::vector<std::string>& markers, const ConditionOptions& opt = ConditionOptions());

// helper: path of the frame f (1-based) of channel c and the TIFF directory to read
struct FrameRef { std::string path; std::size_t directory; };
FrameRef frame_ref(const MovieData& d, std::size_t c, std::size_t f);

} // namespace cme
