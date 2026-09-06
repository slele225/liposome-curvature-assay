#include "cme/condition_data.hpp"
#include "cme/tiff_io.hpp"
#include <algorithm>
#include <cctype>
#include <filesystem>
#include <regex>
#include <set>
#include <stdexcept>
#include <iostream>

namespace fs = std::filesystem;

namespace cme {

namespace {

const char SEP = static_cast<char>(fs::path::preferred_separator);

std::string with_sep(std::string s) {
    if (s.empty() || s.back() != SEP) s.push_back(SEP);
    return s;
}

std::string lower(std::string s) {
    for (char& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return s;
}

// MATLAB dir() order on Windows: case-insensitive lexicographic by name.
bool name_less(const std::string& a, const std::string& b) {
    const std::string la = lower(a), lb = lower(b);
    if (la != lb) return la < lb;
    return a < b;
}

std::vector<std::string> list_visible_dirs(const std::string& d) {
    std::vector<std::string> names;
    std::error_code ec;
    for (const auto& e : fs::directory_iterator(fs::path(d), ec)) {
        std::error_code ec2;
        if (!e.is_directory(ec2)) continue;
        const std::string nm = e.path().filename().string();
        if (!nm.empty() && nm[0] == '.') continue;
        names.push_back(nm);
    }
    std::sort(names.begin(), names.end(), name_less);
    return names;
}

std::vector<std::string> list_visible_files(const std::string& d) {
    std::vector<std::string> names;
    std::error_code ec;
    for (const auto& e : fs::directory_iterator(fs::path(d), ec)) {
        std::error_code ec2;
        if (e.is_directory(ec2)) continue;
        const std::string nm = e.path().filename().string();
        if (!nm.empty() && nm[0] == '.') continue;
        names.push_back(nm);
    }
    std::sort(names.begin(), names.end(), name_less);
    return names;
}

// recursiveDir(d, 2): query path first, then BFS levels 1 and 2 (visible dirs).
std::vector<std::string> recursiveDir2(const std::string& d) {
    std::vector<std::string> p = {with_sep(d)};
    std::vector<std::string> files = p;
    for (int iter = 0; iter < 2 && !files.empty(); ++iter) {
        std::vector<std::string> next;
        for (const auto& q : files) {
            for (const auto& nm : list_visible_dirs(q)) next.push_back(q + nm + SEP);
        }
        p.insert(p.end(), next.begin(), next.end());
        files = next;
    }
    return p;
}

// getDirFromPath
void getDirFromPath(const std::string& dpath, std::string& dirName, std::string& dirPath) {
    std::vector<std::size_t> idx;
    for (std::size_t i = 0; i < dpath.size(); ++i) if (dpath[i] == SEP) idx.push_back(i);
    if (idx.empty()) { dirName = dpath; dirPath.clear(); return; }
    if (idx.back() == dpath.size() - 1) {
        if (idx.size() < 2) { dirName = dpath.substr(0, dpath.size() - 1); dirPath.clear(); return; }
        dirName = dpath.substr(idx[idx.size() - 2] + 1, idx.back() - idx[idx.size() - 2] - 1);
        dirPath = dpath.substr(0, idx[idx.size() - 2] + 1);
    } else {
        dirName = dpath.substr(idx.back() + 1);
        dirPath = dpath.substr(0, idx.back() + 1);
    }
}

} // namespace

std::vector<MovieData> loadConditionData(const std::string& condDirIn, const std::vector<std::string>& chNames,
                                         const std::vector<std::string>& markers, const ConditionOptions& opt) {
    if (chNames.empty()) throw std::invalid_argument("at least one channel name is required");
    const std::string condDir = with_sep(condDirIn);
    if (!fs::is_directory(condDir)) throw std::runtime_error("Condition directory not found: " + condDir);
    std::cout << "Root directory: " << condDir << "\n";

    std::vector<std::string> cellPath = recursiveDir2(condDir);
    std::vector<std::string> cellDirs(cellPath.size()), cellPar(cellPath.size());
    for (std::size_t i = 0; i < cellPath.size(); ++i) getDirFromPath(cellPath[i], cellDirs[i], cellPar[i]);

    // idx = regexpi(cellDirs, MovieSelector, 'once'); (StrictSelector: idx must be 1)
    const std::regex selRe(opt.movieSelector, std::regex::icase);
    std::vector<std::string> selPath, selPar;
    for (std::size_t i = 0; i < cellPath.size(); ++i) {
        std::smatch m;
        if (std::regex_search(cellDirs[i], m, selRe)) {
            if (opt.strictSelector && m.position(0) != 0) continue;
            selPath.push_back(cellPath[i]);
            selPar.push_back(cellPar[i]);
        }
    }
    if (selPath.empty()) throw std::runtime_error("No movies found in: " + condDir);

    // sort by cell number within each parent: sortStringsByToken(s, token, 'post')
    // q = (?<=token)\d+ ; entries without a match are removed; stable sort by number
    {
        std::vector<std::string> parents(selPar);
        std::sort(parents.begin(), parents.end());
        parents.erase(std::unique(parents.begin(), parents.end()), parents.end());
        std::vector<std::string> out(selPath.size());
        std::vector<char> valid(selPath.size(), 0);
        const std::regex tokRe("(" + opt.movieSelector + ")(\\d+)", std::regex::icase);
        for (const auto& par : parents) {
            std::vector<std::size_t> pos;
            for (std::size_t i = 0; i < selPath.size(); ++i) if (selPar[i] == par) pos.push_back(i);
            std::vector<std::pair<double, std::string>> items;
            for (std::size_t i : pos) {
                std::smatch m;
                // regexpi(s, q, 'match', 'once') on the full path
                if (std::regex_search(selPath[i], m, tokRe)) items.emplace_back(std::stod(m[2].str()), selPath[i]);
            }
            std::stable_sort(items.begin(), items.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
            for (std::size_t j = 0; j < items.size() && j < pos.size(); ++j) { out[pos[j]] = items[j].second; valid[pos[j]] = 1; }
        }
        std::vector<std::string> fin;
        for (std::size_t i = 0; i < out.size(); ++i) if (valid[i]) fin.push_back(out[i]);
        selPath = fin;
    }

    const std::size_t nCh = chNames.size();
    for (std::size_t c = 0; c < nCh; ++c) std::cout << "Channel " << (c + 1) << " name: \"" << chNames[c] << "\"\n";

    std::vector<MovieData> data;
    const std::regex dateRe("\\d{6}");
    const std::regex frRe("_(\\d+)?(\\.)?\\d+s");
    const std::regex msRe("_\\d+ms");
    const std::regex tifRe("\\.tif|\\.stk", std::regex::icase);
    const std::regex numRe("\\d+(?=\\.)");

    for (const auto& cp : selPath) {
        MovieData d;
        d.cellPath = cp;
        std::smatch m;
        // date: cell2mat of all 6-digit matches
        {
            std::string all;
            for (auto it = std::sregex_iterator(cp.begin(), cp.end(), dateRe); it != std::sregex_iterator(); ++it) all += it->str();
            d.date = all.empty() ? "000000" : all;
        }
        if (std::regex_search(cp, m, frRe)) {
            std::string s = m.str();
            d.framerate = std::stod(s.substr(1, s.size() - 2));
        } else if (std::regex_search(cp, m, msRe)) {
            std::string s = m.str();
            d.framerate = std::stod(s.substr(1, s.size() - 3)) / 1000.0;
        } else {
            d.framerate = 2.0;
        }

        d.channels.resize(nCh);
        d.framePaths.assign(nCh, {});
        bool allFound = true;
        for (std::size_t c = 0; c < nCh; ++c) {
            d.channels[c] = chNames[c].empty() ? cp : cp + chNames[c] + SEP;
            if (!fs::is_directory(d.channels[c])) {
                throw std::runtime_error("Channel directory not found (MATLAB would open a dialog): " + d.channels[c]);
            }
            std::vector<std::string> tmp;
            for (const auto& nm : list_visible_files(d.channels[c])) {
                if (std::regex_search(nm, tifRe)) tmp.push_back(nm);
            }
            // sort files in case leading zeros are missing: numbers immediately before a '.'
            {
                std::vector<std::string> nums(tmp.size());
                std::set<std::size_t> lens;
                for (std::size_t i = 0; i < tmp.size(); ++i) {
                    // emulate regexp(tmp, '\d+(?=\.)', 'match', 'once'): first digit run followed by '.'
                    std::string found;
                    for (std::size_t p = 0; p < tmp[i].size();) {
                        if (std::isdigit(static_cast<unsigned char>(tmp[i][p]))) {
                            std::size_t q = p;
                            while (q < tmp[i].size() && std::isdigit(static_cast<unsigned char>(tmp[i][q]))) ++q;
                            if (q < tmp[i].size() && tmp[i][q] == '.') { found = tmp[i].substr(p, q - p); break; }
                            // MATLAB regex would backtrack to shorter runs; a shorter run
                            // followed by '.' is impossible unless the char after is '.', so continue
                            p = q;
                        } else {
                            ++p;
                        }
                    }
                    nums[i] = found;
                    lens.insert(found.size());
                }
                if (lens.size() != 1) {
                    std::vector<std::pair<double, std::string>> items;
                    for (std::size_t i = 0; i < tmp.size(); ++i) {
                        double v = nums[i].empty() ? std::numeric_limits<double>::quiet_NaN() : std::stod(nums[i]);
                        items.emplace_back(v, tmp[i]);
                    }
                    // MATLAB sort puts NaN last, stable
                    std::stable_sort(items.begin(), items.end(), [](const auto& a, const auto& b) {
                        if (std::isnan(a.first)) return false;
                        if (std::isnan(b.first)) return true;
                        return a.first < b.first;
                    });
                    for (std::size_t i = 0; i < tmp.size(); ++i) tmp[i] = items[i].second;
                }
            }
            for (const auto& nm : tmp) d.framePaths[c].push_back(d.channels[c] + nm);
            if (tmp.empty()) allFound = false;
        }
        d.source = d.channels[0];
        if (allFound) {
            d.hasFrames = true;
            if (d.framePaths[0].size() == 1) {
                d.singleFile = true;
                for (std::size_t c = 0; c < nCh; ++c) {
                    if (d.framePaths[c].size() != 1) {
                        // MATLAB: cellfun(@(i) i{1}, framePaths) takes the first file of every channel
                        d.framePaths[c].resize(1);
                    }
                }
                TiffInfo info = tiff_info(d.framePaths[0][0]);
                d.imageHeight = info.height;
                d.imageWidth = info.width;
                d.movieLength = info.numDirectories;
            } else {
                d.singleFile = false;
                TiffInfo info = tiff_info(d.framePaths[0][0]);
                d.imageHeight = info.height;
                d.imageWidth = info.width;
                d.movieLength = d.framePaths[0].size();
            }
        } else {
            d.hasFrames = false;
            std::cerr << "Warning: not all channels contain TIFF frames in " << cp << "\n";
        }
        d.markers = markers;
        std::cout << "Loaded: " << cp << "\n";
        data.push_back(d);
    }
    if (opt.ignoreEmptyFolders) {
        std::vector<MovieData> keep;
        for (auto& d : data) if (d.hasFrames) keep.push_back(d);
        data.swap(keep);
    }
    return data;
}

FrameRef frame_ref(const MovieData& d, std::size_t c, std::size_t f) {
    if (d.singleFile) return {d.framePaths[c][0], f};
    return {d.framePaths[c][f - 1], 1};
}

} // namespace cme
