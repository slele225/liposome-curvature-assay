// Parallelisation level of cme_detect.
//
// Either level gives identical output: every unit of work (a movie, a frame,
// a candidate fit, an image column) writes only to its own slots, no
// floating-point reduction is split across threads, and the RNG-consuming
// GMM step is always sequential.
#pragma once

namespace cme {

enum class ParLevel {
    Candidate,   // outer loops sequential; threads over candidate fits and image columns (default)
    Movie,       // threads over movies / frames / sigma-estimation images; inner work sequential (legacy)
};

} // namespace cme
