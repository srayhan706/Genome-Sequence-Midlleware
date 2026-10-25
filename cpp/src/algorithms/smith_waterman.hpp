#pragma once

#include <string>
#include <utility>

namespace genomic {

// Smith-Waterman local alignment (linear gap penalty). For large inputs the
// DP matrix is filled as a tiled anti-diagonal wavefront with OpenMP.
// Memory is O(len(a) * len(b)).
std::pair<std::string, std::string> smith_waterman(
    const std::string& a,
    const std::string& b,
    int match_score = 2,
    int mismatch_penalty = -1,
    int gap_penalty = -1,
    int* out_score = nullptr
);

} // namespace genomic
