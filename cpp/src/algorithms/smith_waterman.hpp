#pragma once

#include <string>
#include <utility>

namespace genomic {

// Smith-Waterman local alignment. OpenMP used for parallel band computation
// if sequence length warrants it.
std::pair<std::string, std::string> smith_waterman(
    const std::string& a,
    const std::string& b,
    int match_score = 2,
    int mismatch_penalty = -1,
    int gap_penalty = -1,
    int* out_score = nullptr
);

} // namespace genomic
