#pragma once

#include <string>
#include <vector>

namespace genomic {

// KMP pattern search. Returns 0-based positions of all matches.
// Large texts are split into chunks of start positions searched in parallel
// with OpenMP (chunks overlap by len(pattern) - 1 so no match is missed).
std::vector<int> kmp_search(const std::string& text, const std::string& pattern);

} // namespace genomic
