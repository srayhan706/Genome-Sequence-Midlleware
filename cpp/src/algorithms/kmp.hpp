#pragma once

#include <string>
#include <vector>

namespace genomic {

// KMP pattern search. Returns 0-based positions of all matches.
// OpenMP parallelizes search over partitioned text chunks.
std::vector<int> kmp_search(const std::string& text, const std::string& pattern);

} // namespace genomic
