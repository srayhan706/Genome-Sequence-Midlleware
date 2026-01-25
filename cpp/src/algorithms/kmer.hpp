#pragma once

#include <string>
#include <unordered_map>
#include <cstdint>

namespace genomic {

// Count k-mers in sequence. OpenMP parallel reduction over chunks.
std::unordered_map<std::string, int64_t> kmer_count(const std::string& sequence, int k);

} // namespace genomic
