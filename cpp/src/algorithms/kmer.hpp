#pragma once

#include <string>
#include <unordered_map>
#include <cstdint>

namespace genomic {

// Largest k supported: k-mers are packed 2 bits per base into a uint64.
constexpr int kMaxKmerK = 32;

// Count k-mers in sequence (case-insensitive; k-mers containing a non-ACGT
// character are skipped). Keys are upper-case. Returns an empty map when
// k is outside [1, kMaxKmerK] or longer than the sequence.
// Large inputs are split across OpenMP threads by start position.
std::unordered_map<std::string, int64_t> kmer_count(const std::string& sequence, int k);

} // namespace genomic
