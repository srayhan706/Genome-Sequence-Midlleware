#pragma once

#include <string>
#include <vector>
#include <unordered_map>
#include <cstdint>

namespace genomic {

enum class TaskType : int32_t {
  kAlign = 1,
  kPattern = 2,
  kKmer = 3,
  kShutdown = 4,
};

struct AlignParams {
  std::string seq_a, seq_b;
  int32_t match_s = 2, mismatch = -1, gap = -1;
};

struct PatternParams {
  std::string text, pattern;
};

struct KmerParams {
  std::string sequence;
  int32_t k = 4;
};

// Serialize request; returns packed bytes (length-prefixed for MPI).
std::string pack_align_request(const AlignParams& p);
std::string pack_pattern_request(const PatternParams& p);
std::string pack_kmer_request(const KmerParams& p);

// Parse request payload (after receiving from MPI).
bool parse_align_request(const std::string& payload, AlignParams& out);
bool parse_pattern_request(const std::string& payload, PatternParams& out);
bool parse_kmer_request(const std::string& payload, KmerParams& out);

// Result serialization (worker -> middleware).
std::string pack_align_result(int score, const std::string& a, const std::string& b);
std::string pack_pattern_result(const std::vector<int>& positions);
std::string pack_kmer_result(const std::unordered_map<std::string, int64_t>& counts);

// Parse result (middleware side).
bool parse_align_result(const std::string& payload, int& score,
                        std::string& aligned_a, std::string& aligned_b);
bool parse_pattern_result(const std::string& payload, std::vector<int>& positions);
bool parse_kmer_result(const std::string& payload,
                       std::unordered_map<std::string, int64_t>& counts);

// MPI tags (middleware <-> workers)
constexpr int kTagType = 0;
constexpr int kTagLength = 1;
constexpr int kTagPayload = 2;
constexpr int kTagResultLength = 10;
constexpr int kTagResult = 11;
constexpr int kRoot = 0;

}  // namespace genomic
