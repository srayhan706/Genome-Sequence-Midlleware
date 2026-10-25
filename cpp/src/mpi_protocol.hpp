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

const char* task_name(TaskType t);

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

// Serialize request payloads (middleware -> worker).
std::string pack_align_request(const AlignParams& p);
std::string pack_pattern_request(const PatternParams& p);
std::string pack_kmer_request(const KmerParams& p);

// Parse request payloads. Return false on truncated/malformed input.
bool parse_align_request(const std::string& payload, AlignParams& out);
bool parse_pattern_request(const std::string& payload, PatternParams& out);
bool parse_kmer_request(const std::string& payload, KmerParams& out);

// Result serialization (worker -> middleware).
std::string pack_align_result(int score, const std::string& a, const std::string& b);
std::string pack_pattern_result(const std::vector<int>& positions);
std::string pack_kmer_result(const std::unordered_map<std::string, int64_t>& counts);

// Parse result (middleware side). Return false on truncated/malformed input.
bool parse_align_result(const std::string& payload, int& score,
                        std::string& aligned_a, std::string& aligned_b);
bool parse_pattern_result(const std::string& payload, std::vector<int>& positions);
bool parse_kmer_result(const std::string& payload,
                       std::unordered_map<std::string, int64_t>& counts);

// MPI framing. Each direction is a single message whose size is discovered
// with MPI_Mprobe/MPI_Get_count:
//   task   (rank 0 -> worker, kTagTask):   [int32 TaskType][payload...]
//   result (worker -> rank 0, kTagResult): [int32 status][int64 compute_us][body...]
// status is 0 on success; otherwise body holds an error message.
constexpr int kTagTask = 1;
constexpr int kTagResult = 2;
constexpr int kRoot = 0;

constexpr size_t kTaskHeaderSize = 4;
constexpr size_t kResultHeaderSize = 12;

std::string frame_task(TaskType type, const std::string& payload);
bool unframe_task(const std::string& msg, TaskType* type, std::string* payload);

std::string frame_result(int32_t status, int64_t compute_us, const std::string& body);
bool unframe_result(const std::string& msg, int32_t* status, int64_t* compute_us,
                    std::string* body);

}  // namespace genomic
