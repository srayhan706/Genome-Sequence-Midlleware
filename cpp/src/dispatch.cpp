#include "dispatch.hpp"
#include "mpi_protocol.hpp"
#include "algorithms/smith_waterman.hpp"
#include "algorithms/kmp.hpp"
#include "algorithms/kmer.hpp"
#include <mpi.h>
#include <mutex>
#include <cstring>
#include <cstdlib>

namespace genomic {

static int g_rank = 0;
static int g_world_size = 0;
static int g_next_worker = 1;
static std::mutex g_dispatch_mutex;

void dispatch_init(int rank, int world_size) {
  g_rank = rank;
  g_world_size = world_size;
  g_next_worker = (world_size > 1) ? 1 : -1;
}

int dispatch_worker_count() {
  return (g_world_size > 1) ? (g_world_size - 1) : 0;
}

static int pick_worker() {
  if (g_world_size <= 1) return -1;
  int w = g_next_worker;
  g_next_worker = (g_next_worker + 1) % g_world_size;
  if (g_next_worker == 0) g_next_worker = 1;
  return w;
}

bool dispatch_align(const AlignParams& params, int* score,
                    std::string* aligned_a, std::string* aligned_b) {
  std::lock_guard<std::mutex> lock(g_dispatch_mutex);
  if (g_world_size > 1) {
    int w = pick_worker();
    std::string payload = pack_align_request(params);
    int32_t type = static_cast<int32_t>(TaskType::kAlign);
    int32_t len = static_cast<int32_t>(payload.size());
    MPI_Send(&type, 1, MPI_INT, w, kTagType, MPI_COMM_WORLD);
    MPI_Send(&len, 1, MPI_INT, w, kTagLength, MPI_COMM_WORLD);
    MPI_Send(payload.data(), len, MPI_BYTE, w, kTagPayload, MPI_COMM_WORLD);

    int32_t res_len = 0;
    MPI_Recv(&res_len, 1, MPI_INT, w, kTagResultLength, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    if (res_len <= 0) return false;
    std::string result(static_cast<size_t>(res_len), '\0');
    MPI_Recv(&result[0], res_len, MPI_BYTE, w, kTagResult, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    return parse_align_result(result, *score, *aligned_a, *aligned_b);
  }

  auto [a, b] = smith_waterman(
      params.seq_a, params.seq_b,
      params.match_s, params.mismatch, params.gap,
      score);
  *aligned_a = a;
  *aligned_b = b;
  return true;
}

bool dispatch_pattern(const PatternParams& params, std::vector<int>* positions) {
  std::lock_guard<std::mutex> lock(g_dispatch_mutex);
  if (g_world_size > 1) {
    int w = pick_worker();
    std::string payload = pack_pattern_request(params);
    int32_t type = static_cast<int32_t>(TaskType::kPattern);
    int32_t len = static_cast<int32_t>(payload.size());
    MPI_Send(&type, 1, MPI_INT, w, kTagType, MPI_COMM_WORLD);
    MPI_Send(&len, 1, MPI_INT, w, kTagLength, MPI_COMM_WORLD);
    MPI_Send(payload.data(), len, MPI_BYTE, w, kTagPayload, MPI_COMM_WORLD);

    int32_t res_len = 0;
    MPI_Recv(&res_len, 1, MPI_INT, w, kTagResultLength, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    if (res_len <= 0) return false;
    std::string result(static_cast<size_t>(res_len), '\0');
    MPI_Recv(&result[0], res_len, MPI_BYTE, w, kTagResult, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    return parse_pattern_result(result, *positions);
  }

  *positions = kmp_search(params.text, params.pattern);
  return true;
}

bool dispatch_kmer(const KmerParams& params,
                   std::unordered_map<std::string, int64_t>* counts) {
  std::lock_guard<std::mutex> lock(g_dispatch_mutex);
  if (g_world_size > 1) {
    int w = pick_worker();
    std::string payload = pack_kmer_request(params);
    int32_t type = static_cast<int32_t>(TaskType::kKmer);
    int32_t len = static_cast<int32_t>(payload.size());
    MPI_Send(&type, 1, MPI_INT, w, kTagType, MPI_COMM_WORLD);
    MPI_Send(&len, 1, MPI_INT, w, kTagLength, MPI_COMM_WORLD);
    MPI_Send(payload.data(), len, MPI_BYTE, w, kTagPayload, MPI_COMM_WORLD);

    int32_t res_len = 0;
    MPI_Recv(&res_len, 1, MPI_INT, w, kTagResultLength, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    if (res_len <= 0) return false;
    std::string result(static_cast<size_t>(res_len), '\0');
    MPI_Recv(&result[0], res_len, MPI_BYTE, w, kTagResult, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    return parse_kmer_result(result, *counts);
  }

  *counts = kmer_count(params.sequence, params.k);
  return true;
}

void dispatch_shutdown_workers() {
  int32_t type = static_cast<int32_t>(TaskType::kShutdown);
  for (int w = 1; w < g_world_size; ++w) {
    MPI_Send(&type, 1, MPI_INT, w, kTagType, MPI_COMM_WORLD);
  }
}

}  // namespace genomic
