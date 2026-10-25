#include "worker.hpp"
#include "algorithms/smith_waterman.hpp"
#include "algorithms/kmp.hpp"
#include "algorithms/kmer.hpp"
#include <mpi.h>
#include <chrono>
#include <string>

namespace genomic {

bool execute_task(TaskType type, const std::string& payload, std::string* result) {
  switch (type) {
    case TaskType::kAlign: {
      AlignParams params;
      if (!parse_align_request(payload, params)) break;
      int score = 0;
      auto [a, b] = smith_waterman(
          params.seq_a, params.seq_b,
          params.match_s, params.mismatch, params.gap,
          &score);
      *result = pack_align_result(score, a, b);
      return true;
    }
    case TaskType::kPattern: {
      PatternParams params;
      if (!parse_pattern_request(payload, params)) break;
      *result = pack_pattern_result(kmp_search(params.text, params.pattern));
      return true;
    }
    case TaskType::kKmer: {
      KmerParams params;
      if (!parse_kmer_request(payload, params)) break;
      *result = pack_kmer_result(kmer_count(params.sequence, params.k));
      return true;
    }
    default:
      *result = "unknown task type";
      return false;
  }
  *result = std::string("malformed ") + task_name(type) + " request";
  return false;
}

void worker_loop(int /*rank*/) {
  std::string msg;
  std::string payload;
  for (;;) {
    // Probe first so the receive buffer can be sized to the incoming task.
    MPI_Message handle;
    MPI_Status status;
    MPI_Mprobe(kRoot, kTagTask, MPI_COMM_WORLD, &handle, &status);
    int len = 0;
    MPI_Get_count(&status, MPI_BYTE, &len);
    msg.resize(static_cast<size_t>(len));
    MPI_Mrecv(msg.data(), len, MPI_BYTE, &handle, MPI_STATUS_IGNORE);

    TaskType type{};
    if (!unframe_task(msg, &type, &payload)) continue;
    if (type == TaskType::kShutdown) break;

    auto t0 = std::chrono::steady_clock::now();
    std::string body;
    bool ok = execute_task(type, payload, &body);
    int64_t us = std::chrono::duration_cast<std::chrono::microseconds>(
        std::chrono::steady_clock::now() - t0).count();

    std::string reply = frame_result(ok ? 0 : 1, us, body);
    MPI_Send(reply.data(), static_cast<int>(reply.size()), MPI_BYTE,
             kRoot, kTagResult, MPI_COMM_WORLD);
  }
}

}  // namespace genomic
