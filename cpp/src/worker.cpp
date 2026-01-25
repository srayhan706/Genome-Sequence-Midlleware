#include "mpi_protocol.hpp"
#include "algorithms/smith_waterman.hpp"
#include "algorithms/kmp.hpp"
#include "algorithms/kmer.hpp"
#include <mpi.h>
#include <string>
#include <cstring>
#include <cstdlib>

namespace genomic {

void worker_loop(int rank, int world_size) {
  MPI_Status status;
  char* buf = nullptr;
  size_t buf_cap = 0;

  for (;;) {
    int32_t task_type = 0;
    MPI_Recv(&task_type, 1, MPI_INT, kRoot, kTagType, MPI_COMM_WORLD, &status);
    if (task_type == static_cast<int32_t>(TaskType::kShutdown)) {
      break;
    }

    int32_t payload_len = 0;
    MPI_Recv(&payload_len, 1, MPI_INT, kRoot, kTagLength, MPI_COMM_WORLD, &status);
    if (payload_len <= 0) continue;

    if (static_cast<size_t>(payload_len) > buf_cap) {
      buf_cap = static_cast<size_t>(payload_len) + 4096;
      buf = static_cast<char*>(std::realloc(buf, buf_cap));
    }
    MPI_Recv(buf, payload_len, MPI_BYTE, kRoot, kTagPayload, MPI_COMM_WORLD, &status);
    std::string payload(buf, buf + payload_len);

    std::string result;

    switch (static_cast<TaskType>(task_type)) {
      case TaskType::kAlign: {
        AlignParams params;
        if (!parse_align_request(payload, params)) break;
        int score = 0;
        auto [a, b] = smith_waterman(
            params.seq_a, params.seq_b,
            params.match_s, params.mismatch, params.gap,
            &score);
        result = pack_align_result(score, a, b);
        break;
      }
      case TaskType::kPattern: {
        PatternParams params;
        if (!parse_pattern_request(payload, params)) break;
        auto positions = kmp_search(params.text, params.pattern);
        result = pack_pattern_result(positions);
        break;
      }
      case TaskType::kKmer: {
        KmerParams params;
        if (!parse_kmer_request(payload, params)) break;
        auto counts = kmer_count(params.sequence, params.k);
        result = pack_kmer_result(counts);
        break;
      }
      default:
        break;
    }

    int32_t res_len = static_cast<int32_t>(result.size());
    MPI_Send(&res_len, 1, MPI_INT, kRoot, kTagResultLength, MPI_COMM_WORLD);
    if (res_len > 0) {
      MPI_Send(result.data(), res_len, MPI_BYTE, kRoot, kTagResult, MPI_COMM_WORLD);
    }
  }

  std::free(buf);
}

}  // namespace genomic
