#pragma once

#include "mpi_protocol.hpp"
#include <string>
#include <vector>

namespace genomic {

// Static facts about one MPI rank, gathered to rank 0 at startup.
struct NodeInfo {
  int rank = 0;
  std::string hostname;
  int omp_threads = 1;
};

// Where and how long a dispatched task ran.
struct DispatchInfo {
  int worker_rank = 0;     // 0 = ran locally on the middleware (no workers)
  double compute_ms = 0;   // time inside the algorithm on the worker
  double total_ms = 0;     // queueing + MPI round trip + compute
};

struct WorkerSnapshot {
  NodeInfo node;
  int64_t jobs = 0;
  int64_t failures = 0;
  int active = 0;
  std::string current_task;
  double compute_ms = 0;
  int64_t bytes_sent = 0;      // middleware -> worker
  int64_t bytes_received = 0;  // worker -> middleware
};

// Call on rank 0 after MPI_Init_thread. nodes is indexed by rank.
// thread_multiple: MPI granted MPI_THREAD_MULTIPLE, so different gRPC threads
// may talk to different workers concurrently; otherwise MPI use is serialized.
void dispatch_init(int world_size, bool thread_multiple, std::vector<NodeInfo> nodes);

int dispatch_worker_count();

// Run a task on the next idle worker (blocking until one is free), or locally
// when there are no workers. On success *result holds the packed result.
// On failure returns false and sets *error.
bool dispatch_task(TaskType type, const std::string& payload,
                   std::string* result, DispatchInfo* info, std::string* error);

// Per-rank statistics. With no workers, contains only rank 0.
std::vector<WorkerSnapshot> dispatch_snapshot();

// Signal workers to shut down (middleware only).
void dispatch_shutdown_workers();

}  // namespace genomic
