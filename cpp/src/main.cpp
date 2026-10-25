#include "dispatch.hpp"
#include "server.hpp"
#include "worker.hpp"
#include <mpi.h>
#include <unistd.h>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

struct RawNodeInfo {
  char hostname[64];
  int32_t omp_threads;
};

// Every rank reports its host and OpenMP width; rank 0 receives them all.
std::vector<genomic::NodeInfo> gather_node_info(int rank, int world_size) {
  RawNodeInfo mine{};
  gethostname(mine.hostname, sizeof(mine.hostname) - 1);
#ifdef _OPENMP
  mine.omp_threads = omp_get_max_threads();
#else
  mine.omp_threads = 1;
#endif
  std::vector<RawNodeInfo> all(rank == 0 ? world_size : 0);
  MPI_Gather(&mine, sizeof(RawNodeInfo), MPI_BYTE,
             all.data(), sizeof(RawNodeInfo), MPI_BYTE, 0, MPI_COMM_WORLD);

  std::vector<genomic::NodeInfo> nodes;
  for (int r = 0; r < static_cast<int>(all.size()); ++r) {
    nodes.push_back({r, all[r].hostname, all[r].omp_threads});
  }
  return nodes;
}

std::string parse_bind(int argc, char** argv) {
  std::string bind = "0.0.0.0:50051";
  for (int i = 1; i < argc; ++i) {
    if (std::strcmp(argv[i], "--bind") == 0 && i + 1 < argc) bind = argv[++i];
  }
  return bind;
}

}  // namespace

int main(int argc, char** argv) {
  // gRPC handlers run on many threads, each of which may talk to a different
  // worker at the same time, so ask for full thread support.
  int provided = MPI_THREAD_SINGLE;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  int rank = 0, world_size = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  auto nodes = gather_node_info(rank, world_size);

  int ret = 0;
  if (rank == 0) {
    bool multiple = provided >= MPI_THREAD_MULTIPLE;
    if (!multiple) {
      std::cerr << "warning: MPI_THREAD_MULTIPLE not available; "
                   "worker dispatch will be serialized" << std::endl;
    }
    genomic::dispatch_init(world_size, multiple, nodes);
    ret = genomic::run_grpc_server(parse_bind(argc, argv), multiple);
    genomic::dispatch_shutdown_workers();
  } else {
    genomic::worker_loop(rank);
  }

  MPI_Finalize();
  return ret;
}
