#include "dispatch.hpp"
#include "server.hpp"
#include "worker.hpp"
#include <mpi.h>
#include <cstdlib>
#include <iostream>
#include <string>

namespace {

int run_middleware(int argc, char** argv);
int run_worker(int rank, int world_size);

}  // namespace

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  int rank = 0, world_size = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  int ret = 0;
  if (rank == 0) {
    ret = run_middleware(argc, argv);
    genomic::dispatch_shutdown_workers();
  } else {
    ret = run_worker(rank, world_size);
  }

  MPI_Finalize();
  return ret;
}

namespace {

int run_middleware(int argc, char** argv) {
  int rank = 0, world_size = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  genomic::dispatch_init(rank, world_size);

  std::string bind = "0.0.0.0:50051";
  for (int i = 1; i < argc; ++i) {
    std::string a = argv[i];
    if (a == "--bind" && i + 1 < argc) {
      bind = argv[++i];
      break;
    }
  }

  return genomic::run_grpc_server(bind);
}

int run_worker(int rank, int world_size) {
  genomic::worker_loop(rank, world_size);
  return 0;
}

}  // namespace
