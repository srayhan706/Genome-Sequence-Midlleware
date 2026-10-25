#pragma once

#include <string>

namespace genomic {

// Run gRPC server (blocks). Returns exit code.
// mpi_thread_multiple is reported through GetStats.
int run_grpc_server(const std::string& bind_address, bool mpi_thread_multiple);

}  // namespace genomic
