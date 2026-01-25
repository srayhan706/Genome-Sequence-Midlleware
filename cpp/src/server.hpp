#pragma once

#include <string>

namespace genomic {

// Run gRPC server (blocks). Returns exit code.
int run_grpc_server(const std::string& bind_address);

}  // namespace genomic
