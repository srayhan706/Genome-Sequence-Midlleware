#pragma once

#include "mpi_protocol.hpp"
#include <string>

namespace genomic {

// Run one task: parse payload, run the algorithm, pack the result.
// On failure returns false and writes an error message into *result.
// Shared by MPI workers and by rank 0 when running without workers.
bool execute_task(TaskType type, const std::string& payload, std::string* result);

// Worker main loop (ranks >= 1): receive tasks from rank 0 until kShutdown.
void worker_loop(int rank);

}  // namespace genomic
