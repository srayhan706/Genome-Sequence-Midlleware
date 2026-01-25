#pragma once

#include "mpi_protocol.hpp"
#include <string>
#include <vector>
#include <unordered_map>

namespace genomic {

// Call after MPI_Init. Uses world_size and rank.
void dispatch_init(int rank, int world_size);

int dispatch_worker_count();

// Run alignment on a worker (or locally if single process). Returns true on success.
bool dispatch_align(const AlignParams& params, int* score,
                    std::string* aligned_a, std::string* aligned_b);

bool dispatch_pattern(const PatternParams& params, std::vector<int>* positions);

bool dispatch_kmer(const KmerParams& params,
                   std::unordered_map<std::string, int64_t>* counts);

// Signal workers to shutdown (middleware only).
void dispatch_shutdown_workers();

}  // namespace genomic
