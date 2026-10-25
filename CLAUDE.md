# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

Course project (CSE-4269 Parallel & Distributed Systems): a C++ middleware that exposes genomic algorithms (Smith–Waterman, KMP pattern search, k-mer counting) over gRPC, distributes jobs to MPI worker processes, and uses OpenMP inside workers. Python CLI clients and a web dashboard (`web/`) talk to it.

## Commands

Dependencies (macOS): `brew install cmake open-mpi grpc protobuf pkg-config libomp`.

```bash
# Build (defaults to Release) -> build/cpp/genomic_middleware, build/cpp/algorithms_test
cmake -S . -B build && cmake --build build -j

# Tests (algorithms vs brute force + MPI wire format; no MPI/gRPC needed)
./build/cpp/algorithms_test          # or: ctest --test-dir build

# Python env + stubs (regenerate stubs whenever proto/genomic.proto changes)
python3 -m venv .venv && source .venv/bin/activate
pip install -r python/requirements.txt
./scripts/generate_grpc_python.sh

# Run: np=1 -> coordinator only (computes locally); np>=2 -> 1 coordinator + (np-1) workers
mpirun -np 5 ./build/cpp/genomic_middleware --bind 0.0.0.0:50051
./scripts/run_web.sh localhost:50051 8080          # dashboard at http://127.0.0.1:8080
python python/client.py --health | --stats | --align A B | --pattern TEXT PAT | --kmer SEQ K | --batch REF SEQ...
```

There is no linter. Verify end-to-end changes by running the middleware and exercising it with `client.py` or the dashboard.

## Architecture

**Single binary, role by MPI rank** (`cpp/src/main.cpp`): `MPI_Init_thread(MPI_THREAD_MULTIPLE)`, every rank's hostname/OpenMP width is `MPI_Gather`ed to rank 0, then rank 0 runs the gRPC server (`server.cpp`) and ranks ≥1 run `worker_loop` (`worker.cpp`).

**Request flow**: gRPC handler validates/normalizes input (strip whitespace, upper-case) → packs request with `pack_*` → `run_job` checks the LRU cache (key = task type byte + packed payload) → `dispatch_task` (`dispatch.cpp`) → worker → result cached and parsed into the proto response along with `JobInfo{worker_rank, compute_ms, total_ms}`.

- **Dispatch**: a FIFO pool of idle worker ranks guarded by a mutex + condvar. A gRPC thread takes one worker, does a blocking send/recv with it, then returns it. The mutex is never held across MPI calls, so concurrent requests run on different workers in parallel. If MPI lacks `THREAD_MULTIPLE`, `g_mpi_serial` serializes the round trips. With one rank, `execute_task` runs locally on rank 0.
- **MPI wire format** (`mpi_protocol.hpp`): one message per direction, sized with `MPI_Mprobe`/`MPI_Get_count`/`MPI_Mrecv`. Task `[int32 TaskType][payload]` on `kTagTask`; result `[int32 status][int64 compute_us][body]` on `kTagResult` (non-zero status → body is an error string). Payloads are custom little-endian binary (`pack_*`/`parse_*`), not protobuf; parsers return false on truncation.
- **`execute_task`** (`worker.cpp`) is the one place that maps a `TaskType` to an algorithm; both workers and local mode use it.
- **Adding an operation**: new `TaskType` + pack/parse pair in `mpi_protocol.*`, a case in `execute_task`, proto RPC/messages, a handler in `GenomicServiceImpl` using `run_job`, and (optionally) a route in `web/gateway.py` + UI tab.
- **`BatchAlignment`** spawns one thread per worker that pulls jobs from an atomic index, so every worker is busy; results are returned in input order.
- **`GetStats`** reads `dispatch_snapshot()` (per-rank jobs, active, compute ms, MPI bytes) plus cache counters; the dashboard polls it every second.
- **CMake targets**: `genomic_core` (algorithms + `mpi_protocol`, no MPI/gRPC, links OpenMP), `genomic_proto` (generated into `build/cpp/gen/`), `genomic_middleware`, `algorithms_test`.

**Algorithms** (`cpp/src/algorithms/`) each have a serial path plus an `#ifdef _OPENMP` parallel path above a size threshold; tests force both paths and compare against brute force:
- Smith–Waterman: 256×256 tiles, anti-diagonal wavefront; best cell picked in a row-major scan afterwards so serial and parallel results are identical.
- KMP: blocks of start positions per thread, each scanning `m-1` chars past its block.
- k-mer: 2-bit packed rolling `uint64` (hence k ≤ 32, `kMaxKmerK`), split by start position.

**Web** (`web/`): `gateway.py` is stdlib-only (`ThreadingHTTPServer`) and imports stubs from `python/`. It serves `web/static/` (vanilla JS + hand-written SVG charts, no CDN) and exposes `/api/{health,stats,align,pattern,kmer,batch,loadtest,sample,fasta}`. gRPC errors map to HTTP (INVALID_ARGUMENT→400, UNAVAILABLE→503).

## Gotchas

- Python stubs (`python/genomic_pb2*.py`) are checked in and must be regenerated after proto edits; C++ stubs are generated at build time.
- `VersionError: mismatched Protobuf Gencode/Runtime` means the stubs and the installed protobuf runtime are different major versions. Use `.venv`, or regenerate the stubs.
- Server input limits live at the top of `server.cpp` (`kMaxAlignCells`, `kMaxSequenceLen`, message size, cache budget). Smith–Waterman memory is O(M·N), which is why alignment is capped by cell count.
- On one machine, MPI ranks × OpenMP threads can oversubscribe cores; use `mpirun -x OMP_NUM_THREADS=…`.
