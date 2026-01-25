# Middleware for Distributed Genomic Sequence Analysis

**CSE-4269 Parallel and Distributed Systems Lab — Project (Roll 41)**

This project implements a middleware-based distributed system for genomic sequence analysis.

## Quick start (how to run)

```bash
# 1. Build C++ middleware (requires: cmake, MPI, grpc, protobuf, pkg-config)
cd genomic-middleware && mkdir -p build && cd build
cmake .. && cmake --build .

# 2. Python client: venv + deps + generate gRPC stubs
cd .. && python3 -m venv .venv && source .venv/bin/activate
pip install -r python/requirements.txt
python -m grpc_tools.protoc -Iproto --python_out=python --grpc_python_out=python proto/genomic.proto

# 3. Start middleware (e.g. 1 coordinator + 2 workers)
mpirun -np 3 ./build/cpp/genomic_middleware --bind 0.0.0.0:50051
# (Leave running. In another terminal, from genomic-middleware/:)

# 4. Run client (use project venv; base Anaconda may have Protobuf 5 vs generated 6)
cd genomic-middleware && source .venv/bin/activate
python python/client.py --health
python python/client.py --align ACGTACGT ACGTTGCA
python python/client.py --kmer ACGTACGTACGTACGT 3
```

The middleware coordinates between client applications and a cluster of worker nodes, handling job submission, data distribution, and result aggregation via **RPC** (gRPC), **MPI**, and **OpenMP**.

## Features

- **RPC (gRPC)**: Client–middleware communication for sequence alignment, pattern search, and k-mer counting.
- **MPI**: Task distribution across worker processes; result aggregation via scatter/gather.
- **OpenMP**: Intra-node parallelism in workers (KMP, k-mer) where applicable.
- **Caching**: Repeated queries served from cache; `cache_hit` reported in responses.
- **Multi-client**: Job scheduling with round-robin dispatch across workers.
- **Streaming**: Simulated sequencing pipeline (`streaming_client.py`) sending chunks for near real-time k-mer stats.

## Project Structure

```
genomic-middleware/
├── proto/
│   └── genomic.proto          # gRPC service definitions
├── cpp/
│   ├── CMakeLists.txt
│   └── src/
│       ├── main.cpp           # MPI rank 0 = middleware, rank ≥1 = workers
│       ├── server.cpp/hpp     # gRPC server, cache, dispatch
│       ├── worker.cpp/hpp     # MPI worker loop (algorithms)
│       ├── dispatch.cpp/hpp   # Middleware→worker MPI dispatch
│       ├── mpi_protocol.*     # Task/result serialization
│       └── algorithms/        # Smith–Waterman, KMP, k-mer (OpenMP)
├── python/
│   ├── client.py              # gRPC client (align, pattern, kmer, health)
│   ├── streaming_client.py    # Simulated streaming k-mer pipeline
│   ├── visualize.py           # Matplotlib/Seaborn charts
│   ├── genomic_pb2.py         # Generated (see below)
│   └── genomic_pb2_grpc.py
├── data/
│   └── sample.fasta
├── scripts/
│   ├── run_middleware.sh
│   └── generate_grpc_python.sh
└── README.md
```

## Requirements

- **C++ build**: CMake ≥3.16, C++17, MPI (OpenMPI/MPICH), Protobuf, gRPC, pkg-config.
- **Python client**: Python 3.9+, `grpcio`, `grpcio-tools`, `protobuf`, `matplotlib`, `seaborn`, `numpy`.

### macOS (Homebrew)

```bash
brew install cmake open-mpi grpc protobuf pkg-config
```

Optional (for OpenMP in workers):

```bash
brew install libomp
```

## How to Build

### 1. C++ middleware (MPI + gRPC)

```bash
cd genomic-middleware
mkdir -p build && cd build
cmake ..
cmake --build .
```

The executable is `build/cpp/genomic_middleware`.

### 2. Python stubs and venv

```bash
cd genomic-middleware
python3 -m venv .venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -r python/requirements.txt
python -m grpc_tools.protoc -Iproto --python_out=python --grpc_python_out=python proto/genomic.proto
```

## How to Run

### 1. Start the middleware

Run with MPI. Use **at least 1** process (middleware only, no workers; compute runs locally). For workers, use **≥2** (e.g. 1 middleware + 2 workers).

**Single process (no workers, local compute):**

```bash
cd genomic-middleware/build
mpirun -np 1 ./cpp/genomic_middleware --bind 0.0.0.0:50051
```

**With workers (1 middleware + 2 workers):**

```bash
mpirun -np 3 ./cpp/genomic_middleware --bind 0.0.0.0:50051
```

Or use the helper script (default: 3 processes, bind `0.0.0.0:50051`):

```bash
./scripts/run_middleware.sh [np] [--bind host:port]
# e.g. ./scripts/run_middleware.sh 5 --bind 0.0.0.0:50051
```

You should see:

```
Genomic middleware listening on 0.0.0.0:50051 (workers: 2)
```

### 2. Run the Python client

Activate the venv, then:

```bash
source .venv/bin/activate
cd python
```

**Health check:**

```bash
python client.py --health
# Health: ok=True, workers=2
```

**Sequence alignment (Smith–Waterman):**

```bash
python client.py --align ACGTACGT ACGTTGCA
```

**Pattern search (KMP):**

```bash
python client.py --pattern "ACGTACGTACGTACGT" ACGT
```

**K-mer count:**

```bash
python client.py --kmer ACGTACGTACGTACGT 3
```

**Streaming (simulated sequencer):**

```bash
python streaming_client.py --chunks 10 --chunk-size 500 --k 4 --delay 0.2
```

**Visualization (k-mer bar chart):**

```bash
python visualize.py --sequence "ACGTACGTACGTACGTACGTACGTACGTACGT" --k 3 --out genomic_analysis.png
```

Use `--host` and `--port` if the middleware is not on `localhost:50051`:

```bash
python client.py --host 192.168.1.10 --port 50051 --health
```

## Batch alignment

The gRPC API supports `BatchAlignment`: either a `reference` sequence plus a list of `sequences` (each aligned to `reference`), or pairs from `sequences` (0–1, 2–3, …). Use the generated `genomic_pb2` / `genomic_pb2_grpc` modules from your own Python (or other) client to call it.

## Sample data

- `data/sample.fasta`: small FASTA for manual testing. Parse locally and send sequences to the client (e.g. `--align`, `--kmer`).

## Performance notes

- **OpenMP**: Used in KMP (large text) and k-mer counting when available. Build with `libomp` and a compiler that supports OpenMP for speedups.
- **Workers**: More MPI processes (workers) improve throughput for batch-style workloads.
- **Caching**: Repeated identical requests return `cache_hit=True` and skip recomputation.

## Troubleshooting

**`VersionError: Detected mismatched Protobuf Gencode/Runtime major versions`**  
The Python stubs were generated with Protobuf 6.x; your environment (e.g. Anaconda base) has Protobuf 5.x.

- **Fix:** Use the project venv (recommended): `source .venv/bin/activate`, then run the client. The venv has matching `grpcio`/`protobuf`.
- **Alternative:** In your current env, upgrade Protobuf: `pip install 'protobuf>=6,<7'` (may conflict with other packages).

## License

For use in CSE-4269 course context only.
