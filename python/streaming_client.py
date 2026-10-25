#!/usr/bin/env python3
"""Simulated streaming client: sends sequence chunks to the middleware as a
   simulated high-throughput sequencing pipeline (near real-time k-mer stats).
"""

import sys
import os
import argparse
import time
import random

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import grpc
import genomic_pb2
import genomic_pb2_grpc

BASES = "ACGT"


def random_chunk(size: int) -> str:
    return "".join(random.choices(BASES, k=size))


def main():
    ap = argparse.ArgumentParser(description="Streaming client (simulated sequencer)")
    ap.add_argument("--host", default="localhost")
    ap.add_argument("--port", type=int, default=50051)
    ap.add_argument("--chunks", type=int, default=10, help="Number of chunks to stream")
    ap.add_argument("--chunk-size", type=int, default=500, help="Bases per chunk")
    ap.add_argument("--k", type=int, default=4, help="K-mer size")
    ap.add_argument("--delay", type=float, default=0.2, help="Seconds between chunks")
    args = ap.parse_args()

    addr = f"{args.host}:{args.port}"
    channel = grpc.insecure_channel(addr)
    stub = genomic_pb2_grpc.GenomicServiceStub(channel)

    print(f"Streaming {args.chunks} chunks (size={args.chunk_size}, k={args.k}) to {addr}")
    total_kmers = 0
    cache_hits = 0
    for i in range(args.chunks):
        chunk = random_chunk(args.chunk_size)
        req = genomic_pb2.KmerCountRequest(sequence=chunk, k=args.k)
        t0 = time.perf_counter()
        r = stub.KmerCount(req)
        dt = time.perf_counter() - t0
        n = sum(r.counts.values())
        total_kmers += n
        if r.cache_hit:
            cache_hits += 1
        print(f"  Chunk {i+1}: {n} k-mers, {len(r.counts)} unique, cache_hit={r.cache_hit}, rank {r.info.worker_rank}, {dt*1000:.1f} ms")
        time.sleep(args.delay)

    print(f"Total k-mers: {total_kmers}, cache hits: {cache_hits}/{args.chunks}")


if __name__ == "__main__":
    main()
