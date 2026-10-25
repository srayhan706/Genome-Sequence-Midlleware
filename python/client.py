#!/usr/bin/env python3
"""Python client for the Genomic Middleware (gRPC)."""

import sys
import os
import argparse

# Allow imports from same dir (genomic_pb2, genomic_pb2_grpc)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import grpc
import genomic_pb2
import genomic_pb2_grpc


def describe(r):
    """Where a response was computed: cache, or which MPI rank and how long."""
    if r.cache_hit:
        return "(cache hit)"
    return f"(rank {r.info.worker_rank}, compute {r.info.compute_ms:.2f} ms)"


def main():
    ap = argparse.ArgumentParser(description="Genomic Middleware client")
    ap.add_argument("--host", default="localhost", help="Middleware host")
    ap.add_argument("--port", type=int, default=50051, help="Middleware port")
    ap.add_argument("--health", action="store_true", help="Health check only")
    ap.add_argument("--align", nargs=2, metavar=("A", "B"), help="Sequence alignment (Smith-Waterman)")
    ap.add_argument("--pattern", nargs=2, metavar=("TEXT", "PATTERN"), help="Pattern search (KMP)")
    ap.add_argument("--kmer", nargs=2, metavar=("SEQUENCE", "K"), help="K-mer count (k as int)")
    ap.add_argument("--batch", nargs="+", metavar=("REF", "SEQ"),
                    help="Batch alignment: align each SEQ against REF in parallel across workers")
    ap.add_argument("--stats", action="store_true", help="Per-worker and cache statistics")
    args = ap.parse_args()

    addr = f"{args.host}:{args.port}"
    channel = grpc.insecure_channel(addr)
    stub = genomic_pb2_grpc.GenomicServiceStub(channel)

    try:
        if args.health:
            r = stub.HealthCheck(genomic_pb2.HealthRequest())
            print(f"Health: ok={r.ok}, workers={r.worker_count}")
            return

        if args.stats:
            r = stub.GetStats(genomic_pb2.StatsRequest())
            print(f"Coordinator {r.coordinator_host}: up {r.uptime_s:.0f}s, requests={r.requests}, "
                  f"in_flight={r.in_flight}, openmp={r.openmp}, mpi_thread_multiple={r.mpi_thread_multiple}")
            lookups = r.cache_hits + r.cache_misses
            rate = f"{100 * r.cache_hits / lookups:.1f}%" if lookups else "-"
            print(f"Cache: {r.cache_entries} entries, {r.cache_bytes / 1024:.1f} KB, hit rate {rate}")
            for w in r.workers:
                print(f"  rank {w.rank} ({w.hostname}, {w.omp_threads} OpenMP threads): "
                      f"jobs={w.jobs}, compute={w.compute_ms:.1f} ms, "
                      f"MPI out/in={w.bytes_sent}/{w.bytes_received} B"
                      + (f", running {w.current_task}" if w.active else ""))
            return

        if args.batch:
            if len(args.batch) < 2:
                ap.error("--batch needs REF and at least one SEQ")
            ref, *seqs = args.batch
            r = stub.BatchAlignment(genomic_pb2.BatchAlignmentRequest(
                reference=ref, sequences=seqs,
                match_score=2, mismatch_penalty=-1, gap_penalty=-1,
            ))
            print(f"Batch of {len(r.results)} alignments in {r.wall_ms:.1f} ms:")
            for i, x in enumerate(r.results):
                print(f"  #{i + 1}: score={x.score} {describe(x)}")
            return

        if args.align:
            a, b = args.align
            req = genomic_pb2.AlignmentRequest(
                sequence_a=a, sequence_b=b,
                match_score=2, mismatch_penalty=-1, gap_penalty=-1,
            )
            r = stub.SequenceAlignment(req)
            print(f"Alignment {describe(r)}:")
            print(f"  Score: {r.score}")
            print(f"  A: {r.aligned_a}")
            print(f"  B: {r.aligned_b}")
            return

        if args.pattern:
            text, pat = args.pattern
            req = genomic_pb2.PatternSearchRequest(text=text, pattern=pat)
            r = stub.PatternSearch(req)
            print(f"Pattern search {describe(r)}: count={r.count}")
            print(f"  Positions: {list(r.positions)}")
            return

        if args.kmer:
            seq, k = args.kmer
            req = genomic_pb2.KmerCountRequest(sequence=seq, k=int(k))
            r = stub.KmerCount(req)
            print(f"K-mer count k={k} {describe(r)}:")
            for kmer, cnt in sorted(r.counts.items(), key=lambda x: -x[1])[:20]:
                print(f"  {kmer}: {cnt}")
            if len(r.counts) > 20:
                print(f"  ... and {len(r.counts) - 20} more")
            return

        ap.print_help()
        sys.exit(0)
    except grpc.RpcError as e:
        print(f"RPC error: {e.code()} - {e.details()}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
