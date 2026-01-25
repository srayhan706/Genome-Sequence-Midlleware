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


def main():
    ap = argparse.ArgumentParser(description="Genomic Middleware client")
    ap.add_argument("--host", default="localhost", help="Middleware host")
    ap.add_argument("--port", type=int, default=50051, help="Middleware port")
    ap.add_argument("--health", action="store_true", help="Health check only")
    ap.add_argument("--align", nargs=2, metavar=("A", "B"), help="Sequence alignment (Smith-Waterman)")
    ap.add_argument("--pattern", nargs=2, metavar=("TEXT", "PATTERN"), help="Pattern search (KMP)")
    ap.add_argument("--kmer", nargs=2, metavar=("SEQUENCE", "K"), help="K-mer count (k as int)")
    args = ap.parse_args()

    addr = f"{args.host}:{args.port}"
    channel = grpc.insecure_channel(addr)
    stub = genomic_pb2_grpc.GenomicServiceStub(channel)

    try:
        if args.health:
            r = stub.HealthCheck(genomic_pb2.HealthRequest())
            print(f"Health: ok={r.ok}, workers={r.worker_count}")
            return

        if args.align:
            a, b = args.align
            req = genomic_pb2.AlignmentRequest(
                sequence_a=a, sequence_b=b,
                match_score=2, mismatch_penalty=-1, gap_penalty=-1,
            )
            r = stub.SequenceAlignment(req)
            print(f"Alignment (cache_hit={r.cache_hit}):")
            print(f"  Score: {r.score}")
            print(f"  A: {r.aligned_a}")
            print(f"  B: {r.aligned_b}")
            return

        if args.pattern:
            text, pat = args.pattern
            req = genomic_pb2.PatternSearchRequest(text=text, pattern=pat)
            r = stub.PatternSearch(req)
            print(f"Pattern search (cache_hit={r.cache_hit}): count={r.count}")
            print(f"  Positions: {list(r.positions)}")
            return

        if args.kmer:
            seq, k = args.kmer
            req = genomic_pb2.KmerCountRequest(sequence=seq, k=int(k))
            r = stub.KmerCount(req)
            print(f"K-mer count (k={k}, cache_hit={r.cache_hit}):")
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
