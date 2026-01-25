#!/usr/bin/env python3
"""Visualize analysis results (k-mer distributions, alignment scores) using Matplotlib/Seaborn."""

import sys
import os
import argparse

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import grpc
import genomic_pb2
import genomic_pb2_grpc
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


def main():
    ap = argparse.ArgumentParser(description="Visualize genomic analysis results")
    ap.add_argument("--host", default="localhost")
    ap.add_argument("--port", type=int, default=50051)
    ap.add_argument("--sequence", default="ACGTACGTACGTACGTACGTACGTACGTACGT", help="Sequence for k-mer analysis")
    ap.add_argument("--k", type=int, default=3)
    ap.add_argument("--out", default="genomic_analysis.png", help="Output plot path")
    args = ap.parse_args()

    addr = f"{args.host}:{args.port}"
    channel = grpc.insecure_channel(addr)
    stub = genomic_pb2_grpc.GenomicServiceStub(channel)

    req = genomic_pb2.KmerCountRequest(sequence=args.sequence, k=args.k)
    r = stub.KmerCount(req)
    kmers = list(r.counts.keys())
    counts = list(r.counts.values())

    sns.set_theme(style="whitegrid")
    fig, ax = plt.subplots(figsize=(10, 5))
    x = np.arange(len(kmers))
    ax.bar(x, counts, color="steelblue", edgecolor="navy", alpha=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(kmers, rotation=45, ha="right")
    ax.set_xlabel(f"{args.k}-mer")
    ax.set_ylabel("Count")
    ax.set_title(f"K-mer distribution (k={args.k}, sequence length={len(args.sequence)})")
    fig.tight_layout()
    fig.savefig(args.out, dpi=150)
    print(f"Saved {args.out}")


if __name__ == "__main__":
    main()
