#!/usr/bin/env python3
"""HTTP/JSON gateway + static web dashboard for the Genomic Middleware.

Browsers can't speak native gRPC, so this process sits between them:

    browser --HTTP/JSON--> gateway.py --gRPC--> middleware (MPI rank 0) --MPI--> workers

Standard library only (plus the grpc stubs in ../python).

    python web/gateway.py --grpc localhost:50051 --port 8080
"""

import argparse
import json
import os
import random
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "python"))

import grpc  # noqa: E402
import genomic_pb2 as pb  # noqa: E402
import genomic_pb2_grpc as pb_grpc  # noqa: E402

STATIC_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "static")
SAMPLE_FASTA = os.path.join(ROOT, "data", "sample.fasta")
MAX_MESSAGE = 64 * 1024 * 1024
MAX_POSITIONS = 5000   # pattern-search positions returned to the browser
MAX_LOADTEST_JOBS = 256

_stub = None
_target = None


class ApiError(Exception):
    def __init__(self, status, message):
        super().__init__(message)
        self.status = status
        self.message = message


def stub():
    return _stub


def job_info(info):
    return {
        "worker": info.worker_rank,
        "compute_ms": round(info.compute_ms, 3),
        "total_ms": round(info.total_ms, 3),
    }


def parse_fasta(text):
    records, name, buf = [], None, []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                records.append({"name": name, "sequence": "".join(buf)})
            name, buf = line[1:].strip() or f"seq{len(records) + 1}", []
        else:
            buf.append(line)
    if name is not None:
        records.append({"name": name, "sequence": "".join(buf)})
    return records


def random_dna(n, rnd=random):
    return "".join(rnd.choice("ACGT") for _ in range(n))


def align_request(body):
    return pb.AlignmentRequest(
        sequence_a=str(body.get("a", "")),
        sequence_b=str(body.get("b", "")),
        match_score=int(body.get("match", 2)),
        mismatch_penalty=int(body.get("mismatch", -1)),
        gap_penalty=int(body.get("gap", -1)),
    )


def align_result(r):
    return {
        "score": r.score,
        "aligned_a": r.aligned_a,
        "aligned_b": r.aligned_b,
        "cache_hit": r.cache_hit,
        **job_info(r.info),
    }


# ---- API handlers: body dict -> response dict ------------------------------

def api_health(_):
    r = stub().HealthCheck(pb.HealthRequest(), timeout=2)
    return {"ok": r.ok, "workers": r.worker_count, "target": _target}


def api_stats(_):
    r = stub().GetStats(pb.StatsRequest(), timeout=2)
    return {
        "target": _target,
        "workers": [
            {
                "rank": w.rank,
                "hostname": w.hostname,
                "omp_threads": w.omp_threads,
                "jobs": w.jobs,
                "failures": w.failures,
                "active": w.active,
                "current_task": w.current_task,
                "compute_ms": round(w.compute_ms, 2),
                "bytes_sent": w.bytes_sent,
                "bytes_received": w.bytes_received,
            }
            for w in r.workers
        ],
        "requests": r.requests,
        "in_flight": r.in_flight,
        "cache": {
            "hits": r.cache_hits,
            "misses": r.cache_misses,
            "entries": r.cache_entries,
            "bytes": r.cache_bytes,
            "capacity": r.cache_capacity_bytes,
        },
        "uptime_s": round(r.uptime_s, 1),
        "openmp": r.openmp,
        "mpi_thread_multiple": r.mpi_thread_multiple,
        "coordinator_host": r.coordinator_host,
    }


def api_align(body):
    return align_result(stub().SequenceAlignment(align_request(body)))


def api_pattern(body):
    r = stub().PatternSearch(pb.PatternSearchRequest(
        text=str(body.get("text", "")), pattern=str(body.get("pattern", ""))))
    positions = list(r.positions)
    return {
        "count": r.count,
        "positions": positions[:MAX_POSITIONS],
        "truncated": len(positions) > MAX_POSITIONS,
        "cache_hit": r.cache_hit,
        **job_info(r.info),
    }


def api_kmer(body):
    k = int(body.get("k", 3))
    top = max(1, min(int(body.get("top", 40)), 5000))
    r = stub().KmerCount(pb.KmerCountRequest(sequence=str(body.get("sequence", "")), k=k))
    items = sorted(r.counts.items(), key=lambda kv: (-kv[1], kv[0]))
    return {
        "k": k,
        "distinct": len(items),
        "total": sum(r.counts.values()),
        "top": [{"kmer": km, "count": c} for km, c in items[:top]],
        "cache_hit": r.cache_hit,
        **job_info(r.info),
    }


def api_batch(body):
    seqs = [str(s) for s in body.get("sequences", [])]
    req = pb.BatchAlignmentRequest(
        sequences=seqs,
        reference=str(body.get("reference", "")),
        match_score=int(body.get("match", 2)),
        mismatch_penalty=int(body.get("mismatch", -1)),
        gap_penalty=int(body.get("gap", -1)),
    )
    r = stub().BatchAlignment(req)
    return {"wall_ms": round(r.wall_ms, 2), "results": [align_result(x) for x in r.results]}


def api_loadtest(body):
    """Fire many concurrent requests (unique inputs, so no cache hits) and
    record when each one started/finished and which worker ran it."""
    op = body.get("op", "align")
    n = max(1, min(int(body.get("requests", 32)), MAX_LOADTEST_JOBS))
    concurrency = max(1, min(int(body.get("concurrency", 8)), 64))
    length = max(10, min(int(body.get("length", 1500)), 4000 if op == "align" else 2_000_000))
    rnd = random.Random()

    if op == "align":
        inputs = [(random_dna(length, rnd), random_dna(length, rnd)) for _ in range(n)]
        call = lambda x: stub().SequenceAlignment(pb.AlignmentRequest(
            sequence_a=x[0], sequence_b=x[1], match_score=2, mismatch_penalty=-1, gap_penalty=-1))
    elif op == "kmer":
        inputs = [random_dna(length, rnd) for _ in range(n)]
        call = lambda x: stub().KmerCount(pb.KmerCountRequest(sequence=x, k=8))
    elif op == "pattern":
        inputs = [random_dna(length, rnd) for _ in range(n)]
        call = lambda x: stub().PatternSearch(pb.PatternSearchRequest(text=x, pattern="ACGTAC"))
    else:
        raise ApiError(400, f"unknown op {op!r}")

    t0 = time.perf_counter()

    def run(i):
        start = time.perf_counter()
        r = call(inputs[i])
        end = time.perf_counter()
        return {
            "i": i,
            "start_ms": round((start - t0) * 1000, 2),
            "end_ms": round((end - t0) * 1000, 2),
            "cache_hit": r.cache_hit,
            **job_info(r.info),
        }

    with ThreadPoolExecutor(concurrency) as ex:
        jobs = list(ex.map(run, range(n)))
    wall_ms = (time.perf_counter() - t0) * 1000
    return {
        "op": op,
        "requests": n,
        "concurrency": concurrency,
        "length": length,
        "wall_ms": round(wall_ms, 2),
        "throughput": round(n / (wall_ms / 1000), 2),
        "jobs": jobs,
    }


def api_sample(_):
    try:
        with open(SAMPLE_FASTA) as f:
            return {"records": parse_fasta(f.read())}
    except OSError:
        return {"records": []}


def api_fasta(body):
    return {"records": parse_fasta(str(body.get("text", "")))}


ROUTES = {
    ("GET", "/api/health"): api_health,
    ("GET", "/api/stats"): api_stats,
    ("GET", "/api/sample"): api_sample,
    ("POST", "/api/align"): api_align,
    ("POST", "/api/pattern"): api_pattern,
    ("POST", "/api/kmer"): api_kmer,
    ("POST", "/api/batch"): api_batch,
    ("POST", "/api/loadtest"): api_loadtest,
    ("POST", "/api/fasta"): api_fasta,
}

GRPC_TO_HTTP = {
    grpc.StatusCode.INVALID_ARGUMENT: 400,
    grpc.StatusCode.UNAVAILABLE: 503,
    grpc.StatusCode.DEADLINE_EXCEEDED: 504,
}


class Handler(SimpleHTTPRequestHandler):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, directory=STATIC_DIR, **kwargs)

    def log_message(self, fmt, *args):
        if not self.path.startswith("/api/stats"):
            sys.stderr.write("%s %s\n" % (self.log_date_time_string(), fmt % args))

    def do_GET(self):
        if self.path.split("?")[0].startswith("/api/"):
            self.handle_api("GET")
        else:
            super().do_GET()

    def do_POST(self):
        self.handle_api("POST")

    def handle_api(self, method):
        path = self.path.split("?")[0]
        handler = ROUTES.get((method, path))
        try:
            if handler is None:
                raise ApiError(404, f"no route {method} {path}")
            body = {}
            if method == "POST":
                length = int(self.headers.get("Content-Length") or 0)
                if length:
                    body = json.loads(self.rfile.read(length))
            self.send_json(200, handler(body))
        except ApiError as e:
            self.send_json(e.status, {"error": e.message})
        except grpc.RpcError as e:
            status = GRPC_TO_HTTP.get(e.code(), 502)
            msg = e.details() or str(e.code())
            if e.code() == grpc.StatusCode.UNAVAILABLE:
                msg = f"middleware not reachable at {_target}"
            self.send_json(status, {"error": msg, "grpc_code": e.code().name})
        except (ValueError, TypeError) as e:
            self.send_json(400, {"error": f"bad request: {e}"})

    def send_json(self, status, obj):
        data = json.dumps(obj).encode()
        self.send_response(status)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(data)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(data)


def main():
    global _stub, _target
    ap = argparse.ArgumentParser(description="Web dashboard for the Genomic Middleware")
    ap.add_argument("--grpc", default="localhost:50051", help="middleware gRPC address")
    ap.add_argument("--host", default="127.0.0.1", help="HTTP bind host")
    ap.add_argument("--port", type=int, default=8080, help="HTTP port")
    args = ap.parse_args()

    _target = args.grpc
    channel = grpc.insecure_channel(_target, options=[
        ("grpc.max_send_message_length", MAX_MESSAGE),
        ("grpc.max_receive_message_length", MAX_MESSAGE),
    ])
    _stub = pb_grpc.GenomicServiceStub(channel)

    server = ThreadingHTTPServer((args.host, args.port), Handler)
    server.daemon_threads = True
    print(f"Dashboard: http://{args.host}:{args.port}  (middleware: {_target})")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass


if __name__ == "__main__":
    main()
