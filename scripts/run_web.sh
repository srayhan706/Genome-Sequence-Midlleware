#!/bin/bash
# Start the web dashboard (HTTP gateway -> gRPC middleware).
# Usage: ./run_web.sh [grpc_target] [http_port]
# Default: grpc_target=localhost:50051, http_port=8080
set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
PY="python3"
[[ -x "$ROOT/.venv/bin/python" ]] && PY="$ROOT/.venv/bin/python"
exec "$PY" "$ROOT/web/gateway.py" --grpc "${1:-localhost:50051}" --port "${2:-8080}"
