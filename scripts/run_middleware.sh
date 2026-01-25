#!/bin/bash
# Run the genomic middleware with MPI.
# Usage: ./run_middleware.sh [np] [--bind host:port]
# Default: np=3 (1 middleware + 2 workers), bind=0.0.0.0:50051
set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
BUILD="${BUILD:-$ROOT/build}"
EXE="$BUILD/cpp/genomic_middleware"
if [[ ! -f "$EXE" ]]; then
  echo "Build first: cd $ROOT && mkdir -p build && cd build && cmake .. && cmake --build ."
  exit 1
fi
NP="${1:-3}"
shift || true
ARGS=("$@")
if [[ ${#ARGS[@]} -eq 0 ]]; then
  ARGS=(--bind "0.0.0.0:50051")
fi
echo "Running: mpirun -np $NP $EXE ${ARGS[*]}"
exec mpirun -np "$NP" "$EXE" "${ARGS[@]}"
