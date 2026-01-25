#!/bin/bash
# Generate Python gRPC stubs from proto
set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
PROTO_DIR="$ROOT/proto"
OUT_DIR="$ROOT/python"
mkdir -p "$OUT_DIR"
python -m grpc_tools.protoc \
  -I"$PROTO_DIR" \
  --python_out="$OUT_DIR" \
  --grpc_python_out="$OUT_DIR" \
  "$PROTO_DIR/genomic.proto"
echo "Generated $OUT_DIR/genomic_pb2.py $OUT_DIR/genomic_pb2_grpc.py"
