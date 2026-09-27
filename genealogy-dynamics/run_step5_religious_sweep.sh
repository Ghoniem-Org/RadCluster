#!/bin/bash
# Step-5 religious Case-1 sweep: shallowest/deepest x alpha 0.0/0.8/0.9.
# Usage: bash run_step5_religious_sweep.sh
set -e
cd "$(dirname "$0")"
OUTDIR="outputs/$(date +%Y%m%d_%H%M%S)_colonial_religious"
mkdir -p "$OUTDIR"
for RULE in shallowest_inheritance deepest_inheritance; do
  for A in 0.0 0.8 0.9; do
    echo "=== $RULE alpha=$A ==="
    GENEALOGY_MAX_DEPTH=14 python3 run_colonial_religious.py "$RULE" "$A" "$OUTDIR"
  done
done
echo "wrote $OUTDIR"
