#!/usr/bin/env bash
# Run the PH-CPF baseline on every example model (or the ones given) and refresh web/data.
#   scripts/run_all.sh                 # all 8 models
#   scripts/run_all.sh metro vase      # a subset
# A model that fails the CPF stage is retried with iter_opt = 1..3 (continuity weight doubled each
# warm-up pass), as the official runme.m does after a barrier failure.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
MATLAB="${MATLAB:-/Applications/MATLAB_R2026a.app/bin/matlab}"
[ -x "$ROOT/build/phcpf_mesher/phcpf_mesher" ] || "$ROOT/scripts/build_mesher.sh"
MODELS=("$@")
if [ ${#MODELS[@]} -eq 0 ]; then
  MODELS=(metro intersection six_s0 vase 3holes_r car_r sphere2con_s0 c_s0)
fi
mkdir -p "$ROOT/results"
for id in "${MODELS[@]}"; do
  echo "=== $id ==="
  "$MATLAB" -batch "cd('$ROOT'); addpath('matlab'); run_with_retries('$id');" 2>&1 \
    | grep -vE "^\[?Warning: Name is nonexistent|rename the function|^\s*$" | tee "$ROOT/results/$id.log"
done
