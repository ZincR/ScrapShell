#!/usr/bin/env bash
# Build tools/phcpf_mesher natively (arm64) against the `scrapshell` conda env (GMP + CMake).
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
ENV_PREFIX="${CONDA_ENV_PREFIX:-$HOME/miniconda3/envs/scrapshell}"
git -C "$ROOT" submodule update --init --depth 1 third_party/Directional
"$ENV_PREFIX/bin/cmake" -S "$ROOT/tools/phcpf_mesher" -B "$ROOT/build/phcpf_mesher" \
  -G Ninja -DCMAKE_MAKE_PROGRAM="$ENV_PREFIX/bin/ninja" \
  -DCMAKE_PREFIX_PATH="$ENV_PREFIX" -DCMAKE_BUILD_TYPE=Release
"$ENV_PREFIX/bin/cmake" --build "$ROOT/build/phcpf_mesher"
echo "built: $ROOT/build/phcpf_mesher/phcpf_mesher"
