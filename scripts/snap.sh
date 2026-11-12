#!/usr/bin/env bash
# Headless-Chrome screenshots of the viewer, one per (model, stage).
#   scripts/snap.sh OUTDIR "metro:1 metro:4 vase:8" [extra query, e.g. "&theme=dark"]
# Stage numbers are 1..8; frame defaults to the last iteration.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
OUTDIR="$1"; SHOTS="$2"; EXTRA="${3:-}"
PORT="${PORT:-8799}"
PY="${PYTHON:-$HOME/miniconda3/envs/scrapshell/bin/python}"
CHROME="${CHROME:-/Applications/Google Chrome.app/Contents/MacOS/Google Chrome}"
mkdir -p "$OUTDIR"
"$ROOT/scripts/build_web.sh" >/dev/null
"$PY" -m http.server "$PORT" -d "$ROOT/web" >/dev/null 2>&1 &
SERVER=$!
trap 'kill $SERVER 2>/dev/null || true' EXIT
sleep 1
PROFILE="$(mktemp -d)"
for shot in $SHOTS; do
  model="${shot%%:*}"; stage="${shot##*:}"
  "$CHROME" --headless=new --user-data-dir="$PROFILE" --hide-scrollbars \
    --use-angle=swiftshader --enable-unsafe-swiftshader --window-size=1440,900 \
    --virtual-time-budget=9000 --screenshot="$OUTDIR/${model}_s${stage}.png" \
    "http://localhost:$PORT/?model=$model&stage=$stage&frame=-1$EXTRA" >/dev/null 2>&1 || echo "failed: $shot"
done
rm -rf "$PROFILE"
ls "$OUTDIR"
