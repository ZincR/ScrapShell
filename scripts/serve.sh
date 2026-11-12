#!/usr/bin/env bash
# Serve the viewer locally: http://localhost:${PORT:-8765}/
# (Browsers refuse to load web/data/*.json from a file:// page, so the viewer needs a server.)
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PY="${PYTHON:-$HOME/miniconda3/envs/scrapshell/bin/python}"
command -v "$PY" >/dev/null 2>&1 || PY=python3
"$ROOT/scripts/build_web.sh" >/dev/null
echo "Viewer: http://localhost:${PORT:-8765}/"
exec "$PY" -m http.server "${PORT:-8765}" -d "$ROOT/web"
