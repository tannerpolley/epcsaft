#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$repo_root"

if ! command -v uv >/dev/null 2>&1; then
    echo "missing required command: uv" >&2
    exit 1
fi

uv sync --python 3.13 --all-packages --no-install-workspace --locked --inexact
bash tools/build-engine.sh
