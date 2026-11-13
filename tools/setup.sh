#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

uv sync --python 3.13 --all-packages --no-install-workspace --locked --inexact
bash data/tools/restore_feos.sh
.venv/bin/python data/tools/validate_packets.py
bash tools/build-engine.sh "$@"
