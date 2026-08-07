#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$repo_root"

for required_command in git uv cmake ninja c++ flock; do
    if ! command -v "$required_command" >/dev/null 2>&1; then
        echo "missing required command: $required_command" >&2
        exit 1
    fi
done

uv sync --python 3.13 --no-install-project --locked --inexact
bash tools/ensure-engine-wheel.sh
