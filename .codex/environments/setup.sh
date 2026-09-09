#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$repo_root"

mkdir -p build/numerical-evidence
if [[ "${EPCSAFT_RECORD_LOCK_HELD:-0}" != 1 ]]; then
    exec {record_lock_fd}>build/numerical-evidence/.lock
    flock "$record_lock_fd"
    export EPCSAFT_RECORD_LOCK_HELD=1
fi

if ! command -v uv >/dev/null 2>&1; then
    echo "missing required command: uv" >&2
    exit 1
fi

uv sync --python 3.13 --all-packages --no-install-workspace --locked
bash tools/ensure-engine-wheel.sh
