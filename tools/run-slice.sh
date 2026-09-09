#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"
if [[ "${EPCSAFT_RECORD_LOCK_HELD:-0}" != 1 ]]; then
    mkdir -p build/numerical-evidence
    exec {record_lock_fd}>build/numerical-evidence/.lock
    flock "$record_lock_fd"
    export EPCSAFT_RECORD_LOCK_HELD=1
fi
slice="${1:-}"
if [[ "$slice" != "khudaida" ]]; then
    echo "usage: $0 khudaida" >&2
    exit 2
fi

EPCSAFT_PREPARE_NATIVE_BUILD=1 bash .codex/environments/setup.sh
.venv/bin/python data/tools/validate_packets.py
.venv/bin/python tools/run-numerical-evidence.py \
    study.reproduction.2026-khudaida --no-build --timeout 600
