#!/usr/bin/env bash
set -euo pipefail

engine_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
workspace_root="$(dirname "$engine_root")"
cd "$workspace_root"
if [[ "${EPCSAFT_RECORD_LOCK_HELD:-0}" != 1 ]]; then
    mkdir -p build/numerical-evidence
    exec {record_lock_fd}>build/numerical-evidence/.lock
    flock "$record_lock_fd"
    export EPCSAFT_RECORD_LOCK_HELD=1
fi

# Keep the supported direct native route on the measured fast native subset.
# Each registered case is run through the shared wheel/identity recorder;
# full equilibrium, boundary, and model-floor suites remain explicit owners.
if (($# == 0)); then
    names=(contributions properties diagnostics specification)
else
    names=("$@")
fi
for name in "${names[@]}"; do
    case "${name,,}" in
        contributions|properties|diagnostics|specification) ;;
        *)
            echo "unsupported direct native suite: $name (slow suites require an explicit owner command)" >&2
            exit 2
            ;;
    esac
done

EPCSAFT_PREPARE_NATIVE_BUILD=1 bash .codex/environments/setup.sh
for name in "${names[@]}"; do
    .venv/bin/python tools/run-numerical-evidence.py \
        --native-suite "$name" --explicit --no-build --timeout 600
done
