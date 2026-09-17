#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

if [[ "${1:-}" != "khudaida" ]]; then
    echo "usage: $0 khudaida" >&2
    exit 2
fi
EPCSAFT_PREPARE_NATIVE_BUILD=1 bash .codex/environments/setup.sh
.venv/bin/python data/tools/validate_packets.py
candidate_dir="$(mktemp -d /tmp/epcsaft-khudaida.XXXXXX)"
if .venv/bin/python analyses/2026-khudaida-idaes-settler/scripts/run_fixed_settler.py \
    --output-dir "$candidate_dir"; then
    case "$candidate_dir" in
        /tmp/epcsaft-khudaida.*) rm -rf -- "$candidate_dir" ;;
        *) echo "refusing to remove unexpected candidate path: $candidate_dir" >&2; exit 1 ;;
    esac
else
    echo "Khudaida candidate retained at $candidate_dir" >&2
    exit 1
fi
