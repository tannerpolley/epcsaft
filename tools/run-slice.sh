#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

if [[ "${1:-}" != "khudaida" ]]; then
    echo "usage: $0 khudaida" >&2
    exit 2
fi
bash .codex/environments/setup.sh
.venv/bin/python data/tools/validate_packets.py
candidate_dir="$(mktemp -d /tmp/epcsaft-khudaida.XXXXXX)"
qualification_dir="$(mktemp -d /tmp/epcsaft-khudaida-qualification.XXXXXX)"
if .venv/bin/python analyses/2026-khudaida-idaes-settler/scripts/run_fixed_settler.py \
    --output-dir "$candidate_dir" && \
    .venv/bin/python analyses/2026-khudaida-idaes-settler/scripts/qualify.py \
    --results-dir "$candidate_dir" --output-dir "$qualification_dir"; then
    find "$candidate_dir" "$qualification_dir" -depth -delete
else
    echo "Khudaida candidates retained at $candidate_dir and $qualification_dir" >&2
    exit 1
fi
