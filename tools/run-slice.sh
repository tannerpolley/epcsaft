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
(cd analyses/reproduction/2026-khudaida && ../../../.venv/bin/python -m pytest -q tests/test_phase_equilibrium_campaign.py)
