#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"
slice="${1:-}"
if [[ "$slice" != "khudaida" ]]; then
    echo "usage: $0 khudaida" >&2
    exit 2
fi

bash .codex/environments/setup.sh
.venv/bin/python data/tools/validate_packets.py
(cd validation && ../.venv/bin/python -m pytest -q -m campaign \
    campaigns/2026-khudaida/tests/test_phase_equilibrium_campaign.py)
