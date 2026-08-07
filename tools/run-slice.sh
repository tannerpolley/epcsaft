#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"
slice="${1:-}"
if [[ "$slice" != "khudaida" && "$slice" != "figiel" ]]; then
    echo "usage: $0 khudaida|figiel" >&2
    exit 2
fi

bash .codex/environments/setup.sh
.venv/bin/python data/tools/validate_packets.py
if [[ "$slice" == "khudaida" ]]; then
    (cd validation && ../.venv/bin/python -m pytest -q -m campaign \
        campaigns/2026-khudaida/tests/test_khudaida_phase_equilibrium.py)
else
    (cd engine && ../.venv/bin/python -m pytest -q \
        tests/regression/numerical/test_generic_fit.py)
fi
