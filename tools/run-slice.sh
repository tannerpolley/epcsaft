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
    (cd engine && ../.venv/bin/python -m pytest -q -m campaign \
        tests/equilibrium/test_khudaida_held2_numerical_regression.py)
    (cd engine && ../.venv/bin/python -m pytest -q \
        tests/integration/test_scientific_slice.py)
else
    (cd engine && ../.venv/bin/python -m pytest -q -m campaign \
        tests/regression/test_figiel_aqueous_kij.py \
        tests/regression/test_figiel_water_factor.py \
        tests/regression/test_born_tracer.py)
    (cd validation && ../.venv/bin/python -m pytest -q \
        campaigns/2025-figiel/tests/test_figiel_regression_campaigns.py)
fi
