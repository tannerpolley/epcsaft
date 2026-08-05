#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"
bash .codex/environments/setup.sh

.venv/bin/python data/tools/validate_packets.py
(cd engine && ../.venv/bin/python -m pytest -q)
(cd engine && ../.venv/bin/ruff check --select E4,E7,E9,F src tests)
(cd engine && ../.venv/bin/mypy)
(cd validation && ../.venv/bin/python -m pytest -q -m 'not campaign')
(cd validation && ../.venv/bin/python -m pytest -q -m campaign)
(cd lab && ../.venv/bin/python -m pytest -q reference_oracles)
python3 governance/tools/management_preflight.py check
python3 governance/tools/artifact_store.py check
bash engine/tools/check-native.sh
bash engine/tools/check-native-sanitizers.sh
