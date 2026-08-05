#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"
base_ref="${1:-HEAD~1}"
if ! git rev-parse --verify "$base_ref" >/dev/null 2>&1; then
    base_ref="HEAD"
fi
mapfile -t changed < <(
    {
        git diff --name-only "$base_ref" HEAD
        git diff --name-only
        git diff --cached --name-only
        git ls-files --others --exclude-standard
    } | sort -u
)
if ((${#changed[@]} == 0)); then
    echo "no changed paths"
    exit 0
fi

paths="$(printf '%s\n' "${changed[@]}")"
if grep -Eq '^(AGENTS.md|CONTEXT.md|README.md|pyproject.toml|uv.lock|\.codex/|tools/)' <<<"$paths"; then
    exec bash tools/check-all.sh
fi
bash .codex/environments/setup.sh
if grep -q '^data/' <<<"$paths"; then
    .venv/bin/python data/tools/validate_packets.py
fi
if grep -q '^engine/' <<<"$paths"; then
    (cd engine && ../.venv/bin/python -m pytest -q)
    (cd engine && ../.venv/bin/ruff check --select E4,E7,E9,F src tests)
    (cd engine && ../.venv/bin/mypy)
fi
if grep -q '^validation/' <<<"$paths"; then
    (cd validation && ../.venv/bin/python -m pytest -q -m 'not campaign')
    if grep -q '^validation/campaigns/' <<<"$paths"; then
        (cd validation && ../.venv/bin/python -m pytest -q -m campaign)
    fi
fi
if grep -q '^lab/' <<<"$paths"; then
    (cd lab && ../.venv/bin/python -m pytest -q reference_oracles)
fi
if grep -q '^governance/' <<<"$paths"; then
    python3 governance/tools/management_preflight.py check
    python3 governance/tools/artifact_store.py check
fi
