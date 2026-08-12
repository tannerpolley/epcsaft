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
if grep -q '^data/' <<<"$paths"; then
    python3 data/tools/validate_packets.py
fi
if grep -Eq '^(engine/(CMakeLists\.txt|LICENSE|README\.md|pyproject\.toml|cpp/|src/|data/|tests/)|validation/|pyproject\.toml$|uv\.lock$|\.codex/)' <<<"$paths"; then
    bash .codex/environments/setup.sh
fi
if grep -q '^engine/' <<<"$paths"; then
    engine_all=false
    engine_eos=false
    engine_equilibrium=false
    engine_regression=false
    while IFS= read -r path; do
        case "$path" in
            engine/tests/eos/*) engine_eos=true ;;
            engine/tests/equilibrium/*) engine_equilibrium=true ;;
            engine/tests/regression/*) engine_regression=true ;;
            engine/cpp/equilibrium/*|engine/src/epcsaft/equilibrium/*)
                engine_equilibrium=true
                ;;
            engine/data/reference/*)
                engine_eos=true
                engine_equilibrium=true
                ;;
            engine/src/epcsaft/regression/*) engine_regression=true ;;
            engine/AGENTS.md|engine/ARCHITECTURE.yaml|engine/CONTEXT.md|engine/README.md|engine/CMakeLists.txt|engine/pyproject.toml|engine/cpp/*|engine/src/*|engine/data/*|engine/tests/*)
                engine_all=true
                ;;
        esac
    done <<<"$paths"

    engine_test_paths=()
    if [[ "$engine_all" == true ]]; then
        engine_test_paths=(tests)
    else
        [[ "$engine_eos" == true ]] && engine_test_paths+=(tests/eos)
        [[ "$engine_equilibrium" == true ]] && engine_test_paths+=(tests/equilibrium)
        [[ "$engine_regression" == true ]] && engine_test_paths+=(tests/regression)
    fi
    if ((${#engine_test_paths[@]})); then
        (cd engine && ../.venv/bin/python -m pytest -q "${engine_test_paths[@]}")
    fi
    if grep -Eq '^engine/(src/|tests/)' <<<"$paths"; then
        (cd engine && ../.venv/bin/ruff check --select E4,E7,E9,F src tests)
    fi
    if grep -Eq '^engine/(src/|pyproject\.toml$)' <<<"$paths"; then
        (cd engine && ../.venv/bin/mypy)
    fi
    if grep -q '^engine/tools/' <<<"$paths"; then
        bash -n engine/tools/*.sh
    fi
fi
if grep -q '^validation/' <<<"$paths"; then
    (cd validation && ../.venv/bin/python -m pytest -q -m 'not campaign')
    if grep -q '^validation/campaigns/' <<<"$paths"; then
        echo "Validation campaign replay is explicit; run the changed campaign or tools/run-slice.sh when its scientific claim requires it"
    fi
fi
if grep -q '^governance/' <<<"$paths"; then
    python3 governance/tools/management_preflight.py check
    python3 governance/tools/artifact_store.py check
fi
