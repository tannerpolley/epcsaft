#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"
if [[ "${EPCSAFT_RECORD_LOCK_HELD:-0}" != 1 ]]; then
    mkdir -p build/numerical-evidence
    exec {record_lock_fd}>build/numerical-evidence/.lock
    flock "$record_lock_fd"
    export EPCSAFT_RECORD_LOCK_HELD=1
fi

run_recorded() {
    local test_id="$1"
    shift
    python3 tools/record-command.py --test-id "$test_id" -- "$@"
}

base_ref="${1:-HEAD~1}"
source_only="${EPCSAFT_CHECK_SOURCE_ONLY:-0}"
if [[ "$source_only" != 0 && "$source_only" != 1 ]]; then
    echo "EPCSAFT_CHECK_SOURCE_ONLY must be 0 or 1" >&2
    exit 2
fi
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
# The full gate owns every shared check; decide that before running any here.
if [[ "$source_only" == 0 ]] && grep -Eq '^(AGENTS.md|CONTEXT.md|README.md|pyproject.toml|uv.lock|\.codex/|tools/)' <<<"$paths"; then
    exec bash tools/check-all.sh
fi
run_recorded tooling.check-all.roadmap python3 tools/check-roadmap.py
run_recorded tooling.check-all.budgets python3 tools/diffstat.py --budgets
# Documentation readiness must not trigger wheel installation or numerical studies.
docs_only=true
for path in "${changed[@]}"; do
    case "$path" in
        docs/*.md|engine/docs/*.md|docs/*.json|docs/*.csv|engine/docs/equations_registry.yaml|AGENTS.md|engine/AGENTS.md|README.md|CONTEXT.md|engine/README.md|engine/CONTEXT.md|.gitattributes) ;;
        *) docs_only=false ;;
    esac
done
run_recorded tooling.check-all.equation-registry python3 engine/scripts/docs/sync_equation_registry.py --check
if [[ "$docs_only" == true ]]; then
    bash -n tools/check-affected.sh tools/check-all.sh
    echo "documentation readiness passed; no wheel build or numerical tests ran"
    exit 0
fi

architecture_inputs=false
prepare_native_build=0
if grep -Eq '^(analyses/|engine/)' <<<"$paths"; then
    run_recorded tooling.check-all.analysis-layout .venv/bin/python tools/check-analysis-layout.py
fi
if grep -Eq '^(engine/(ARCHITECTURE\.yaml|meson\.build|native/|src/|tests/native/)|docs/scientific/|tools/(check-engine-architecture\.py|verify-installed-wheel\.py))' <<<"$paths"; then
    architecture_inputs=true
    prepare_native_build=1
fi
if grep -q '^data/' <<<"$paths"; then
    run_recorded tooling.check-all.packets .venv/bin/python data/tools/validate_packets.py
fi
if [[ "$source_only" == 0 ]] && grep -Eq '^(engine/(ARCHITECTURE\.yaml|meson\.build|build/subprojects/|LICENSE|README\.md|pyproject\.toml|native/|src/|data/|tests/)|docs/scientific/|pyproject\.toml$|uv\.lock$|\.codex/)' <<<"$paths"; then
    run_recorded tooling.check-all.environment env EPCSAFT_PREPARE_NATIVE_BUILD="$prepare_native_build" bash .codex/environments/setup.sh
fi
if [[ "$source_only" == 0 && "$architecture_inputs" == true ]]; then
    build_dir="$(<build/environment-wheel/meson-build-dir)"
    mapfile -t engine_wheels < <(find build/environment-wheel -maxdepth 1 -type f -name 'epcsaft-*.whl' -print)
    if ((${#engine_wheels[@]} != 1)); then
        echo "expected exactly one Engine wheel, found ${#engine_wheels[@]}" >&2
        exit 1
    fi
    run_recorded tooling.check-all.architecture .venv/bin/python tools/check-engine-architecture.py --build-dir "$build_dir" --wheel "${engine_wheels[0]}"
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
            engine/native/equilibrium/*|engine/src/epcsaft/equilibrium/*)
                engine_equilibrium=true
                ;;
            engine/data/reference/*)
                engine_eos=true
                engine_equilibrium=true
                ;;
            engine/src/epcsaft/regression/*) engine_regression=true ;;
            engine/AGENTS.md|engine/ARCHITECTURE.yaml|engine/CONTEXT.md|engine/README.md|engine/meson.build|engine/build/subprojects/*|engine/pyproject.toml|engine/native/*|engine/src/*|engine/data/*|engine/tests/*)
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
    if [[ "$source_only" == 0 ]] && ((${#engine_test_paths[@]})); then
        run_recorded tooling.check-all.recorded-evidence .venv/bin/python tools/run-numerical-evidence.py --automatic --no-build
    fi
    if grep -Eq '^engine/(src/|tests/)' <<<"$paths"; then
        run_recorded tooling.check-all.ruff bash -c 'cd "$1/engine" && ../.venv/bin/ruff check --select E4,E7,E9,F src tests' _ "$repo_root"
    fi
    if grep -Eq '^engine/(src/|pyproject\.toml$)' <<<"$paths"; then
        run_recorded tooling.check-all.mypy bash -c 'cd "$1/engine" && ../.venv/bin/mypy' _ "$repo_root"
    fi
    if grep -q '^engine/tools/' <<<"$paths"; then
        bash -n engine/tools/*.sh
    fi
fi
if grep -q '^analyses/' <<<"$paths"; then
    echo "Analysis replay is explicit; run only the changed scientific claim or tools/run-slice.sh khudaida"
fi
if [[ "$source_only" == 1 ]]; then
    echo "source-only affected checks passed; no wheel build or installed-artifact tests ran"
fi
