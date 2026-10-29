#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

# Select from every change since the base: by default the merge base with origin/main.
base_ref="${1:-$(git merge-base origin/main HEAD 2>/dev/null || git merge-base main HEAD 2>/dev/null || true)}"
if ! git rev-parse --verify "$base_ref" >/dev/null 2>&1; then
    [[ -z "${1:-}" ]] || { echo "check-affected: base '$1' does not resolve to a commit" >&2; exit 2; }
    base_ref="HEAD"
fi
base="$base_ref ($(git rev-parse --short "$base_ref"))"
# Runs one step and prints its elapsed time, keeping the step's exit status.
timed() {
    local name="$1" start="$EPOCHREALTIME" status=0
    shift
    "$@" || status=$?
    awk -v n="$name" -v s="$start" -v e="$EPOCHREALTIME" 'BEGIN { printf "[%s: %.1f s]\n", n, e - s }'
    return "$status"
}
echo "check-affected base: $base"
mapfile -t changed < <(
    {
        git diff --name-only "$base_ref" HEAD
        git diff --name-only
        git diff --cached --name-only
        git ls-files --others --exclude-standard
    } | sort -u
)
if ((${#changed[@]} == 0)); then
    echo "affected checks passed: no changed paths since base $base"
    exit 0
fi

notebooks=()
for path in "${changed[@]}"; do
    [[ "$path" == analyses/*/notebook.qmd || "$path" == analyses/*/README.md ]] && notebooks+=("$path")
done
((${#notebooks[@]} == 0)) || timed notebook-numbers python3 tools/check-notebook-numbers.py --base "$base_ref" "${notebooks[@]}"

paths="$(printf '%s\n' "${changed[@]}")"
timed project-map python3 tools/project-map.py
if grep -Eq '^(docs/|engine/docs/|AGENTS.md|README.md|ARCHITECTURE\.yaml)' <<<"$paths"; then
    timed equation-registry python3 engine/scripts/docs/sync_equation_registry.py --check
fi
if grep -Eq '^(analyses/|docs/|engine/|tools/|ARCHITECTURE\.yaml)' <<<"$paths"; then
    timed analysis-layout python3 tools/check-analysis-layout.py
fi

declarations_affected=false
analysis_checkers='^tools/(check-analysis-declarations\.py|analysis-declarations-excluded\.tsv|check-notebook-numbers\.py|tests/test-analysis-checks\.py)'
if grep -Eq "^(analyses/|data/|engine/|pyproject\.toml$|uv\.lock$|\.codex/)|$analysis_checkers" <<<"$paths"; then
    declarations_affected=true
fi
check_declarations() {
    if [[ "$declarations_affected" == true ]]; then
        if grep -Eq "$analysis_checkers" <<<"$paths"; then
            # Runs the full declaration check, then the checkers' known-failure examples.
            timed declarations .venv/bin/python tools/tests/test-analysis-checks.py
        else
            timed declarations .venv/bin/python tools/check-analysis-declarations.py
        fi
    fi
}

check_tools() {
    grep -Eq '^(engine/)?tools/' <<<"$paths" || return 0
    for script in $(grep -ls '^#!/usr/bin/env bash' tools/* engine/tools/*); do bash -n "$script"; done
    if grep -q '^tools/project-map\.py' <<<"$paths"; then python3 tools/project-map.py --self-test; fi
    if grep -Eq '^tools/(time-key-solves|tests/test-time-key-solves-gate\.py)$' <<<"$paths"; then
        python3 tools/tests/test-time-key-solves-gate.py
    fi
}

docs_only=true
while IFS= read -r path; do
    case "$path" in
        docs/*|engine/docs/*|AGENTS.md|README.md|engine/README.md|.gitattributes) ;;
        *) docs_only=false ;;
    esac
done <<<"$paths"
if [[ "$docs_only" == true ]]; then
    for script in tools/check-affected.sh tools/check-all.sh; do bash -n "$script"; done
    echo "affected checks passed: documentation readiness only, no Engine build or numerical tests, base $base"
    exit 0
fi

slow_tests=()
if grep -Eq '^engine/(native/|tests/native/|build/subprojects/|meson\.build)' <<<"$paths"; then
    # The slow solved-state programs run when their layer or any layer below it changes
    # (owners from ARCHITECTURE.yaml): a property or EOS change moves every solved state.
    owners="$(printf '%s\n' "${changed[@]}" | PYTHONPATH=tools python3 -c '
import importlib, sys
m = importlib.import_module("project-map")
rules = m.declaration()["rules"]
print(*{m.resolve(p, rules)[0] for p in sys.stdin.read().splitlines()})')"
    # A program's *_checks_shared.hpp reaches only that program, so the header rules skip the others'.
    if [[ " $owners " =~ \ engine\.(core|contributions|properties|equilibrium)\  ]] \
        || grep -Eq '^engine/tests/native/((equilibrium|boundary)_|.*\.hpp$)' <<<"$(grep -Ev '/(properties|regression)_checks_shared\.hpp$' <<<"$paths")"; then
        slow_tests+=('equilibrium/Automatic*' boundary)
    fi
    # The 16 binary128 reactive-caloric and reactive-cp cases and the reacting Cp fit (Explicit)
    # run only when the code they verify changes, or on request:
    # run-native-tests.sh 'equilibrium/Explicit/*' 'regression/Explicit/*'.
    if grep -Eq '^engine/(native/equilibrium/(implicit_sensitivity\.cpp|request_compile\.cpp|problem\.hpp)|native/properties/(thermochemistry\.|infinite_dilution_reference\.cpp)|native/regression/observations\.(cpp|hpp)|tests/native/((equilibrium_)?reactive_|quad_base\.hpp|binary128_support\.hpp))' <<<"$paths"; then
        slow_tests+=('equilibrium/Explicit/reactive-*' 'regression/Explicit/reactive-cp')
    fi
    # The envelope traces (about 100 s, Explicit) run when the tracer or its checks change.
    if grep -Eq '^engine/(native/equilibrium/(phase_boundary|implicit_sensitivity)\.cpp|tests/native/(equilibrium_critical_envelope_checks\.cpp|caloric_specification_shared\.hpp))' <<<"$paths"; then
        slow_tests+=('equilibrium/Explicit/envelope-*')
    fi
    if [[ " $owners " =~ \ engine\.(core|contributions|properties|equilibrium|regression)\  ]] \
        || grep -Eq '^engine/tests/native/(regression_|.*\.hpp$)' <<<"$(grep -Ev '/(equilibrium|properties)_checks_shared\.hpp$' <<<"$paths")"; then
        slow_tests+=(regression_fit)
    fi
    if grep -Eq '^engine/(native/contributions/association\.hpp|tests/native/association_comparison_checks\.cpp)' <<<"$paths"; then
        slow_tests+=(association-comparison)
    fi
fi

if grep -Eq '^(pyproject\.toml|uv\.lock|\.codex/|tools/(check-all|setup)\.sh)' <<<"$paths"; then
    timed check-all bash tools/check-all.sh
    check_declarations
    check_tools
    if ((${#slow_tests[@]})); then
        timed native-slow bash engine/tools/run-native-tests.sh "${slow_tests[@]}"
    fi
    echo "affected checks passed (base $base)"
    exit 0
fi

if [[ "$declarations_affected" == true ]] || grep -Eq '^(ARCHITECTURE\.yaml|engine/(meson\.build|build/subprojects/|LICENSE|pyproject\.toml|native/|src/|data/|tests/)|docs/scientific/|pyproject\.toml$|uv\.lock$|\.codex/|tools/(check-engine-architecture\.py|build-engine\.sh))' <<<"$paths"; then
    timed setup bash tools/setup.sh
fi
timed budgets python3 tools/diffstat.py --budgets

check_declarations
if grep -Eq '^(ARCHITECTURE\.yaml|engine/(meson\.build|build/subprojects/|native/|src/|tests/native/)|tools/(check-engine-architecture\.py|build-engine\.sh))' <<<"$paths"; then
    timed engine-architecture .venv/bin/python tools/check-engine-architecture.py --build-dir "$(<build/environment-wheel/meson-build-dir)"
fi

if grep -q '^engine/' <<<"$paths"; then
    test_paths=()
    grep -Eq '^engine/(native/|src/|tests/|data/|build/subprojects/|meson\.build|pyproject\.toml)' <<<"$paths" && test_paths=(tests)
    if ((${#test_paths[@]})); then
        mkdir -p build/test-results
        xml=build/test-results/pytest.xml s=0
        timed pytest env -C engine ../.venv/bin/python -m pytest -q "${test_paths[@]}" --junitxml="../$xml.tmp" || s=$?
        mv "$xml.tmp" "$xml" && python3 tools/report.py --stamp "$xml"
        ((s == 0)) || exit "$s"
    fi
    if grep -Eq '^engine/(native/|tests/native/|build/subprojects/|meson\.build)' <<<"$paths"; then
        # run-native-tests.sh compiles the programs it runs; slow cases run two at a time by default.
        timed native-default bash engine/tools/run-native-tests.sh
        if ((${#slow_tests[@]})); then
            timed native-slow bash engine/tools/run-native-tests.sh "${slow_tests[@]}"
        fi
    fi
    if grep -Eq '^engine/(src/|tests/)' <<<"$paths"; then
        timed ruff env -C engine ../.venv/bin/ruff check --select E4,E7,E9,F src tests
    fi
    if grep -Eq '^engine/(src/|tests/.*\.py$|pyproject\.toml$)' <<<"$paths"; then
        timed mypy env -C engine ../.venv/bin/mypy
    fi
fi

check_tools

grep -q '^analyses/' <<<"$paths" && bash -n tools/run-slice.sh

if grep -q '^data/' <<<"$paths"; then
    timed packets .venv/bin/python data/tools/validate_packets.py
fi

echo "affected checks passed (base $base)"
