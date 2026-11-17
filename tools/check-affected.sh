#!/usr/bin/env bash
# The native gate: native setup, readiness, the source rules and the default native suites on every change, plus the
# checks the changed paths select; then the remaining gates of docs/agents/workflow.md "Gates" those paths need.
#   check-affected.sh [base]   paths changed since base (default: the merge base with origin/main) plus uncommitted
#   check-affected.sh --all    no Git history and no path selection (source archives, the container check)
#   check-affected.sh --python [base|--all]   the lead's batch: editable extension, Python suite, linked-image rules (#58)
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"
# The report reads whatever records exist, pass or fail; it never changes this gate's status.
trap 'python3 tools/report.py || true' EXIT

# Runs one step and prints its elapsed time, keeping the step's exit status.
timed() {
    local name="$1" start="$EPOCHREALTIME" status=0
    shift
    "$@" || status=$?
    awk -v n="$name" -v s="$start" -v e="$EPOCHREALTIME" 'BEGIN { printf "[%s: %.1f s]\n", n, e - s }'
    return "$status"
}

python_batch=false
[[ "${1:-}" != --python ]] || { python_batch=true; shift; }
changed=() owners=""
if [[ "${1:-}" == --all ]]; then
    echo "check-affected: all fast checks, no path selection"
else
    base_ref="${1:-$(git merge-base origin/main HEAD 2>/dev/null || git merge-base main HEAD 2>/dev/null || true)}"
    if ! git rev-parse --verify "$base_ref" >/dev/null 2>&1; then
        [[ -z "${1:-}" ]] || { echo "check-affected: base '$1' does not resolve to a commit" >&2; exit 2; }
        base_ref="HEAD"
    fi
    echo "check-affected base: $base_ref ($(git rev-parse --short "$base_ref"))"
    mapfile -t changed < <(
        {
            git diff --name-only "$base_ref" HEAD
            git diff --name-only
            git diff --cached --name-only
            git ls-files --others --exclude-standard
        } | sort -u
    )
    # Layer owners from ARCHITECTURE.yaml.
    owners="$(printf '%s\n' "${changed[@]}" | PYTHONPATH=tools python3 -c '
import importlib, sys
m = importlib.import_module("project-map")
rules = m.declaration()["rules"]
print(*{m.resolve(p, rules)[0] for p in sys.stdin.read().splitlines() if p})')"
fi
paths="$(printf '%s\n' "${changed[@]}")"
# A tools/ change runs every tools/tests/* script once, in check_tools; the steps that also cover one skip it.
tools_changed=false
grep -q '^tools/' <<<"$paths" && tools_changed=true

notebooks=()
for path in "${changed[@]}"; do
    [[ "$path" == analyses/*/notebook.qmd || "$path" == analyses/*/README.md ]] && notebooks+=("$path")
done
((${#notebooks[@]} == 0)) || timed notebook-numbers python3 tools/check-notebook-numbers.py --base "$base_ref" "${notebooks[@]}"

# Setup first: the budget count formats Python with the locked .venv Ruff. Only the batch builds the extension.
setup_args=(--native)
[[ "$python_batch" == false ]] || setup_args=()
timed setup bash tools/setup.sh "${setup_args[@]}"
timed readiness bash -e -c '
python3 engine/tools/sync_equation_registry.py --check
python3 tools/diffstat.py --budgets
python3 tools/project-map.py
[[ "$1" == true ]] || python3 tools/tests/test-time-key-solves-gate.py
.venv/bin/python data/tools/validate_packets.py' readiness "$tools_changed"
if [[ -d analyses ]] && grep -Eq '^(analyses/|docs/|engine/docs/|engine/src/|engine/tests/.*\.py$|tools/|ARCHITECTURE\.yaml)' <<<"$paths"; then
    timed analysis-layout python3 tools/check-analysis-layout.py
fi

declarations_affected=false
analysis_checkers='^tools/(check-analysis-declarations\.py|analysis-declarations-excluded\.tsv|check-notebook-numbers\.py|tests/test-analysis-checks\.py)'
if grep -Eq "^(analyses/|data/|engine/|pyproject\.toml$|uv\.lock$|\.codex/)|$analysis_checkers" <<<"$paths"; then
    declarations_affected=true
fi
check_declarations() {
    # With tools changed, tools/tests/test-analysis-checks.py in check_tools runs the full check.
    if [[ -d analyses && "$declarations_affected" == true && "$tools_changed" == false ]]; then
        timed declarations .venv/bin/python tools/check-analysis-declarations.py
    fi
}

check_tools() {
    grep -Eq '^(engine/)?tools/' <<<"$paths" || return 0
    for script in $(grep -ls '^#!/usr/bin/env bash' tools/* engine/tools/*); do bash -n "$script"; done
    if grep -q '^tools/project-map\.py' <<<"$paths"; then python3 tools/project-map.py --self-test; fi
    [[ "$tools_changed" == true ]] || return 0
    for test in tools/tests/*; do
        [[ -d analyses || "$test" != tools/tests/test-analysis-checks.py ]] || continue
        if [[ "$test" == *.sh ]]; then timed "$test" bash "$test"; else timed "$test" .venv/bin/python "$test"; fi
    done
}

check_declarations
check_tools
while read -r script; do [[ ! -f "$script" ]] || bash -n "$script"; done < <(grep -E '^analyses/.*\.sh$' <<<"$paths")

slow_tests=()
if grep -Eq '^engine/(native/|tests/native/|build/subprojects/|meson\.build)' <<<"$paths"; then
    # The slow solved-state programs run when their layer or any layer below it changes:
    # a property or EOS change moves every solved state.
    # A program's *_checks_shared.hpp reaches only that program, so the header rules skip the others'.
    if [[ " $owners " =~ \ engine\.(core|contributions|properties|equilibrium)\  ]] \
        || grep -Eq '^engine/tests/native/((equilibrium|boundary)_|.*\.hpp$)' <<<"$(grep -Ev '/(properties|regression)_checks_shared\.hpp$' <<<"$paths")"; then
        slow_tests+=('equilibrium/Automatic*' boundary)
    fi
    # The Explicit cases run here only when the code they verify changes.
    if grep -Eq '^engine/(native/equilibrium/(implicit_sensitivity\.cpp|request_compile\.cpp|row_evaluation\.cpp|problem\.hpp)|native/properties/(thermochemistry\.|infinite_dilution_reference\.cpp)|native/regression/observations\.(cpp|hpp)|tests/native/((equilibrium_)?reactive_|quad_base\.hpp|binary128_support\.hpp))' <<<"$paths"; then
        slow_tests+=('equilibrium/Explicit/reactive-*' 'regression/Explicit/reactive-cp')
    fi
    if grep -Eq '^engine/(native/equilibrium/(phase_boundary|implicit_sensitivity)\.cpp|tests/native/(equilibrium_critical_envelope_checks\.cpp|caloric_specification_shared\.hpp))' <<<"$paths"; then
        slow_tests+=('equilibrium/Explicit/envelope-*')
    fi
    if [[ " $owners " =~ \ engine\.(core|contributions|properties|equilibrium|regression)\  ]] \
        || grep -Eq '^engine/tests/native/(regression_|.*\.hpp$)' <<<"$(grep -Ev '/(equilibrium|properties)_checks_shared\.hpp$' <<<"$paths")"; then
        slow_tests+=(regression_fit)
    fi
fi

# The native suites share one Meson build directory, so they run in sequence, alongside the Python checks.
(
    timed native-default bash engine/tools/run-native-tests.sh
    ((${#slow_tests[@]} == 0)) || timed native-slow bash engine/tools/run-native-tests.sh "${slow_tests[@]}"
) &
native_pid=$!

failed=0
mkdir -p build/test-results
# JUnit goes to a temporary file first, so a reader never sees a partial report.
xml=build/test-results/pytest.xml
if [[ "$python_batch" == true ]]; then
    timed pytest env -C engine ../.venv/bin/python -m pytest -q --junitxml="../$xml.tmp" || failed=1
    mv "$xml.tmp" "$xml" && python3 tools/report.py --stamp "$xml"
    # The committed native fixture is generated from these inputs through the extension this batch built (#62).
    if grep -Eq '^(tools/fixtures/mea-f1(-identity)?\.json(\.gz)?|engine/src/epcsaft/parameters\.py|engine/tests/native/fixtures/mea_capped_addition\.(py|inc))$' <<<"$paths" || [[ "${1:-}" == --all ]]; then
        timed mea-fixture-drift bash -c '.venv/bin/python engine/tests/native/fixtures/mea_capped_addition.py \
            "$(.venv/bin/python -c "import epcsaft._core as c; print(c.__file__)")" build/mea_capped_addition.inc &&
            cmp build/mea_capped_addition.inc engine/tests/native/fixtures/mea_capped_addition.inc' || failed=1
    fi
fi
timed ruff env -C engine ../.venv/bin/ruff check src tests || failed=1
timed mypy env -C engine ../.venv/bin/mypy || failed=1
architecture_args=()
[[ "$python_batch" == false ]] || architecture_args=(--build-dir "$(<build/environment-wheel/meson-build-dir)")
timed engine-architecture .venv/bin/python tools/check-engine-architecture.py "${architecture_args[@]}" || failed=1
wait "$native_pid" || failed=1
((failed == 0)) || exit 1

if [[ "${1:-}" == --all ]]; then
    echo "native gate passed; other gates and when they run: docs/agents/workflow.md \"Gates\""
    exit 0
fi
echo "native gate passed (base $base_ref)"
# One reader of docs/agents/workflow.md "Gates": the remaining rows these paths select.
gates=()
[[ " $owners " =~ \ engine\.(core|contributions|properties|equilibrium)\  ]] \
    && gates+=("lead's batch: tools/check-equilibrium-corpora <prepared-A> <prepared-B> --out build/<dir>")
[[ " $owners " =~ \ engine\.(equilibrium|regression)\  ]] \
    && gates+=("bash engine/tools/run-native-tests.sh 'equilibrium/Explicit/*' 'regression/Explicit/*'")
if grep -Eq '^engine/(native|src)/' <<<"$paths"; then
    gates+=("lead's batch: tools/time-key-solves <prepared-A> <prepared-B> --pairs 1 --out build/<dir>, and an A=A control (two workers from A)"
        "lead's batch, if fit scheduling, caching, optimizer work, observation composition or derivatives, reference computations, equilibrium solve work or density algorithms changed: add --fit --control <A=A-summary.json> (A=A run with --fit from the previous qualified batch environment)")
fi
grep -Eq '^engine/.*\.(cpp|hpp|h)$' <<<"$paths" \
    && gates+=("lead's batch: agent-heavy tools/check-pr-native --pr N --head <full-sha> --base <previous-qualified-sha>")
[[ "$python_batch" == true ]] || echo "lead's batch: bash tools/check-affected.sh --python <previous-qualified-sha> (Python suite, linked-image rules)"
echo "remaining gates (docs/agents/workflow.md \"Gates\"):"
((${#gates[@]})) && printf '  %s\n' "${gates[@]}" || echo "  none"
