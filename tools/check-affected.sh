#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

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

notebooks=()
for path in "${changed[@]}"; do
    [[ "$path" == analyses/*/notebook.qmd ]] && notebooks+=("$path")
done
((${#notebooks[@]} == 0)) || python3 tools/check-notebook-numbers.py "${notebooks[@]}"

paths="$(printf '%s\n' "${changed[@]}")"
python3 tools/diffstat.py --budgets
python3 tools/project-map.py
if grep -Eq '^(docs/|engine/docs/|AGENTS.md|README.md|ARCHITECTURE\.yaml)' <<<"$paths"; then
    python3 engine/scripts/docs/sync_equation_registry.py --check
fi
if grep -Eq '^(analyses/|docs/|engine/|tools/|ARCHITECTURE\.yaml)' <<<"$paths"; then
    python3 tools/check-analysis-layout.py
fi

declarations_affected=false
if grep -Eq '^(analyses/|data/|engine/|pyproject\.toml$|uv\.lock$|\.codex/|tools/(check-analysis-declarations\.py|analysis-declarations-excluded\.tsv|tests/test-analysis-checks\.py))' <<<"$paths"; then
    declarations_affected=true
fi
check_declarations() {
    if [[ "$source_only" == 0 && "$declarations_affected" == true ]]; then
        .venv/bin/python tools/check-analysis-declarations.py
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
    bash -n tools/check-affected.sh tools/check-all.sh
    echo "documentation readiness passed; no wheel build or numerical tests ran"
    exit 0
fi

slow_tests=()
if [[ "$source_only" == 0 ]] && grep -Eq '^engine/(native/|tests/native/|build/subprojects/|meson\.build)' <<<"$paths"; then
    # The slow solved-state programs run when their layer or any layer below it changes
    # (owners from ARCHITECTURE.yaml): a property or EOS change moves every solved state.
    owners="$(printf '%s\n' "${changed[@]}" | PYTHONPATH=tools python3 -c '
import importlib, sys
m = importlib.import_module("project-map")
rules = m.declaration()["rules"]
print(*{m.resolve(p, rules)[0] for p in sys.stdin.read().splitlines()})')"
    if [[ " $owners " =~ \ engine\.(core|contributions|properties|equilibrium)\  ]] || grep -Eq '^engine/tests/native/(equilibrium|boundary)_' <<<"$paths"; then
        slow_tests+=('equilibrium/Automatic*' boundary)
    fi
    # The 16 binary128 reactive-caloric and reactive-cp cases (Explicit) run only when
    # the code they verify changes, or on request: run-native-tests.sh 'equilibrium/Explicit/*'.
    if grep -Eq '^engine/(native/equilibrium/implicit_sensitivity\.cpp|native/properties/thermochemistry\.|tests/native/(reactive_|quad_base\.hpp))' <<<"$paths"; then
        slow_tests+=('equilibrium/Explicit/*')
    fi
    # The envelope traces (about 100 s, Explicit) run when the tracer or its checks change.
    if grep -Eq '^engine/(native/equilibrium/phase_boundary\.cpp|tests/native/critical_envelope_checks\.hpp)' <<<"$paths"; then
        slow_tests+=('equilibrium/Explicit/envelope-*')
    fi
    if [[ " $owners " =~ \ engine\.(core|contributions|properties|equilibrium|regression)\  ]] || grep -q '^engine/tests/native/regression_checks' <<<"$paths"; then
        slow_tests+=(regression_fit)
    fi
fi

if [[ "$source_only" == 0 ]] && grep -Eq '^(AGENTS.md|README.md|pyproject.toml|uv.lock|\.codex/|tools/check-all\.sh)' <<<"$paths"; then
    bash tools/check-all.sh
    check_declarations
    if ((${#slow_tests[@]})); then
        MESON_TESTTHREADS="${MESON_TESTTHREADS:-2}" EPCSAFT_NATIVE_BUILD_READY=1 bash engine/tools/run-native-tests.sh "${slow_tests[@]}"
    fi
    exit 0
fi

prepare_native_build=0
if [[ "$source_only" == 0 ]] && grep -Eq '^(ARCHITECTURE\.yaml|engine/(meson\.build|build/subprojects/|native/|src/|tests/native/)|tools/(check-engine-architecture\.py|verify-installed-wheel\.py))' <<<"$paths"; then
    prepare_native_build=1
fi
if [[ "$source_only" == 0 ]] && { [[ "$declarations_affected" == true ]] || grep -Eq '^(ARCHITECTURE\.yaml|engine/(meson\.build|build/subprojects/|LICENSE|README\.md|pyproject\.toml|native/|src/|data/|tests/)|docs/scientific/|pyproject\.toml$|uv\.lock$|\.codex/)' <<<"$paths"; }; then
    EPCSAFT_PREPARE_NATIVE_BUILD="$prepare_native_build" bash .codex/environments/setup.sh
fi

check_declarations
if [[ "$prepare_native_build" == 1 ]]; then
    build_dir="$(<build/environment-wheel/meson-build-dir)"
    mapfile -t engine_wheels < <(find build/environment-wheel -maxdepth 1 -type f -name 'epcsaft-*.whl' -print)
    if ((${#engine_wheels[@]} != 1)); then
        echo "expected exactly one Engine wheel, found ${#engine_wheels[@]}" >&2
        exit 1
    fi
    .venv/bin/python tools/check-engine-architecture.py --build-dir "$build_dir" --wheel "${engine_wheels[0]}"
fi

if grep -q '^engine/' <<<"$paths"; then
    test_paths=()
    grep -Eq '^engine/(native/|src/|tests/|data/|build/subprojects/|meson\.build|pyproject\.toml)' <<<"$paths" && test_paths=(tests)
    if [[ "$source_only" == 0 ]] && ((${#test_paths[@]})); then
        mkdir -p build/test-results
        (cd engine && ../.venv/bin/python -m pytest -q "${test_paths[@]}" --junitxml=../build/test-results/pytest.xml)
    fi
if [[ "$source_only" == 0 ]] && grep -Eq '^engine/(native/|tests/native/|build/subprojects/|meson\.build)' <<<"$paths"; then
        # run-native-tests.sh compiles the programs it runs; slow cases run two at a time by default.
        EPCSAFT_NATIVE_BUILD_READY=1 bash engine/tools/run-native-tests.sh
        if ((${#slow_tests[@]})); then
            MESON_TESTTHREADS="${MESON_TESTTHREADS:-2}" EPCSAFT_NATIVE_BUILD_READY=1 bash engine/tools/run-native-tests.sh "${slow_tests[@]}"
        fi
    fi
    grep -Eq '^engine/(src/|tests/)' <<<"$paths" && (cd engine && ../.venv/bin/ruff check --select E4,E7,E9,F src tests)
    grep -Eq '^engine/(src/|pyproject\.toml$)' <<<"$paths" && (cd engine && ../.venv/bin/mypy)
    grep -q '^engine/tools/' <<<"$paths" && bash -n engine/tools/*.sh
fi

if grep -q '^analyses/' <<<"$paths"; then
    bash -n tools/run-slice.sh
    if [[ "$source_only" == 1 ]]; then
        echo "source-only mode skipped installed-artifact analysis replay"
    fi
fi

if grep -q '^data/' <<<"$paths"; then
    .venv/bin/python data/tools/validate_packets.py
fi
