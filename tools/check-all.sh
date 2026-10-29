#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"
# The report reads whatever records exist, pass or fail; it never changes this gate's status.
trap 'python3 tools/report.py || true' EXIT

# Setup first: the budget count formats Python with the locked .venv Ruff.
/usr/bin/time -f 'engine-build: %e s' bash tools/setup.sh
build_dir="$(<build/environment-wheel/meson-build-dir)"

/usr/bin/time -f 'readiness: %e s' bash -e -c '
python3 engine/scripts/docs/sync_equation_registry.py --check
python3 tools/diffstat.py --budgets
python3 tools/project-map.py
python3 tools/tests/test-time-key-solves-gate.py
'
mkdir -p build/test-results

failed=0
run_check() {
    local name="$1"
    shift
    printf '\n[%s]\n' "$name"
    if ! /usr/bin/time -f "$name: %e s" "$@"; then
        printf '[%s FAILED]\n' "$name" >&2
        failed=1
    fi
}

# Native tests run alongside the Python checks.
(
    run_check native-tests engine/tools/run-native-tests.sh
    exit "$failed"
) &
native_pid=$!

run_check packets .venv/bin/python data/tools/validate_packets.py
# JUnit goes to a temporary file first, so a reader never sees a partial report.
run_check python-tests bash -c 'cd "$1/engine" && s=0 && xml=../build/test-results/pytest.xml
    ../.venv/bin/python -m pytest -q --junitxml="$xml.tmp" || s=$?
    mv "$xml.tmp" "$xml" && python3 ../tools/report.py --stamp "$xml"; exit "$s"' _ "$repo_root"
run_check ruff bash -c 'cd "$1/engine" && ../.venv/bin/ruff check --select E4,E7,E9,F src tests' _ "$repo_root"
run_check mypy bash -c 'cd "$1/engine" && ../.venv/bin/mypy' _ "$repo_root"
run_check architecture .venv/bin/python tools/check-engine-architecture.py --build-dir "$build_dir"

wait "$native_pid" || failed=1
if ((failed)); then
    exit 1
fi
echo "fast checks passed"
echo "other gates and when they run: docs/agents/workflow.md \"Gates\""
