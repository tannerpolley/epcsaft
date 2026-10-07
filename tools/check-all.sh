#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

/usr/bin/time -f 'readiness: %e s' bash -e -c '
python3 engine/scripts/docs/sync_equation_registry.py --check
python3 tools/diffstat.py --budgets
python3 tools/project-map.py
'

/usr/bin/time -f 'wheel-preparation: %e s' env EPCSAFT_PREPARE_NATIVE_BUILD=1 bash .codex/environments/setup.sh
build_dir="$(<build/environment-wheel/meson-build-dir)"
mapfile -t engine_wheels < <(find build/environment-wheel -maxdepth 1 -type f -name 'epcsaft-*.whl' -print)
if ((${#engine_wheels[@]} != 1)); then
    echo "expected exactly one Engine wheel, found ${#engine_wheels[@]}" >&2
    exit 1
fi
wheel="${engine_wheels[0]}"
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

# One native process overlaps one checking process; numerical libraries stay single-threaded.
(
    run_check native-tests env MESON_TESTTHREADS=1 EPCSAFT_NATIVE_BUILD_READY=1 engine/tools/run-native-tests.sh
    exit "$failed"
) &
native_pid=$!

run_check packets .venv/bin/python data/tools/validate_packets.py
run_check python-tests bash -c \
    'cd "$1/engine" && ../.venv/bin/python -m pytest -q --junitxml=../build/test-results/pytest.xml' \
    _ "$repo_root"
run_check ruff bash -c 'cd "$1/engine" && ../.venv/bin/ruff check --select E4,E7,E9,F src tests' _ "$repo_root"
run_check mypy bash -c 'cd "$1/engine" && ../.venv/bin/mypy' _ "$repo_root"
run_check architecture .venv/bin/python tools/check-engine-architecture.py --build-dir "$build_dir" --wheel "$wheel"

wait "$native_pid" || failed=1
if ((failed)); then
    exit 1
fi
echo "fast installed-wheel checks passed"
echo "slow Engine evidence, analysis reproductions, clang-tidy, and sanitizers run only by explicit owner command"
