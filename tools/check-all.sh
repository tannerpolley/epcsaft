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

python3 tools/record-command.py --test-id tooling.check-all.equation-registry -- \
    python3 engine/scripts/docs/sync_equation_registry.py --check
python3 tools/record-command.py --test-id tooling.check-all.roadmap -- \
    python3 tools/check-roadmap.py
python3 tools/record-command.py --test-id tooling.check-all.budgets -- \
    python3 tools/diffstat.py --budgets
python3 tools/record-command.py --test-id tooling.check-all.project-map -- \
    python3 tools/project-map.py

# Build and install the immutable wheel before any test process starts.
EPCSAFT_PREPARE_NATIVE_BUILD=1 python3 tools/record-command.py \
    --test-id tooling.check-all.environment -- bash .codex/environments/setup.sh

build_dir="$(<build/environment-wheel/meson-build-dir)"
mapfile -t engine_wheels < <(find build/environment-wheel -maxdepth 1 -type f -name 'epcsaft-*.whl' -print)
if ((${#engine_wheels[@]} != 1)); then
    echo "expected exactly one Engine wheel, found ${#engine_wheels[@]}" >&2
    exit 1
fi
wheel="${engine_wheels[0]}"

failed=0
run_check() {
    local name="$1"
    shift
    printf '\n[%s]\n' "$name"
    if ! /usr/bin/time -f 'wall_seconds=%e' "$@"; then
        printf '[%s FAILED]\n' "$name" >&2
        failed=1
    fi
}

run_check packets python3 tools/record-command.py --test-id tooling.check-all.packets -- \
    .venv/bin/python data/tools/validate_packets.py
run_check recorded-evidence python3 tools/record-command.py \
    --test-id tooling.check-all.recorded-evidence -- \
    .venv/bin/python tools/run-numerical-evidence.py --automatic --no-build
run_check ruff python3 tools/record-command.py --test-id tooling.check-all.ruff -- \
    bash -c 'cd "$1/engine" && ../.venv/bin/ruff check --select E4,E7,E9,F src tests' _ "$repo_root"
run_check mypy python3 tools/record-command.py --test-id tooling.check-all.mypy -- \
    bash -c 'cd "$1/engine" && ../.venv/bin/mypy' _ "$repo_root"
run_check analysis-layout python3 tools/record-command.py \
    --test-id tooling.check-all.analysis-layout -- .venv/bin/python tools/check-analysis-layout.py
run_check architecture python3 tools/record-command.py --test-id tooling.check-all.architecture -- \
    .venv/bin/python tools/check-engine-architecture.py --build-dir "$build_dir" --wheel "$wheel"

if ((failed)); then
    exit 1
fi

echo "fast installed-wheel checks passed"
echo "slow Engine evidence, analysis reproductions, clang-tidy, and sanitizers run only by explicit owner command"
