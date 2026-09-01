#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

# Build and install the immutable wheel before any test process starts.
EPCSAFT_PREPARE_NATIVE_BUILD=1 bash .codex/environments/setup.sh

build_dir="$(<build/environment-wheel/meson-build-dir)"
mapfile -t engine_wheels < <(find build/environment-wheel -maxdepth 1 -type f -name 'epcsaft-*.whl' -print)
if ((${#engine_wheels[@]} != 1)); then
    echo "expected exactly one Engine wheel, found ${#engine_wheels[@]}" >&2
    exit 1
fi
wheel="${engine_wheels[0]}"

run_root="$(mktemp -d)"
cleanup() {
    rm -rf -- "$run_root"
}
trap cleanup EXIT

names=()
pids=()
active=0
max_workers=4
declare -A exit_codes=()
run_check() {
    local name="$1"
    shift
    names+=("$name")
    (/usr/bin/time -f 'wall_seconds=%e' "$@") >"$run_root/$name.log" 2>&1 &
    pids+=("$!")
    ((active += 1))
    if ((active == max_workers)); then
        local finished_pid
        if wait -n -p finished_pid; then
            exit_codes["$finished_pid"]=0
        else
            exit_codes["$finished_pid"]=$?
        fi
        ((active -= 1))
    fi
}

run_check packets .venv/bin/python data/tools/validate_packets.py
run_check engine-eos bash -c 'cd "$1/engine" && ../.venv/bin/python -m pytest -q -m "not slow and not expensive" tests/eos' _ "$repo_root"
run_check engine-equilibrium bash -c 'cd "$1/engine" && ../.venv/bin/python -m pytest -q -m "not slow and not expensive" tests/equilibrium' _ "$repo_root"
run_check engine-regression bash -c 'cd "$1/engine" && ../.venv/bin/python -m pytest -q -m "not slow and not expensive" tests/regression' _ "$repo_root"
run_check ruff bash -c 'cd "$1/engine" && ../.venv/bin/ruff check --select E4,E7,E9,F src tests' _ "$repo_root"
run_check mypy bash -c 'cd "$1/engine" && ../.venv/bin/mypy' _ "$repo_root"
run_check analysis-layout .venv/bin/python tools/check-analysis-layout.py
run_check architecture .venv/bin/python tools/check-engine-architecture.py --build-dir "$build_dir" --wheel "$wheel"

failed=0
for index in "${!pids[@]}"; do
    name="${names[$index]}"
    pid="${pids[$index]}"
    status="${exit_codes[$pid]-}"
    if [[ -z "$status" ]]; then
        if wait "$pid"; then
            status=0
        else
            status=$?
        fi
    fi
    if ((status == 0)); then
        printf '\n[%s]\n' "$name"
        cat "$run_root/$name.log"
    else
        printf '\n[%s FAILED]\n' "$name" >&2
        cat "$run_root/$name.log" >&2
        failed=1
    fi
done

if ((failed)); then
    exit 1
fi

echo "fast installed-wheel checks passed"
echo "slow Engine evidence, analysis reproductions, clang-tidy, and sanitizers run only by explicit owner command"
