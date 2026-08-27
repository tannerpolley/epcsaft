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
run_check() {
    local name="$1"
    shift
    names+=("$name")
    ("$@") >"$run_root/$name.log" 2>&1 &
    pids+=("$!")
}

run_check packets .venv/bin/python data/tools/validate_packets.py
run_check engine-eos bash -c 'cd "$1/engine" && ../.venv/bin/python -m pytest -q -m "not slow and not expensive" tests/eos' _ "$repo_root"
run_check engine-equilibrium-compiler bash -c 'cd "$1/engine" && ../.venv/bin/python -m pytest -q -m "not slow and not expensive" tests/equilibrium/test_chemical_compilation.py' _ "$repo_root"
run_check engine-equilibrium-other bash -c 'cd "$1/engine" && ../.venv/bin/python -m pytest -q -m "not slow and not expensive" tests/equilibrium --ignore=tests/equilibrium/test_chemical_compilation.py' _ "$repo_root"
run_check engine-regression bash -c 'cd "$1/engine" && ../.venv/bin/python -m pytest -q -m "not slow and not expensive" tests/regression' _ "$repo_root"
run_check ruff bash -c 'cd "$1/engine" && ../.venv/bin/ruff check --select E4,E7,E9,F src tests' _ "$repo_root"
run_check mypy bash -c 'cd "$1/engine" && ../.venv/bin/mypy' _ "$repo_root"
run_check validation-fast bash -c 'cd "$1/validation" && ../.venv/bin/python -m pytest -q -m "not campaign and not slow and not expensive"' _ "$repo_root"
run_check management python3 governance/tools/management_preflight.py check
run_check artifacts python3 governance/tools/artifact_store.py check
run_check architecture .venv/bin/python tools/check-engine-architecture.py --build-dir "$build_dir" --wheel "$wheel"

failed=0
for index in "${!pids[@]}"; do
    name="${names[$index]}"
    if wait "${pids[$index]}"; then
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
echo "slow Engine evidence, Validation campaigns, clang-tidy, and sanitizers run only by explicit owner command"
