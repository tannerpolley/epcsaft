#!/usr/bin/env bash
set -euo pipefail

engine_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
workspace_root="$(dirname "$engine_root")"
cd "$engine_root"

for required_command in cmake ninja clang++-18 clang-tidy-18 run-clang-tidy-18; do
    if ! command -v "$required_command" >/dev/null 2>&1; then
        echo "missing required command: $required_command" >&2
        exit 1
    fi
done
if [[ ! -x "$workspace_root/.venv/bin/python" ]]; then
    echo "missing root environment; run bash .codex/environments/setup.sh" >&2
    exit 1
fi

build_root="$workspace_root/build/native-quality"
pybind11_dir="$($workspace_root/.venv/bin/python -m pybind11 --cmakedir)"
cmake -S "$engine_root" -B "$build_root" -G Ninja \
    -DCMAKE_BUILD_TYPE=Debug \
    -DCMAKE_CXX_COMPILER=clang++-18 \
    -DPython_EXECUTABLE="$workspace_root/.venv/bin/python" \
    -Dpybind11_DIR="$pybind11_dir" \
    -DSKBUILD_PROJECT_NAME=epcsaft \
    -DSKBUILD_PROJECT_VERSION=0.2.0
cmake --build "$build_root" --target \
    epcsaft_core \
    epcsaft_engine_equilibrium_native \
    epcsaft_engine_regression_native

run-clang-tidy-18 \
    -clang-tidy-binary clang-tidy-18 \
    -p "$build_root" \
    -quiet \
    '(src/epcsaft/_native|src/epcsaft/regression/native|cpp/equilibrium)/.*\.cpp$'
