#!/usr/bin/env bash
set -euo pipefail

repo_root="$(git rev-parse --show-toplevel)"
cd "$repo_root"

for required_command in uv cmake ninja clang++-18; do
    if ! command -v "$required_command" >/dev/null 2>&1; then
        echo "missing required command: $required_command" >&2
        exit 1
    fi
done
if [[ ! -x .venv/bin/python ]]; then
    echo "missing Engine environment; run .codex/environments/setup.sh" >&2
    exit 1
fi

run_root="$repo_root/build/native-sanitizers"
wheelhouse="$run_root/wheelhouse"
cmake -E make_directory "$run_root"
cmake -E remove_directory "$wheelhouse"
cmake -E remove_directory "$run_root/venv"
cmake -E make_directory "$wheelhouse"

CMAKE_ARGS="-DEPCSAFT_ENABLE_SANITIZERS=ON -DCMAKE_CXX_COMPILER=clang++-18 -DCMAKE_INTERPROCEDURAL_OPTIMIZATION=OFF" \
    .venv/bin/python -m build --wheel --outdir "$wheelhouse" \
        --config-setting="build-dir=$run_root/build" \
        --config-setting="cmake.build-type=Debug"

mapfile -t wheels < <(find "$wheelhouse" -maxdepth 1 -type f -name 'epcsaft-*.whl')
if ((${#wheels[@]} != 1)); then
    echo "expected exactly one sanitizer wheel, found ${#wheels[@]}" >&2
    exit 1
fi

uv venv --python 3.13 "$run_root/venv"
uv pip install --python "$run_root/venv/bin/python" pytest "${wheels[0]}"

asan_runtime="$(clang++-18 -print-file-name=libclang_rt.asan-x86_64.so)"
if [[ ! -f "$asan_runtime" ]]; then
    echo "clang ASan runtime not found: $asan_runtime" >&2
    exit 1
fi

LD_PRELOAD="$asan_runtime" \
ASAN_OPTIONS="detect_leaks=0:halt_on_error=1" \
UBSAN_OPTIONS="halt_on_error=1:print_stacktrace=1" \
PYTHONNOUSERSITE=1 \
    "$run_root/venv/bin/python" -m pytest -q \
        tests/integration/test_scientific_slice.py \
        tests/equilibrium/test_saturation.py
