#!/usr/bin/env bash
set -euo pipefail

engine_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
workspace_root="$(dirname "$engine_root")"
cd "$engine_root"

for required_command in uv cmake ninja clang++-18; do
    if ! command -v "$required_command" >/dev/null 2>&1; then
        echo "missing required command: $required_command" >&2
        exit 1
    fi
done
if [[ ! -x "$workspace_root/.venv/bin/python" ]]; then
    echo "missing root environment; run bash .codex/environments/setup.sh" >&2
    exit 1
fi

run_root="$workspace_root/build/native-sanitizers"
wheelhouse="$run_root/wheelhouse"
cleanup() {
    cmake -E remove_directory "$run_root"
}
trap cleanup EXIT
cmake -E make_directory "$run_root"
cmake -E remove_directory "$wheelhouse"
cmake -E remove_directory "$run_root/venv"
cmake -E make_directory "$wheelhouse"

CMAKE_ARGS="-DEPCSAFT_ENABLE_SANITIZERS=ON -DCMAKE_CXX_COMPILER=clang++-18 -DCMAKE_INTERPROCEDURAL_OPTIMIZATION=OFF" \
    "$workspace_root/.venv/bin/python" -m build --wheel --outdir "$wheelhouse" \
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
ubsan_runtime="$(clang++-18 -print-file-name=libclang_rt.ubsan_standalone-x86_64.so)"
if [[ ! -f "$ubsan_runtime" ]]; then
    echo "clang UBSan runtime not found: $ubsan_runtime" >&2
    exit 1
fi
sanitizer_preload="$asan_runtime:$ubsan_runtime"

LD_PRELOAD="$sanitizer_preload" \
UBSAN_OPTIONS="halt_on_error=1:print_stacktrace=1" \
PYTHONNOUSERSITE=1 \
    "$run_root/venv/bin/python" "$workspace_root/tools/verify-installed-wheel.py" \
        --wheel "${wheels[0]}"

LD_PRELOAD="$sanitizer_preload" \
ASAN_OPTIONS="detect_leaks=0:halt_on_error=1" \
UBSAN_OPTIONS="halt_on_error=1:print_stacktrace=1" \
PYTHONNOUSERSITE=1 \
    "$run_root/venv/bin/python" -m pytest -q \
        tests/eos/numerical/test_pure_saturation.py \
        tests/regression/numerical/test_generic_fit.py \
        tests/equilibrium/test_phase_equilibrium.py::test_declared_single_phase_solves_role_and_hints_do_not_change_root
