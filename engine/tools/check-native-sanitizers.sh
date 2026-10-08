#!/usr/bin/env bash
set -euo pipefail

engine_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
workspace_root="$(dirname "$engine_root")"
cd "$engine_root"
export PATH="$workspace_root/.venv/bin:$PATH"

for required_command in uv meson cmake ninja clang-18 clang++-18 file nm; do
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
# The sanitizer Meson tree ($run_root/build) stays between runs; the wheel and its venv do not.
rm -rf -- "$wheelhouse" "$run_root/venv"
mkdir -p "$wheelhouse"

CC=clang-18 CXX=clang++-18 \
    uv build --python "$workspace_root/.venv/bin/python" \
        --wheel --no-build-isolation --out-dir "$wheelhouse" \
        --config-setting="build-dir=$run_root/build" \
        --config-setting="setup-args=-Dbuildtype=debug" \
        --config-setting="setup-args=-Db_ndebug=false" \
        --config-setting="setup-args=-Db_lundef=false" \
        --config-setting="setup-args=-Db_sanitize=address,undefined" \
        --config-setting="setup-args=-Db_lto=false" \
        --config-setting="setup-args=-Dstrip=false"

if grep -Fq -- '-DNDEBUG' "$run_root/build/compile_commands.json"; then
    echo "sanitizer compile commands unexpectedly define NDEBUG" >&2
    exit 1
fi
echo "sanitizer compile commands: NDEBUG absent"

mapfile -t wheels < <(find "$wheelhouse" -maxdepth 1 -type f -name 'epcsaft-*.whl')
if ((${#wheels[@]} != 1)); then
    echo "expected exactly one sanitizer wheel, found ${#wheels[@]}" >&2
    exit 1
fi

uv venv --python 3.13 "$run_root/venv"
uv export --locked --all-packages --no-emit-workspace --no-hashes \
    --output-file "$run_root/constraints.txt" >/dev/null
uv pip install --python "$run_root/venv/bin/python" \
    --constraint "$run_root/constraints.txt" pytest "${wheels[0]}"
uv pip check --python "$run_root/venv/bin/python"
mapfile -t installed_extensions < <(
    find "$run_root/venv/lib" -type f \
        -name '_core*.so'
)
if ((${#installed_extensions[@]} != 1)); then
    echo "expected one installed sanitizer extension, found ${#installed_extensions[@]}" >&2
    exit 1
fi
for extension in "${installed_extensions[@]}"; do
    if ! file "$extension" | grep -Fq 'not stripped'; then
        echo "sanitizer extension lacks useful symbols: $extension" >&2
        exit 1
    fi
    nm --defined-only "$extension" >"$run_root/defined-symbols.txt"
    if [[ ! -s "$run_root/defined-symbols.txt" ]]; then
        echo "sanitizer extension has no defined symbols: $extension" >&2
        exit 1
    fi
done
echo "sanitizer extensions: installed with symbols retained"

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

# Leak detection stays off in both steps: CPython never unloads _core, so its pybind11 type and enum registrations
# (every record under pybind11_init__core) are reported as leaks at exit. Address and undefined-behaviour errors halt.
LD_PRELOAD="$sanitizer_preload" \
ASAN_OPTIONS="detect_leaks=0:halt_on_error=1" \
UBSAN_OPTIONS="halt_on_error=1:print_stacktrace=1" \
PYTHONNOUSERSITE=1 \
    "$run_root/venv/bin/python" "$workspace_root/tools/verify-installed-wheel.py" \
        --wheel "${wheels[0]}"

LD_PRELOAD="$sanitizer_preload" \
ASAN_OPTIONS="detect_leaks=0:halt_on_error=1" \
UBSAN_OPTIONS="halt_on_error=1:print_stacktrace=1" \
PYTHONNOUSERSITE=1 \
EPCSAFT_PYTEST_RESULTS="$run_root/pytest-results.json" \
    "$run_root/venv/bin/python" -m pytest -q \
        tests/eos/test_activity.py \
        tests/eos/test_diagnostics.py \
        tests/equilibrium/test_problem.py
