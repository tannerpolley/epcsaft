#!/usr/bin/env bash
set -euo pipefail

[[ "${1:-}" == --prepared && -d "${2:-}" ]] || {
    echo 'use tools/check-pr-native --pr N --head <full-sha> --base <full-sha> to prepare the shared tree' >&2; exit 2;
}
run_root="$(realpath "$2")"
results="${EPCSAFT_GATE_RESULTS:-$run_root}"
workspace_root="$run_root/source"
cd "$workspace_root/engine"
export PATH="$run_root/venv/bin:$PATH"
export EPCSAFT_SOURCE_ROOT="$workspace_root" EPCSAFT_BUILD_DIR="$run_root/meson"

if grep -Fq -- '-DNDEBUG' "$run_root/meson/compile_commands.json"; then
    echo "sanitizer compile commands unexpectedly define NDEBUG" >&2
    exit 1
fi
echo "sanitizer compile commands: NDEBUG absent"

uv pip check --python "$run_root/venv/bin/python"
for library in "$run_root/meson"/_core*.so "$run_root/meson/libepcsaft_native.so" "$run_root/meson/libepcsaft_asl.so"; do
    file "$library" | grep -Fq 'not stripped' || { echo "stripped or missing library: $library" >&2; exit 1; }
    nm --defined-only "$library" >"$results/$(basename "$library").symbols"
    [[ -s "$results/$(basename "$library").symbols" ]] || { echo "no symbols: $library" >&2; exit 1; }
    nm -D "$library" | grep '__asan_init' > /dev/null || { echo "uninstrumented library: $library" >&2; exit 1; }
    nm -D "$library" | grep '__ubsan_handle' > /dev/null || { echo "missing UBSan: $library" >&2; exit 1; }
done
echo "sanitizer extension, native library and ASL: instrumentation and symbols retained"

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

# Leak detection stays off in both steps. Address and undefined-behaviour errors halt.
LD_PRELOAD="$sanitizer_preload" \
ASAN_OPTIONS="detect_leaks=0:halt_on_error=1" \
UBSAN_OPTIONS="halt_on_error=1:print_stacktrace=1" \
PYTHONNOUSERSITE=1 \
    "$run_root/venv/bin/python" "$workspace_root/tools/engine_identity.py" >"$results/prepared-identity.json"

LD_PRELOAD="$sanitizer_preload" \
ASAN_OPTIONS="detect_leaks=0:halt_on_error=1" \
UBSAN_OPTIONS="halt_on_error=1:print_stacktrace=1" \
PYTHONNOUSERSITE=1 \
EPCSAFT_PYTEST_RESULTS="$results/pytest-results.json" \
    "$run_root/venv/bin/python" -m pytest -q --junitxml="$results/sanitizers.xml" \
        tests/eos/test_activity.py \
        tests/eos/test_diagnostics.py \
        tests/eos/test_state.py::test_density_anchor_outside_domain_keeps_native_convergence_failure \
        tests/equilibrium/test_problem.py
