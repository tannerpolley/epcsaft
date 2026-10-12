#!/usr/bin/env bash
set -euo pipefail

engine_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
workspace_root="$(dirname "$engine_root")"
cd "$engine_root"
export PATH="$workspace_root/.venv/bin:$PATH"

# --changed [base] checks only the native .cpp files changed since base (default: the merge base with
# origin/main), for review rounds; with no argument every native translation unit is checked.
files=('engine/(native|tests/native)/.*\.cpp$')
if [[ "${1:-}" == --changed ]]; then
    base="${2:-$(git merge-base origin/main HEAD)}"
    mapfile -t files < <({
        git diff --name-only --diff-filter=d "$base" -- 'native/*.cpp' 'tests/native/*.cpp'
        git ls-files --full-name --others --exclude-standard -- 'native/*.cpp' 'tests/native/*.cpp'
    } | sort -u | sed 's/\./\\./g; s/^/\//; s/$/$/')
    ((${#files[@]})) || { echo "no native .cpp file changed since $base"; exit 0; }
fi

for required_command in meson cmake ninja clang-18 clang++-18 clang-tidy-18 run-clang-tidy-18; do
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
# The checkout's build-job count sets the compile and clang-tidy parallelism.
jobs="$(cat "$workspace_root/build/environment-wheel/build-jobs" 2>/dev/null || echo 2)"
extension_suffix="$($workspace_root/.venv/bin/python -c 'import sysconfig; print(sysconfig.get_config_var("EXT_SUFFIX"))')"
# The tree is always the Clang debug build, so a warm tree reconfigures instead of rebuilding.
setup_args=()
if [[ -f "$build_root/meson-private/coredata.dat" ]]; then
    setup_args+=(--reconfigure)
elif [[ -d "$build_root" ]]; then
    rm -rf -- "$build_root"
fi
# ccache as tools/locked-ninja sets it up; debug info names the build directory, so a new worktree
# reuses another's objects only when ccache leaves the directory out of its hash.
export CCACHE_BASEDIR="$workspace_root" CCACHE_NOHASHDIR=1
CC="ccache clang-18" CXX="ccache clang++-18" meson setup \
    "${setup_args[@]}" \
    --buildtype=debug \
    "$build_root" \
    "$engine_root"
meson compile -C "$build_root" -j "$jobs" \
    "_core$extension_suffix"

# Clang ships no quadmath.h; the binary128 test units take GCC's after Clang's own headers.
run-clang-tidy-18 \
    -clang-tidy-binary clang-tidy-18 \
    -p "$build_root" \
    -extra-arg="-idirafter$(gcc -print-file-name=include)" \
    -j "$jobs" \
    -quiet \
    "${files[@]}"
