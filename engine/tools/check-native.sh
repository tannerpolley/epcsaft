#!/usr/bin/env bash
set -euo pipefail

engine_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
workspace_root="$(dirname "$engine_root")"
cd "$engine_root"
export PATH="$workspace_root/.venv/bin:$PATH"

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
extension_suffix="$($workspace_root/.venv/bin/python -c 'import sysconfig; print(sysconfig.get_config_var("EXT_SUFFIX"))')"
setup_args=()
if [[ -f "$build_root/meson-private/coredata.dat" ]]; then
    setup_args+=(--wipe)
elif [[ -d "$build_root" ]]; then
    rm -rf -- "$build_root"
fi
CC=clang-18 CXX=clang++-18 meson setup \
    "${setup_args[@]}" \
    --buildtype=debug \
    "$build_root" \
    "$engine_root"
meson compile -C "$build_root" \
    "_core$extension_suffix"

run-clang-tidy-18 \
    -clang-tidy-binary clang-tidy-18 \
    -p "$build_root" \
    -quiet \
    '(src/epcsaft/_native|src/epcsaft/regression/native|cpp/equilibrium)/.*\.cpp$'
