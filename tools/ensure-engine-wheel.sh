#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

python="$repo_root/.venv/bin/python"
if [[ ! -x "$python" ]]; then
    echo "missing root environment: $python" >&2
    exit 1
fi

source_files=(
    engine/CMakeLists.txt
    engine/LICENSE
    engine/README.md
    engine/pyproject.toml
    uv.lock
)
while IFS= read -r -d '' path; do
    source_files+=("$path")
done < <(find engine/cpp engine/src -type f -print0 | sort -z)

source_sha256="$({
    for path in "${source_files[@]}"; do
        sha256sum "$path"
    done
} | sha256sum | cut -d ' ' -f 1)"

runtime_id="$("$python" - <<'PY'
import platform
import sys
import sysconfig

print("-".join((sys.implementation.cache_tag, sysconfig.get_platform(), platform.machine())))
PY
)"
toolchain_sha256="$({
    "$python" --version
    cmake --version | head -n 1
    c++ --version | head -n 1
} | sha256sum | cut -d ' ' -f 1)"

wheelhouse="$repo_root/build/environment-wheel"
stamp="$wheelhouse/source.sha256"
cache_root="${XDG_CACHE_HOME:-$HOME/.cache}/epcsaft/wheels"
cache_key="${source_sha256}-${runtime_id}-${toolchain_sha256}"
cache_entry="$cache_root/$cache_key"
cmake -E make_directory "$cache_root"
exec {wheelhouse_lock_fd}>"$cache_root/wheelhouse.lock"
flock "$wheelhouse_lock_fd"
mapfile -t wheels < <(find "$wheelhouse" -maxdepth 1 -type f -name 'epcsaft-*.whl' -print 2>/dev/null)
wheel_current=false
if ((${#wheels[@]} == 1)) && [[ -f "$stamp" ]] && [[ "$(<"$stamp")" == "$source_sha256" ]]; then
    wheel_current=true
fi

if [[ "$wheel_current" == true ]]; then
    wheel="${wheels[0]}"
    if verification_output="$("$python" tools/verify-installed-wheel.py --wheel "$wheel" 2>&1)"; then
        printf '%s\n' "$verification_output"
        echo "Engine wheel is current; native rebuild skipped"
        exit 0
    fi
fi

wheel_origin="local"
if [[ "$wheel_current" != true ]]; then
    exec {cache_lock_fd}>"$cache_root/$cache_key.lock"
    flock "$cache_lock_fd"

    cmake -E remove_directory "$wheelhouse"
    cmake -E make_directory "$wheelhouse"

    mapfile -t cached_wheels < <(
        find "$cache_entry" -maxdepth 1 -type f -name 'epcsaft-*.whl' -print 2>/dev/null
    )
    cache_valid=false
    if ((${#cached_wheels[@]} == 1)) && [[ -f "$cache_entry/wheel.sha256" ]]; then
        cached_sha256="$(sha256sum "${cached_wheels[0]}" | cut -d ' ' -f 1)"
        if [[ "$(<"$cache_entry/wheel.sha256")" == "$cached_sha256" ]]; then
            cache_valid=true
        fi
    fi

    if [[ "$cache_valid" == true ]]; then
        cmake -E copy "${cached_wheels[0]}" "$wheelhouse/"
        wheel_origin="cache"
    else
        "$python" -m build --wheel --no-isolation --outdir "$wheelhouse" engine
        mapfile -t built_wheels < <(
            find "$wheelhouse" -maxdepth 1 -type f -name 'epcsaft-*.whl' -print
        )
        if ((${#built_wheels[@]} != 1)); then
            echo "expected exactly one Engine wheel, found ${#built_wheels[@]}" >&2
            exit 1
        fi
        cmake -E remove_directory "$cache_entry"
        cmake -E make_directory "$cache_entry"
        cmake -E copy "${built_wheels[0]}" "$cache_entry/"
        sha256sum "${built_wheels[0]}" | cut -d ' ' -f 1 >"$cache_entry/wheel.sha256"
        wheel_origin="rebuilt"
    fi
    mapfile -t wheels < <(find "$wheelhouse" -maxdepth 1 -type f -name 'epcsaft-*.whl' -print)
    if ((${#wheels[@]} != 1)); then
        echo "expected exactly one Engine wheel, found ${#wheels[@]}" >&2
        exit 1
    fi
fi

wheel="${wheels[0]}"
uv pip install --python "$python" --reinstall-package epcsaft "$wheel"
"$python" tools/verify-installed-wheel.py --wheel "$wheel"
printf '%s\n' "$source_sha256" >"$stamp"

if [[ -n "${cache_lock_fd:-}" ]]; then
    flock -u "$cache_lock_fd"
fi
flock -u "$wheelhouse_lock_fd"

case "$wheel_origin" in
    cache) echo "Engine wheel restored from shared cache for source $source_sha256" ;;
    rebuilt) echo "Engine wheel rebuilt and cached for source $source_sha256" ;;
    local) echo "Engine wheel reinstalled after local verification failed" ;;
esac
