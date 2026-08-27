#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

python="$repo_root/.venv/bin/python"
if [[ ! -x "$python" ]]; then
    echo "missing root environment: $python" >&2
    exit 1
fi
export PATH="$repo_root/.venv/bin:$PATH"
read -r -a cc_command <<<"${CC:-cc}"
read -r -a cxx_command <<<"${CXX:-c++}"
for required_command in uv meson cmake ninja patchelf pkg-config sha256sum readlink flock xargs \
    "${cc_command[0]}" "${cxx_command[0]}"; do
    if ! command -v "$required_command" >/dev/null 2>&1; then
        echo "missing required command: $required_command" >&2
        exit 1
    fi
done

check_native_dependency() {
    local package="$1"
    local minimum="$2"
    if ! pkg-config --exists "$package"; then
        echo "missing native dependency: $package >= $minimum" >&2
        exit 1
    fi
    if ! pkg-config --atleast-version "$minimum" "$package"; then
        echo "unsupported native dependency: $package $(pkg-config --modversion "$package") (need >= $minimum)" >&2
        exit 1
    fi
}

check_native_dependency ipopt 3.11.9
check_native_dependency eigen3 3.4

fingerprint_tree() {
    local label="$1"
    local root="$2"
    if [[ ! -d "$root" ]] || [[ -z "$(find "$root" -type f -print -quit)" ]]; then
        echo "cannot fingerprint $label tree: $root" >&2
        return 1
    fi
    printf '%s=%s\n' "$label" "$(
        cd "$root"
        find . -type f -print0 | sort -z | xargs -0 sha256sum | sha256sum | cut -d ' ' -f 1
    )"
}

source_files=(
    tools/ensure-engine-wheel.sh
    engine/LICENSE
    engine/README.md
    engine/meson.build
    engine/pyproject.toml
)
while IFS= read -r -d '' path; do
    source_files+=("$path")
done < <(find engine/cpp engine/src -type f -print0 | sort -z)
while IFS= read -r -d '' path; do
    source_files+=("$path")
done < <(
    {
        find engine/build/subprojects -maxdepth 1 -type f -name '*.wrap' -print0
        find engine/build/subprojects/packagefiles -type f -print0
    } | sort -z
)

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
toolchain_fingerprint="$({
    "$python" - <<'PY'
import importlib.metadata
import platform
import sys
import sysconfig

for item in (
    sys.version,
    sys.implementation.cache_tag,
    sysconfig.get_config_var("SOABI"),
    sysconfig.get_config_var("EXT_SUFFIX"),
    sysconfig.get_platform(),
    platform.machine(),
    importlib.metadata.version("meson-python"),
    importlib.metadata.version("pybind11"),
):
    print(item)
PY
    uv --version
    meson --version
    cmake --version | sed -n '1p'
    ninja --version
    patchelf --version
    if command -v ccache >/dev/null 2>&1; then
        ccache_path="$(readlink -f "$(command -v ccache)")"
        printf 'ccache=%s\n' "$(sha256sum "$ccache_path" | cut -d ' ' -f 1)"
        ccache --version | sed -n '1p'
    fi
    for compiler in "${cc_command[0]}" "${cxx_command[0]}"; do
        compiler_path="$(readlink -f "$(command -v "$compiler")")"
        printf '%s=%s\n' "$compiler" "$(sha256sum "$compiler_path" | cut -d ' ' -f 1)"
    done
    printf 'CC-command=%q ' "${cc_command[@]}"
    printf '\nCXX-command=%q ' "${cxx_command[@]}"
    printf '\n'
    "${cc_command[@]}" --version | sed -n '1p'
    "${cc_command[@]}" -dumpmachine
    "${cc_command[@]}" -dumpfullversion -dumpversion
    "${cxx_command[@]}" --version | sed -n '1p'
    "${cxx_command[@]}" -dumpmachine
    "${cxx_command[@]}" -dumpfullversion -dumpversion
    linker_name="$("${cxx_command[@]}" -print-prog-name=ld)"
    linker_path="$(command -v "$linker_name")"
    linker_path="$(readlink -f "$linker_path")"
    [[ -x "$linker_path" ]] || { echo "cannot fingerprint linker: $linker_path" >&2; exit 1; }
    printf 'linker=%s\n' "$(sha256sum "$linker_path" | cut -d ' ' -f 1)"
    "$linker_path" --version | sed -n '1p'
    for variable in CC CXX AR LD CPPFLAGS CFLAGS CXXFLAGS LDFLAGS CCACHE_CONFIGPATH \
        CCACHE_DIR CCACHE_DISABLE CCACHE_NODIRECT CCACHE_PREFIX CCACHE_SLOPPINESS \
        PKG_CONFIG_PATH PKG_CONFIG_LIBDIR PKG_CONFIG_SYSROOT_DIR \
        CMAKE_PREFIX_PATH MESON_ARGS SOURCE_DATE_EPOCH; do
        printf '%s=%s\n' "$variable" "${!variable-}"
    done
    for package in ipopt eigen3; do
        pc_file="$(pkg-config --variable=pcfiledir "$package")/$package.pc"
        [[ -f "$pc_file" ]] || { echo "cannot fingerprint $package metadata: $pc_file" >&2; exit 1; }
        printf '%s-version=%s\n' "$package" "$(pkg-config --modversion "$package")"
        printf '%s-flags=%s\n' "$package" "$(pkg-config --cflags --libs "$package")"
        printf '%s-metadata=%s\n' "$package" "$(sha256sum "$pc_file" | cut -d ' ' -f 1)"
    done
    fingerprint_tree eigen-headers "$(pkg-config --variable=prefix eigen3)/include/eigen3/Eigen"
    fingerprint_tree ipopt-headers "$(pkg-config --variable=includedir ipopt)"
    for library in ipopt; do
        library_path="$("${cxx_command[@]}" -print-file-name="lib$library.so")"
        [[ "$library_path" != "lib$library.so" && -f "$library_path" ]] || {
            echo "cannot fingerprint native library: lib$library.so" >&2
            exit 1
        }
        library_path="$(readlink -f "$library_path")"
        printf '%s-library=%s:%s\n' \
            "$library" "$(basename "$library_path")" \
            "$(sha256sum "$library_path" | cut -d ' ' -f 1)"
    done
} | sha256sum | cut -d ' ' -f 1)"
build_fingerprint="$(
    printf 'source=%s\ntoolchain=%s\n' "$source_sha256" "$toolchain_fingerprint" \
        | sha256sum | cut -d ' ' -f 1
)"

wheelhouse="$repo_root/build/environment-wheel"
stamp="$wheelhouse/build.sha256"
wheel_hash_stamp="$wheelhouse/wheel.sha256"
cache_root="${XDG_CACHE_HOME:-$HOME/.cache}/epcsaft/wheels"
cache_key="$build_fingerprint"
cache_entry="$cache_root/$cache_key"
build_dir="$repo_root/build/engine-meson/$runtime_id-$toolchain_fingerprint"
build_dir_stamp="$wheelhouse/meson-build-dir"
configured_build_stamp="$build_dir/epcsaft-build.sha256"
mkdir -p "$repo_root/build" "$cache_root"
exec {wheelhouse_lock_fd}>"$repo_root/build/environment-wheel.lock"
flock "$wheelhouse_lock_fd"

ensure_configured_build() {
    if [[ ! -f "$build_dir/meson-info/meson-info.json" ]]; then
        meson setup --wrap-mode=nodownload "$build_dir" engine
    fi
    if [[ ! -f "$configured_build_stamp" ]] \
        || [[ "$(<"$configured_build_stamp")" != "$build_fingerprint" ]]; then
        meson compile -C "$build_dir"
        printf '%s\n' "$build_fingerprint" >"$configured_build_stamp"
    fi
    printf '%s\n' "$build_dir" >"$build_dir_stamp"
}

record_rebuilt_configured_build() {
    printf '%s\n' "$build_fingerprint" >"$configured_build_stamp"
    printf '%s\n' "$build_dir" >"$build_dir_stamp"
}

mapfile -t wheels < <(find "$wheelhouse" -maxdepth 1 -type f -name 'epcsaft-*.whl' -print 2>/dev/null)
wheel_current=false
if ((${#wheels[@]} == 1)) \
    && [[ -f "$stamp" && -f "$wheel_hash_stamp" ]] \
    && [[ "$(<"$stamp")" == "$build_fingerprint" ]] \
    && [[ "$(<"$wheel_hash_stamp")" == "$(sha256sum "${wheels[0]}" | cut -d ' ' -f 1)" ]]; then
    wheel_current=true
fi

if [[ "$wheel_current" == true ]]; then
    wheel="${wheels[0]}"
    if verification_output="$("$python" tools/verify-installed-wheel.py --wheel "$wheel" 2>&1)"; then
        ensure_configured_build
        uv pip check --python "$python"
        printf '%s\n' "$verification_output"
        echo "Engine wheel is current; native rebuild skipped"
        exit 0
    fi
fi

wheel_origin="local"
if [[ "$wheel_current" != true ]]; then
    exec {cache_lock_fd}>"$cache_root/$cache_key.lock"
    flock "$cache_lock_fd"

    rm -rf -- "$wheelhouse"
    mkdir -p "$wheelhouse"

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
        cp "${cached_wheels[0]}" "$wheelhouse/"
        wheel_origin="cache"
    else
        uv build \
            --python "$python" \
            --wheel \
            --no-build-isolation \
            --out-dir "$wheelhouse" \
            --config-setting="build-dir=$build_dir" \
            engine
        mapfile -t built_wheels < <(
            find "$wheelhouse" -maxdepth 1 -type f -name 'epcsaft-*.whl' -print
        )
        if ((${#built_wheels[@]} != 1)); then
            echo "expected exactly one Engine wheel, found ${#built_wheels[@]}" >&2
            exit 1
        fi
        rm -rf -- "$cache_entry"
        mkdir -p "$cache_entry"
        cp "${built_wheels[0]}" "$cache_entry/"
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
uv pip install --python "$python" --no-deps --reinstall-package epcsaft "$wheel"
uv pip check --python "$python"
"$python" tools/verify-installed-wheel.py --wheel "$wheel"
if [[ "$wheel_origin" == rebuilt ]]; then
    record_rebuilt_configured_build
else
    ensure_configured_build
fi
printf '%s\n' "$build_fingerprint" >"$stamp"
sha256sum "$wheel" | cut -d ' ' -f 1 >"$wheel_hash_stamp"

case "$wheel_origin" in
    cache) echo "Engine wheel restored from shared cache for build $build_fingerprint" ;;
    rebuilt) echo "Engine wheel rebuilt and cached for build $build_fingerprint" ;;
    local) echo "Engine wheel installed from the current checkout wheel" ;;
esac
