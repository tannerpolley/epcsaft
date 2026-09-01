#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

prepare_native_build="${EPCSAFT_PREPARE_NATIVE_BUILD:-0}"
if [[ "$prepare_native_build" != 0 && "$prepare_native_build" != 1 ]]; then
    echo "EPCSAFT_PREPARE_NATIVE_BUILD must be 0 or 1" >&2
    exit 2
fi

python="$repo_root/.venv/bin/python"
if [[ ! -x "$python" ]]; then
    echo "missing root environment: $python" >&2
    exit 1
fi
export PATH="$repo_root/.venv/bin:$PATH"
read -r -a cc_command <<<"${CC:-cc}"
read -r -a cxx_command <<<"${CXX:-c++}"
for required_command in uv meson cmake ninja patchelf pkg-config sha256sum readlink flock xargs \
    dpkg dpkg-query ldd mktemp \
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

ceres_probe_dir="$(mktemp -d)"
trap 'rm -rf -- "$ceres_probe_dir"' EXIT
printf '%s\n' \
    'cmake_minimum_required(VERSION 3.25)' \
    'project(epcsaft_ceres_identity LANGUAGES CXX)' \
    'find_package(Ceres 2.2 CONFIG REQUIRED)' \
    'if(Ceres_VERSION VERSION_LESS 2.2 OR NOT Ceres_VERSION VERSION_LESS 2.3)' \
    '  message(FATAL_ERROR "Ceres 2.2.x is required, found ${Ceres_VERSION}")' \
    'endif()' \
    'get_target_property(ceres_includes Ceres::ceres INTERFACE_INCLUDE_DIRECTORIES)' \
    'get_target_property(ceres_links Ceres::ceres INTERFACE_LINK_LIBRARIES)' \
    'get_target_property(ceres_features Ceres::ceres INTERFACE_COMPILE_FEATURES)' \
    'get_target_property(ceres_configurations Ceres::ceres IMPORTED_CONFIGURATIONS)' \
    'set(ceres_locations "")' \
    'foreach(configuration IN LISTS ceres_configurations)' \
    '  get_target_property(location Ceres::ceres IMPORTED_LOCATION_${configuration})' \
    '  list(APPEND ceres_locations "${location}")' \
    'endforeach()' \
    'file(WRITE "${CMAKE_BINARY_DIR}/ceres-identity.txt"' \
    '  "version=${Ceres_VERSION}\nconfig=${Ceres_DIR}\nincludes=${ceres_includes}\n"' \
    '  "links=${ceres_links}\nfeatures=${ceres_features}\n"' \
    '  "configurations=${ceres_configurations}\nlocations=${ceres_locations}\n")' \
    >"$ceres_probe_dir/CMakeLists.txt"
cmake -S "$ceres_probe_dir" -B "$ceres_probe_dir/build" >/dev/null
ceres_identity="$ceres_probe_dir/build/ceres-identity.txt"
[[ -f "$ceres_identity" ]] || { echo "Ceres identity probe did not produce metadata" >&2; exit 1; }
ceres_version="$(sed -n 's/^version=//p' "$ceres_identity")"
ceres_config_dir="$(sed -n 's/^config=//p' "$ceres_identity")"
ceres_include_dirs="$(sed -n 's/^includes=//p' "$ceres_identity")"
ceres_locations="$(sed -n 's/^locations=//p' "$ceres_identity")"
ceres_library=""
while IFS= read -r candidate; do
    if [[ -f "$candidate" ]]; then
        ceres_library="$(readlink -f "$candidate")"
        break
    fi
done < <(tr ';' '\n' <<<"$ceres_locations")
[[ -n "$ceres_library" ]] || { echo "Ceres::ceres has no resolved imported library" >&2; exit 1; }
ceres_include_root=""
while IFS= read -r candidate; do
    if [[ -f "$candidate/ceres/version.h" ]]; then
        ceres_include_root="$(readlink -f "$candidate")"
        break
    fi
done < <(tr ';' '\n' <<<"$ceres_include_dirs")
[[ -n "$ceres_include_root" ]] || { echo "Ceres::ceres has no resolved header tree" >&2; exit 1; }
[[ -d "$ceres_config_dir" ]] || { echo "Ceres package metadata is unavailable" >&2; exit 1; }
for package in libceres-dev libceres4t64; do
    dpkg-query -W -f='${Status}\n' "$package" 2>/dev/null | grep -qx 'install ok installed' || {
        echo "missing Ubuntu Ceres package: $package" >&2
        exit 1
    }
done
ceres_dev_version="$(dpkg-query -W -f='${Version}' libceres-dev)"
ceres_runtime_version="$(dpkg-query -W -f='${Version}' libceres4t64)"
ceres_architecture="$(dpkg --print-architecture)"
if [[ "$ceres_dev_version" != "$ceres_runtime_version" ]] \
    || [[ "$(dpkg-query -W -f='${Architecture}' libceres-dev)" != "$ceres_architecture" ]] \
    || [[ "$(dpkg-query -W -f='${Architecture}' libceres4t64)" != "$ceres_architecture" ]]; then
    echo "Ubuntu Ceres development/runtime package identity is inconsistent" >&2
    exit 1
fi
dpkg-query -S "$ceres_config_dir/CeresConfig.cmake" | grep -Eq '^libceres-dev(:[^:]+)?:' || {
    echo "Ceres package metadata is not owned by libceres-dev" >&2
    exit 1
}
dpkg-query -S "$ceres_include_root/ceres/version.h" | grep -Eq '^libceres-dev(:[^:]+)?:' || {
    echo "Ceres headers are not owned by libceres-dev" >&2
    exit 1
}
dpkg-query -S "$ceres_library" | grep -Eq '^libceres4t64(:[^:]+)?:' || {
    echo "Ceres imported library is not owned by libceres4t64" >&2
    exit 1
}
ceres_runtime_links="$(ldd "$ceres_library")"
if grep -q '=> not found' <<<"$ceres_runtime_links"; then
    echo "Ceres has unresolved runtime libraries" >&2
    exit 1
fi

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
    printf 'ceres-imported-target=%s\n' "$(sha256sum "$ceres_identity" | cut -d ' ' -f 1)"
    printf 'ceres-version=%s\n' "$ceres_version"
    printf 'ceres-packages=%s\n' "$({
        dpkg-query -W -f='${Package}=${Version}\n' libceres-dev libceres4t64
    } | sort | tr '\n' ' ')"
    fingerprint_tree ceres-cmake-metadata "$ceres_config_dir"
    fingerprint_tree ceres-headers "$ceres_include_root/ceres"
    printf 'ceres-library=%s:%s\n' "$(basename "$ceres_library")" \
        "$(sha256sum "$ceres_library" | cut -d ' ' -f 1)"
    while IFS= read -r library_path; do
        library_path="$(readlink -f "$library_path")"
        printf 'ceres-runtime-library=%s:%s\n' "$(basename "$library_path")" \
            "$(sha256sum "$library_path" | cut -d ' ' -f 1)"
    done < <(awk '/=> \/[^ ]+/{print $3}' <<<"$ceres_runtime_links" | sort -u)
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
        meson setup "$build_dir" engine
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
        if [[ "$prepare_native_build" == 1 ]]; then
            ensure_configured_build
        fi
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
elif [[ "$prepare_native_build" == 1 ]]; then
    ensure_configured_build
fi
printf '%s\n' "$build_fingerprint" >"$stamp"
sha256sum "$wheel" | cut -d ' ' -f 1 >"$wheel_hash_stamp"

case "$wheel_origin" in
    cache) echo "Engine wheel restored from shared cache for build $build_fingerprint" ;;
    rebuilt) echo "Engine wheel rebuilt and cached for build $build_fingerprint" ;;
    local) echo "Engine wheel installed from the current checkout wheel" ;;
esac
