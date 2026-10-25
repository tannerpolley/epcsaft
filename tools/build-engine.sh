#!/usr/bin/env bash
# Build the Engine from this checkout into build/engine-meson/shared, the Meson tree that
# imports, native tests, tools/native-run and wheels share.
#   (no argument)    install engine/ editable into .venv: the next `import epcsaft` rebuilds
#                    changed native sources, and Python edits need no rebuild
#   --wheel DIR      build one wheel into DIR (new or empty) and print its path: application
#                    pins (tools/refresh-main-wheel); comparisons use prepared environments
#   --install-wheel  build a wheel into build/environment-wheel/, install it non-editably and
#                    record wheel.sha256 and preparation.json, for an analysis that records a
#                    wheel; rerunning setup restores the editable install
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
python="$PWD/.venv/bin/python"
[[ -x "$python" ]] || { echo "missing root environment: run bash .codex/environments/setup.sh" >&2; exit 1; }
export PATH="$PWD/.venv/bin:$PATH" NINJA="$PWD/tools/locked-ninja"
# locked-ninja locks build/engine-meson/ninja.lock, so the directory must exist before any Ninja run.
mkdir -p build/engine-meson

# Every CPU, but no more compile jobs than memory holds: a cold build needs about 1 GiB per job
# and the binding files about 1.4 GiB (#319, #373). Memory is free RAM or, if tighter, the
# headroom under the nearest cgroup limit (agent sessions run under one). EPCSAFT_BUILD_JOBS overrides.
cpus="$(env -u OMP_NUM_THREADS nproc)"
memory_kib=$(awk '/^MemAvailable:/ {print $2}' /proc/meminfo)
group=/sys/fs/cgroup$(cut -d: -f3 /proc/self/cgroup)
while [[ $group == /sys/fs/cgroup/* ]]; do
    for limit_file in memory.high memory.max; do
        limit=$(cat "$group/$limit_file" 2>/dev/null) || continue
        [[ $limit == max ]] && continue
        headroom=$(( (limit - $(cat "$group/memory.current")) / 1024 ))
        ((headroom < memory_kib)) && memory_kib=$headroom
        break
    done
    group=$(dirname "$group")
done
memory_jobs=$(( memory_kib / (1536 * 1024) ))
build_jobs="${EPCSAFT_BUILD_JOBS:-$(( memory_jobs < 1 ? 1 : memory_jobs < cpus ? memory_jobs : cpus ))}"
[[ "$build_jobs" =~ ^[1-9][0-9]*$ ]] || { echo "EPCSAFT_BUILD_JOBS must be a positive integer" >&2; exit 2; }

# HiGHS compiles once per toolchain into a shared static build (#319). engine/meson.build
# falls back to the subproject without it; these options repeat the subproject's.
read -r -a cc_command <<<"${CC:-cc}"
read -r -a cxx_command <<<"${CXX:-c++}"
highs_wrap=engine/build/subprojects/highs.wrap
highs_source="$(sed -n 's/^directory = //p' "$highs_wrap")"
highs_options=(--buildtype=release -Db_ndebug=if-release -Ddefault_library=static -Ddo_install=false
    -Dfast_build=true -Dhighsint64=false -Duse_zlib=disabled -Dwith_pybind11=false -Dwith_tests=false)
highs_key="$({
    cat "$highs_wrap"
    printf '%s\n' "${highs_options[*]}" "${cc_command[*]}" "${cxx_command[*]}" \
        "${CPPFLAGS-}" "${CFLAGS-}" "${CXXFLAGS-}" "${LDFLAGS-}"
    meson --version
    for compiler in "${cc_command[0]}" "${cxx_command[0]}"; do
        sha256sum "$(readlink -f "$(command -v "$compiler")")"
    done
} | sha256sum | cut -c 1-16)"
highs_prefix="${XDG_CACHE_HOME:-$HOME/.cache}/epcsaft/${highs_source,,}-$highs_key"
highs_library="$highs_prefix/build/highs/libhighs.a"
mkdir -p "$(dirname "$highs_prefix")"
exec {highs_lock_fd}>"$highs_prefix.lock"
flock "$highs_lock_fd"
if [[ ! -f "$highs_prefix/pkgconfig/epcsaft-highs.pc" || ! -f "$highs_library" ]]; then
    rm -rf -- "$highs_prefix"
    meson subprojects download --sourcedir engine highs >&2
    mkdir -p "$highs_prefix/pkgconfig"
    cp -r "engine/build/subprojects/$highs_source" "$highs_prefix/src"
    meson setup "${highs_options[@]}" "$highs_prefix/build" "$highs_prefix/src" >&2
    meson compile -C "$highs_prefix/build" -j "$build_jobs" >&2
    s="$highs_prefix/src"
    printf '%s\n' 'Name: epcsaft-highs' 'Description: HiGHS for the Engine build' \
        "Version: ${highs_source#HiGHS-}" "Libs: $highs_library" \
        "Cflags: -I$s/extern -I$s/extern/catch -I$s/extern/pdqsort -I$(dirname "$highs_library") -I$s/highs" \
        >"$highs_prefix/pkgconfig/epcsaft-highs.pc"
fi
flock -u "$highs_lock_fd"

# The HiGHS key hashes the compiler binaries, flags and Meson version. A tree built under another
# key is wiped, so a compiler or HiGHS change cannot link objects Ninja would not rebuild.
build_dir="$PWD/build/engine-meson/shared"
[[ "$(cat "$build_dir/toolchain-key" 2>/dev/null)" == "$highs_key" ]] || rm -rf -- "$build_dir"
records=build/environment-wheel
mkdir -p "$records"
printf '%s\n' "$build_dir" >"$records/meson-build-dir"
printf '%s\n' "$build_jobs" >"$records/build-jobs"
config=(-C "build-dir=$build_dir" -C "setup-args=-Dpkg_config_path=$highs_prefix/pkgconfig"
    -C "compile-args=-j$build_jobs")

build_wheel() {
    # Commit-time member timestamps: the same commit and bytes give the same wheel SHA-256,
    # which tools/refresh-main-wheel reads as "unchanged" (while the bundled apt packages are unchanged).
    # auditwheel copies the Ipopt, MUMPS and Ceres closure into epcsaft.libs/, so the wheel runs on
    # glibc 2.39 or newer without them (#576); HSL loads at run time and never enters the wheel.
    local epoch="${SOURCE_DATE_EPOCH:-$(git log -1 --format=%ct)}"
    rm -rf -- "$1/unrepaired"
    SOURCE_DATE_EPOCH="$epoch" uv build --python "$python" --wheel --no-build-isolation --out-dir "$1/unrepaired" "${config[@]}" engine >&2
    SOURCE_DATE_EPOCH="$epoch" auditwheel repair --plat manylinux_2_39_x86_64 -w "$1" "$1"/unrepaired/epcsaft-*.whl >&2
    rm -rf -- "$1/unrepaired"
    find "$1" -maxdepth 1 -name 'epcsaft-*.whl'
}

case "${1:-}" in
"")
    # A wheel record would misidentify the editable build.
    rm -f "$records"/*.whl "$records/wheel.sha256" "$records/preparation.json"
    if "$python" - "$PWD/engine" <<'PY' 2>/dev/null && [[ -f "$build_dir/build.ninja" ]]; then
import json, sys
from importlib.metadata import distribution
url = json.loads(distribution("epcsaft").read_text("direct_url.json") or "{}")
sys.exit(not (url.get("dir_info", {}).get("editable") and url.get("url") == "file://" + sys.argv[1]))
PY
        meson compile -C "$build_dir" -j "$build_jobs"
    else
        uv pip install --python "$python" --no-deps --no-build-isolation --reinstall-package epcsaft \
            "${config[@]}" -C editable-verbose=true -e engine
    fi
    bash tools/prepare-serena-cpp.sh
    ;;
--wheel)
    (($# == 2)) && [[ ! -d "$2" || -z "$(ls -A "$2")" ]] || { echo "usage: $0 --wheel <new or empty dir>" >&2; exit 2; }
    build_wheel "$2"
    ;;
--install-wheel)
    rm -f "$records"/*.whl "$records/wheel.sha256" "$records/preparation.json"
    wheel="$(build_wheel "$records")"
    uv pip install --python "$python" --no-deps --reinstall-package epcsaft "$wheel"
    "$python" tools/verify-installed-wheel.py --wheel "$wheel"
    sha256sum "$wheel" | cut -d ' ' -f 1 >"$records/wheel.sha256"
    printf '{"build_dir": "%s", "wheel": "%s", "wheel_sha256": "%s"}\n' \
        build/engine-meson/shared "$wheel" "$(<"$records/wheel.sha256")" >"$records/preparation.json"
    ;;
*)
    echo "usage: $0 [--wheel <dir> | --install-wheel]" >&2
    exit 2
    ;;
esac
printf '%s\n' "$highs_key" >"$build_dir/toolchain-key"
printf '%s\n' "$python" >"$build_dir/python"
