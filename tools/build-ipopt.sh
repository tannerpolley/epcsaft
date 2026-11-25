#!/usr/bin/env bash
# Linux source builds with the sequential Debian/Ubuntu MUMPS libraries.
set -euo pipefail
jobs=$1
version=3.14.20
digest=43bddd6fa793b1694aa94d6129fe3e4a8d452d97a84ef5f2ff721c8047c75605
# Print the search path for both the release and instrumented Engine builds.
if pkg-config --exact-version="$version" ipopt; then
    printf '%s\n' "${PKG_CONFIG_PATH:-}"
    exit 0
fi
key="$({
    sha256sum < "$0"
    for compiler in "${CC:-cc}" "${CXX:-c++}" "${FC:-gfortran}"; do
        read -ra command <<<"$compiler"
        "${command[@]}" --version
        printf '%s\n' "$compiler"
    done
    printf '%s\n' "${CPPFLAGS-}" "${CFLAGS-}" "${CXXFLAGS-}" "${FCFLAGS-}" "${LDFLAGS-}"
} | sha256sum | cut -c 1-16)"
prefix="${XDG_CACHE_HOME:-$HOME/.cache}/epcsaft/ipopt-$version-$key"
mkdir -p "$(dirname "$prefix")"
exec {lock}>"$prefix.lock"
flock "$lock"
printf '%s\n' "$prefix/lib/pkgconfig${PKG_CONFIG_PATH:+:$PKG_CONFIG_PATH}"
[[ -f "$prefix/lib/pkgconfig/ipopt.pc" && -f "$prefix/lib/libipopt.so.$version" ]] && exit 0
exec >&2
work=$(mktemp -d)
trap 'rm -rf -- "$work"' EXIT
cd "$work"
curl --fail --location --retry 3 --output ipopt.tar.gz "https://codeload.github.com/coin-or/Ipopt/tar.gz/refs/tags/releases/$version"
echo "$digest  ipopt.tar.gz" | sha256sum --check
tar -xzf ipopt.tar.gz --strip-components=1
./configure --prefix="$prefix" --enable-shared --disable-static --disable-java --disable-sipopt \
    --disable-mpiinit --disable-pardisomkl --without-hsl --without-spral --without-asl \
    --with-precision=double --with-intsize=32 --with-mumps-cflags='-I/usr/include -I/usr/include/mumps_seq' \
    --with-mumps-lflags=-ldmumps_seq --with-lapack-lflags='-llapack -lblas'
make -j"$jobs"
make install
