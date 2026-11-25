#!/usr/bin/env bash
# Linux source builds with the sequential Debian/Ubuntu MUMPS libraries.
set -euo pipefail
prefix=$1
jobs=$2
version=3.14.20
digest=43bddd6fa793b1694aa94d6129fe3e4a8d452d97a84ef5f2ff721c8047c75605
mkdir -p "$(dirname "$prefix")"
exec {lock}>"$prefix.lock"
flock "$lock"
[[ -f "$prefix/lib/pkgconfig/ipopt.pc" && -f "$prefix/lib/libipopt.so.$version" ]] && exit 0
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
