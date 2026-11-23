#!/usr/bin/env bash
set -euo pipefail
source "$RECIPE_DIR/build-env.sh"
python="$BUILD_PREFIX/bin/python"
if [[ $target_platform == osx-* ]]; then
    find "$SRC_DIR/ipopt-source" \( -name config.guess -o -name config.sub \) -exec sh -c 'cp "$BUILD_PREFIX/share/gnuconfig/$(basename "$1")" "$1"' _ {} \;
fi
mkdir -p "$BUILD_PREFIX/epcsaft-ipopt-build"
cd "$BUILD_PREFIX/epcsaft-ipopt-build"
CFLAGS="${CFLAGS:-} -std=gnu11" \
CXXFLAGS="${CXXFLAGS:-} -std=gnu++11 -DIPOPT_BUILD" \
CPPFLAGS="${CPPFLAGS:-} -I$PREFIX/include" \
LDFLAGS="${LDFLAGS:-} -L$PREFIX/lib -Wl,-rpath,$PREFIX/lib" \
"$SRC_DIR/ipopt-source/configure" \
    --prefix="$IPOPT_PREFIX" --disable-static --enable-dependency-linking \
    --disable-linear-solver-loader --without-hsl --without-asl \
    --without-pardiso --without-wsmp \
    --with-blas=-lblas --with-lapack=-llapack \
    --with-mumps-incdir="$PREFIX/include/mumps_seq" \
    --with-mumps-lib="$(pkg-config --libs dmumps_seq)"
jobs=$("$python" "$RECIPE_DIR/build-jobs.py" ipopt)
make -j"$jobs"
make install
# The native cache can be restored under a different work directory.
"$python" - "$IPOPT_PREFIX/lib/pkgconfig/ipopt.pc" <<'PYTHON'
from pathlib import Path
import re, sys
pc = Path(sys.argv[1])
pc.write_text(re.sub(r"(?m)^prefix=.*$", r"prefix=${pcfiledir}/../..", pc.read_text()))
PYTHON
if [[ $target_platform == osx-* ]]; then
    ipopt_library=$(basename "$(otool -D "$IPOPT_PREFIX/lib/libipopt.dylib" | tail -n 1)")
    install_name_tool -id "@rpath/$ipopt_library" "$IPOPT_PREFIX/lib/$ipopt_library"
fi
cd "$SRC_DIR"
cp -a engine/build/subprojects/packagefiles/cppad/. engine/build/subprojects/CppAD-20250000.3/
build_dir="$BUILD_PREFIX/epcsaft-native-build"
stage=engine/conda-native
if [[ ${EPCSAFT_QUALIFICATION:-} ]]; then
    "$python" "$RECIPE_DIR/qualification/prepare.py"
fi
meson setup "$build_dir" engine --wrap-mode=nodownload -Depcsaft_build=native
jobs=$("$python" "$RECIPE_DIR/build-jobs.py" native-ninja)
meson compile -C "$build_dir" -j"$jobs" -v
mkdir -p "$stage/cppad"
cp "$build_dir"/libepcsaft_*"$SHLIB_EXT" "$stage/"
cp "$build_dir/build/subprojects/CppAD-20250000.3/cppad/configure.hpp" "$stage/cppad/"
if [[ ${EPCSAFT_QUALIFICATION:-} ]]; then
    "$python" "$RECIPE_DIR/qualification/retain.py" --compile "$build_dir" "$jobs"
fi
