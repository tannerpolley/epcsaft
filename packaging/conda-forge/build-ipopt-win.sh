#!/usr/bin/env bash
# Sourced by run_autotools_clang_conda_build; PREFIX is the private install prefix.
set -euo pipefail
: "${LIBRARY_PREFIX:?}" "${SRC_DIR:?}" "${PREFIX:?}" "${BUILD_PREFIX:?}"
[[ ${target_platform:?} == win-64 ]]
declare -F patch_libtool >/dev/null
export EPCSAFT_BUILD_JOBS=${EPCSAFT_BUILD_JOBS:-6}
source_dir=$(cygpath -u "$SRC_DIR")/ipopt-source
private_prefix=$(cygpath -m "$PREFIX")
library_prefix=$(cygpath -m "$LIBRARY_PREFIX")
evidence="$private_prefix/share/epcsaft-ipopt-probe"
mkdir -p "$evidence" "$source_dir/../ipopt-build"
for library in dmumps openblas; do
    test -f "$library_prefix/lib/$library.lib"
    printf '%s\n' "$library_prefix/lib/$library.lib"
done | tee "$evidence/link-inputs.txt"
FC="$(cygpath -u "$BUILD_PREFIX")/Library/bin/flang.exe"
export FC
test -x "$FC"
export F77="$FC"
# Autoconf 2.59 misreads Flang's MSVC verbose link output as GNU libraries.
# Supply real libraries, while leaving its mixed-language link probes enabled.
build_library_prefix=$(cygpath -m "$BUILD_PREFIX")/Library/lib
export FLIBS="-L$build_library_prefix -lflang -lflangrti -lompstub"
export FFLAGS="-O2 $FLIBS -lflangmain -lmsvcrt -Xlinker /subsystem:console"
export CFLAGS="$CFLAGS -DIPOPT_BUILD -fno-fast-math -ffp-contract=off"
export CXXFLAGS="$CXXFLAGS -DIPOPT_BUILD -fno-fast-math -ffp-contract=off -Wno-reserved-user-defined-literal"
# Preserve forward slashes in configure-generated include flags.
export CYGPATH_W="cygpath -m"
export CPPFLAGS="$CPPFLAGS -I$library_prefix/include"
export LDFLAGS="$LDFLAGS -L$library_prefix/lib"
{
    "$CC" --version
    "$CXX" --version
    "$F77" --version
} > "$evidence/compilers.txt" 2>&1
cd "$source_dir/../ipopt-build"
build_dir=$PWD
# Keep failed configure/link probes as well as successful configuration headers.
trap 'find "$build_dir" -name config.log -exec cp --parents {} "$evidence" \;' EXIT
# Native Windows make cannot resolve MSYS absolute /d source paths.
../ipopt-source/configure --prefix="$private_prefix" \
    --disable-static --enable-dependency-linking \
    --disable-linear-solver-loader --without-hsl --without-asl \
    --without-pardiso --without-wsmp \
    --with-blas=-lopenblas --with-lapack=-lopenblas \
    --with-mumps-incdir="$library_prefix/include/mumps_seq" \
    --with-mumps-lib=-ldmumps 2>&1 | tee "$evidence/configure.txt"
# The conda wrapper supplies COFF C++ exports. Rename the DLL before linking:
# renaming an installed DLL alone leaves its old identity in the import library.
while IFS= read -r libtool; do
    (
        cd "$(dirname "$libtool")"
        patch_libtool
        sed -i 's/^soname_spec=.*/soname_spec="epcsaft_ipopt_3119.dll"/; s/^library_names_spec=.*/library_names_spec="epcsaft_ipopt_3119.dll"/' libtool
    )
done < <(find "$PWD" -name libtool -type f)
nice -n 10 make -j"$EPCSAFT_BUILD_JOBS" 2>&1 | tee "$evidence/make.txt"
make install 2>&1 | tee "$evidence/install.txt"
mkdir -p "$private_prefix/lib/pkgconfig" "$private_prefix/bin"
dll=$(find "$private_prefix" -name epcsaft_ipopt_3119.dll -type f -print -quit)
test -n "$dll"
[[ $dll == "$private_prefix/bin/epcsaft_ipopt_3119.dll" ]] || cp "$dll" "$private_prefix/bin/"
implib=$(find "$PWD/Ipopt" -name '*ipopt.dll.lib' -type f -print -quit)
test -n "$implib"
cp "$implib" "$private_prefix/lib/epcsaft_ipopt_3119.lib"
cat > "$private_prefix/lib/pkgconfig/ipopt.pc" <<PC
prefix=$private_prefix
libdir=\${prefix}/lib
includedir=\${prefix}/include/coin
Name: Ipopt
Description: Private Ipopt 3.11.9 with MUMPS
Version: 3.11.9
Libs: -L\${libdir} -lepcsaft_ipopt_3119
Cflags: -I\${includedir}
PC
find . \( -name 'config*.h' -o -name libtool -o -name '*.def' \) -exec cp --parents {} "$evidence" \;
sha256sum "$private_prefix/bin/epcsaft_ipopt_3119.dll" \
    "$private_prefix/lib/epcsaft_ipopt_3119.lib" > "$evidence/sha256.txt"
llvm-readobj --coff-imports --coff-exports "$private_prefix/bin/epcsaft_ipopt_3119.dll" > "$evidence/dll.txt"
