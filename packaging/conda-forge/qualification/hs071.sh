#!/usr/bin/env bash
# Sourced by run_autotools_clang_conda_build after the private Ipopt build.
set -eo pipefail
private_prefix=$(cygpath -m "$PREFIX")
cd "$(cygpath -u "$SRC_DIR")"
"$BUILD_PREFIX/python.exe" packaging/conda-forge/qualification/ipopt-prepare.py
example=ipopt-source/Ipopt/examples/hs071_cpp
for main in "$example/hs071_main.cpp" packaging/conda-forge/qualification/internal_caller.cpp; do
    name=hs071
    [[ $main == packaging/conda-forge/qualification/internal_caller.cpp ]] && name=internal_caller
    nice -n 10 "$CXX" $CXXFLAGS -DHAVE_CONFIG_H -I"$private_prefix/include/coin" -I"$example" \
        "$main" "$example/hs071_nlp.cpp" -L"$private_prefix/lib" \
        -lepcsaft_ipopt_3119 -lpsapi -o "$private_prefix/bin/$name.exe"
done
