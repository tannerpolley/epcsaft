#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../../../.."
leaf=analyses/2026-errors-in-variables/scripts/diagnosis
probe_build="$PWD/build/ionic-reference-allowance"
mkdir -p "$probe_build"
patch --silent -o "$probe_build/binary128_support_probe.hpp" \
    engine/tests/native/binary128_support.hpp "$leaf/binary128-reference.patch"
nice -n 10 c++ -DHAVE_CSTDDEF -std=c++17 -O2 -I/usr/include/coin \
    -Iengine/native -Iengine/tests/native -I/usr/include/eigen3 -I"$probe_build" \
    -Iengine/build/subprojects/CppAD-20250000.3/include "$leaf/ladder.cpp" \
    -Lbuild/engine-meson/shared -Wl,-rpath,"$PWD/build/engine-meson/shared" \
    -lepcsaft_native -lquadmath -o "$probe_build/ladder"
exec nice -n 10 "$probe_build/ladder" "$probe_build"
