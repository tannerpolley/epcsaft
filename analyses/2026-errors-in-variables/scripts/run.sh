#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../../.."
export EPCSAFT_BUILD_JOBS=6 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
leaf=analyses/2026-errors-in-variables
if [[ ! -x build/issue-510/campaign-run || "$leaf/scripts/campaign.cpp" -nt build/issue-510/campaign-run || engine/tests/native/binary128_support.hpp -nt build/issue-510/campaign-run ]]; then
    nice -n 10 c++ -DEPCSAFT_EIV -DHAVE_CSTDDEF -std=c++17 -O2 -I/usr/include/coin -Iengine/native -Iengine/tests/native -Ibuild/issue-510 -I/usr/include/eigen3 -Iengine/build/subprojects/CppAD-20250000.3/include "$leaf/scripts/campaign.cpp" -Lbuild/engine-meson/shared -Wl,-rpath,"$PWD/build/engine-meson/shared" -lepcsaft_native -lquadmath -o build/issue-510/campaign-run
fi
destination="$leaf"
[[ "$1" != prepare && "$1" != bootstrap ]] || destination="$leaf/inputs"
exec nice -n 10 build/issue-510/campaign-run "$1" "$destination"
