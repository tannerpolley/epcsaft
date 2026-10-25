#!/usr/bin/env bash
set -euo pipefail
build="$PWD/build/issue-789/2026-khudaida"
mkdir -p "$build"
nb="$(.venv/bin/python -c 'print(__import__("nanobind").include_dir())')"
py="$(.venv/bin/python -c 'print(__import__("sysconfig").get_path("include"))')"
nice -n 10 c++ -shared -fPIC -std=c++17 -fvisibility=hidden \
    -Iengine/native -Iengine/build/subprojects/CppAD-20250000.3/include -Iengine/build/subprojects/CppAD-20250000.3 \
    -Ibuild/engine-meson/shared/build/subprojects/CppAD-20250000.3 -I/usr/include/eigen3 -I"$nb" -I"$nb/../ext/robin_map/include" -I"$py" \
    analyses/reproduction/2026-khudaida/scripts/figure2.cpp "$nb/../src/nb_combined.cpp" \
    -Lbuild/engine-meson/shared -Wl,-rpath,"$PWD/build/engine-meson/shared" -lepcsaft_native -o "$build/figure2.so"
PYTHONPATH="$build${PYTHONPATH:+:$PYTHONPATH}" nice -n 10 .venv/bin/python analyses/reproduction/2026-khudaida/scripts/reproduce.py
