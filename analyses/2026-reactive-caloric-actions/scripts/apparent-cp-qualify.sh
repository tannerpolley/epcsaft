#!/usr/bin/env bash
set -euo pipefail
native_dir=${EPCSAFT_NATIVE_CHECKS:-$(<build/environment-wheel/meson-build-dir)}
raw=${EPCSAFT_CALORIC_RAW:-build/current-calorics/slice4}
for selector in reactive-cp-F1 reactive-cp-F3 reactive-caloric-F2; do
    "$native_dir/equilibrium_checks" "$selector" | tee "$raw/qualification-$selector.log"
done
for selector in reactive-cp addition-enthalpy reactive-caloric-F6-addition ionic-addition; do
    "$native_dir/regression_checks" "$selector" | tee "$raw/qualification-$selector.log"
done
"${EPCSAFT_PYTHON:-.venv/bin/python}" -m pytest engine/tests/equilibrium/test_apparent_heat_capacity.py -q | tee "$raw/qualification-installed-wheel.log"
