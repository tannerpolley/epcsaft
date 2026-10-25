#!/usr/bin/env bash
# Compile and run the independent row references with tools/native-driver from the repository
# root; each source's compile and run logs go to the output directory.
set -euo pipefail
out=${1:?usage: references.sh <output directory>}
src=$(dirname "${BASH_SOURCE[0]}")
mkdir -p "$out"
for job in properties_checks:surface-properties equilibrium_checks:surface-equilibrium regression_checks:prediction-reference \
           regression_checks:fit-row-references properties_checks:surface-state-references; do
 program=${job%%:*} name=${job#*:}
 printf 'independent reference %s\n' "$name"
 driver=$(tools/native-driver "$program" "$src/$name.cpp" 2> "$out/$name-compile.log")
 "$driver" > "$out/$name.log" 2>&1
done
