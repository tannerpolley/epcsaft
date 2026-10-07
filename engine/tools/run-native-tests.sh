#!/usr/bin/env bash
set -euo pipefail

engine_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
workspace_root="$(dirname "$engine_root")"
cd "$workspace_root"

if [[ "${EPCSAFT_NATIVE_BUILD_READY:-0}" != 1 ]]; then
    EPCSAFT_PREPARE_NATIVE_BUILD=1 bash .codex/environments/setup.sh
fi
build_dir="$(<build/environment-wheel/meson-build-dir)"
tests=("$@")
if ((${#tests[@]} == 0)); then
    tests=(contributions properties diagnostics specification surface_example regression)
fi
# Test programs are not built by default (#319); compile only the ones the selected tests run.
mapfile -t programs < <(.venv/bin/meson introspect "$build_dir" --tests | .venv/bin/python -c '
import fnmatch, json, os, sys
tests = json.load(sys.stdin)
print(*sorted({os.path.basename(t["cmd"][0]) for t in tests if any(fnmatch.fnmatch(t["name"], p) for p in sys.argv[1:])}), sep="\n")' "${tests[@]}")
.venv/bin/meson compile -C "$build_dir" -j "$(<build/environment-wheel/build-jobs)" "${programs[@]}"
set +e
logbase="testlog-${BASHPID}-$(date +%s%N)"
.venv/bin/meson test -C "$build_dir" "${tests[@]}" --print-errorlogs --no-rebuild --logbase "$logbase"
meson_status=$?
set -e
mkdir -p build/test-results
report="$(mktemp build/test-results/native.xml.XXXXXX)"
meson_report="$build_dir/meson-logs/${logbase}.junit.xml"
if [[ -f "$meson_report" ]]; then
    cp "$meson_report" "$report"
else
    cat >"$report" <<EOF
<?xml version="1.0" encoding="UTF-8"?>
<testsuites tests="1" errors="1" failures="0" skipped="0"><testsuite name="native" tests="1" errors="1"><testcase name="meson invocation"><error message="Meson did not produce a JUnit report for this invocation"/></testcase></testsuite></testsuites>
EOF
fi
mv "$report" build/test-results/native.xml
rm -f "$build_dir/meson-logs/${logbase}."*
exit "$meson_status"
