#!/usr/bin/env bash
# Runs each argument as one bash command, all at the same time, from the checkout root; each logs to
# <n>.log in a new directory under build/checks/. Prints each command's exit status and wall time as
# it ends and exits 1 if any failed. Commands that share a Meson build directory must be chained inside one argument.
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
(($#)) || { echo "usage: $0 '<command>' ['<command>' ...]" >&2; exit 2; }
mkdir -p build/checks
log_dir="$(mktemp -d build/checks/run-XXXXXX)"
commands=("$@") pids=()
for i in "${!commands[@]}"; do
    (
        start="$EPOCHREALTIME"
        bash -c "${commands[i]}" >"$log_dir/$i.log" 2>&1
        status=$?
        elapsed="$(awk -v s="$start" -v e="$EPOCHREALTIME" 'BEGIN { printf "%.1f", e - s }')"
        printf '[exit %d, %s s] %s (%s)\n' "$status" "$elapsed" "${commands[i]}" "$log_dir/$i.log"
        exit "$status"
    ) &
    pids+=("$!")
done
failed=0
for pid in "${pids[@]}"; do wait "$pid" || failed=1; done
exit "$failed"
