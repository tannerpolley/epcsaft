#!/usr/bin/env bash
# Run "$@" on physical core 15/31 with every other user slice fenced off it; restore the slices on exit.
set -euo pipefail
slices=(app.slice background.slice session.slice agents.slice)
restore() { for s in "${slices[@]}"; do systemctl --user set-property --runtime "$s" AllowedCPUs=; done; }
trap restore EXIT
for s in "${slices[@]}"; do systemctl --user set-property --runtime "$s" AllowedCPUs=0-14,16-30; done
systemd-run --user --wait --pipe --collect --same-dir --unit="timing-reserved-$$" --slice=timing.slice -p AllowedCPUs=15,31 -- "$@"
