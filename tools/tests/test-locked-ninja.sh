#!/usr/bin/env bash
set -euo pipefail
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd -P)"
lock="$root/build/engine-meson/ninja.lock"
mkdir -p "$(dirname "$lock")"
set +e
flock "$lock" env EPCSAFT_NINJA_LOCKED=/another-checkout timeout 1 "$root/tools/locked-ninja" --version
status=$?
set -e
[[ $status == 124 ]]
flock "$lock" env EPCSAFT_NINJA_LOCKED="$root" timeout 1 "$root/tools/locked-ninja" --version
