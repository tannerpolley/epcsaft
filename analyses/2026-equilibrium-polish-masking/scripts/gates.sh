#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../../.."
bash tools/check-affected.sh origin/main
bash tools/check-all.sh
