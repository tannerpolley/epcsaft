#!/usr/bin/env bash
# Fixed-settler candidates and their qualification, kept under the caller's output directory.
# Needs the current non-editable wheel in .venv: bash tools/build-engine.sh --install-wheel
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "usage: $0 OUTPUT_DIR" >&2
    exit 2
fi
out="$(realpath -m "$1")"
cd "$(dirname "${BASH_SOURCE[0]}")/../../.."
.venv/bin/python data/tools/validate_packets.py
.venv/bin/python analyses/2026-khudaida-idaes-settler/scripts/run_fixed_settler.py --output-dir "$out/candidates"
.venv/bin/python analyses/2026-khudaida-idaes-settler/scripts/qualify.py \
    --results-dir "$out/candidates" --output-dir "$out/qualification"
