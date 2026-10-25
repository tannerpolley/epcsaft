#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../../.."
out=build/issue-272/resumed
.venv/bin/python analyses/2026-equilibrium-polish-masking/scripts/compare_polish.py
tools/time-key-solves "$out"/main-wheel/*.whl "$out"/candidate-wheel/*.whl --out analyses/2026-equilibrium-polish-masking/results/resumed/timing
bash analyses/2026-equilibrium-polish-masking/scripts/gates.sh
