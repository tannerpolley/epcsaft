#!/usr/bin/env bash
# Qualification build: the shipped build plus the qualification programs and development records.
set -euo pipefail
"$PYTHON" packaging/conda-forge/qualification/prepare.py
EPCSAFT_BUILD_DIR="$SRC_DIR/engine/build-unix" bash "$RECIPE_DIR/build.sh"
"$PYTHON" packaging/conda-forge/qualification/retain.py engine/build-unix
