#!/usr/bin/env bash
set -euo pipefail
# Source the build to retain the private Ipopt pkg-config path for qualification.
source "$RECIPE_DIR/build.sh"
"$PYTHON" "$RECIPE_DIR/qualification/retain.py" "$jobs"
