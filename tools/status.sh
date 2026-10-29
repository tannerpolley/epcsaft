#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"
if git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
    git status --short --branch
fi

if [[ ! -x .venv/bin/python ]]; then
    echo "environment: missing (run bash tools/setup.sh)"
    exit 0
fi
.venv/bin/python -c 'import epcsaft; print("epcsaft provenance: %s %s" % epcsaft.provenance())'
