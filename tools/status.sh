#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"
git status --short --branch

if [[ ! -x .venv/bin/python ]]; then
    echo "environment: missing (run bash .codex/environments/setup.sh)"
    exit 0
fi
.venv/bin/python -c 'import epcsaft; print("epcsaft wheel: %s\nsha256: %s" % epcsaft.installed_wheel())'
