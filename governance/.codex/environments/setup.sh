#!/usr/bin/env bash
set -euo pipefail

repo_root="$(git rev-parse --show-toplevel)"
cd "$repo_root"

for required_command in git gh python3 uv; do
    if ! command -v "$required_command" >/dev/null 2>&1; then
        echo "missing required command: $required_command" >&2
        exit 1
    fi
done

if [[ ! -x .venv/bin/python ]]; then
    uv venv --python 3.13 .venv
fi

uv pip install --python .venv/bin/python pytest
.venv/bin/python tools/management_preflight.py check
