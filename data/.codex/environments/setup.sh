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

mapfile -t runtime_distributions < <(
    .venv/bin/python - <<'PY'
from importlib import metadata
import re

for distribution in metadata.distributions():
    name = distribution.metadata.get("Name", "")
    if re.sub(r"[-_.]+", "-", name).lower().startswith("epcsaft"):
        print(name)
PY
)
if ((${#runtime_distributions[@]})); then
    uv pip uninstall --python .venv/bin/python "${runtime_distributions[@]}"
fi

.venv/bin/python tools/validate_packets.py
.venv/bin/python -m unittest discover -s tests -v
