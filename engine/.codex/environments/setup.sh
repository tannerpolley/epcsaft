#!/usr/bin/env bash
set -euo pipefail

repo_root="$(git rev-parse --show-toplevel)"
cd "$repo_root"

for required_command in uv cmake ninja c++; do
    if ! command -v "$required_command" >/dev/null 2>&1; then
        echo "missing required command: $required_command" >&2
        exit 1
    fi
done

if [[ ! -x .venv/bin/python ]]; then
    uv venv --python 3.13 .venv
fi

mapfile -t stale_distributions < <(
    .venv/bin/python - <<'PY'
from importlib import metadata
import re

for distribution in metadata.distributions():
    name = distribution.metadata.get("Name", "")
    normalized = re.sub(r"[-_.]+", "-", name).lower()
    if normalized.startswith("epcsaft") and normalized != "epcsaft":
        print(name)
PY
)
if ((${#stale_distributions[@]})); then
    uv pip uninstall --python .venv/bin/python "${stale_distributions[@]}"
fi

uv pip install --python .venv/bin/python \
    --reinstall-package epcsaft \
    --editable . \
    build \
    mypy \
    numpy \
    pytest \
    ruff

.venv/bin/python - <<'PY'
from importlib import metadata
from pathlib import Path

import epcsaft
from epcsaft import equilibrium, regression

repo_root = Path.cwd().resolve()
for module in (epcsaft, equilibrium, regression):
    module_path = Path(module.__file__).resolve()
    if repo_root not in module_path.parents:
        raise SystemExit(f"{module.__name__} is not installed from this worktree: {module_path}")
if metadata.version("epcsaft") != "0.2.0.dev0":
    raise SystemExit(f"unexpected epcsaft version: {metadata.version('epcsaft')}")

unexpected = []
for distribution in metadata.distributions():
    name = distribution.metadata.get("Name", "")
    normalized = name.lower().replace("_", "-").replace(".", "-")
    if normalized.startswith("epcsaft") and normalized != "epcsaft":
        unexpected.append(name)
if unexpected:
    raise SystemExit(f"unexpected additional ePC-SAFT distributions: {sorted(unexpected)}")
PY
