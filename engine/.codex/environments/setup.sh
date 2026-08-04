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

uv pip uninstall --python .venv/bin/python \
    epcsaft-equilibrium \
    epcsaft-regression

uv pip install --python .venv/bin/python \
    --reinstall-package epcsaft \
    --editable . \
    build \
    mypy \
    pytest \
    ruff

.venv/bin/python - <<'PY'
from importlib import metadata
from pathlib import Path

import epcsaft

repo_root = Path.cwd().resolve()
module_root = Path(epcsaft.__file__).resolve().parent
if repo_root not in module_root.parents:
    raise SystemExit(f"epcsaft is not installed from this worktree: {module_root}")
if metadata.version("epcsaft") != "0.2.0.dev0":
    raise SystemExit(f"unexpected epcsaft version: {metadata.version('epcsaft')}")

for forbidden in ("epcsaft-equilibrium", "epcsaft-regression"):
    try:
        metadata.distribution(forbidden)
    except metadata.PackageNotFoundError:
        continue
    raise SystemExit(f"forbidden sibling distribution is installed: {forbidden}")
PY
