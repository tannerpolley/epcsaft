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

stale_build_root=""
while IFS= read -r -d '' cmake_cache; do
    cached_root="$(sed -n 's/^CMAKE_HOME_DIRECTORY:INTERNAL=//p' "$cmake_cache")"
    if [[
        -n "$cached_root"
        && "$cached_root" != "$repo_root"
        && "$cached_root" != "$repo_root/"*
    ]]; then
        stale_build_root="$cached_root"
        break
    fi
done < <(find build -type f -name CMakeCache.txt -print0 2>/dev/null)
if [[ -n "$stale_build_root" ]]; then
    echo "removing stale CMake cache from moved checkout: $stale_build_root"
    cmake -E remove_directory "$repo_root/build"
fi

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
    build \
    mypy \
    numpy \
    pytest \
    ruff

bash .codex/environments/install-worktree-wheel.sh

.venv/bin/python - <<'PY'
from importlib import metadata
import json
from pathlib import Path
from urllib.parse import unquote, urlparse

import epcsaft
from epcsaft import equilibrium, regression

distribution = metadata.distribution("epcsaft")
direct_url = json.loads(distribution.read_text("direct_url.json") or "{}")
if direct_url.get("dir_info", {}).get("editable", False):
    raise SystemExit("editable Engine install is forbidden")
installed_url = direct_url.get("url", "")
installed_wheel = Path(unquote(urlparse(installed_url).path)).resolve()
if not installed_wheel.is_file() or installed_wheel.suffix != ".whl":
    raise SystemExit(f"Engine is not installed from a worktree wheel: {installed_url}")
installed_root = Path(distribution.locate_file("")).resolve()
for module in (epcsaft, equilibrium, regression):
    module_path = Path(module.__file__).resolve()
    if installed_root not in module_path.parents:
        raise SystemExit(f"{module.__name__} is outside the installed wheel: {module_path}")
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
