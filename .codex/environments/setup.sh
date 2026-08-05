#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$repo_root"

for required_command in git uv cmake ninja c++; do
    if ! command -v "$required_command" >/dev/null 2>&1; then
        echo "missing required command: $required_command" >&2
        exit 1
    fi
done

uv sync --python 3.13 --no-install-project

wheelhouse="$repo_root/build/environment-wheel"
cmake -E remove_directory "$wheelhouse"
cmake -E make_directory "$wheelhouse"
.venv/bin/python -m build --wheel --outdir "$wheelhouse" engine

mapfile -t wheels < <(find "$wheelhouse" -maxdepth 1 -type f -name 'epcsaft-*.whl' -print)
if ((${#wheels[@]} != 1)); then
    echo "expected exactly one Engine wheel, found ${#wheels[@]}" >&2
    exit 1
fi
uv pip install --python .venv/bin/python --reinstall "${wheels[0]}"

.venv/bin/python - "$repo_root" "${wheels[0]}" <<'PY'
from importlib import metadata
import json
from pathlib import Path
import re
import sys
from urllib.parse import unquote, urlparse

root = Path(sys.argv[1]).resolve()
expected_wheel = Path(sys.argv[2]).resolve()
runtime = []
for distribution in metadata.distributions():
    name = distribution.metadata.get("Name", "")
    normalized = re.sub(r"[-_.]+", "-", name).lower()
    if normalized.startswith("epcsaft"):
        runtime.append((name, distribution))
if [name for name, _ in runtime] != ["epcsaft"]:
    raise SystemExit(f"expected only epcsaft, found {[name for name, _ in runtime]}")
distribution = runtime[0][1]
direct_url = json.loads(distribution.read_text("direct_url.json") or "{}")
if direct_url.get("dir_info", {}).get("editable", False):
    raise SystemExit("editable Engine install is forbidden")
installed = Path(unquote(urlparse(direct_url.get("url", "")).path)).resolve()
if installed != expected_wheel or root / "build" / "environment-wheel" not in installed.parents:
    raise SystemExit(f"unexpected Engine install source: {installed}")
import epcsaft
from epcsaft import equilibrium, regression
installed_root = Path(distribution.locate_file("")).resolve()
for module in (epcsaft, equilibrium, regression):
    if installed_root not in Path(module.__file__).resolve().parents:
        raise SystemExit(f"{module.__name__} is outside the installed wheel")
print(f"epcsaft {metadata.version('epcsaft')} installed from {installed}")
PY
