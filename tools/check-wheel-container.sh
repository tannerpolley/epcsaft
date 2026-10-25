#!/usr/bin/env bash
# Check one Engine wheel where no host library exists: a digest-pinned python:3.13-slim container
# (Debian 13) installs it with the locked test environment beside this commit's tracked tree, then runs
# verify-installed-wheel.py and the Python suite. Fails first if any member names HSL: Ipopt loads it
# at run time, never from the wheel.
#   bash tools/check-wheel-container.sh <wheel>
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
wheel="$(realpath "$1")"
python3 - "$wheel" <<'PY'
import re, sys, zipfile
hsl = [name for name in zipfile.ZipFile(sys.argv[1]).namelist() if re.search("hsl|ma27|ma57", name, re.I)]
sys.exit(f"HSL members in wheel: {hsl}" if hsl else None)
PY
# The suite's locked test environment (pytest, pyomo, ...), installed by the host's uv.
requirements="$(mktemp)"
trap 'rm -f "$requirements"' EXIT
uv export --locked --no-emit-workspace --no-hashes --all-packages --quiet >"$requirements"
git archive HEAD | docker run --rm -i -v "$wheel:/wheel/${wheel##*/}:ro" -v "$requirements:/requirements.txt:ro" \
    -v "$(command -v uv):/usr/local/bin/uv:ro" \
    python:3.13-slim@sha256:70729b46c69b4f1e97c4822c1af3df53a1476cf5ddc6c087c0c10bc3a5678c2f bash -euo pipefail -c '
mkdir /src && tar -x -C /src && cd /src/engine
uv pip install -q --system -r /requirements.txt /wheel/*.whl
python ../tools/verify-installed-wheel.py --wheel /wheel/*.whl
python -m pytest -q -p no:cacheprovider'
