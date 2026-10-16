#!/usr/bin/env bash
# Check one Engine wheel where no host library exists: a digest-pinned python:3.13-slim container
# (Debian 13) installs it beside this commit's tracked tree and runs verify-installed-wheel.py and the
# Python suite. Fails first if any member names HSL: Ipopt loads it at run time, never from the wheel.
#   bash tools/check-wheel-container.sh <wheel>
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
wheel="$(realpath "$1")"
python3 - "$wheel" <<'PY'
import re, sys, zipfile
hsl = [name for name in zipfile.ZipFile(sys.argv[1]).namelist() if re.search("hsl|ma27|ma57", name, re.I)]
sys.exit(f"HSL members in wheel: {hsl}" if hsl else None)
PY
pytest="$(python3 -c "import tomllib; print(next(p['version'] for p in tomllib.load(open('uv.lock', 'rb'))['package'] if p['name'] == 'pytest'))")"
git archive HEAD | docker run --rm -i -v "$wheel:/wheel/${wheel##*/}:ro" \
    python:3.13-slim@sha256:70729b46c69b4f1e97c4822c1af3df53a1476cf5ddc6c087c0c10bc3a5678c2f bash -euo pipefail -c '
mkdir /src && tar -x -C /src && cd /src/engine
pip install -q --root-user-action=ignore "pytest==$1" /wheel/*.whl
python ../tools/verify-installed-wheel.py --wheel /wheel/*.whl
python -m pytest -q -p no:cacheprovider' _ "$pytest"
