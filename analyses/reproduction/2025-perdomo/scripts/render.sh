#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../../../.."
mkdir -p build/issue-729/render
cp analyses/reproduction/2025-perdomo/notebook.qmd build/issue-729/render/notebook.qmd
.venv/bin/python -m ipykernel install --prefix "$PWD/build/issue-729/jupyter" --name python3
JUPYTER_PATH="$PWD/build/issue-729/jupyter/share/jupyter" \
QUARTO_PYTHON="$PWD/.venv/bin/python" quarto render \
  build/issue-729/render/notebook.qmd --to html --execute --execute-dir "$PWD" "$@"
.venv/bin/python analyses/reproduction/2025-perdomo/scripts/finalize-render.py
