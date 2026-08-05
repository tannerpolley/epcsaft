#!/usr/bin/env bash
set -euo pipefail

repo_root="$(git rev-parse --show-toplevel)"
cd "$repo_root"

if [[ ! -x .venv/bin/python ]]; then
    echo "missing Engine environment; run .codex/environments/setup.sh" >&2
    exit 1
fi

wheelhouse="$repo_root/build/environment-wheel"
cmake -E remove_directory "$wheelhouse"
.venv/bin/python -m build --wheel --outdir "$wheelhouse"

mapfile -t worktree_wheels < <(find "$wheelhouse" -maxdepth 1 -type f -name 'epcsaft-*.whl')
if ((${#worktree_wheels[@]} != 1)); then
    echo "expected exactly one Engine worktree wheel, found ${#worktree_wheels[@]}" >&2
    exit 1
fi

uv pip install --python .venv/bin/python --reinstall "${worktree_wheels[0]}"
