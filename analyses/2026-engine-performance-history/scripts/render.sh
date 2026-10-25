#!/usr/bin/env bash
set -euo pipefail
repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
cd "$repo_root"
presentation_python="$repo_root/build/issue-549/presentation/bin/python"
export QUARTO_PYTHON="$presentation_python"
export QUARTO_JUPYTER_JULIA=""

"$presentation_python" - <<'PY'
import importlib.util
assert importlib.util.find_spec('epcsaft') is None, 'rendering environment must not contain the Engine'
PY
"$presentation_python" analyses/2026-engine-performance-history/scripts/plot.py
"$presentation_python" - <<'PY'
import hashlib, json, shutil
from pathlib import Path
p = Path('analyses/2026-engine-performance-history')
sources = [p/'notebook.qmd', *sorted((p/'results').glob('*')), *sorted((p/'figures').glob('*.svg'))]
hashes = {str(f): hashlib.sha256(f.read_bytes()).hexdigest() for f in sources if f.is_file()}
Path('build/issue-549/render-inputs.json').write_text(json.dumps(hashes, indent=2)+'\n')
render = Path('build/issue-549/render')
render.mkdir(parents=True, exist_ok=True)
for name in ('notebook.qmd', 'README.md'):
    shutil.copy2(p/name, render/name)
for name in ('results', 'figures', 'inputs'):
    shutil.copytree(p/name, render/name, dirs_exist_ok=True)
PY
quarto render "$repo_root/build/issue-549/render/notebook.qmd" --to html --execute
"$presentation_python" - <<'PY'
import hashlib, json, shutil, subprocess
from pathlib import Path
hashes = json.loads(Path('build/issue-549/render-inputs.json').read_text())
assert all(hashlib.sha256(Path(p).read_bytes()).hexdigest() == h for p, h in hashes.items()), 'render changed retained inputs'
rendered = Path('build/issue-549/render/notebook.html')
assert rendered.stat().st_size > 10000
generated_ignore = Path('analyses/2026-engine-performance-history/.gitignore')
if generated_ignore.exists() and generated_ignore.read_text() == '/.quarto/\n**/*.quarto_ipynb\n':
    generated_ignore.unlink()
main = Path(subprocess.check_output(['git', 'rev-parse', '--git-common-dir'], text=True).strip()).resolve().parent
destination = main/'build/issue-549/notebook.html'
destination.parent.mkdir(parents=True, exist_ok=True)
shutil.copy2(rendered, destination)
print(destination)
PY
