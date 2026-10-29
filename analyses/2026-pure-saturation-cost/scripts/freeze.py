"""Reuse the campaign image freezer for the selected clean prepared checkout."""
import ast
import importlib.util
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[3]
path = ROOT / 'analyses/2026-audit-timing-campaign/scripts/freeze.py'
spec = importlib.util.spec_from_file_location('image_freezer', path)
f = importlib.util.module_from_spec(spec)
spec.loader.exec_module(f)
f.ROOT, f.HERE = Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve()
(f.HERE / 'results').mkdir(parents=True, exist_ok=True)
# This calculation has no fits; retain the original image freeze and omit its fit launcher.
tree = ast.parse(path.read_text())
fn = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'main')
end = next(i for i, n in enumerate(fn.body) if isinstance(n, ast.Assign) and any(isinstance(t, ast.Name) and t.id == 'launcher' for t in n.targets))
fn.body = fn.body[:end]
exec(compile(ast.fix_missing_locations(ast.Module(body=[fn], type_ignores=[])), str(path), 'exec'), f.__dict__)
f.main()
