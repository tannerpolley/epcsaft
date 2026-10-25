"""Use the existing persistent FeOs worker unchanged, with this leaf's declarations."""
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
ROOT = HERE.parents[1]
spec = importlib.util.spec_from_file_location('release_feos', ROOT / 'analyses/2026-audit-timing-campaign/scripts/feos_worker.py')
w = importlib.util.module_from_spec(spec)
spec.loader.exec_module(w)
w.INPUTS = json.loads((HERE / 'inputs/inputs.json').read_text())
w.T1 = {r['id']: r for r in w.INPUTS['t1']}
if __name__ == '__main__':
    w.main()
