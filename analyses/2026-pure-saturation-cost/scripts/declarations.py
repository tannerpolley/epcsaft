"""Declare the unchanged C1-300 request without a numerical solve."""
import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]


def declarations():
    spec = importlib.util.spec_from_file_location('cost_inputs', ROOT / 'analyses/2026-competitor-cross-check/scripts/engine.py')
    cc = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cc)
    inputs = json.loads((ROOT / 'analyses/2026-competitor-cross-check/inputs/cases.json').read_text())
    request = next(r for r in inputs['requests'] if r['id'] == 'C1-300')
    return [cc.prepare(inputs, request, declare_only=True)[:2]]
