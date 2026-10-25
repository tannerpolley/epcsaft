"""Retain #210 numerical comparisons, with each row stored once."""
import argparse
import csv
import hashlib
import json
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
STUDY = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--zero-inventory', action='store_true', help='retain #231 rows separately from the existing #210 evidence')
parser.add_argument('--stem', help='retain a replay separately from historical results')
parser.add_argument('--base', help='explicit implementation base for a replay')
parser.add_argument('--baseline-probe', type=Path)
parser.add_argument('--baseline-preparation', type=Path)
parser.add_argument('--candidate-probe', type=Path)
parser.add_argument('--failure-log', type=Path, action='append', default=[])
parser.add_argument('--artifact', type=Path, action='append', default=[], help='retain supporting gate/output identities without duplicating numerical rows')
parser.add_argument('--preparation', type=Path, help='preparation identity for the retained numerical execution')
parser.add_argument('logs', nargs='+', type=Path)
args = parser.parse_args()
rows = []
for path in args.logs:
    for line in path.read_text().splitlines():
        if not line.startswith('EPCSAFT_COMPARE_JSON '):
            continue
        item = json.loads(line.split(' ', 1)[1])
        key = item['comparison_id']
        prefixes = ('inventory-face/', 'neutral-electrostatics/', 'zero-inventory/', 'incipient/') if args.zero_inventory else (
            'addition-', 'addition/', 'temperature-parameter', 'enthalpy-parameter', 'enthalpy-reaction',
            'reference-calorics/parameter/', 'addition-enthalpy', 'reference-parameter-calorics/')
        if not any(token in key for token in prefixes):
            continue
        criterion = item.get('criterion', {})
        if not isinstance(criterion, dict):
            criterion = {name: item.get(name, 0) for name in ('atol', 'rtol')}
        model, reference = item.get('model'), item.get('reference')
        # #231 judgments and comparison arithmetic are native; Python only retains them.
        error = None if args.zero_inventory else abs(model - reference) if model is not None and reference is not None else None
        atol, rtol = criterion.get('atol', 0), criterion.get('rtol', 0)
        rows.append(dict(comparison=key, quantity=item.get('quantity', ''), unit=item.get('unit', ''),
            model=model, reference=reference, absolute_error=error, atol=atol, rtol=rtol,
            passed=item.get('verdict', None if item.get('judgment') == 'Observation only' else item.get('judgment') == 'Criterion satisfied'),
            inputs=json.dumps({key: value for key, value in item.get('inputs', {}).items() if key != 'gradient'}, separators=(',', ':')),
            criterion=json.dumps(item.get('rule', item.get('criterion_rule', item.get('criterion', {}))), separators=(',', ':'))))
if not rows:
    raise SystemExit('No #210 numerical rows found')
results = STUDY / 'results'
results.mkdir(exist_ok=True)
stem = args.stem or ('zero-inventory' if args.zero_inventory else 'comparisons')
with (results / f'{stem}.csv').open('w', newline='') as out:
    writer = csv.DictWriter(out, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
preparation = json.loads((args.preparation or ROOT / 'build/environment-wheel/preparation.json').read_text())
summary = dict(base=args.base or ('a0adb1cd5ff5f4ee32b3608c78a02c756b06ebd9' if args.zero_inventory else '3f5ee19f'), wheel_sha256=preparation['wheel_sha256'],
    logs=[str(path.relative_to(ROOT) if path.is_relative_to(ROOT) else path) for path in args.logs],
    rows=len(rows), failed=[row['comparison'] for row in rows if row['passed'] is False],
    evidence='Numerical verification on manufactured Engine cases; no physical validation or parameter adoption')
if args.zero_inventory:
    summary['native_binary_sha256'] = hashlib.sha256((ROOT / 'build/native-drivers/zero_inventory/driver').read_bytes()).hexdigest()
    summary['log_sha256'] = {str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in args.logs}
    summary['artifacts'] = {str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in args.artifact}
    if args.baseline_preparation:
        summary['baseline_preparation'] = json.loads(args.baseline_preparation.read_text())
    summary['inputs'] = {path.name: hashlib.sha256(path.read_bytes()).hexdigest()
                         for path in (ROOT / 'engine/data/reference/mea-zero-inventory').glob('*.json')}
    summary['source_head'] = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()
    summary['runtime_patch_sha256'] = hashlib.sha256(subprocess.check_output(
        ['git', 'diff', summary['base'], '--', 'engine/native', 'engine/src'], cwd=ROOT)).hexdigest()
    summary['native_test_patch_sha256'] = hashlib.sha256(subprocess.check_output(
        ['git', 'diff', summary['base'], '--', 'engine/tests/native'], cwd=ROOT)).hexdigest()
    summary['probes'] = {name: json.loads(path.read_text()) for name, path in
                         (('baseline', args.baseline_probe), ('candidate', args.candidate_probe)) if path is not None}
    summary['failures'] = {str(path): '\n'.join(line for line in path.read_text().splitlines()
                          if not line.startswith('EPCSAFT_COMPARE_JSON ')) for path in args.failure_log}
    summary['evidence'] = 'Numerical verification of exact inventory faces and neutral electrostatic restrictions; no calorimetry validation, fit or parameter adoption'
(results / (f'{stem}.json' if args.stem else 'zero-inventory.json' if args.zero_inventory else 'summary.json')).write_text(json.dumps(summary, indent=2) + '\n')
print(json.dumps(summary, indent=2))
