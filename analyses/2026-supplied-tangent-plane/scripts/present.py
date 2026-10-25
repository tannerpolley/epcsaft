"""Present retained #300 verification rows; no Engine calls."""
import csv
import json
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
RESULTS = HERE / 'results'
rows = list(csv.DictReader((RESULTS/'qualification.csv').open()))
groups = defaultdict(list)
for row in rows:
    groups[row['item']].append(row)
lines = ['| Installed-wheel check | Count | Maximum absolute error | Maximum error / criterion | Passed |',
         '|---|---:|---:|---:|---|']
for name, items in groups.items():
    error = max(float(r['error']) for r in items)
    ratio = max(float(r['error'])/float(r['criterion']) if float(r['criterion']) else 0. for r in items)
    lines.append(f'| {name} | {len(items)} | {error:.17g} | {ratio:.17g} | {all(r["passed"] == "True" for r in items)} |')
native = list(csv.DictReader((RESULTS/'native-qualification.csv').open()))
lines += ['', '| Native fixture check | Count | Maximum absolute error | Maximum error / criterion | Passed |',
          '|---|---:|---:|---:|---|']
groups = defaultdict(list)
for row in native:
    name = row['comparison_id'].split('/')[2]
    if name in ('public', 'corrected'): name += '/'+row['comparison_id'].split('/')[-1]
    groups[name].append(row)
for name, items in groups.items():
    error = max(float(r['error']) for r in items)
    ratio = max(float(r['error'])/float(r['criterion']) if float(r['criterion']) else 0. for r in items)
    lines.append(f'| {name} | {len(items)} | {error:.17g} | {ratio:.17g} | {all(r["passed"] == "True" for r in items)} |')
(RESULTS/'qualification-table.md').write_text('\n'.join(lines)+'\n')
provenance = json.loads((RESULTS/'provenance.json').read_text())
key = json.loads((RESULTS/'key-solves.json').read_text())
lines = [f'All **{provenance["default_rows"]:,} retained comparison rows** match main bitwise; explicitly false matches omission. This includes candidate fields, trial order, raw chart/residual/composition arrays, optimizer values, quantity maps, donor/split arrays and work counts.',
         '', f'Main wheel SHA-256: `{provenance["main_wheel_sha256"]}`.',
         f'Candidate wheel SHA-256: `{provenance["candidate_wheel_sha256"]}`.', '',
         '| Declared solve | Main / candidate work (iterations / evaluations / Hessians) | Maximum relative state difference | Paired median time change (%) | Paired range (%) |',
         '|---|---|---:|---:|---|']
for name in ('case-6/figure-2-point-7/cold','case-7/flash-5918.5286-Pa/cold'):
    c = key['cases'][name]
    work = c['work_A'] if 'work_A' in c else c
    lines.append(f'| {name} | {work} | {c["max_relative_state_difference"]:.17g} | {c["median_change_percent"]:.17g} | {c["range_percent"]} |')
(RESULTS/'default-table.md').write_text('\n'.join(lines)+'\n')
