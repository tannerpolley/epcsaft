"""Present the retained declared references without running thermodynamics."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import statistics

import costs

HERE = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--qualification', type=Path, default=HERE/'results/qualification-8470b49a')
parser.add_argument('--output', type=Path, default=HERE/'results')
args = parser.parse_args()
args.output.mkdir(parents=True, exist_ok=True)
rows = list(csv.DictReader((HERE/'results/before/declared.csv').open()))
cases = dict.fromkeys(row['case'] for row in rows if row['quantity'] == 'objective')
lines = ['| Declared reference | G/RT (mol) | Iterations | Evaluations | Hessian calls | Maximum classified residual |',
         '|---|---:|---:|---:|---:|---:|']
for case in cases:
    values = {row['quantity']: row['value'] for row in rows if row['case'] == case}
    residual = max(abs(float(v)) for k, v in values.items() if k.startswith('classified_residuals/'))
    fields = [case, *(values[key] for key in ('objective', 'iterations', 'evaluations', 'hessian_calls')), f'{residual:.17g}']
    lines.append('| '+' | '.join(fields)+' |')
baseline = '\n'.join(lines)+'\n'
canonical = HERE/'results/baseline-summary.md'
if args.output != canonical.parent and canonical.read_text() == baseline:
    (args.output/'baseline-link.json').write_text(json.dumps(dict(path=os.path.relpath(canonical, args.output),
        sha256=hashlib.sha256(canonical.read_bytes()).hexdigest(), qualification_directory=str(args.qualification),
        comparison='byte-identical generated table linked once'), indent=2)+'\n')
else:
    (args.output/'baseline-summary.md').write_text(baseline)

lines = ['| Case and operation | Median of three (s) | Driver calls | Trial calls | Density roots | Density evaluations | Iterations | Evaluations | Hessian calls |',
         '|---|---:|---:|---:|---:|---:|---:|---:|---:|']
for name in ('neutral', 'feos', 'ionic', 'permittivity', 'reacting', 'gibbs-scaling'):
    path = args.qualification/f'{name}.json'
    if not path.exists():
        continue
    data = json.loads(path.read_text())
    observed = costs.read(path, data)
    for case, role in dict.fromkeys((row['case'], row['role']) for row in observed):
        samples = [row for row in observed if (row['case'], row['role']) == (case, role)]
        assert len(samples) == 3
        counts = [{key: row[key] for key in ('iterations', 'evaluations', 'hessian_calls')} for row in samples]
        assert counts[0] == counts[1] == counts[2], 'work differs across timing samples; retain its variation before presentation'
        work = samples[0].get('native_work', samples[0].get('trace_work'))
        assert all(row.get('native_work', row.get('trace_work')) == work for row in samples)
        calls = [str(int(work.get(key, 0))) for key in ('driver_invocations', 'trial_invocations', 'root_calls', 'root_evaluations')] if work else [str(samples[0]['driver_finalizations']), '0', '—', '—']
        label = case+' / '+role+(' (from failed qualification run)' if data['failure'] is not None else '')
        fields = [label, repr(statistics.median(row['seconds'] for row in samples)), *calls,
                  *(str(int(counts[0][key])) for key in ('iterations', 'evaluations', 'hessian_calls'))]
        lines.append('| '+' | '.join(fields)+' |')
(args.output/'cost-summary.md').write_text('\n'.join(lines)+'\n')
