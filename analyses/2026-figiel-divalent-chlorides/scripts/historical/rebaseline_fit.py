import csv
import json
import subprocess
import sys
from pathlib import Path

root = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(root/'analyses/2026-figiel-divalent-chlorides/scripts'))
import stage2_fit as S
from epcsaft import _core

original = S.OUT
baseline = original/'stage1-rebaseline'
output = original/'provisional-stage2'
output.mkdir(exist_ok=True)
S.OUT = S.F.OUT = output
source = S.M.assembled()
fixed = S.M.assembled(document=S.BASE)
_, _, obs, meta, groups = S.F.declarations(source, 'all')
new = _core.observations_at(fixed._native, [], [], [], obs, [0.]*len(obs))
old = S.F.read(original/'stage2-combined-stage1-preservation.csv')
assert len(new) == len(old) == 948
rows = []
for v, a, b in zip(meta, old, new):
    for key in v:
        if key != 'observed':
            assert str(v[key]) == a[key], (key, v[key], a[key])
    assert abs(v['observed']-float(a['observed'])) <= 1e-12*max(abs(float(a['observed'])), 1.)
    rows.append(v | dict(original_reference=float(a['reference']), reference=b.observed,
        original_change=b.observed-float(a['reference']),
        tolerance=1e-12*abs(b.observed) if b.observed else 1e-12))
with (baseline/'stage1-reference.csv').open('w') as stream:
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator='\n')
    writer.writeheader()
    writer.writerows(rows)
record = dict(authority='Accountable lead decision in the 2026-10-06 build assignment: re-baseline current-wheel predictions; preserve /1 and 1e-12 rule.',
    reason='PR #273 density.cpp::refine retains the pressure-closest endpoint; reproduced bisection improves closure in every moved row, with individual Phi identity errors improving or worsening as retained in diagnosis.json.',
    reference='stage1-reference.csv', rows=len(rows),
    wheel_sha256=S.digest(next((root/'build/environment-wheel').glob('*.whl'))),
    source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=root,text=True).strip(),
    merged_origin_main=subprocess.check_output(['git','rev-parse','origin/main'],cwd=root,text=True).strip(),
    original_reference=str((original/'stage2-combined-stage1-preservation.csv').relative_to(root)),
    original_reference_sha256=S.digest(original/'stage2-combined-stage1-preservation.csv'),
    original_reference_wheel_sha256=json.loads((original/'stage2-summary.json').read_text())['wheel_sha256'],
    reference_sha256=S.digest(baseline/'stage1-reference.csv'), stage1_parameters_sha256=S.digest(S.SOURCE),
    density_correction='https://github.com/tannerpolley/ePC-SAFT/pull/273',
    diagnosis='diagnosis.json', diagnosis_sha256=S.digest(baseline/'diagnosis.json'),
    preservation_rule='abs(candidate-reference) <= 1e-12*abs(reference) for nonzero reference, otherwise <=1e-12.',
    evidence_type='Numerical preservation reference; not absolute 1e-12 accuracy or physical validation.',
    maximum_original_change=max(abs(v['original_change']) for v in rows),
    original_preservation_failures=sum(abs(v['original_change'])>float(a['tolerance']) for v,a in zip(rows,old)))
(baseline/'rebaseline.json').write_text(json.dumps(record,indent=2)+'\n')
print('Re-baseline retained', record, flush=True)
S.guards(S.BASE, 'baseline')
for salt in S.SALTS:
    S.fit(salt, 'ceres')
    S.fit(salt, 'minimax')
    S.qualify(salt)
doc = S.extended()
for salt in S.SALTS:
    report = json.loads((output/f'stage2-{salt}-qualification.json').read_text())
    selected = json.loads((output/f'stage2-{salt}-minimax.json').read_text())['results'][report['selected']]
    assert selected['minimax']['feasible'] and selected['minimax']['stationary']
    doc = S.candidate(salt, selected['physical'], doc)
doc = S.packet_document(doc)
S.save('combined-candidate', doc)
for salt in S.SALTS:
    S.qualify(salt, doc)
S.guards(doc, 'combined')
stored = S.F.read(output/'stage2-combined-stage1-preservation.csv')
for a,b in zip(rows, stored):
    assert abs(float(b['stored'])-a['reference']) <= a['tolerance']
print('Re-baselined preservation passed',len(stored), flush=True)
