import csv
import json
import sys
from pathlib import Path
root = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(root/'analyses/2026-figiel-divalent-chlorides/scripts'))
import stage2_fit as S
from epcsaft import _core
reference = S.OUT
output = reference/'current-wheel'
output.mkdir(exist_ok=True)
S.OUT = S.F.OUT = output
source = S.M.assembled()
fixed = S.M.assembled(document=S.BASE)
_,_,obs,meta,_ = S.F.declarations(source, 'all')
new = _core.observations_at(fixed._native, [], [], [], obs, [0.]*len(obs))
old = S.F.read(reference/'stage2-combined-stage1-preservation.csv')
assert len(new)==len(old)==948
rows = [v|dict(reference=float(a['reference']),current=b.observed,error=abs(b.observed-float(a['reference'])),
    tolerance=float(a['tolerance']),passed=abs(b.observed-float(a['reference']))<=float(a['tolerance']))
    for v,a,b in zip(meta,old,new)]
S.write('historical-stage1-preservation', rows)
print('Historical stage1',len(rows),'failures',sum(not v['passed'] for v in rows),'maximum change',max(v['error'] for v in rows),flush=True)
assert all(v['passed'] for v in rows), 'Current wheel changes historical stage-1 predictions beyond 1e-12; stop before fitting or packet.'
S.guards(S.BASE,'baseline')
for salt in S.SALTS:
    S.fit(salt,'ceres')
    S.fit(salt,'minimax')
    S.qualify(salt)
S.finish()
S.diagnose()
