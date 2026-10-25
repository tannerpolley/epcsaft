import json, sys
sys.path.insert(0, 'analyses/2026-head-to-head-comparisons/scripts')
import h2h_engine as h
from epcsaft import Mixture, Parameters
model = Mixture(Parameters.from_mapping(h.khudaida_document('matched'), components=h.KH_IDS))
rows = json.load(open('build/h2h/probe-engine-terms.json')); out = []
terms = ('hard_chain','dispersion','association','debye_huckel','born')
for r in rows:
    x, rho = r['x'], r['rho']; V = 1 / rho; d = {t: [] for t in terms}
    for i in range(5):
        vals = []
        for hstep in (1e-6, -1e-6):
            n = list(x); n[i] += hstep; N = sum(n)
            s = model.state(293.15, rho=N / V, x=[v / N for v in n])
            vals.append({t: N * getattr(s, t) for t in terms})
        for t in terms: d[t].append((vals[0][t] - vals[1][t]) / 2e-6)
    out.append(d)
json.dump(out, open('build/h2h/probe-engine-mu-terms.json', 'w'))
