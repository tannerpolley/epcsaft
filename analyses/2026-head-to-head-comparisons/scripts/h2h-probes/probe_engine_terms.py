import json, sys
sys.path.insert(0, 'analyses/2026-head-to-head-comparisons/scripts')
import h2h_engine as h
from epcsaft import Mixture, Parameters
e = json.load(open('build/h2h/engine-h2h.json'))
model = Mixture(Parameters.from_mapping(h.khudaida_document('matched'), components=h.KH_IDS))
rows = []
for r in e['khudaida']['matched'][:1]:
    for x, rho in zip(r['x'], r['rho']):
        s = model.state(293.15, rho=rho, x=x)
        rows.append(dict(x=x, rho=rho, **{k: getattr(s, k) for k in ('hard_chain','dispersion','association','debye_huckel','born','bulk_relative_permittivity','residual_helmholtz','pressure')}, lnphi=list(s.log_fugacity_coefficient)))
json.dump(rows, open('build/h2h/probe-engine-terms.json','w'), indent=1); print(json.dumps(rows, indent=1))
