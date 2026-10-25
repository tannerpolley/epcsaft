import json, sys
sys.path.insert(0, 'analyses/2026-head-to-head-comparisons/scripts')
import h2h_engine as h
from epcsaft import Mixture, Parameters
model = Mixture(Parameters.from_mapping(h.khudaida_document('matched'), components=h.KH_IDS))
base = json.load(open('build/h2h/probe-engine-terms.json'))[0]
states = []
for i, dx in ((None, 0), (0, 0.05), (1, 0.05), (2, 0.05)):
    n = list(base['x'])
    if i is not None: n[i] += dx
    N = sum(n); x = [v / N for v in n]; rho = base['rho'] * N
    s = model.state(293.15, rho=rho, x=x)
    states.append(dict(x=x, rho=rho, debye_huckel=s.debye_huckel, born=s.born, eps=s.bulk_relative_permittivity))
json.dump(states, open('build/h2h/probe-perturbed.json', 'w')); print(json.dumps(states))
