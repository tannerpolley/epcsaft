import json, sys
from epcsaft import Mixture, Parameters
from epcsaft import _core
from epcsaft.equilibrium import Problem, Phase, Pinned, Free, Amounts, solve_equilibrium
P = "data/packets/gross-2001-methane-ethane/3/parameters/parameters.json"
params = Parameters.from_json(P, components=("methane", "ethane"))
model = Mixture(params)
def problem(n_ethane):
    return Problem(phases=[Phase("L", kind="liquid"), Phase("V", kind="vapor", amount=Pinned(0.2))],
                   T=180.0, P=Free(), feed=Amounts({"methane": 0.5, "ethane": n_ethane}))
def declarations():
    return model, problem(0.5)
if '--declare-only' in sys.argv:
    declarations()
    raise SystemExit(0)
base = solve_equilibrium(model, problem(0.5))
assert base.success, base.message
compiled = _core.compile_problem(model._model, problem(0.5)._native(model))
C = 2; A = 0
lever = [0.0] * (2 + C + A); lever[2 + 1] = 1.0   # ethane feed column
d = _core.sensitivity(compiled, base.coordinates, lever)
p_index = compiled.primal_size - 1                 # P is the last free intensive
analytic = d[p_index]
h = 1e-4
up = solve_equilibrium(model, problem(0.5 + h), base); dn = solve_equilibrium(model, problem(0.5 - h), base)
assert up.success and dn.success
central = (up.pressure - dn.pressure) / (2 * h)
print(json.dumps({"case": "methane-ethane 180 K, V pinned 0.2, P free, dP/dn_ethane",
                  "P_pa": base.pressure, "analytic_pa_per_mol": analytic, "central_pa_per_mol": central,
                  "relative_difference": abs(analytic - central) / abs(central)}, indent=1))
