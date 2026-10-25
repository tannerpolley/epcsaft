import copy
import json
import math
import sys
from pathlib import Path
sys.path.insert(0, str(Path.cwd() / 'engine'))
from tests.native_inputs import ETHANE, parameters
from epcsaft import Mixture, Parameters, _core
from epcsaft.equilibrium import *
from read_capture import b
cases = [(name.replace('/', '-'), model, problem) for name, model, problem in b.declared_cases() if name in ('case-1/two-liquid', 'case-6/figure-2-point-7', 'case-7/flash-5918.5286-Pa')]
# Frozen copy of the #172 masking declaration in engine/tests/equilibrium/test_trace_reliability.py
# (test_trace_inventory_below_the_composition_box, inventory 1e-12 mol) as run at 9fdfe7ce.
model = Mixture(parameters((ETHANE,) * 3, [("k_ij", pair, 0.0) for pair in ((0, 1), (0, 2), (1, 2))]))
temperature, k = 350.0, 3.0
a = math.log(k) - 1200 / temperature - 0.4 * math.log(temperature / 300) - 0.002 * temperature
law = ReactionLogPolynomial(a, 1200, 0.4, 0.002, 300, 250, 450, EOS_STANDARD_STATE_ID)
problem = Problem(phases=[Phase("V", kind="vapor")], T=temperature, P=2e5, row_tolerance=1e-10,
                  feed=Amounts({"s0": 1.0, "s1": 1e-12, "s2": 0.0}),
                  reactions=[Reaction({"s1": -1.0, "s2": 1.0}, correlation=law)])
cases.append(('masking', model, problem))
name, model, problem = cases[0]
problem = copy.copy(problem); problem.phases = Discover(); problem.topology_declared = False
cases.append(('discovery', model, problem))
site = json.loads(Path('analyses/2026-native-diagnostic-trace/inputs/site-trace-input.json').read_text())
inputs = site['parameters']
site_parameters = parameters(tuple(zip((v * 1000 for v in inputs['molar_mass_kg_per_mol']),
    inputs['segment_count'], inputs['sigma_constant'], inputs['dispersion_energy_over_k'])),
    [('k_ij', (0, 1), inputs['k_ij'][1]),
     ('association_energy_over_k', (0, 1), inputs['association_energy_over_k'][1]),
     ('association_volume', (0, 1), inputs['association_volume'][1])],
    sites=(inputs['site_component'], inputs['site_multiplicity']))
site_model = Mixture(site_parameters)


def declarations():
    return (cases, site_model)
