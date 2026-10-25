"""Reproduce the #187 cold free-temperature branch-window failure; no solver changes."""
import json
import math
import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "engine"))
from tests.regression.test_fit import neutral_binary, bubble
from epcsaft import Mixture
from epcsaft.equilibrium import Free, Reaction, ReactionLogPolynomial, EOS_STANDARD_STATE_ID, solve_equilibrium

p = neutral_binary()
model = Mixture(p)
problem = bubble(180.0)


def declarations():
    return model, problem


if '--declare-only' in sys.argv:
    raise SystemExit(0)
base = solve_equilibrium(model, problem)
x = base.mole_fractions[:2]
state = model.state(180.0, rho=base.molar_densities[0], x=x)
L = math.log(x[1]/x[0]) + state.residual_chemical_potential_over_rt[1] - state.residual_chemical_potential_over_rt[0]
a = L - 300.0/180.0 - 0.5*math.log(180.0/200.0) - 0.001*180.0
problem.reactions = [Reaction({"s0": -1.0, "s1": 1.0}, name="iso", correlation=
    ReactionLogPolynomial(a, 300.0, 0.5, 0.001, 200.0, 150.0, 220.0, EOS_STANDARD_STATE_ID))]
reacting = solve_equilibrium(model, problem)
problem.T = Free(180.0)
problem.P = reacting.pressure

def record(result):
    return {"success": result.success, "solver_status": result.solver_status,
            "message": result.message, "topology": str(result.topology_event),
            "T_K": result.temperature, "P_Pa": result.pressure,
            "x": result.mole_fractions, "packing": result.packing_fractions,
            "residual": result.convergence.residual}

cold = solve_equilibrium(model, problem)
for q, phase in enumerate(problem.phases):
    phase.composition_guess = tuple(reacting.mole_fractions[q*2:(q+1)*2])
    phase.packing_guess = reacting.packing_fractions[q]
seeded = solve_equilibrium(model, problem)
output = {
    "source": "engine/tests/regression/test_fit.py:neutral_binary and bubble",
    "inputs": {"molar_mass_g_mol": [16.043,30.07], "m": [1.,1.6069],
               "sigma_angstrom": [3.7039,3.5206], "epsilon_K": [150.03,191.42],
               "kij": 0., "feed_mol": [.5,.5], "T_guess_K": 180.,
               "P_fixed_Pa": reacting.pressure, "reaction": "iso",
               "stoichiometry": [-1.,1.], "a": a, "b_K": 300., "c": .5,
               "d_per_K": .001, "Tref_K": 200., "Tmin_K": 150., "Tmax_K": 220.,
               "standard_state_id": EOS_STANDARD_STATE_ID,
               "phases": ["liquid, free amount", "vapor, pinned zero amount"],
               "cold_composition_guesses": [], "cold_packing_guesses": [None,None],
               "thermochemistry": None, "source_reference": None,
               "solver": "unchanged Engine defaults"},
    "J1": record(reacting), "cold": record(cold), "declared_J1": record(seeded)}
(ROOT / "analyses/2026-reaction-constant-fitting/results/cold-start.json").write_text(
    json.dumps(output, indent=2) + "\n")
