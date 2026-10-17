"""The #204 free-temperature incipient vapor must recover the retained J1 state."""
import json
from pathlib import Path
from typing import Any

from epcsaft import Mixture
from epcsaft.equilibrium import Free, Reaction, ReactionLogPolynomial, solve_equilibrium

from tests.native_inputs import bubble, parameters, same


def test_cold_free_temperature_vapor_recovers_j1(evidence: Any) -> None:
    case = json.loads((Path(__file__).resolve().parents[3] / "analyses" /
                       "2026-reaction-constant-fitting/results/cold-start.json").read_text())
    inputs, reference = case["inputs"], case["J1"]
    model = Mixture(parameters(list(zip(inputs["molar_mass_g_mol"], inputs["m"],
                                        inputs["sigma_angstrom"], inputs["epsilon_K"])),
                               pairs=[("k_ij", (0, 1), inputs["kij"])]))
    problem = bubble(inputs["T_guess_K"])
    problem.T = Free(inputs["T_guess_K"])
    problem.P = inputs["P_fixed_Pa"]
    problem.reactions = (Reaction(dict(zip(model.component_ids, inputs["stoichiometry"])),
                                 name=inputs["reaction"], correlation=ReactionLogPolynomial(
                                     inputs["a"], inputs["b_K"], inputs["c"], inputs["d_per_K"],
                                     inputs["Tref_K"], inputs["Tmin_K"], inputs["Tmax_K"],
                                     inputs["standard_state_id"])),)
    result = solve_equilibrium(model, problem)
    assert result.success, result.message
    assert result.requested_tolerance_met
    # Existing J1 comparisons: 1e-10 absolute for x, relative for T/P/eta.
    for quantity, values, expected in (("T_K", [result.temperature], reference["T_K"]),
                                     ("P_Pa", [result.pressure], reference["P_Pa"]),
                                     ("x", result.mole_fractions, reference["x"]),
                                     ("packing", result.packing_fractions, reference["packing"])):
        for index, (value, target) in enumerate(zip(values, expected, strict=True)):
            same(evidence, f"cold-vapor/{quantity}/{index}", value, target,
                 atol=1e-10 if quantity == "x" else 0.0,
                 rtol=0.0 if quantity == "x" else 1e-10)
