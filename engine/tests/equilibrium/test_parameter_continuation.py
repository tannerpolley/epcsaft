"""Identical pseudo-species cancel residual potentials in the two exchange laws (#234)."""
import json
from typing import Any

from epcsaft import ActiveParameter, Mixture, ParameterFamily, Parameters, _core, regression
from epcsaft import equilibrium as eq
from tests.native_inputs import ETHANE, parameters, same
from tests.packets import parameters_json


def test_dispersion_path_preserves_exchange_amounts_and_compositions(evidence: Any) -> None:
    p = parameters((ETHANE,) * 6, pairs=[("k_ij", (i, j), 0.0) for i in range(6) for j in range(i + 1, 6)])
    active = [ActiveParameter(ParameterFamily.DispersionEnergyOverK, (i,)) for i in range(6)]
    model = Mixture(p, active_parameters=active)
    problem = eq.Problem(T=200.0, P=5e6,
        phases=[eq.Phase("aqueous", kind="liquid", support=("s0", "s2", "s4"), composition_guess=(10/13, 1/13, 2/13)),
                eq.Phase("organic", kind="liquid", support=("s1", "s3", "s5"), composition_guess=(10/13, 2/13, 1/13))],
        feed=eq.Amounts({"s0": 1.0, "s1": 1.0, "s2": 0.3, "s4": 0.3}),
        reactions=[eq.Reaction({"s2": -1.0, "s3": 1.0}, name="exchange_a", equilibrium_constant=2.0),
                   eq.Reaction({"s4": -1.0, "s5": 1.0}, name="exchange_b", equilibrium_constant=0.5)])
    result = eq.continue_equilibrium(model, problem, problem,
        eq.Continuation(initial_step=0.25, minimum_step=0.001, max_steps=20, max_retries_per_step=8),
        parameter_targets=(201.42,) * 6)
    assert result.success, result.message
    expected = (10/13, 0.0, 1/13, 0.0, 2/13, 0.0, 0.0, 10/13, 0.0, 2/13, 0.0, 1/13)
    for index, state in enumerate(result.states):
        for phase, amount in enumerate(state.amounts):
            same(evidence, f"parameter-path/{index}/amount/{phase}", amount, 1.3, atol=1e-9, rtol=0.0)
        for species, value in enumerate(state.mole_fractions):
            same(evidence, f"parameter-path/{index}/x/{species}", value, expected[species], atol=1e-9, rtol=0.0)
        same(evidence, f"parameter-path/{index}/rows", max(map(abs, state.residuals)), 0.0, atol=1e-10, rtol=0.0)
        assert max(map(abs, state.classified_residuals)) <= 1e-10
    compiled = eq.compile_problem(model, problem)
    observables = [_core.SolvedStateObservable(_core.SolvedStateObservableKind.PhaseAmount, phase) for phase in (0, 1)]
    observables += [_core.SolvedStateObservable(_core.SolvedStateObservableKind.PhaseMoleFraction, phase, species)
                    for phase, support in enumerate(((0, 2, 4), (1, 3, 5))) for species in support]
    action = _core.solved_state_first_actions(compiled, result.states[0], observables,
        [_core.SolvedStateActionDirection(0.0, 0.0, (0.0,) * 6, (10.0,) * 6)])[0]
    for index, item in enumerate(action.results):
        same(evidence, f"parameter-path/tangent/{index}", item.action, 0.0, atol=1e-9, rtol=0.0)


def test_engine_readme_parameter_path_reaches_direct_solve() -> None:
    """engine/README.md water/isobutanol record, k_ij and isobutanol eps/k active, fixed step 0.25. Each step's
    L1 density moves 1.03 %, beyond the 0.00995 ln-density window around the previous accepted density;
    the density at the predicted coordinates, which the exact tangent misses by O(h^2), also witnesses the branch."""
    mapping = json.loads(parameters_json("khudaida-2026-figure-2-electrolyte-lle", 7).read_text())

    def model(theta: tuple[float, ...] | None = None) -> Mixture:
        if theta is not None:
            pair = next(p for p in mapping["pairs"] if p["coefficients"][0]["identity"] == "pair/isobutanol/water/k_ij")
            pair["coefficients"][0]["value"]["magnitude"] = theta[0]
            iso = next(c for c in mapping["components"] if c["component_id"] == "isobutanol")
            next(k for k in iso["coefficients"] if k["identity"].endswith("dispersion_energy_over_k"))["value"]["magnitude"] = theta[1]
        p = Parameters.from_mapping(mapping, components=("water", "isobutanol"))
        return Mixture(p, active_parameters=[regression.active_parameter(p, "k_ij", ("isobutanol", "water")),
                                             regression.active_parameter(p, "dispersion_energy_over_k", "isobutanol")])

    both = ("water", "isobutanol")
    problem = eq.Problem(phases=[eq.Phase("L1", kind="liquid", support=both, composition_guess=(0.46, 0.54)),
                                 eq.Phase("L2", kind="liquid", support=both, composition_guess=(0.99, 0.01))],
                         T=298.15, P=1.0e5, feed=eq.Amounts({"water": 0.7, "isobutanol": 0.3}))

    def amounts(result: eq.EquilibriumResult) -> list[float]:
        return [result.amounts[p] * result.mole_fractions[2 * p + c] for p in range(2) for c in range(2)]

    base = model()
    theta0 = tuple(base.active_parameter_values)
    start = eq.solve_equilibrium(base, problem)
    for target in ((-0.0021, 257.13), (-0.0021, 258.13)):  # the k_ij shift +0.005 alone, then with eps/k +1 K
        path = eq.continue_equilibrium(base, problem, problem, eq.Continuation(
            initial_step=0.25, minimum_step=0.25, max_steps=4, max_retries_per_step=1), start, parameter_targets=target)
        assert path.success and path.accepted_steps == 4, path.message
        direct = eq.solve_equilibrium(model(target), problem)
        assert direct.success, direct.message
        assert max(abs(a / b - 1.0) for a, b in zip(amounts(path.states[-1]), amounts(direct), strict=True)) <= 1e-12
    compiled = _core.compile_problem(base._model, problem._native(base))
    tangent = _core.sensitivity(compiled, start.coordinates, [0.0, 0.0, 0.0, 0.0, *(b - a for a, b in zip(theta0, target))])

    def predictor_error(s: float) -> float:
        y = eq.solve_equilibrium(model(tuple(a + s * (b - a) for a, b in zip(theta0, target))), problem).coordinates
        return max(abs(v - v0 - s * t) for v, v0, t in zip(y, start.coordinates, tangent, strict=True))

    assert 3.5 < predictor_error(0.25) / predictor_error(0.125) < 4.5  # second order in the step
