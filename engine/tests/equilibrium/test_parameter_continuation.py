"""Identical pseudo-species cancel residual potentials in the two exchange laws (#234)."""
from typing import Any

from epcsaft import ActiveParameter, Mixture, ParameterFamily, _core
from epcsaft import equilibrium as eq
from tests.native_inputs import ETHANE, parameters, same


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
