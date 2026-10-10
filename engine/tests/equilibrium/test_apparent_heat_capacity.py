"""Installed-wheel transport of EqID: equilibrium_isobaric_heat_capacity."""
import pytest
from epcsaft import Mixture, ActiveParameter, ParameterFamily, regression
from epcsaft.equilibrium import SolvedStateObservable, SolvedStateObservableKind as Kind
from epcsaft.equilibrium import SolvedStateActionRequest, SolvedStateActionDirection
from epcsaft.equilibrium import compile_problem, solve_equilibrium, solved_state_actions, SolvedStateActionStatus
from tests.native_inputs import A_RECORD, dimer_problem, records, twin_ethane


@pytest.mark.parametrize("reference_pressure", [None, 5.0e6])
def test_apparent_heat_capacity_transport(reference_pressure: float | None) -> None:
    parameters = twin_ethane(A_RECORD, []).parameters
    model = Mixture(parameters, active_parameters=[ActiveParameter(ParameterFamily.DispersionEnergyOverK, [0])],
                    thermochemistry=records(A_RECORD, []))
    problem = dimer_problem(260.0, reference_pressure)
    observation = regression.observation(parameters, "equilibrium_isobaric_heat_capacity", problem,
                                        observed=1.0, form="difference", scale=1.0, mass_basis=False)
    prediction = regression.predict(model, observation)
    direction = SolvedStateActionDirection(0.0, 0.0, [0.0, 0.0], [1.0])
    state = solve_equilibrium(model, problem)
    action = solved_state_actions(compile_problem(model, problem), state, SolvedStateActionRequest(
        [SolvedStateObservable(Kind.EquilibriumIsobaricHeatCapacity)], [direction])).results[0]
    columns = regression.differentiate(model, observation)
    assert prediction.value is not None and action.value == prediction.value
    assert prediction.unit == action.unit == "J/(mol-feed K)"
    assert action.action is not None and columns.value[0] == action.action


def test_apparent_heat_capacity_preserves_value_on_support_insertion() -> None:
    model = twin_ethane(A_RECORD, A_RECORD)
    problem = dimer_problem(260.0, None, feed={"s0": 2.0, "s1": 0.0})
    problem.reactions = ()
    state = solve_equilibrium(model, problem)
    compiled = compile_problem(model, problem)
    observable = SolvedStateObservable(Kind.EquilibriumIsobaricHeatCapacity)
    inserted = solved_state_actions(compiled, state, SolvedStateActionRequest(
        [observable], [SolvedStateActionDirection(0.0, 0.0, [0.0, 1.0], [])])).results[0]
    regular = solved_state_actions(compiled, state, SolvedStateActionRequest(
        [observable], [SolvedStateActionDirection(0.0, 0.0, [0.0, 0.0], [])])).results[0]
    assert inserted.status == SolvedStateActionStatus.SupportChange and inserted.action is None
    assert inserted.value is not None and inserted.value == regular.value
