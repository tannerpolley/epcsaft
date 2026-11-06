"""Internal-energy and volume inputs and exact energy action on a homogeneous vapor."""
import pytest
from epcsaft.equilibrium import (Amounts, Free, Phase, Problem, SolvedStateActionDirection,
    SolvedStateActionRequest, SolvedStateObservable, SolvedStateObservableKind, compile_problem,
    problem_from_mapping, solve_equilibrium, solved_state_actions)
from tests.native_inputs import A_RECORD, twin_ethane


def test_internal_energy_volume_and_temperature_action() -> None:
    model = twin_ethane(A_RECORD, A_RECORD)
    state = model.state(300.0, rho=100.0, x=[1.0, 0.0])
    assert state.total_internal_energy is not None
    assert state.total_isochoric_heat_capacity is not None
    problem = problem_from_mapping({
        "T": {"free": True, "guess": 280.0}, "P": {"free": True, "guess": 1e5}, "internal_energy_j": state.total_internal_energy,
        "volume_m3": 0.01, "feed": {"basis": "amounts", "values": {"s0": 1.0, "s1": 0.0}},
        "phases": [{"name": "V", "kind": "vapor"}], "row_tolerance": 1e-11,
    })
    result = solve_equilibrium(model, problem)
    assert result.temperature == pytest.approx(300.0, rel=1e-8), result.message
    assert result.pressure == pytest.approx(state.pressure, rel=1e-8), result.message
    request = SolvedStateActionRequest([SolvedStateObservable(SolvedStateObservableKind.Temperature)],
        [SolvedStateActionDirection(0.0, 0.0, [0.0, 0.0], [], internal_energy_j=1.0)])
    batch = solved_state_actions(compile_problem(model, problem), result, request)
    assert batch.results[0].action == pytest.approx(1.0 / state.total_isochoric_heat_capacity, rel=1e-8)
    with pytest.raises(ValueError, match="UV requires"):
        solve_equilibrium(model, Problem(T=Free(280.0), P=Free(1e5), internal_energy_j=state.total_internal_energy,
            feed=Amounts({"s0": 1.0, "s1": 0.0}), phases=[Phase("V", kind="vapor")]))
