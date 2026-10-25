"""Installed-wheel PH/PS transport on a one-mole declared neutral phase."""
import pytest
from typing import Any
from epcsaft.equilibrium import Amounts, Free, Phase, Problem, solve_equilibrium
from epcsaft.equilibrium import SolvedStateActionDirection, SolvedStateActionRequest, SolvedStateObservable, SolvedStateObservableKind, compile_problem, solved_state_actions
from tests.native_inputs import A_RECORD, twin_ethane

@pytest.mark.parametrize("field,quantity", [("enthalpy_j", "total_enthalpy"), ("entropy_j_per_k", "total_entropy")])
def test_caloric_transport(field: str, quantity: str) -> None:
    model = twin_ethane(A_RECORD, A_RECORD)
    # At 300 K and 1 MPa the sole root, 434.83 mol/m3, lies below the 6,337 mol/m3 negative-slope sample: vapor (#423).
    state = model.state(300.0, P=1.0e6, x=[1.0, 0.0], phase="vapor")
    target = getattr(state, quantity)
    problem = Problem(T=Free(280.0), P=1.0e6, feed=Amounts({"s0": 1.0, "s1": 0.0}), phases=[Phase("V", kind="vapor")], row_tolerance=1e-11, **{field: target})
    result = solve_equilibrium(model, problem)
    assert result.temperature == pytest.approx(300.0, rel=1e-8), result.message
    controls: dict[str, Any] = {field: 1.0}
    direction = SolvedStateActionDirection(0.0, 0.0, [0.0, 0.0], [], **controls)
    batch = solved_state_actions(compile_problem(model, problem), result, SolvedStateActionRequest([SolvedStateObservable(SolvedStateObservableKind.Temperature)], [direction]))
    # One mole at fixed P: dT/dH = 1/c_p and dT/dS = T/c_p.
    assert state.total_isobaric_heat_capacity is not None
    expected = (1.0 if field == "enthalpy_j" else 300.0) / state.total_isobaric_heat_capacity
    assert batch.results[0].action == pytest.approx(expected, rel=1e-10)
