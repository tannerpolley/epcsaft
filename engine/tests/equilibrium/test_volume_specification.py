"""Total-volume declaration and exact target action on a neutral homogeneous vapor."""
import pytest
from epcsaft import Mixture, Parameters, _core
from epcsaft.equilibrium import (Amounts, Free, Phase, Problem, SolvedStateActionDirection, SolvedStateActionRequest,
    SolvedStateObservable, SolvedStateObservableKind, compile_problem, solve_equilibrium, solved_state_actions)
from tests.packets import parameters_json


def test_total_volume_and_pressure_action() -> None:
    model = Mixture(Parameters.from_json(parameters_json("gross-2001-methane-ethane", 3), components=("methane", "ethane")))
    problem = Problem(T=220.0, P=Free(), volume_m3=0.002, feed=Amounts({"methane": 0.5, "ethane": 0.5}),
                      phases=[Phase("V", kind="vapor")], row_tolerance=1e-11)
    result = solve_equilibrium(model, problem)
    # One mole in 0.002 m3 is the state at rho = 500 mol/m3; the native check owns the reference pressure.
    state = model.state(220.0, rho=500.0, x=[0.5, 0.5])
    assert result.pressure == pytest.approx(state.pressure, rel=1e-8), result.message
    compiled = compile_problem(model, problem)
    volume_row = next(row for row in compiled.rows if row.eqid == "extensive_volume_specification")
    assert volume_row.kind == _core.RowKind.Volume
    direction = SolvedStateActionDirection(0.0, 0.0, [0.0, 0.0], [], volume_m3=0.002)
    batch = solved_state_actions(compiled, result,
        SolvedStateActionRequest([SolvedStateObservable(SolvedStateObservableKind.Pressure)], [direction]))
    assert batch.results[0].action == pytest.approx(-500 * state.pressure_density_derivative, rel=1e-10)
