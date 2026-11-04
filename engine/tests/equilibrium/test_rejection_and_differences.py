"""#976: rejected declared splits carry one verdict; native central differences verify exact derivatives."""
from functools import cache

import pytest

from epcsaft import Mixture, Parameters, PressureActionDirection, PropertyObservable, Pure, Stop
from epcsaft._core import EquilibriumResult, StabilityReport
from epcsaft.equilibrium import (
    Amounts, Phase, Problem, SolvedStateActionDirection, SolvedStateObservable, SolvedStateObservableKind, bubble_point,
    compile_problem, jacobian_table, solve_equilibrium, solved_state_difference,
)

from tests.packets import PACKETS

UNSUPPORTED_CLAIMS = ("stable phase", "does not exist", "exceeds", "no bubble point exists")


@cache
def methane_ethane() -> Mixture:
    return Mixture(Parameters.from_packets(PACKETS, [("gross-2001-methane-ethane/3", ["methane", "ethane"])]))


def flash(T: float, P: float) -> Problem:
    return Problem(phases=[Phase("L", kind="liquid"), Phase("V", kind="vapor")], T=T, P=P,
                   feed=Amounts({"methane": 0.3, "ethane": 0.7}))


def charged(report: StabilityReport) -> tuple[float, ...]:
    """Driver plus native tangent-plane work the #987 trial search was charged."""
    work = report.work
    return tuple(work.get(driver, 0.0) + work.get(native, 0.0) for driver, native in (
        ("iterations", "trial_tangent_plane_iterations"), ("evaluations", "trial_tangent_plane_evaluations"),
        ("hessian_calls", "trial_tangent_plane_hessian_evaluations")))


def work(result: EquilibriumResult) -> tuple[int, int, int]:
    return result.iterations, result.evaluations, result.hessian_calls


def test_single_phase_flash_verdict() -> None:
    """F1: the 1 bar vapor feed's declared L+V split coalesces; the verdict says so and claims nothing more. The finite
    trial search finds no descent, so recovery (512 iterations at A) is not run: work is the first run plus the search."""
    problem = flash(250.0, 1e5)
    result = solve_equilibrium(methane_ethane(), problem)
    assert not result.success and result.solver_status == "SUCCESS"
    assert result.convergence.stop == Stop.infeasible and result.convergence.cause == "phases_coalesced"
    assert result.message.startswith("declared phases L and V coalesced at T = 250 K, P = 100000 Pa (phase distance ")
    assert not any(claim in result.message for claim in UNSUPPORTED_CLAIMS)
    assert ("tangent-plane search before recovery from single-phase V1, 7 trials: no witness in the recorded trial set; "
            "unresolved trials 0; not a global stability certificate: a split outside this finite trial set is not "
            "excluded; recovery not run") in result.message and "; recovery: " not in result.message
    report = result.recovery_stability
    assert report is not None
    assert report.completed and len(report.trials) == 7
    assert all(t.evidence == "no_resolved_descent" and t.optimizer.success for t in report.trials)
    first = compile_problem(methane_ethane(), problem).solve()
    assert work(result) == tuple(a + b for a, b in zip(work(first), charged(report)))


def test_recoverable_flash_keeps_recovery() -> None:
    """F7: at 240 K and 4.535 MPa the search finds a resolved descending split, so the unchanged recovery runs; the
    recovered state is the one recovery reached before #987 and only the search's work is added."""
    problem = Problem(phases=[Phase("L", kind="liquid"), Phase("V", kind="vapor")], T=240.0, P=4535437.839536398,
                      feed=Amounts({"methane": 0.7, "ethane": 0.3}))
    result = solve_equilibrium(methane_ethane(), problem)
    report = result.recovery_stability
    assert report is not None and result.success and any(t.resolved_descent for t in report.trials)
    assert list(result.amounts) == [0.08578566729375602, 0.914214332706244]  # before #987, 71 / 144 / 92 work
    assert work(result) == tuple(a + b for a, b in zip((71, 144, 92), charged(report)))


def test_bubble_point_above_coalescence_verdict() -> None:
    """F3: a coalesced bubble point names what was observed and leaves the critical limit untested."""
    problem = bubble_point(liquid=Phase("L", kind="liquid"), T=250.0, feed=Amounts({"methane": 0.9, "ethane": 0.1}))
    result = solve_equilibrium(methane_ethane(), problem)
    assert not result.success and result.convergence.stop == Stop.infeasible
    assert result.message.startswith("no distinct bubble point was found at T = 250 K: declared phases L and V coalesced")
    assert "one possible, untested explanation" in result.message
    assert not any(claim in result.message for claim in UNSUPPORTED_CLAIMS)


def test_jacobian_difference_and_labels() -> None:
    """F5/F10: every residual-Jacobian entry at a solved flash, including exact zeros, matches central differences."""
    problem = flash(250.0, 3.75e6)
    compiled, result = compile_problem(methane_ethane(), problem), solve_equilibrium(methane_ethane(), problem)
    checks = compiled.jacobian_difference(result.coordinates, 1e-6, 1e-6, 10.0)
    assert {check.status for check in checks} == {"passed"}
    assert any(check.exact == 0.0 and check.ratio == 0.0 for check in checks)
    table = jacobian_table(compiled, result.coordinates)
    assert table["columns"] == ["ln x[L,component 0]", "ln x[L,component 1]", "n[L] / mol", "ln eta[L]",
                                "ln x[V,component 0]", "ln x[V,component 1]", "n[V] / mol", "ln eta[V]",
                                "lambda[moiety 0]", "lambda[moiety 1]"]
    assert [len(row) for row in table["values"]] == [len(result.coordinates)] * len(table["rows"])


def test_directional_differences() -> None:
    """F8/F9: dH_res/dT of liquid water and dP_bubble/dT of methane/ethane against native central differences."""
    water = Pure(Parameters.from_packets(PACKETS, [("held-cameretti-sadowski-2008/4", ["water"])]))
    property_checks = water.pressure_action_difference(
        298.15, 1e5, phase="liquid", observables=[PropertyObservable.ResidualEnthalpy],
        direction=PressureActionDirection(temperature=1.0, composition=[0.0]), step=1e-3, tau=1e-8, kappa=10.0)
    problem = bubble_point(liquid=Phase("L", kind="liquid"), T=230.0, feed=Amounts({"methane": 0.3, "ethane": 0.7}))
    pressure = SolvedStateObservable(SolvedStateObservableKind.Pressure, 0, 0)
    state_checks = solved_state_difference(compile_problem(methane_ethane(), problem), solve_equilibrium(methane_ethane(), problem),
                                           [pressure], SolvedStateActionDirection(1.0, 0.0, [0.0, 0.0], []), 1e-3, 1e-8, 10.0)
    assert [check.status for check in property_checks + state_checks] == ["passed", "passed"]
    with pytest.raises(ValueError):
        water.pressure_action_difference(298.15, 1e5, phase="liquid", observables=[PropertyObservable.ResidualEnthalpy],
                                         direction=PressureActionDirection(temperature=1.0, composition=[0.0]),
                                         step=0.0, tau=1e-8, kappa=10.0)


def test_failed_neighbor_is_a_status() -> None:
    """A 70 K step from the 230 K bubble point fails its 300 K neighbor; the check reports it instead of crashing."""
    problem = bubble_point(liquid=Phase("L", kind="liquid"), T=230.0, feed=Amounts({"methane": 0.3, "ethane": 0.7}))
    pressure = SolvedStateObservable(SolvedStateObservableKind.Pressure, 0, 0)
    [check] = solved_state_difference(compile_problem(methane_ethane(), problem), solve_equilibrium(methane_ethane(), problem),
                                      [pressure], SolvedStateActionDirection(1.0, 0.0, [0.0, 0.0], []), 70.0, 1e-8, 10.0)
    assert check.status == "neighbor_unavailable" and check.cause
