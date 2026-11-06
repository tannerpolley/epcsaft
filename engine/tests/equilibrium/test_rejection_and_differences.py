"""#976: rejected declared splits carry one verdict; native central differences verify exact derivatives."""
import math
from functools import cache

import pytest

from epcsaft import (
    DensityActionDirection, Mixture, ParameterError, Parameters, PressureActionDirection, PropertyObservable, Pure, Stop,
    ThermochemistryPacket, _core,
)
from epcsaft._core import EquilibriumResult, EquilibriumStart, EquilibriumWorkLimits, StabilityReport
from epcsaft.equilibrium import (
    Amounts, Phase, Problem, SolvedStateActionDirection, SolvedStateActionRequest, SolvedStateObservable, SolvedStateObservableKind, bubble_point,
    compile_problem, equilibrium_numerical_estimates, jacobian_table, solve_equilibrium, solve_equilibrium_starts,
    solved_state_difference,
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
    recovered state includes the equilibrated tangent. The search work and one additional eighth-step evaluation in the
    subsequently refused Newton attempt are added."""
    problem = Problem(phases=[Phase("L", kind="liquid"), Phase("V", kind="vapor")], T=240.0, P=4535437.839536398,
                      feed=Amounts({"methane": 0.7, "ethane": 0.3}))
    result = solve_equilibrium(methane_ethane(), problem)
    report = result.recovery_stability
    assert report is not None and result.success and any(t.resolved_descent for t in report.trials)
    assert list(result.amounts) == [0.08578566729375332, 0.9142143327062465]  # equilibrated tangent; 71 / 144 / 92 recovery work
    assert work(result) == tuple(a + b for a, b in zip((71, 145, 92), charged(report)))


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


def test_polish_trial_beyond_close_packing_is_a_domain_refusal() -> None:
    """#1000: N-(2-aminoethyl)ethanolamine (Esper 2023) saturation at 150 K; the polish's damped trial at 1.4e162
    mol/m³ reached the site solve. No retained warning may carry a density at or above close packing."""
    if not (PACKETS / "feos-esper-2023/1/packet.toml").is_file():
        pytest.skip("CC BY-NC Esper packet absent; restore it with bash data/tools/restore_feos.sh")
    model = Pure(Parameters.from_packets(PACKETS, [("feos-esper-2023/1", ["n-(2-aminoethyl)ethanolamine"])]))
    result = solve_equilibrium(model, model.saturation(150.0))
    close_packed = model.molar_density_upper_bound(150.0)
    # Site-solve warnings report (T, molar density, x).
    densities = [r.first_conditioning_warning.controls[1] for r in result.rank_records
                 if r.basis == "mass_action_site_fractions" and r.first_conditioning_warning is not None]
    assert all(0.0 < d < close_packed for d in densities), densities


def test_dilute_vapor_estimate_resolves_with_pressure_in_pascal() -> None:
    """#1007: at N-(2-aminoethyl)ethanolamine's (Esper 2023) 150 K saturation, P = 1.4e-16 Pa, the estimate once
    called the nonsingular rows rank 2. #1008: with the pressure scale h_P = P its pressure radius is 4.2e-13 P, not
    the 11 P of a 1 Pa floor. At 300 K a zero or negative pressure coordinate refuses before it becomes a scale."""
    if not (PACKETS / "feos-esper-2023/1/packet.toml").is_file():
        pytest.skip("CC BY-NC Esper packet absent; restore it with bash data/tools/restore_feos.sh")
    model = Pure(Parameters.from_packets(PACKETS, [("feos-esper-2023/1", ["n-(2-aminoethyl)ethanolamine"])]))
    result = solve_equilibrium(model, model.saturation(150.0))
    assert result.success, result.message
    estimate = equilibrium_numerical_estimates(model, result)
    assert estimate.cause == "" and estimate.rank == estimate.dimension
    assert all(e.available and e.resolved and math.isfinite(e.absolute) for e in (estimate.pressure, *estimate.densities))
    assert result.pressure is not None and estimate.pressure.absolute <= 1e-11 * result.pressure
    compiled = compile_problem(model, model.saturation(150.0))
    tangent = _core.sensitivity(compiled, result.coordinates, [1.0, 0.0, 0.0])
    pressure_index = compiled.coordinate_labels().index("P / Pa")
    request = SolvedStateActionRequest([SolvedStateObservable(SolvedStateObservableKind.Pressure, 0, 0)],
                                       [SolvedStateActionDirection(1.0, 0.0, [0.0], [])])
    action = _core.solved_state_actions(compiled, result, request)
    assert action.diagnostics.rank == estimate.rank
    derivative = action.results[0].action
    assert derivative is not None
    assert tangent[pressure_index] == pytest.approx(derivative, rel=1e-10, abs=0.0)
    [check] = solved_state_difference(compiled, result, request.observables, request.directions[0], 0.001, 1e-6, 10.0)
    assert check.status == "passed"
    problem = model.saturation(300.0)
    state, compiled = solve_equilibrium(model, problem), compile_problem(model, problem)
    k = compiled.coordinate_labels().index("P / Pa")
    for P in (0.0, -1.0):
        x = [P if i == k else q for i, q in enumerate(state.coordinates)]
        assert compiled.numerical_estimates(x, 300.0, P, state.amounts, state.mole_fractions, state.molar_densities,
                                            state.packing_fractions).cause == "nonpositive_or_nonfinite_pressure"


def test_ideal_action_beyond_close_packing_refuses() -> None:
    """#1000: an ideal-only action at twice toluene's close-packing density (HDA PC-SAFT, NASA9 at 300 K) refuses
    like its residual action instead of returning an ideal entropy."""
    model = Pure(Parameters.from_json(PACKETS / "hda-pc-saft/1/parameters/parameters.json", components=["toluene"]),
                 thermochemistry=ThermochemistryPacket(
                     PACKETS / "nasa-glenn-2002-hda-species/1/parameters/nasa9-records.json").select(["toluene"]))
    with pytest.raises(ParameterError, match="outside the fluid domain"):
        model.density_actions(300.0, 2.0 * model.molar_density_upper_bound(300.0),
                              observables=(PropertyObservable.IdealEntropy,), directions=(DensityActionDirection(0.0, 1.0),))


@pytest.mark.parametrize(("T", "packing"), [(150.0, None), (220.0, None), (300.0, 1e-200)])
def test_dilute_vapor_saturation_and_unrepresentable_refusal(T: float, packing: float | None) -> None:
    """#1004: N-(2-aminoethyl)ethanolamine (Esper 2023) saturation with vapor packing below 1e-12 meets
    P_sat = P0 phi_L(T, P0), P0 = 1e-3 Pa, within 1e-6 (the dropped Poynting and vapor terms are below 3e-10).
    A supplied vapor packing of 1e-200 at 300 K is not representable and must return no state."""
    if not (PACKETS / "feos-esper-2023/1/packet.toml").is_file():
        pytest.skip("CC BY-NC Esper packet absent; restore it with bash data/tools/restore_feos.sh")
    model = Pure(Parameters.from_packets(PACKETS, [("feos-esper-2023/1", ["n-(2-aminoethyl)ethanolamine"])]))
    if packing is None:
        result = solve_equilibrium(model, model.saturation(T))
        assert result.success, result.message
        log_phi = model.state(T, P=1e-3, phase="liquid").log_fugacity_coefficient
        assert log_phi is not None and log_phi[0] is not None
        reference = 1e-3 * math.exp(log_phi[0])
        assert result.pressure is not None and abs(result.pressure / reference - 1.0) <= 1e-6
        return
    out = solve_equilibrium_starts(model, Problem(phases=[Phase("V", kind="vapor")], T=T, P=1e-192,
                                                  feed=Amounts({model.component_ids[0]: 1.0})),
                                   [EquilibriumStart(compositions=[[1.0]], packing_fractions=[packing])],
                                   EquilibriumWorkLimits(attempts=1))
    [refused] = out.results
    assert out.selected is None and refused is not None and not refused.success
    assert refused.pressure is None and not list(refused.molar_densities) and not list(refused.packing_fractions)
