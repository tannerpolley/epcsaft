"""The ASL equilibrium record on the engine/README.md water/isobutanol record (298.15 K, 1e5 Pa, feed (0.7, 0.3) mol) with
k_ij(isobutanol, water) and isobutanol eps/k active, q = (T, P, u_water, u_isobutanol, k_ij, eps/k). F1: values and exact actions
equal the direct Engine (bitwise at theta0 and for a record with no active parameters; 1e-12 relative at theta*). F2: the
columns match central differences (verification only). F3: two temperatures recover theta* through Ipopt. F4: where the
liquids merge the callback names the cause, and the next calls recover. M1-M3 repeat F1-F3 on a declared vapor-liquid
pair: the packet-3 MEA/CO2/water flash (charged, reacting liquid; neutral vapor support) at 318.15 K and 5918.5286 Pa.
R1-R3 move T and P on the reference MEA declaration, whose ln K(T) laws, reaction references, neutral reference and model
thermochemistry the record carries. E1-E4 check its system enthalpy callback: the direct Engine, central differences, an
adiabatic flash through Ipopt, and the record without thermochemistry."""
from __future__ import annotations

import dataclasses
import functools
import json
import shutil
import subprocess
from pathlib import Path
from typing import Any

import numpy as np
import pytest
from pyomo.environ import ConcreteModel, Constraint, ExternalFunction, Var, check_optimal_termination, value  # type: ignore[import-untyped]

from epcsaft import Mixture, Parameters, _core, regression
from epcsaft import equilibrium as eq
from tests.equilibrium.test_declared_starts import _mea_flash
from tests.equilibrium.test_trace_reliability import _MEA, _MEA_REACTIONS, _mea_mapping, _reference_mea
from tests.packets import parameters_json

NAMES, G0, G1, TP0 = ("water", "isobutanol"), (0.7, 0.3), (0.7, 0.285), (298.15, 1.0e5)
THETA0, THETA_STAR, THETA_MERGED = (-0.0071, 257.13), (-0.0021, 258.13), (-0.0371, 257.13)
OPTIONS = eq.Continuation(initial_step=0.05, minimum_step=1e-3, max_steps=200, max_retries_per_step=5)
FUNCTION, ENTHALPY, SOLID = "epcsaft_equilibrium_phase_component_amount", "epcsaft_equilibrium_system_enthalpy", "epcsaft_equilibrium_solid_amount"


def mixture(theta: Any = None, active: bool = True) -> Mixture:
    mapping = json.loads(parameters_json("khudaida-2026-figure-2-electrolyte-lle", 7).read_text())
    if theta is not None:  # the edited packet: k_ij and isobutanol eps/k at theta
        pair = next(p for p in mapping["pairs"] if p["coefficients"][0]["identity"] == "pair/isobutanol/water/k_ij")
        pair["coefficients"][0]["value"]["magnitude"] = theta[0]
        iso = next(c for c in mapping["components"] if c["component_id"] == "isobutanol")
        next(k for k in iso["coefficients"] if k["identity"].endswith("dispersion_energy_over_k"))["value"]["magnitude"] = theta[1]
    parameters = Parameters.from_mapping(mapping, components=NAMES)
    return Mixture(parameters, active_parameters=[regression.active_parameter(parameters, "k_ij", ("isobutanol", "water")),
                   regression.active_parameter(parameters, "dispersion_energy_over_k", "isobutanol")] if active else None)


def problem(feed: Any, T: float = 298.15) -> eq.Problem:
    return eq.Problem(phases=[eq.Phase("L1", kind="liquid", support=NAMES, composition_guess=(0.46, 0.54)),
                              eq.Phase("L2", kind="liquid", support=NAMES, composition_guess=(0.99, 0.01))],
                      T=T, P=1.0e5, feed=eq.Amounts(dict(zip(NAMES, feed, strict=True))))


def record(model: Mixture, T: float = 298.15) -> str:
    return eq.asl_record(model, problem(G0, T), OPTIONS, [0.0, 0.0], [[1.0, 0.0], [0.0, 1.0]], list(G0), 1, 0.1)


def fgh(text: str, q: Any, row: int | None, order: int, components: int = 2, phases: int = 2, fixed: Any = None) -> tuple[float, Any, Any]:
    """Row p*components + c at q, solid s at row phases*components + s, or the system enthalpy for row None; the Hessian
    unpacked from the lower triangle i + j(j+1)/2. `fixed` marks the real arguments ASL passes as constants (dig)."""
    from epcsaft.idaes import PropertyExternalFunction, asl_library
    solid = row is not None and row >= phases * components
    callback = PropertyExternalFunction(library=str(asl_library()), function=ENTHALPY if row is None else SOLID if solid else FUNCTION)
    strings = [] if row is None else [str(row - phases * components)] if solid else [str(row // components), str(row % components)]
    marks = None if fixed is None else [*fixed, *[False] * (1 + len(strings))]
    v, g, h = callback.evaluate_fgh([*q, text, *strings], fixed=marks, fgh=order)
    n = len(q)
    hessian = np.array([[h[min(i, j) + max(i, j) * (max(i, j) + 1) // 2] for j in range(n)] for i in range(n)]) if order == 2 else None
    return float(v), np.array(g[:n]) if order else None, hessian


def native(model: Mixture, feed: Any, target: eq.Problem | None = None, enthalpy: bool = False, matrix: Any = None) -> list[list[float]]:
    """Per row: value, first actions, second actions (i <= j, packed order) over (T, P, feed, theta), from solve_equilibrium
    and solved_state_actions; the feed controls move the components by `matrix`, else the first len(feed) components. The
    rows are the phase-component amounts then the solid amounts, or the system enthalpy alone."""
    target = target or problem(feed)
    central = eq.solve_equilibrium(model, target)
    compiled = _core.compile_problem(model._model, target._native(model))
    F, S = central.phase_count, len(central.solid_amounts)
    k, n, m = len(feed), len(central.mole_fractions) // F, len(model._model.active_parameter_values)
    e = np.eye(2 + k + m).tolist()
    d = [eq.SolvedStateActionDirection(c[0], c[1], (np.array(matrix) @ c[2:2 + k]).tolist() if matrix else c[2:2 + k] + [0.0] * (n - k),
                                       c[2 + k:]) for c in e]
    observables = [eq.SolvedStateObservable(eq.SolvedStateObservableKind.SystemEnthalpy)] if enthalpy else [
        eq.SolvedStateObservable(eq.SolvedStateObservableKind.PhaseComponentAmount, p, c) for p in range(F) for c in range(n)] + [
        eq.SolvedStateObservable(eq.SolvedStateObservableKind.SolidAmount, s) for s in range(S)]

    def batch(*directions: Any) -> Any:
        return _core.solved_state_actions(compiled, central, eq.SolvedStateActionRequest(observables, list(directions))).results

    zero = batch(eq.SolvedStateActionDirection(0.0, 0.0, [0.0] * n, [0.0] * m))
    first = [batch(d[i]) for i in range(2 + k + m)]
    second = [batch(d[i], d[j]) for j in range(2 + k + m) for i in range(j + 1)]
    return [[zero[r].value, *(b[r].action for b in first), *(b[r].action for b in second)] for r in range(len(observables))]


def callback_row(text: str, q: Any, row: int | None, components: int = 2, phases: int = 2, fixed: Any = None) -> list[float]:
    v, g, h = fgh(text, q, row, 2, components, phases, fixed)
    return [v, *g, *(h[i, j] for j in range(len(q)) for i in range(j + 1))]


def test_record_equals_direct_engine() -> None:
    """F1: at theta0 (G0) and with no active parameters (G0, G1) every row's value, gradient and Hessian equal the direct
    Engine bitwise, and each order returns the same value; at theta* (G0, G1) they match the edited packet to 1e-12."""
    model, inactive = mixture(), mixture(active=False)
    assert model._model.active_parameter_values == list(THETA0)
    cases = [(record(model), [*TP0, *G0, *THETA0], native(model, G0), True)]
    cases += [(record(inactive), [*TP0, *G], native(inactive, G), True) for G in (G0, G1)]
    cases += [(record(model), [*TP0, *G, *THETA_STAR], native(mixture(THETA_STAR), G), False) for G in (G0, G1)]
    for text, q, want, bitwise in cases:
        values = [[fgh(text, q, row, order)[0] for row in range(4)] for order in (0, 1)]  # before any Hessian at q
        got = [callback_row(text, q, row) for row in range(4)]
        assert values == [[r[0] for r in got]] * 2
        for row in range(4):
            if bitwise:
                assert [x.hex() for x in got[row]] == [float(x).hex() for x in want[row]], (q, row)
            else:
                assert got[row] == pytest.approx(want[row], rel=1e-12, abs=0.0), (q, row)


def test_constant_arguments_keep_the_other_columns() -> None:
    """A1: at a state the callback first enters with arguments marked constant, the value and every first and second
    action among the other arguments equal the direct Engine bitwise, and no entry of a constant argument is written."""
    from epcsaft.idaes import PropertyExternalFunction, asl_library
    model, pairs = mixture(active=False), [(i, j) for j in range(4) for i in range(j + 1)]
    for G, fixed in (((0.7, 0.295), [True, True, False, False]), ((0.7, 0.29), [False, False, True, False])):
        callback = PropertyExternalFunction(library=str(asl_library()), function=FUNCTION)
        _, g, h = callback._evaluate([*TP0, *G, record(model), "0", "0"], [*fixed, False, False, False], 2)
        assert {g[i] for i in range(4) if fixed[i]} | {h[i + j * (j + 1) // 2] for i, j in pairs if fixed[i] or fixed[j]} == {0.0}, G
        want = native(model, G)
        keep = [True, *(not f for f in fixed), *(not (fixed[i] or fixed[j]) for i, j in pairs)]
        for row in range(4):
            got = callback_row(record(model), [*TP0, *G], row, fixed=fixed)
            assert [x.hex() for x, k in zip(got, keep, strict=True) if k] == [float(x).hex() for x, k in zip(want[row], keep, strict=True) if k], (G, row)


def test_columns_match_central_differences() -> None:
    """F2, verification only: feed (h 1e-5 mol), k_ij (h 1e-5) and eps/k (h 1e-3 K) at h x {4, 1, 0.25}: every gradient
    entry (tau 1e-6) and every Hessian entry with a parameter (tau 1e-5) meet |D - A| <= tau|A| + 10 eps (|f+| + |f-|)/2h."""
    text, eps, q = record(mixture()), float(np.finfo(float).eps), np.array([*TP0, *G0, *THETA0])
    exact = [fgh(text, q, row, 2) for row in range(4)]
    for b, h0 in ((2, 1e-5), (3, 1e-5), (4, 1e-5), (5, 1e-3)):
        for h in (4 * h0, h0, 0.25 * h0):
            plus, minus = ([fgh(text, q + s * h * np.eye(6)[b], row, 1) for row in range(4)] for s in (1, -1))
            for row, (_, gradient, hessian) in enumerate(exact):
                pairs = [(gradient[b], plus[row][0], minus[row][0], 1e-6)]
                pairs += [(hessian[a, b], plus[row][1][a], minus[row][1][a], 1e-5) for a in range(2, 6) if b >= 4 or a >= 4]
                for exact_value, fp, fm, tau in pairs:
                    assert abs((fp - fm) / (2 * h) - exact_value) <= tau * abs(exact_value) + 10 * eps * (abs(fp) + abs(fm)) / (2 * h), (row, b, h)


def test_parameter_estimation_ipopt_round_trip(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """F3: the L2 isobutanol amounts at theta* (298.15 K and 348.15 K, feed G0) recover theta* from theta0 through Ipopt;
    the external AMPL Ipopt exits with code 0 on the same problem."""
    from idaes.core.solvers import get_solver  # type: ignore[import-untyped]
    from epcsaft.idaes import asl_library, solve_with_diagnostics
    solver, executable = get_solver(), shutil.which("ipopt") or str(Path.home() / ".idaes/bin/ipopt")
    if not solver.available(exception_flag=False) or not Path(executable).is_file():
        pytest.skip("Ipopt executable unavailable")
    library = str(asl_library())
    monkeypatch.setenv("AMPLFUNC", library)
    model, targets = mixture(), (0.002356016921554457, 0.004444355639825189)
    m = ConcreteModel()
    m.theta = Var([0, 1], initialize=dict(enumerate(THETA0)))
    m.amount = ExternalFunction(library=library, function=FUNCTION)
    for k, (T, target) in enumerate(zip((298.15, 348.15), targets, strict=True)):
        m.add_component(f"target_{k}", Constraint(expr=m.amount(T, TP0[1], *G0, m.theta[0], m.theta[1], record(model, T), "1", "1") / target == 1))
    m.write(str(tmp_path / "estimate.nl"), format="nl")
    (tmp_path / "ipopt.opt").write_text("linear_solver mumps\n")
    external = subprocess.run([executable, "estimate.nl", "-AMPL"], cwd=tmp_path, capture_output=True, text=True, check=False)
    assert external.returncode == 0, external.stdout + external.stderr
    assert "Optimal Solution Found" in (tmp_path / "estimate.sol").read_text()
    log = tmp_path / "ipopt.log"
    result = solve_with_diagnostics(solver, m, logfile=log)
    assert "epcsaft_equilibrium failure" not in log.read_text()
    assert check_optimal_termination(result)
    for got, want in zip(m.theta.values(), THETA_STAR, strict=True):
        assert abs(value(got) - want) <= 1e-6 * max(abs(want), 1.0)


def test_merged_liquids_name_the_cause_and_recover() -> None:
    """F4: at k_ij -0.0371 the liquids merge, so the direct solve and the parameter leg fail and the value-only call is NaN
    with the named cause; the following theta* and theta0 calls equal their earlier values bitwise."""
    from epcsaft.idaes import PropertyExternalFunction, asl_library
    model = mixture()
    text = record(model)
    star, start = [callback_row(text, [*TP0, *G0, *THETA_STAR], row) for row in range(4)], native(model, G0)
    callback = PropertyExternalFunction(library=str(asl_library()), function=FUNCTION)
    with pytest.raises(RuntimeError) as failure:
        callback.evaluate([*TP0, *G0, *THETA_MERGED, text, "1", "1"])
    assert str(failure.value).startswith("epcsaft_equilibrium:phase_boundary_not_differentiable:")
    assert "parameter continuation: phase_disappeared" in str(failure.value)
    assert [[x.hex() for x in callback_row(text, [*TP0, *G0, *THETA_STAR], row)] for row in range(4)] == [[x.hex() for x in r] for r in star]
    assert [[x.hex() for x in callback_row(text, [*TP0, *G0, *THETA0], row)] for row in range(4)] == [[float(x).hex() for x in r] for r in start]



# M: the packet-3 MEA flash with the vapor support CO2, MEA and water; feed controls u = (CO2, MEA, water) mol.
# Reaction targets are affinities -R T ln K: K of R3, exp(-948.2312), underflows. U_STAR is CO2 x 1.05.
MEA_T, MEA_P, MEA_LN_K, C = 318.15, 5918.5286, (-605.348, -587.8, -948.2312, -25.1014, -298.6415), len(_MEA)
U0 = (.25557798094241313 * 9.7292821, 9.7292821, 76.9751925)
U_STAR = (2.61091978889402, *U0[1:])


def mea_model() -> Mixture:
    return Mixture(Parameters.from_mapping(_mea_mapping(parameters_json("mea-co2-h2o-nine-species-estimation", 3)), components=_MEA))


def mea_problem(u: Any, vapor_support: Any = _MEA[:3]) -> eq.Problem:
    reactions = [eq.Reaction(nu, affinity_j_per_mol=-_core.GAS_CONSTANT_J_PER_MOL_K * MEA_T * a) for nu, a in zip(_MEA_REACTIONS, MEA_LN_K, strict=True)]
    return dataclasses.replace(_mea_flash(MEA_P), phases=[eq.Phase("L", kind="liquid", support=_MEA), eq.Phase("V", kind="vapor", support=vapor_support)],
                               reactions=reactions, feed=eq.Amounts(dict(zip(_MEA[:3], u, strict=True))), row_tolerance=1e-12)


def mea_record(model: Mixture, vapor_support: Any = _MEA[:3]) -> str:
    matrix = np.eye(C, 3).tolist()  # the 9 x 3 identity on CO2, MEA and water
    options = eq.Continuation(initial_step=0.25, minimum_step=0.01, max_steps=40, max_retries_per_step=2)
    return eq.asl_record(model, mea_problem(U0, vapor_support), options, [0.0] * C, matrix, list(U0), 1, 0.03)


def test_vapor_liquid_record_equals_direct_engine() -> None:
    """M1: at u0 and u* all 18 values equal amounts[p] x[p][c] of the direct solve to 1e-12 relative (the six vapor ions
    exactly 0), and each gradient and packed Hessian equals solved_state_actions on the same state to 1e-12 relative."""
    model = mea_model()
    text = mea_record(model)
    assert "6:liquid," in text and "5:vapor," in text
    for u in (U0, U_STAR):
        direct = eq.solve_equilibrium(model, mea_problem(u))
        assert direct.success and direct.phase_count == 2
        amounts = [direct.amounts[p] * direct.mole_fractions[p * C + c] for p in range(2) for c in range(C)]
        want = native(model, u, mea_problem(u))
        for row in range(2 * C):
            got = callback_row(text, [MEA_T, MEA_P, *u], row, C)
            assert got[0] == pytest.approx(amounts[row], rel=1e-12, abs=0.0), (u, row)
            assert got[1:] == pytest.approx(want[row][1:], rel=1e-12, abs=0.0), (u, row)
        assert [callback_row(text, [MEA_T, MEA_P, *u], C + c, C)[0] for c in range(3, C)] == [0.0] * 6


def test_vapor_liquid_columns_match_central_differences() -> None:
    """M2, verification only: each feed control at h = 1e-4 u0 x {1, 1/4, 1/16}; every gradient entry meets
    |D - A| <= 1e-6 |A| + 10 max(eps, row_tolerance) (|f+| + |f-|)/2h, since a re-solve reproduces amounts to its row tolerance."""
    text, noise, q = mea_record(mea_model()), 10 * max(float(np.finfo(float).eps), 1e-12), np.array([MEA_T, MEA_P, *U0])
    exact = [fgh(text, q, row, 1, C)[1] for row in range(2 * C)]
    for b in range(2, 5):
        for h in (1e-4 * U0[b - 2] * s for s in (1.0, 0.25, 0.0625)):
            plus, minus = ([fgh(text, q + s * h * np.eye(5)[b], row, 0, C)[0] for row in range(2 * C)] for s in (1, -1))
            for row, gradient in enumerate(exact):
                fp, fm = plus[row], minus[row]
                assert abs((fp - fm) / (2 * h) - gradient[b]) <= 1e-6 * abs(gradient[b]) + noise * (abs(fp) + abs(fm)) / (2 * h), (row, b, h)


def test_vapor_liquid_ipopt_round_trip(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """M3: with u_MEA and u_water fixed, Ipopt frees u_CO2 in (2.3, 2.9) from u0 and recovers u*_CO2 to 1e-6 relative
    from the liquid MEAH+ amount at u*; the external AMPL Ipopt exits with code 0 and reports no callback failure."""
    from idaes.core.solvers import get_solver
    from epcsaft.idaes import asl_library, solve_with_diagnostics
    solver, executable = get_solver(), shutil.which("ipopt") or str(Path.home() / ".idaes/bin/ipopt")
    if not solver.available(exception_flag=False) or not Path(executable).is_file():
        pytest.skip("Ipopt executable unavailable")
    library = str(asl_library())
    monkeypatch.setenv("AMPLFUNC", library)
    m = ConcreteModel()
    m.u = Var(initialize=U0[0], bounds=(2.3, 2.9))
    m.amount = ExternalFunction(library=library, function=FUNCTION)
    m.target = Constraint(expr=m.amount(MEA_T, MEA_P, m.u, *U0[1:], mea_record(mea_model()), "0", "3") == 6.302671442747203)
    m.write(str(tmp_path / "roundtrip.nl"), format="nl")
    (tmp_path / "ipopt.opt").write_text("linear_solver mumps\n")
    external = subprocess.run([executable, "roundtrip.nl", "-AMPL"], cwd=tmp_path, capture_output=True, text=True, check=False)
    assert external.returncode == 0 and "epcsaft_equilibrium failure" not in external.stderr, external.stdout + external.stderr
    assert "Optimal Solution Found" in (tmp_path / "roundtrip.sol").read_text()
    log = tmp_path / "ipopt.log"
    assert check_optimal_termination(solve_with_diagnostics(solver, m, logfile=log))
    assert "epcsaft_equilibrium failure" not in log.read_text()
    assert abs(value(m.u) - U_STAR[0]) <= 1e-6 * U_STAR[0]


def test_vapor_support_refusals() -> None:
    """A vapor support with a charged component is refused when the record is built and when it is read; two vapors and
    the previous record version are refused."""
    from epcsaft.idaes import PropertyExternalFunction, asl_library
    model, contract = mea_model(), "invalid ASL affine feed/phase contract"
    with pytest.raises(_core.ParameterError, match=contract):
        mea_record(model, _MEA[:4])
    text = mea_record(model)
    vapor = "1:V,5:vapor,1:3,1:0,1:1,1:2,"
    assert vapor in text
    callback = PropertyExternalFunction(library=str(asl_library()), function=FUNCTION)
    for edited, cause in ((text.replace(vapor, "1:V,5:vapor,1:4,1:0,1:1,1:2,1:3,"), contract),
                          (text.replace("1:L,6:liquid,", "1:L,5:vapor,"), contract),
                          (text.replace("equilibrium.v8", "equilibrium.v7"), "")):
        with pytest.raises(RuntimeError, match="epcsaft_equilibrium:.*" + cause):
            callback.evaluate([MEA_T, MEA_P, *U0, edited, "0", "3"])


# R: the reference MEA declaration (engine/data/reference/mea-nominally-unloaded, with its ideal calorics unless stated)
# with a declared liquid + vapor pair at 353.15 K and 45 kPa; q = (T, P, u_CO2, u_MEA, u_water), the loading of the
# Hilliard (2008) Cp row. Both phases guess the anchor's own solved compositions.
R_Q0 = (353.15, 45000.0, 0.358, 1.0, 7.929776437398854)


@functools.cache
def reference_mea(thermochemistry: bool = True) -> tuple[Mixture, Any]:
    """The model and target(q), the declared problem at q."""
    model, reactions, reference, guess, _ = _reference_mea(thermochemistry)

    def target(q: Any, liquid: Any = guess, vapor: Any = (0.01, 0.002, 0.988)) -> eq.Problem:
        return eq.Problem(T=q[0], P=q[1], row_tolerance=1e-10, topology_declared=True, reactions=reactions, neutral_reference=reference,
                          feed=eq.Amounts(dict(zip(_MEA, (*q[2:], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0), strict=True))),
                          phases=[eq.Phase("L", kind="liquid", support=_MEA, amount=eq.Free(), composition_guess=list(liquid)),
                                  eq.Phase("V", kind="vapor", support=_MEA[:3], amount=eq.Free(), composition_guess=list(vapor))])
    x = eq.solve_equilibrium(model, target(R_Q0)).mole_fractions
    return model, lambda q: target(q, x[:C], x[C:C + 3])


def reference_record(thermochemistry: bool = True) -> str:
    model, target = reference_mea(thermochemistry)
    options = eq.Continuation(initial_step=0.25, minimum_step=1e-3, max_steps=64, max_retries_per_step=8)
    return eq.asl_record(model, target(R_Q0), options, [0.0] * C, np.eye(C, 3).tolist(), list(R_Q0[2:]), 1, 0.01)


def test_reference_continuation_fallback() -> None:
    model, target = reference_mea()
    q = (353.15, 41000.0, *R_Q0[2:])
    assert not eq.solve_equilibrium(model, target(q)).success
    options = eq.Continuation(initial_step=0.25, minimum_step=1e-3, max_steps=64, max_retries_per_step=8)
    path = eq.continue_equilibrium(model, target(R_Q0), target(q), options, eq.solve_equilibrium(model, target(R_Q0)))
    assert path.success, path.message
    x = path.states[-1].mole_fractions
    seeded = target(q)
    seeded.phases[0].composition_guess = list(x[:C])
    seeded.phases[1].composition_guess = list(x[C:C + 3])
    direct = eq.solve_equilibrium(model, seeded)
    assert direct.success, direct.message
    want = native(model, q[2:], seeded)
    text = reference_record()
    for row in range(2 * C):
        assert callback_row(text, q, row, C) == pytest.approx(want[row], rel=1e-8, abs=1e-10), row


def test_reference_record_equals_direct_engine() -> None:
    """R1: at the anchor every row's value, gradient and packed Hessian over q equal solve_equilibrium plus
    solved_state_actions on the declared problem bitwise; at (353.65 K, 45500 Pa), and there with u_CO2 0.36 mol, they
    equal the direct solve of the declared target to 1e-12 relative."""
    model, target = reference_mea()
    text = reference_record()
    for q, bitwise in ((R_Q0, True), ((353.65, 45500.0, *R_Q0[2:]), False), ((353.65, 45500.0, 0.36, *R_Q0[3:]), False)):
        want = native(model, q[2:], target(q))
        for row in range(2 * C):
            got = callback_row(text, q, row, C)
            if bitwise:
                assert [x.hex() for x in got] == [float(x).hex() for x in want[row]], (q, row)
            else:
                assert got == pytest.approx(want[row], rel=1e-12, abs=0.0), (q, row)


def central_differences(rows: list[int | None], feed_step: float) -> None:
    """Verification only: T at {1e-3, 2.5e-4, 6.25e-5} K, P at {0.45, 0.1125, 0.028125} Pa and feed j at feed_step u_j
    x {1, 1/4, 1/16} mol. Each value difference meets |D - A| <= 1e-6 |A| + 10 eps (|f+| + |f-|)/2h against its gradient
    entry, and each first-action difference meets it with tau 1e-5 against its Hessian entry, so every second action,
    including T-T, T-P and P-P, is checked from both of its columns."""
    text, eps, q = reference_record(), float(np.finfo(float).eps), np.array(R_Q0)
    exact = [fgh(text, q, row, 2, C) for row in rows]
    for b, h0 in enumerate((1e-3, 0.45, *(feed_step * u for u in R_Q0[2:]))):
        for h in (h0, h0 / 4, h0 / 16):
            plus, minus = ([fgh(text, q + s * h * np.eye(5)[b], row, 1, C) for row in rows] for s in (1, -1))
            for k, (_, gradient, hessian) in enumerate(exact):
                pairs = [(gradient[b], plus[k][0], minus[k][0], 1e-6)]
                pairs += [(hessian[a, b], plus[k][1][a], minus[k][1][a], 1e-5) for a in range(5)]
                for exact_value, fp, fm, tau in pairs:
                    assert abs((fp - fm) / (2 * h) - exact_value) <= tau * abs(exact_value) + 10 * eps * (abs(fp) + abs(fm)) / (2 * h), (rows[k], b, h)


def test_reference_columns_match_central_differences() -> None:
    """R2: every amount row with feed steps 1e-5 u_j x {16, 4, 1} mol; at 2.5e-6 u_water the re-solve noise in the P
    action (1e-13 relative) exceeds the 10 eps allowance on the small P-u_water entries."""
    central_differences(list(range(2 * C)), 1.6e-4)

def test_reference_temperature_pressure_ipopt_round_trip(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """R3: with u fixed at the anchor, Ipopt (tol 1e-10, "Optimal Solution Found" only) frees T in [353.15, 354.15] K and
    P in [44000, 46000] Pa from the anchor and recovers (353.65 K, 45500 Pa) from the vapor CO2 and liquid carbonate amounts
    there, each relative to its target, within |J*^-1| (tol + 1e-12) rounded up: 6e-9 K and 1.1e-5 Pa. The external AMPL
    Ipopt exits with code 0 and reports no callback failure."""
    from idaes.core.solvers import get_solver
    from epcsaft.idaes import asl_library, solve_with_diagnostics
    solver, executable = get_solver(), shutil.which("ipopt") or str(Path.home() / ".idaes/bin/ipopt")
    if not solver.available(exception_flag=False) or not Path(executable).is_file():
        pytest.skip("Ipopt executable unavailable")
    library = str(asl_library())
    monkeypatch.setenv("AMPLFUNC", library)
    text = reference_record()
    m = ConcreteModel()
    m.T, m.P = Var(initialize=353.15, bounds=(353.15, 354.15)), Var(initialize=45000.0, bounds=(44000.0, 46000.0))
    m.amount = ExternalFunction(library=library, function=FUNCTION)
    m.co2 = Constraint(expr=m.amount(m.T, m.P, *R_Q0[2:], text, "1", "0") / 0.03306804468947244 == 1)
    m.carbonate = Constraint(expr=m.amount(m.T, m.P, *R_Q0[2:], text, "0", "6") / 0.0006570839846201848 == 1)
    m.write(str(tmp_path / "roundtrip.nl"), format="nl")
    (tmp_path / "ipopt.opt").write_text("linear_solver mumps\ntol 1e-10\nacceptable_iter 0\n")
    external = subprocess.run([executable, "roundtrip.nl", "-AMPL"], cwd=tmp_path, capture_output=True, text=True, check=False)
    assert external.returncode == 0 and "epcsaft_equilibrium failure" not in external.stderr, external.stdout + external.stderr
    assert "Optimal Solution Found" in (tmp_path / "roundtrip.sol").read_text()
    solver.options["tol"] = 1e-10
    log = tmp_path / "ipopt.log"
    assert check_optimal_termination(solve_with_diagnostics(solver, m, logfile=log))
    assert "Optimal Solution Found" in log.read_text() and "epcsaft_equilibrium failure" not in log.read_text()  # not "Acceptable"
    assert abs(value(m.T) - 353.65) <= 6e-9 and abs(value(m.P) - 45500.0) <= 1.1e-5, (value(m.T), value(m.P))


H_ANCHOR = float.fromhex("-0x1.28afb4ecf2405p+21")  # -2430454.6156964325 J, the direct native SystemEnthalpy at R_Q0


def test_reference_enthalpy_equals_direct_engine() -> None:
    """E1: the system enthalpy (J) with its gradient and packed Hessian over q equals solve_equilibrium plus
    solved_state_actions on the declared problem bitwise at the anchor, where H = -2430454.6156964325 J, and to 1e-12
    relative at (353.65 K, 45500 Pa) and there with u_CO2 0.36 mol; each order returns the same value."""
    model, target = reference_mea()
    text = reference_record()
    for q, bitwise in ((R_Q0, True), ((353.65, 45500.0, *R_Q0[2:]), False), ((353.65, 45500.0, 0.36, *R_Q0[3:]), False)):
        values = [fgh(text, q, None, order, C)[0] for order in (0, 1)]  # before any Hessian at q
        got, want = callback_row(text, q, None, C), native(model, q[2:], target(q), enthalpy=True)[0]
        assert values == [got[0]] * 2
        if bitwise:
            assert got[0] == H_ANCHOR and [x.hex() for x in got] == [float(x).hex() for x in want]
        else:
            assert got == pytest.approx(want, rel=1e-12, abs=0.0), q


def test_reference_enthalpy_matches_central_differences() -> None:
    """E2: the system enthalpy with feed steps 1e-5 u_j x {4, 1, 0.25} mol."""
    central_differences([None], 4e-5)


def test_enthalpy_parameter_columns_match_central_differences() -> None:
    from epcsaft import IdealPolynomial
    from tests.native_inputs import record as thermal_record
    base = mixture()
    form = IdealPolynomial([30.0, 0.0], 300.0, 1000.0, 50.0)
    model = Mixture(base.parameters, thermochemistry=thermal_record(form, form),
                    active_parameters=[regression.active_parameter(base.parameters, "k_ij", ("isobutanol", "water")),
                                       regression.active_parameter(base.parameters, "dispersion_energy_over_k", "isobutanol")])
    text, q, eps = record(model), np.array([*TP0, *G0, *THETA0]), float(np.finfo(float).eps)
    _, gradient, hessian = fgh(text, q, None, 2)
    for b, h0 in ((4, 1e-5), (5, 1e-3)):
        for h in (4 * h0, h0, 0.25 * h0):
            plus, minus = (fgh(text, q + s * h * np.eye(6)[b], None, 1) for s in (1, -1))
            pairs = [(gradient[b], plus[0], minus[0], 1e-6)]
            pairs += [(hessian[a, b], plus[1][a], minus[1][a], 1e-5) for a in range(6)]
            for exact, fp, fm, tau in pairs:
                assert abs((fp - fm) / (2 * h) - exact) <= tau * abs(exact) + 10 * eps * (abs(fp) + abs(fm)) / (2 * h), (b, h)


def test_reference_adiabatic_flash_ipopt(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """E3: at 45 kPa and the anchor feed, Ipopt (tol 1e-10, "Optimal Solution Found" only) frees T in [353.15, 354.15] K
    from 353.15 K with H/H* = 1, H* = -2399045.378695539 J the direct native enthalpy at 353.65 K, and y equal to the vapor
    CO2 amount, so both callbacks read one state per q. With dH/dT / H* = -0.0295 1/K both rows stay unscaled, and the
    bounds (tol + 1e-12)/|dg/dT| and |da/dT| dT + tol + 1e-12 a* round up to 3.5e-9 K and 1.9e-10 mol about
    (353.65 K, a* = 0.040325995174790685 mol). The external AMPL Ipopt exits with code 0 and reports no callback failure."""
    from idaes.core.solvers import get_solver
    from epcsaft.idaes import asl_library, solve_with_diagnostics
    solver, executable = get_solver(), shutil.which("ipopt") or str(Path.home() / ".idaes/bin/ipopt")
    if not solver.available(exception_flag=False) or not Path(executable).is_file():
        pytest.skip("Ipopt executable unavailable")
    library = str(asl_library())
    monkeypatch.setenv("AMPLFUNC", library)
    text = reference_record()
    m = ConcreteModel()
    m.T, m.y = Var(initialize=353.15, bounds=(353.15, 354.15)), Var(initialize=0.0)
    m.enthalpy = ExternalFunction(library=library, function=ENTHALPY)
    m.amount = ExternalFunction(library=library, function=FUNCTION)
    m.energy = Constraint(expr=m.enthalpy(m.T, R_Q0[1], *R_Q0[2:], text) / float.fromhex("-0x1.24da2b079186ep+21") == 1)
    m.co2 = Constraint(expr=m.y == m.amount(m.T, R_Q0[1], *R_Q0[2:], text, "1", "0"))
    m.write(str(tmp_path / "adiabatic.nl"), format="nl")
    (tmp_path / "ipopt.opt").write_text("linear_solver mumps\ntol 1e-10\nacceptable_iter 0\n")
    external = subprocess.run([executable, "adiabatic.nl", "-AMPL"], cwd=tmp_path, capture_output=True, text=True, check=False)
    assert external.returncode == 0 and "epcsaft_equilibrium failure" not in external.stderr, external.stdout + external.stderr
    assert "Optimal Solution Found" in (tmp_path / "adiabatic.sol").read_text()
    solver.options["tol"] = 1e-10
    log = tmp_path / "ipopt.log"
    assert check_optimal_termination(solve_with_diagnostics(solver, m, logfile=log))
    assert "Optimal Solution Found" in log.read_text() and "epcsaft_equilibrium failure" not in log.read_text()  # not "Acceptable"
    assert abs(value(m.T) - 353.65) <= 3.5e-9 and abs(value(m.y) - 0.040325995174790685) <= 1.9e-10, (value(m.T), value(m.y))


def test_record_without_thermochemistry_refuses_only_enthalpy() -> None:
    """E4: the record with its thermochemistry block replaced by the absent flag is the record of the model without
    thermochemistry. There the enthalpy callback names missing_thermochemistry at orders 0, 1 and 2, and after each refusal
    the amount rows at that order equal those of the record with thermochemistry bitwise."""
    with_block, without = reference_record(), reference_record(False)
    start, end = with_block.index("1:1,17:standard_pressure,"), with_block.index("20:continuation_options,")
    assert with_block[:start] + "1:0," + with_block[end:] == without
    for order in (0, 1, 2):
        with pytest.raises(RuntimeError) as failure:
            fgh(without, R_Q0, None, order, C)
        assert str(failure.value).startswith("epcsaft_equilibrium:missing_thermochemistry:system.enthalpy:"), str(failure.value)
        assert "reference calorics need record anchors" in str(failure.value)
        for row in range(2 * C):
            got, want = (np.concatenate([np.ravel(x) for x in fgh(t, R_Q0, row, order, C) if x is not None]) for t in (without, with_block))
            assert got.tobytes() == want.tobytes(), (order, row)


# H: Held 2014 strategy 2 water/Na+/Cl- beside declared-present halite (tests/equilibrium/test_pure_solid.py declarations)
# at 298.15 K and 101325 Pa with 1 kg of water and 10 mol NaCl; q = (T, P, u_water, u_NaCl, eps/k[Na+], sigma[Cl-],
# k_ij[Cl-, water]). Rows 0-2 are the liquid amounts, row 3 the halite amount in mol of NaCl.
H_Q0, H_MATRIX = (298.15, 101325.0, 55.508435061791985, 10.0, 230.0, 2.756, -0.25), [[1.0, 0.0], [0.0, 1.0], [0.0, 1.0]]


def halite_problem(q: Any, amount: Any = None) -> eq.Problem:
    from tests.equilibrium import test_pure_solid as hs
    target = hs.problem(10.0, amount or eq.Free(), q[0], q[1])
    target.phases = [eq.Phase("L", kind="liquid", support=hs.SPECIES)]
    target.feed = eq.Amounts(dict(zip(hs.SPECIES, (q[2], q[3], q[3]), strict=True)))
    return target


def halite(shift: tuple[str, float] = ("", 0.0)) -> tuple[Mixture, str]:
    from tests.equilibrium import test_pure_solid as hs
    model = hs.held(shift, tuple(item for _, item, _ in hs.PARAMETERS))
    return model, eq.asl_record(model, halite_problem(H_Q0), OPTIONS, [0.0] * 3, H_MATRIX, list(H_Q0[2:4]), 0, 0.0)


def test_halite_record_equals_direct_engine() -> None:
    """H1: the saturated molality and halite amount; every row's value, gradient and Hessian equal the direct Engine
    bitwise at theta0 and to 1e-12 relative at k_ij + 0.01, and each order returns the same value."""
    from tests.equilibrium import test_pure_solid as hs
    model, text = halite()
    assert model._model.active_parameter_values == list(H_Q0[4:])
    sodium, solid = (fgh(text, H_Q0, row, 0, 3, 1)[0] for row in (1, 3))
    assert sodium / (H_Q0[2] * hs.MW) == pytest.approx(6.29890460910619, rel=1e-12, abs=0.0)
    assert solid == pytest.approx(3.7010953908938036, rel=1e-12, abs=0.0) and abs(solid + sodium - H_Q0[3]) <= 1e-13
    star = (*H_Q0[:6], -0.24)
    for q, mixture_at, bitwise in ((H_Q0, model, True), (star, halite(("pair/chloride-anion/water/k_ij", 0.01))[0], False)):
        want = native(mixture_at, q[2:4], halite_problem(q), matrix=H_MATRIX)
        values = [[fgh(text, q, row, order, 3, 1)[0] for row in range(4)] for order in (0, 1)]
        got = [callback_row(text, q, row, 3, 1) for row in range(4)]
        assert values == [[r[0] for r in got]] * 2
        for row in range(4):
            if bitwise:
                assert [x.hex() for x in got[row]] == [float(x).hex() for x in want[row]], (q, row)
            else:
                assert got[row] == pytest.approx(want[row], rel=1e-12, abs=0.0), (q, row)


def test_halite_saturation_condition() -> None:
    """H2: at the callback's liquid the property record's solid row closes sum nu (ln x + ln phi) = l_s to 1e-12."""
    from tests.eos.test_asl_property_record import halite as property_halite, property_fgh
    _, text = halite()
    _, row_text, _, _ = property_halite()
    amounts = np.array([fgh(text, H_Q0, row, 0, 3, 1)[0] for row in range(3)])
    q = [*H_Q0[:2], *(amounts / amounts.sum()), *H_Q0[4:]]
    ln_x = np.log(q[3]) + np.log(q[4])
    residual = ln_x + property_fgh(row_text, q, 3, 0)[0] + property_fgh(row_text, q, 4, 0)[0] - property_fgh(row_text, q, 5, 0)[0]
    assert abs(residual) <= 1e-12, residual


def test_halite_columns_match_central_differences() -> None:
    """H3, H4, verification only: T at {4, 1, 0.25} 1e-3 K, P at {4, 1, 0.25} 1e3 Pa, u at {4, 1, 0.25} 1e-5 mol and theta at
    {4, 1, 0.25} 1e-5 max(1, |theta|). Every gradient entry (tau 1e-6) and every Hessian entry against the difference of the
    exact gradient (tau 1e-5) meet |D - A| <= tau |A| + 10 eps (|f+| + |f-|)/2h. The halite column is the design's native
    one; dn_s/du_w = -m_s M_w and dn_s/du_NaCl = 1 (the saturated liquid does not move with the feed), so the u_NaCl
    Hessian row is zero and the liquid sodium row is the negative of the halite row elsewhere. The liquid water row is u_w
    exactly: its columns close that balance to 1e-14 of the largest halite entry of the same order. The u_NaCl Hessian column
    is that analytic zero, so its differences carry only re-solve noise (25 ulp of the gradient): 1e-13 replaces 10 eps there."""
    from tests.equilibrium import test_pure_solid as hs
    _, text = halite()
    eps, q = float(np.finfo(float).eps), np.array(H_Q0)
    exact = [fgh(text, q, row, 2, 3, 1) for row in range(4)]
    solid = exact[3][1]
    assert solid == pytest.approx([-1.430830593594595e-02, -2.1295896304825715e-09, -0.11347653022633869, 1.0,
                                   8.92286053703689e-03, 2.304237756233167, -4.849349365690943], rel=1e-12, abs=0.0)
    assert abs(solid[2] + exact[1][0] / (H_Q0[2] * hs.MW) * hs.MW) <= 1e-15 and abs(solid[3] - 1.0) <= 1e-15
    assert not exact[1][2][3].any() and not exact[3][2][3].any()
    water = exact[0]
    # The two identities below hold exactly in the equations and to one ulp of their own scale in binary64. Bitwise
    # equality was rounding luck: Ipopt's unrelaxed bounds alone move each by one ulp (the same Engine reproduces it with
    # an options file setting only bound_relax_factor 0).
    assert abs(water[0] - H_Q0[2]) <= eps * H_Q0[2] and np.abs(water[1] - np.eye(7)[2]).max() <= 1e-14 * np.abs(solid).max()
    assert np.abs(water[2]).max() <= 1e-14 * np.abs(exact[3][2]).max()
    assert np.delete(exact[1][1] + solid, 3) == pytest.approx(0.0, abs=1e-15)
    assert np.delete(exact[1][2] + exact[3][2], 3) == pytest.approx(0.0, abs=eps * np.abs(exact[3][2]).max())
    for b, h0 in enumerate((1e-3, 1e3, 1e-5, 1e-5, *(1e-5 * max(1.0, abs(t)) for t in H_Q0[4:]))):
        for h in (4 * h0, h0, 0.25 * h0):
            plus, minus = ([fgh(text, q + s * h * np.eye(7)[b], row, 1, 3, 1) for row in range(4)] for s in (1, -1))
            for row, (_, gradient, hessian) in enumerate(exact[1:], 1):
                pairs = [(gradient[b], plus[row][0], minus[row][0], 1e-6)]
                pairs += [(hessian[a, b], plus[row][1][a], minus[row][1][a], 1e-5) for a in range(7)]
                noise = 1e-13 if b == 3 else 10 * eps
                for exact_value, fp, fm, tau in pairs:
                    assert abs((fp - fm) / (2 * h) - exact_value) <= tau * abs(exact_value) + noise * (abs(fp) + abs(fm)) / (2 * h), (row, b, h)


def test_halite_onset_and_interval_name_the_cause_and_recover() -> None:
    """H5: below the solubility (6.0 mol NaCl) halite is absent, and above its law (308.16 K) the solve is refused; each
    names its cause, and the next call at q0 returns the q0 value bitwise. At 6.30 mol halite is present."""
    from epcsaft.idaes import PropertyExternalFunction, asl_library
    _, text = halite()
    callback = PropertyExternalFunction(library=str(asl_library()), function=SOLID)
    start = callback.evaluate([*H_Q0, text, "0"])
    for q, cause in (((*H_Q0[:3], 6.0, *H_Q0[4:]), "phase_boundary_not_differentiable:solid[0].amount: solid halite is absent at the solved state"),
                     ((308.16, *H_Q0[1:]), "outside_interval:solid[0].amount: solid halite: temperature outside reaction correlation")):
        with pytest.raises(RuntimeError) as failure:
            callback.evaluate([*q, text, "0"])
        assert str(failure.value).startswith("epcsaft_equilibrium:" + cause), str(failure.value)
        assert callback.evaluate([*H_Q0, text, "0"]).hex() == start.hex()
    assert callback.evaluate([*H_Q0[:3], 6.30, *H_Q0[4:], text, "0"]) == pytest.approx(1.0953908938e-3, rel=1e-10, abs=0.0)


def test_halite_temperature_ipopt(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """H6: with u and theta fixed, Ipopt (exact Hessian, tol 1e-12, "Optimal Solution Found" only) frees T in [288.15,
    308.15] K from 298.15 K and recovers 303.15 K from the halite amount there, the native 3.633567004371989 mol: dn_s/dT =
    -0.0127 mol/K at 303.15 K, so 1e-8 K is 1.27e-10 mol and the unscaled amount residual must be at most 1e-11 mol. The
    external AMPL Ipopt returns the same temperature and reports no callback failure."""
    from idaes.core.solvers import get_solver
    from pyomo.opt import ReaderFactory  # type: ignore[import-untyped]
    from epcsaft.idaes import PropertyExternalFunction, asl_library, solve_with_diagnostics
    solver, executable = get_solver(), shutil.which("ipopt") or str(Path.home() / ".idaes/bin/ipopt")
    if not solver.available(exception_flag=False) or not Path(executable).is_file():
        pytest.skip("Ipopt executable unavailable")
    library, target = str(asl_library()), 3.633567004371989
    monkeypatch.setenv("AMPLFUNC", library)
    _, text = halite()
    m = ConcreteModel()
    m.T = Var(initialize=298.15, bounds=(288.15, 308.15))
    m.solid = ExternalFunction(library=library, function=SOLID)
    m.amount = Constraint(expr=m.solid(m.T, *H_Q0[1:], text, "0") == target)
    m.write(str(tmp_path / "halite.nl"), format="nl")
    (tmp_path / "ipopt.opt").write_text("linear_solver mumps\ntol 1e-12\nacceptable_iter 0\n")
    external = subprocess.run([executable, "halite.nl", "-AMPL"], cwd=tmp_path, capture_output=True, text=True, check=False)
    assert external.returncode == 0 and "epcsaft_equilibrium failure" not in external.stderr, external.stdout + external.stderr
    sol = ReaderFactory("sol")(str(tmp_path / "halite.sol"))
    assert "Optimal Solution Found" in str(sol.solver.message)
    external_T = sol.solution(0).variable["v0"]["Value"]
    solver.options["tol"] = 1e-12
    log = tmp_path / "ipopt.log"
    result = solve_with_diagnostics(solver, m, logfile=log)
    assert check_optimal_termination(result) and result.solver.error_rc == 0
    assert "Optimal Solution Found" in log.read_text() and "epcsaft_equilibrium failure" not in log.read_text()
    residual = PropertyExternalFunction(library=library, function=SOLID).evaluate([value(m.T), *H_Q0[1:], text, "0"]) - target
    print(f"H6 T={value(m.T)!r} K dT={value(m.T) - 303.15:.3e} K residual={residual:.3e} mol external T={external_T!r} K")
    assert abs(value(m.T) - 303.15) <= 1e-8 and abs(residual) <= 1e-11 and abs(external_T - 303.15) <= 1e-8


def test_halite_record_refusals() -> None:
    """H7: a Pinned(0) and a scalar-K halite are refused when the record is built; with halite present the system enthalpy
    is refused at every order naming the solid, and in a native batch it loses its value and action while the amounts keep
    theirs."""
    from tests.equilibrium import test_pure_solid as hs
    model, text = halite()
    with pytest.raises(_core.ParameterError, match="invalid ASL affine feed/phase contract"):
        eq.asl_record(model, halite_problem(H_Q0, eq.Pinned(0.0)), OPTIONS, [0.0] * 3, H_MATRIX, list(H_Q0[2:4]), 0, 0.0)
    scalar = halite_problem(H_Q0)
    scalar.solids = [dataclasses.replace(scalar.solids[0], correlation=None, equilibrium_constant=37.9)]
    with pytest.raises(ValueError, match=r"ASL solid rows require a ln K\(T\) law"):
        eq.asl_record(model, scalar, OPTIONS, [0.0] * 3, H_MATRIX, list(H_Q0[2:4]), 0, 0.0)
    for order in (0, 1, 2):
        with pytest.raises(RuntimeError) as failure:
            fgh(text, H_Q0, None, order, 3, 1)
        assert str(failure.value).startswith("epcsaft_equilibrium:missing_thermochemistry:system.enthalpy: solid halite has no enthalpy record")
    target = halite_problem(H_Q0)
    central = eq.solve_equilibrium(model, target)
    kinds = (eq.SolvedStateObservableKind.SystemEnthalpy, eq.SolvedStateObservableKind.SolidAmount)
    request = eq.SolvedStateActionRequest([eq.SolvedStateObservable(k) for k in kinds],
                                          [eq.SolvedStateActionDirection(0.0, 0.0, [1.0, 0.0, 0.0], [0.0] * 3)])
    enthalpy, solid = _core.solved_state_actions(_core.compile_problem(model._model, target._native(model)), central, request).results
    assert (enthalpy.status, enthalpy.value, enthalpy.action) == (_core.SolvedStateActionStatus.MissingThermochemistry, None, None)
    assert solid.value == central.solid_amounts[0] and solid.action == pytest.approx(-hs.MW * 6.29890460910619, rel=1e-12)


def raw(function: str, arguments: list[Any], order: int) -> Any:
    """The wrapper's raw call with every argument variable and a buffer for each derivative order through `order`."""
    from epcsaft.idaes import PropertyExternalFunction, asl_library
    return PropertyExternalFunction(library=str(asl_library()), function=function)._evaluate(arguments, [False] * len(arguments), order)


def test_first_order_names_compute_no_hessian() -> None:
    """Q1: each `_first_order` name, at a state it enters first, returns the direct Engine's value and gradient bitwise and
    NaN in every Hessian entry; the exact amount name called next at that state returns the direct Engine's Hessian."""
    (model, (h_model, h_text)), (r_model, target) = (mixture(), halite()), reference_mea()
    cases = [(FUNCTION, record(model), [*TP0, *G0, *THETA0], [*TP0, *G1, *THETA0], ["0", "1"], native(model, G0)[1]),
             (ENTHALPY, reference_record(), [*R_Q0], [353.65, 45500.0, *R_Q0[2:]], [], native(r_model, R_Q0[2:], target(R_Q0), enthalpy=True)[0]),
             (SOLID, h_text, [*H_Q0], [*H_Q0[:6], -0.24], ["0"], native(h_model, H_Q0[2:4], halite_problem(H_Q0), matrix=H_MATRIX)[3])]
    for function, text, q, away, strings, want in cases:
        n = len(q)
        raw(function, [*away, text, *strings], 0)  # the record now holds another state
        v, g, h = raw(function + "_first_order", [*q, text, *strings], 2)
        assert [x.hex() for x in (v, *g[:n])] == [float(x).hex() for x in want[:1 + n]] and np.isnan(h[:n * (n + 1) // 2]).all(), function
        if function == FUNCTION:
            v, g, h = raw(function, [*q, text, *strings], 2)
            assert [x.hex() for x in (v, *g[:n], *h[:n * (n + 1) // 2])] == [float(x).hex() for x in want]
