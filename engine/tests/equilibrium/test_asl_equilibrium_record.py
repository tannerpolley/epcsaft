"""The ASL equilibrium record on the engine/README.md water/isobutanol record (298.15 K, 1e5 Pa, feed (0.7, 0.3) mol) with
k_ij(isobutanol, water) and isobutanol eps/k active, q = (u_water, u_isobutanol, k_ij, eps/k). F1: values and exact actions
equal the direct Engine (bitwise at theta0 and for a record with no active parameters; 1e-12 relative at theta*). F2: the
columns match central differences (verification only). F3: two temperatures recover theta* through Ipopt. F4: where the
liquids merge the callback names the cause, and the next calls recover. M1-M3 repeat F1-F3 on a declared vapor-liquid
pair: the packet-3 MEA/CO2/water flash (charged, reacting liquid; neutral vapor support) at 318.15 K and 5918.5286 Pa."""
from __future__ import annotations

import dataclasses
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
from tests.equilibrium.test_trace_reliability import _MEA, _MEA_REACTIONS, _mea_mapping
from tests.packets import parameters_json

NAMES, G0, G1 = ("water", "isobutanol"), (0.7, 0.3), (0.7, 0.285)
THETA0, THETA_STAR, THETA_MERGED = (-0.0071, 257.13), (-0.0021, 258.13), (-0.0371, 257.13)
OPTIONS = eq.Continuation(initial_step=0.05, minimum_step=1e-3, max_steps=200, max_retries_per_step=5)
FUNCTION = "epcsaft_equilibrium_phase_component_amount"


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


def fgh(text: str, q: Any, row: int, order: int, components: int = 2) -> tuple[float, Any, Any]:
    """Row p*components + c at q; the Hessian unpacked from the lower triangle i + j(j+1)/2."""
    from epcsaft.idaes import PropertyExternalFunction, asl_library
    callback = PropertyExternalFunction(library=str(asl_library()), function=FUNCTION)
    v, g, h = callback.evaluate_fgh([*q, text, str(row // components), str(row % components)], fgh=order)
    n = len(q)
    hessian = np.array([[h[min(i, j) + max(i, j) * (max(i, j) + 1) // 2] for j in range(n)] for i in range(n)]) if order == 2 else None
    return float(v), np.array(g[:n]) if order else None, hessian


def native(model: Mixture, feed: Any, target: eq.Problem | None = None) -> list[list[float]]:
    """Per row: value, first actions, second actions (i <= j, packed order), from solve_equilibrium and solved_state_actions;
    the feed controls are the first len(feed) components."""
    target = target or problem(feed)
    central = eq.solve_equilibrium(model, target)
    compiled = _core.compile_problem(model._model, target._native(model))
    k, n, m = len(feed), len(central.mole_fractions) // 2, len(model._model.active_parameter_values)
    e = np.eye(k + m).tolist()
    d = [eq.SolvedStateActionDirection(0.0, 0.0, column[:k] + [0.0] * (n - k), column[k:]) for column in e]
    observables = [eq.SolvedStateObservable(eq.SolvedStateObservableKind.PhaseComponentAmount, p, c) for p in range(2) for c in range(n)]

    def batch(*directions: Any) -> Any:
        return _core.solved_state_actions(compiled, central, eq.SolvedStateActionRequest(observables, list(directions))).results

    zero = batch(eq.SolvedStateActionDirection(0.0, 0.0, [0.0] * n, [0.0] * m))
    first = [batch(d[i]) for i in range(k + m)]
    second = [batch(d[i], d[j]) for j in range(k + m) for i in range(j + 1)]
    return [[zero[r].value, *(b[r].action for b in first), *(b[r].action for b in second)] for r in range(2 * n)]


def callback_row(text: str, q: Any, row: int, components: int = 2) -> list[float]:
    v, g, h = fgh(text, q, row, 2, components)
    return [v, *g, *(h[i, j] for j in range(len(q)) for i in range(j + 1))]


def test_record_equals_direct_engine() -> None:
    """F1: at theta0 (G0) and with no active parameters (G0, G1) every row's value, gradient and Hessian equal the direct
    Engine bitwise, and each order returns the same value; at theta* (G0, G1) they match the edited packet to 1e-12."""
    model, inactive = mixture(), mixture(active=False)
    assert model._model.active_parameter_values == list(THETA0)
    cases = [(record(model), [*G0, *THETA0], native(model, G0), True)]
    cases += [(record(inactive), list(G), native(inactive, G), True) for G in (G0, G1)]
    cases += [(record(model), [*G, *THETA_STAR], native(mixture(THETA_STAR), G), False) for G in (G0, G1)]
    for text, q, want, bitwise in cases:
        values = [[fgh(text, q, row, order)[0] for row in range(4)] for order in (0, 1)]  # before any Hessian at q
        got = [callback_row(text, q, row) for row in range(4)]
        assert values == [[r[0] for r in got]] * 2
        for row in range(4):
            if bitwise:
                assert [x.hex() for x in got[row]] == [float(x).hex() for x in want[row]], (q, row)
            else:
                assert got[row] == pytest.approx(want[row], rel=1e-12, abs=0.0), (q, row)


def test_columns_match_central_differences() -> None:
    """F2, verification only: feed (h 1e-5 mol), k_ij (h 1e-5) and eps/k (h 1e-3 K) at h x {4, 1, 0.25}: every gradient
    entry (tau 1e-6) and every Hessian entry with a parameter (tau 1e-5) meet |D - A| <= tau|A| + 10 eps (|f+| + |f-|)/2h."""
    text, eps, q = record(mixture()), float(np.finfo(float).eps), np.array([*G0, *THETA0])
    exact = [fgh(text, q, row, 2) for row in range(4)]
    for b, h0 in ((0, 1e-5), (1, 1e-5), (2, 1e-5), (3, 1e-3)):
        for h in (4 * h0, h0, 0.25 * h0):
            plus, minus = ([fgh(text, q + s * h * np.eye(4)[b], row, 1) for row in range(4)] for s in (1, -1))
            for row, (_, gradient, hessian) in enumerate(exact):
                pairs = [(gradient[b], plus[row][0], minus[row][0], 1e-6)]
                pairs += [(hessian[a, b], plus[row][1][a], minus[row][1][a], 1e-5) for a in range(4) if b >= 2 or a >= 2]
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
        m.add_component(f"target_{k}", Constraint(expr=m.amount(*G0, m.theta[0], m.theta[1], record(model, T), "1", "1") / target == 1))
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
    star, start = [callback_row(text, [*G0, *THETA_STAR], row) for row in range(4)], native(model, G0)
    callback = PropertyExternalFunction(library=str(asl_library()), function=FUNCTION)
    with pytest.raises(RuntimeError) as failure:
        callback.evaluate([*G0, *THETA_MERGED, text, "1", "1"])
    assert str(failure.value).startswith("epcsaft_equilibrium:phase_boundary_not_differentiable:")
    assert "parameter continuation: phase_disappeared" in str(failure.value)
    assert [[x.hex() for x in callback_row(text, [*G0, *THETA_STAR], row)] for row in range(4)] == [[x.hex() for x in r] for r in star]
    assert [[x.hex() for x in callback_row(text, [*G0, *THETA0], row)] for row in range(4)] == [[float(x).hex() for x in r] for r in start]



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
            got = callback_row(text, u, row, C)
            assert got[0] == pytest.approx(amounts[row], rel=1e-12, abs=0.0), (u, row)
            assert got[1:] == pytest.approx(want[row][1:], rel=1e-12, abs=0.0), (u, row)
        assert [callback_row(text, u, C + c, C)[0] for c in range(3, C)] == [0.0] * 6


def test_vapor_liquid_columns_match_central_differences() -> None:
    """M2, verification only: each feed control at h = 1e-4 u0 x {1, 1/4, 1/16}; every gradient entry meets
    |D - A| <= 1e-6 |A| + 10 max(eps, row_tolerance) (|f+| + |f-|)/2h, since a re-solve reproduces amounts to its row tolerance."""
    text, noise, q = mea_record(mea_model()), 10 * max(float(np.finfo(float).eps), 1e-12), np.array(U0)
    exact = [fgh(text, q, row, 1, C)[1] for row in range(2 * C)]
    for b in range(3):
        for h in (1e-4 * U0[b] * s for s in (1.0, 0.25, 0.0625)):
            plus, minus = ([fgh(text, q + s * h * np.eye(3)[b], row, 0, C)[0] for row in range(2 * C)] for s in (1, -1))
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
    m.target = Constraint(expr=m.amount(m.u, *U0[1:], mea_record(mea_model()), "0", "3") == 6.302671442747203)
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
                          (text.replace("equilibrium.v5", "equilibrium.v4"), "")):
        with pytest.raises(RuntimeError, match="epcsaft_equilibrium:.*" + cause):
            callback.evaluate([*U0, edited, "0", "3"])
