"""The ASL equilibrium record on the engine/README.md water/isobutanol record (298.15 K, 1e5 Pa, feed (0.7, 0.3) mol) with
k_ij(isobutanol, water) and isobutanol eps/k active, q = (u_water, u_isobutanol, k_ij, eps/k). F1: values and exact actions
equal the direct Engine (bitwise at theta0 and for a record with no active parameters; 1e-12 relative at theta*). F2: the
columns match central differences (verification only). F3: two temperatures recover theta* through Ipopt. F4: where the
liquids merge the callback names the cause, and the next calls recover."""
from __future__ import annotations

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


def fgh(text: str, q: Any, row: int, order: int) -> tuple[float, Any, Any]:
    """Row p*2 + c at q; the Hessian unpacked from the lower triangle i + j(j+1)/2."""
    from epcsaft.idaes import PropertyExternalFunction, asl_library
    callback = PropertyExternalFunction(library=str(asl_library()), function=FUNCTION)
    v, g, h = callback.evaluate_fgh([*q, text, str(row // 2), str(row % 2)], fgh=order)
    n = len(q)
    hessian = np.array([[h[min(i, j) + max(i, j) * (max(i, j) + 1) // 2] for j in range(n)] for i in range(n)]) if order == 2 else None
    return float(v), np.array(g[:n]) if order else None, hessian


def native(model: Mixture, feed: Any) -> list[list[float]]:
    """Per row: value, first actions, second actions (i <= j, packed order), from solve_equilibrium and solved_state_actions."""
    target = problem(feed)
    central = eq.solve_equilibrium(model, target)
    compiled = _core.compile_problem(model._model, target._native(model))
    m = len(model._model.active_parameter_values)
    e = np.eye(2 + m).tolist()
    d = [eq.SolvedStateActionDirection(0.0, 0.0, column[:2], column[2:]) for column in e]
    observables = [eq.SolvedStateObservable(eq.SolvedStateObservableKind.PhaseComponentAmount, p, c) for p in range(2) for c in range(2)]

    def batch(*directions: Any) -> Any:
        return _core.solved_state_actions(compiled, central, eq.SolvedStateActionRequest(observables, list(directions))).results

    zero = batch(eq.SolvedStateActionDirection(0.0, 0.0, [0.0, 0.0], [0.0] * m))
    first = [batch(d[i]) for i in range(2 + m)]
    second = [batch(d[i], d[j]) for j in range(2 + m) for i in range(j + 1)]
    return [[zero[r].value, *(b[r].action for b in first), *(b[r].action for b in second)] for r in range(4)]


def callback_row(text: str, q: Any, row: int) -> list[float]:
    v, g, h = fgh(text, q, row, 2)
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

