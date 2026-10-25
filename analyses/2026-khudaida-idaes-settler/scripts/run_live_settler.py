"""Installed-wheel qualification of the live fixed-T/P Khudaida ASL settler (#76)."""

import argparse
import ctypes
import importlib.util
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import epcsaft
import numpy as np
import run_fixed_settler as fixed
from epcsaft import equilibrium as eq
from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from pyomo.contrib.pynumero.interfaces.ampl_nlp import AslNLP
from pyomo.core.base.external import _ARGLIST
from pyomo.environ import (
    ConcreteModel,
    Constraint,
    ExternalFunction,
    Objective,
    Suffix,
    Var,
    check_optimal_termination,
    value,
)
from pyomo.environ import units as u

spec = importlib.util.spec_from_file_location("neutral_asl", fixed.REPO / "analyses/2026-idaes-asl-adapter/scripts/run.py")
neutral = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = neutral
spec.loader.exec_module(neutral)
F0 = np.array((0.7198946513021484, 0.03377833329223333, 0.21864518731791907, 0.02768182808769922))
OPTIONS = eq.Continuation(initial_step=0.25, minimum_step=0.25, max_steps=4, max_retries_per_step=1)
NAME = "epcsaft_equilibrium_phase_component_amount"
FLOW = u.mol / u.s


def feed(F):
    return (*F[:3], F[3], F[3])


def encode_record(model, problem, options, anchor):
    return eq.asl_record(
        model, problem, options, [0.0] * 5, np.array([feed(row) for row in np.eye(4)]).T.tolist(), list(anchor), 0, fixed.GAP
    )


class Reference:
    def __init__(self):
        self.model = epcsaft.Mixture(epcsaft.Parameters.from_json(fixed.PARAMETERS, components=fixed.COMPONENTS))
        self.problem = fixed.engine_problem(feed(F0))
        self.anchor = eq.solve_equilibrium(self.model, self.problem)
        self.record = encode_record(self.model, self.problem, OPTIONS, F0)
        self.observables = [
            eq.SolvedStateObservable(eq.SolvedStateObservableKind.PhaseComponentAmount, p, i) for p in range(2) for i in range(5)
        ]
        self.directions = [eq.SolvedStateActionDirection(0.0, 0.0, feed(row), []) for row in np.eye(4)]
        _, _, self.library = neutral.installed_identity()
        self.external = ExternalFunction(library=str(self.library), function=NAME, units=FLOW, arg_units=[FLOW] * 4 + [None] * 3)

    def public(self, F, order=2):
        problem = fixed.engine_problem(feed(F))
        continuation = (
            None if np.array_equal(F, F0) else eq.continue_equilibrium(self.model, self.problem, problem, OPTIONS, start_result=self.anchor)
        )
        if continuation is not None and (not continuation.success or continuation.stop != eq.ContinuationStop.reached_target):
            raise RuntimeError(f"continuation: {continuation.stop}: {continuation.message}")
        central = self.anchor if continuation is None else continuation.states[-1]
        _, diagnostics, failures = fixed.phase_data(self.model, central, feed(F))
        if failures or central.mole_fractions[0] - central.mole_fractions[5] <= fixed.GAP:
            raise RuntimeError(f"branch/certificate: {failures}, {diagnostics}")
        compiled = eq.compile_problem(self.model, problem)
        y, g, h, statuses = np.zeros(10), np.zeros((10, 4)), np.zeros((10, 4, 4)), []
        pairs = [(-1,)] if order == 0 else [(i,) for i in range(4)]
        if order == 2:
            pairs += [(i, j) for j in range(4) for i in range(j + 1)]
        for indices in pairs:
            directions = (
                [eq.SolvedStateActionDirection(0.0, 0.0, [0.0] * 5, [])] if indices == (-1,) else [self.directions[i] for i in indices]
            )
            batch = eq.solved_state_actions(compiled, central, eq.SolvedStateActionRequest(self.observables, directions))
            for row, item in enumerate(batch.results):
                statuses.append(item.status.name)
                if (
                    item.status != eq.SolvedStateActionStatus.Available
                    or item.value is None
                    or item.action is None
                    or not np.isfinite([item.value, item.action]).all()
                ):
                    raise RuntimeError(f"unavailable E5 action: {item.status}")
                y[row] = item.value
                if indices == (-1,):
                    continue
                if len(indices) == 1:
                    g[row, indices[0]] = item.action
                else:
                    i, j = indices
                    h[row, i, j] = h[row, j, i] = item.action
        assert batch.diagnostics.rank == batch.diagnostics.state_dimension
        return (y, g, h), {
            "rank": batch.diagnostics.rank,
            "dimension": batch.diagnostics.state_dimension,
            "statuses": sorted(set(statuses)),
            "classified_residual": diagnostics["classified"],
            "reduced_hessian_minimum": diagnostics["hessian"],
            "steps": 0 if continuation is None else continuation.accepted_steps,
            "requested_tolerance_met": central.requested_tolerance_met,
            "water_gap": central.mole_fractions[0] - central.mole_fractions[5],
        }

    def asl(self, F, order=2, rows=range(10)):
        y, g, h = np.zeros(10), np.zeros((10, 4)), np.zeros((10, 4, 4))
        for row in rows:
            y[row], gradient, packed = self.external.evaluate_fgh([*F, self.record, str(row // 5), str(row % 5)], fgh=order)
            if order >= 1:
                g[row] = gradient[:4]
            if order == 2:
                h[row] = [[packed[min(i, j) + max(i, j) * (max(i, j) + 1) // 2] for j in range(4)] for i in range(4)]
        if not all(np.isfinite(array).all() for array in (y, g, h)):
            raise RuntimeError("nonfinite ASL response")
        return y, g, h


def stage(reference, start, optimize):
    F = F0 * np.array((1.0, 1.0, *start))
    (y, _, _), _ = reference.public(F, 0)
    phase_data, _, failures = fixed.phase_data(reference.model, reference.anchor, feed(F0))
    assert not failures
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.settler = fixed.KhudaidaSettler(
        dynamic=False, has_holdup=False, inlet_flows=feed(F0), phase_data=[(p["label"], p) for p in phase_data]
    )
    b = m.fs.settler
    for name in ("aqueous_flow_mapping", "organic_flow_mapping", "mapping_factor", "accepted_phase_amount", "accepted_phase_mole_fraction"):
        b.del_component(name)
    m.si, m.ss = Var(bounds=(0.85, 1.15), initialize=start[0]), Var(bounds=(0.85, 1.15), initialize=start[1])
    if not optimize:
        m.si.fix()
        m.ss.fix()
    for i in (2, 3, 4):
        b.flow_mol_in[fixed.COMPONENTS[i]].unfix()
        b.flow_mol_in[fixed.COMPONENTS[i]].set_value(feed(F)[i])
    b.feed_links = Constraint(
        (2, 3, 4), rule=lambda b, i: b.flow_mol_in[fixed.COMPONENTS[i]] == (F0[2] * m.si if i == 2 else F0[3] * m.ss) * FLOW
    )
    m.response = ExternalFunction(library=str(reference.library), function=NAME, units=FLOW, arg_units=[FLOW] * 4 + [None] * 3)
    flows = [b.flow_mol_aqueous, b.flow_mol_organic]
    b.live_outlets = Constraint(
        range(2),
        range(5),
        rule=lambda b, p, i: (
            (
                flows[p][fixed.COMPONENTS[i]]
                - m.response(*[b.flow_mol_in[c] for c in fixed.COMPONENTS[:4]], reference.record, str(p), str(i))
            )
            / FLOW
            == 0
        ),
    )
    for row in range(10):
        flows[row // 5][fixed.COMPONENTS[row % 5]].set_value(y[row])
    m.recovery = Constraint(expr=0.965 - b.flow_mol_organic["ethanol"] / b.flow_mol_in["ethanol"] <= 0)
    m.objective = Objective(expr=m.si + m.ss)
    for name in ("dual", "ipopt_zL_out", "ipopt_zU_out"):
        setattr(m, name, Suffix(direction=Suffix.IMPORT))
    assert degrees_of_freedom(m) == (2 if optimize else 0)
    return m


def nl_checks(m, path, reference, compare=True):
    os.environ["AMPLFUNC"] = str(reference.library)
    nlp = AslNLP(str(path))
    columns, rows = neutral.named_file_rows(path)
    variables = [m.find_component(name) for name in columns]
    constraints = [m.find_component(name) for name in rows[: nlp.n_constraints()]]
    x = np.array([value(v) for v in variables])
    lam = -np.array([m.dual.get(c, 0.0) for c in constraints])
    lower, upper = nlp.primals_lb(), nlp.primals_ub()
    zl, zu = np.array([m.ipopt_zL_out.get(v, 0.0) for v in variables]), -np.array([m.ipopt_zU_out.get(v, 0.0) for v in variables])
    nlp.set_primals(x)
    nlp.set_duals(lam)
    nlp.set_obj_factor(1.0)
    c, jac = nlp.evaluate_constraints(), nlp.evaluate_jacobian().toarray()
    cl, cu = nlp.constraints_lb(), nlp.constraints_ub()
    masks = (np.isfinite(lower), np.isfinite(upper), np.isfinite(cl), np.isfinite(cu))
    comp = np.r_[zl[masks[0]] * (x - lower)[masks[0]], zu[masks[1]] * (upper - x)[masks[1]], lam[masks[3]] * (cu - c)[masks[3]]]
    metrics = {
        "stationarity": float(np.max(abs(nlp.evaluate_grad_objective() + jac.T @ lam - zl + zu))),
        "feasibility": float(
            max(
                np.max(np.maximum(cl - c, 0)),
                np.max(np.maximum(c - cu, 0)),
                np.max(np.maximum(lower - x, 0)),
                np.max(np.maximum(x - upper, 0)),
            )
        ),
        "multiplier_sign": float(
            max(0.0, np.max(-zl), np.max(-zu), *[-lam[i] for i in range(len(c)) if cl[i] < cu[i] and np.isfinite(cu[i])])
        ),
        "complementarity": float(np.max(abs(comp))),
    }
    metrics.update(
        {"lambda": dict(zip(rows[: len(lam)], lam.tolist())), "zL": dict(zip(columns, zl.tolist())), "zU": dict(zip(columns, zu.tolist()))}
    )
    if compare:
        F = np.array([value(m.fs.settler.flow_mol_in[i]) for i in fixed.COMPONENTS[:4]])
        (y, g, h), _ = reference.public(F)
        D = np.diag(F0[2:4])
        P = np.zeros((len(columns), 2))
        for name, axis, factor in (
            ("si", 0, 1.0),
            ("ss", 1, 1.0),
            ("fs.settler.flow_mol_in[isobutanol]", 0, F0[2]),
            ("fs.settler.flow_mol_in[sodium-cation]", 1, F0[3]),
            ("fs.settler.flow_mol_in[chloride-anion]", 1, F0[3]),
        ):
            P[columns.index(name), axis] = factor
        weighted = np.zeros((2, 2))
        errors = []
        for row in range(10):
            index = rows.index(f"fs.settler.live_outlets[{row // 5},{row % 5}]")
            outlet = m.fs.settler.flow_mol_aqueous if row < 5 else m.fs.settler.flow_mol_organic
            expected = np.zeros(len(columns))
            expected[columns.index(outlet[fixed.COMPONENTS[row % 5]].name)] = 1.0
            expected[[columns.index(f"fs.settler.flow_mol_in[{i}]") for i in fixed.COMPONENTS[2:4]]] = -g[row, 2:4]
            errors.append(float(np.max(abs(jac[index] - expected))))
            errors.extend(
                (
                    abs(c[index] - (value(outlet[fixed.COMPONENTS[row % 5]]) - y[row])),
                    float(np.linalg.norm(jac[index] @ P + g[row, 2:4] @ D, ord=np.inf)),
                )
            )
            weighted -= lam[index] * D @ h[row, 2:4, 2:4] @ D
            dual = np.zeros(len(c))
            dual[index] = 1.0
            nlp.set_duals(dual)
            nlp.set_obj_factor(0.0)
            errors.append(float(np.linalg.norm(P.T @ nlp.evaluate_hessian_lag().toarray() @ P + D @ h[row, 2:4, 2:4] @ D, ord=np.inf)))
        nlp.set_duals(lam)
        nlp.set_obj_factor(1.0)
        errors.append(float(np.linalg.norm(P.T @ nlp.evaluate_hessian_lag().toarray() @ P - weighted, ord=np.inf)))
        metrics["nl_transport_max_abs"] = max(errors)
        assert max(errors) <= 1e-9
    assert all(
        metrics[k] <= tolerance
        for k, tolerance in (("stationarity", 1e-7), ("feasibility", 1e-8), ("multiplier_sign", 1e-8), ("complementarity", 1e-7))
    ), metrics
    return metrics


def solve(reference, start, output, ipopt, optimize=True):
    output.mkdir(parents=True, exist_ok=True)
    m = stage(reference, start, optimize)
    path = output / "stage.nl"
    m.write(str(path), format="nl", io_options=neutral.NL_IO_OPTIONS)
    options = {**neutral.SOLVER_OPTIONS, "print_user_options": "yes", "file_print_level": 5, "output_file": str(output / "ipopt.log")}
    if optimize:
        options.update(neutral.DERIVATIVE_TEST_OPTIONS)
    option_file = output / "ipopt.opt"
    option_file.write_text("".join(f"{key} {v}\n" for key, v in options.items()))
    solver = get_solver(solver_options={"option_file_name": str(option_file)})
    solver.set_executable(str(ipopt))
    os.environ["AMPLFUNC"] = str(reference.library)
    started = time.monotonic()
    result = solver.solve(str(path), load_solutions=False, suffixes=["dual", "ipopt_zL_out", "ipopt_zU_out"])
    assert check_optimal_termination(result), str(result.solver)
    columns, rows = neutral.named_file_rows(path)
    for index, name in enumerate(columns):
        v = m.find_component(name)
        record = result.solution[0].variable[f"v{index}"]
        v.set_value(record["Value"])
        for suffix in ("ipopt_zL_out", "ipopt_zU_out"):
            getattr(m, suffix)[v] = record.get(suffix, 0.0)
    for index, name in enumerate(rows[:-1]):
        m.dual[m.find_component(name)] = result.solution[0].constraint[f"c{index}"].get("Dual", 0.0)
    text = (output / "ipopt.log").read_text()
    applied = {m[0]: m[1] for m in neutral.re.findall(r"^\s*(\w+)\s*=\s*(\S+)\s+yes\s*$", text, neutral.re.M)}
    for key, expected_option in options.items():
        assert key in applied and (
            float(applied[key]) == expected_option if isinstance(expected_option, (int, float)) else applied[key] == expected_option
        ), (key, applied)
    if optimize:
        assert "No errors detected by derivative checker." in text
        for quantity in ("Number of Lagrangian Hessian evaluations", "Number of equality constraint Jacobian evaluations"):
            assert any(quantity in line and int(line.rsplit(None, 1)[-1]) > 0 for line in text.splitlines())
    b = m.fs.settler
    y = np.array([value(v[i]) for v in (b.flow_mol_aqueous, b.flow_mol_organic) for i in fixed.COMPONENTS])
    metrics = {
        "si": value(m.si),
        "ss": value(m.ss),
        "objective": value(m.objective),
        "recovery": y[6] / F0[1],
        "material": max(abs(value(v)) for v in b.component_balance_residual.values()),
        "charge": max(abs(value(v)) for v in b.charge_residual.values()),
        "elapsed_s": time.monotonic() - started,
        "flows": y.tolist(),
        "termination": str(result.solver.termination_condition),
        "nl_sha256": fixed.sha256(path),
        "options": options,
        "option_file_sha256": fixed.sha256(option_file),
        **nl_checks(m, path, reference, optimize),
    }
    assert metrics["material"] <= 1e-10 and metrics["charge"] <= 1e-10 and metrics["recovery"] >= 0.965 - 1e-8
    return metrics


def derivatives(reference, points, output):
    transport, ladders, endpoints, certificates = [], [], [], {}
    directions = [*np.eye(4), *[np.array(d) / np.sqrt(2) for d in ((1, 1, 0, 0), (1, -1, 0, 0), (0, 0, 1, 1), (0, 0, 1, -1))]]
    for label, F in points.items():
        public, certificates[label] = reference.public(F)
        actual = reference.asl(F)
        for order, (expected, observed) in enumerate(zip(public, actual)):
            for index in np.ndindex(expected.shape):
                error = abs(float(expected[index] - observed[index]))
                limit = 1e-11 + 1e-9 * max(abs(float(expected[index])), abs(float(observed[index])))
                transport.append(
                    {
                        "case": label,
                        "order": order,
                        "index": str(index),
                        "public": float(expected[index]),
                        "asl": float(observed[index]),
                        "error": error,
                        "limit": limit,
                    }
                )
                assert error <= limit, transport[-1]
        if label not in ("anchor", "optimum"):
            continue
        for number, direction in enumerate(directions):
            dF = F * direction
            estimates, hessians, admitted = [], [], []
            for step in neutral.STEPS:
                pair = []
                for sign in (1, -1):
                    try:
                        values_, certificate = reference.public(F + sign * step * dF, 1)
                        pair.append(values_)
                    except RuntimeError as error:
                        certificate = {"error": str(error)}
                    endpoints.append(
                        {"case": label, "direction": number, "step": step, "sign": sign, "certificate": json.dumps(certificate)}
                    )
                if len(pair) != 2:
                    continue
                plus, minus = pair
                estimates.append((plus[0] - minus[0]) / (2 * step))
                hessians.append((plus[1] - minus[1]) / (2 * step))
                admitted.append(step)
            assert len(admitted) >= 3, (label, number, admitted)
            for row in range(10):
                for order, values_, expected, absolute, relative in (
                    (1, [v[row : row + 1] for v in estimates], np.atleast_1d(public[1][row] @ dF), 1e-8, 1e-6),
                    (2, [v[row] for v in hessians], public[2][row] @ dF, 1e-7, 1e-4),
                ):
                    check = neutral.ladder_summary(values_, expected, tuple(admitted), absolute, relative)
                    for k, step in enumerate(admitted):
                        ladders.append(
                            {
                                "case": label,
                                "row": row,
                                "direction": number,
                                "order": order,
                                "step": step,
                                "error": check["errors"][k],
                                "limit": check["targets"][k],
                                "mode": check["mode"],
                                "reference": json.dumps(expected.tolist()),
                                "estimate": json.dumps(values_[k].tolist()),
                                "scale": max(1.0, float(np.max(abs(expected))), float(np.max(abs(values_[k])))),
                                "scaled_error": check["errors"][k] / max(1.0, float(np.max(abs(expected))), float(np.max(abs(values_[k])))),
                                "adjacent_support": check["adjacent_support"],
                                "floor_support": check["floor_support"],
                                "floor": check["floor"],
                                "floor_target": check["floor_target"],
                            }
                        )
    for name, rows in (("transport", transport), ("ladders", ladders), ("endpoints", endpoints)):
        fixed.atomic_csv(output / f"{name}.csv", rows, tuple(rows[0]))
    return certificates


def failure_probe(name):
    r = Reference()
    args = [*F0, r.record, "0", "0"]
    accepted = r.public(F0)[0]
    if name == "determinism":
        assert np.allclose(r.asl(F0, 0)[0], accepted[0], atol=1e-11, rtol=1e-9)
        r.asl(F0)
        r.asl(F0 * (1.0, 1.0, 0.95, 0.95))
        assert all(np.allclose(a, b, atol=1e-11, rtol=1e-9) for a, b in zip(accepted, r.asl(F0, rows=reversed(range(10)))))
        print('PROBE_JSON={"case":"determinism","cold_warm_reordered_fgh":true}')
        return
    cold = name in ("version", "truncated", "feed_map", "phase_swap", "support", "boundary", "unavailable")
    if not cold:
        r.asl(F0)
    if name == "version":
        args[4] = r.record.replace("equilibrium.v3", "equilibrium.v2")
    elif name == "truncated":
        args[4] = r.record[:-1]
    elif name == "feed_map":
        assert "11:feed_offset,1:5,6:0x0p+0," in r.record
        args[4] = r.record.replace("11:feed_offset,1:5,6:0x0p+0,", "11:feed_offset,1:5,6:0x1p+0,", 1)
    elif name == "record_change":
        args[4] = r.record + "0:,"
    elif name == "selector":
        args[5] = "2"
    elif name == "negative":
        args[0] = -1.0
    elif name == "nonfinite":
        args[0] = float("nan")
    elif name in ("boundary", "unavailable"):
        F = F0.copy()
        if name == "boundary":
            # Ethanol and isobutanol 0.001 mol: the single liquid has positive TPD against salt-free and tested electrolyte
            # trials (#461; ionic trials not certified), so the isobutanol-rich phase vanishes.
            F[1:3] = 0.001
        else:
            # Salt only: the solve returns "rank-deficient declared equations" without a state (#461). The cause, two liquids
            # of pure NaCl with undetermined amounts, is inferred, not inspected.
            F[:3] = 0.0
        problem = fixed.engine_problem(feed(F))
        central = eq.solve_equilibrium(r.model, problem)
        if name == "unavailable":
            batch = eq.solved_state_actions(
                eq.compile_problem(r.model, problem), central, eq.SolvedStateActionRequest(r.observables, [r.directions[0]])
            )
            assert not central.success and not central.requested_tolerance_met and all(
                item.status == eq.SolvedStateActionStatus.CentralSolveUnavailable and item.action is None for item in batch.results
            )
        else:
            assert central.topology_event.name == "PhaseDisappeared"
        args = [*F, encode_record(r.model, problem, OPTIONS, F), "0", "0"]
    elif name == "support":
        support = "7:phase-1,6:liquid,1:5,1:0,1:1,1:2,1:3,1:4,1:1,"
        assert support in r.record
        args[4] = r.record.replace(support, support[:-4] + "1:0,", 1)
    elif name == "phase_swap":
        phases = list(r.problem.phases)
        phases.reverse()
        problem = eq.Problem(
            phases=phases,
            T=fixed.TEMPERATURE,
            P=fixed.PRESSURE,
            feed=eq.Amounts(dict(zip(fixed.COMPONENTS, feed(F0)))),
            row_tolerance=fixed.ROW_TOLERANCE,
        )
        args[4] = encode_record(r.model, problem, OPTIONS, F0)
    r.external.load_library()
    packed = _ARGLIST(args, fgh=2)
    result = r.external._known_functions[NAME][0](ctypes.byref(packed))
    message = packed.Errmsg.decode() if packed.Errmsg else ""
    failed = np.isnan(result) and all(np.isnan(packed.derivs[i]) for i in range(4)) and all(np.isnan(packed.hes[i]) for i in range(10))
    expected = {
        "negative": "outside_interval:phase[0].amount[0]",
        "nonfinite": "nonfinite:phase[0].amount[0]",
        "selector": "defining_equations_not_satisfied:request",
        "support": "defining_equations_not_satisfied:phase[0].amount[0]",
        "phase_swap": "phase_boundary_not_differentiable:phase[0].amount[0]",
        "boundary": "phase_boundary_not_differentiable:phase[0].amount[0]",
        "unavailable": "central_solve_unavailable:phase[0].amount[0]",
    }.get(name, "defining_equations_not_satisfied:phase[0].amount[0]")
    assert failed and message == '"epcsaft_equilibrium:' + expected, (name, result, message)
    assert all(np.allclose(a, b, atol=1e-11, rtol=1e-9) for a, b in zip(accepted, r.asl(F0)))
    print("PROBE_JSON=" + json.dumps({"case": name, "error": message, "nonfinite_outputs": failed, "cache_recovery_fgh": True}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=fixed.CAMPAIGN / "results/live")
    parser.add_argument("--ipopt", type=Path, default=Path("/home/tnnrpolley21/.idaes/bin/ipopt"))
    parser.add_argument("--probe")
    args = parser.parse_args()
    if args.probe:
        failure_probe(args.probe)
        return
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    r = Reference()
    _, _, library = neutral.installed_identity()
    solver = get_solver()
    solver.set_executable(str(args.ipopt))
    dependencies = neutral.shared_library_identity(library)
    outer = subprocess.check_output(["ldd", str(args.ipopt)], text=True)
    exports = subprocess.check_output(["nm", "-D", "--defined-only", str(args.ipopt)], text=True)
    assert not neutral.re.search(r"IpoptApplicationFactory|_Z\w*5Ipopt", exports), "outer Ipopt exports collide with Engine"
    assert all(name in dependencies["needed"] for name in neutral.re.findall(r"libipopt\.so\.\d+", outer)), "incompatible shared Ipopt ABI"
    identity = {
        **fixed.campaign_identity(solver),
        "library": str(library),
        "library_sha256": fixed.sha256(library),
        "dependencies": dependencies,
        "source_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=fixed.REPO, text=True).strip(),
        "source_dirty": bool(subprocess.check_output(["git", "status", "--porcelain"], cwd=fixed.REPO, text=True).strip()),
        "source_hashes": {
            str(p.resolve().relative_to(fixed.REPO)): fixed.sha256(p)
            for p in (
                fixed.REPO / "engine/native/bindings/asl_adapter.cpp",
                fixed.REPO / "engine/src/epcsaft/parameters.py",
                Path(__file__),
                fixed.REPO / "engine/native/bindings/equilibrium_record.hpp",
                fixed.REPO / "engine/native/bindings/equilibrium_bindings.cpp",
                fixed.REPO / "engine/src/epcsaft/equilibrium/__init__.py",
                Path(fixed.__file__),
                Path(neutral.__file__),
            )
        },
        "python": sys.version,
        "hardware": neutral.platform.uname()._asdict(),
        "threads": {k: os.environ.get(k) for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")},
        "record_sha256": __import__("hashlib").sha256(r.record.encode()).hexdigest(),
        "time_basis_s": 1.0,
        "writer_options": neutral.NL_IO_OPTIONS,
    }
    identity["solver"].pop("options")  # Per-stage option files own the live settings, not the fixed-feed helper.
    fixed.atomic_json(output / "identity.json", identity)
    (output / "record.netstrings").write_text(r.record, encoding="utf-8")
    summary = {"status": "running", "identity": identity}
    try:
        baseline = solve(r, (1.0, 1.0), output / "baseline", args.ipopt, False)
        with (fixed.CAMPAIGN / "results/fixed-feed/stage-results.csv").open() as stream:
            source = [row for row in __import__("csv").DictReader(stream) if row["point_id"] == "4"]
        expected = [
            float(row["phase_amount_mol"]) * float(row[name]) * sum(feed(F0))
            for phase in fixed.PHASES
            for row in source
            if row["phase"] == phase
            for name in ("x_water", "x_ethanol", "x_isobutanol", "x_sodium_cation", "x_chloride_anion")
        ]
        assert np.max(abs(np.array(baseline["flows"]) - expected)) <= 1e-10
        cases = [solve(r, (s, s), output / f"start-{s}", args.ipopt) for s in (0.9, 1.0, 1.1)]
        decisions = np.array([[c[k] for k in ("si", "ss", "objective")] for c in cases])
        assert np.max(np.ptp(decisions, axis=0)) <= 1e-7 and all(0.85 < c[k] < 1.15 for c in cases for k in ("si", "ss"))
        points = {
            "anchor": F0,
            **{f"start-{s}": F0 * (1.0, 1.0, s, s) for s in (0.9, 1.0, 1.1)},
            "optimum": F0 * (1.0, 1.0, cases[0]["si"], cases[0]["ss"]),
        }
        certificates = derivatives(r, points, output)
        failures = []
        for name in (
            "determinism",
            "version",
            "truncated",
            "feed_map",
            "record_change",
            "selector",
            "negative",
            "nonfinite",
            "phase_swap",
            "support",
            "boundary",
            "unavailable",
        ):
            probe = subprocess.run([sys.executable, __file__, "--probe", name], capture_output=True, text=True, timeout=60, check=False)
            (output / f"failure-{name}.log").write_text(probe.stdout + probe.stderr)
            assert probe.returncode == 0, (name, probe.stdout, probe.stderr)
            failures.append(json.loads(next(line[11:] for line in probe.stdout.splitlines() if line.startswith("PROBE_JSON="))))
        summary.update(
            status="accepted",
            baseline=baseline,
            optimization=cases,
            certificates=certificates,
            failures=failures,
            limits="Declared local two-liquid branch, fixed T/P and stated feed box; no discovery, global optimum, calorics or reactive Lithium qualification.",
        )
    except Exception as error:
        summary.update(status="failed", error=repr(error))
        raise
    finally:
        fixed.atomic_json(output / "summary.json", summary)


if __name__ == "__main__":
    main()
