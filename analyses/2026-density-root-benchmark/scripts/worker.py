"""Persistent worker in the selected checkout's installed Python environment."""

import json
import math
import os
import subprocess
import sys
import time
from pathlib import Path

from declarations import HERE, cross_check, fit_declaration, key_solves, load, models
from epcsaft import (
    PressureActionDirection,
    PropertyObservable,
    _core,
    provenance,
    trace,
)
from epcsaft import (
    equilibrium as eq,
)
from epcsaft import (
    regression as r,
)
from verification import application_check, density_check
from verification import density_values as state_values

DOCUMENT = json.loads((HERE / "inputs/cases.json").read_text())
CASES = {q["id"]: q for q in DOCUMENT["cases"]}
MODELS = models()
KEY_TOOL, KEYS = key_solves()
CROSS, INPUTS = cross_check()
REQUESTS = {q["id"]: q for q in INPUTS["requests"]}


def physical(result):
    if not result.success or result.requested_tolerance_met is not True:
        raise RuntimeError(result.message)
    values = {k: [getattr(result, k)] if k in ('temperature', 'pressure') else list(getattr(result, k)) for k in KEY_TOOL.PHYSICAL}
    assert all(v and all(map(math.isfinite, v)) for v in values.values())
    return values


def application(q):
    """Return one callable, with setup and accepted starts outside the timing boundary."""
    op = q["operation"]
    if op in ("C1", "C2", "bubble", "flash", "stability"):
        if op in ("C1", "C2"):
            model, problem, _, _ = CROSS.prepare(
                INPUTS, REQUESTS["C1-300" if op == "C1" else "C2-0.5"]
            )
        elif op == "bubble":
            model, problem = load(
                "analyses/2026-binary-methane-ethane-bubble-point/scripts/collect_results.py",
                "worker_bubble",
            ).collect_engine(declare_only=True)
        else:
            model, problem = load(
                "analyses/2026-binary-methane-ethane-tp-flash/scripts/performance.py",
                "worker_flash",
            ).declarations()
        if "T_K" in q:
            problem.T = q["T_K"]
        if op != "stability":
            return lambda: [physical(eq.solve_equilibrium(model, problem))]
        central = eq.solve_equilibrium(model, problem)
        physical(central)
        compiled = eq.compile_problem(model, problem)

        def stability():
            report = eq.stability_test(compiled, central)
            assert report.completed, report.message
            stability.native_work = {
                "solve_density_calls": int(report.work["root_calls"]),
                "pressure_evaluations": int(report.work["root_evaluations"]),
                "energy_evaluations": int(report.work["energy_evaluations"]),
            }
            return {
                "central": [physical(central)],
                "standard": report.standard_over_rt,
                "plane": report.plane_over_rt,
                "compositions": [t.composition for t in report.trials],
                "tpd_over_rt": [t.quantities.get("tpd_over_rt") for t in report.trials],
            }

        return stability
    if op == "key":
        name = q["id"]
        if name in KEY_TOOL.CASES[:4]:
            base, mode = name.rsplit("/", 1)
            model, problem = KEYS[base]
            start = None
            if mode == "accepted":
                start = eq.EquilibriumStart()
                start.coordinates = DOCUMENT["accepted_starts"][base]

            def ordinary():
                result = (
                    eq.solve_equilibrium(model, problem)
                    if start is None
                    else eq.solve_equilibrium_starts(
                        model, problem, [start], eq.EquilibriumWorkLimits()
                    ).results[0]
                )
                if result is None:
                    raise RuntimeError("accepted start returned no state")
                return [physical(result)]

            return ordinary
        calls = KEYS[name] if name == "khudaida/sequential-8" else [KEYS[name]]

        def sequence():
            previous = None
            rows = []
            for model, problem in calls:
                previous = eq.solve_equilibrium(model, problem, previous)
                rows.append(physical(previous))
            return rows

        return sequence
    if op.startswith("actions-"):
        cid = op.removeprefix("actions-")
        density_case = CASES["A-C2-liquid" if cid == "C2" else "B-C3-1-mol-kg"]

        def actions():
            batch = MODELS[cid].pressure_actions(
                q.get("T_K", density_case["T_K"]),
                density_case["P_Pa"],
                x=density_case["x"],
                phase="liquid",
                observables=[PropertyObservable.LogFugacityCoefficient],
                directions=[
                    PressureActionDirection(
                        pressure=1, composition=[0] * len(density_case["x"])
                    )
                ],
            )
            assert batch.density_status == _core.PropertyActionStatus.Available, (
                batch.density_message
            )
            assert all(
                item.status == _core.PropertyActionStatus.Available
                for item in batch.results
            )
            return {
                "density": [batch.molar_density],
                "values": [item.value for item in batch.results],
                "actions": [item.action for item in batch.results],
            }

        return actions
    if op in ("objective", "fit"):
        parameters, coordinates, observations, starts, controls = fit_declaration()
        controls.workers = 1
        chart = r.CoordinateSet()
        chart.coordinates = coordinates
        optimizer = chart.start(parameters._native, observations)
        if op == "objective":

            def objective():
                result = _core.evaluate_observations(
                    parameters._native, coordinates, optimizer, [], observations
                )
                assert all(s == r.ObservationStatus.Available for s in result.statuses)
                return {
                    "predictions": result.predictions,
                    "prediction_gradient": result.gradient,
                    "residual_optimizer_gradient": result.residual_optimizer_gradient,
                }

            return objective

        def fit():
            results = r.fit(
                parameters,
                coordinates,
                observations,
                starts=starts[:1],
                controls=controls,
                shared_iterations=controls.maximum_iterations,
                shared_seconds=controls.maximum_elapsed_time_seconds,
            )
            assert len(results) == 1
            assert all(result.usable for result in results), [
                (v.status, v.message) for v in results
            ]
            return [
                {
                    "physical": result.physical,
                    "weighted_residuals": result.weighted_residuals,
                    "cost": [result.final_cost],
                }
                for result in results
            ]

        return fit
    raise ValueError(op)


def work(captured):
    # Ancestor-dependent density counters are diagnostics; total evaluate spans are authoritative.
    roots = sum(s.name == "solve_density" for s in captured.spans)
    evaluations = 0
    for span in captured.spans:
        if span.name != "evaluate":
            continue
        parent = span.parent
        while parent is not None:
            if captured.spans[parent].name == "solve_density":
                evaluations += 1
                break
            parent = captured.spans[parent].parent
    return {
        "solve_density_calls": roots,
        "pressure_evaluations": evaluations,
        "eos_evaluations": sum(s.name == "evaluate" for s in captured.spans),
        "energy_evaluations": 0,
    }


def density_values(result, fixed=False):
    values = state_values(result, fixed)
    if not fixed:
        values.update(
            reported_criterion=result.density_convergence.criterion.value,
            reported_floor=result.density_convergence.floor or 0,
        )
    return values


def calls(q):
    """Fresh models/declarations outside both timing and work capture, with matching histories."""
    global MODELS, KEYS
    MODELS = models()
    if q["tier"] in ("A", "B"):
        model = MODELS[q["model"]]

        def pair(state):
            returned_density = [state["oracle_rho_mol_m3"]]

            def pressure():
                return model.state(
                    state["T_K"], P=state["P_Pa"], x=state["x"], phase=state["branch"]
                )

            def fixed():
                return model.state(state["T_K"], rho=returned_density[0], x=state["x"])

            fixed.set_density = lambda value: returned_density.__setitem__(0, value)
            return pressure, fixed

        return [pair(state) for state in q["states"]]
    KEYS = key_solves()[1]
    return [(application(q | state), None) for state in q["states"]]


def measure(q, repeated, references):
    targets = calls(q)
    indices = [0] * q["N"] if repeated else [i % len(targets) for i in range(q["N"])]
    records, errors = [], []
    density_case = q["tier"] in ("A", "B")

    def check(index, kind, value):
        try:
            if density_case:
                density_check(
                    q["states"][index],
                    value,
                    DOCUMENT["gas_constant_J_mol_K"],
                    fixed=kind == "fixed",
                )
            else:
                application_check(value, references[q["id"]][index])
        except (AssertionError, TypeError, ValueError) as error:
            errors.append(str(error))

    def invoke(index, kind, stage, timed=False):
        call = targets[index][kind == "fixed"]
        begun = time.perf_counter_ns() if timed else None
        result = call()
        elapsed = (time.perf_counter_ns() - begun) / 1000 if timed else None
        # Timer has stopped before any result extraction or numerical acceptance.
        value = (
            density_values(result, fixed=kind == "fixed") if density_case else result
        )
        if density_case and kind == "pressure":
            targets[index][1].set_density(value["rho_mol_m3"])
        check(index, kind, value)
        record = dict(state_index=index, kind=kind, stage=stage, values=value)
        if timed:
            record["us"] = elapsed
        records.append(record)

    kinds = ["pressure", "fixed"] if density_case else ["application"]
    for _ in range(q["warmups"]):
        for kind in kinds:
            invoke(0, kind, "warmup")
    for kind in kinds:
        for index in indices:
            invoke(index, kind, "timed", timed=True)

    # A fresh replay with the same warm-up and call history measures the timed mode's work.
    targets = calls(q)
    for _ in range(q["warmups"]):
        for kind in kinds:
            invoke(0, kind, "work-warmup")
    for kind in kinds:
        for index in indices:
            with trace(diagnostics=True) as captured:
                call = targets[index][kind == "fixed"]
                result = call()
            value = (
                density_values(result, fixed=kind == "fixed")
                if density_case
                else result
            )
            if density_case and kind == "pressure":
                targets[index][1].set_density(value["rho_mol_m3"])
            check(index, kind, value)
            counts = work(captured)
            native = getattr(call, "native_work", {})
            if native:
                counts.update(native)
            records.append(
                dict(
                    state_index=index,
                    kind=kind,
                    stage="work",
                    values=value,
                    counts=counts,
                )
            )
    return dict(
        mode="repeated identical input (rewards caching)" if repeated else q["mode"],
        records=records,
        errors=errors,
    )


def run(q, rep, references):
    primary = measure(q, False, references)
    secondary = measure(q, True, references) if len(q["states"]) > 1 else None
    return dict(
        replica=rep, pid=os.getpid(), primary=primary, repeated_identical=secondary
    )


def main():
    assert os.getpriority(os.PRIO_PROCESS, 0) == 10, "timing worker must have nice 10"
    if "--reference" in sys.argv:
        old = __import__("csv").DictReader((HERE / "results/main-0225e03b.csv").open())
        old_values = {
            row["id"]: json.loads(row["values"]) for row in old if row["tier"] == "C"
        }
        references = {}
        for q in DOCUMENT["cases"]:
            if q["tier"] == "C":
                print("application reference", q["id"], flush=True)
                references[q["id"]] = [call() for call, _ in calls(q)]
                application_check(references[q["id"]][0], old_values[q["id"]])
        (HERE / "results/application-reference.json").write_text(
            json.dumps(references, indent=2) + "\n"
        )
        return
    references = json.loads((HERE / "results/application-reference.json").read_text())
    reply = os.fdopen(os.dup(1), "w")
    os.dup2(2, 1)
    sys.path.insert(0, str(HERE.parents[1] / "tools"))
    from engine_identity import loaded_identity
    loaded = loaded_identity()
    identity, detail = loaded["native_library_sha256"], loaded
    source_tree = subprocess.check_output(
        ["git", "rev-parse", f"{loaded['source_revision']}:engine"],
        cwd=loaded["source_root"], text=True).strip()
    print(
        json.dumps(
            dict(
                engine_identity=str(identity),
                engine_detail=detail,
                engine_tree=source_tree,
                nice=os.getpriority(os.PRIO_PROCESS, 0),
                pid=os.getpid(),
            )
        ),
        file=reply,
        flush=True,
    )
    for line in sys.stdin:
        request = json.loads(line)
        try:
            answer = run(CASES[request["id"]], request["rep"], references)
        except Exception as error:  # noqa: BLE001 -- preserve native failures on the response pipe
            answer = {"error_type": type(error).__name__, "error": str(error)}
        print(json.dumps(answer, allow_nan=False), file=reply, flush=True)


if __name__ == "__main__":
    main()
