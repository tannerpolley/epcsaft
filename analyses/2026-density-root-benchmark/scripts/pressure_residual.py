"""Native nested density scans and direct pressure-residual verification (#686)."""

import argparse
import csv
import hashlib
import json
import math
import os
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import prepare
from declarations import HERE, models
from epcsaft import _core, provenance
from verification import (
    DENSITY_RELATIVE_LIMIT,
    PRESSURE_RESIDUAL_LIMIT,
    density_check,
    density_mutations,
    density_values,
)

OUTPUT = HERE / "results/main-8298c602-pressure-residual"
WHEEL_SHA256 = "8897a0caab9f2508ed642f42c8b5226957fbe7035ed7bfdad71e4166798785b3"


def initialize():
    global MODELS
    MODELS = models()


def calculate(job):
    case, state = job
    model = MODELS[case["model"]]
    # Scanned references belong to the run declaration, not the frozen state.
    state = {k: v for k, v in state.items() if k != "scanned_references"}
    record = dict(case=case["id"], state=state, scans=[], returns=[], failures=[])
    T, P, x = state["T_K"], state["P_Pa"], state["x"]
    for intervals in (prepare.INTERVALS, 2 * prepare.INTERVALS):
        scan = dict(intervals=intervals, roots=[], spinodals=[], unavailable_samples=[])
        record["scans"].append(scan)
        try:
            native = []
            prepare.enumerate_roots(model, state, intervals=intervals, native_result=native)
            result, lo, hi, tolerance = native[0]
            scan.update(
                lower=lo,
                upper=hi,
                density_tolerance=tolerance,
                log_tolerance=prepare.LOG_TOLERANCE,
                max_evaluations=3 * intervals,
                max_seconds=120,
                stop=result.stop,
                evaluations=result.evaluations,
                sample_count=len(result.samples),
                unavailable_samples=[prepare.diagnostic_point(p) for p in result.samples if p.cause],
                unresolved=[prepare.diagnostic_point(p) for p in result.unresolved],
            )
            for name, points in (("roots", result.pressure_roots), ("spinodals", result.spinodals)):
                for point in points:
                    row = dict(
                        prepare.diagnostic_point(point),
                        endpoints=[],
                        bracket_log_width=math.log(point.upper / point.lower),
                    )
                    scan[name].append(row)
                    endpoints = model.density_diagnostics(
                        T,
                        P,
                        x=x,
                        nodes=[
                            point.lower,
                            point.upper if point.upper > point.lower else math.nextafter(point.lower, math.inf),
                        ],
                        density_tolerance=tolerance,
                        log_tolerance=prepare.LOG_TOLERANCE,
                        max_evaluations=64,
                        max_seconds=120,
                        adaptive=False,
                        finite_upper=hi,
                        evaluation_order=_core.DerivativeOrder.PressureSlope,
                    )
                    row.update(
                        endpoints=[prepare.diagnostic_point(p) for p in endpoints.samples], endpoint_stop=endpoints.stop
                    )
        except Exception as error:
            record["failures"].append(
                dict(operation=f"scan-{intervals}", cause=str(error), exception=type(error).__name__)
            )
    returned = None
    for mode in ("pressure", "returned-density", "frozen-density"):
        try:
            if mode == "returned-density" and returned is None:
                raise RuntimeError("pressure request unavailable; returned density does not exist")
            args = (
                dict(P=P, phase=state["branch"])
                if mode == "pressure"
                else dict(rho=returned if mode == "returned-density" else state["oracle_rho_mol_m3"])
            )
            value = density_values(model.state(T, x=x, **args), fixed=mode != "pressure")
            if mode == "pressure":
                returned = value["rho_mol_m3"]
            record["returns"].append(dict(mode=mode, **value))
        except Exception as error:
            record["failures"].append(dict(operation=mode, cause=str(error), exception=type(error).__name__))
    return record


def check(records, document):
    frozen = json.loads((HERE / "inputs/cases.json").read_text())
    cases = [c for c in frozen["cases"] if c["tier"] != "C"]
    assert len(cases) == len(document["cases"]) == 22
    old = list(csv.DictReader((HERE / "results/oracle.csv").open()))
    residuals, frozen_deltas, scan_deltas, nested, brackets, reproduction = [], [], [], [], [], []
    jobs = [(c, s) for c in document["cases"] for s in c["states"]]
    assert len(records) == len(jobs) == 704
    assert document["pressure_criterion"] == dict(
        reduced_residual_limit=PRESSURE_RESIDUAL_LIMIT, density_relative_limit=DENSITY_RELATIVE_LIMIT
    )
    for case, before in zip(document["cases"], cases, strict=True):
        assert {k: v for k, v in case.items() if k != "states"} == {k: v for k, v in before.items() if k != "states"}
        assert len(case["states"]) == len(before["states"]) == len({s["P_Pa"] for s in case["states"]}) == 32
        for state, original in zip(case["states"], before["states"], strict=True):
            assert {k: v for k, v in state.items() if k != "scanned_references"} == {
                k: v for k, v in original.items() if k != "main_pressure_floor"
            }
    for (case, state), record in zip(jobs, records, strict=True):
        identity = (case["id"], state["index"])
        assert record["case"] == case["id"] and record["state"] == {
            k: v for k, v in state.items() if k != "scanned_references"
        }, identity
        assert not record["failures"], (identity, record["failures"])
        assert [s["intervals"] for s in record["scans"]] == [prepare.INTERVALS, 2 * prepare.INTERVALS], identity
        for scan in record["scans"]:
            assert (
                scan["stop"] == "resolved_crossings" and not scan["unavailable_samples"] and not scan["unresolved"]
            ), identity
            assert scan["sample_count"] == scan["intervals"] + 1, identity
            for name in ("roots", "spinodals"):
                for point in scan[name]:
                    endpoints = point["endpoints"]
                    operands = [
                        p["pressure"] - state["P_Pa"] if name == "roots" else p["pressure_density_derivative"]
                        for p in endpoints
                    ]
                    assert len(endpoints) == 2 and not any(p["cause"] for p in endpoints), identity
                    assert all(math.isfinite(v) for v in operands) and min(operands) <= 0 <= max(operands), identity
                    assert all(
                        math.isfinite(point[k])
                        for k in ("molar_density", "pressure", "pressure_density_derivative", "residual")
                    ), identity
                    assert point["bracket_log_width"] == math.log(point["upper"] / point["lower"]), identity
                    if point["stop"] == "bracket_tolerance":
                        assert (
                            point["bracket_log_width"] <= scan["log_tolerance"]
                            and point["upper"] - point["lower"] <= scan["density_tolerance"]
                        ), identity
                    else:
                        assert point["stop"] == "exact_numerical_zero" and point["residual"] == 0 and 0 in operands, (
                            identity
                        )
                    brackets.append((name, point))
            stable = [p for p in scan["roots"] if p["mechanical_sign"] > 0]
            assert len(stable) == state["oracle_stable_root_count"], identity
            scan["selected_rho_mol_m3"] = stable[0 if state["branch"] == "vapor" else -1]["molar_density"]
        coarse, fine = record["scans"]
        if state["index"] == 0:
            previous = [r for r in old if r["id"] == case["id"] and int(r["state_index"]) == 0]
            for before, density, delta in prepare.compare_roots(
                previous, [p["molar_density"] for p in coarse["roots"]]
            ):
                reproduction.append(
                    dict(
                        case=case["id"],
                        root=before["root"],
                        previous_density=float(before["rho_mol_m3"]),
                        density=density,
                        relative_difference=delta,
                    )
                )
            if case["id"] == "A-propane-vapor":
                assert len(coarse["roots"]) == 3 and len(coarse["spinodals"]) == 2, identity
        for name in ("roots", "spinodals"):
            for a, b in zip(coarse[name], fine[name], strict=True):
                delta = abs(a["molar_density"] / b["molar_density"] - 1)
                assert delta <= DENSITY_RELATIVE_LIMIT and a["mechanical_sign"] == b["mechanical_sign"], identity
                if name == "spinodals":
                    slopes = [[e["pressure_density_derivative"] for e in p["endpoints"]] for p in (a, b)]
                    assert (slopes[0][-1] > slopes[0][0]) == (slopes[1][-1] > slopes[1][0]), identity
                a["nested_relative_difference"] = delta
                nested.append(delta)
        references = [{k: s[k] for k in ("intervals", "selected_rho_mol_m3")} for s in record["scans"]]
        assert state["scanned_references"] == references, identity
        assert [v["mode"] for v in record["returns"]] == ["pressure", "returned-density", "frozen-density"], identity
        assert record["returns"][1]["rho_mol_m3"] == record["returns"][0]["rho_mol_m3"], identity
        assert record["returns"][2]["rho_mol_m3"] == state["oracle_rho_mol_m3"], identity
        for value in record["returns"]:
            frozen_deltas.append(
                density_check(state, value, document["gas_constant_J_mol_K"], fixed=value["mode"] != "pressure")
            )
            delta = max(abs(value["rho_mol_m3"] / s["selected_rho_mol_m3"] - 1) for s in record["scans"])
            assert delta <= DENSITY_RELATIVE_LIMIT, identity
            scan_deltas.append(delta)
            residuals.append(
                abs(value["pressure_Pa"] - state["P_Pa"])
                / (value["rho_mol_m3"] * document["gas_constant_J_mol_K"] * state["T_K"])
            )
    rejected = 0
    for field, changed in density_mutations(records[0]["returns"][0]):
        try:
            density_check(document["cases"][0]["states"][0], changed, document["gas_constant_J_mol_K"])
        except AssertionError:
            rejected += 1
            continue
        raise AssertionError(f"accepted mutation: {field}")
    assert len(reproduction) == 56
    summary = dict(
        requests=len(jobs),
        returns=len(residuals),
        oracle_comparisons=len(reproduction),
        scans=2 * len(jobs),
        failures=[],
        mutation_rejections=rejected,
        pressure_limit=PRESSURE_RESIDUAL_LIMIT,
        density_limit=DENSITY_RELATIVE_LIMIT,
        max_reduced_pressure_residual=max(residuals),
        max_frozen_density_difference=max(frozen_deltas),
        max_scan_density_difference=max(scan_deltas),
        max_oracle_density_difference=max(p["relative_difference"] for p in reproduction),
        max_nested_density_difference=max(nested),
        max_bracket_log_width=max(p["bracket_log_width"] for _, p in brackets if p["stop"] == "bracket_tolerance"),
        max_exact_zero_log_width=max(
            p["bracket_log_width"] for _, p in brackets if p["stop"] == "exact_numerical_zero"
        ),
        exact_numerical_zeros=sum(p["stop"] == "exact_numerical_zero" for _, p in brackets),
        pressure_crossings=sum(name == "roots" for name, _ in brackets),
        spinodals=sum(name == "spinodals" for name, _ in brackets),
    )
    return summary, reproduction


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="check retained values without calculating or writing")
    args = parser.parse_args()
    document = json.loads((OUTPUT / "cases.json").read_text())
    identity = json.loads((OUTPUT / "provenance.json").read_text())
    for name, key in (("inputs/cases.json", "input_sha256"), ("results/oracle.csv", "frozen_oracle_sha256")):
        assert hashlib.sha256((HERE / name).read_bytes()).hexdigest() == identity[key]
    paths = [OUTPUT / f"{c['id']}-state-{s['index']}.json" for c in document["cases"] for s in c["states"]]
    if args.check:
        records = [json.loads(p.read_text()) for p in paths]
    else:
        assert os.getpriority(os.PRIO_PROCESS, 0) == 10
        wheel, digest = provenance()
        assert digest == WHEEL_SHA256 == hashlib.sha256(Path(wheel).read_bytes()).hexdigest()
        assert [str(wheel), digest] == document["actual_provenance"] == identity["actual_provenance"]
        jobs = [(c, s) for c in document["cases"] for s in c["states"]]
        with ProcessPoolExecutor(max_workers=os.cpu_count(), initializer=initialize) as pool:
            records = []
            for record in pool.map(calculate, jobs, chunksize=1):
                records.append(record)
                print(record["case"], record["state"]["index"], record["failures"], flush=True)
    # Compare the complete records, including binary64 signs, rather than absorb any change.
    summary, reproduction = check(records, document)
    for path, record in zip(paths, records, strict=True):
        assert json.dumps(record, sort_keys=True) == json.dumps(json.loads(path.read_text()), sort_keys=True), path
    expected = json.loads((OUTPUT / "summary.json").read_text())
    assert summary == {k: v for k, v in expected.items() if k != "retained_files"}
    assert list(csv.DictReader((OUTPUT / "oracle-reproduction.csv").open())) == [
        {k: str(v) for k, v in row.items()} for row in reproduction
    ]
    print(json.dumps(summary), flush=True)
    print(f"{len(records)} canonical records checked; retained evidence unchanged", flush=True)


if __name__ == "__main__":
    main()
