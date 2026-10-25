#!/usr/bin/env python3
"""Retained numerical evidence and the independent review's result-mutation checks."""

import ast
import copy
import csv
import hashlib
import json
import statistics
from pathlib import Path
from types import SimpleNamespace

from run import validate
from verification import density_mutations

HERE = Path(__file__).resolve().parents[1]


def main():
    document = json.loads((HERE / "inputs/cases.json").read_text())
    references = json.loads((HERE / "results/application-reference.json").read_text())
    samples = json.loads((HERE / "results/main-samples.json").read_text())
    cases = {q["id"]: q for q in document["cases"]}
    checked = 0
    for id, records in samples["cases"].items():
        ok, _, n, errors = validate(cases[id], records, document, references)
        assert ok, (id, errors)
        checked += n
        for replica in records:
            for mode in (replica["primary"], replica["repeated_identical"]):
                if mode is None:
                    continue
                kinds = (
                    ("pressure", "fixed")
                    if cases[id]["tier"] != "C"
                    else ("application",)
                )
                for kind in kinds:
                    timed = [
                        r
                        for r in mode["records"]
                        if r["stage"] == "timed" and r["kind"] == kind
                    ]
                    work = [
                        r
                        for r in mode["records"]
                        if r["stage"] == "work" and r["kind"] == kind
                    ]
                    assert len(timed) == len(work) == cases[id]["N"]
                    assert [r["state_index"] for r in timed] == [
                        r["state_index"] for r in work
                    ]
                    if kind == "fixed":
                        pressure = [
                            r
                            for r in mode["records"]
                            if r["stage"] == "timed" and r["kind"] == "pressure"
                        ]
                        assert all(
                            p["values"]["rho_mol_m3"] == f["values"]["rho_mol_m3"]
                            for p, f in zip(pressure, timed, strict=True)
                        )

                    if mode is replica["primary"] and len(cases[id]["states"]) > 1:
                        assert len({r["state_index"] for r in timed}) == 32

    id = "A-propane-vapor"
    control = next(
        r["values"]
        for r in samples["cases"][id][0]["primary"]["records"]
        if r["stage"] == "timed" and r["kind"] == "pressure"
    )
    for field, mutation in density_mutations(control):
        changed = copy.deepcopy(samples["cases"][id])
        # Mutate only a timed return, leaving all warm-ups and traced returns correct.
        result = next(
            r["values"] for r in changed[0]["primary"]["records"] if r["stage"] == "timed" and r["kind"] == "pressure"
        )
        result.update(mutation)
        result["reported_criterion"] = result["reported_floor"] = 100
        assert not validate(cases[id], changed, document, references)[0], field
    changed = copy.deepcopy(samples["cases"]["C1-saturation"])
    result = next(
        r["values"] for r in changed[0]["primary"]["records"] if r["stage"] == "timed"
    )
    result[0]["temperature"][0] *= 2
    assert not validate(cases["C1-saturation"], changed, document, references)[0]

    # Execute the actual counting function without importing or calculating with the Engine.
    tree = ast.parse((HERE / "scripts/worker.py").read_text())
    function = next(
        node
        for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name == "work"
    )
    namespace = {}
    exec(
        compile(ast.Module(body=[function], type_ignores=[]), "worker.py:work", "exec"),
        namespace,
    )
    spans = [
        SimpleNamespace(name="solve_density", parent=None),
        SimpleNamespace(name="evaluate", parent=0),
        SimpleNamespace(name="evaluate", parent=None),
    ]
    before = namespace["work"](SimpleNamespace(spans=spans))
    spans[0].name = "renamed_density_solver"
    after = namespace["work"](SimpleNamespace(spans=spans))
    assert before["eos_evaluations"] == after["eos_evaluations"] == 2
    assert before["pressure_evaluations"] == 1 and after["pressure_evaluations"] == 0
    roots = list(csv.DictReader((HERE / "results/oracle-reproduction.csv").open()))
    assert (
        len(roots) == 56
        and max(float(r["relative_difference"]) for r in roots) <= 1e-10
    )
    for q in cases.values():
        if q["tier"] != "C":
            assert len({s["P_Pa"] for s in q["states"]}) == 32
            assert (
                len(
                    {
                        (s["oracle_branch"], s["oracle_stable_root_count"])
                        for s in q["states"]
                    }
                )
                == 1
            )
    competitors = json.loads((HERE / "results/competitor-samples.json").read_text())
    assert (
        competitors["protocol"]["cases_sha256"]
        == hashlib.sha256((HERE / "inputs/cases.json").read_bytes()).hexdigest()
    )
    identities = competitors["processes"]
    assert len({p["pid"] for p in identities}) == 6
    assert all(p["nice"] == 10 and p["threadpools"] == 1 for p in identities)
    summaries = {
        (r["id"], r["library"]): r
        for r in csv.DictReader((HERE / "results/competitors.csv").open())
    }
    competitor_returns = 0
    for row in competitors["cases"]:
        timed = row["samples"]
        if not timed:
            assert row["status"] == "unrepresentable"
            continue
        assert len(timed) == 192
        primary = [s for s in timed if s["mode"] == "unique states (one pass)"]
        assert len(primary) == 96
        for replica in range(3):
            assert {
                s["state_index"] for s in primary if s["replica"] == replica
            } == set(range(32))
        deltas = []
        for s in timed:
            if s["rho_mol_m3"] is None:
                assert s["status"] == "unrepresentable"
                continue
            state = cases[row["id"]]["states"][s["state_index"]]
            assert s["oracle_rho_mol_m3"] == state["oracle_rho_mol_m3"]
            delta = abs(s["rho_mol_m3"] / state["oracle_rho_mol_m3"] - 1)
            assert delta == s["relative_density_difference"]
            assert s["accuracy_within_1e10"] == (delta <= 1e-10)
            deltas.append(delta)
            competitor_returns += 1
        summary = summaries[row["id"], row["library"]]
        assert float(summary["median_us"]) == statistics.median(
            s["elapsed_us"] for s in primary
        )
        if deltas:
            assert float(summary["max_relative_density_difference"]) == max(deltas)
            assert (summary["accuracy_within_1e10"] == "true") == (max(deltas) <= 1e-10)
    print(
        f"{competitor_returns} competitor densities verified on the frozen grid in six processes"
    )
    print(
        f"{checked} retained returns verified; timed density/slope/branch/count/pressure and application mutations rejected; ancestor rename preserves total EOS count; 56 original roots reproduced"
    )


if __name__ == "__main__":
    main()
