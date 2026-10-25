"""Held-out comparison (#575): the frozen #556 compare.py, by path, on every stage-A optimum within the 10% cost band.

Each run pairs the unfitted held-2014-strategy2/2 (first slot) with one optimum (second slot), so the #556 transfer
table reads as the optimum's improvement over the settled transcription. Complete #556 outputs stay under build/;
the leaf retains verdicts, nominal bucket metrics, reference checks and the lowest-cost optimum's nominal residuals.

Owner amendment (2026-10-16): each Dinane-Mounir bucket is judged, in every frozen realization, against the larger of
the frozen limit and the Rard-Miller/Dinane-Mounir source disagreement on the same rows and realization (amend below,
conversion in consistency.py). Every other criterion is compare.py's frozen verdict.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import consistency  # noqa: E402
import fit  # noqa: E402
from consistency import DINANE_MOUNIR, RARD_MILLER  # noqa: E402
from packets import fitted_document  # noqa: E402

HERE, ROOT = fit.HERE, fit.ROOT
UNFITTED = ROOT / "data/packets/held-2014-strategy2/2"
HELD_OUT = ("robinson-bower-1966", DINANE_MOUNIR)
SLOT = {"Held": "unfitted", "Figiel": "candidate"}
BUCKETS = ("I<=3", "I>3")
LIMITS = json.loads((ROOT / fit.DECLARATION["held_out"]["specification"]).read_text())["primary_limits"]


def read(path: Path) -> list[dict]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def write(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def agreement_limits(compare, scenarios: list[dict]) -> dict:
    """Per (realization, Dinane-Mounir bucket): the source disagreement and the amended limit. A realization perturbs
    only its own source, the other stays nominal, as in compare.py's effective realization."""
    frozen = LIMITS["aw"]
    nominal = {s: compare.realize(s, scenarios[0]) for s in (RARD_MILLER, DINANE_MOUNIR)}
    buckets = consistency.nominal_buckets(nominal[DINANE_MOUNIR])
    limits = {}
    for scenario in scenarios:
        rows = {s: compare.realize(s, scenario) if scenario["source"] in ("all", s) else nominal[s] for s in nominal}
        differences = consistency.implied(rows[RARD_MILLER], rows[DINANE_MOUNIR], buckets, compare.MW_SOURCE)
        for bucket in BUCKETS:
            d = [o["difference"] for o in differences if o["I_bucket"] == bucket]
            rms, largest = math.sqrt(sum(x * x for x in d) / len(d)), max(map(abs, d))
            limits[scenario["realization_id"], bucket] = {
                "realization_id": scenario["realization_id"], "I_bucket": bucket, "rows": len(d),
                "disagreement_RMS": rms, "disagreement_maximum": largest, "limit_RMS": max(frozen["RMS"], rms),
                "limit_maximum": max(frozen["maximum_absolute_error"], largest)}
    return limits


def amend(results: Path, scenarios: list[dict], limits: dict) -> tuple[list[dict], dict[str, set], list[dict]]:
    """compare.py's verdict rows, Dinane-Mounir buckets re-judged per realization against the amended limits, each
    slot's realizations failing any held-out bucket, and one detail row per failing Dinane-Mounir realization (the leaf
    keeps the candidate's). A verdict row lists its failing realizations when there are at most 20."""
    metrics = {(r["realization_id"], r["source"], r["I_bucket"], r["declaration"]): r for r in read(results / "bucket-metrics.csv")}
    coverage = read(results / "realization-coverage.csv")
    rows, failed, details = [], {"Held": set(), "Figiel": set()}, []
    for v in read(results / "verdicts.csv"):
        frozen = LIMITS[v["quantity"]]
        limit = (frozen["RMS"], frozen["maximum_absolute_error"])
        nominal_passed, passing, robust = v["nominal_passed"], v["passing_realizations"], v["robust_passed"]
        failing = [r["realization_id"] for r in coverage if r["source"] == v["source"] and r["I_bucket"] == v["I_bucket"]
                   and r[v["declaration"] + "_passed"] != "True"]
        if v["source"] == DINANE_MOUNIR:
            failing = []
            for scenario in scenarios:
                effective = scenario["realization_id"] if scenario["source"] in ("all", DINANE_MOUNIR) else "r0000"
                m = metrics[effective, DINANE_MOUNIR, v["I_bucket"], v["declaration"]]
                bound = limits[scenario["realization_id"], v["I_bucket"]]
                if not (int(m["native_failure_count"]) == 0 and m["RMS"] != "" and float(m["RMS"]) <= bound["limit_RMS"]
                        and float(m["maximum_absolute_error"]) <= bound["limit_maximum"]):
                    failing.append(scenario["realization_id"])
                    details.append({"declaration": v["declaration"], "I_bucket": v["I_bucket"], "effective_source_realization": effective,
                                    **{k: m[k] for k in ("RMS", "maximum_absolute_error", "native_failure_count")}, **bound})
            limit = (limits["r0000", v["I_bucket"]]["limit_RMS"], limits["r0000", v["I_bucket"]]["limit_maximum"])
            nominal_passed, passing, robust = "r0000" not in failing, len(scenarios) - len(failing), not failing
        if v["source"] in HELD_OUT:
            failed[v["declaration"]].update(failing)
        rows.append({**v, "limit_RMS": limit[0], "limit_maximum": limit[1], "nominal_passed": nominal_passed,
                     "passing_realizations": passing, "robust_passed": robust,
                     "frozen_passing_realizations": v["passing_realizations"], "frozen_robust_passed": v["robust_passed"],
                     "failing_realizations": ";".join(failing) if len(failing) <= 20 else f"{len(failing)} realizations"})
    return rows, failed, details


def evaluate(compare, results: Path, scratch: Path, document: dict, scenarios: list[dict], limits: dict) -> tuple[list[dict], set]:
    shutil.rmtree(scratch, ignore_errors=True)
    (scratch / "packet/parameters").mkdir(parents=True)
    (scratch / "packet/parameters/parameters.json").write_text(json.dumps(document, indent=2))
    compare.PACKETS = {"Held": str(UNFITTED), "Figiel": str(scratch / "packet")}  # compare.py names its two slots
    compare.RESULTS = scratch / "results"
    sys.argv = ["compare.py", "--workers", str(os.cpu_count())]
    compare.main()
    kept = results / "held-out" / scratch.name
    kept.mkdir(parents=True, exist_ok=True)
    verdicts, failed, details = amend(compare.RESULTS, scenarios, limits)
    for name, rows in (("verdicts.csv", verdicts), ("dinane-mounir-failures.csv", [d for d in details if d["declaration"] == "Figiel"]),
                       ("nominal-buckets.csv", [r for r in read(compare.RESULTS / "bucket-metrics.csv") if r["realization_id"] == "r0000"]),
                       ("binary-reference-checks.csv", read(compare.RESULTS / "binary-reference-checks.csv")),
                       ("nominal-residuals.csv", read(compare.RESULTS / "nominal-residuals.csv"))):
        if not rows:
            continue
        write(scratch / name, [{**r, "declaration": SLOT[r["declaration"]]} for r in rows])
        if name != "nominal-residuals.csv":
            shutil.copy(scratch / name, kept / name)
    shutil.copy(compare.RESULTS / "summary.json", kept / "summary.json")
    return [r for r in read(kept / "verdicts.csv") if r["declaration"] == "candidate"], failed["Figiel"]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", default=fit.DECLARATION["run"], help="results/<run> holding fit-H.json and fit-F.json")
    run = parser.parse_args().run
    results, scratch = HERE / "results" / run, ROOT / "build/issue-575/held-out" / run
    compare = fit.comparison()
    scenarios = compare.read(compare.INPUTS / "frozen-realizations.csv")
    limits = agreement_limits(compare, scenarios)
    write(HERE / "results/source-agreement-limits.csv", list(limits.values()))
    selection = {"rule": fit.DECLARATION["optimum_selection"], "held_out_amendment": fit.DECLARATION["held_out"]["amendment"],
                 "forms": {}}
    for form, spec in fit.DECLARATION["forms"].items():
        record = json.loads((results / f"fit-{form}.json").read_text())
        labels = record["starts"][0]["labels"]
        verdicts = []
        for k, optimum in enumerate(o for o in record["optima"] if o["within_band"]):
            label = f"{form}-optimum-{k + 1}"
            rows, failing = evaluate(compare, results, scratch / label, fitted_document(spec, labels, optimum["physical"], optimum["starts"][0], label),
                            scenarios, limits)
            held = [r for r in rows if r["source"] in HELD_OUT]
            verdicts.append({"label": label, "starts": optimum["starts"], "final_cost": optimum["final_cost"],
                             "physical": optimum["physical"], "qualified": all(r["robust_passed"] == "True" for r in held),
                             "realizations_passing_every_bucket": len(scenarios) - len(failing),
                             "held_out": [{"source": r["source"], "I_bucket": r["I_bucket"], "RMS": float(r["RMS"]),
                                           "maximum_absolute_error": float(r["maximum_absolute_error"]), "bias": float(r["bias"]),
                                           "limit_RMS": float(r["limit_RMS"]), "limit_maximum": float(r["limit_maximum"]),
                                           "native_failure_count": int(r["native_failure_count"]),
                                           "passing_realizations": int(r["passing_realizations"]),
                                           "robust_passed": r["robust_passed"] == "True",
                                           "frozen_robust_passed": r["frozen_robust_passed"] == "True",
                                           "failing_realizations": r["failing_realizations"]}
                                          for r in held]})
        best = verdicts[0]
        nominal = read(scratch / best["label"] / "nominal-residuals.csv")
        write(results / f"nominal-residuals-{form}.csv", [r for r in nominal if r["declaration"] == "candidate"])
        write(results / "nominal-residuals-unfitted.csv", [r for r in nominal if r["declaration"] == "unfitted"])
        selection["forms"][form] = {**best, "labels": labels, "start": best["starts"][0], "optima": verdicts,
                                    "verdict_differs_between_optima": len({v["qualified"] for v in verdicts}) > 1}
        print(form, [(v["label"], v["qualified"], v["realizations_passing_every_bucket"],
                      [(b["source"][:5], b["I_bucket"], b["passing_realizations"], b["failing_realizations"]) for b in v["held_out"]])
                     for v in verdicts], flush=True)
    (results / "selection.json").write_text(json.dumps(selection, indent=2) + "\n")


if __name__ == "__main__":
    main()
