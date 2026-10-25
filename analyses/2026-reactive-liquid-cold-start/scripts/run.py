"""Cold starts of the reacting MEA liquid (#154, #157): reproduction, independent-route agreement, sweep, timing.

Usage: run.py {states,sweep,timing} LABEL. LABEL names the installed wheel (before: 684b213a, after: this
branch); every JSON records its SHA-256. Application inputs are read-only clones (cases.py).
"""
from __future__ import annotations

import csv
import hashlib
import importlib.metadata
import json
import os
import random
import sys
from pathlib import Path
from time import perf_counter
from urllib.parse import unquote, urlparse

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cases  # noqa: E402
from epcsaft import equilibrium as q  # noqa: E402

RESULTS = Path(__file__).resolve().parents[1] / "results"
MODE, LABEL = sys.argv[1], sys.argv[2]


def wheel() -> str:
    url = json.loads(importlib.metadata.distribution("epcsaft").read_text("direct_url.json"))["url"]
    return hashlib.sha256(Path(unquote(urlparse(url).path)).read_bytes()).hexdigest()


PROVENANCE = {"wheel_sha256": wheel(), "label": LABEL, "omp_num_threads": os.environ.get("OMP_NUM_THREADS")}


def solve(model, problem):
    clock = perf_counter()
    try:
        result = q.solve_equilibrium(model, problem)
    except Exception as error:  # noqa: BLE001 - an Engine exception is a failed attempt
        return None, perf_counter() - clock, str(error)
    return result, perf_counter() - clock, result.message


def row(prefix: str, result, seconds: float, message: str) -> dict:
    ok = bool(result and result.success)
    return {f"{prefix}_ok": ok, f"{prefix}_iterations": result.convergence.iterations if result else None,
            f"{prefix}_residual": result.convergence.residual if result else None, f"{prefix}_s": seconds,
            f"{prefix}_message": "" if ok else message[:120]}


def write(name: str, rows: list[dict], summary: dict) -> None:
    RESULTS.mkdir(exist_ok=True)
    with (RESULTS / f"{name}-{LABEL}.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(dict.fromkeys(k for r in rows for k in r)), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    (RESULTS / f"{name}-{LABEL}.json").write_text(json.dumps({**PROVENANCE, **summary}, indent=1) + "\n")


def continuation(model, start, target):
    """Native continuation from a converged start problem to the target (the absorber's fallback route)."""
    options = q.Continuation(initial_step=0.25, minimum_step=1.0 / 64, max_steps=64, max_retries_per_step=4)
    path = q.continue_equilibrium(model, start, target, options)
    return path.states[-1] if path.success else None


def states() -> None:
    """#157 rows cold and by continuation from a converged neighbour; #154 states cold at the pre-refit and
    refit-C parameters and from the other parameter set's solution as the liquid guess (MEA's anchor route)."""
    rows = []
    model, problem, _ = cases.absorber_liquid()
    anchor = (318.15, 111200.0, (0.2 * 9.7292821, 9.7292821, 76.9751925))  # a converged absorber state
    for temperature, pressure, feed in cases.ABSORBER_STATES:
        cold = solve(model, problem(temperature, pressure, feed))
        item = {"issue": 157, "state": f"{temperature} K {pressure} Pa CO2 {feed[0]}", **row("cold", *cold)}
        other = continuation(model, problem(anchor[0], pressure, anchor[2]), problem(temperature, pressure, feed))
        item["other_route"] = "continuation"
        item["other_ok"] = other is not None
        if cold[0] is not None and cold[0].success and other is not None:
            item["max_relative_composition_difference"] = cases.relative_composition_difference(
                cold[0].mole_fractions, other.mole_fractions)
        rows.append(item)
    timing, models = cases.mea_evaluation()
    for identity in cases.MEA_STATES:
        cold = {name: solve(m, timing.problem(identity)) for name, m in models.items()}
        for name, other_name in (("pre_refit", "refit_c"), ("refit_c", "pre_refit")):
            item = {"issue": 154, "state": f"{identity} {name}", **row("cold", *cold[name])}
            seed = cold[other_name][0]
            item["other_route"] = f"anchored on {other_name} solution"
            other = None
            if seed is not None and seed.success:
                other = solve(models[name], timing.problem(identity, timing.anchor(identity, seed)))[0]
            item["other_ok"] = bool(other and other.success)
            if cold[name][0] is not None and cold[name][0].success and item["other_ok"]:
                item["max_relative_composition_difference"] = cases.relative_composition_difference(
                    cold[name][0].mole_fractions, other.mole_fractions)
            rows.append(item)
    for item in rows:
        print(item["state"], item["cold_ok"], item["cold_iterations"], item.get("max_relative_composition_difference"), flush=True)
    write("states", rows, {"failures": [r["state"] for r in rows if not r["cold_ok"]],
                           "worst_route_difference": max((r.get("max_relative_composition_difference", 0.0) for r in rows), default=0.0),
                           "criterion_relative": 1e-8})


def sweep() -> None:
    """#157's random absorber sweep: 310-362 K, 1.05-1.12 bar, CO2/MEA loading 0.03-0.56, seeds 1 and 2."""
    model, problem, _ = cases.absorber_liquid()
    rows = []
    for seed, count in ((1, 150), (2, 200)):
        rng = random.Random(seed)
        for _ in range(count):
            temperature, pressure, loading = rng.uniform(310, 362), rng.uniform(1.05e5, 1.12e5), rng.uniform(0.03, 0.56)
            feed = (loading * 9.7292821, 9.7292821, 76.9751925)
            result, seconds, message = solve(model, problem(temperature, pressure, feed))
            rows.append({"seed": seed, "temperature_k": temperature, "pressure_pa": pressure, "loading": loading,
                         **row("cold", result, seconds, message),
                         "reduced_hessian_minimum": result.reduced_hessian_minimum if result else None})
    iterations = sorted(r["cold_iterations"] for r in rows if r["cold_iterations"] is not None)
    failures = [r for r in rows if not r["cold_ok"]]
    print("failures", len(failures), "median iterations", iterations[len(iterations) // 2], flush=True)
    write("sweep", rows, {"states": len(rows), "failures": len(failures),
                          "failure_messages": sorted({r["cold_message"] for r in failures}),
                          "iterations": {"median": iterations[len(iterations) // 2],
                                         "p90": iterations[int(0.9 * (len(iterations) - 1))], "max": iterations[-1]}})


def timing() -> None:
    """#152's 123-state MEA packet protocol (its timing table), written here under LABEL."""
    module, _ = cases.mea_evaluation()
    module.RESULTS, module.LABEL = RESULTS, LABEL
    module.timing()


if __name__ == "__main__":
    {"states": states, "sweep": sweep, "timing": timing}[MODE]()
