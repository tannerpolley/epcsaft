"""MEA reactive refit through the Engine fitter (#152): state timing, exact derivatives, Jacobian cost, refit C.

Usage: run.py MEA_CHECKOUT {timing,derivatives,jacobian,fit} [LABEL]. MEA's production request mapping
(`shared_evaluation._problem_from_request`), state packet and selected reactions are read from MEA_CHECKOUT;
the parameter record is its committed HEAD version (the pre-refit record refit C started from). The Engine is
the installed non-editable wheel; every result records its SHA-256 and the MEA commit.
"""
from __future__ import annotations

import copy
import csv
import hashlib
import importlib.metadata
import json
import os
import statistics
import subprocess
import sys
from pathlib import Path
from time import perf_counter
from urllib.parse import unquote, urlparse

MEA = Path(sys.argv[1]).resolve()
MODE = sys.argv[2]
LABEL = sys.argv[3] if len(sys.argv) > 3 else "after"
sys.path[:0] = [str(MEA / "analyses/mea_parameter_bundle/scripts"), str(MEA / "src")]
import shared_evaluation as shared  # noqa: E402
import epcsaft  # noqa: E402
from epcsaft import equilibrium, regression  # noqa: E402

RESULTS = Path(__file__).resolve().parents[1] / "results"
BUNDLE = MEA / "analyses/mea_parameter_bundle"
RECORD_PATH = "analyses/mea_parameter_bundle/results/selected-current-best-parameters.json"
RECORD_SHA256 = "868a501831b87e95dedf18ce40e9e7ac949f7c6a4aaf137f717cc493ecfcb7be"  # refit C's start record
REFIT_C = json.loads((BUNDLE / "calibration-misfit/refit-C-converged.json").read_text())
SLOPE = "/k_ij/reciprocal_temperature_slope"
IDS = list(shared.COMPONENT_IDS)


def git(*args: str) -> str:
    return subprocess.run(["git", "-C", str(MEA), *args], capture_output=True, text=True, check=True).stdout


def wheel() -> str:
    url = json.loads(importlib.metadata.distribution("epcsaft").read_text("direct_url.json"))["url"]
    return hashlib.sha256(Path(unquote(urlparse(url).path)).read_bytes()).hexdigest()


shared.ENGINE_WHEEL_SHA256 = WHEEL = wheel()
MEA_COMMIT = git("rev-parse", "HEAD").strip()
RECORD_TEXT = git("show", f"HEAD:{RECORD_PATH}")
assert hashlib.sha256(RECORD_TEXT.encode()).hexdigest() == RECORD_SHA256, "MEA HEAD record is not refit C's start"
REACTIONS = shared._selected_reactions()
PACKET = {o["identity"]: o for o in shared.load_state_packet()["observations"]}
PROVENANCE = {"wheel_sha256": WHEEL, "mea_commit": MEA_COMMIT, "record_sha256": RECORD_SHA256,
              "state_packet_sha256": shared.sha256(shared.STATE_PACKET), "omp_num_threads": os.environ.get("OMP_NUM_THREADS")}


def mapping(values: dict[str, float] | None = None) -> dict:
    """The pre-refit record with MEA's runtime defaults, plus named overrides; a missing k_ij 1/T slope node
    is added at value 0 in the form the record uses for CO2-MEA (as calibration-misfit/probe.py does)."""
    result = json.loads(RECORD_TEXT)
    for family in result.get("model_families", ()):
        if family.get("kind") == "electrolyte" and family.get("choice") == "born":
            for key, value in shared.MODEL_RUNTIME_DEFAULTS.items():
                family.setdefault(key, value)
    have = set(shared.parameter_values(result))
    for ident in values or {}:
        if ident.endswith(SLOPE) and ident not in have:
            _, a, b, _, _ = ident.split("/")
            pair = next(p for p in result["pairs"] if {p["component_id_a"], p["component_id_b"]} == {a, b})
            node = copy.deepcopy(next(c for p in result["pairs"] for c in p["coefficients"]
                                      if c["family"] == "k_ij_reciprocal_temperature_slope"))
            node.update({"identity": ident, "value": {"magnitude": 0.0, "unit": "kelvin"}})
            pair["coefficients"].append(node)
    return shared.with_parameter_values(result, values) if values else result


def problem(identity: str, anchor=None):
    return shared._problem_from_request(shared.corrected_request(PACKET[identity]["request"], REACTIONS), anchor)


def solve(model, prob) -> tuple[object | None, float]:
    clock = perf_counter()
    try:
        result = equilibrium.solve_equilibrium(model, prob)
    except Exception:  # noqa: BLE001 - an Engine exception is a failed attempt
        result = None
    return result, perf_counter() - clock


def anchor(identity: str, result) -> shared.Anchor:
    """MEA's anchor (shared.liquid_anchor): liquid pressure, mole fractions and molar volume of a solve."""
    width = len(IDS)
    return shared.Anchor(round(PACKET[identity]["request"]["temperature"]["value"] - 273.15), 0.0,
                         float(result.pressure), tuple(float(x) for x in result.mole_fractions[:width]),
                         1.0 / float(result.molar_densities[0]))


def summary(values: list[float]) -> dict[str, float]:
    ordered = sorted(values)
    return {"n": len(values), "median": statistics.median(values), "p90": ordered[int(0.9 * (len(values) - 1))],
            "max": ordered[-1], "sum": sum(values)}


def timing() -> None:
    """MEA probe protocol per packet state, in process: base cold solve; the refit-C start trial anchored from
    it; a cold trial fallback when the anchor fails. `child_s` repeats the anchored trial through MEA's own
    forked-child path (fork, solve, snapshot with its central action batch, JSON pipe)."""
    base = epcsaft.Mixture(epcsaft.Parameters.from_mapping(mapping()))
    trial = epcsaft.Mixture(epcsaft.Parameters.from_mapping(mapping(dict(zip(REFIT_C["parameters"], REFIT_C["start"])))))
    rows = []
    for identity, item in PACKET.items():
        cold, cold_s = solve(base, problem(identity))
        row = {"identity": identity, "temperature_c": round(item["request"]["temperature"]["value"] - 273.15),
               "kind": "pco2" if item["request"]["pressure"]["role"] == "solved" else "speciation",
               "base_cold_s": cold_s, "base_cold_ok": bool(cold and cold.success),
               "base_cold_iterations": cold.convergence.iterations if cold else None,
               "base_cold_pressure": cold.pressure if cold and cold.success else None}
        seed = anchor(identity, cold) if cold and cold.success else None
        warm, warm_s = solve(trial, problem(identity, seed))
        row |= {"trial_anchor_s": warm_s, "trial_anchor_ok": bool(warm and warm.success),
                "trial_anchor_iterations": warm.convergence.iterations if warm else None}
        fallback_s = 0.0
        if not row["trial_anchor_ok"]:
            fallback, fallback_s = solve(trial, problem(identity))
            row["trial_cold_ok"] = bool(fallback and fallback.success)
        row["trial_protocol_s"] = warm_s + fallback_s
        clock = perf_counter()
        try:
            snapshot = shared._solve_in_child(trial, problem(identity, seed), 120.0,
                                              request=shared.corrected_request(item["request"], REACTIONS))
            row["child_ok"] = snapshot.status == "evaluated"
        except Exception:  # noqa: BLE001
            row["child_ok"] = False
        row["child_s"] = perf_counter() - clock
        rows.append(row)
        print(identity, *(round(row[k], 2) for k in ("base_cold_s", "trial_protocol_s", "child_s")), flush=True)
    RESULTS.mkdir(exist_ok=True)
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with (RESULTS / f"timing-{LABEL}.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    stats = {key: summary([row[key] for row in rows]) for key in ("base_cold_s", "trial_anchor_s", "trial_protocol_s", "child_s")}
    failures = {key: [row["identity"] for row in rows if not row.get(key, True)]
                for key in ("base_cold_ok", "trial_anchor_ok", "trial_cold_ok", "child_ok")}
    (RESULTS / f"timing-{LABEL}.json").write_text(json.dumps({**PROVENANCE, "summary": stats, "failures": failures}, indent=1) + "\n")


# ---- the refit-C problem through epcsaft.regression

SIGMA_LN_P = 0.3  # refit.py: pCO2 ln(pred/obs)/0.3; species (pred - obs)/(0.1 obs + 0.001)
HCO3_POOL = ("bicarbonate-anion", "carbonate-anion")  # NMR "HCO3-" is the fast-exchange carbon pool


def start_parameters():
    return epcsaft.Parameters.from_mapping(mapping(dict(zip(REFIT_C["parameters"], REFIT_C["start"]))))


def coordinates(parameters, count: int = 5) -> list:
    """Refit C's first `count` coordinates at its start, scaled by its sweep steps (scipy x_scale), inside its bounds."""
    result = []
    for ident, x0, lower, upper in list(zip(REFIT_C["parameters"], REFIT_C["start"], *REFIT_C["bounds"]))[:count]:
        _, a, b, *rest = ident.split("/")
        family = "k_ij_reciprocal_temperature_slope" if ident.endswith(SLOPE) else "k_ij"
        result.append(regression.coordinate(parameters, [(family, (a, b), 1)], origin=x0, scale=10.0 if rest[-1] != "k_ij" else 0.01,
                                            bounds=(lower, upper)))
    return result


def active_model(parameters):
    items = regression.CoordinateSet()
    items.coordinates = coordinates(parameters)
    return epcsaft.Mixture(parameters, active_parameters=items.active())


def observations(parameters, identity: str, prob) -> list[tuple[str, object, float]]:
    """(target identity, observation, weight) for one packet state, exactly as refit.py scores it."""
    outputs = {o["identity"]: o for o in PACKET[identity]["request"]["outputs"]}
    names = [phase.name for phase in prob.phases]
    result = []
    for target in PACKET[identity]["targets"]:
        output = outputs[target["prediction_identity"]]
        phase = names.index(output["phase_identity"])
        support = prob.phases[phase].support or IDS
        species = HCO3_POOL if target["identity"].endswith("::HCO3-") else tuple(
            name for name, c in zip(support, output["coefficients"], strict=True) if c)
        assert all(c in (0.0, 1.0) for c in output["coefficients"]), output
        if target["prediction_identity"] == "co2-partial-pressure":
            item = regression.observation(parameters, "partial_pressure", prob, observed=target["observed"], form="log_ratio",
                                          mass_basis=False, phase=output["phase_identity"], species=species)
            result.append((target["identity"], item, 1.0 / SIGMA_LN_P**2))
        else:
            item = regression.observation(parameters, "mole_fraction", prob, observed=target["observed"], form="difference",
                                          scale=0.1 * target["observed"] + 0.001, mass_basis=False,
                                          phase=output["phase_identity"], species=species)
            result.append((target["identity"], item, 1.0))
    return result


def richardson(value, center: float, h: float) -> tuple[float, float]:
    """Centered re-solve differences at h and h/2: Richardson value and resolution |D(h/2) - D(h)|/3."""
    coarse = (value(center + h) - value(center - h)) / (2 * h)
    fine = (value(center + h / 2) - value(center - h / 2)) / h
    return (4 * fine - coarse) / 3, abs(fine - coarse) / 3


def derivatives() -> None:
    """Exact prediction gradients (one batched Jacobian row) against centered re-solves of rebuilt models,
    for every refit-C coordinate; and the batched row against separate single-direction action requests."""
    base = epcsaft.Parameters.from_mapping(mapping())
    start = start_parameters()
    model = active_model(start)
    steps = [1e-4, 1e-4, 1e-4, 1e-4, 0.5]
    rows = []
    for identity in ("vle_obs_0131", "Bottinger2008_state_050", "Bottinger2008_state_059", "vle_obs_0224"):
        seed_result, _ = solve(epcsaft.Mixture(base), problem(identity))
        prob = problem(identity, anchor(identity, seed_result))
        compiled = equilibrium.compile_problem(model, prob)
        solved = equilibrium.solve_equilibrium(model, prob)
        for target, item, _ in observations(start, identity, prob):
            gradient = regression.differentiate(model, item)
            assert gradient.status == regression.ObservationStatus.Available, gradient.message
            # The same row from separate single-direction requests (the pre-#152 per-parameter path).
            kind = equilibrium.SolvedStateObservableKind
            observable = [equilibrium.SolvedStateObservable(kind.PhaseMoleFraction, item.phase, s) for s in item.species]
            partial = item.kind == regression.ObservationKind.PartialPressure
            observable += [equilibrium.SolvedStateObservable(kind.Pressure, 0, 0)] if partial else []
            for a, (ident, h) in enumerate(zip(REFIT_C["parameters"], steps)):
                unit = [0.0] * len(steps)
                unit[a] = 1.0
                direction = equilibrium.SolvedStateActionDirection(0.0, 0.0, [0.0] * len(IDS), unit)
                actions = [r.action for r in equilibrium.solved_state_actions(
                    compiled, solved, equilibrium.SolvedStateActionRequest(observable, [direction])).results]
                fractions = sum(solved.mole_fractions[item.phase * len(IDS) + s] for s in item.species)
                separate = (actions[-1] * fractions + solved.pressure * sum(actions[:-1])) if partial else sum(actions)
                x0 = REFIT_C["start"][a]

                def value(theta: float, a: int = a, ident: str = ident, target_item=item) -> float:
                    shifted = epcsaft.Mixture(epcsaft.Parameters.from_mapping(mapping(
                        dict(zip(REFIT_C["parameters"], REFIT_C["start"])) | {ident: theta})))
                    result = regression.predict(shifted, target_item)
                    assert result.status == regression.ObservationStatus.Available, result.message
                    return result.value

                reference, resolution = richardson(value, x0, h)
                exact = gradient.value[a]
                rows.append({"identity": identity, "target": target, "kind": "partial_pressure" if partial else "mole_fraction",
                             "parameter": ident, "temperature_k": round(solved.temperature, 2), "step": h,
                             "exact": exact, "resolve": reference, "resolution": resolution,
                             "relative_error": abs(exact - reference) / abs(reference) if reference else abs(exact),
                             "separate": separate, "batched_vs_separate": abs(exact - separate) / max(abs(exact), 1e-300)})
                print(identity, target, ident.split("/")[-1], f"{exact:.10e} {reference:.10e} {rows[-1]['relative_error']:.2e}", flush=True)
    write_rows("derivatives", rows, {"rows": len(rows), "criterion_relative": 1e-6,
                                     "worst_relative_error": max(row["relative_error"] for row in rows)})


def steps() -> None:
    """Step study of the weakest-sensitivity row (vle_obs_0224 pCO2 against the carbamate-MEAH+ k_ij): a
    Richardson reference limited by re-solve noise at small h converges onto the exact action as h grows."""
    start = start_parameters()
    identity, a = "vle_obs_0224", 3
    seed, _ = solve(epcsaft.Mixture(epcsaft.Parameters.from_mapping(mapping())), problem(identity))
    _, item, _ = observations(start, identity, problem(identity, anchor(identity, seed)))[0]
    exact = regression.differentiate(active_model(start), item).value[a]
    ident, x0 = REFIT_C["parameters"][a], REFIT_C["start"][a]

    def value(theta: float) -> float:
        shifted = dict(zip(REFIT_C["parameters"], REFIT_C["start"])) | {ident: theta}
        return regression.predict(epcsaft.Mixture(epcsaft.Parameters.from_mapping(mapping(shifted))), item).value

    rows = []
    for h in (1e-4, 5e-4, 1e-3, 2e-3):
        reference, resolution = richardson(value, x0, h)
        rows.append({"identity": identity, "parameter": ident, "step": h, "exact": exact, "resolve": reference,
                     "relative_error": abs(exact - reference) / abs(reference), "resolution": resolution})
    write_rows("step-study", rows, {"rows": len(rows)})


def write_rows(name: str, rows: list[dict], extra: dict | None = None) -> None:
    RESULTS.mkdir(exist_ok=True)
    with (RESULTS / f"{name}.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(dict.fromkeys(k for row in rows for k in row)), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    if extra is not None:
        (RESULTS / f"{name}.json").write_text(json.dumps({**PROVENANCE, **extra}, indent=1) + "\n")


def calibration(base_model) -> list[tuple[str, object]]:
    """Refit C's calibration states (every packet state except 80 degC) that the pre-refit record solves,
    each declared from MEA's anchor of that solve, as refit.py's warm-started probe evaluates them."""
    result = []
    for identity, item in PACKET.items():
        if round(item["request"]["temperature"]["value"] - 273.15) == 80:
            continue
        cold, _ = solve(base_model, problem(identity))
        if cold is not None and cold.success:
            result.append((identity, problem(identity, anchor(identity, cold))))
    return result


def jacobian() -> None:
    """Jacobian cost per calibration state: the pre-#152 regression path (compile, one action request per
    parameter) against one first-action batch on the solve's compiled problem; the four k_ij, and on a wheel
    with the slope family also the fifth coordinate."""
    start = start_parameters()
    families = 5 if hasattr(epcsaft.ParameterFamily, "KijReciprocalTemperatureSlope") else 4
    items = regression.CoordinateSet()
    items.coordinates = coordinates(start, families)
    model = epcsaft.Mixture(start, active_parameters=items.active())
    kind = equilibrium.SolvedStateObservableKind
    rows = []
    for identity, prob in calibration(epcsaft.Mixture(epcsaft.Parameters.from_mapping(mapping()))):
        solved = equilibrium.solve_equilibrium(model, prob)
        if not solved.success:
            rows.append({"identity": identity, "solved": False})
            continue
        observable = [equilibrium.SolvedStateObservable(kind.PhaseMoleFraction, 0, i) for i in range(len(IDS))]
        observable += [equilibrium.SolvedStateObservable(kind.Pressure, 0, 0)]
        for width in sorted({4, families}):
            directions = [equilibrium.SolvedStateActionDirection(0.0, 0.0, [0.0] * len(IDS), [float(k == a) for k in range(families)])
                          for a in range(width)]
            clock = perf_counter()
            compiled = equilibrium.compile_problem(model, prob)
            separate = [equilibrium.solved_state_actions(compiled, solved, equilibrium.SolvedStateActionRequest(observable, [d]))
                        for d in directions]
            separate_s = perf_counter() - clock
            clock = perf_counter()
            batched = equilibrium.solved_state_first_actions(compiled, solved, observable, directions)
            batched_s = perf_counter() - clock
            pairs = [(x.action, y.action) for b, c in zip(batched, separate) for x, y in zip(b.results, c.results)]
            ok = all(x is not None and y is not None for x, y in pairs)
            rows.append({"identity": identity, "solved": True, "directions": width, "separate_s": separate_s, "batched_s": batched_s,
                         "available": ok, "max_relative_difference": max(abs(x - y) / max(abs(x), abs(y), 1e-300) for x, y in pairs) if ok else None})
            print(identity, width, round(separate_s, 3), round(batched_s, 3), rows[-1]["max_relative_difference"], flush=True)
    done = [r for r in rows if r.get("available")]
    stats = {str(w): {k: summary([r[k] for r in done if r["directions"] == w]) for k in ("separate_s", "batched_s")}
             for w in sorted({r["directions"] for r in done})}
    write_rows(f"jacobian-{LABEL}", rows, {"summary": stats, "unsolved_or_unavailable": [r["identity"] for r in rows if not r.get("available")]})


def fit() -> None:
    """Refit C through epcsaft.regression.fit: the five coordinates, the same scaled residuals and weights,
    capped at refit C's six iterations."""
    clock = perf_counter()
    states = calibration(epcsaft.Mixture(epcsaft.Parameters.from_mapping(mapping())))
    setup_s = perf_counter() - clock
    start = start_parameters()
    rows = [row for identity, prob in states for row in observations(start, identity, prob)]
    controls = regression.FitControls(maximum_iterations=6, maximum_elapsed_time_seconds=36000.0)
    clock = perf_counter()
    result = regression.fit(start, coordinates(start), [item for _, item, _ in rows], weights=[w for _, _, w in rows], controls=controls)
    fit_s = perf_counter() - clock
    print(result.status, result.message, result.initial_cost, result.final_cost, result.iterations, fit_s, flush=True)
    write_rows("refit-c-residuals", [{"target": target, "prediction": y, "weighted_residual": r}
                                     for (target, _, _), y, r in zip(rows, result.predictions, result.weighted_residuals)],
               {"status": result.status, "termination_message": result.message, "residuals": len(rows),
                "states": len(states), "initial_cost": result.initial_cost, "final_cost": result.final_cost,
                "scipy_cost": REFIT_C["cost"], "scipy_residuals": REFIT_C["n_residuals"], "iterations": result.iterations,
                "iteration_seconds": list(result.iteration_seconds), "iteration_costs": list(result.iteration_costs),
                "residual_evaluations": result.residual_evaluations, "jacobian_evaluations": result.jacobian_evaluations,
                "trial_failures": [{"trial": f.trial, "observation": f.observation, "status": f.status.name, "message": f.message,
                                    "jacobian": f.jacobian} for f in result.failures],
                "parameters": REFIT_C["parameters"], "start": REFIT_C["start"], "fitted": list(result.physical),
                "scipy_fitted": REFIT_C["fitted"], "active_bounds": list(result.covariance.active_bounds),
                "covariance_status": result.covariance.status, "setup_s": setup_s, "fit_s": fit_s})


if __name__ == "__main__":
    {"timing": timing, "derivatives": derivatives, "jacobian": jacobian, "fit": fit, "steps": steps}[MODE]()
