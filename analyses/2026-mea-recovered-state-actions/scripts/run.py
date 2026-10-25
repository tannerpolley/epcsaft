"""Consumed actions on recovered MEA states and the anchored/continuation sweep (#131).

Usage: run.py MEA_CHECKOUT {actions,sweep,heat}. MEA inputs and its production request
mapping are read from MEA_CHECKOUT; the Engine is the installed non-editable wheel.
"""
from __future__ import annotations

import csv
import hashlib
import importlib.metadata
import json
import subprocess
import sys
from pathlib import Path
from time import perf_counter
from urllib.parse import unquote, urlparse

MEA = Path(sys.argv[1]).resolve()
MODE = sys.argv[2]
sys.path[:0] = [str(MEA / "analyses/mea_parameter_bundle/scripts"), str(MEA / "src")]
import numpy as np
import shared_evaluation as shared
import verify_reference_calorics as vrc
import verify_temperature_reference as vt
from epcsaft import ActiveParameter, GAS_CONSTANT_J_PER_MOL_K as R, ParameterFamily, PropertyObservable, equilibrium

RESULTS = Path(__file__).resolve().parents[1] / "results"
CO2 = shared.COMPONENT_IDS.index("carbon-dioxide")
WATER = shared.COMPONENT_IDS.index("water")
NEUTRAL = [shared.COMPONENT_IDS.index(c) for c in ("carbon-dioxide", "monoethanolamine", "water")]
WIDTH = len(shared.COMPONENT_IDS)
# Frozen protocol (#131 body): steps (h, h/2), consumer criterion, resolution limit. The k_ij
# steps are #142's (dimensionless k_ij).
STEPS = {"T": (0.04, 0.02), "P": (200.0, 100.0), "feed": (2e-4, 1e-4), "kij": (1e-3, 5e-4)}
ACTION_RTOL, RESOLUTION_RTOL, BRANCH_RTOL = 1e-5, 1e-3, 1e-8
# #61 state criteria, frozen before the run. ROW_TOL is the Engine's default caller row
# criterion (EQUILIBRIUM_ROW_RESIDUAL) applied to every unscaled row, trace species included;
# the independent interphase check re-evaluates each phase at its solved (T, rho, x). Its phase
# pressure is checked as the Engine's pressure row, |P_phase - P|/(rho_phase R T) <= ROW_TOL (#331).
ROW_TOL, BALANCE_RTOL, FUGACITY_TOL, IDENTITY_RTOL = 1e-10, 1e-10, 1e-9, 1e-6
TRACE_X = 1e-6
CONTINUATION = equilibrium.Continuation(initial_step=0.25, minimum_step=1 / 64, max_steps=64, max_retries_per_step=4)
# Declared with the frozen protocol on every continuation request: the 1% anchored density window bounds each step.
FINE = equilibrium.Continuation(initial_step=0.25, minimum_step=1 / 512, max_steps=512, max_retries_per_step=8)
Kind = equilibrium.SolvedStateObservableKind


def wheel() -> tuple[str, str]:
    url = json.loads(importlib.metadata.distribution("epcsaft").read_text("direct_url.json"))["url"]
    path = Path(unquote(urlparse(url).path))
    return path.name, hashlib.sha256(path.read_bytes()).hexdigest()


WHEEL_NAME, shared.ENGINE_WHEEL_SHA256 = wheel()
shared.verify_wheel()
MODEL = shared.epcsaft.Mixture(shared.load_parameters())
KIJ = "pair/carbon-dioxide/water/k_ij"
MAPPING = shared.parameter_mapping()
KIJ_VALUE = shared.parameter_values(MAPPING)[KIJ]
MODEL_KIJ = shared.epcsaft.Mixture(shared.load_parameters(),
                                   active_parameters=[ActiveParameter(ParameterFamily.Kij, (CO2 * WIDTH + WATER, WATER * WIDTH + CO2))])
# MEA's neutral ideal-gas records; ions are reaction-completed by the Engine (EqID reference_species_calorics).
MODEL_H = shared.epcsaft.Mixture(shared.load_parameters(), thermochemistry=vrc.record())


def kij_model(delta: float):
    mapping = shared.with_parameter_values(MAPPING, {KIJ: KIJ_VALUE + delta})
    return shared.epcsaft.Mixture(shared.epcsaft.Parameters.from_mapping(mapping))
REACTIONS = shared._selected_reactions()
PACKET = {o["identity"]: o["request"] for o in shared.load_state_packet()["observations"]}
REQUESTS: dict[str, dict] = {}


def request(identity: str, temperature: float | None = None) -> dict:
    base = shared.corrected_request(PACKET[identity], REACTIONS)
    return base if temperature is None else vt.request_at(base, temperature)


def loading(req: dict) -> float:
    return float(req["reaction_system"]["feed_amounts_mol"][CO2])


def accept(problem, result, req: dict) -> tuple[bool, str]:
    """The cold-sweep acceptance: Engine success plus MEA balance, charge and pressure checks."""
    if not result.success:
        return False, shared._stop_code(result)[0]
    snapshot = shared._snapshot_from_result(MODEL, problem, req, result)
    return snapshot["status"] == "evaluated", snapshot["failure_code"]


def criteria(problem, result, req: dict) -> dict:
    """#61 intensive checks on an accepted state (values recorded; pass flags against the frozen criteria)."""
    system = req["reaction_system"]
    nu = np.asarray(system["reaction_matrix"], dtype=float)
    feed = np.asarray(system["feed_amounts_mol"], dtype=float)
    x = np.asarray(result.mole_fractions, dtype=float).reshape(result.phase_count, WIDTH)
    total = np.asarray(result.amounts) @ x
    # Species/moiety balance: n_out - n_feed must be a reaction extent (every conserved moiety).
    extent, *_ = np.linalg.lstsq(nu.T, total - feed, rcond=None)
    balance = float(np.max(np.abs(nu.T @ extent - (total - feed))) / max(1.0, float(np.max(feed))))
    charge = abs(float(np.asarray(system["charges"], dtype=float) @ total)) / float(feed.sum())
    raw = np.asarray(result.residuals, dtype=float)
    rows = result.rows
    intensive = [i for i, r in enumerate(rows) if r.kind.name in ("Chemical", "Pressure")]
    trace = [i for i in intensive if rows[i].kind.name == "Chemical" and x[rows[i].phase, rows[i].index] < TRACE_X]
    liquid = vt.phase_index(problem, "liquid")
    chem = {rows[i].index: raw[i] for i in intensive if rows[i].kind.name == "Chemical" and rows[i].phase == liquid}
    affinity = [abs(sum(v * chem[j] for j, v in enumerate(r) if v)) / sum(map(abs, r)) for r in nu]
    states = [MODEL.state(result.temperature, rho=result.molar_densities[p], x=list(x[p]))
              for p in range(result.phase_count)]
    out = {"balance_rel": balance, "charge_rel": charge,
           "raw_row_max": float(np.max(np.abs(raw[intensive]))),
           "trace_row_max": float(np.max(np.abs(raw[trace]))) if trace else 0.0,
           "trace_species": ";".join(sorted({shared.COMPONENT_IDS[rows[i].index] for i in trace})),
           "affinity_max": max(affinity), "min_liquid_x": float(np.min(x[liquid])),
           "eta_liquid": result.packing_fractions[liquid],
           "dPdrho_min": min(s.pressure_density_derivative for s in states),
           "reduced_pressure_residual_liquid": abs(states[liquid].pressure - result.pressure)
           / (result.molar_densities[liquid] * R * result.temperature),
           "reduced_hessian_minimum": result.reduced_hessian_minimum}
    mu = [np.array([np.nan if v is None else v for v in s.residual_chemical_potential_over_rt]) + np.log(np.maximum(x[p], 1e-300) * s.molar_density)
          for p, s in enumerate(states)]
    if result.phase_count > 1:
        vapor = vt.phase_index(problem, "vapor")
        out["fugacity_max"] = float(max(abs(mu[liquid][i] - mu[vapor][i]) for i in NEUTRAL))
        out["eta_vapor"] = result.packing_fractions[vapor]
        out["reduced_pressure_residual_vapor"] = (abs(states[vapor].pressure - result.pressure)
                                                  / (result.molar_densities[vapor] * R * result.temperature))
        out["vapor_support_ok"] = bool(np.all(x[vapor][NEUTRAL] > 0) and np.sum(np.delete(x[vapor], NEUTRAL)) == 0)
    out["criteria_ok"] = bool(
        balance <= BALANCE_RTOL and charge <= BALANCE_RTOL and out["raw_row_max"] <= ROW_TOL
        and out["affinity_max"] <= ROW_TOL and out["min_liquid_x"] > 0 and out["dPdrho_min"] > 0
        and out["reduced_pressure_residual_liquid"] <= ROW_TOL and out["eta_liquid"] > 0.2
        and (result.phase_count == 1 or (out["fugacity_max"] <= FUGACITY_TOL and out["eta_vapor"] < 0.01
                                         and out["reduced_pressure_residual_vapor"] <= ROW_TOL
                                         and out["vapor_support_ok"])))
    return out


def solve(problem, req: dict) -> dict:
    clock = perf_counter()
    try:
        result = equilibrium.solve_equilibrium(MODEL, problem)
    except Exception as exc:  # noqa: BLE001 - an Engine exception is a failed attempt, never a pass
        return {"problem": problem, "result": None, "ok": False, "code": type(exc).__name__,
                "wall_s": perf_counter() - clock, "work": (0, 0)}
    ok, code = accept(problem, result, req)
    return {"problem": problem, "result": result, "ok": ok, "code": code, "wall_s": perf_counter() - clock,
            "work": (result.iterations, result.evaluations)}


def cold(req: dict) -> dict:
    """MEA cold path: the packet start; a free-pressure bubble seeds itself natively from its speciated liquid (#152)."""
    outcome = solve(shared._problem_from_request(req), req)
    outcome["path"] = "cold-packet-start"
    return outcome


def observables(problem, calorics: bool = False) -> list[tuple[str, object]]:
    liquid = vt.phase_index(problem, "liquid")
    if calorics:  # what MEA's heat and Cp stages consume on a liquid-only TP state
        return [("H_L", equilibrium.SolvedStateObservable(Kind.PhaseProperty, liquid,
                                                          property=PropertyObservable.TotalEnthalpy)),
                ("n_L", equilibrium.SolvedStateObservable(Kind.PhaseAmount, liquid)),
                ("lna_L:CO2", equilibrium.SolvedStateObservable(Kind.PhaseLogActivity, liquid, CO2))]
    items = [(f"x_L:{c}", equilibrium.SolvedStateObservable(Kind.PhaseMoleFraction, liquid, i))
             for i, c in enumerate(shared.COMPONENT_IDS)]
    items += [("lna_L:CO2", equilibrium.SolvedStateObservable(Kind.PhaseLogActivity, liquid, CO2)),
              ("n_L", equilibrium.SolvedStateObservable(Kind.PhaseAmount, liquid)),
              ("V_L", equilibrium.SolvedStateObservable(Kind.PhaseVolume, liquid))]
    items += [(f"mures_L:{c}", equilibrium.SolvedStateObservable(Kind.PhaseResidualChemicalPotentialOverRt, liquid, i))
              for i, c in enumerate(shared.COMPONENT_IDS)]
    if len(problem.phases) > 1:
        vapor = vt.phase_index(problem, "vapor")
        items += [("P", equilibrium.SolvedStateObservable(Kind.Pressure)),
                  ("y:CO2", equilibrium.SolvedStateObservable(Kind.PhaseMoleFraction, vapor, CO2))]
    return items


def direction(name: str) -> equilibrium.SolvedStateActionDirection:
    feed = [0.0] * WIDTH
    if name == "feed":
        feed[CO2] = 1.0
    return equilibrium.SolvedStateActionDirection(float(name == "T"), float(name == "P"), feed,
                                                  [1.0] if name == "kij" else [])


def values(model, problem, result, items) -> list[float]:
    batch = equilibrium.solved_state_first_actions(
        equilibrium.compile_problem(model, problem), result, [o for _, o in items], [direction("none")])[0]
    return [item.value for item in batch.results]


def perturbed(req: dict, name: str, delta: float) -> dict:
    temperature = float(req["temperature"]["value"]) + (delta if name == "T" else 0.0)
    moved = vt.request_at(req, temperature)
    if name == "P":
        moved["pressure"]["value"] = float(moved["pressure"]["value"]) + delta
    if name == "feed":
        moved["reaction_system"]["feed_amounts_mol"][CO2] += delta
    return moved


def action_rows(label: str, path: str, req: dict, central: dict, names: tuple[str, ...],
                model=MODEL, calorics: bool = False) -> list[dict]:
    """Engine A1 against external Richardson re-solves; k_ij re-solves rebuild the model at k_ij +- h."""
    problem, result = central["problem"], central["result"]
    warm = vt.anchor(problem, result)
    items = observables(problem, calorics)
    batches = [equilibrium.solved_state_first_actions(
        equilibrium.compile_problem(MODEL_KIJ if name == "kij" else model, problem), result,
        [o for _, o in items], [direction(name)])[0] for name in names]
    at = lambda name, delta: kij_model(delta) if name == "kij" else model  # noqa: E731
    rows = []
    for name, batch in zip(names, batches, strict=True):
        differences, scheme = [], "central"
        for step in STEPS[name]:
            lower = values(at(name, -step), *vt.solve(at(name, -step), perturbed(req, name, -step), warm), items)
            try:
                upper = values(at(name, step), *vt.solve(at(name, step), perturbed(req, name, step), warm), items)
                differences.append([(u - v) / (2 * step) for u, v in zip(upper, lower, strict=True)])
            except shared.epcsaft.ParameterError as exc:
                # The state sits on a source correlation's upper temperature bound: use the
                # second-order backward difference (same h^2 error order, same Richardson).
                if "outside reaction correlation" not in str(exc):
                    raise
                scheme = "backward"
                centre = values(model, problem, result, items)
                lower2 = values(model, *vt.solve(model, perturbed(req, name, -2 * step), warm), items)
                differences.append([(3 * c - 4 * v + w) / (2 * step)
                                    for c, v, w in zip(centre, lower, lower2, strict=True)])
        for (observable, _), item, coarse, fine in zip(items, batch.results, *differences, strict=True):
            reference = (4 * fine - coarse) / 3
            resolution = abs(fine - coarse) / 3
            error = None if item.action is None else abs(item.action - reference)
            rows.append({
                "state": label, "recovery": path, "direction": name, "scheme": scheme, "observable": observable,
                "unit": item.unit, "status": item.status.name, "value": item.value, "action": item.action,
                "reference": reference, "resolution": resolution, "abs_error": error,
                "criterion": ACTION_RTOL * abs(reference) + resolution,
                "resolved": resolution <= RESOLUTION_RTOL * abs(reference) or resolution == 0.0,
                "passed": error is not None and item.status.name == "Available"
                and error <= ACTION_RTOL * abs(reference) + resolution,
            })
    return rows


def gibbs_helmholtz(label: str, temperature: float, rows: list[dict]) -> dict:
    """Heat of CO2 absorption h_CO2^ig - d(nH)/dn_CO2 against RT^2 d ln a_CO2/dT + RT (MEA's heat identity)."""
    get = lambda direction, observable: next(  # noqa: E731
        r for r in rows if r["state"] == label and r["direction"] == direction and r["observable"] == observable)
    h, n = get("feed", "H_L"), get("feed", "n_L")
    heat = vrc.shomate_enthalpy(temperature, vrc.CO2_SHOMATE, -393510.0) - (n["value"] * h["action"] + h["value"] * n["action"])
    identity = R * temperature**2 * get("T", "lna_L:CO2")["action"] + R * temperature
    th, tn = get("T", "H_L"), get("T", "n_L")
    return {"state": label, "heat_kJ_per_mol": heat / 1e3, "identity_rel": abs(heat - identity) / abs(identity),
            # closed-system d(nh)/dT per mole of liquid; the molar dh/dT alone depends on the enthalpy gauge
            "cp_equilibrium_J_per_mol_K": (tn["value"] * th["action"] + th["value"] * tn["action"]) / tn["value"]}


def liquid_caloric_state(req: dict, central: dict) -> dict:
    """The liquid-only TP state at the recovered state's T and solved P, solved with the thermal record."""
    problem, result = central["problem"], central["result"]
    liquid = vt.request_at(req, float(req["temperature"]["value"]), result.pressure, liquid_only=True)
    problem, result = vt.solve(MODEL_H, liquid, vt.anchor(problem, result))
    return {"problem": problem, "result": result}, liquid


def run_actions() -> dict:
    states = []
    base = request("Bottinger2008_state_050")
    b050 = cold(base)
    states.append(("B050", "cold packet start (#142 state)", base, b050, b050))
    sentinel = cold(vt.request_at(base, 345.0))
    assert sentinel["ok"], sentinel["code"]
    target = vt.request_at(base, 360.0)
    unanchored = cold(target)
    anchored = solve(shared._problem_from_request(target, vt.anchor(sentinel["problem"], sentinel["result"])[0]), target)
    states.append(("R360", "345 K liquid anchor", target, anchored, unanchored))
    for identity in ("vle_obs_0128", "vle_obs_0145", "vle_obs_0227"):
        req = request(identity)
        packet = cold(req)
        states.append((identity, "packet start (#152 native speciated seed)", req, packet, packet))
    rows, summary, identities = [], [], []
    for label, path, req, recovered, reference in states:
        entry = {"state": label, "recovery": path, "recovered_ok": recovered["ok"], "recovered_code": recovered["code"],
                 "without_recovery_ok": reference["ok"], "without_recovery_code": reference["code"]}
        if recovered["ok"]:
            result = recovered["result"]
            entry.update(pressure_pa=result.pressure, work=recovered["work"], wall_s=round(recovered["wall_s"], 3),
                         classified_residual_max=max(map(abs, result.classified_residuals)),
                         **criteria(recovered["problem"], result, req))
            names = ("T", "P", "feed", "kij") if req["pressure"]["role"] == "fixed" else ("T", "feed", "kij")
            rows += action_rows(label, path, req, recovered, names)
            caloric, liquid = liquid_caloric_state(req, recovered)
            rows += action_rows(f"{label}:liquid-TP", "caloric", liquid, caloric, ("T", "feed"), MODEL_H, calorics=True)
            identities.append(gibbs_helmholtz(f"{label}:liquid-TP", float(req["temperature"]["value"]), rows))
        summary.append(entry)
        print(json.dumps(entry, default=str), flush=True)
    write_csv("actions.csv", rows)
    by_direction = {}
    for row in rows:
        entry = by_direction.setdefault(row["direction"] + (":caloric" if row["recovery"] == "caloric" else ""),
                                        {"rows": 0, "passed": 0, "max_rel_error": 0.0})
        entry["rows"] += 1
        entry["passed"] += bool(row["passed"])
        if row["abs_error"] is not None and row["reference"]:
            entry["max_rel_error"] = max(entry["max_rel_error"], row["abs_error"] / abs(row["reference"]))
    return {"states": summary, "actions": len(rows), "passed": sum(r["passed"] for r in rows),
            "unresolved_references": sum(not r["resolved"] for r in rows), "by_direction": by_direction,
            "gibbs_helmholtz": identities,
            "identity_passed": sum(i["identity_rel"] <= IDENTITY_RTOL for i in identities),
            "max_error_over_criterion": max((r["abs_error"] / r["criterion"] for r in rows
                                             if r["abs_error"] is not None and r["criterion"] > 0), default=None)}


def branch_match(outcome: dict, reference: dict) -> bool | None:
    if not (outcome["ok"] and reference["ok"]):
        return None
    a, b = outcome["result"], reference["result"]
    pco2 = lambda r, o: r.pressure * r.mole_fractions[vt.phase_index(o["problem"], "vapor") * WIDTH + CO2]
    return (abs(a.pressure - b.pressure) <= BRANCH_RTOL * b.pressure
            and abs(pco2(a, outcome) - pco2(b, reference)) <= BRANCH_RTOL * pco2(b, reference))


def sweep_row(identity: str, path: str, source: str | None, outcome: dict, reference: dict, extra: dict | None = None) -> dict:
    result = outcome.get("result")
    return {"identity": identity, "path": path, "neighbour": source or "", "ok": outcome["ok"],
            "failure_code": "" if outcome["ok"] else outcome["code"], "same_branch": branch_match(outcome, reference),
            "iterations": outcome["work"][0], "evaluations": outcome["work"][1], "wall_s": round(outcome["wall_s"], 3),
            "pressure_pa": result.pressure if outcome["ok"] else None, **(extra or {}),
            **(criteria(outcome["problem"], result, REQUESTS[identity]) if outcome["ok"] else {})}


def continuation(start: dict, req: dict, options=CONTINUATION) -> dict:
    clock = perf_counter()
    target = shared._problem_from_request(req)
    try:
        path = equilibrium.continue_equilibrium(MODEL, start["problem"], target, options, start["result"])
    except Exception as exc:  # noqa: BLE001 - recorded as a failed continuation
        return {"problem": target, "result": None, "ok": False, "code": f"{type(exc).__name__}: {exc}"[:120],
                "wall_s": perf_counter() - clock, "work": (0, 0), "extra": {}}
    result = path.states[-1] if path.success else None
    ok, code = accept(target, result, req) if result is not None else (False, path.stop.name)
    work = sum(s.iterations for s in path.states), sum(s.evaluations for s in path.states)
    return {"problem": target, "result": result, "ok": ok, "code": code, "wall_s": perf_counter() - clock, "work": work,
            "extra": {"continuation_stop": path.stop.name, "accepted_steps": path.accepted_steps, "attempts": path.attempts,
                      "cause": "" if path.cause_stop is None else path.cause_stop.name}}


def run_sweep() -> dict:
    identities = [i for i, r in PACKET.items() if r["pressure"]["role"] == "solved"]
    REQUESTS.update({i: request(i) for i in identities})
    requests = REQUESTS
    reference = {}
    rows = []
    for identity in identities:
        reference[identity] = cold(requests[identity])
        rows.append(sweep_row(identity, reference[identity]["path"], None, reference[identity], reference[identity]))
        print(identity, "cold", reference[identity]["ok"], reference[identity]["work"], flush=True)
    # Neighbours are chosen from the declared temperature and loading alone; an unsolved neighbour is recorded.
    declared = {i: (round(float(requests[i]["temperature"]["value"]), 2), loading(requests[i])) for i in identities}
    unsolved = {"ok": False, "code": "neighbour_unsolved", "wall_s": 0.0, "work": (0, 0), "extra": {}}
    for identity in identities:
        temperature, load = declared[identity]
        others = [(i, a) for i, a in declared.items() if i != identity]
        same = sorted((a for a in others if a[1][0] == temperature), key=lambda a: abs(a[1][1] - load))
        cross = sorted((a for a in others if a[1][0] != temperature),
                       key=lambda a: (abs(a[1][0] - temperature), abs(a[1][1] - load)))
        for path, pool in (("same-temperature-anchor", same), ("cross-temperature-anchor", cross)):
            if pool:
                source = pool[0][0]
                outcome = unsolved if not reference[source]["ok"] else solve(shared._problem_from_request(
                    requests[identity], vt.anchor(reference[source]["problem"], reference[source]["result"])[0]),
                    requests[identity])
                rows.append(sweep_row(identity, path, source, outcome, reference[identity]))
        if same:
            source = same[0][0]
            for path, options in (("feed-continuation", CONTINUATION), ("feed-continuation-fine", FINE)):
                outcome = continuation(reference[source], requests[identity], options) if reference[source]["ok"] else unsolved
                rows.append(sweep_row(identity, path, source, outcome, reference[identity], outcome["extra"]))
        print(identity, [(r["path"], r["ok"], r["same_branch"]) for r in rows if r["identity"] == identity], flush=True)
    write_csv("sweep.csv", rows)
    counts: dict[str, dict[str, int]] = {}
    for row in rows:
        entry = counts.setdefault(row["path"], {"attempts": 0, "accepted": 0, "same_branch": 0, "iterations": 0,
                                                "criteria_ok": 0, "max": {}})
        entry["attempts"] += 1
        entry["accepted"] += bool(row["ok"])
        entry["criteria_ok"] += bool(row.get("criteria_ok"))
        for key in ("balance_rel", "charge_rel", "raw_row_max", "trace_row_max", "affinity_max", "fugacity_max",
                    "reduced_pressure_residual_liquid", "reduced_pressure_residual_vapor"):
            if row.get(key) is not None:
                entry["max"][key] = max(entry["max"].get(key, 0.0), row[key])
        entry["same_branch"] += bool(row["same_branch"])
        entry["iterations"] += int(row["iterations"])
        if not row["ok"]:
            entry[row["failure_code"]] = entry.get(row["failure_code"], 0) + 1
    return {"requests": len(identities), "paths": counts}


def run_heat() -> dict:
    """MEA's rebuilt calorimetry stage on this wheel; paired with MEA's pinned-wheel comparison."""
    import evaluate_direct_absorption_heat as heat
    model, fingerprint = heat.thermal_model()
    observations = [r for r in heat.read_rows(heat.OBSERVATIONS) if r["regression_eligible"] == "true"]
    _, _, comparison = heat.score(shared._selected_reactions(), model, fingerprint, observations,
                                  shared.EvaluationLimits(60.0, 3600.0, perf_counter() + 3600.0), "engine-61")
    pinned = {r["record_id"]: r for r in heat.read_rows(heat.RESULTS / "current-selected-direct-enthalpy-comparison.csv")}
    shared_rows = [r for r in comparison if r["status"] == "evaluated" and pinned[r["record_id"]]["status"] == "evaluated"]
    write_csv("heat.csv", comparison)
    templates, cold = heat.pressure_templates(heat.load_packet()), {}
    for loading in (0.003, 0.037):  # the endpoints MEA's pinned wheel cannot solve
        req = shared.corrected_request(heat.heat_request(templates, 40, loading), shared._selected_reactions())
        result = equilibrium.solve_equilibrium(MODEL, shared._problem_from_request(req))
        cold[f"40C-{loading}"] = {"success": result.success, "iterations": result.iterations,
                                   "pressure_pa": result.pressure if result.success else None}
    return {**heat.metrics(comparison), "cold_packet_start": cold,
            "by_campaign_partition": {p: heat.metrics([r for r in comparison if r["campaign_partition"] == p])
                                      for p in sorted({r["campaign_partition"] for r in comparison})},
            "pinned_wheel": json.loads((heat.RESULTS / "current-selected-direct-enthalpy-summary.json").read_text())[
                "engine_wheel_sha256"], "paired_with_pinned": len(shared_rows),
            "max_abs_difference_from_pinned_kj_per_mol": max(
                abs(r["predicted_heat_release_kj_per_mol_CO2"]
                    - float(pinned[r["record_id"]]["predicted_heat_release_kj_per_mol_CO2"])) for r in shared_rows)}


def write_csv(name: str, rows: list[dict]) -> None:
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with (RESULTS / name).open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    RESULTS.mkdir(exist_ok=True)
    started = perf_counter()
    summary = {"actions": run_actions, "sweep": run_sweep, "heat": run_heat}[MODE]()
    manifest_path = RESULTS / "summary.json"
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
    manifest[MODE] = {**summary, "wall_s": round(perf_counter() - started, 1)}
    manifest.pop("provenance", None)
    manifest[MODE]["provenance"] = {"wheel": WHEEL_NAME, "wheel_sha256": shared.ENGINE_WHEEL_SHA256,
                              "mea_commit": subprocess.check_output(["git", "-C", str(MEA), "rev-parse", "HEAD"], text=True).strip(),
                              "parameters_sha256": shared.sha256(shared.PARAMETERS),
                              "state_packet_sha256": shared.sha256(shared.STATE_PACKET)}
    manifest_path.write_text(json.dumps(manifest, indent=1, default=str) + "\n", encoding="utf-8")
