"""Polar (PCP-SAFT, own parameters) vs nonpolar (PC-SAFT, own parameters) on the same data (#28).

Run from the repository root with the installed wheel:
    .venv/bin/python analyses/2026-polar-dipole/scripts/compare.py
Pure fluids: inputs/reference-saturation.csv (reference-EOS saturation rows, alternating roles) and
inputs/dmso-saturation.csv (measured DMSO rows; a row carries p_sat, rho_liq or both).
Polar parameters: data/packets/gross-vrabec-2006/1 (GV2006 Table 2). Nonpolar parameters: the published
Gross-Sadowski 2001 set where one exists (packet gross-sadowski-2001/3), otherwise a documented refit of
(m, sigma, eps/k) on the calibration rows only. A secondary refit of both models on the same calibration
rows removes the in-sample advantage. Writes results/pure-rows.csv, results/pure-summary.csv and
results/summary.json.
"""
from __future__ import annotations

import copy
import csv
import json
import math
import statistics
from pathlib import Path

from epcsaft import Parameters, Pure, regression
from epcsaft.equilibrium import Amounts, Free, Phase, Pinned, Problem

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parents[1]
GV = ROOT / "data/packets/gross-vrabec-2006/1/parameters/parameters.json"
GS = ROOT / "data/packets/gross-sadowski-2001/3/parameters/parameters.json"
FAMILY = {"segment_count": "segment_count", "segment_diameter": "segment_diameter", "dispersion_energy_over_k": "dispersion_energy_over_k"}
SCALES = {"segment_count": (0.5, (1.0, 6.0)), "segment_diameter": (0.2, (2.5, 4.5)), "dispersion_energy_over_k": (20.0, (100.0, 450.0))}
# GV2006 Table 2 AAD-% (PCP-SAFT, PC-SAFT) for psat and rho_liq.
GV_AAD = {"acetone": ((0.55, 1.61), (1.31, 2.94)), "dimethyl-ether": ((0.24, 0.25), (0.53, 0.80)), "diethyl-ether": ((0.66, 0.68), (1.01, 1.04)),
          "dimethyl-sulfoxide": ((0.46, 1.45), (0.33, 0.58))}
GS_ID = {"dimethyl-ether": "dimethyl-ether", "diethyl-ether": "diethyl-ether"}
# Secondary, added after inspection: the PCP-SAFT fit ranges (acetone: Kleiner and Sadowski 2007 Table 3,
# p. 15547, 200-508 K; ethers and DMSO: GV2006 Table 2). The predeclared all-row disposition stays primary.
SOURCE_RANGE = {"acetone": (200.0, 508.0), "dimethyl-ether": (200.0, 400.0), "diethyl-ether": (220.0, 466.0), "dimethyl-sulfoxide": (291.0, 519.0)}
U_KEY = {"p_rel_error": "u_p_sat_rel", "rho_rel_error": "u_rho_liq_rel"}


def document(component: str, polar: bool, values: dict[str, float] | None = None) -> dict:
    from_gv = polar or component not in GS_ID  # a nonpolar refit starts from the GV values without dipoles
    source = json.loads((GV if from_gv else GS).read_text())
    record = copy.deepcopy(next(c for c in source["components"] if c["component_id"] == (component if from_gv else GS_ID[component])))
    record["component_id"] = component
    if not polar:
        record["coefficients"] = [c for c in record["coefficients"] if c["family"] not in ("dipole_moment", "dipole_count")]
    for coefficient in record["coefficients"]:
        if values and coefficient["family"] in values:
            coefficient["value"]["magnitude"] = values[coefficient["family"]]
    families = [f for f in json.loads(GV.read_text())["model_families"] if polar or f["kind"] != "polar"]
    return {"schema": "epcsaft.parameters", "schema_version": 1, "components": [record], "pairs": [], "model_families": families}


def values(doc: dict) -> dict[str, float]:
    return {c["family"]: c["value"]["magnitude"] for c in doc["components"][0]["coefficients"] if c["family"] in FAMILY}


def saturation(component: str, T: float, guess: float) -> Problem:
    return Problem(phases=[Phase("L", kind="liquid", composition_guess=[1.0]), Phase("V", kind="vapor", amount=Pinned(0.0), composition_guess=[1.0])],
                   T=T, P=Free(guess), feed=Amounts({component: 1.0}), row_tolerance=1.0e-10)


def observations(parameters: Parameters, component: str, rows: list[dict]) -> list[regression.Observation]:
    result = []
    for row in rows:
        T, guess = float(row["T_K"]), row["p_guess"]
        if row["p_sat_Pa"]:
            p = float(row["p_sat_Pa"])
            result.append(regression.observation(parameters, "pressure", saturation(component, T, guess), observed=p, form="log_ratio", mass_basis=False))
        if row["rho_liq_mol_m3"]:
            rho = float(row["rho_liq_mol_m3"])
            result.append(regression.observation(parameters, "phase_density", saturation(component, T, guess), observed=rho, phase="L", form="log_ratio", mass_basis=False))
    return result


def refit(component: str, polar: bool, start: dict, rows: list[dict], *, declare_only=False) -> tuple[dict, dict]:
    # Multistart over eps/k only: a nonpolar start at the polar eps/k can have T_c below the calibration range.
    parameters = Parameters.from_mapping(start, components=[component])
    origins = values(start)
    coordinates = [regression.coordinate(parameters, [(FAMILY[f], component, 1)], origin=origins[f], scale=SCALES[f][0], bounds=SCALES[f][1])
                   for f in FAMILY]
    if declare_only:
        return parameters, coordinates, observations(parameters, component, rows)
    factors = (1.0, 1.05, 1.1)
    chart = regression.CoordinateSet(); chart.coordinates = coordinates
    starts = [chart.start(Parameters.from_mapping(document(component,polar,
        {f: factor*origins[f] if f=="dispersion_energy_over_k" else origins[f] for f in FAMILY}),components=[component])._native,
        observations(parameters,component,rows)) for factor in factors]
    controls = regression.FitControls(maximum_iterations=100, maximum_elapsed_time_seconds=600.0)
    results = regression.fit(parameters, coordinates, observations(parameters, component, rows), starts=starts, controls=controls,
                             shared_iterations=len(starts) * controls.maximum_iterations,
                             shared_seconds=len(starts) * controls.maximum_elapsed_time_seconds)
    if len(results) != len(starts):
        raise RuntimeError("native multistart exhausted its budget before all declared starts")
    attempts = list(zip(factors, results))
    usable = [(f, r) for f, r in attempts if r.usable]
    factor, result = min(usable, key=lambda item: item[1].final_cost) if usable else attempts[0]
    fitted = dict(zip(FAMILY, result.physical))
    cov = result.covariance
    n = len(result.physical)
    correlation = None
    if cov.status == "available":
        diag = [cov.physical[k * n + k] for k in range(n)]
        correlation = [[cov.physical[i * n + j] / math.sqrt(diag[i] * diag[j]) for j in range(n)] for i in range(n)]
    sv = cov.singular_values
    return document(component, polar, fitted), {
        "coordinates": {f: {"origin": origins[f], "scale": SCALES[f][0], "bounds": SCALES[f][1]} for f in FAMILY},
        "optimizer_starts": starts,
        "starts": [{"eps_factor": f, "status": r.status, "message": r.message, "usable": r.usable,
                    "final_cost": r.final_cost if r.usable else None, "sse": 2.0 * r.final_cost if r.usable else None,
                    "fitted": r.physical if r.usable else None, "trial_failures": r.trial_failures,
                    "failures": [{"observation": e.observation, "status": e.status.name, "message": e.message} for e in r.failures],
                    "covariance_status": r.covariance.status, "rank": r.covariance.rank,
                    "singular_values": r.covariance.singular_values, "active_bounds": r.covariance.active_bounds,
                    "physical_covariance": r.covariance.physical, "variance_factor": r.covariance.variance_factor,
                    "projected_gradient_norm": r.covariance.projected_gradient_norm} for f, r in attempts],
        "selected_eps_factor": factor, "status": result.status, "usable": result.usable, "fitted": fitted, "final_cost": result.final_cost, "sse": 2.0 * result.final_cost,
        "trial_failures": result.trial_failures, "covariance_status": cov.status, "rank": cov.rank, "singular_values": sv,
        "condition_number": (max(sv) / min(sv)) if sv and min(sv) > 0 else None, "correlation": correlation,
        "active_bounds": cov.active_bounds, "physical_covariance": cov.physical,
        "variance_factor": cov.variance_factor, "covariance_assumption": cov.assumption}


def predict(doc: dict, component: str, row: dict) -> tuple[float | None, float | None, str]:
    parameters = Parameters.from_mapping(doc, components=[component])
    model = Pure(parameters)
    T, p = float(row["T_K"]), row["p_guess"]
    pressure = regression.predict(model, regression.observation(parameters, "pressure", saturation(component, T, p), observed=p, form="log_ratio", mass_basis=False))
    if pressure.status != regression.ObservationStatus.Available or pressure.value is None:
        return None, None, f"{pressure.status.name}: {pressure.message}"
    density = regression.predict(model, regression.observation(parameters, "phase_density", saturation(component, T, p), observed=1.0, phase="L", form="log_ratio", mass_basis=False))
    if density.status != regression.ObservationStatus.Available:
        return pressure.value, None, f"{density.status.name}: {density.message}"
    return pressure.value, density.value, "Available"


def subset(rows: list[dict], key: str, role: str | None = None, span: tuple[float, float] | None = None, stated: bool = False) -> list[dict]:
    return [r for r in rows if r[key] is not None and (role is None or r["role"] == role)
            and (span is None or span[0] <= r["T_K"] <= span[1]) and (not stated or r[U_KEY[key]] is not None)]


def aad(rows: list[dict], key: str, role: str | None = None, **kw) -> float | None:
    errors = [abs(r[key]) for r in subset(rows, key, role, **kw)]
    return 100.0 * sum(errors) / len(errors) if errors else None


def bound(a_rows: list[dict], b_rows: list[dict], key: str, **kw) -> float | None:
    """Largest change (%) of the AAD gap a-b when every reference value moves within its stated uncertainty.

    Row-wise worst case of (|e_a + d| - |e_b + d|) - (|e_a| - |e_b|) over |d| <= u; the expression is piecewise
    linear in d, so the maximum sits at d = +-u or at a kink d = -e. First order in d (relative errors).
    """
    chosen, other = subset(a_rows, key, **kw), {r["T_K"]: r for r in subset(b_rows, key, **kw)}
    if not chosen or any(r[U_KEY[key]] is None or r["T_K"] not in other for r in chosen):
        return None
    total = 0.0
    for r in chosen:
        ea, eb, u = r[key], other[r["T_K"]][key], r[U_KEY[key]]
        total += max(abs(abs(ea + d) - abs(eb + d) - abs(ea) + abs(eb)) for d in (-u, u, -ea, -eb) if abs(d) <= u)
    return 100.0 * total / len(chosen)


def input_rows(component):
    reference = [row for name in ("reference-saturation.csv", "dmso-saturation.csv") for row in csv.DictReader((HERE / "inputs" / name).open())]
    data = [r for r in reference if r["component_id"] == component]
    # Initial pressure for rows without a p_sat observation: ln p linear in 1/T through the fluid's p_sat rows.
    slope, intercept = statistics.linear_regression([1.0 / float(r["T_K"]) for r in data if r["p_sat_Pa"]],
                                                    [math.log(float(r["p_sat_Pa"])) for r in data if r["p_sat_Pa"]])
    for r in data:
        r["p_guess"] = float(r["p_sat_Pa"]) if r["p_sat_Pa"] else math.exp(intercept + slope / float(r["T_K"]))
    return data


def main() -> None:
    rows, summary, fits = [], [], {}
    for component in GV_AAD:
        data = input_rows(component)
        calibration = [r for r in data if r["role"] == "calibration"]
        models = {"polar/published": document(component, True)}
        if component in GS_ID:
            models["nonpolar/published"] = document(component, False)
        else:
            models["nonpolar/refit"], fits[f"{component}/nonpolar/refit"] = refit(component, False, document(component, False), calibration)
        models["polar/refit"], fits[f"{component}/polar/refit"] = refit(component, True, document(component, True), calibration)
        if component in GS_ID:
            models["nonpolar/refit"], fits[f"{component}/nonpolar/refit"] = refit(component, False, document(component, False), calibration)
        for label, doc in models.items():
            model_rows = []
            for row in data:
                p, rho, status = predict(doc, component, row)
                p_ref = float(row["p_sat_Pa"]) if row["p_sat_Pa"] else None
                rho_ref = float(row["rho_liq_mol_m3"]) if row["rho_liq_mol_m3"] else None
                model_rows.append({"component_id": component, "model": label, "role": row["role"], "T_K": float(row["T_K"]),
                                   "p_ref_Pa": p_ref, "p_model_Pa": p, "rho_ref_mol_m3": rho_ref, "rho_model_mol_m3": rho, "status": status,
                                   "p_rel_error": None if p is None or p_ref is None else p / p_ref - 1.0,
                                   "rho_rel_error": None if rho is None or rho_ref is None else rho / rho_ref - 1.0,
                                   **{u: float(row[u]) if row[u] else None for u in U_KEY.values()}})
            rows += model_rows
            summary.append({"component_id": component, "model": label, **values(doc),
                            "AAD_psat_all": aad(model_rows, "p_rel_error"), "AAD_rho_all": aad(model_rows, "rho_rel_error"),
                            "AAD_psat_validation": aad(model_rows, "p_rel_error", "validation"),
                            "AAD_rho_validation": aad(model_rows, "rho_rel_error", "validation"),
                            "AAD_psat_source_range": aad(model_rows, "p_rel_error", span=SOURCE_RANGE[component]),
                            "AAD_rho_source_range": aad(model_rows, "rho_rel_error", span=SOURCE_RANGE[component]),
                            "AAD_psat_stated_u": aad(model_rows, "p_rel_error", stated=True),
                            "AAD_rho_stated_u": aad(model_rows, "rho_rel_error", stated=True),
                            "failures": sum(r["status"] != "Available" for r in model_rows)})
    results = HERE / "results"
    for name, table in (("pure-rows.csv", rows), ("pure-summary.csv", summary)):
        with (results / name).open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(table[0]), lineterminator="\n")
            writer.writeheader()
            writer.writerows(table)
    # Predeclared disposition (#28 design item 5 as revised): the gap test requires the published PCP-minus-PC sign
    # with at least half its magnitude, using the nonpolar model's own (published or validated refit) set. GV2006's
    # KDB/DIPPR rows are not retained, so absolute AADs are against the named reference, not GV's data. A gap is
    # resolved by the reference when it exceeds the worst-case change a within-uncertainty reference error can cause.
    reproduction = []
    for component, gv in GV_AAD.items():
        polar = next(i for i in summary if i["component_id"] == component and i["model"] == "polar/published")
        nonpolar = next(i for i in summary if i["component_id"] == component and i["model"] == ("nonpolar/published" if component in GS_ID else "nonpolar/refit"))
        polar_rows = [r for r in rows if r["component_id"] == component and r["model"] == "polar/published"]
        nonpolar_rows = [r for r in rows if r["component_id"] == component and r["model"] == nonpolar["model"]]
        for k, (key, name, kw) in enumerate((("AAD_psat_all", "psat", {}), ("AAD_rho_all", "rho_liq", {}),
                                             ("AAD_psat_source_range", "psat, source range (secondary)", {"span": SOURCE_RANGE[component]}),
                                             ("AAD_rho_source_range", "rho_liq, source range (secondary)", {"span": SOURCE_RANGE[component]}),
                                             ("AAD_psat_stated_u", "psat, stated-uncertainty rows (secondary)", {"stated": True}),
                                             ("AAD_rho_stated_u", "rho_liq, stated-uncertainty rows (secondary)", {"stated": True}))):
            k %= 2
            if polar[key] is None or nonpolar[key] is None:
                continue
            ours, theirs = polar[key] - nonpolar[key], gv[k][0] - gv[k][1]
            limit = bound(polar_rows, nonpolar_rows, ("p_rel_error", "rho_rel_error")[k], **kw)
            reproduction.append({"component_id": component, "observable": name, "gv_pcp": gv[k][0], "gv_pc": gv[k][1],
                                 "ours_pcp": polar[key], "ours_pc_own": nonpolar[key], "gv_gap": theirs, "our_gap": ours,
                                 "reference_bound": limit, "gap_resolved": None if limit is None else abs(ours) > limit,
                                 "gap_criterion": "pass" if ours * theirs > 0 and abs(ours) >= 0.5 * abs(theirs) else "fail"})
    preparation = json.loads((ROOT / "build/environment-wheel/preparation.json").read_text())
    out = json.loads((results / "summary.json").read_text())
    if "binary_kij_fits" in out:
        out.setdefault("binary_wheel_sha256", out["wheel_sha256"])
    out.update({"wheel_sha256": preparation["wheel_sha256"], "gv2006_table2_aad": GV_AAD, "fits": fits,
                "pure_summary": summary, "reproduction": reproduction})
    (results / "summary.json").write_text(json.dumps(out, indent=1) + "\n")
    for item in summary:
        print(item["component_id"], item["model"], {k: round(v, 3) if isinstance(v, float) else v for k, v in item.items() if k.startswith("AAD") or k == "failures"})
    for item in reproduction:
        print(item)


def declarations():
    return [
        refit(c, polar, document(c, polar), input_rows(c), declare_only=True)
        for c in GV_AAD
        for polar in (False, True)
    ]


if __name__ == "__main__":
    main()
