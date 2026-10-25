"""Native Cai MEA-water topology selection with the frozen pressure-grouped scoring."""

from __future__ import annotations

import csv
import hashlib
import importlib.metadata
import importlib.util
import json
import math
import statistics
import subprocess
from copy import deepcopy
from pathlib import Path
from typing import Any
from urllib.parse import unquote, urlparse

import epcsaft
from epcsaft import regression, _core

HERE = Path(__file__).resolve().parents[1] / "results/analysis/neutral-mea-water"
CAMPAIGN = Path(__file__).resolve().parents[1]
ROOT = Path(__file__).resolve().parents[3]
PARAMETERS = (
    CAMPAIGN / "inputs/data/mea-co2-h2o-nine-species-estimation/1/parameters/parameters.json"
)
CAI = (
    CAMPAIGN
    / "inputs/source/library/observations/vapor_liquid_equilibrium/Cai_1996_MEA_H2O_VLE.csv"
)
RESULT = HERE / "results/neutral-mea-water.json"
PLOT_DATA = HERE / "figures/output/topology-selection.csv"
ENGINE_SUBJECT = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, check=True, capture_output=True, text=True).stdout.strip()
spec = importlib.util.spec_from_file_location("cai_fit", ROOT / "analyses/2026-wolbach-sandler-requalification/scripts/mea_water.py")
cai_fit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cai_fit)
DOMAIN = "mea-candidate-293-15-to-393-15-k"
COORDINATE = "pair/monoethanolamine/water/k_ij"
SCHEMES = {"2B": (1, 1), "3B": (1, 2), "4C": (2, 2)}
BOUNDS = (-0.5, 0.5)
FROZEN_ORIGIN = -0.052
SCALES = {"temperature_k": 1.0, "vapor_water_mole_fraction": 0.01}
ALLOWED_ROLES = {"binary_model_selection", "binary_training"}


def sha256(path: Path) -> str:
    return f"sha256:{hashlib.sha256(path.read_bytes()).hexdigest()}"


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def installed_wheel() -> Path:
    direct_url = json.loads(
        importlib.metadata.distribution("epcsaft").read_text("direct_url.json") or "{}"
    ).get("url", "")
    wheel = Path(unquote(urlparse(direct_url).path))
    if (
        not direct_url.startswith("file://")
        or wheel.suffix != ".whl"
        or not wheel.is_file()
    ):
        raise RuntimeError(
            "epcsaft must be installed from one local non-editable wheel"
        )
    return wheel


def cai_rows() -> list[dict[str, str]]:
    with CAI.open(newline="", encoding="utf-8-sig") as stream:
        source_rows = list(csv.DictReader(stream))
    unexpected = {
        row["role"]
        for row in source_rows
        if row["role"] not in ALLOWED_ROLES | {"pure_endpoint_context"}
    }
    if unexpected:
        raise RuntimeError(
            f"Cai admission contains forbidden roles: {sorted(unexpected)}"
        )
    rows = [row for row in source_rows if row["role"] in ALLOWED_ROLES]
    if len(rows) != 25 or sorted({float(row["pressure_pa"]) for row in rows}) != [
        66660.0,
        101330.0,
    ]:
        raise RuntimeError("the frozen 25-row, two-pressure Cai admission changed")
    return rows


def scheme_parameters(source: dict[str, Any], scheme: str) -> epcsaft.Parameters:
    mapping = deepcopy(source)
    next(domain for domain in mapping["domains"] if domain["domain_id"] == DOMAIN)[
        "temperature_max"
    ] = {"magnitude": 460.0, "unit": "kelvin"}
    # The retained loader excluded like-charge dispersion; both selected species are neutral.
    next(f for f in mapping["model_families"] if f["kind"] == "electrolyte")["ionic_dispersion"] = "exclude-same-sign-ion-pairs"
    sites = mapping["topology"]["sites"]
    mea_sites = [site for site in sites if site["component_id"] == "monoethanolamine"]
    if {(site["site_id"], site["site_role"]) for site in mea_sites} != {
        ("a", "donor"),
        ("b", "acceptor"),
    }:
        raise RuntimeError("retained MEA donor/acceptor identities changed")
    multiplicities = dict(zip(("a", "b"), SCHEMES[scheme], strict=True))
    for site in mea_sites:
        site["multiplicity"] = multiplicities[site["site_id"]]

    edges = mapping["topology"]["edges"]
    self_edges = {
        edge["endpoint_a"]["component_id"]: edge
        for edge in edges
        if edge["endpoint_a"]["component_id"] == edge["endpoint_b"]["component_id"]
        and edge["endpoint_a"]["component_id"] in {"monoethanolamine", "water"}
    }
    if set(self_edges) != {"monoethanolamine", "water"}:
        raise RuntimeError("retained MEA/water self-association inputs changed")
    retained = [
        edge
        for edge in edges
        if {edge["endpoint_a"]["component_id"], edge["endpoint_b"]["component_id"]}
        != {"monoethanolamine", "water"}
    ]
    provenance = sorted(
        self_edges["monoethanolamine"]["source"]["provenance"]
        + self_edges["water"]["source"]["provenance"],
        key=lambda item: (item["source_id"], item["locator"], item["domain_id"]),
    )
    inputs = {
        "energy_over_k": [
            self_edges[component]["energy_over_k"]["identity"]
            for component in ("monoethanolamine", "water")
        ],
        "volume": [
            self_edges[component]["volume"]["identity"]
            for component in ("monoethanolamine", "water")
        ],
    }
    candidate = {
        field: self_edges["monoethanolamine"][field]
        for field in ("candidate_domain_id", "qualification", "source_sha256")
    }
    for mea_site, water_site in (("a", "b"), ("b", "a")):
        retained.append(
            {
                **candidate,
                "endpoint_a": {"component_id": "monoethanolamine", "site_id": mea_site},
                "endpoint_b": {"component_id": "water", "site_id": water_site},
                "source": {
                    "kind": "combining-rule",
                    "rule_id": "arithmetic-energy-geometric-volume",
                    "inputs": inputs,
                    "provenance": provenance,
                },
            }
        )
    mapping["topology"]["edges"] = retained
    return epcsaft.Parameters.from_mapping(
        mapping, components=("water", "monoethanolamine")
    )


def evaluate(parameters, rows, value):
    coordinate = regression.coordinate(parameters, [("k_ij", cai_fit.COMPONENTS, 1)],
                                       origin=value, scale=1.0, bounds=BOUNDS)
    items = cai_fit.observations(parameters, rows)
    evaluated = _core.evaluate_observations(parameters._native, [coordinate], [0.0], [], items)
    residuals, derivatives, records = [], [], []
    for i, row in enumerate(rows):
        available = all(status == regression.ObservationStatus.Available for status in evaluated.statuses[2*i:2*i+2])
        record = {"identity": row["row_id"], "attempted": True,
                  "outcome": "evaluated" if available else "unavailable",
                  "statuses": [str(x) for x in evaluated.statuses[2*i:2*i+2]]}
        if available:
            temperature, vapor = evaluated.predictions[2*i:2*i+2]
            pair = (temperature-float(row["temperature_k"]),
                    (vapor-float(row["water_vapor_mole_fraction"]))/0.01)
            residuals.extend(pair)
            derivatives.extend((evaluated.optimizer_gradient[2*i], evaluated.optimizer_gradient[2*i+1]/0.01))
            record.update(pressure_pa=float(row["pressure_pa"]),
                          liquid_water_mole_fraction=float(row["water_liquid_mole_fraction"]),
                          observed_temperature_k=float(row["temperature_k"]), predicted_temperature_k=temperature,
                          observed_vapor_water_mole_fraction=float(row["water_vapor_mole_fraction"]),
                          predicted_vapor_water_mole_fraction=vapor,
                          scaled_temperature_residual=pair[0], scaled_vapor_water_residual=pair[1])
        else:
            chart = regression.evaluate_coordinates(parameters, [coordinate], [0.0], [])
            model = epcsaft.Mixture(epcsaft.Parameters(chart.parameters, parameters.component_ids, parameters.site_ids))
            # Read the owner's cause for failed rows; no retry or substitute residual.
            record["diagnostics"] = [regression.predict(model, item).message for item in items[2*i:2*i+2]]
        records.append(record)
    return residuals, derivatives, records


def fit(parameters, rows, origin):
    print(f"fit multiplicities={parameters._document.site_multiplicity} rows={len(rows)}", flush=True)
    # Diagnosis: function stopping at 1e-12 missed the frozen scaled-gradient criterion.
    controls = regression.FitControls(maximum_iterations=100, maximum_elapsed_time_seconds=600.0,
                                      function_tolerance=1e-15, gradient_tolerance=1e-12)
    result = cai_fit.fit_rows(parameters, rows, origin, BOUNDS, controls)
    value = result.physical[0] if result.physical else origin
    residuals, derivatives, records = evaluate(parameters, rows, value)
    phase_success = len(residuals) == 2*len(rows)
    gradient = math.fsum(r*d for r,d in zip(residuals, derivatives, strict=True)) if phase_success else None
    curvature = math.fsum(d*d for d in derivatives) if phase_success else None
    gradient_pass = phase_success and abs(gradient) <= 1e-8*max(1.0, curvature)
    usable = result.usable and gradient_pass
    standard_error = math.sqrt(result.covariance.optimizer[0] * result.covariance.variance_factor) if result.covariance.status == "available" else None
    confidence = [value-1.96*standard_error, value+1.96*standard_error] if standard_error is not None else None
    bounds_pass = BOUNDS[0] < value < BOUNDS[1] and not result.covariance.active_bounds
    uncertainty_pass = (standard_error is not None and math.isfinite(standard_error) and standard_error > 0
                        and BOUNDS[0] < confidence[0] < confidence[1] < BOUNDS[1])
    print(f"k={value:.12g} cost={result.final_cost:.12g} gradient={gradient} curvature={curvature} gradient_pass={gradient_pass}", flush=True)
    return {"usable": usable, "converged": result.status == "converged" and gradient_pass,
            "phase_success": phase_success, "value": value, "iterations": result.iterations,
            "objective": result.final_cost if result.usable else None,
            "scaled_rmse": math.sqrt(2*result.final_cost/(2*len(rows))) if result.usable else None,
            "failures": len(rows)-len(residuals)//2, "records": records,
            "native": cai_fit.fit_numbers(result),
            "exact_gradient": gradient, "exact_curvature": curvature, "gradient_pass": gradient_pass,
            "bounds": {"interval": list(BOUNDS), "active": bool(result.covariance.active_bounds), "pass": bounds_pass},
            "uncertainty": {"method": "local-linearized-covariance-with-residual-variance",
                            "degrees_of_freedom": 2*len(rows)-1, "standard_error": standard_error,
                            "confidence_level": 0.95, "confidence_interval": confidence, "pass": uncertainty_pass,
                            "limitation": "conditional on the selected topology, frozen residual scales, and local linearization"}}


def compare_topology(source, rows, result):
    old = json.loads(cai_fit.historical(RESULT))
    table, diagnostics = [], []
    old_candidates = {c["scheme"]: c for c in old["candidates"]}
    for candidate in result["candidates"]:
        scheme = candidate["scheme"]
        prior = old_candidates[scheme]
        parameters = scheme_parameters(source, scheme)
        pairs = [("complete", rows, prior["complete_refit"], candidate["complete_refit"])]
        for before, after in zip(prior["folds"], candidate["folds"], strict=True):
            pressure = before["held_pressure_pa"]
            training = [r for r in rows if float(r["pressure_pa"]) != pressure]
            pairs.append((f"training_hold_{pressure}", training, before["fit"], after["fit"]))
            held = [r for r in rows if float(r["pressure_pa"]) == pressure]
            old_held_residuals,_,old_held_native = evaluate(parameters, held, before["fit"]["value"])
            table.append(cai_fit.comparison(scheme,"all","held_scaled_mean_square",f"held_{pressure}_native_at_historical_coordinate",
                                           before["held_scaled_mean_square"],
                                           math.fsum(r*r for r in old_held_residuals)/len(old_held_residuals),
                                           "objective" if before["held_scaled_mean_square"] is not None else "secondary"))
            for b,replayed,a in zip(before["held_records"], old_held_native, after["held_records"], strict=True):
                for quantity,key in (("temperature_k", "predicted_temperature_k"), ("vapor_water", "predicted_vapor_water_mole_fraction")):
                    table.extend((cai_fit.comparison(scheme,b["identity"],quantity,f"held_{pressure}_native_at_historical_coordinate",b.get(key),replayed.get(key)),
                                  cai_fit.comparison(scheme,b["identity"],quantity,f"held_{pressure}",b.get(key),a.get(key))))
            table.append(cai_fit.comparison(scheme,"all","held_scaled_mean_square",f"held_{pressure}",
                                           before["held_scaled_mean_square"],after["held_scaled_mean_square"],
                                           "objective" if before["held_scaled_mean_square"] is not None else "secondary"))
        for case, fit_rows, before, after in pairs:
            old_residuals,old_columns,replayed = evaluate(parameters, fit_rows, before["value"])
            if len(old_residuals) != 2*len(fit_rows):
                raise RuntimeError(f"historical coordinate replay unavailable: {scheme} {case}")
            old_native_sse = math.fsum(r*r for r in old_residuals)
            if before.get("objective") is not None:
                table.append(cai_fit.comparison(scheme,"all","sse",f"{case}_native_at_historical_coordinate",
                                               2*before["objective"],old_native_sse,"objective"))
            table.append(cai_fit.comparison(scheme,"all","exact_gradient",f"{case}_native_at_historical_coordinate",
                                           None,math.fsum(r*d for r,d in zip(old_residuals,old_columns,strict=True)),"secondary"))
            coordinate = regression.coordinate(parameters, [("k_ij", cai_fit.COMPONENTS, 1)], origin=after["value"], scale=1.0, bounds=BOUNDS)
            chart = regression.evaluate_coordinates(parameters, [coordinate], [0.0], [])
            model = epcsaft.Mixture(epcsaft.Parameters(chart.parameters, parameters.component_ids, parameters.site_ids))
            for row in fit_rows:
                cai_fit.bubble(model, row, diagnostics)
                diagnostics[-1].update(case=scheme, stage=case, k_ij=after["value"])
            for old_record, replay, new_record in zip(before["records"], replayed, after["records"], strict=True):
                for quantity,key in (("temperature_k", "predicted_temperature_k"), ("vapor_water", "predicted_vapor_water_mole_fraction")):
                    table.extend((cai_fit.comparison(scheme,old_record["identity"],quantity,f"{case}_native_at_historical_coordinate",old_record.get(key),replay.get(key)),
                                  cai_fit.comparison(scheme,old_record["identity"],quantity,f"{case}_native_optimum",old_record.get(key),new_record.get(key))))
            table.append(cai_fit.comparison(scheme,"all","k_ij",case,before["value"],after["value"],"secondary"))
            if before.get("objective") is not None:
                table.append(cai_fit.comparison(scheme,"all","sse",case,2*before["objective"],2*after["objective"],"objective"))
            old_uncertainty, new_uncertainty = before.get("uncertainty"), after.get("uncertainty")
            if old_uncertainty and new_uncertainty:
                table.append(cai_fit.comparison(scheme,"all","standard_error",case,old_uncertainty["standard_error"],new_uncertainty["standard_error"],tolerance=1e-7))
        table.append(cai_fit.comparison(scheme,"all","grouped_mean_score","selection",prior["grouped_mean_score"],candidate["grouped_mean_score"],
                                       "objective" if prior["grouped_mean_score"] is not None else "secondary"))
    with (RESULT.parent / "batch1-native-rows.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(diagnostics[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(diagnostics)
    cai_fit.write_comparison(RESULT.parent / "batch1-before-after.csv", table)
    result["comparison"] = {"base_commit": cai_fit.REFERENCE_COMMIT, "historical_engine_commit": old["engine_subject_commit"],
                            "old_wheel_sha256": old["wheel_sha256"], "old_selection": old["selection"],
                            "table": "batch1-before-after.csv", "failed_comparisons": sum(r["passes"] is False for r in table),
                            "claim_changes": ["4C now passes all fit gates and has a grouped score; its historical complete and first-fold fits stopped at the 350-460 K bracket. The retained witness in batch1-diagnosis.json lies above that bracket. 2B remains the only one-SE topology. This is candidate extrapolation, not physical validation."],
                            "score_comparison_limits": "Held-pressure scores can change when training coordinates satisfy the same stopping criterion at different points. Native-at-historical-coordinate rows separate state/residual changes from parameter stopping differences; all comparison misses remain in the table.",
                            "historical_incomplete_fits": [f"{c['scheme']} {name}" for c in old["candidates"] for name,f in [("complete",c["complete_refit"])]+[(f"hold_{x['held_pressure_pa']}",x["fit"]) for x in c["folds"]] if not f["usable"]]}


def main() -> None:
    wheel = installed_wheel()
    source = read_json(PARAMETERS)
    packet_temperature_max_k = float(
        next(domain for domain in source["domains"] if domain["domain_id"] == DOMAIN)[
            "temperature_max"
        ]["magnitude"]
    )
    if packet_temperature_max_k != 393.15:
        raise RuntimeError("the frozen packet candidate-domain maximum changed")
    rows = cai_rows()
    pressures = sorted({float(row["pressure_pa"]) for row in rows})
    candidates: list[dict[str, Any]] = []
    for scheme in SCHEMES:
        parameters = scheme_parameters(source, scheme)
        origin = FROZEN_ORIGIN
        folds = []
        for held_pressure in pressures:
            training = [
                row for row in rows if float(row["pressure_pa"]) != held_pressure
            ]
            held = [row for row in rows if float(row["pressure_pa"]) == held_pressure]
            fitted = fit(parameters, training, origin)
            held_residuals, _, held_records = evaluate(
                parameters, held, fitted["value"]
            )
            folds.append(
                {
                    "held_pressure_pa": held_pressure,
                    "fit": fitted,
                    "held_scaled_mean_square": (
                        math.fsum(r * r for r in held_residuals) / len(held_residuals)
                        if len(held_residuals) == 2 * len(held)
                        else None
                    ),
                    "held_records": held_records,
                }
            )
        complete = fit(parameters, rows, origin)
        phase_success = complete["usable"] and all(
            fold["fit"]["usable"]
            and fold["held_scaled_mean_square"] is not None
            and all(record["outcome"] == "evaluated" for record in fold["held_records"])
            for fold in folds
        )
        scores = [fold["held_scaled_mean_square"] for fold in folds]
        gates = {
            "phase_success": phase_success,
            "optimizer_convergence": complete.get("converged", False)
            and all(fold["fit"].get("converged", False) for fold in folds),
            "bounds": complete.get("bounds", {}).get("pass", False),
            "uncertainty": (complete.get("uncertainty") or {}).get("pass", False),
        }
        candidates.append(
            {
                "scheme": scheme,
                "association": {"mea_site_multiplicities": dict(zip(("a", "b"), SCHEMES[scheme], strict=True)),
                                "combining_rule": "arithmetic-energy-geometric-volume", "source_sha256": sha256(PARAMETERS)},
                "folds": folds,
                "grouped_mean_score": statistics.fmean(scores)
                if all(s is not None for s in scores)
                else None,
                "grouped_standard_error": (
                    statistics.stdev(scores) / math.sqrt(len(scores))
                    if all(s is not None for s in scores)
                    else None
                ),
                "complete_refit": complete,
                "gates": gates,
                "eligible": all(gates.values()),
            }
        )

    eligible = [candidate for candidate in candidates if candidate["eligible"]]
    if not eligible:
        raise RuntimeError(
            "no topology passed phase-success, bounds, and uncertainty gates"
        )
    best = min(eligible, key=lambda candidate: candidate["grouped_mean_score"])
    threshold = best["grouped_mean_score"] + best["grouped_standard_error"]
    within_one_se = [
        candidate
        for candidate in eligible
        if candidate["grouped_mean_score"] <= threshold
    ]
    selected = min(
        within_one_se, key=lambda candidate: tuple(SCHEMES).index(candidate["scheme"])
    )
    result = {
        "schema": "epcsaft.validation.neutral-mea-water",
        "schema_version": 1,
        "classification": "MODERN_REFIT",
        "engine_subject_commit": ENGINE_SUBJECT,
        "engine_tree": subprocess.run(
            ["git", "rev-parse", "HEAD:engine"],
            cwd=ROOT,
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip(),
        "wheel_sha256": sha256(wheel),
        "inputs": {
            "parameters": {
                "path": str(PARAMETERS.relative_to(ROOT)),
                "sha256": sha256(PARAMETERS),
            },
            "cai": {"path": str(CAI.relative_to(ROOT)), "sha256": sha256(CAI)},
            "source": "Cai et al. (1996), Tables 2-3, DOI 10.1021/je960118o",
            "roles_used": ["binary_model_selection", "binary_training"],
            "reactive_or_holdout_rows_used": 0,
        },
        "domain": {
            "row_count": len(rows),
            "temperature_k": [
                min(float(row["temperature_k"]) for row in rows),
                max(float(row["temperature_k"]) for row in rows),
            ],
            "pressure_pa": pressures,
            "selected_coordinate_evidence_range": "the admitted Cai observation range only",
            "packet_candidate_temperature_max_k": packet_temperature_max_k,
            "historical_execution_candidate_temperature_max_k": 460.0,
            "native_free_temperature": "unbounded driver solve; predictions above 393.15 K are candidate extrapolation; historical 350-460 K bracket is a comparison limit",
            "historical_temperature_bracket_k": [350.0, 460.0],
            "override_reason": "the admitted Cai rows extend above the provisional packet's 393.15 K candidate range",
            "qualification": "candidate_extrapolation above 393.15 K; no claim outside the admitted Cai observations",
        },
        "fit": {
            "coordinate": COORDINATE,
            "bounds": list(BOUNDS),
            "origin": FROZEN_ORIGIN,
            "origin_role": "frozen pre-selection numerical start; not a fitted result",
            "native_controls": {"function_tolerance": 1e-15, "parameter_tolerance": "native default (1e-10; historical step cutoff was 1e-9)", "gradient_tolerance": 1e-12,
                                "reason": "diagnosed premature function stopping; frozen exact scaled-gradient criterion unchanged"},
            "residual_scales": SCALES,
            "groups": "leave-one-pressure-level-out",
            "selection": "simplest eligible topology within the best grouped mean plus its standard error",
        },
        "candidates": candidates,
        "selection": {
            "status": "selected_without_one_se_ambiguity"
            if len(within_one_se) == 1
            else "selected_by_simplicity_within_one_se_ambiguity",
            "selected_scheme": selected["scheme"],
            "selected_k_mea_water": selected["complete_refit"]["value"],
            "one_se_threshold": threshold,
            "within_one_se": [candidate["scheme"] for candidate in within_one_se],
            "falsification_condition": "reopen if a topology outside the retained one-SE set passes all gates under the frozen Cai rows and scoring, or if the selected fit loses phase success, interior bounds, or finite interior uncertainty",
        },
        "claim_boundary": "Cai MEA-water VLE only; the fitted k_ij is correlation evidence, not reactive-mixture or holdout validation",
    }
    compare_topology(source, rows, result)
    write_json(RESULT, result)
    PLOT_DATA.parent.mkdir(parents=True, exist_ok=True)
    with PLOT_DATA.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=(
                "scheme",
                "grouped_mean_score",
                "grouped_standard_error",
                "eligible",
                "selected",
            ),
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(
            {
                "scheme": candidate["scheme"],
                "grouped_mean_score": candidate["grouped_mean_score"],
                "grouped_standard_error": candidate["grouped_standard_error"],
                "eligible": str(candidate["eligible"]).lower(),
                "selected": str(candidate is selected).lower(),
            }
            for candidate in candidates
        )
    print(json.dumps(result["selection"], indent=2, sort_keys=True))


def declarations():
    from epcsaft import regression as reg

    return [
        (
            p,
            reg.coordinate(
                p, [('k_ij', cai_fit.COMPONENTS, 1)], origin=FROZEN_ORIGIN, scale=1.0, bounds=BOUNDS
            ),
            cai_fit.observations(p, cai_rows()),
        )
        for p in (scheme_parameters(read_json(PARAMETERS), s) for s in SCHEMES)
    ]


if __name__ == "__main__":
    main()
