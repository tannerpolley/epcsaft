"""Native Ceres fitting of the 39 retained Kiepe CO2-water pressures."""

from __future__ import annotations

import csv
import hashlib
import json
import math
from importlib.metadata import distribution, version
from itertools import pairwise
from pathlib import Path
from typing import Any
from urllib.parse import unquote, urlparse

import epcsaft
import numpy as np
from epcsaft import equilibrium, regression

HERE = Path(__file__).resolve().parents[1] / "results/analysis/physical-co2-water"
ROOT = HERE.parents[4]
CAMPAIGN = HERE.parents[2]
PARAMETERS = (
    CAMPAIGN
    / "inputs/data/mea-co2-h2o-nine-species-estimation/1/parameters/parameters.json"
)
OBSERVATIONS = (
    CAMPAIGN
    / "inputs/source/library/observations/physical_co2_water/kiepe_2002_table_1.csv"
)
CANDIDATE_DOMAIN = "mea-candidate-293-15-to-393-15-k"
K_ID = "pair/carbon-dioxide/water/k_ij"
SLOPE_ID = "pair/carbon-dioxide/water/k_ij/reciprocal_temperature_slope"
T_REF_K = 353.11
BOUNDS = np.array(((-1.0, -1000.0), (1.0, 1000.0)))
SCALES = np.array((1.0, 1000.0))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path: Path, value: Any) -> None:
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def normalize_candidate_provenance(
    value: Any,
    path: str = "$",
    mutations: list[dict[str, str]] | None = None,
) -> list[dict[str, str]]:
    """Apply the campaign's established candidate-domain replay qualification."""
    mutations = [] if mutations is None else mutations
    if isinstance(value, dict):
        if (
            value.get("candidate_domain_id") == CANDIDATE_DOMAIN
            and "provenance" in value
        ):
            original = value["provenance"]["domain_id"]
            mutations.append(
                {
                    "json_path": f"{path}.provenance.domain_id",
                    "from": original,
                    "to": CANDIDATE_DOMAIN,
                }
            )
            value["provenance"]["domain_id"] = CANDIDATE_DOMAIN
        for key, child in value.items():
            normalize_candidate_provenance(child, f"{path}.{key}", mutations)
    elif isinstance(value, list):
        for index, child in enumerate(value):
            normalize_candidate_provenance(child, f"{path}[{index}]", mutations)
    return mutations


def load_model() -> tuple[epcsaft.Parameters, list[dict[str, str]]]:
    mapping = json.loads(PARAMETERS.read_text(encoding="utf-8"))
    component = next(
        item
        for item in mapping["components"]
        if item["component_id"] == "carbon-dioxide"
    )
    segment_count = next(
        item for item in component["coefficients"] if item["family"] == "segment_count"
    )
    if segment_count["value"] != {"magnitude": 2.0729, "unit": "dimensionless"}:
        raise RuntimeError("the corrected CO2 segment count is not exactly 2.0729")

    domain = next(
        item for item in mapping["domains"] if item["domain_id"] == CANDIDATE_DOMAIN
    )
    domain["temperature_max"] = {"magnitude": 393.17, "unit": "kelvin"}
    provenance_mutations = normalize_candidate_provenance(mapping)
    # Current admission requires this declaration; the selected neutral binary
    # has no ionic pairs and the native builder removes its ionic model choices.
    electrolyte = next(
        item for item in mapping["model_families"] if item["kind"] == "electrolyte"
    )
    electrolyte["ionic_dispersion"] = "exclude-same-sign-ion-pairs"
    pair = next(
        item
        for item in mapping["pairs"]
        if {item["component_id_a"], item["component_id_b"]}
        == {"carbon-dioxide", "water"}
    )
    source = next(item for item in pair["coefficients"] if item["identity"] == K_ID)
    pair["coefficients"].append(
        {
            "candidate_domain_id": source["candidate_domain_id"],
            "family": "k_ij_reciprocal_temperature_slope",
            "identity": SLOPE_ID,
            "provenance": {
                **source["provenance"],
                "locator": source["provenance"]["locator"]
                + f"; Issue 115 diagnostic extension at T_ref={T_REF_K} K",
            },
            "qualification": source["qualification"],
            "reference_temperature": {"magnitude": T_REF_K, "unit": "kelvin"},
            "source_sha256": source["source_sha256"],
            "value": {"magnitude": 0.0, "unit": "kelvin"},
        }
    )
    selected = epcsaft.Parameters.from_mapping(
        mapping, components=("carbon-dioxide", "water")
    )
    return selected, provenance_mutations


def load_rows() -> list[dict[str, Any]]:
    with OBSERVATIONS.open(newline="", encoding="utf-8-sig") as stream:
        rows = [row for row in csv.DictReader(stream) if row["included"] == "true"]
    return [
        {
            "row_id": row["row_id"],
            "temperature_k": float(row["temperature_k"]),
            "liquid_co2_mole_fraction": float(row["liquid_co2_mole_fraction"]),
            "observed_pressure_pa": float(row["observed_pressure_kpa"]) * 1000.0,
        }
        for row in rows
    ]


def coordinates(parameters: epcsaft.Parameters, count: int):
    return [
        regression.coordinate(
            parameters,
            [(family, ("carbon-dioxide", "water"), 1)],
            origin=0.0,
            scale=float(SCALES[index]),
            bounds=tuple(BOUNDS[:, index]),
        )
        for index, family in enumerate(
            ("k_ij", "k_ij_reciprocal_temperature_slope")[:count]
        )
    ]


def observations(parameters: epcsaft.Parameters, rows: list[dict[str, Any]]):
    return [
        regression.observation(
            parameters,
            "pressure",
            equilibrium.Problem(
                T=row["temperature_k"],
                P=equilibrium.Free(row["observed_pressure_pa"]),
                phases=[
                    equilibrium.Phase(
                        "liquid",
                        kind="liquid",
                        composition_guess=(
                            row["liquid_co2_mole_fraction"],
                            1.0 - row["liquid_co2_mole_fraction"],
                        ),
                    ),
                    equilibrium.Phase(
                        "vapor",
                        kind="vapor",
                        amount=equilibrium.Pinned(0.0),
                        composition_guess=(0.99, 0.01),
                    ),
                ],
                feed=equilibrium.Amounts(
                    {
                        "carbon-dioxide": row["liquid_co2_mole_fraction"],
                        "water": 1.0 - row["liquid_co2_mole_fraction"],
                    }
                ),
            ),
            observed=row["observed_pressure_pa"],
            form="log_ratio",
            mass_basis=False,
        )
        for row in rows
    ]


def unit_coordinates(parameters, count):
    chart = coordinates(parameters,count)
    for c in chart: c.scale = 1.0
    return chart


def physical_start(parameters, chart, values):
    point = regression.evaluate_coordinates(parameters,unit_coordinates(parameters,len(chart)),list(values),())
    native = regression.CoordinateSet(); native.coordinates = chart
    return native.start(point.parameters)


def make_evaluator(parameters, rows, identities):
    chart = coordinates(parameters, len(identities))
    items = observations(parameters, rows)

    def evaluate(values, diagnostics=True):
        point = regression.evaluate_coordinates(
            parameters, unit_coordinates(parameters,len(chart)), list(values), ()
        )
        selected = epcsaft.Parameters(
            point.parameters, parameters.component_ids, parameters.site_ids
        )
        trial = epcsaft.Mixture(selected, active_parameters=point.active)
        records, residuals, jacobian = [], [], []
        for row, item in zip(rows, items, strict=True):
            prediction = regression.predict(trial, item)
            if prediction.status != regression.ObservationStatus.Available:
                raise RuntimeError(f"{row['row_id']}: {prediction.message}")
            residuals.append(prediction.residual)
            if not diagnostics:
                continue
            derivative = regression.differentiate(trial, item)
            if derivative.status != regression.ObservationStatus.Available:
                raise RuntimeError(f"{row['row_id']}: {derivative.message}")
            jacobian.append(np.asarray(derivative.value) / prediction.value)
            records.append(
                {
                    **row,
                    "predicted_pressure_pa": prediction.value,
                    "log_pressure_residual": prediction.residual,
                    "failure": None,
                    "equilibrium_residual": prediction.convergence.residual,
                    "requested_tolerance_met": prediction.requested_tolerance_met,
                    "solver_message": prediction.message,
                    "phase_declaration": "liquid/free + vapor/pinned-zero",
                }
            )
        return np.asarray(residuals), np.asarray(jacobian), records

    return evaluate


def fit(parameters, rows, identities, origin, retained):
    chart = coordinates(parameters, len(identities))
    controls = regression.FitControls(maximum_iterations=40, maximum_elapsed_time_seconds=600.0, function_tolerance=1.0e-11,
                                      parameter_tolerance=1.0e-11, gradient_tolerance=1.0e-11)
    attempts = regression.fit(
        parameters,
        chart,
        observations(parameters, rows),
        weights=[1.0] * len(rows),
        holds=(),
        controls=controls,
        starts=[
            physical_start(parameters,chart,start) for start in (origin, retained)
        ],
        shared_iterations=40,
        shared_seconds=600.0,
    )
    usable = [
        attempt
        for attempt in attempts
        if attempt.usable
        and attempt.status == "converged"
        and attempt.covariance.projected_gradient_norm <= 1e-3
    ]
    if not usable:
        raise RuntimeError(
            f"No accepted native fit: {[(a.status, a.message) for a in attempts]}"
        )
    solved = min(usable, key=lambda attempt: attempt.final_cost)
    optimality = solved.covariance.projected_gradient_norm
    if not solved.usable or solved.status != "converged" or optimality > 1.0e-3:
        raise RuntimeError(
            f"Ceres fit gate failed: {solved.status}, {solved.message}, gradient={optimality}; failures={[(f.observation, f.message) for f in solved.failures]}"
        )
    residuals = np.asarray(solved.residuals)
    identity = unit_coordinates(parameters,len(chart))
    native = regression.CoordinateSet(); native.coordinates = identity
    point = epcsaft._core.evaluate_observations(solved.parameters,identity,native.start(solved.parameters),[],observations(parameters,rows))
    jacobian = np.asarray(point.residual_optimizer_gradient).reshape(len(rows),len(chart))
    records = [
        {
            **row,
            "predicted_pressure_pa": value,
            "log_pressure_residual": residual,
            "failure": None,
        }
        for row, value, residual in zip(
            rows, solved.predictions, residuals, strict=True
        )
    ]
    return {
        "attempts": [
            {
                "start_physical": list(start),
                "values": a.physical,
                "sse": 2 * a.final_cost,
                "status": a.status,
                "usable": a.usable,
                "projected_gradient": a.covariance.projected_gradient_norm,
                "message": a.message,
                "iterations": a.iterations,
                "trial_failures": [
                    {
                        "observation": f.observation,
                        "message": f.message,
                        "status": str(f.status),
                        "optimizer": f.optimizer,
                    }
                    for f in a.failures
                ],
            }
            for start, a in zip((origin, retained), attempts, strict=False)
        ],
        "values": np.asarray(solved.physical),
        "residuals": residuals,
        "jacobian": jacobian,
        "records": records,
        "sse": 2.0 * solved.final_cost,
        "rank": solved.covariance.rank,
        "scaled_singular_values": list(solved.covariance.singular_values),
        "nfev": solved.residual_evaluations,
        "success": solved.usable,
        "solver_status": solved.status,
        "solver_message": solved.message,
        "optimality_gate": 1.0e-3,
        "optimality": optimality,
        "active_mask": [
            int(index in solved.covariance.active_bounds) for index in range(len(chart))
        ],
        "covariance": (
            np.asarray(solved.covariance.physical) * solved.covariance.variance_factor
        ).tolist(),
        "covariance_status": solved.covariance.status,
        "variance_factor": solved.covariance.variance_factor,
        "trial_failures": [
            {
                "trial": f.trial,
                "observation": f.observation,
                "optimizer": f.optimizer,
                "status": str(f.status),
                "message": f.message,
                "jacobian": f.jacobian,
            }
            for f in solved.failures
        ],
    }


def profile(
    evaluate,
    optimum: np.ndarray,
    coordinate: int,
    sse_min: float,
    residual_count: int,
    covariance: list[float],
) -> dict[str, Any]:
    parameter_count = len(optimum)
    variance = sse_min / (residual_count - parameter_count)
    target = sse_min + 3.841458820694124 * variance
    covariance = np.asarray(covariance).reshape(parameter_count, parameter_count)
    standard_error = math.sqrt(float(covariance[coordinate, coordinate]))
    low, high = BOUNDS[:, coordinate]
    grid = np.unique(
        np.clip(
            optimum[coordinate] + standard_error * np.linspace(-3.0, 3.0, 9),
            low,
            high,
        )
    )
    points = []
    nuisance = 1 - coordinate if parameter_count == 2 else None
    for fixed in grid:
        print(f"profile coordinate {coordinate}: {fixed:.10g}", flush=True)
        values = optimum.copy()
        if nuisance is not None:
            nuisance_value = optimum[nuisance] + (
                covariance[nuisance, coordinate]
                / covariance[coordinate, coordinate]
                * (fixed - optimum[coordinate])
            )
            values[nuisance] = float(
                np.clip(nuisance_value, BOUNDS[0, nuisance], BOUNDS[1, nuisance])
            )
        values[coordinate] = fixed
        residuals = evaluate(values, diagnostics=False)[0]
        points.append(
            {
                "coordinate_value": float(fixed),
                "nuisance_value": float(values[nuisance])
                if nuisance is not None
                else None,
                "sse": float(residuals @ residuals),
                "delta_sse": float(residuals @ residuals - sse_min),
            }
        )

    crossings: list[float] = []
    for left, right in pairwise(points):
        a = left["sse"] - target
        b = right["sse"] - target
        if a == 0.0:
            crossings.append(left["coordinate_value"])
        elif a * b < 0.0:
            fraction = -a / (b - a)
            crossings.append(
                left["coordinate_value"]
                + fraction * (right["coordinate_value"] - left["coordinate_value"])
            )
    return {
        "standard_error": standard_error,
        "target_sse": target,
        "threshold_basis": (
            "nonlinear one-dimensional profile grid; chi-square(1,0.95) times SSE/(n-p)"
            if nuisance is None
            else "local quadratic profile; chi-square(1,0.95) times SSE/(n-p)"
        ),
        "crossings": crossings,
        "points": points,
    }


def metrics(records: list[dict[str, Any]]) -> dict[str, float | int]:
    residuals = [
        row["log_pressure_residual"] for row in records if row["failure"] is None
    ]
    return {
        "row_count": len(records),
        "phase_failure_count": sum(row["failure"] is not None for row in records),
        "rmse_log_pressure": math.sqrt(
            sum(value * value for value in residuals) / len(residuals)
        ),
        "mean_log_bias": sum(residuals) / len(residuals),
        "max_absolute_log_error": max(abs(value) for value in residuals),
    }


def grouped_metrics(records: list[dict[str, Any]]) -> dict[str, Any]:
    groupers = {
        "temperature_k": lambda row: f"{row['temperature_k']:.2f}",
        "observed_pressure": lambda row: (
            "lt_1_mpa"
            if row["observed_pressure_pa"] < 1.0e6
            else "1_to_5_mpa"
            if row["observed_pressure_pa"] < 5.0e6
            else "ge_5_mpa"
        ),
        "liquid_co2_mole_fraction": lambda row: (
            "lt_0.001"
            if row["liquid_co2_mole_fraction"] < 0.001
            else "0.001_to_0.01"
            if row["liquid_co2_mole_fraction"] < 0.01
            else "ge_0.01"
        ),
    }
    result: dict[str, Any] = {}
    for name, grouper in groupers.items():
        groups: dict[str, list[dict[str, Any]]] = {}
        for row in records:
            groups.setdefault(grouper(row), []).append(row)
        result[name] = {key: metrics(value) for key, value in sorted(groups.items())}
    return result


def aicc(sse: float, n: int, p: int) -> float:
    return n * math.log(sse / n) + 2.0 * p + 2.0 * p * (p + 1) / (n - p - 1)


def reproduce_profile_failure(parameters, row):
    values = (0.02, -28.0)
    point = regression.evaluate_coordinates(
        parameters, unit_coordinates(parameters, 2), list(values), ()
    )
    trial = epcsaft.Mixture(
        epcsaft.Parameters(
            point.parameters, parameters.component_ids, parameters.site_ids
        )
    )
    predicted = regression.predict(trial, observations(parameters, [row])[0])
    return {
        "row_id": row["row_id"],
        "values": {K_ID: values[0], SLOPE_ID: values[1]},
        "status": str(predicted.status),
        "message": predicted.message,
        "predicted_pressure_pa": predicted.value,
        "residual": predicted.residual,
    }


def compare(
    old, before, old_evaluations, parameters, rows, constant, temperature, result
):
    comparisons, columns, summaries = [], [], []
    names = (
        ("corrected_m_fixed_k_zero", "baseline_pressure_pa", None),
        ("constant_k", "constant_k_pressure_pa", constant),
        ("reciprocal_temperature_k", "temperature_k_pressure_pa", temperature),
    )
    for name, column, fitted in names:
        residuals, jacobian, records = old_evaluations[name]
        stages = [("historical_coordinates", residuals, jacobian, records)]
        if fitted is not None:
            ids = (K_ID,) if name == "constant_k" else (K_ID, SLOPE_ID)
            cold = make_evaluator(parameters, rows, ids)(fitted["values"])
            stages.append(("native_optimum", *cold))
        for stage, residuals, jacobian, records in stages:
            for prior, row, derivative in zip(before, records, jacobian, strict=True):
                pressure_old = float(prior[column])
                pressure_new = row["predicted_pressure_pa"]
                relative_error = abs(pressure_new / pressure_old - 1.0)
                residual_old = math.log(pressure_old / row["observed_pressure_pa"])
                comparisons.append(
                    {
                        "model": name,
                        "stage": stage,
                        "row_id": row["row_id"],
                        "old_pressure_pa": pressure_old,
                        "new_pressure_pa": pressure_new,
                        "relative_pressure_difference": relative_error,
                        "pressure_limit": 1e-7,
                        "pressure_pass": relative_error <= 1e-7,
                        "old_log_residual": residual_old,
                        "new_log_residual": row["log_pressure_residual"],
                        "absolute_residual_difference": abs(
                            row["log_pressure_residual"] - residual_old
                        ),
                        "residual_limit": 1e-7,
                        "residual_pass": abs(
                            row["log_pressure_residual"] - residual_old
                        )
                        <= 1e-7,
                        "equilibrium_residual": row["equilibrium_residual"],
                        "requested_tolerance_met": row["requested_tolerance_met"],
                        "phase_declaration": row["phase_declaration"],
                        "solver_message": row["solver_message"],
                    }
                )
                columns.append(
                    {
                        "model": name,
                        "stage": stage,
                        "row_id": row["row_id"],
                        "dlnP_dk": derivative[0],
                        "dlnP_db_per_k": derivative[1] if len(derivative) == 2 else "",
                    }
                )
            old_sse = old["models"][name].get(
                "sse_log_pressure",
                sum(
                    math.log(float(r[column]) / float(r["observed_pressure_pa"])) ** 2
                    for r in before
                ),
            )
            new_sse = (
                fitted["sse"]
                if stage == "native_optimum"
                else float(residuals @ residuals)
            )
            summaries.append(
                {
                    "model": name,
                    "stage": stage,
                    "old_sse": old_sse,
                    "new_sse": new_sse,
                    "sse_increase_limit": 1e-8 * max(1.0, old_sse),
                    "sse_pass": new_sse - old_sse <= 1e-8 * max(1.0, old_sse),
                }
            )
    for filename, table in (
        ("batch2-comparison.csv", comparisons),
        ("batch2-columns.csv", columns),
        ("batch2-objectives.csv", summaries),
    ):
        with (HERE / filename).open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=table[0], lineterminator="\n")
            writer.writeheader()
            writer.writerows(table)
    result["conversion_comparison"] = {
        "historical_base": "f4075c0f",
        "before_files": [
            "batch2-before-result.json",
            "batch2-before-predictions.csv",
            "batch2-before-profiles.csv",
        ],
        "pressure_relative_limit": 1e-7,
        "residual_absolute_limit": 1e-7,
        "objective_rule": "SSE=2*final_cost; increase <=1e-8*max(1,old SSE)",
        "pressure_misses": sum(not row["pressure_pass"] for row in comparisons),
        "residual_misses": sum(not row["residual_pass"] for row in comparisons),
        "objective_misses": sum(not row["sse_pass"] for row in summaries),
        "holdouts": "none in the retained analysis",
    }
    print(result["conversion_comparison"], flush=True)
    if any(
        result["conversion_comparison"][key]
        for key in ("pressure_misses", "residual_misses", "objective_misses")
    ):
        raise RuntimeError(
            "Conversion comparison failed; retain the comparison tables for diagnosis"
        )


def main() -> None:
    wheel = (
        ROOT / "build/environment-wheel/epcsaft-1.0.0-cp313-cp313-linux_x86_64.whl"
    )
    if not wheel.is_file():
        raise RuntimeError("current checkout wheel is missing")
    installed = distribution("epcsaft")
    site_packages = Path(installed.locate_file("")).resolve()
    module_path = Path(epcsaft.__file__).resolve()
    direct_url_path = (
        site_packages / f"epcsaft-{installed.version}.dist-info/direct_url.json"
    )
    direct_url = json.loads(direct_url_path.read_text(encoding="utf-8"))["url"]
    direct_wheel = Path(unquote(urlparse(direct_url).path)).resolve()
    if not module_path.is_relative_to(site_packages) or direct_wheel != wheel:
        raise RuntimeError("imported epcsaft is not bound to the hashed wheel")
    if (
        sha256(PARAMETERS)
        != "12e380a27eb615e25bae580f8c1ec852020879917b08fa16a5821ae76311550d"
    ):
        raise RuntimeError("materialized parameter document hash changed")
    if (
        sha256(OBSERVATIONS)
        != "ec28659dc696643c49bca771b65c79c7852b4c2e24296c9aa79b555ffdc7b9bd"
    ):
        raise RuntimeError("Kiepe input hash changed")

    parameters, provenance_mutations = load_model()
    rows = load_rows()
    if len(rows) != 39:
        raise RuntimeError(f"expected 39 retained Kiepe rows, found {len(rows)}")

    old = json.loads((HERE / "batch2-before-result.json").read_text())
    before = list(csv.DictReader((HERE / "batch2-before-predictions.csv").open()))
    old_evaluations = {}
    for name, identities in (
        ("corrected_m_fixed_k_zero", (K_ID,)),
        ("constant_k", (K_ID,)),
        ("reciprocal_temperature_k", (K_ID, SLOPE_ID)),
    ):
        print(f"historical coordinates: {name}", flush=True)
        old_evaluations[name] = make_evaluator(parameters, rows, identities)(
            tuple(old["models"][name]["values"][identity] for identity in identities)
        )
    constant_evaluate = make_evaluator(parameters, rows, (K_ID,))
    baseline_records = constant_evaluate((0.0,))[2]
    print("native constant-k fit", flush=True)
    constant = fit(
        parameters,
        rows,
        (K_ID,),
        (0.02,),
        (old["models"]["constant_k"]["values"][K_ID],),
    )
    constant["profile"] = profile(
        constant_evaluate,
        constant["values"],
        0,
        constant["sse"],
        len(rows),
        constant["covariance"],
    )
    temperature_evaluate = make_evaluator(parameters, rows, (K_ID, SLOPE_ID))
    print("native reciprocal-temperature fit", flush=True)
    temperature = fit(
        parameters,
        rows,
        (K_ID, SLOPE_ID),
        (constant["values"][0], -30.0),
        tuple(
            old["models"]["reciprocal_temperature_k"]["values"][i]
            for i in (K_ID, SLOPE_ID)
        ),
    )
    temperature["profiles"] = {
        K_ID: profile(
            temperature_evaluate,
            temperature["values"],
            0,
            temperature["sse"],
            len(rows),
            temperature["covariance"],
        ),
        SLOPE_ID: profile(
            temperature_evaluate,
            temperature["values"],
            1,
            temperature["sse"],
            len(rows),
            temperature["covariance"],
        ),
    }
    profile_failure = reproduce_profile_failure(
        parameters,
        next(row for row in rows if row["row_id"] == "kiepe-313.20-13"),
    )

    model_metrics = {
        "corrected_m_fixed_k_zero": metrics(baseline_records),
        "constant_k": metrics(constant["records"]),
        "reciprocal_temperature_k": metrics(temperature["records"]),
    }
    constant_aicc = aicc(constant["sse"], len(rows), 1)
    temperature_aicc = aicc(temperature["sse"], len(rows), 2)
    selected = "constant_k"

    result = {
        "schema_version": 1,
        "engine_version": version("epcsaft"),
        "engine_import": str(module_path.relative_to(ROOT)),
        "installed_wheel": str(wheel.relative_to(ROOT)),
        "installed_wheel_sha256": sha256(wheel),
        "installed_distribution_direct_url_wheel": str(direct_wheel.relative_to(ROOT)),
        "parameter_document": str(PARAMETERS.relative_to(ROOT)),
        "parameter_document_sha256": sha256(PARAMETERS),
        "observations": str(OBSERVATIONS.relative_to(ROOT)),
        "observations_sha256": sha256(OBSERVATIONS),
        "co2_segment_count": 2.0729,
        "historical_segment_count_not_executed": 2.079,
        "candidate_provenance_normalization": {
            "status": "analysis-only; physical values and topology unchanged",
            "mutation_count": len(provenance_mutations),
            "field_mutations": provenance_mutations,
            "resolved_candidate_temperature_k": [293.15, 393.17],
        },
        "temperature_extension": "393.15 K active ceiling extended to 393.17 K only to execute the six retained top-temperature rows",
    }
    result["native_declarations"] = {
        "observation_row_count": len(rows),
        "phases": [
            {
                "name": "liquid",
                "kind": "liquid",
                "amount": "free",
                "composition_guess": "row liquid feed",
            },
            {
                "name": "vapor",
                "kind": "vapor",
                "amount": 0.0,
                "composition_guess": [0.99, 0.01],
            },
        ],
        "feed_basis": "1 mol, carbon-dioxide = row liquid mole fraction, water = 1-x",
        "pressure": "free, guess = row observed pressure Pa",
        "residual": "EquilibriumPressure LogRatio, unit weights, holds=()",
        "coordinates": {
            "origins": [0.0, 0.0],
            "scales": SCALES.tolist(),
            "bounds": BOUNDS.tolist(),
            "reference_temperature_k": T_REF_K,
        },
        "starts_physical": {
            "constant_k": [
                [0.02],
                list(old["models"]["constant_k"]["values"].values()),
            ],
            "reciprocal_temperature_k": [
                [float(constant["values"][0]), -30.0],
                list(old["models"]["reciprocal_temperature_k"]["values"].values()),
            ],
        },
        "controls": {
            "function_tolerance": 1e-11,
            "gradient_tolerance": 1e-11,
            "parameter_tolerance": 1e-11,
            "maximum_iterations": 40,
            "shared_iterations": 40,
            "shared_seconds": 600.0,
        },
        "ionic_dispersion_admission": "exclude-same-sign-ion-pairs added only for current packet admission; native selection removes ionic choices from the uncharged binary",
    }
    write_json(HERE / "input-result.json", result)

    result = {
        "schema_version": 1,
        "status": "completed_with_nonlinear_profile_gate_unresolved"
        if all(value["phase_failure_count"] == 0 for value in model_metrics.values())
        else "failed_rows",
        "selected_model": selected,
        "selection_rule": "minimum supported complexity; reciprocal-temperature candidate cannot be promoted without completed two-sided nonlinear profiles",
        "bounds": {K_ID: [-1.0, 1.0], SLOPE_ID: [-1000.0, 1000.0]},
        "reference_temperature_k": T_REF_K,
        "effective_k_definition": "k(T) = k_ref + slope_K * (1/T - 1/T_ref)",
        "models": {
            "corrected_m_fixed_k_zero": {
                "values": {K_ID: 0.0},
                "metrics": model_metrics["corrected_m_fixed_k_zero"],
            },
            "constant_k": {
                "scientific_status": "best_supported_parameterization",
                "values": {K_ID: float(constant["values"][0])},
                "metrics": model_metrics["constant_k"],
                "metrics_by_axis": grouped_metrics(constant["records"]),
                "sse_log_pressure": constant["sse"],
                "aicc": constant_aicc,
                "rank": constant["rank"],
                "scaled_singular_values": constant["scaled_singular_values"],
                "active_mask": constant["active_mask"],
                "nfev": constant["nfev"],
                "solver_success": constant["success"],
                "solver_status": constant["solver_status"],
                "solver_message": constant["solver_message"],
                "optimality_gate": constant["optimality_gate"],
                "optimality": constant["optimality"],
                "nonlinear_profile": constant["profile"],
            },
            "reciprocal_temperature_k": {
                "scientific_status": "candidate_profile_gate_unresolved",
                "values": {
                    K_ID: float(temperature["values"][0]),
                    SLOPE_ID: float(temperature["values"][1]),
                },
                "metrics": model_metrics["reciprocal_temperature_k"],
                "metrics_by_axis": grouped_metrics(temperature["records"]),
                "sse_log_pressure": temperature["sse"],
                "aicc": temperature_aicc,
                "rank": temperature["rank"],
                "scaled_singular_values": temperature["scaled_singular_values"],
                "active_mask": temperature["active_mask"],
                "nfev": temperature["nfev"],
                "solver_success": temperature["success"],
                "solver_status": temperature["solver_status"],
                "solver_message": temperature["solver_message"],
                "optimality_gate": temperature["optimality_gate"],
                "optimality": temperature["optimality"],
                "local_quadratic_uncertainty": temperature["profiles"],
                "nonlinear_profile_gate": {
                    "status": "unresolved",
                    "current_trial": profile_failure,
                    "historical_failure_trial": old["models"][
                        "reciprocal_temperature_k"
                    ]["nonlinear_profile_gate"]["reproduced_failure_trial"],
                    "scope": "two-sided nuisance-refitted nonlinear profiles were not completed",
                    "interpretation": "local quadratic crossings are not finite nonlinear profile bounds and are not used for promotion",
                },
            },
        },
        "supported_data_domain": {
            "temperature_k": [
                min(row["temperature_k"] for row in rows),
                max(row["temperature_k"] for row in rows),
            ],
            "observed_pressure_pa": [
                min(row["observed_pressure_pa"] for row in rows),
                max(row["observed_pressure_pa"] for row in rows),
            ],
            "liquid_co2_mole_fraction": [
                min(row["liquid_co2_mole_fraction"] for row in rows),
                max(row["liquid_co2_mole_fraction"] for row in rows),
            ],
            "active_domain_caveat": "393.17 K rows are 0.02 K above the active candidate ceiling and remain boundary evidence only",
        },
    }
    for name, fitted in (
        ("constant_k", constant),
        ("reciprocal_temperature_k", temperature),
    ):
        result["models"][name]["conditional_covariance"] = {
            "physical": fitted["covariance"],
            "status": fitted["covariance_status"],
            "variance_factor": fitted["variance_factor"],
            "basis": "native exact columns, SSE/(n-rank)",
        }
        result["models"][name]["native_starts"] = fitted["attempts"]
        result["models"][name]["trial_failures"] = fitted["trial_failures"]
        result["models"][name]["solver"] = "native Ceres"
    result["claim_changes"] = []
    if profile_failure["predicted_pressure_pa"] is not None:
        result["claim_changes"].append(
            "The historical multiple_unresolved_coexistence_branches trial at k=0.02, b=-28 K now returns a locally converged declared-phase state; this does not establish global branch uniqueness. The two-sided nuisance-refitted nonlinear profile remains uncompleted."
        )
    compare(
        old, before, old_evaluations, parameters, rows, constant, temperature, result
    )
    write_json(HERE / "result.json", result)

    with (HERE / "predictions.csv").open("w", newline="", encoding="utf-8") as stream:
        fields = [
            "row_id",
            "temperature_k",
            "liquid_co2_mole_fraction",
            "observed_pressure_pa",
            "baseline_pressure_pa",
            "constant_k_pressure_pa",
            "temperature_k_pressure_pa",
            "temperature_k_log_residual",
            "phase_failure",
        ]
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for baseline, fixed, variable in zip(
            baseline_records, constant["records"], temperature["records"], strict=True
        ):
            writer.writerow(
                {
                    "row_id": variable["row_id"],
                    "temperature_k": f"{variable['temperature_k']:.2f}",
                    "liquid_co2_mole_fraction": f"{variable['liquid_co2_mole_fraction']:.8g}",
                    "observed_pressure_pa": f"{variable['observed_pressure_pa']:.10g}",
                    "baseline_pressure_pa": f"{baseline['predicted_pressure_pa']:.10g}",
                    "constant_k_pressure_pa": f"{fixed['predicted_pressure_pa']:.10g}",
                    "temperature_k_pressure_pa": f"{variable['predicted_pressure_pa']:.10g}",
                    "temperature_k_log_residual": f"{variable['log_pressure_residual']:.12g}",
                    "phase_failure": variable["failure"] or "",
                }
            )

    with (HERE / "profiles.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=(
                "identity",
                "coordinate_value",
                "nuisance_value",
                "sse",
                "delta_sse",
                "target_sse",
            ),
            lineterminator="\n",
        )
        writer.writeheader()
        for identity, evidence in temperature["profiles"].items():
            for point in evidence["points"]:
                writer.writerow(
                    {
                        "identity": identity,
                        **point,
                        "target_sse": evidence["target_sse"],
                    }
                )

    print(
        {
            name: {"values": model["values"], "sse": model["sse_log_pressure"]}
            for name, model in result["models"].items()
            if "sse_log_pressure" in model
        },
        flush=True,
    )


def declarations():
    return [(p, coordinates(p, 2), observations(p, load_rows())) for p in [load_model()[0]]]


if __name__ == "__main__":
    main()
