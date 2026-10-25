"""Replay the retained propane saturation fit through the installed wheel (#49 R4).

Training rows: inputs/source/propane-baseline.csv (Glos 2004 Table 2 at 210, 270, 330 K; saturation
pressure and saturated-liquid mass density, log-ratio residuals, unit weights). Coordinates, starts and
bounds: inputs/propane-baseline-v3.fit.json (three families, four multistart points). Holdout and stress rows:
the Glos 2004 saturation table retained by analyses/reproduction/2004-glos (roles held_out and stress),
compared as MODERN_REFIT evidence with their expanded uncertainties; no accuracy bar is claimed.
Run from the repository root: .venv/bin/python analyses/2026-pure-saturation-regression/scripts/replay_propane_fit.py
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import sys
from datetime import date
from pathlib import Path

from epcsaft import Parameters, Pure, regression
from epcsaft.equilibrium import Amounts, Free, Phase, Pinned, Problem

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parents[1]
FIT = HERE / "inputs/propane-baseline-v3.fit.json"
TRAINING = HERE / "inputs/source/propane-baseline.csv"
PARAMETERS = ROOT / "data/packets/gross-2001-propane/3/parameters/parameters.json"
GLOS = ROOT / "analyses/reproduction/2004-glos/inputs/source/glos-2004-propane-saturation.csv"
FAMILY = {"dispersion_energy_over_k": "dispersion_energy_over_k", "segment_count": "segment_count", "segment_diameter": "segment_diameter"}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def saturation(T: float, guess: float) -> Problem:
    return Problem(phases=[Phase("L", kind="liquid", composition_guess=[1.0]), Phase("V", kind="vapor", amount=Pinned(0.0), composition_guess=[1.0])],
                   T=T, P=Free(guess), feed=Amounts({"propane": 1.0}), row_tolerance=1.0e-10)


def observations(parameters: Parameters, rows: list[dict[str, str]]) -> list[regression.Observation]:
    result = []
    for row in rows:
        T, p, rho = float(row["temperature_k"]), float(row["saturation_pressure"]), float(row["saturated_liquid_mass_density"])
        result.append(regression.observation(parameters, "pressure", saturation(T, p), observed=p, form="log_ratio", mass_basis=False))
        result.append(regression.observation(parameters, "phase_density", saturation(T, p), observed=rho, phase="L", form="log_ratio", mass_basis=True))
    return result


def evaluate(model: Pure, parameters: Parameters, T: float, p_source: float) -> tuple[float | None, float | None, str]:
    p = regression.predict(model, regression.observation(parameters, "pressure", saturation(T, p_source), observed=p_source, form="log_ratio", mass_basis=False))
    if p.status != regression.ObservationStatus.Available or p.value is None:
        return None, None, f"{p.status.name}: {p.message}"
    rho = regression.predict(model, regression.observation(parameters, "phase_density", saturation(T, p_source), observed=1.0, phase="L", form="log_ratio", mass_basis=True))
    if rho.status != regression.ObservationStatus.Available:
        return p.value, None, f"{rho.status.name}: {rho.message}"
    return p.value, rho.value, p.status.name


def main() -> None:
    request = json.loads(FIT.read_text())
    digest = request["parameters"]["sha256"].removeprefix("sha256:")
    if digest != sha(PARAMETERS):
        sys.exit(f"parameter document digest mismatch: {digest} != {sha(PARAMETERS)}")
    document = Parameters.from_json(PARAMETERS, components=["propane"])
    training = list(csv.DictReader(TRAINING.open()))
    glos = list(csv.DictReader(GLOS.open()))
    identities = request["active"]["identities"]
    starts = request["solver"]["multistart"]["starts"]
    coordinates = []
    for identity in identities:
        spec = request["active"]["coordinates"][identity]
        coordinates.append(regression.coordinate(document, [(FAMILY[identity.rsplit("/", 1)[1]], "propane", 1)], origin=spec["origin"],
                                                 scale=spec["scale"], bounds=tuple(spec["bounds"])))
    import copy
    chart = regression.CoordinateSet(); chart.coordinates = coordinates
    if '--declare-only' in sys.argv:
        observations(document, training)
        return
    optimizer_starts = []
    for start in starts:
        declared = copy.deepcopy(json.loads(PARAMETERS.read_text()))
        magnitudes = dict(zip(identities,start["values"]))
        for component in declared["components"]:
            for coefficient in component["coefficients"]:
                if coefficient["identity"] in magnitudes: coefficient["value"]["magnitude"] = magnitudes[coefficient["identity"]]
        optimizer_starts.append(chart.start(Parameters.from_mapping(declared,components=["propane"])._native,observations(document,training)))
    controls = regression.FitControls(maximum_iterations=100, maximum_elapsed_time_seconds=600.0)
    results = regression.fit(document, coordinates, observations(document, training), starts=optimizer_starts, controls=controls,
                             shared_iterations=len(starts) * controls.maximum_iterations,
                             shared_seconds=len(starts) * controls.maximum_elapsed_time_seconds)
    if len(results) != len(starts):
        raise RuntimeError("native multistart exhausted its budget before all declared starts")
    fits = list(zip((start["identity"] for start in starts), results))
    reference = fits[0][1]
    rows, failures = [], []
    for label, result in fits:
        if not result.usable:
            failures.append({"start": label, "status": result.status, "message": result.message,
                             "trial_failures": [(item.observation, item.status.name, item.message) for item in result.failures]})
    if not reference.usable:
        raise RuntimeError(f"first declared propane start has no usable reference: {reference.status}: {reference.message}")
    model = Pure(regression.fitted_parameters(document, reference))
    packet = Pure(document)
    for row in glos:
        T, p_src, rho_src = float(row["T_K"]), float(row["p_sat_Pa"]), float(row["rho_sat_liq_kg_m3"])
        p_fit, rho_fit, status = evaluate(model, document, T, p_src)
        p_doc, rho_doc, doc_status = evaluate(packet, document, T, p_src)
        rows.append({"row_id": row["row_id"], "role": row["role"], "T_K": T, "p_source_Pa": p_src, "p_expanded_uncertainty_Pa": row["p_sat_expanded_uncertainty_Pa"],
                     "rho_source_kg_m3": rho_src, "rho_expanded_uncertainty_kg_m3": row["rho_sat_liq_expanded_uncertainty_kg_m3"],
                     "p_fitted_Pa": p_fit, "rho_fitted_kg_m3": rho_fit, "fitted_status": status,
                     "p_fitted_relative_error": None if p_fit is None else p_fit / p_src - 1.0,
                     "rho_fitted_relative_error": None if rho_fit is None else rho_fit / rho_src - 1.0,
                     "p_packet_Pa": p_doc, "rho_packet_kg_m3": rho_doc, "packet_status": doc_status,
                     "p_packet_relative_error": None if p_doc is None else p_doc / p_src - 1.0,
                     "rho_packet_relative_error": None if rho_doc is None else rho_doc / rho_src - 1.0})
        if p_fit is None or rho_fit is None or p_doc is None or rho_doc is None:
            failures.append({"row_id": row["row_id"], "fitted_status": status, "packet_status": doc_status})
    stamp = date.today().isoformat()
    out_csv = HERE / f"results/propane-replay-{stamp}.csv"
    with out_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    preparation = json.loads((ROOT / "build/environment-wheel/preparation.json").read_text())
    variance = reference.covariance.variance_factor
    physical_std = None
    if reference.covariance.status == "available" and variance is not None:
        n = len(reference.physical)
        physical_std = [math.sqrt(variance * reference.covariance.physical[k * n + k]) for k in range(n)]
    summary = {
        "date": stamp, "class": "MODERN_REFIT", "wheel_sha256": preparation["wheel_sha256"],
        "sources": {"training_csv": str(TRAINING.relative_to(ROOT)), "training_sha256": sha(TRAINING), "fit_request_sha256": sha(FIT),
                    "parameters_sha256": sha(PARAMETERS), "holdout_csv": str(GLOS.relative_to(ROOT)), "holdout_sha256": sha(GLOS)},
        "objective": "sum of ln(y_pred / y_obs)^2 over 3 pressures and 3 saturated-liquid mass densities, unit weights",
        "reference_start": starts[0]["identity"], "optimizer_starts": optimizer_starts,
        "coordinates": {identity: request["active"]["coordinates"][identity] for identity in identities},
        "starts": [{"start": label, "values": start["values"], "status": result.status, "iterations": result.iterations,
                    "residual_evaluations": result.residual_evaluations, "jacobian_evaluations": result.jacobian_evaluations,
                    "trial_failures": result.trial_failures, "initial_cost": result.initial_cost, "final_cost": result.final_cost if result.usable else None,
                    "sse": 2.0 * result.final_cost if result.usable else None,
                    "covariance_status": result.covariance.status, "rank": result.covariance.rank,
                    "singular_values": result.covariance.singular_values, "active_bounds": result.covariance.active_bounds,
                    "physical_covariance": result.covariance.physical, "variance_factor": result.covariance.variance_factor,
                    "projected_gradient_norm": result.covariance.projected_gradient_norm,
                    "fitted": result.physical if result.usable else None} for (label, result), start in zip(fits, starts)],
        "multistart_spread": {"epsilon_K": max(r.physical[0] for _, r in fits if r.usable) - min(r.physical[0] for _, r in fits if r.usable),
                              "m": max(r.physical[1] for _, r in fits if r.usable) - min(r.physical[1] for _, r in fits if r.usable),
                              "sigma_angstrom": max(r.physical[2] for _, r in fits if r.usable) - min(r.physical[2] for _, r in fits if r.usable)},
        "fitted": {"identities": identities, "values": reference.physical, "units": ["K", "1", "angstrom"],
                   "packet_values": [208.11, 2.002, 3.6184], "historical_fitted_2026_07": [209.27852823846877, 1.9835426308533168, 3.627861545037305],
                   "historical_note": "the retained July result fitted eleven Glos rows from the packet start with another optimizer build; agreement is descriptive"},
        "covariance": {"status": reference.covariance.status, "assumption": reference.covariance.assumption, "rank": reference.covariance.rank,
                       "singular_values": reference.covariance.singular_values, "active_bounds": reference.covariance.active_bounds,
                       "variance_factor": variance, "physical": reference.covariance.physical,
                       "physical_std_scaled_by_variance_factor": physical_std,
                       "meaning": "rank and the conditional covariance describe local conditioning of the least-squares problem at the fitted point under the stated assumption; six residuals and three parameters leave three degrees of freedom, so the variance factor is a rough scale, not a validated uncertainty; neither establishes structural or global identifiability nor predictive accuracy"},
        "training_residuals": {"weighted": reference.weighted_residuals, "statuses": [item.name for item in reference.statuses]},
        "holdout": {"rows": len(rows), "csv": str(out_csv.relative_to(ROOT)),
                    "role_note": "roles come from the 2004-glos packet; its 150 K 'training' row is outside this replay's training set (210, 270, 330 K) and is a genuine holdout here",
                    "max_abs_relative_error_held_out_pressure": max(abs(r["p_fitted_relative_error"]) for r in rows if r["role"] == "held_out" and r["p_fitted_relative_error"] is not None),
                    "max_abs_held_out_pressure_error_in_expanded_uncertainties": max(abs(r["p_fitted_Pa"] - r["p_source_Pa"]) / float(r["p_expanded_uncertainty_Pa"]) for r in rows if r["role"] == "held_out" and r["p_fitted_Pa"] is not None and float(r["p_expanded_uncertainty_Pa"]) > 0),
                    "max_abs_relative_error_held_out_density": max(abs(r["rho_fitted_relative_error"]) for r in rows if r["role"] == "held_out" and r["rho_fitted_relative_error"] is not None),
                    "max_abs_held_out_density_error_in_expanded_uncertainties": max(abs(r["rho_fitted_kg_m3"] - r["rho_source_kg_m3"]) / float(r["rho_expanded_uncertainty_kg_m3"]) for r in rows if r["role"] == "held_out" and r["rho_fitted_kg_m3"] is not None and float(r["rho_expanded_uncertainty_kg_m3"]) > 0)},
        "failures": failures,
        "claim_limits": ["the converged multistart spread measures optimizer reproducibility across those starts; rank and the conditional covariance (physical_std_scaled_by_variance_factor) describe local conditioning under the stated assumptions; neither establishes structural or global identifiability of the three parameters nor predictive accuracy",
                         "unit weights on log-ratio residuals treat every training value as equally uncertain in relative terms; the source expanded uncertainties were not used as weights and the covariance assumption therefore means unit variance in log space",
                         "held-out rows are descriptive comparisons against Glos 2004 with their expanded uncertainties; no accuracy bar or promotion",
                         "the stress row at 110 K lies outside the fitted domain and is retained whatever its error",
                         "no held-out predictive claim transfers to other fluids, mixtures or properties"],
    }
    (HERE / f"results/propane-replay-{stamp}.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps({k: summary[k] for k in ("fitted", "multistart_spread", "holdout", "failures")}, indent=1))
    print("covariance", summary["covariance"]["status"], summary["covariance"]["singular_values"], physical_std)


def declarations():
    return main()


if __name__ == "__main__":
    main()
