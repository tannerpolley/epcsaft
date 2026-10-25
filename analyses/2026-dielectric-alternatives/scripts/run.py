"""Issue 27: common-data comparison of the admitted permittivity alternatives.

Run from the repository root with the installed wheel, after selections.py:
    .venv/bin/python analyses/2026-dielectric-alternatives/scripts/run.py

Six variants of the Figiel 2025 water/Na+/Cl- model at 298.15 K and 1 bar:
F0/D0/P0 keep every Figiel parameter and exchange only epsilon_r
(ion-fraction suppression, density-temperature, dipolar perturbation);
F1/D1/P1 refit one shared water-Na+ = water-Cl- k_ij to Hamer-Wu NaCl MIAC
with m <= 2 mol/kg (Rueben's relative loss). Rows with m > 2 validate.
"""

from __future__ import annotations

import copy
import csv
import hashlib
import json
import math
from importlib import metadata
from pathlib import Path

from epcsaft import Mixture, Parameters, _core, regression

ROOT = Path(__file__).resolve().parents[3]
LEAF = ROOT / "analyses/2026-dielectric-alternatives"
HAMER_WU = ROOT / "data/packets/hamer-wu-1972-aqueous-alkali-halides/1/observations/hamer-wu-1972-aqueous-alkali-halides.csv"
LEDGER = ROOT / "analyses/reproduction/2025-figiel/inputs/source/figiel-2025-regression-target-ledger.csv"
NIST = ROOT / "analyses/2026-shared-derivative-qualification/inputs/nist-srd69-solvent-isobars.tsv"
T, P = 298.15, 1.0e5
COMPONENTS = ("water", "sodium-cation", "chloride-anion")
BUDGET = regression.FitControls(maximum_iterations=100, maximum_elapsed_time_seconds=600.0)
FORMS = {"F": "ion-fraction-suppression", "D": "density-temperature", "P": "dipolar-perturbation"}
CALIBRATION_MAX_MOLALITY = 2.0
FIGIEL_K = -0.3  # figiel-yu-held-2025/3 water-Na+ and water-Cl- k_ij (equal)


def selection(form: str, k: float) -> Mixture:
    document = copy.deepcopy(json.loads((LEAF / f"inputs/selections/{FORMS[form]}.json").read_text()))
    for pair in document["pairs"]:
        if "water" in (pair["component_id_a"], pair["component_id_b"]):
            next(c for c in pair["coefficients"] if c["family"] == "k_ij")["value"]["magnitude"] = k
    return Mixture(Parameters.from_mapping(document, components=COMPONENTS))


def water_molar_mass() -> float:
    document = json.loads((LEAF / "inputs/selections/ion-fraction-suppression.json").read_text())
    water = next(c for c in document["components"] if c["component_id"] == "water")
    return water["fixed"]["molar_mass"]["value"]["magnitude"] * 1.0e-3


def coordinate(parameters: Parameters):
    return regression.coordinate(parameters, [("k_ij", target, 1) for target in (('water', ion) for ion in COMPONENTS[1:])],
                                 origin=FIGIEL_K, scale=1.0, bounds=(-0.8, 0.4))


def electrolyte(model: Mixture, kind: str, *, observed: float, m: float = 0.0,
                x_ion: float = 0.0, species: str | None = None):
    return regression.electrolyte(
        model.parameters, kind, T, P,
        dict(zip(COMPONENTS, (1.0 - x_ion, x_ion / 2, x_ion / 2))),
        ions=COMPONENTS[1:], stoichiometry=(1, 1), molality=m, species=species,
        observed=observed, form="difference", scale=observed if kind == "mean_ionic_activity_coefficient" else 1.0)


def prediction(model: Mixture, item):
    result = regression.predict(model, item)
    if result.status != regression.ObservationStatus.Available:
        raise RuntimeError(f"{result.status.name}: {result.message}")
    return result


def observations() -> dict[str, list[dict]]:
    with HAMER_WU.open() as stream:
        miac_rows = [{"m": float(r["molality_mol_kg"]), "value": float(r["gamma_pm_m"])}
                     for r in csv.DictReader(stream) if r["salt"] == "NaCl"]
    with LEDGER.open() as stream:
        ledger = list(csv.DictReader(stream))
    dielectric = [{"x_ion": float(r["composition_value"]), "value": float(r["value"]), "source": r["underlying_reference"]}
                  for r in ledger if r["quantity"] == "relative_dielectric_constant" and r["species"] == "NaCl"]
    solvation = []
    for ion, component in (("Na+", "sodium-cation"), ("Cl-", "chloride-anion")):
        rows = [r for r in ledger if r["quantity"] == "single_ion_solvation_gibbs_energy" and r["species"] == ion
                and r["solvent"] == "water"]
        values = [float(r["value"]) for r in rows if r["role"] == "underlying_literature_value"]
        average = next(float(r["value"]) for r in rows if r["role"].startswith("reported_average"))
        solvation.append({"ion": ion, "component": component, "value": average,
                          "half_spread": (max(values) - min(values)) / 2})
    return {"miac": miac_rows, "dielectric": dielectric, "solvation": solvation}


def attempt(calculation):
    """Return (value, status); a failed calculation is retained as a row, never dropped."""
    try:
        return calculation(), "ok"
    except Exception as error:  # every Engine refusal is evidence for the comparison
        return None, f"failed: {type(error).__name__}: {error}"


def refit(form: str, rows: list[dict]) -> dict:
    model = selection(form, FIGIEL_K)
    items = [electrolyte(model, "mean_ionic_activity_coefficient", observed=r["value"], m=r["m"]) for r in rows]
    weights = [1.0 if r["m"] <= CALIBRATION_MAX_MOLALITY else 0.0 for r in rows]
    holdout = regression.CheckGroup()
    holdout.name, holdout.rows = "MIAC above 2 mol/kg", [i for i, w in enumerate(weights) if w == 0]
    holdout.metric, holdout.report_only = regression.CheckMetric.AbsoluteDifference, True
    result, status = attempt(lambda: regression.fit(
        model.parameters, [coordinate(model.parameters)], items, weights=weights,
        checks=[holdout], holds=(), starts=[[0.0]], controls=BUDGET)[0])
    if result is None:
        return {"form": form, "calibration_rows": sum(w > 0 for w in weights), "k_water_ion": None,
                "sigma_k": None, "relative_loss": None, "at_bound": None, "status": status}
    covariance_available = result.covariance.status == "available" and result.covariance.variance_factor is not None
    k = result.physical[0] if result.physical else None
    print(f"{form}: k={k}, relative SSE={2 * result.final_cost:.12g}", flush=True)
    return {"form": form, "calibration_rows": sum(w > 0 for w in weights), "k_water_ion": k,
            "sigma_k": math.sqrt(result.covariance.physical[0] * result.covariance.variance_factor) if covariance_available else None,
            "uncertainty_status": result.covariance.status, "relative_loss": 2 * result.final_cost,
            "at_bound": bool(result.covariance.active_bounds),
            "status": "ok" if result.usable else f"failed: {result.status}: {result.message}",
            "rank": result.covariance.rank, "singular_value": result.covariance.singular_values[0] if result.covariance.singular_values else None,
            "variance_factor": result.covariance.variance_factor,
            "projected_gradient_norm": result.covariance.projected_gradient_norm,
            "fit_status": result.status, "fit_message": result.message,
            "iterations": result.iterations, "trial_failures": result.trial_failures,
            "failures": [{"trial": f.trial, "row": f.observation, "optimizer": f.optimizer,
                          "status": f.status.name, "message": f.message} for f in result.failures],
            "holdout_check": {"maximum": result.checks[0].maximum, "mean": result.checks[0].mean,
                              "failed_rows": result.checks[0].failed_rows} if result.checks else None,
            "covariance_assumption": result.covariance.assumption,
            "native_failed_fit_predictions": result.predictions if not result.usable else [],
            "native_failed_fit_residuals": result.residuals if not result.usable else []}


def rows_for(variant: str, form: str, k: float, sigma_k: float | None, data: dict, water_density: float) -> list[dict]:
    fitted_miac = variant == "F0" or variant.endswith("1")
    model = selection(form, k)
    base = {"variant": variant, "form": FORMS[form], "k_water_ion": k, "sigma_k": sigma_k}
    rows = []
    chart = regression.CoordinateSet()
    chart.coordinates = [coordinate(model.parameters)]
    for r in data["miac"]:
        item = electrolyte(model, "mean_ionic_activity_coefficient", observed=r["value"], m=r["m"])
        result, miac_status = attempt(lambda: prediction(model, item))
        g = result.value if result else None
        # Native Difference column; the independent k perturbation below is unit-scale only.
        gradient, derivative_status = attempt(lambda: _core.evaluate_observations(
            model.parameters._native,chart.coordinates,chart.start(model.parameters._native),[],[item])) if sigma_k and result else (None, "not_requested")
        if gradient is not None and any(status != regression.ObservationStatus.Available for status in gradient.statuses):
            derivative_status = 'failed: ' + '; '.join(status.name for status in gradient.statuses)
        slope = gradient.residual_optimizer_gradient[0] if gradient is not None and derivative_status == "ok" else None
        difference, difference_error, difference_tolerance = None, None, None
        verification_status = "not_requested"
        if slope is not None:
            step = 1.0e-4
            shifted, verification_status = attempt(lambda: (
                prediction(selection(form, k + step), item).value,
                prediction(selection(form, k - step), item).value))
            if shifted is not None:
                difference = (shifted[0] - shifted[1]) / (2 * step) / r["value"]
                difference_error = abs(slope - difference)
                difference_tolerance = 2e-5 * abs(slope) + 1e-7
                if difference_error > difference_tolerance:
                    verification_status = "failed: exact tied derivative differs from central difference"
        fitted = fitted_miac and (variant == "F0" or r["m"] <= CALIBRATION_MAX_MOLALITY)
        rows.append({**base, "property": "miac", "species": "NaCl", "condition": "molality_mol_kg",
                     "condition_value": r["m"], "reference_kind": "Hamer-Wu 1972", "observed": r["value"],
                     "model": g, "fitted_to_row": fitted, "relative_error_k_slope": slope,
                     "reference_uncertainty": None, "branch": "reference-path",
                     "status": miac_status, "relative_error_k_slope_difference": difference,
                     "relative_error_k_slope_difference_error": difference_error,
                     "relative_error_k_slope_difference_tolerance": difference_tolerance,
                     "derivative_status": derivative_status, "verification_status": verification_status})
    for r in data["dielectric"]:
        item = electrolyte(model, "bulk_relative_permittivity", observed=r["value"], x_ion=r["x_ion"])
        result, status = attempt(lambda: prediction(model, item))
        rows.append({**base, "property": "relative_permittivity", "species": "NaCl", "condition": "x_ion",
                     "condition_value": r["x_ion"], "reference_kind": r["source"], "observed": r["value"],
                     "model": result.value if result else None, "fitted_to_row": form == "F",
                     "relative_error_k_slope": None, "reference_uncertainty": None,
                     "branch": "liquid" if result else "", "status": status})
    for r in data["solvation"]:
        item = electrolyte(model, "solvation_gibbs_energy", observed=r["value"] * 1000, species=r["component"])
        result, status = attempt(lambda: prediction(model, item))
        rows.append({**base, "property": "solvation_gibbs_kJ_mol", "species": r["ion"], "condition": "x_ion",
                     "condition_value": 0.0, "reference_kind": "Figiel SI S5 average", "observed": r["value"],
                     "model": result.value / 1000.0 if result else None, "fitted_to_row": form == "F",
                     "relative_error_k_slope": None, "reference_uncertainty": r["half_spread"],
                     "branch": "reference-path", "status": status})
    mw = water_molar_mass()
    for m in (0.0, 1.0, 2.0, 4.0, 6.0):
        x_ion = 2 * m * mw / (1 + 2 * m * mw)
        item = regression.density(model.parameters, T, P,
                                  dict(zip(COMPONENTS, (1 - x_ion, x_ion / 2, x_ion / 2))),
                                  phase="liquid", observed=water_density, form="difference", scale=1.0, mass_basis=False)
        result, status = attempt(lambda: prediction(model, item))
        rows.append({**base, "property": "liquid_molar_density_mol_m3", "species": "NaCl", "condition": "molality_mol_kg",
                     "condition_value": m, "reference_kind": "NIST SRD 69 IAPWS-95" if m == 0.0 else "F0 model",
                     "observed": water_density if m == 0.0 else None,
                     "model": result.value if result else None, "fitted_to_row": False,
                     "relative_error_k_slope": None, "reference_uncertainty": None,
                     "branch": "liquid" if result else "", "status": status,
                     "native_density_branch": result.branch.name if result else "",
                     "stable_root_count": result.stable_root_count if result else None})
    for row in rows:
        known = row["observed"] is not None and row["model"] is not None
        row["error"] = row["model"] - row["observed"] if known else None
        row["relative_error"] = row["error"] / abs(row["observed"]) if known else None
    return rows


def summarize(rows: list[dict]) -> dict:
    out: dict = {}
    for row in rows:
        partition = "fitted" if row["fitted_to_row"] else "unfitted"
        if row["property"] == "liquid_molar_density_mol_m3":
            partition = "vs_iapws_pure_water" if row["condition_value"] == 0.0 else "difference_from_same_fit_F"
        elif row["property"] == "miac":
            partition += "_m_le_2" if row["condition_value"] <= CALIBRATION_MAX_MOLALITY else "_m_gt_2"
        out.setdefault(f"{row['variant']}/{row['property']}/{partition}", []).append(row)
    metrics = {}
    for key, group in sorted(out.items()):
        v = [r for r in group if r["status"] == "ok" and r["relative_error"] is not None]
        if not v:
            metrics[key] = {"n": 0, "failed": len(group)}
            continue
        slopes = [r["relative_error_k_slope"] for r in v]
        metrics[key] = {
            "n": len(v), "failed": sum(r["status"] != "ok" for r in group),
            "mard": sum(abs(r["relative_error"]) for r in v) / len(v),
            "max_ard": max(abs(r["relative_error"]) for r in v),
            "mean_abs_error": sum(abs(r["error"]) for r in v) / len(v),
            # First-order propagation of sigma_k through the compared metric (MARD) itself.
            "mard_sigma_from_k": (v[0]["sigma_k"] * abs(sum(math.copysign(1.0, r["relative_error"]) * d
                                                             for r, d in zip(v, slopes))) / len(v)
                                  if all(d is not None and r["derivative_status"] == "ok" and r["verification_status"] == "ok"
                                         for r, d in zip(v, slopes)) else None)}
    return metrics


def main() -> None:
    data = observations()
    water_density = next(float(r["mass_density_kg_m3"]) for r in csv.DictReader(NIST.open(), delimiter="\t")
                         if r["solvent"] == "water" and float(r["T_K"]) == T) / water_molar_mass()
    fits = [refit(form, data["miac"]) for form in FORMS]
    rows = [row for form in FORMS for row in rows_for(f"{form}0", form, FIGIEL_K, None, data, water_density)]
    rows += [row for fit in fits if fit["status"] == "ok" for row in
             rows_for(f"{fit['form']}1", fit["form"], fit["k_water_ion"], fit["sigma_k"], data, water_density)]
    reference = {(r["variant"][1], r["condition_value"]): r["model"] for r in rows
                 if r["property"] == "liquid_molar_density_mol_m3" and r["variant"][0] == "F"}
    for row in rows:
        if row["property"] == "liquid_molar_density_mol_m3" and row["condition_value"] > 0.0:
            row["observed"] = reference.get((row["variant"][1], row["condition_value"]))
            row["reference_kind"] = f"F{row['variant'][1]} model"
            if row["observed"] is not None and row["model"] is not None:
                row["error"] = row["model"] - row["observed"]
                row["relative_error"] = row["error"] / row["observed"]
    results = LEAF / "results"
    results.mkdir(exist_ok=True)
    with (results / "comparison.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(dict.fromkeys(k for row in rows for k in row)), lineterminator="\n")
        writer.writeheader()
        writer.writerows({k: ("" if v is None else v) for k, v in row.items()} for row in rows)
    with (results / "refit.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=[k for k in dict.fromkeys(k for fit in fits for k in fit)
                                                   if k not in ("failures", "holdout_check", "covariance_assumption",
                                                                "native_failed_fit_predictions", "native_failed_fit_residuals")], lineterminator="\n")
        writer.writeheader()
        writer.writerows({k: v for k, v in fit.items() if k in writer.fieldnames} for fit in fits)
    wheel = next(p for p in (metadata.distribution("epcsaft").files or []) if p.name.endswith(".so"))
    inputs = [HAMER_WU, LEDGER, NIST, *sorted((LEAF / "inputs/selections").glob("*.json")),
              ROOT / "data/packets/dielectric-response-records/1/parameters/dielectric-records.csv"]
    summary = {
        "question": "Do the admitted density-temperature and dipolar-perturbation epsilon_r forms, paired with the "
                    "Figiel 2025 water/Na+/Cl- base, describe NaCl(aq) MIAC, dielectric decrement, solvation and "
                    "density better than Figiel's ion-fraction suppression at 298.15 K and 1 bar?",
        "wheel_sha256": (ROOT / "build/environment-wheel/wheel.sha256").read_text().strip(),
        "method": {"owner": "Engine Ceres", "objective": "sum(relative MIAC residual squared) = 2 * final_cost",
                   "origin": FIGIEL_K, "scale": 1.0, "bounds": [-0.8, 0.4], "optimizer_starts": [[0.0]],
                   "holds": [], "holdout_weights": 0.0, "uncertainty": "conditional SSE/(21-1), exact tied columns"},
        "engine_version": metadata.version("epcsaft"),
        "extension_sha256": hashlib.sha256(Path(wheel.locate()).read_bytes()).hexdigest(),
        "inputs": {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
        "calibration": f"shared water-ion k_ij, Hamer-Wu NaCl m <= {CALIBRATION_MAX_MOLALITY} mol/kg, relative loss",
        "refits": fits,
        "metrics": summarize(rows),
        "sole_root_states": sorted({(r["variant"], r["property"]) for r in rows if r["branch"] == "sole-root"}),
        "failures": [{k: r[k] for k in ("variant", "property", "condition_value", "status")}
                     for r in rows if r["status"] != "ok"] + [f for f in fits if f["status"] != "ok"],
        "derivative_failures": [{k: r[k] for k in ("variant", "condition_value", "derivative_status", "verification_status")}
                                for r in rows if r.get("derivative_status", "").startswith("failed")
                                or r.get("verification_status", "").startswith("failed")],
        "limits": [
            "Hybrids: Figiel ion, Born and water parameters were fitted with ion-fraction suppression; "
            "only one shared water-ion k_ij is refitted, so D1/P1 are not source models.",
            "Hamer-Wu MIAC uncertainty is not retained; propagated bands cover k_ij only.",
            "Dielectric rows are Figiel Fig. 2 digitizations without pointwise uncertainty; static-permittivity "
            "data can include kinetic depolarization (Rueben 2024). They calibrated Figiel's C_epsilon.",
            "No retained NaCl(aq) solution-density reference: solution density is a model difference only.",
            "298.15 K and 1 bar only; no temperature, mixed-solvent or other-salt claim.",
        ],
    }
    (results / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    assert not summary["failures"] and not summary["derivative_failures"], "failed calculations retained in summary.json"


def declarations():
    return [
        (
            model,
            coordinate(model.parameters),
            [
                electrolyte(model, kind, observed=1)
                for kind in (
                    'mean_ionic_activity_coefficient',
                    'osmotic_coefficient',
                    'bulk_relative_permittivity',
                )
            ],
        )
        for model in [selection(form, FIGIEL_K) for form in FORMS]
    ]


if __name__ == "__main__":
    main()
