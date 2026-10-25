"""MEA-water bubble points against Cai (1996), packet v3 against v4, and the v4 k_ij refit.

v3 resolves the MEA/water cross-association volume without the size factor; v4 restores
it (Fakouri Baygi 2015 Eq. 6). Each admitted Cai row fixes P and the liquid composition;
the Engine solves the bubble temperature with a pinned zero-amount vapor. The objective
and scales follow the retained neutral replay (1 K, 0.01 vapor water mole fraction).
"""

import copy
import csv
import hashlib
import importlib.metadata
import json
import math
import io
import subprocess
from pathlib import Path

import epcsaft
from epcsaft.equilibrium import Free, MoleFractions, Phase, Pinned, Problem, solve_equilibrium
from epcsaft import regression

REFERENCE_COMMIT = "f4075c0f"
ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parents[1]
CAI = ROOT / "analyses/2026-mea-parameter-estimation/inputs/source/library/observations/vapor_liquid_equilibrium/Cai_1996_MEA_H2O_VLE.csv"
PACKETS = {v: ROOT / f"data/packets/mea-co2-h2o-nine-species-estimation/{v}/parameters/parameters.json" for v in (3, 4)}
COMPONENTS = ("water", "monoethanolamine")
ADMITTED = {"binary_model_selection", "binary_training"}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def with_kij(document: dict, value: float) -> epcsaft.Mixture:
    mapping = copy.deepcopy(document)
    pair = next(p for p in mapping["pairs"] if {p["component_id_a"], p["component_id_b"]} == set(COMPONENTS))
    next(c for c in pair["coefficients"] if c["family"] == "k_ij")["value"]["magnitude"] = value
    return epcsaft.Mixture(epcsaft.Parameters.from_mapping(mapping, components=COMPONENTS))


def bubble_problem(row):
    water = float(row["water_liquid_mole_fraction"])
    return Problem(phases=[Phase("L", kind="liquid"), Phase("V", kind="vapor", amount=Pinned(0.0))],
                   T=Free(float(row["temperature_k"])), P=float(row["pressure_pa"]),
                   feed=MoleFractions({"water": water, "monoethanolamine": 1.0-water}))


def bubble(model: epcsaft.Mixture, row: dict, diagnostics=None) -> tuple[float, float]:
    result = solve_equilibrium(model, bubble_problem(row))
    if not result.success:
        raise RuntimeError(f"{row['row_id']}: {result.message}")
    if diagnostics is not None:
        diagnostics.append({"row": row["row_id"], "temperature_k": result.temperature,
                            "pressure_pa": result.pressure,
                            "liquid_water": float(row["water_liquid_mole_fraction"]),
                            "observed_temperature_k": float(row["temperature_k"]),
                            "observed_vapor_water": float(row["water_vapor_mole_fraction"]),
                            "vapor_water": result.mole_fractions[2],
                            "scaled_temperature_residual": result.temperature-float(row["temperature_k"]),
                            "scaled_vapor_water_residual": (result.mole_fractions[2]-float(row["water_vapor_mole_fraction"]))/0.01,
                            "liquid_density_mol_m3": result.molar_densities[0],
                            "vapor_density_mol_m3": result.molar_densities[1],
                            "maximum_row_residual": max(abs(x) for x in result.residuals),
                            "requested_tolerance_met": result.requested_tolerance_met,
                            "jacobian_rank": result.jacobian_rank, "unknown_count": result.unknown_count,
                            "topology_event": str(result.topology_event), "solver_status": result.solver_status})
    return result.temperature, result.mole_fractions[2]


def observations(parameters, rows):
    items = []
    for row in rows:
        problem = bubble_problem(row)
        items.extend((
            regression.observation(parameters, "temperature", problem, observed=float(row["temperature_k"]),
                                   form="difference", scale=1.0, mass_basis=False),
            regression.observation(parameters, "mole_fraction", problem,
                                   observed=float(row["water_vapor_mole_fraction"]), phase="V", species=("water",),
                                   form="difference", scale=0.01, mass_basis=False)))
    return items


def fit_rows(parameters, rows, origin, bounds, controls):
    coordinate = regression.coordinate(parameters, [("k_ij", COMPONENTS, 1)], origin=origin, scale=1.0, bounds=bounds)
    return regression.fit(parameters, [coordinate], observations(parameters, rows), weights=[1.0] * (2 * len(rows)), holds=(), controls=controls)


def fit_numbers(result):
    return {"k_ij": result.physical[0] if result.physical else None,
            "mean_scaled_square": 2 * result.final_cost / len(result.residuals) if result.residuals else None,
            "sse": 2 * result.final_cost, "status": result.status, "message": result.message,
            "usable": result.usable, "iterations": result.iterations,
            "rank": result.covariance.rank, "singular_values": result.covariance.singular_values,
            "active_bounds": result.covariance.active_bounds, "covariance_status": result.covariance.status,
            "variance_factor": result.covariance.variance_factor,
            "optimizer_covariance": result.covariance.optimizer,
            "projected_gradient_norm": result.covariance.projected_gradient_norm,
            "failures": [{"trial": f.trial, "observation": f.observation, "optimizer": f.optimizer,
                          "status": str(f.status), "message": f.message, "jacobian": f.jacobian} for f in result.failures]}


def objective(model: epcsaft.Mixture, rows: list[dict]) -> float:
    predictions = [regression.predict(model, item) for item in observations(model.parameters, rows)]
    if any(p.status != regression.ObservationStatus.Available for p in predictions):
        raise RuntimeError([(str(p.status), p.message) for p in predictions if p.status != regression.ObservationStatus.Available])
    return math.fsum(p.residual ** 2 for p in predictions) / len(predictions)


def historical(path):
    """Frozen pre-conversion retained values, never regenerated by a Python optimizer."""
    return subprocess.run(["git", "show", f"{REFERENCE_COMMIT}:{path.relative_to(ROOT)}"],
                          cwd=ROOT, check=True, capture_output=True, text=True).stdout


def comparison(case, row, quantity, stage, old, new, rule="absolute", tolerance=None):
    difference = new-old if old is not None and new is not None else None
    if tolerance is None:
        tolerance = 1e-6 if quantity == "temperature_k" else 1e-7
    if rule == "objective":
        tolerance = 1e-8*max(1.0, old)
        passed = difference <= tolerance if difference is not None else None
    elif rule == "secondary":
        tolerance, passed = None, None
    else:
        passed = abs(difference) <= tolerance if difference is not None else None
    return {"case": case, "row": row, "quantity": quantity, "stage": stage,
            "old": old, "new": new, "difference": difference, "rule": rule,
            "tolerance": tolerance, "passes": passed}


def write_comparison(output, table):
    with output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(table[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(table)


def compare_cai(output, rows, documents, fits, summary):
    old = json.loads(historical(output / "mea-water-cai.json"))
    old_rows = list(csv.DictReader(io.StringIO(historical(output / "mea-water-cai.csv"))))
    table, diagnostics = [], []
    for version, document in documents.items():
        case = f"v{version}"
        model = with_kij(document, old["packet_k_ij"])
        for row, prior in zip(rows, old_rows, strict=True):
            temperature, vapor = bubble(model, row, diagnostics)
            diagnostics[-1].update(case=case, stage="packet_coordinate", k_ij=old["packet_k_ij"])
            table.extend((comparison(case, row["row_id"], "temperature_k", "packet_coordinate", float(prior[f"temperature_v{version}_k"]), temperature),
                          comparison(case, row["row_id"], "vapor_water", "packet_coordinate", float(prior[f"vapor_water_v{version}"]), vapor)))
        table.append(comparison(case, "all", "mean_scaled_square", "packet_coordinate",
                                old["at_packet_k_ij"][case]["mean_scaled_square"], summary["at_packet_k_ij"][case]["mean_scaled_square"], "objective"))
        if version not in fits:
            continue
        prior = old["refit"][case] if "refit" in old else old["v5_refit"]
        result = fits[version]
        prior_model = with_kij(document, prior["k_ij"])
        old_coordinate_predictions = [regression.predict(prior_model, item) for item in observations(prior_model.parameters, rows)]
        if any(p.status != regression.ObservationStatus.Available for p in old_coordinate_predictions):
            raise RuntimeError("historical optimum replay has unavailable observations")
        for row in rows:
            bubble(prior_model, row, diagnostics)
            diagnostics[-1].update(case=case, stage="native_at_historical_optimum", k_ij=prior["k_ij"])
        table.extend((comparison(case, "all", "mean_scaled_square", "native_at_historical_optimum", prior["mean_scaled_square"],
                                 math.fsum(p.residual**2 for p in old_coordinate_predictions)/len(old_coordinate_predictions), "objective"),
                      comparison(case, "all", "mean_scaled_square", "native_optimum", prior["mean_scaled_square"], 2*result.final_cost/len(result.residuals), "objective"),
                      comparison(case, "all", "k_ij", "native_optimum", prior["k_ij"], result.physical[0], "secondary")))
        for i,row in enumerate(rows):
            for j,quantity in enumerate(("temperature_k", "vapor_water")):
                # Old optimum row values were not retained; identify this native replay explicitly.
                table.append(comparison(case, row["row_id"], quantity, "historical_coordinate_native_vs_native_optimum",
                                        old_coordinate_predictions[2*i+j].value, result.predictions[2*i+j]))
        fitted_model = epcsaft.Mixture(regression.fitted_parameters(model.parameters, result))
        for row in rows:
            bubble(fitted_model, row, diagnostics)
            diagnostics[-1].update(case=case, stage="native_optimum", k_ij=result.physical[0])
    with (output / "batch1-native-rows.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(diagnostics[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(diagnostics)
    write_comparison(output / "batch1-before-after.csv", table)
    summary["comparison"] = {"base_commit": REFERENCE_COMMIT, "old_wheel_sha256": old["identity"]["wheel_sha256"],
                             "old_refit": old.get("refit", old.get("v5_refit")),
                             "table": "batch1-before-after.csv", "failed_comparisons": sum(r["passes"] is False for r in table),
                             "old_optimum_rows": "native replay at the historical coordinate; original optimum rows were not retained"}


def main() -> None:
    rows = [r for r in csv.DictReader(CAI.open(encoding="utf-8-sig")) if r["role"] in ADMITTED]
    assert len(rows) == 25
    documents = {v: json.loads(p.read_text()) for v, p in PACKETS.items()}
    kij = next(c["value"]["magnitude"] for p in documents[3]["pairs"] if {p["component_id_a"], p["component_id_b"]} == set(COMPONENTS)
               for c in p["coefficients"] if c["family"] == "k_ij")
    table = []
    for row in rows:
        entry = {"row_id": row["row_id"], "pressure_pa": float(row["pressure_pa"]),
                 "liquid_water": float(row["water_liquid_mole_fraction"]), "temperature_cai_k": float(row["temperature_k"]),
                 "vapor_water_cai": float(row["water_vapor_mole_fraction"])}
        for version, document in documents.items():
            entry[f"temperature_v{version}_k"], entry[f"vapor_water_v{version}"] = bubble(with_kij(document, kij), row)
        table.append(entry)
    refits = {}
    fits = {}
    budget = regression.FitControls(maximum_iterations=100, maximum_elapsed_time_seconds=600.0)
    for version, document in documents.items():
        best = fit_rows(with_kij(document, kij).parameters, rows, kij, (kij - 0.03, kij + 0.03), budget)
        if not best.usable:
            raise RuntimeError(best.message)
        refits[f"v{version}"] = fit_numbers(best)
        fits[version] = best
    output = HERE / "results"
    with (output / "mea-water-cai.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(table[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(table)
    wheel = Path(json.loads(importlib.metadata.distribution("epcsaft").read_text("direct_url.json"))["url"].removeprefix("file://"))
    summary = {
        "packet_k_ij": kij,
        "at_packet_k_ij": {f"v{v}": {"mean_scaled_square": objective(with_kij(d, kij), rows),
                                     "max_abs_temperature_change_vs_v3_k": max(abs(e[f"temperature_v{v}_k"] - e["temperature_v3_k"]) for e in table),
                                     "max_abs_vapor_water_change_vs_v3": max(abs(e[f"vapor_water_v{v}"] - e["vapor_water_v3"]) for e in table)}
                           for v, d in documents.items()},
        "refit": refits,
        "identity": {"wheel_sha256": sha(wheel), "packets": {f"v{v}": sha(p) for v, p in PACKETS.items()}, "cai": sha(CAI)},
        "claim_limits": "25 admitted Cai 1996 rows, 66.66 and 101.33 kPa, 2B MEA; one-coordinate least squares with "
                        "the retained scales; no uncertainty interval or topology reselection.",
    }
    compare_cai(output, rows, documents, fits, summary)
    (output / "mea-water-cai.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps(summary, indent=1))


def declarations():
    from epcsaft import Parameters

    return [
        observations(
            Parameters.from_json(p, components=COMPONENTS),
            list(csv.DictReader(CAI.open(encoding='utf-8-sig'))),
        )
        for p in PACKETS.values()
    ]


if __name__ == "__main__":
    main()
