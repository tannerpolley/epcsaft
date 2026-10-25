"""Acetone-n-pentane VLE: PCP-SAFT (own parameters) vs PC-SAFT (own parameters) (#28).

Run after compare.py: .venv/bin/python analyses/2026-polar-dipole/scripts/binary.py
Data: inputs/gv2006-fig15-acetone-pentane-digitized.csv (Campbell et al. 1986 experimental symbols
digitized from GV2006 Fig. 15). Liquid points give bubble pressures, vapor points dew pressures.
Acetone: GV2006 Table 2 (polar) or the results/summary.json nonpolar refit; n-pentane: Gross-Sadowski 2001
in both. Each model is run with k_ij = 0 and with one k_ij calibrated at 372.7 K and validated at
397.7 and 422.6 K. Unconverged points are retained as failure rows and excluded from AAD.
The printed Campbell et al. 1986 Table I rows (inputs/campbell1986-table1.csv) are then evaluated as bubble
points at their measured T and x: mean |y_model - y_data| of n-pentane and bubble-P AAD, next to the paper's
own four-suffix Margules correlation (eq 8, fitted per isotherm to P-x only by Barker's method, Table III p. 427;
Table I "calculated" columns), whose y is a Gibbs-Duhem prediction and so the scatter floor of the data.
"""
from __future__ import annotations

import copy
import csv
import json
from pathlib import Path

from epcsaft import Mixture, Parameters, regression
from epcsaft.equilibrium import Amounts, Free, Phase, Pinned, Problem

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parents[1]
GV = ROOT / "data/packets/gross-vrabec-2006/1/parameters/parameters.json"
GS = ROOT / "data/packets/gross-sadowski-2001/3/parameters/parameters.json"
DATA = HERE / "inputs/gv2006-fig15-acetone-pentane-digitized.csv"
PRINTED = HERE / "inputs/campbell1986-table1.csv"
CALIBRATION_T = 372.7
NAMES = ("acetone", "n-pentane")


def document(polar: bool, kij: float) -> dict:
    gv = json.loads(GV.read_text())
    acetone = copy.deepcopy(next(c for c in gv["components"] if c["component_id"] == "acetone"))
    pentane = copy.deepcopy(next(c for c in json.loads(GS.read_text())["components"] if c["component_id"] == "pentane"))
    pentane["component_id"] = "n-pentane"
    if polar:
        pentane["coefficients"] += [{"identity": f"n-pentane/{f}", "family": f, "value": {"magnitude": 0.0, "unit": u}}
                                    for f, u in (("dipole_moment", "debye"), ("dipole_count", "dimensionless"))]
    else:
        fitted = json.loads((HERE / "results/summary.json").read_text())["fits"]["acetone/nonpolar/refit"]["fitted"]
        acetone["coefficients"] = [c for c in acetone["coefficients"] if c["family"] not in ("dipole_moment", "dipole_count")]
        for c in acetone["coefficients"]:
            c["value"]["magnitude"] = fitted.get(c["family"], c["value"]["magnitude"])
    return {"schema": "epcsaft.parameters", "schema_version": 1, "components": [acetone, pentane],
            "pairs": [{"component_id_a": "acetone", "component_id_b": "n-pentane", "coefficients": [
                {"identity": "pair/acetone/n-pentane/k_ij", "family": "k_ij", "value": {"magnitude": kij, "unit": "dimensionless"}}]}],
            "model_families": [f for f in gv["model_families"] if polar or f["kind"] != "polar"]}


def point(row: dict) -> Problem:
    T, z, P = float(row["T_K"]), float(row["x_or_y_acetone"]), float(row["P_Pa"])
    liquid = row["phase"] == "liquid"
    return Problem(phases=[Phase("L", kind="liquid", amount=Free() if liquid else Pinned(0.0)),
                           Phase("V", kind="vapor", amount=Pinned(0.0) if liquid else Free())],
                   T=T, P=Free(P), feed=Amounts({"acetone": z, "n-pentane": 1.0 - z}), row_tolerance=1.0e-10)


def calibrate(polar: bool, rows: list[dict]) -> tuple[float | None, dict]:
    parameters = Parameters.from_mapping(document(polar, 0.0), components=NAMES)
    observations = [regression.observation(parameters, "pressure", point(r), observed=float(r["P_Pa"]), form="log_ratio", mass_basis=False)
                    for r in rows]
    coordinate = regression.coordinate(parameters, [("k_ij", ("acetone", "n-pentane"), 1)], origin=0.0, scale=0.05, bounds=(-0.2, 0.3))
    controls = regression.FitControls(maximum_iterations=100, maximum_elapsed_time_seconds=600.0)
    result = regression.fit(parameters, [coordinate], observations, controls=controls)
    return (result.physical[0] if result.usable else None), {"status": result.status, "final_cost": result.final_cost,
                                                             "trial_failures": result.trial_failures}


def paired(calibrated: dict[str, float]) -> tuple[list[dict], list[dict]]:
    """Bubble points at the printed Table I (T, x) rows: vapor n-pentane fraction and pressure."""
    data = [r for r in csv.DictReader(PRINTED.open()) if 0.0 < float(r["x_pentane"]) < 1.0]
    rows = []
    for polar in (True, False):
        label = "polar/GV2006" if polar else "nonpolar/refit"
        for variant, value in (("k_ij=0", 0.0), ("k_ij calibrated", calibrated[label])):
            parameters = Parameters.from_mapping(document(polar, value), components=NAMES)
            model = Mixture(parameters)
            for r in data:
                problem = point({"T_K": r["T_K"], "phase": "liquid", "x_or_y_acetone": 1.0 - float(r["x_pentane"]), "P_Pa": 1.0e3 * float(r["P_kPa"])})
                P = regression.predict(model, regression.observation(parameters, "pressure", problem, observed=1.0e3 * float(r["P_kPa"]),
                                                                     form="log_ratio", mass_basis=False))
                y = regression.predict(model, regression.observation(parameters, "mole_fraction", problem, observed=float(r["y_pentane"]),
                                                                     form="difference", scale=0.005, mass_basis=False, phase="V", species=["n-pentane"]))
                ok = all(q.status == regression.ObservationStatus.Available and q.value is not None for q in (P, y))
                rows.append({"model": label, "variant": variant, "k_ij": value, "T_K": float(r["T_K"]), "x_pentane": float(r["x_pentane"]),
                             "y_data": float(r["y_pentane"]), "y_model": y.value if ok else None, "P_data_kPa": float(r["P_kPa"]),
                             "P_model_kPa": P.value / 1.0e3 if ok else None, "status": "Available" if ok else f"{P.status.name}/{y.status.name}",
                             "role": "calibration" if abs(float(r["T_K"]) - CALIBRATION_T) < 0.5 else "validation"})
    rows += [{"model": "Margules/Campbell1986", "variant": "four-suffix eq 8 (Barker P-x fit) per isotherm", "k_ij": None, "T_K": float(r["T_K"]),
              "x_pentane": float(r["x_pentane"]), "y_data": float(r["y_pentane"]), "y_model": float(r["y_calc_pentane"]),
              "P_data_kPa": float(r["P_kPa"]), "P_model_kPa": float(r["P_calc_kPa"]), "status": "Available",
              "role": "calibration" if abs(float(r["T_K"]) - CALIBRATION_T) < 0.5 else "validation"} for r in data]
    summary = []
    for key in dict.fromkeys((r["model"], r["variant"], r["k_ij"]) for r in rows):
        for role in ("calibration", "validation"):
            group = [r for r in rows if (r["model"], r["variant"], r["k_ij"]) == key and r["role"] == role]
            done = [r for r in group if r["y_model"] is not None]
            summary.append({"model": key[0], "variant": key[1], "k_ij": key[2], "role": role, "points": len(done), "failures": len(group) - len(done),
                            "mean_abs_dy": sum(abs(r["y_model"] - r["y_data"]) for r in done) / len(done),
                            "max_abs_dy": max(abs(r["y_model"] - r["y_data"]) for r in done),
                            "AAD_P_percent": 100.0 * sum(abs(r["P_model_kPa"] / r["P_data_kPa"] - 1.0) for r in done) / len(done)})
    return rows, summary


def main() -> None:
    data = list(csv.DictReader(DATA.open()))
    calibration = [r for r in data if abs(float(r["T_K"]) - CALIBRATION_T) < 0.5]
    rows, summary, fits, calibrated = [], [], {}, {}
    for polar in (True, False):
        label = "polar/GV2006" if polar else "nonpolar/refit"
        kij, fits[label] = calibrate(polar, calibration)
        calibrated[label] = kij
        for variant, value in (("k_ij=0", 0.0), ("k_ij calibrated", kij)):
            if value is None:
                continue
            parameters = Parameters.from_mapping(document(polar, value), components=NAMES)
            model = Mixture(parameters)
            model_rows = []
            for r in data:
                prediction = regression.predict(model, regression.observation(parameters, "pressure", point(r), observed=float(r["P_Pa"]),
                                                                              form="log_ratio", mass_basis=False))
                ok = prediction.status == regression.ObservationStatus.Available and prediction.value is not None
                model_rows.append({"model": label, "variant": variant, "k_ij": value, "T_K": float(r["T_K"]), "phase": r["phase"],
                                   "x_or_y_acetone": float(r["x_or_y_acetone"]), "P_data_Pa": float(r["P_Pa"]),
                                   "P_model_Pa": prediction.value if ok else None, "status": prediction.status.name if ok else f"{prediction.status.name}: {prediction.message}",
                                   "role": "calibration" if abs(float(r["T_K"]) - CALIBRATION_T) < 0.5 else "validation",
                                   "P_rel_error": prediction.value / float(r["P_Pa"]) - 1.0 if ok else None})
            rows += model_rows
            for role in ("calibration", "validation"):
                for phase in ("liquid", "vapor"):
                    errors = [abs(m["P_rel_error"]) for m in model_rows if m["role"] == role and m["phase"] == phase and m["P_rel_error"] is not None]
                    summary.append({"model": label, "variant": variant, "k_ij": value, "role": role, "phase": phase,
                                    "points": len(errors), "failures": sum(1 for m in model_rows if m["role"] == role and m["phase"] == phase and m["P_rel_error"] is None),
                                    "AAD_P_percent": 100.0 * sum(errors) / len(errors) if errors else None})
    xy_rows, xy_summary = paired(calibrated)
    for name, table in (("binary-rows.csv", rows), ("binary-summary.csv", summary), ("binary-xy-rows.csv", xy_rows),
                        ("binary-xy-summary.csv", xy_summary)):
        with (HERE / "results" / name).open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(table[0]), lineterminator="\n")
            writer.writeheader()
            writer.writerows(table)
    out = json.loads((HERE / "results/summary.json").read_text())
    out["binary_kij_fits"] = fits
    out["binary_summary"] = summary
    out["binary_xy_summary"] = xy_summary
    (HERE / "results/summary.json").write_text(json.dumps(out, indent=1) + "\n")
    for item in summary + xy_summary:
        print(item)


def declarations():
    return [
        (
            Parameters.from_mapping(document(polar, 0)),
            [point(r) for r in list(csv.DictReader(DATA.open(encoding='utf-8-sig')))],
        )
        for polar in (False, True)
    ]


if __name__ == "__main__":
    main()
