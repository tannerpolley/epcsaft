"""Declare Baygi Table 2 pure MEA states; retain native sigma-cubed solutions and figure comparisons."""

import csv
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import time

import epcsaft
from epcsaft.equilibrium import solve_equilibrium

HERE = Path(__file__).resolve().parents[1]
ROOT = HERE.parents[2]
SOURCE = HERE / "inputs/source"
MW = 0.061083  # kg/mol; inherited comparison mass conversion, not Table 2.
TOLERANCES = {"row": 1e-10, "density_kg_m3": 1.5, "figure2_pressure_kPa": 1.1,
              "figure3_temperature_K": 1.5, "table2_percentage_point": 0.005}
SITES = {"2B": (1, 1), "3B": (2, 1), "4C": (2, 2)}


def read_csv(name):
    with (SOURCE / name).open() as stream:
        return list(csv.DictReader(stream))


def parameter_mapping(scheme, table):
    values = {r["quantity"]: r for r in table if r["source_row"] == f"MEA {scheme}"}
    def entry(quantity):
        r = values[quantity]
        return {"identity": f"MEA/{scheme}/{quantity}",
                "value": {"magnitude": float(r["value"]),
                          "unit": {"1": "dimensionless", "K": "kelvin"}.get(r["unit"], r["unit"])}}
    component = {"component_id": "monoethanolamine", "fixed": {
        "charge_number": {"value": {"magnitude": 0, "unit": "elementary-charge"}},
        "molar_mass": {"value": {"magnitude": MW, "unit": "kilogram / mole"}}}, "coefficients": [
        {"family": q, **entry(q)} for q in ("segment_count", "segment_diameter", "dispersion_energy_over_k")]}
    return {"schema": "epcsaft.parameters", "schema_version": 1, "components": [component],
            "model_families": [{"kind": k, "choice": v} for k, v in (
                ("base", "pc-saft"), ("association", "general-site"), ("electrolyte", "none"),
                ("permittivity", "none"))],
            "topology": {"sites": [{"component_id": "monoethanolamine", "site_id": s, "multiplicity": n}
                                    for s, n in zip(("donor", "acceptor"), SITES[scheme])],
                         "edges": [{"endpoint_a": {"component_id": "monoethanolamine", "site_id": "donor"},
                                    "endpoint_b": {"component_id": "monoethanolamine", "site_id": "acceptor"},
                                    "source": {"kind": "explicit"}, "volume": entry("association_volume"),
                                    "energy_over_k": entry("association_energy_over_k")}]}}


def paths(rows, series):
    grouped = {}
    for r in rows:
        if r["series"] == series:
            grouped.setdefault(r["path_index"], []).append((float(r["x"]), float(r["y"])))
    return grouped.values()


def intersects(point, curves, tolerances):
    """Intersect the segment parameter intervals allowed by BOTH coordinate tolerances."""
    for curve in curves:
        for start, end in zip(curve, curve[1:]):
            low, high = 0.0, 1.0
            for p, a, b, tol in zip(point, start, end, tolerances):
                delta = b - a
                if delta == 0:
                    if abs(p - a) > tol:
                        high = -1.0
                else:
                    bounds = sorted(((p - tol - a) / delta, (p + tol - a) / delta))
                    low, high = max(low, bounds[0]), min(high, bounds[1])
            if low <= high:
                return True
    return False


def source_at(curves, coordinate, axis):
    """Linear drawing-coordinate readout, without endpoint extrapolation."""
    for curve in curves:
        for a, b in zip(curve, curve[1:]):
            if min(a[axis], b[axis]) <= coordinate <= max(a[axis], b[axis]) and a[axis] != b[axis]:
                fraction = (coordinate - a[axis]) / (b[axis] - a[axis])
                return a[1 - axis] + fraction * (b[1 - axis] - a[1 - axis])
    return None


def declarations():
    mappings = {scheme: parameter_mapping(scheme, read_csv("baygi-table-2.csv")) for scheme in SITES}
    models = {scheme: epcsaft.Pure(epcsaft.Parameters.from_mapping(mapping, components=["monoethanolamine"]))
              for scheme, mapping in mappings.items()}
    return models, mappings


def main():
    if "EPCSAFT_HSL_LIBRARY" in os.environ:
        raise RuntimeError("Leave EPCSAFT_HSL_LIBRARY unset")
    assert Path(sys.prefix).resolve() == ROOT / ".venv", sys.prefix
    # A rectangle's corner can be outside a diagonal segment's allowed band.
    assert not intersects((0, 1), [[(0, 0), (1, 1)]], (0.1, 0.1))
    assert intersects((0.5, 0.6), [[(0, 0), (1, 1)]], (0.1, 0.1))
    table = read_csv("baygi-table-2.csv")
    figures = {n: read_csv(f"baygi-figure-0{n}.csv") for n in (2, 3)}
    models, mappings = declarations()
    records, diagnostics = [], []
    started = time.perf_counter()
    for scheme, model in models.items():
        for i in range(15):
            temperature = 303.15 + 10 * i
            problem = model.saturation(temperature)
            problem.row_tolerance = TOLERANCES["row"]
            result = solve_equilibrium(model, problem)
            densities = result.molar_densities
            valid = (result.success and result.requested_tolerance_met and result.pressure > 0
                     and len(densities) == 2 and densities[0] > densities[1] > 0)
            density = densities[0] * MW if len(densities) == 2 else math.nan
            row = {"scheme": scheme, "temperature_K": temperature, "pressure_Pa": result.pressure,
                   "rhoL_mol_m3": densities[0] if densities else math.nan,
                   "rhoV_mol_m3": densities[1] if len(densities) == 2 else math.nan,
                   "rhoL_kg_m3": density, "association_convention": "sigma^3",
                   "native_verified": valid, "requested_tolerance_met": result.requested_tolerance_met,
                   "native_success": result.success, "residual": result.convergence.residual,
                   "floor": result.convergence.floor, "stop": str(result.convergence.stop),
                   "iterations": result.iterations, "evaluations": result.evaluations}
            for figure in (2, 3):
                y = result.pressure / 1000 if figure == 2 else temperature
                tol = TOLERANCES["figure2_pressure_kPa" if figure == 2 else "figure3_temperature_K"]
                curves = list(paths(figures[figure], f"MEA-{scheme}"))
                row[f"figure{figure}"] = ("reproduced" if intersects(
                    (density, y), curves, (TOLERANCES["density_kg_m3"], tol))
                                         else "not reproduced") if valid else "undetermined"
                row[f"figure{figure}_cause"] = (
                    "" if row[f"figure{figure}"] == "reproduced" else "undetermined")
                source_value = source_at(curves, density if figure == 2 else temperature, 0 if figure == 2 else 1)
                row[f"figure{figure}_source_readout"] = source_value
                row[f"figure{figure}_signed_deviation"] = (
                    (y if figure == 2 else density) - source_value if source_value is not None and valid else None)
                dippr = list(paths(figures[figure], "DIPPR"))
                row[f"dippr_figure{figure}"] = ("reproduced" if intersects(
                    (density, y), dippr, (TOLERANCES["density_kg_m3"], tol))
                    else "not reproduced") if valid else "undetermined"
                reference = source_at(dippr, density if figure == 2 else temperature, 0 if figure == 2 else 1)
                row[f"dippr_figure{figure}_reference"] = reference
                row[f"dippr_figure{figure}_signed_deviation"] = (
                    (y if figure == 2 else density) - reference if reference is not None and valid else None)
            records.append(row)
            convergence = result.convergence
            diagnostics.append({"scheme": scheme, "temperature_K": temperature, "message": result.message,
                "solver_status": result.solver_status, "amounts_mol": result.amounts,
                "mole_fractions": result.mole_fractions, "residuals": result.residuals,
                "classified_residuals": result.classified_residuals,
                "rows": [{"eqid": r.eqid, "phase": r.phase, "index": r.index} for r in result.rows],
                "convergence": {k: getattr(convergence, k) for k in (
                    "residual", "floor", "amplification", "contraction", "iterations", "evaluations")},
                "criterion": None if convergence.criterion is None else convergence.criterion.value,
                "stop": str(convergence.stop), "linear_solver": result.linear_solver})
            print(scheme, temperature, result.pressure, density, valid, row["stop"], flush=True)
    elapsed = time.perf_counter() - started
    output = HERE / "results/current-pure-saturation"
    with output.with_suffix(".csv").open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)
    with (HERE / "results/baygi-pure-saturation.csv").open() as stream:
        historical = {(r["scheme"], float(r["temperature_k"])): r for r in csv.DictReader(stream)}
    history_checks = []
    for scheme in SITES:
        selected = [r for r in records if r["scheme"] == scheme]
        changes = {key: max(abs(r[key] / float(historical[(scheme, r["temperature_K"])][old]) - 1)
                           for r in selected) for key, old in (
                               ("pressure_Pa", "saturation_pressure_pa"), ("rhoL_mol_m3", "liquid_density_mol_m3"))}
        history_checks.append({"scheme": scheme, "maximum_relative_change_from_historical_sigma3": changes})
    metadata = {"class": "PUBLISHED_TUPLE_PROPERTY_REPLAY", "association_convention": "sigma^3",
        "engine_provenance": epcsaft.provenance(), "engine_module": epcsaft.__file__,
        "interpreter": sys.executable, "base_main": subprocess.check_output(
            ["git", "rev-parse", "origin/main"], cwd=ROOT, text=True).strip(),
        "input_sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in SOURCE.iterdir() if p.is_file()},
        "producer_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "molar_mass_kg_mol": MW, "tolerances": TOLERANCES, "parameter_mappings": mappings,
        "figure_readout_units": {"figure2": "kPa at native liquid density", "figure3": "kg/m3 at native T"},
        "readout_unavailable": "native coordinate outside source path; no extrapolation",
        "evidence": "numerical verification against paper model and empirical correlation curves; no raw measurements",
        "exact_d3_falsifier": {"scheme": "2B", "molar_density_mol_m3": 15000,
                               "temperature_K": [303.15, 373.15, 443.15], "ratio_relative_tolerance": 1e-12},
        "wall_seconds": elapsed, "states": len(records), "verified_states": sum(r["native_verified"] for r in records),
        "historical_sigma3_comparison": history_checks,
        "table2": [{"scheme": r["association_scheme"], "quantity": r["quantity"],
                    "reported_percent": float(r["value"]), "outcome": "undetermined",
                    "cause": "author sampling grid unknown; density Eq. 10 grouping unresolved"}
                   for r in table if r["source_row"].startswith("MEA ") and r["quantity"].endswith("AAD")],
        "diagnostics": diagnostics}
    output.with_suffix(".json").write_text(json.dumps(metadata, indent=2, allow_nan=False) + "\n")
    print(f"wall_seconds={elapsed:.6f}", flush=True)


if __name__ == "__main__":
    main()
