"""Fit the retained IAPWS-95 targets with and without total cp rows on the installed wheel."""
import csv
import json
import time
import sys
from pathlib import Path

from epcsaft import Parameters, Pure, regression
from epcsaft.eos import IdealCorrelation, IdealInterval, IdealPolynomial, ThermochemistryRecord
from epcsaft.equilibrium import Amounts, Free, Phase, Pinned, Problem

SET = Path(__file__).resolve().parents[1]
ROOT = SET.parents[1]
source = json.loads((ROOT / "data/packets/held-cameretti-sadowski-2008/3/parameters/parameters.json").read_text())
parameters = Parameters.from_mapping(source, components=["water"])
thermal = json.loads((SET / "inputs/ideal-gas-polynomial.json").read_text())
record = ThermochemistryRecord(thermal["standard_pressure_Pa"], [[IdealInterval(
    thermal["temperature_min_K"], thermal["temperature_max_K"], True, True, IdealCorrelation(
        IdealPolynomial(thermal["coefficients"], thermal["reference_temperature_K"],
                        thermal["enthalpy_at_reference_J_mol"], thermal["entropy_at_reference_J_mol_K"]),
        thermal["standard_pressure_Pa"]))]])
component = next(c for c in source["components"] if c["component_id"] == "water")
origins = {c["family"]: c["value"]["magnitude"] for c in component["coefficients"]}
sigma = next(c for c in source["correlations"] if c.get("component_id") == "water" and c["family"] == "segment_diameter")
origins["segment_diameter"] = sigma["constant"]["value"]["magnitude"]
edge = next(e for e in source["topology"]["edges"] if e["endpoint_a"]["component_id"] == "water")
origins["association_energy_over_k"] = edge["energy_over_k"]["value"]["magnitude"]
origins["association_volume"] = edge["volume"]["value"]["magnitude"]
chart = [("segment_count", 0.1, (0.8, 2.0)), ("segment_diameter", 0.05, (2.4, 3.2)),
         ("dispersion_energy_over_k", 10, (200, 500)),
         ("association_energy_over_k", 100, (1500, 3500)), ("association_volume", 0.005, (0.005, 0.2))]
coordinates = [regression.coordinate(parameters, [(family, (("water", "a"), ("water", "b")) if family.startswith("association") else "water", 1)],
    origin=origins[family], scale=scale, bounds=bounds) for family, scale, bounds in chart]
with (SET / "inputs/iapws95-water.csv").open() as stream:
    data = [{k: float(v) for k, v in row.items()} for row in csv.DictReader(stream)]
observations, quantities, temperatures = [], [], []
for row in data:
    T, P = row["T_K"], row["P_Pa"]
    problem = Problem(T=T, P=Free(row["P_sat_Pa"]), feed=Amounts({"water": 1}), phases=[
        Phase("L", kind="liquid", amount=Free(), composition_guess=(1.0,)),
        Phase("V", kind="vapor", amount=Pinned(0), composition_guess=(1.0,))], row_tolerance=1e-10)
    observations.extend([
        regression.density(parameters, T, P, {"water": 1}, phase="liquid", observed=row["rho_mol_m3"],
                           form="log_ratio", mass_basis=False),
        regression.isobaric_heat_capacity(parameters, T, P, {"water": 1}, phase="liquid",
                                         observed=row["cp_J_mol_K"], form="log_ratio"),
        regression.observation(parameters, "pressure", problem, observed=row["P_sat_Pa"],
                               form="log_ratio", mass_basis=False)])
    quantities.extend(["rho_mol_m3", "cp_J_mol_K", "P_sat_Pa"])
    temperatures.extend([T] * 3)


def declarations():
    return parameters, observations


if '--declare-only' in sys.argv:
    raise SystemExit(0)

def predictions(model):
    values = [regression.predict(model, observation) for observation in observations]
    assert all(value.status == regression.ObservationStatus.Available for value in values), values
    return [value.value for value in values]


rows, shifts = [], []
summary = {"base": "d47d0475", "wheel_sha256": json.loads(
    (ROOT / "build/environment-wheel/preparation.json").read_text())["wheel_sha256"],
    "source_class": "evaluated IAPWS-95 reference, numerical demonstration only", "fits": {}}
before = predictions(Pure(parameters, thermochemistry=record))
controls = regression.FitControls(maximum_iterations=100, maximum_elapsed_time_seconds=600.0)
for case, indices in [("with_cp", list(range(27))), ("without_cp", [i for i in range(27) if i % 3 != 1])]:
    start = time.perf_counter()
    result = regression.fit(parameters, coordinates, [observations[i] for i in indices], thermochemistry=record, controls=controls)
    wall = time.perf_counter() - start
    summary["fits"][case] = dict(status=result.status, termination=result.message,
        usable=result.usable, trial_failures=result.trial_failures, cost=result.final_cost,
        iterations=result.iterations, jacobian_evaluations=result.jacobian_evaluations,
        wall_seconds=wall, wall_seconds_per_jacobian=wall/result.jacobian_evaluations,
        covariance=result.covariance.status, rank=result.covariance.rank,
        singular_values=result.covariance.singular_values, active_bounds=result.covariance.active_bounds)
    assert result.status == "converged" and result.usable, summary["fits"][case]
    fitted = Pure(regression.fitted_parameters(parameters, result), thermochemistry=record)
    after = predictions(fitted)
    for i, (initial, final) in enumerate(zip(before, after)):
        reference = observations[i].observed
        rows.append(dict(case=case, quantity=quantities[i], T_K=temperatures[i],
            reference=reference, before=initial, after=final,
            before_deviation_percent=100*(initial/reference-1), after_deviation_percent=100*(final/reference-1)))
    for k, (family, scale, bounds) in enumerate(chart):
        shifts.append(dict(case=case, family=family, origin=origins[family], fitted=result.physical[k],
            shift_percent=100*(result.physical[k]/origins[family]-1), scale=scale,
            lower_bound=bounds[0], upper_bound=bounds[1]))
    print(case, summary["fits"][case], flush=True)
    # Timing at the converged model: direct physical Jacobian work on all rows of this fit.
    a_set = regression.CoordinateSet()
    a_set.coordinates = coordinates
    active_model = Pure(regression.fitted_parameters(parameters, result), thermochemistry=record,
                        active_parameters=a_set.active())
    start = time.perf_counter()
    gradients = [regression.differentiate(active_model, observations[i]) for i in indices]
    elapsed = time.perf_counter() - start
    assert all(g.status == regression.ObservationStatus.Available for g in gradients)
    summary["fits"][case]["final_full_jacobian_seconds"] = elapsed

metrics = []
for case in ("with_cp", "without_cp"):
    for quantity in quantities[:3]:
        selected = [r for r in rows if r["case"] == case and r["quantity"] == quantity]
        metrics.append(dict(case=case, quantity=quantity,
            before_AAD_percent=sum(abs(r["before_deviation_percent"]) for r in selected)/len(selected),
            before_maximum_percent=max(abs(r["before_deviation_percent"]) for r in selected),
            after_AAD_percent=sum(abs(r["after_deviation_percent"]) for r in selected)/len(selected),
            after_maximum_percent=max(abs(r["after_deviation_percent"]) for r in selected)))
for name, values in [("predictions.csv", rows), ("parameter-shifts.csv", shifts), ("deviations.csv", metrics)]:
    with (SET / "results" / name).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(values[0]))
        writer.writeheader()
        writer.writerows(values)
(SET / "results/fit-summary.json").write_text(json.dumps(summary, indent=2)+"\n")
print(metrics)
