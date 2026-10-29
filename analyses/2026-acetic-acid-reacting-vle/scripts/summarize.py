"""Reduce retained native states to pressure errors, controls and interval checks."""
import csv
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def write_csv(name, rows):
    with (ROOT / "results" / name).open("w") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    states = json.loads((ROOT / "results/states.json").read_text())
    reference = json.loads((ROOT / "results/reference-refinement.json").read_text())
    config = json.loads((ROOT / "inputs/reactions.json").read_text())
    measured = [r for r in csv.DictReader((ROOT / "inputs/measurements.csv").open()) if r["selected"] == "True"]
    pressure, controls, intervals = [], [], []
    for m in measured:
        number, xw = int(m["row"]), float(m["x_water"])
        group = [r for r in states if r["measurement_row"] == number]
        grid = sorted((r for r in group if r["variant"] == "primary" and not r["tighter"]), key=lambda r: r["x_water"])
        assert len(grid) == config["interval_grid_points"]
        def central(variant):
            return next(r for r in group if r["variant"] == variant and not r["tighter"] and abs(r["x_water"] - xw) < 1e-14)
        p, no, acid, original = (central(v) for v in ("primary", "no-dissociation", "acid-only", "water-only-permittivity"))
        tight = next(r for r in group if r["tighter"])
        error = (p["P_Pa"] - float(m["P_Pa"])) / 1000
        pressure.append(dict(row=number, x_water=xw, P_measured_kPa=float(m["P_Pa"]) / 1000,
            U95_P_kPa=float(m["U95_P_Pa"]) / 1000, P_calculated_kPa=p["P_Pa"] / 1000,
            signed_error_kPa=error, signed_relative_error_percent=100 * error / (float(m["P_Pa"]) / 1000)))
        omission = max(abs(r["P_Pa"] - next(a["P_Pa"] for a in group if a["variant"] == "acid-only" and a["x_water"] == r["x_water"])) for r in grid)
        controls.append(dict(row=number, P_no_dissociation_kPa=no["P_Pa"] / 1000,
            P_original_water_only_kPa=original["P_Pa"] / 1000,
            dissociation_effect_Pa=p["P_Pa"] - no["P_Pa"],
            original_dissociation_effect_Pa=original["P_Pa"] - no["P_Pa"],
            autoprotolysis_effect_Pa=p["P_Pa"] - acid["P_Pa"], autoprotolysis_grid_max_abs_Pa=omission,
            acid_permittivity_effect_Pa=p["P_Pa"] - original["P_Pa"],
            tighter_solve_abs_effect_Pa=abs(tight["P_Pa"] - p["P_Pa"]),
            reference_refinement_max_abs_Pa=next(r["maximum_absolute_pressure_effect_Pa"] for r in reference["comparisons"] if r["row"] == number),
            native_numerical_pressure_estimate_Pa=max(r["numerical_estimate"]["pressure"]["absolute"] for r in grid)))
        slopes = []
        for r in grid:
            assert r["checks_met"]
            est = r["numerical_estimate"]["pressure"]
            assert est["available"] and est["resolved"] and est["absolute"] < config["numerical_pressure_criterion_Pa"]
            a = r["pressure_actions"][0]
            assert a["status"] == "SolvedStateActionStatus.Available"
            slopes.append(a["action"])
        # The retained samples and their exact slopes qualify a monotone smooth branch numerically;
        # they are not interval arithmetic or a proof excluding unsampled extrema.
        monotone = min(slopes) > 0 and all(b["P_Pa"] > a["P_Pa"] for a, b in zip(grid, grid[1:]))
        assert monotone
        allowance = config["numerical_pressure_criterion_Pa"]
        lower, upper = grid[0]["P_Pa"] - allowance, grid[-1]["P_Pa"] + allowance
        ml, mu = float(m["P_Pa"]) - float(m["U95_P_Pa"]), float(m["P_Pa"]) + float(m["U95_P_Pa"])
        intervals.append(dict(row=number, x_low=grid[0]["x_water"], x_high=grid[-1]["x_water"],
            P_calculated_low_kPa=lower / 1000, P_calculated_high_kPa=upper / 1000,
            P_measured_low_kPa=ml / 1000, P_measured_high_kPa=mu / 1000,
            numerical_allowance_Pa=allowance, min_dP_dx_Pa=min(slopes), max_dP_dx_Pa=max(slopes),
            sampled_monotone=monotone, disjoint=lower > mu or upper < ml,
            separation_kPa=max(lower - mu, ml - upper, 0) / 1000))
    write_csv("pressures.csv", pressure)
    write_csv("controls.csv", controls)
    write_csv("uncertainty.csv", intervals)
    failure_record = dict(attempted=len(states) + reference["attempted"],
        failures=[r for r in states if not r.get("checks_met", False)] + reference["failures"])
    (ROOT / "results/failures.json").write_text(json.dumps(failure_record, indent=2) + "\n")
    summary = dict(temperature_K=config["temperature_K"], measurements=len(pressure), attempted=failure_record["attempted"],
        numerical_failures=len(failure_record["failures"]),
        mean_signed_error_kPa=sum(r["signed_error_kPa"] for r in pressure) / len(pressure),
        RMSE_kPa=math.sqrt(sum(r["signed_error_kPa"] ** 2 for r in pressure) / len(pressure)),
        mean_absolute_relative_error_percent=sum(abs(r["signed_relative_error_percent"]) for r in pressure) / len(pressure),
        maximum_absolute_relative_error_percent=max(abs(r["signed_relative_error_percent"]) for r in pressure),
        interval_disjoint_count=sum(r["disjoint"] for r in intervals),
        max_row_residual=max(r["row_residual"] for r in states),
        max_charge_error_mol=max(abs(r["charge_error_mol"]) for r in states),
        max_inventory_error_mol=max(abs(r[k]) for r in states for k in ("water_inventory_error_mol", "acid_inventory_error_mol", "hydrogen_inventory_error_mol")),
        max_relative_pressure_closure=max(r["pressure_closure_relative"] for r in states),
        max_neutral_log_fugacity_difference=max(abs(x) for r in states for x in r["neutral_log_fugacity_difference"]),
        max_reaction_log_residual=max(abs(x) for r in states for x in r["reaction_log_residuals"]),
        calibration_label=config["calibration_label"],
        uncertainty_label="Composition intervals at fixed temperature and fixed parameters; no joint confidence region; sampled monotone branch, not certified interval arithmetic.",
        reference_error_label="Analytic infinite-dilution insertion on accepted main. Native local higher-precision reference replay and pressure/site propagation included in pressure estimate. Independent finite-ladder comparison and native pressure re-solves retained separately. Neither local error estimates nor finite refinement differences are rigorous global bounds.")
    (ROOT / "results/summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
