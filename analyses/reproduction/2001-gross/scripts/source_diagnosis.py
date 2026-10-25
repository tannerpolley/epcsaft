"""Retain native pressure and temperature-uncertainty probes of source selections."""

import csv
import json
import time
from run import CALIBRATION, HERE, TOL, calculate, eq, model, write_rows


def main():
    started = time.perf_counter()
    rows = list(csv.DictReader((HERE / "results/figure-comparisons.csv").open()))
    case = json.loads((HERE / "inputs/digitized/ternary-pixels.json").read_text())["11"]
    data = [
        r
        for r in csv.DictReader((HERE / "inputs/digitized/ternary-values.csv").open())
        if r["figure"] == "11" and r["point"] == "2"
    ]
    x = [float(r["liquid_mole_fraction"]) for r in data]
    y = [float(r["vapor_mole_fraction"]) for r in data]
    m = model(case["components"], interactions=case["pairs"])
    out = []
    for p in (4800000.0, 3450000.0):
        problem = eq.Problem(
            T=case["temperature_k"],
            P=p,
            phases=[
                eq.Phase("L", kind="liquid", composition_guess=x),
                eq.Phase("V", kind="vapor", composition_guess=y),
            ],
            feed=eq.MoleFractions(
                dict(zip(case["components"], [(a + b) / 2 for a, b in zip(x, y)]))
            ),
            row_tolerance=TOL["equilibrium_maximum_absolute_row_residual"],
        )
        r = eq.solve_equilibrium(m, problem)
        assert r.success, r.message
        closure = max(
            abs(
                m.state(
                    case["temperature_k"],
                    rho=r.molar_densities[i],
                    x=r.mole_fractions[3 * i : 3 * i + 3],
                ).pressure
                / p
                - 1
            )
            for i in (0, 1)
        )
        assert (
            max(map(abs, r.residuals))
            <= TOL["equilibrium_maximum_absolute_row_residual"]
        )
        assert closure <= TOL["phase_pressure_relative_closure"]
        out.append(
            dict(
                figure="11",
                probe="printed pressure" if p == 4800000 else "only pressure changed",
                temperature_k=case["temperature_k"],
                pressure_pa=p,
                selected_liquid=json.dumps(x),
                selected_vapor=json.dumps(y),
                predicted_liquid=json.dumps(r.mole_fractions[:3]),
                predicted_vapor=json.dumps(r.mole_fractions[3:]),
                maximum_composition_absolute_error=max(
                    abs(a - b) for a, b in zip(x + y, r.mole_fractions)
                ),
                max_abs_residual=max(map(abs, r.residuals)),
                phase_pressure_relative_closure=closure,
                evidence="reproduced native probe",
            )
        )
    problem = eq.Problem(
        T=case["temperature_k"],
        P=eq.Free(3450000.0),
        phases=[
            eq.Phase("L", kind="liquid", composition_guess=x),
            eq.Phase("V", kind="vapor", amount=eq.Pinned(0), composition_guess=y),
        ],
        feed=eq.MoleFractions(dict(zip(case["components"], x))),
        row_tolerance=TOL["equilibrium_maximum_absolute_row_residual"],
    )
    r = eq.solve_equilibrium(m, problem)
    assert r.success, r.message
    out.append(
        dict(
            figure="11",
            probe="bubble pressure at selected liquid",
            temperature_k=case["temperature_k"],
            pressure_pa=r.pressure,
            selected_liquid=json.dumps(x),
            selected_vapor=json.dumps(y),
            max_abs_residual=max(map(abs, r.residuals)),
            evidence="reproduced native probe",
        )
    )
    out += [
        dict(
            figure="11",
            probe="alternative endpoint pairing",
            evidence="reviewer supplied diagnosis",
            conclusion="does not reduce the discrepancy",
        ),
        dict(
            figure="11",
            probe="interpret coordinates as mass fractions",
            evidence="reviewer supplied diagnosis",
            conclusion="does not reduce the discrepancy",
        ),
    ]
    write_rows("source-diagnosis.csv", out)
    review = list(
        csv.DictReader((HERE / "inputs/digitized/model-point-review.csv").open())
    )
    old = []
    x0, x1, xv0, xv1, y0, y1, yv0, yv1 = CALIBRATION["03"]["axis"]
    for r in review:
        if r["figure"] == "03":
            old.append(
                dict(
                    figure="03",
                    component=r["component"],
                    reference="model",
                    kind="density",
                    temperature_c="0",
                    x=xv0 + (float(r["old_pixel_x"]) - x0) * (xv1 - xv0) / (x1 - x0),
                    y=yv0 + (float(r["old_pixel_y"]) - y0) * (yv1 - yv0) / (y1 - y0),
                    y_uncertainty=TOL["digitization_pixel_uncertainty"]
                    * abs((yv1 - yv0) / (y1 - y0)),
                )
            )
    uncertainty = []
    for basis, source in [
        ("revised source point", rows),
        ("reviewer original source point", old),
    ]:
        for name in ["methane", "propylene"]:
            row = next(
                r
                for r in source
                if r["figure"] == "03"
                and r["component"] == name
                and r["reference"] == "model"
            )
            t = float(row["y"])
            dt = float(row["y_uncertainty"])
            values = []
            for temp in (t - dt, t + dt):
                r = calculate({**row, "y": temp})
                assert r["success"] and r["numerical_tolerance_met"]
                values.append(r["rho_liquid_kg_m3"])
            uncertainty.append(
                dict(
                    component=name,
                    basis=basis,
                    temperature_k=t,
                    temperature_uncertainty_k=dt,
                    pixel_uncertainty=TOL["digitization_pixel_uncertainty"],
                    source_density_kg_m3=row["x"],
                    minimum_liquid_density_kg_m3=min(values),
                    maximum_liquid_density_kg_m3=max(values),
                )
            )
    write_rows("density-temperature-uncertainty.csv", uncertainty)
    (HERE / "results/source-diagnosis-run.json").write_text(
        json.dumps({"wall_seconds": time.perf_counter() - started}, indent=2) + "\n"
    )


if __name__ == "__main__":
    main()
