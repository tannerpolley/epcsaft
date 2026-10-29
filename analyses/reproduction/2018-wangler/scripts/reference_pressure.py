"""Replay the superseded one-bar standards to isolate Eqs. 20–21's correction."""

import csv
import json
import runpy
from pathlib import Path
import epcsaft
from epcsaft import equilibrium

LEAF = Path(__file__).resolve().parents[1]
DECLARATION = runpy.run_path(str(LEAF / "figures/figure-14/scripts/generate_data.py"))


def main():
    rows = []
    fitted = json.loads(
        (
            LEAF / "figures/figure-10/output/wangler-2018-mdeah-osmotic-summary.json"
        ).read_text()
    )["lanes"]
    for lane, dispersion in DECLARATION["DISPERSION"].items():
        model = epcsaft.Mixture(
            DECLARATION["build_parameters"](
                303.0, fitted[lane]["k_ij_303_k"], dispersion
            )
        )
        for source in csv.DictReader(DECLARATION["SOURCE"].open()):
            alpha = float(source["co2_loading_mol_per_mol_mdea"])
            observed = float(source["experimental_total_pressure_bar"])
            problem = DECLARATION["problem"](alpha, observed * 1e5)
            for reaction in problem.reactions:
                reaction.reference = equilibrium.ReactionReference(
                    equilibrium.ReferenceSourceBasis.RawMoleFractionInfiniteDilution,
                    reference_pressure_pa=1e5,
                )
            result = equilibrium.solve_equilibrium(model, problem)
            assert result.success, result.message
            rows.append(
                dict(
                    lane=lane,
                    alpha=alpha,
                    observed_bar=observed,
                    reference_pressure_pa=1e5,
                    pressure_bar=result.pressure / 1e5,
                    maximum_residual=max(map(abs, result.residuals)),
                )
            )
    with (LEAF / "results/reference-pressure-baseline.csv").open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    for lane in fitted:
        selected = [r for r in rows if r["lane"] == lane]
        print(
            lane,
            100
            * sum(abs(r["pressure_bar"] / r["observed_bar"] - 1) for r in selected)
            / len(selected),
        )


if __name__ == "__main__":
    main()
