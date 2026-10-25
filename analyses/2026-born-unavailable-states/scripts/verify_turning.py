"""Re-evaluate selected continued roots directly in the declared Born form."""

import argparse
from pathlib import Path
import csv
import json

import epcsaft
from epcsaft import equilibrium as eq
import numpy as np

from reproduce import HERE, archived, load_case, changed, problem_from_request

parser = argparse.ArgumentParser()
parser.add_argument(
    "--output-dir", type=Path, default=HERE / "results/current-main/turning"
)
args = parser.parse_args()
args.output_dir.mkdir(parents=True, exist_ok=True)
rows = []
corrections = []
for case, folder, indices in (
    (1, "turning", (0, 1, 23, 51, 64)),
    (2, "turning", (0, 1, 6, 7)),
    (2, "turning-refined", (9, 64)),
):
    points = list(csv.DictReader(archived()[f"{folder}/case{case}-coordinates.csv"].splitlines()))
    for index in indices:
        record = points[index]
        x = np.array([float(record[f"coordinate_{i}"]) for i in range(20)])
        q = float(record["Born_coefficient"])
        value = q if case == 1 else 1 / (1 - q)
        base, request, anchor = load_case(case)
        model = epcsaft.Mixture(
            epcsaft.Parameters.from_mapping(changed(base, case, value))
        )
        problem = problem_from_request(request, anchor)
        compiled = eq.compile_problem(model, problem)
        # #353 changed the reference limit. Archived roots supply guesses, not current roots.
        start = eq.EquilibriumStart()
        start.coordinates = list(x)
        result = eq.solve_equilibrium_starts(model, problem, [start], eq.EquilibriumWorkLimits()).results[0]
        assert result.success and result.requested_tolerance_met, result.message
        corrections.append(dict(case=case, path=folder, point=index, old_coordinates=list(x),
            old_rows=list(compiled.evaluate(x).values), coordinates=list(result.coordinates), iterations=result.iterations))
        x = np.array(result.coordinates)
        values = compiled.evaluate(x).values
        maximum = max(abs(v) for v in values)
        assert maximum <= 1e-10, (case, folder, index, maximum)
        for phase, (a, b) in enumerate(((0, 9), (11, 14))):
            fraction = np.exp(x[a:b])
            fraction /= fraction.sum()
            full = list(fraction) + ([0.0] * 6 if phase else [])
            state = model.state(
                request["temperature"]["value"],
                P=float(x[15]),
                x=full,
                phase="vapor" if phase else "liquid",
            )
            rho = state.molar_density
            rows.append(
                (
                    case,
                    folder,
                    index,
                    value,
                    maximum,
                    phase,
                    rho,
                    state.pressure_density_derivative,
                    state.born,
                    state.debye_huckel,
                    state.bulk_relative_permittivity,
                )
            )
            print(rows[-1], flush=True)
with (args.output_dir / "turning-verification.csv").open("w") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(
        (
            "case",
            "path",
            "point",
            "parameter",
            "row_max",
            "phase",
            "rho_mol_m3",
            "dP_drho",
            "born_over_RT",
            "debye_huckel_over_RT",
            "bulk_relative_permittivity",
        )
    )
    w.writerows(rows)
(args.output_dir / "turning-corrections.json").write_text(json.dumps({
    "engine_provenance": list(epcsaft.provenance()), "points": corrections}, indent=2)+"\n")
