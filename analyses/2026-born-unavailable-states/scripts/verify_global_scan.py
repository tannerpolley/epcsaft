"""Historical carbon-omitted chart replay; see results/current-main/diagnosis.md.

The 026b3031/e9b1aed5 wheel comparison identifies the #385 reference change.
Current main rejects these old rows; no native held-composition corrector exists.
"""

import argparse
from pathlib import Path
import csv
import time
import json

import epcsaft
from epcsaft import equilibrium as eq
import numpy as np

from reproduce import (
    HERE,
    changed,
    load_case,
    problem_from_request,
    table,
)

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--output-dir", type=Path, default=HERE / "results/diagnostics/native/global-scan"
)
args = parser.parse_args()
args.output_dir.mkdir(parents=True, exist_ok=True)
started = time.process_time()
policy = {t.name: t.value for t in epcsaft.tolerances()}
checks, phase_checks = [], []
for case in (1, 2):
    base, request, anchor = load_case(case)
    value = 0 if case == 1 else 2
    mapping = changed(base, case, value)
    model = epcsaft.Mixture(epcsaft.Parameters.from_mapping(mapping))
    compiled = eq.compile_problem(model, problem_from_request(request, anchor))
    for folder in sorted((HERE / "results/global-scan").iterdir()):
        if not folder.is_dir() or not (
            folder.name == f"case{case}"
            or folder.name.startswith(f"case{case}-completion")
        ):
            continue
        attempts = {
            int(r["attempt"]): r
            for r in csv.DictReader((folder / "attempts.csv").open())
        }
        rows = {}
        for row in csv.DictReader((folder / "rows.csv").open()):
            rows.setdefault(int(row["attempt"]), []).append(float(row["residual"]))
        points = {
            int(r["attempt"]): np.array(
                [float(r[f"coordinate_{i}"]) for i in range(20)]
            )
            for r in csv.DictReader((folder / "coordinates.csv").open())
        }
        for index, old in rows.items():
            a = attempts[index]
            assert len(old) == 20
            assert max(abs(v) for v in old) == float(a["original_row_max"])
            assert old[16] == float(a["carbon_residual"])
            chart_max = max(abs(v) for i, v in enumerate(old) if i != 16)
            assert chart_max == float(a["chart_row_max"])
            if a["chart_criterion_met"] != "True":
                continue
            x = points[index]
            raw = np.array(compiled.evaluate(x).values)
            error = float(max(abs(raw - old)))
            assert error < 1e-10, (folder.name, index, error)
            maximum = float(max(abs(raw)))
            checks.append((case, folder.name, index, error, maximum, float(raw[16])))
            if maximum > 1e-10:
                continue
            fractions = [np.exp(x[:9]), np.r_[np.exp(x[11:14]), np.zeros(6)]]
            fractions = [f / sum(f) for f in fractions]
            distance = max(abs(x[10] - x[14]), max(abs(fractions[0] - fractions[1])))
            topology_ok = (
                x[9] / sum(request["reaction_system"]["feed_amounts_mol"])
                > policy["PHASE_DISAPPEARANCE_AMOUNT"]
                and distance > policy["PHASE_COINCIDENCE_DISTANCE"]
            )
            for phase in range(2):
                state = model.state(
                    request["temperature"]["value"],
                    P=float(x[15]),
                    x=fractions[phase],
                    phase="vapor" if phase else "liquid",
                )
                rho = state.molar_density
                branch_ok = state.pressure_density_derivative > 0
                phase_checks.append(
                    (
                        case,
                        folder.name,
                        index,
                        phase,
                        maximum,
                        rho,
                        state.pressure_density_derivative,
                        topology_ok,
                        branch_ok,
                        state.born,
                        state.debye_huckel,
                        state.bulk_relative_permittivity,
                    )
                )
table(
    args.output_dir,
    "row-verification.csv",
    (
        "case",
        "scan",
        "attempt",
        "reevaluation_error",
        "original_row_max",
        "carbon_residual",
    ),
    checks,
)
table(
    args.output_dir,
    "phase-verification.csv",
    (
        "case",
        "scan",
        "attempt",
        "phase",
        "original_row_max",
        "rho_mol_m3",
        "dP_drho_Pa_m3_per_mol",
        "topology_ok",
        "density_branch_ok",
        "born_over_RT",
        "debye_huckel_over_RT",
        "bulk_relative_permittivity",
    ),
    phase_checks,
)
(args.output_dir / "verification-summary.json").write_text(
    json.dumps(
        {
            "CPU_s": time.process_time() - started,
            "native_chart_points_rechecked": len(checks),
            "root_phase_checks": len(phase_checks),
            "criterion": 1e-10,
        },
        indent=2,
    )
    + "\n"
)
print(
    "verified chart points",
    len(checks),
    "root phase checks",
    len(phase_checks),
    flush=True,
)
