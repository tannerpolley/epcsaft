"""Inspect the diagnosed dense start and the repaired cold start on one prepared Engine."""
import csv
import argparse
import json
import math
import os
import shlex
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
SET = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "engine"))
from tests.regression.test_fit import bubble, neutral_binary
from epcsaft import GAS_CONSTANT_J_PER_MOL_K as R, Mixture, trace, tolerances, _core
from epcsaft.equilibrium import Free, Reaction, ReactionLogPolynomial, solve_equilibrium

case = json.loads((ROOT / "analyses/2026-reaction-constant-fitting/results/cold-start.json").read_text())
i = case["inputs"]
model = Mixture(neutral_binary())
p = bubble(i["T_guess_K"])
p.T, p.P = Free(i["T_guess_K"]), i["P_fixed_Pa"]
p.reactions = [Reaction({"s0": -1.0, "s1": 1.0}, name=i["reaction"], correlation=ReactionLogPolynomial(
    i["a"], i["b_K"], i["c"], i["d_per_K"], i["Tref_K"], i["Tmin_K"], i["Tmax_K"], i["standard_state_id"]))]


def declarations():
    return model, p


if '--declare-only' in sys.argv:
    raise SystemExit(0)
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output-dir', type=Path, default=SET / 'results')
OUT = parser.parse_args().output_dir
OUT.mkdir(parents=True, exist_ok=True)
wheel, = (ROOT / 'build/environment-wheel').glob('*.whl')
subprocess.run([str(ROOT / ".venv/bin/python"), str(ROOT / "tools/verify-installed-wheel.py"),
                "--wheel", str(wheel)], check=True)
with trace() as captured:
    solved = solve_equilibrium(model, p)
assert solved.success, solved.message
assert solved.requested_tolerance_met
# Reuse the prepared native library and its exact compile/link flags for the diagnostic corrector.
build = Path((ROOT / "build/environment-wheel/meson-build-dir").read_text().strip())
scratch = ROOT / "build/issue204"
scratch.mkdir(parents=True, exist_ok=True)
commands = json.loads((build / "compile_commands.json").read_text())
row = next(r for r in commands if r["file"].endswith("/tests/native/boundary_checks.cpp"))
args = shlex.split(row["command"])
args[args.index("-o") + 1] = str(scratch / "density_probe.o")
args[-1] = str(SET / "scripts/density_probe.cpp")
args += ["-I" + str(ROOT / "engine/tests/native")]
subprocess.run(args, cwd=build, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
link = subprocess.check_output(["ninja", "-C", str(build), "-t", "commands", "boundary_checks"], text=True).splitlines()[-1]
args = shlex.split(link)
args[args.index("-o") + 1] = str(scratch / "density_probe")
args = [str(scratch / "density_probe.o") if a == "boundary_checks.p/tests_native_boundary_checks.cpp.o" else a for a in args]
subprocess.run(args, cwd=build, check=True)
probe = subprocess.run([str(scratch / "density_probe")] + [str(i[k]) for k in
    ("T_guess_K", "P_fixed_Pa", "a", "b_K", "c", "d_per_K", "Tref_K", "Tmin_K", "Tmax_K")],
    text=True, capture_output=True,
    env={**os.environ, 'LD_LIBRARY_PATH': str(build)+os.pathsep+os.environ.get('LD_LIBRARY_PATH', '')})
(OUT / 'density-probe.stdout').write_text(probe.stdout)
(OUT / 'density-probe.stderr').write_text(probe.stderr)
assert probe.returncode == 2 and 'phases_coincide:' in probe.stderr, (probe.returncode, probe.stderr)
assert 'FINAL,' not in probe.stdout
rows = []
for line in probe.stdout.splitlines():
    if line.startswith("SEED,"):
        kind, phase, *values = line.split(",")
        rows.append(["dense_seed", int(phase), *map(float, values)])
assert [row[:2] for row in rows] == [["dense_seed", 0], ["dense_seed", 1]]
for q in range(2):
    rows.append(["cold_fixed", q, solved.temperature, solved.pressure, solved.molar_densities[q],
                 solved.packing_fractions[q], *solved.mole_fractions[2*q:2*q+2]])
for row in rows:
    _, _, T, P, rho, eta, *x = row
    selected = model.state(T, P=P, x=x, phase="liquid" if row[1] == 0 else "vapor")
    anchor = selected.molar_density
    window = next(t.value for t in tolerances() if t.name == "DENSITY_ANCHOR_LOG_WINDOW")
    row.extend([model.state(T, rho=rho, x=x).pressure_density_derivative, anchor,
                anchor * math.exp(-window), anchor * math.exp(window), abs(math.log(rho/anchor))])
with (OUT / "densities.csv").open("w") as f:
    writer = csv.writer(f)
    writer.writerow(["case", "phase", "T_K", "P_Pa", "rho_mol_m3", "eta", "x0", "x1", "P_rho_Pa_m3_mol",
                     "anchor_mol_m3", "window_lower_mol_m3", "window_upper_mol_m3", "log_displacement"])
    writer.writerows(rows)

references = list(csv.DictReader((SET / "inputs/pressure-roots.csv").open()))
root_rows = []
pressure_tolerance = next(t.value for t in tolerances() if t.name == "DENSITY_ROOT_RESIDUAL")
for name, T, x in [("feed", i["T_guess_K"], [0.5, 0.5]),
                   ("dense_corrector_V", 201.10122071679316, [0.7819713975823668, 0.2180286024176333]),
                   ("J1_V", case["J1"]["T_K"][0], case["J1"]["x"][2:])]:
    upper = model.molar_density_upper_bound(T, x=x) * (1-1e-10)
    diagnostic = model.density_diagnostics(
        T, p.P, x=x, nodes=np.geomspace(1e-3, upper, 512), density_tolerance=1e-9,
        log_tolerance=1e-14, max_evaluations=100000, max_seconds=60.0, adaptive=True,
        finite_upper=0.0, evaluation_order=_core.DerivativeOrder.PressureSlope,
    )
    assert diagnostic.stop == "resolved_crossings", diagnostic.stop
    for kind, records in (("pressure_root", diagnostic.pressure_roots), ("spinodal", diagnostic.spinodals)):
        expected = sorted((r for r in references if r["case"] == name and r["kind"] == kind),
                          key=lambda r: float(r["rho_mol_m3"]))
        assert len(records) == len(expected)
        for record, reference in zip(records, expected):
            error = abs(record.molar_density - float(reference["rho_mol_m3"]))
            assert error <= 1e-9, (name, kind, error)
            eta = float(reference["eta"])
            eta_radius = abs(eta) * 1e-9 / float(reference["rho_mol_m3"])
            assert abs(record.packing_fraction - eta) <= eta_radius + math.ulp(eta) + math.ulp(record.packing_fraction)
            assert abs(record.pressure - float(reference["P_Pa"])) / (R * T * record.molar_density) <= pressure_tolerance
            assert record.upper - record.lower <= 1e-9 or record.stop == "exact_numerical_zero"
            if kind == "spinodal":
                assert abs(record.residual) <= 1e-10
            else:
                assert (record.pressure_density_derivative > 0) == (float(reference["P_rho_Pa_m3_mol"]) > 0)
            root_rows.append([name, kind, T, *x, record.molar_density, record.packing_fraction,
                              record.pressure, record.pressure_density_derivative,
                              record.lower, record.upper, record.residual, record.stop, reference["rho_mol_m3"], error])
with (OUT / "pressure-roots.csv").open("w") as f:
    writer = csv.writer(f)
    writer.writerow(["case", "kind", "T_K", "x0", "x1", "rho_mol_m3", "eta", "P_Pa", "P_rho_Pa_m3_mol", "bracket_lower_mol_m3", "bracket_upper_mol_m3", "residual_Pa_or_slope_over_RT", "stop", "reference_rho_mol_m3", "absolute_density_error_mol_m3"])
    writer.writerows(root_rows)
summary = {"input": "analyses/2026-reaction-constant-fitting/results/cold-start.json",
           "reference_commit": "bea7655d",
           "source_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
           "wheel_sha256": (ROOT / "build/environment-wheel/wheel.sha256").read_text().strip(),
           "cold_residual": solved.convergence.residual, "cold_message": solved.message,
           "equimolar_probe": {"returncode": probe.returncode, "message": probe.stderr.strip(),
                               "state_available": False},
           "seed_densities": [e.values[0] for e in captured.events if e.name == "branch_density"],
           "seed_source": "shared native initializer; dense diagnostic declares the former K=1 compositions",
           "evidence_kind": "numerical verification; no physical validation"}
(OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary))
