"""Bounded installed-wheel comparison; imports and model setup are reported separately."""
import argparse
import csv
import hashlib
import json
import os
import resource
import statistics
import time
from pathlib import Path

IMPORT_START = time.perf_counter_ns()
from epcsaft import Mixture, Parameters
from epcsaft.equilibrium import Amounts, Phase, Problem, solve_equilibrium

IMPORT_MS = (time.perf_counter_ns() - IMPORT_START) / 1e6

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parents[1] / "results"
PACKET = ROOT / "data/packets/gross-2001-methane-ethane/3/parameters/parameters.json"
COMPONENTS = ("methane", "ethane")


def declare_inputs():
    return Problem(phases=[Phase("L", kind="liquid"), Phase("V", kind="vapor")],
                          T=180., P=1e6, feed=Amounts({"methane": .5, "ethane": .5}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("label")
    parser.add_argument("--samples", type=int, default=15)
    parser.add_argument("--baseline", type=Path, default=OUT / "performance-before.csv")
    args = parser.parse_args()
    setup_start = time.perf_counter_ns()
    model = Mixture(Parameters.from_json(PACKET, components=COMPONENTS))
    setup_ms = (time.perf_counter_ns() - setup_start) / 1e6
    first_start = time.perf_counter_ns()
    first = solve_equilibrium(model, declare_inputs())
    first_solve_ms = (time.perf_counter_ns() - first_start) / 1e6
    cases = [(180., 1e6, .5, 1., False), (179., 1e6, .5, 1., False),
             (181., 1e6, .5, 1., False), (180., .8e6, .5, 1., False),
             (180., 1.2e6, .5, 1., False), (180., 1e6, .4, 1., False),
             (180., 1e6, .6, 1., False), (180., 1e6, .5, 10., False),
             (180., 1e6, .5, 1., True)]
    rows = []
    root_errors = []
    reference = None
    for case, (temperature, pressure, z, scale, reverse) in enumerate(cases):
        cold_result = None
        for mode in ("cold", "previous"):
            walls, cpus = [], []
            for sample in range(-2, args.samples):
                start, cpu = time.perf_counter_ns(), time.process_time_ns()
                phases = [Phase("L", kind="liquid"), Phase("V", kind="vapor")]
                if reverse:
                    phases.reverse()
                problem = Problem(phases=phases, T=temperature, P=pressure,
                                  feed=Amounts({"methane": scale*z, "ethane": scale*(1-z)}))
                # Same-state for case 0; genuinely changed controls for cases 1..7.
                # Reversed phase order uses its own result because seeds are positional.
                seed = (reference if not reverse else cold_result) if mode == "previous" else None
                result = solve_equilibrium(model, problem, seed)
                wall_ms = (time.perf_counter_ns() - start) / 1e6
                cpu_ms = (time.process_time_ns() - cpu) / 1e6
                if sample < 0:
                    continue
                liq, vap = (1, 0) if reverse else (0, 1)
                rows.append({"case": case, "mode": mode, "sample": sample, "T_K": temperature,
                    "P_Pa": pressure, "z_methane": z, "feed_mol": scale, "reverse": reverse,
                    "wall_ms": wall_ms, "cpu_ms": cpu_ms, "success": result.success,
                    "iterations": result.iterations, "evaluations": result.evaluations,
                    "hessian_calls": result.hessian_calls, "residual": result.convergence.residual,
                    "beta": result.amounts[vap]/scale, "x_methane": result.mole_fractions[liq*2],
                    "y_methane": result.mole_fractions[vap*2], "rho_L": result.molar_densities[liq],
                    "rho_V": result.molar_densities[vap]})
                walls.append(wall_ms)
                cpus.append(cpu_ms)
            if mode == "cold":
                cold_result = result
                if case == 0:
                    reference = result
            print(case, mode, result.success, result.iterations,
                  "wall_ms", statistics.median(walls), "cpu_ms", statistics.median(cpus), flush=True)
            for phase, kind in ((liq, "liquid"), (vap, "vapor")):
                root = model.state(temperature, P=pressure,
                    x=result.mole_fractions[2*phase:2*phase+2], phase=kind)
                root_errors.append(abs(root.molar_density / result.molar_densities[phase] - 1))
    with (OUT / f"performance-{args.label}.csv").open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    wheel = next((ROOT / "build/environment-wheel").glob("epcsaft-*.whl"))
    metadata = {"label": args.label, "warmups": 2, "samples": args.samples, "import_ms": IMPORT_MS,
        "model_setup_ms": setup_ms, "first_solve_ms": first_solve_ms, "first_solve_success": first.success,
        "peak_rss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        "wheel_sha256": hashlib.sha256(wheel.read_bytes()).hexdigest(),
        "parameter_sha256": hashlib.sha256(PACKET.read_bytes()).hexdigest(),
        "affinity": sorted(os.sched_getaffinity(0)), "nice": os.getpriority(os.PRIO_PROCESS, 0),
        "maximum_independent_density_root_relative_error": max(root_errors),
        "threads": {k: os.environ.get(k) for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")},
        "limits": "One binary interior-flash neighborhood; warm process, cold nonlinear starts. No external reruns."}
    (OUT / f"performance-{args.label}.json").write_text(json.dumps(metadata, indent=2) + "\n")
    assert all(row["success"] for row in rows), "Inspect retained failed cases"
    assert max(root_errors) < 1e-7, "Returned phase does not match its requested density branch"
    with (OUT / "values.csv").open() as handle:
        oracles = list(csv.DictReader(handle))
    for oracle in oracles:
        for key, field, tolerance in (("beta", "beta_vapor", 1e-8),
                ("x_methane", "x_methane", 1e-8), ("y_methane", "y_methane", 1e-8),
                ("rho_L", "rho_liquid_mol_m3", 1e-7), ("rho_V", "rho_vapor_mol_m3", 1e-7)):
            if oracle[field]:
                expected = float(oracle[field])
                assert abs(rows[0][key] - expected) <= tolerance * max(1., abs(expected))
    if args.label != "before":
        with args.baseline.open() as handle:
            previous = {(int(row["case"]), row["mode"]): row for row in csv.DictReader(handle)}
        for row in rows:
            old = previous[row["case"], row["mode"]]
            for key in ("beta", "x_methane", "y_methane", "rho_L", "rho_V"):
                tolerance = 1e-7 if key.startswith("rho") else 1e-8
                assert abs(row[key] - float(old[key])) <= tolerance * max(1., abs(float(old[key])))


def declarations():
    return (Mixture(Parameters.from_json(PACKET, components=COMPONENTS)), declare_inputs())


if __name__ == "__main__":
    main()
