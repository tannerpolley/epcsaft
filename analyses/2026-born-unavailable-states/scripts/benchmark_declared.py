"""Best-of-three declared solves for #240; frozen #223 inputs, no search."""
import argparse
import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import time

from epcsaft import equilibrium as eq

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parents[1]
source = ROOT / "analyses/2026-tpd-phase-discovery/scripts/baseline.py"
spec = importlib.util.spec_from_file_location("baseline_223", source)
baseline = importlib.util.module_from_spec(spec)
spec.loader.exec_module(baseline)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--label", required=True)
    parser.add_argument("--wheel", type=Path, required=True)
    parser.add_argument("--runs", type=int, choices=(1, 3), default=3)
    args = parser.parse_args()
    out = HERE / "results/native-parameter/performance"
    out.mkdir(parents=True, exist_ok=True)
    cases = {"case-6/figure-2-point-7", "case-7/flash-5918.5286-Pa"}
    rows, values, inputs = [], {}, {}
    for name, model, problem in baseline.declared_cases():
        if name not in cases:
            continue
        inputs[name] = baseline.problem_input(problem)
        for run in range(args.runs):
            start, cpu_start = time.perf_counter(), time.process_time()
            result = eq.solve_equilibrium(model, problem)
            seconds, cpu_seconds = time.perf_counter() - start, time.process_time() - cpu_start
            assert result.success and result.requested_tolerance_met, result.message
            rows.append((name, run, seconds, cpu_seconds, result.iterations, result.evaluations,
                         result.hessian_calls, max(map(abs, result.classified_residuals))))
            if name not in values:
                values[name] = {key: getattr(result, key) for key in (
                    "amounts", "mole_fractions", "molar_densities")} | {key: [getattr(result, key)] for key in ("pressure", "temperature")}
            print(args.label, name, run, seconds, result.iterations, flush=True)
    assert inputs.keys() == cases
    with (out / f"{args.label}-runs.csv").open("w") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(("case", "run", "wall_s", "cpu_s", "iterations", "evaluations", "hessian_calls", "classified_residual"))
        writer.writerows(rows)
    (out / f"{args.label}.json").write_text(json.dumps({
        "wheel": args.wheel.name, "wheel_sha256": hashlib.sha256(args.wheel.read_bytes()).hexdigest(),
        "baseline_script_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "inputs": inputs, "states": values,
        "packet_sha256": {str(p.relative_to(ROOT)): baseline.digest(p) for p in (baseline.KH_PACKET, baseline.MEA_PACKET)},
        "reaction_inputs": baseline.INPUTS["case-7/reactions"],
    }, indent=2) + "\n")


if __name__ == "__main__":
    main()
