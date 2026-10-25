"""Owner-requested #238 cost control on two frozen default declared solves."""
import argparse
import csv
import importlib.metadata
import json
import os
from pathlib import Path
import time
import zipfile

import baseline as b
from epcsaft import equilibrium as eq

CASES = ("case-6/figure-2-point-7", "case-7/flash-5918.5286-Pa")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--label", required=True)
    parser.add_argument("--source-commit", required=True)
    parser.add_argument("--wheel", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--retain-failed", action="store_true")
    args = parser.parse_args()
    distribution = importlib.metadata.distribution("epcsaft")
    with zipfile.ZipFile(args.wheel) as wheel:
        members = [name for name in wheel.namelist() if name.startswith("epcsaft/") and not name.endswith("/")]
        assert all(distribution.locate_file(name).read_bytes() == wheel.read(name) for name in members), "installed wheel mismatch"
    args.output.mkdir(parents=True, exist_ok=True)
    rows = []
    for case, model, problem in b.declared_cases():
        if case not in CASES:
            continue
        b.INPUTS[case] = dict(components=model.component_ids, **b.problem_input(problem))
        for run in range(3):
            begun = time.perf_counter()
            result = eq.solve_equilibrium(model, problem)  # default ordinary single-start call
            seconds = time.perf_counter() - begun
            rows.append(dict(label=args.label, case=case, run=run, seconds=seconds,
                             iterations=result.iterations, evaluations=result.evaluations,
                             hessian_calls=result.hessian_calls, objective=result.objective,
                             raw_rows=max(map(abs, result.residuals)),
                             classified_rows=max(map(abs, result.classified_residuals)),
                             success=result.success, criterion_met=result.requested_tolerance_met,
                             phase_count=result.phase_count))
            b.state(f"{case}/run-{run}", result)
            b.record(f"{case}/run-{run}", "message", result.message)
            print(rows[-1], flush=True)
            if not args.retain_failed:
                assert result.success and result.requested_tolerance_met is True, result.message
    with (args.output / f"{args.label}-cost.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    b.csv_write(args.output / f"{args.label}-states.csv", b.ROWS)
    metadata = dict(source_commit=args.source_commit, wheel_sha256=b.digest(args.wheel),
                    inputs=b.INPUTS, input_files={str(path.relative_to(b.ROOT)): b.digest(path)
                        for path in (b.KH_PACKET, b.MEA_PACKET,
                            b.KH / "scripts/run_fixed_settler.py",
                            b.ROOT / "engine/tests/equilibrium/test_trace_reliability.py")},
                    thread_environment={key: os.environ.get(key) for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS")},
                    nice=os.getpriority(os.PRIO_PROCESS, 0),
                    evaluation_meaning="Engine evaluations counts Ipopt objective callbacks plus exact polish evaluations; no independent objective-count instrument added.",
                    limits="Best of three on a shared host; wall time does not isolate sub-percent overhead.")
    (args.output / f"{args.label}-provenance.json").write_text(json.dumps(metadata, indent=2) + "\n")


if __name__ == "__main__":
    main()
