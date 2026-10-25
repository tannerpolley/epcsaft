"""Caller-declared water-reaction extents; one native ordered attempt controller."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
import epcsaft
from epcsaft import equilibrium as eq
from reproduce import HERE, archived, load_case, problem_from_request, record


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=HERE / "results/starts/native")
    parser.add_argument("--case", type=int, choices=(1, 2), required=True)
    parser.add_argument("--seconds", type=float)
    parser.add_argument("--attempt-seconds", type=float)
    parser.add_argument("--extents", type=float, nargs="+", default=(.005, .02, .05, .1, .2, .3))
    args = parser.parse_args()
    mapping, request, anchor = load_case(args.case)
    anchor_phases = json.loads((HERE / "inputs" / ("case1-shell-only-1-0-anchor.json" if args.case == 1
                                                else "case2-eps-ion-2-anchor.json")).read_text())["phases"]
    stable_value = "0.330000" if args.case == 1 else "2.031000"
    states = list(csv.DictReader(archived()[f"refined/case{args.case}-continuation-{stable_value}-state.csv"].splitlines()))
    liquid = np.array([float(r["mole_fraction"]) for r in states if r["phase"] == "0"])
    direction = np.array(request["reaction_system"]["reaction_matrix"][0])
    problem = problem_from_request(request, anchor)
    problem.row_tolerance = 1e-10
    starts = []
    for extent in args.extents:  # balanced physical input declaration only
        guess = liquid + extent * direction
        assert min(guess) > 0 and abs(sum(guess) - 1) < 1e-12
        seed = eq.EquilibriumStart()
        seed.compositions = [list(guess), list(problem.phases[1].composition_guess)]
        seed.packing_fractions = [p["packing_fraction"] for p in anchor_phases]
        seed.pressure = anchor["pressure_pa"]
        starts.append(seed)
    model = epcsaft.Mixture(epcsaft.Parameters.from_mapping(mapping))
    limits = eq.EquilibriumWorkLimits()
    limits.seconds = args.seconds
    limits.attempt_seconds = args.attempt_seconds
    with epcsaft.trace() as capture:
        solved = eq.solve_equilibrium_starts(model, problem, starts, limits)
    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    for stale in out.glob(f"case{args.case}-water-extent-*"):
        stale.unlink()  # this native replay owns these generated attempt files
    iteration_spans = {int(event.values[0]): event.span for event in capture.events
                       if event.name == "declared_start_index"}
    for i, result in enumerate(solved.results):
        if result is not None:
            record(out, f"case{args.case}-water-extent-{args.extents[i]:g}", args.case,
                   0 if args.case == 1 else 2, result, capture, solved.attempt_seconds[i], "",
                   iteration_span=iteration_spans[i])
        print(args.case, args.extents[i], solved.causes[i], flush=True)
    with (out / f"case{args.case}-attempt-work.csv").open("w", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(("index", "extent", "seconds", "iterations", "evaluations", "hessian_calls",
                         "phase_count", "success", "state_available", "raw_rows", "classified_rows", "native_termination", "cause"))
        for i, result in enumerate(solved.results):
            writer.writerow((i, args.extents[i], solved.attempt_seconds[i],
                result.iterations if result else 0, result.evaluations if result else 0,
                result.hessian_calls if result else 0, result.phase_count if result else 0,
                bool(result and result.success),
                bool(result and result.success and result.temperature is not None and result.pressure is not None and result.amounts and result.mole_fractions),
                max(map(abs, result.residuals)) if result and result.residuals else None,
                max(map(abs, result.classified_residuals)) if result and result.classified_residuals else None,
                result.solver_status if result else None, solved.causes[i]))
    wheel = next((HERE.parents[1] / "build/environment-wheel").glob("epcsaft-*.whl"))
    (out / f"case{args.case}-summary.json").write_text(json.dumps({"selected": solved.selected,
        "extents": args.extents, "causes": solved.causes, "iterations": solved.iterations,
        "seconds": solved.seconds, "attempt_seconds": solved.attempt_seconds,
        "overrun_seconds": solved.overrun_seconds, "budget_exhausted": solved.budget_exhausted,
        "limits": {"seconds": args.seconds, "attempt_seconds": args.attempt_seconds},
        "wheel_sha256": hashlib.sha256(wheel.read_bytes()).hexdigest()}, indent=2) + "\n")


if __name__ == "__main__":
    main()
