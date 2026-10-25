"""Best-of-three ordinary solves at the frozen #223 ionic LLE and MEA flash inputs."""

import argparse
import csv
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "analyses/2026-tpd-phase-discovery/scripts"))
import baseline as b


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("label", choices=("main", "branch"))
    parser.add_argument("--sample", type=int, choices=(0, 1, 2))
    args = parser.parse_args()
    out = ROOT / "analyses/2026-born-unavailable-states/results/inventory-bounds"
    out.mkdir(parents=True, exist_ok=True)
    wanted = {"case-6/figure-2-point-7", "case-7/flash-5918.5286-Pa"}
    rows, inputs, states = [], {}, {}
    for case, model, problem in b.declared_cases():
        if case not in wanted:
            continue
        inputs[case] = dict(components=model.component_ids, **b.problem_input(problem))
        for sample in range(3) if args.sample is None else [args.sample]:  # Cost samples, not search.
            start = time.perf_counter()
            result = b.eq.solve_equilibrium(model, problem)
            elapsed = time.perf_counter() - start
            assert result.success and result.requested_tolerance_met, result.message
            rows.append((case, sample, elapsed, result.iterations, result.evaluations,
                         result.hessian_calls, result.convergence.residual, result.objective))
            states[f"{case}/{sample}"] = dict(amounts=result.amounts,
                mole_fractions=result.mole_fractions, residuals=result.residuals)
    append = args.sample is not None and args.sample > 0
    with (out / f"cost-{args.label}.csv").open("a" if append else "w", newline="") as stream:
        writer = csv.writer(stream)
        if not append:
            writer.writerow(("case", "sample", "seconds", "iterations", "evaluations",
                             "hessian_calls", "residual", "objective"))
        writer.writerows(rows)
    direct = json.loads(importlib.metadata.distribution("epcsaft").read_text("direct_url.json"))
    wheel = Path(direct["url"].removeprefix("file://"))
    provenance = dict(wheel_sha256=hashlib.sha256(wheel.read_bytes()).hexdigest(),
        threads={key: os.environ.get(key) for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS")},
        niceness=os.getpriority(os.PRIO_PROCESS, 0), inputs=inputs,
        reactions=b.INPUTS["case-7/reactions"], states=states,
        packet_sha256={str(p.relative_to(ROOT)): b.digest(p) for p in (b.KH_PACKET, b.MEA_PACKET)})
    identity = out / f"cost-{args.label}.json"
    if append:
        prior = json.loads(identity.read_text())
        assert prior["wheel_sha256"] == provenance["wheel_sha256"]
        provenance["states"] = {**prior["states"], **states}
    identity.write_text(json.dumps(provenance, separators=(",", ":")) + "\n")
    for case in sorted(wanted):
        samples = [r for r in rows if r[0] == case]
        print(args.label, case, "best seconds", min(r[2] for r in samples),
              "iterations", [r[3] for r in samples], flush=True)


if __name__ == "__main__":
    main()
