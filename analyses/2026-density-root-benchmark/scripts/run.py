#!/usr/bin/env python3
"""Run the identical declarations against an explicitly selected installed Engine."""

import argparse
import csv
import hashlib
import json
import math
import os
import statistics
import subprocess
import sys
import time
from pathlib import Path

from verification import application_check, density_check

HERE = Path(__file__).resolve().parents[1]
ROOT = HERE.parents[1]
ENV = os.environ | {
    key: "1"
    for key in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "BLIS_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
        "NUMEXPR_NUM_THREADS",
        "JULIA_NUM_THREADS",
        "RAYON_NUM_THREADS",
    )
}
ENV.update(PYTHONNOUSERSITE="1", PYTHONPATH="")


def csv_rows(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def write_csv(path, rows):
    assert rows
    with path.open("w", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=list(dict.fromkeys(k for row in rows for k in row)),
            lineterminator="\n",
        )
        w.writeheader()
        w.writerows(rows)


def percentile(samples, p):
    values = sorted(samples)
    index = (len(values) - 1) * p
    lo = math.floor(index)
    hi = math.ceil(index)
    return values[lo] + (values[hi] - values[lo]) * (index - lo)


def validate(q, records, document, references):
    """All warm-up, timed and matching work-replay returns have identical acceptance."""
    maximum = 0.0
    failures = []
    checked = 0
    for replica in records:
        for mode in (replica["primary"], replica["repeated_identical"]):
            if mode is None:
                continue
            failures.extend(mode["errors"])
            for record in mode["records"]:
                checked += 1
                try:
                    if q["tier"] in ("A", "B"):
                        delta = density_check(
                            q["states"][record["state_index"]],
                            record["values"],
                            document["gas_constant_J_mol_K"],
                            fixed=record["kind"] == "fixed",
                        )
                    else:
                        delta = application_check(
                            record["values"], references[q["id"]][record["state_index"]]
                        )
                    maximum = max(maximum, delta)
                except (AssertionError, TypeError, ValueError) as error:
                    failures.append(str(error))
    return not failures, maximum, checked, failures


def metrics(mode, density):
    records = mode["records"]
    kinds = ["pressure", "fixed"] if density else ["application"]
    timed = {
        kind: [r["us"] for r in records if r["stage"] == "timed" and r["kind"] == kind]
        for kind in kinds
    }
    primary = (
        [p - f for p, f in zip(timed["pressure"], timed["fixed"], strict=True)]
        if density
        else timed["application"]
    )
    work_records = [
        r for r in records if r["stage"] == "work" and r["kind"] == kinds[0]
    ]
    assert len(work_records) == len(primary)
    return primary, timed, work_records


def summary(rows):
    lines = [
        "| Case | Primary mode | N/process | Correct | Primary median µs | p10–p90 µs | Raw P / fixed µs | EOS evaluations/call |",
        "|---|---|---:|:---:|---:|---:|---:|---:|",
    ]
    for r in rows:
        raw = (
            f"{r['pressure_call_median_us']:.3f} / {r['fixed_call_median_us']:.3f}"
            if r["tier"] != "C"
            else "—"
        )
        lines.append(
            f"| {r['id']} | {r['mode']} | {r['N']} | {r['correct']} | {r['primary_median_us']:.3f} | {r['primary_p10_us']:.3f}–{r['primary_p90_us']:.3f} | {raw} | {r['eos_evaluations']:.2f} |"
        )
    return "\n".join(lines) + "\n"


def compare(a, b):
    left = {r["id"]: r for r in csv_rows(a)}
    right = {r["id"]: r for r in csv_rows(b)}
    assert left.keys() == right.keys(), "case sets differ"
    print(
        "| Case | B/A primary time | B/A raw pressure/application time | B/A total EOS evaluations |"
    )
    print("|---|---:|---:|---:|")
    for id, x in left.items():
        y = right[id]
        assert x["case_sha256"] == y["case_sha256"], "case declarations differ"
        assert x["correct"] == y["correct"] == "True", "compare requires correct rows"

        def ratio(k):
            return (
                f"{float(y[k]) / float(x[k]):.4f}"
                if float(x[k])
                else ("1.0000" if not float(y[k]) else "undefined")
            )

        raw = "pressure_call_median_us" if x["tier"] != "C" else "primary_median_us"
        print(
            f"| {id} | {ratio('primary_median_us')} | {ratio(raw)} | {ratio('eos_evaluations')} |"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--python", type=Path)
    parser.add_argument("--label")
    parser.add_argument("--compare", nargs=2, type=Path, metavar=("A.csv", "B.csv"))
    parser.add_argument(
        "--record-main",
        action="store_true",
        help="regenerate reference with clean, identical pinned main Engine sources",
    )
    args = parser.parse_args()
    if args.compare:
        compare(*args.compare)
        return
    if not args.python or not args.label:
        parser.error("--python and --label are required")
    if Path(args.label).name != args.label or args.label in (".", ".."):
        parser.error("--label must be one directory name")
    sys.path.insert(0, str(ROOT / "tools"))
    from engine_identity import prepared_python
    python, worker_env = prepared_python(args.python, ENV)
    prepared = json.loads(subprocess.check_output([str(python), str(ROOT / "tools/engine_identity.py")], env=worker_env))
    if not prepared["source_revision"] or not prepared["meson_configuration_sha256"] or not prepared["lock_sha256"]:
        parser.error("selected Engine must identify its prepared source, Meson build, and lock")
    document = json.loads((HERE / "inputs/cases.json").read_text())
    sha = hashlib.sha256((HERE / "inputs/cases.json").read_bytes()).hexdigest()
    references = json.loads((HERE / "results/application-reference.json").read_text())
    competitors = csv_rows(HERE / "results/competitors.csv")
    out = ROOT / "build/density-benchmark" / args.label
    out.mkdir(parents=True, exist_ok=True)
    processes = []
    logs = []
    rows = []
    failures = []
    samples = {}
    if os.getpriority(os.PRIO_PROCESS, 0) > 10:
        parser.error("inherited niceness exceeds 10")
    begun = time.perf_counter()
    try:
        for replica in range(3):
            log = (out / f"process-{replica}.log").open("w")
            logs.append(log)
            processes.append(
                subprocess.Popen(
                    [
                        "nice",
                        "-n",
                        str(10 - os.getpriority(os.PRIO_PROCESS, 0)),
                        str(python),
                        str(HERE / "scripts/worker.py"),
                    ],
                    env=worker_env,
                    stdin=subprocess.PIPE,
                    stdout=subprocess.PIPE,
                    stderr=log,
                    text=True,
                    bufsize=1,
                )
            )
        identities = [json.loads(process.stdout.readline()) for process in processes]
        assert (
            len(
                {
                    (v["engine_identity"], str(v["engine_detail"]), v["engine_tree"])
                    for v in identities
                }
            )
            == 1
        ), "worker Engines differ"
        assert identities[0]["engine_detail"] == prepared, "worker image differs from prepared Engine"
        if args.record_main:
            assert (
                prepared["source_diff_sha256"] == hashlib.sha256(b"").hexdigest()
                and identities[0]["engine_tree"] == document["reference_engine_tree"]
            ), "reference requires clean pinned Engine tree"
        (out / "environment.json").write_text(json.dumps(identities, indent=2) + "\n")
        for q in document["cases"]:
            print(q["id"], flush=True)
            records = []
            for replica in [(len(rows) + i) % 3 for i in range(3)]:
                process = processes[replica]
                process.stdin.write(json.dumps(dict(id=q["id"], rep=replica)) + "\n")
                process.stdin.flush()
                line = process.stdout.readline()
                if not line:
                    raise RuntimeError(
                        f"worker {replica} stopped; see {out}/process-{replica}.log"
                    )
                record = json.loads(line)
                if "error" in record:
                    raise RuntimeError(
                        f"{q['id']}: {record['error_type']}: {record['error']}"
                    )
                records.append(record)
            samples[q["id"]] = records
            (out / (q["id"].replace("/", "_") + ".json")).write_text(
                json.dumps(records) + "\n"
            )
            ok, delta, checked, errors = validate(q, records, document, references)
            density = q["tier"] != "C"
            primary = []
            raw = {}
            counts = []
            secondary = []
            for record in records:
                times, timed, work_records = metrics(record["primary"], density)
                primary.extend(times)
                counts.extend(r["counts"] for r in work_records)
                for kind, values in timed.items():
                    raw.setdefault(kind, []).extend(values)
                if record["repeated_identical"] is not None:
                    secondary.extend(metrics(record["repeated_identical"], density)[0])
            assert len(primary) == 3 * q["N"]
            row = dict(
                id=q["id"],
                tier=q["tier"],
                mode=q["mode"],
                N=q["N"],
                processes=3,
                warmups=q["warmups"],
                primary_metric="marginal pressure request"
                if density
                else "raw application call",
                correct=ok,
                relative_difference=delta,
                checked_returns=checked,
                primary_median_us=statistics.median(primary),
                primary_p10_us=percentile(primary, 0.1),
                primary_p90_us=percentile(primary, 0.9),
                repeated_identical_label="repeated identical input (rewards caching)",
                repeated_identical_median_us=statistics.median(secondary or primary),
                work_mode="fresh matching-history replay of primary",
                **{k: statistics.mean(c[k] for c in counts) for k in counts[0]},
                values=json.dumps(
                    next(
                        r["values"]
                        for r in records[0]["primary"]["records"]
                        if r["stage"] == "timed"
                    ),
                    separators=(",", ":"),
                    sort_keys=True,
                ),
                case_sha256=sha,
                reference_sha=document["main_sha"],
                python=str(python),
                engine_provenance=json.dumps(
                    [identities[0]["engine_identity"], identities[0]["engine_detail"]]
                ),
                engine_tree=identities[0]["engine_tree"],
                errors=json.dumps(errors),
            )
            if density:
                row["pressure_call_median_us"] = statistics.median(raw["pressure"])
                row["pressure_call_p10_us"] = percentile(raw["pressure"], 0.1)
                row["pressure_call_p90_us"] = percentile(raw["pressure"], 0.9)
                row["fixed_call_median_us"] = statistics.median(raw["fixed"])
                row["fixed_call_p10_us"] = percentile(raw["fixed"], 0.1)
                row["fixed_call_p90_us"] = percentile(raw["fixed"], 0.9)
                fixed_counts = [
                    r["counts"]["eos_evaluations"]
                    for rec in records
                    for r in rec["primary"]["records"]
                    if r["stage"] == "work" and r["kind"] == "fixed"
                ]
                row["fixed_eos_evaluations"] = statistics.mean(fixed_counts)
                for lib in ("FeOs", "Clapeyron"):
                    c = next(
                        v
                        for v in competitors
                        if v["id"] == q["id"] and v["library"] == lib
                    )
                    for k in (
                        "status",
                        "rho_mol_m3",
                        "median_us",
                        "p10_us",
                        "p90_us",
                        "max_relative_density_difference",
                        "accuracy_within_1e10",
                    ):
                        row[lib + "_" + k] = c.get(k, "")
            rows.append(row)
            write_csv(out / "rows.csv", rows)
            if not ok:
                failures.append(q["id"])
    finally:
        for process in processes:
            process.stdin.close()
        for process in processes:
            process.wait()
        for log in logs:
            log.close()
    table = summary(rows)
    (out / "summary.md").write_text(
        table
        + f"\nElapsed: {time.perf_counter() - begun:.3f} s. Raw competitor times compare with raw Engine pressure calls, not marginal subtraction.\n"
    )
    print(table)
    if failures:
        raise RuntimeError(f"correctness failures: {failures}")
    if args.record_main:
        write_csv(HERE / f"results/main-{document['main_sha'][:8]}.csv", rows)
        (HERE / "results/main-samples.json").write_text(
            json.dumps(
                dict(environment=identities, cases=samples), separators=(",", ":")
            )
            + "\n"
        )


if __name__ == "__main__":
    main()
