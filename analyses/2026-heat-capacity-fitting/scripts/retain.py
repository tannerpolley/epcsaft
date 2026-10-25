"""Retain the native derivative and recovery comparisons as numerical tables."""
import csv
import json
from pathlib import Path

SET = Path(__file__).resolve().parents[1]
ROOT = SET.parents[1]
comparisons = []
for name in ("derivatives", "recovery"):
    path = ROOT / "build" / f"issue205-{name}.log"
    for line in path.read_text().splitlines():
        if line.startswith("EPCSAFT_COMPARE_JSON "):
            row = json.loads(line.split(" ", 1)[1])
            if row["comparison_id"].startswith(("cp-a/", "cp-b/")):
                comparisons.append(row)
rows, identities, recovery, weak = [], [], [], []
summary = {"base": "d47d0475", "wheel_sha256": json.loads(
    (ROOT / "build/environment-wheel/preparation.json").read_text())["wheel_sha256"]}
for row in comparisons:
    key, inputs = row["comparison_id"], row["inputs"]
    error = abs(row["model"] - row["reference"])
    if key.endswith(("/fine", "/coarse", "/richardson")):
        rows.append(dict(comparison=key, exact=row["model"], finite_difference=row["reference"],
                         absolute_error=error, criterion=row["atol"], signal=inputs["fd_signal"],
                         step_native=inputs["step"], passed=row["verdict"]))
        if key in ("cp-a/a2-298/2/fine", "cp-a/a3/2/fine"):
            weak.append(dict(comparison=key, signal=inputs["fd_signal"],
                fine_minus_coarse=inputs["fine_minus_coarse"], richardson_correction=inputs["richardson_correction"],
                exact_minus_richardson_abs=inputs["richardson_error"],
                fine_agreement_error=error, relative_criterion=row["atol"], passed=row["verdict"]))
    elif key.startswith("cp-b/value/"):
        recovery.append(dict(coordinate=int(key.rsplit("/", 1)[1]), recovered=row["model"], known=row["reference"],
            relative_error=error/abs(row["reference"]), relative_criterion=row["rtol"], passed=row["verdict"]))
    elif key == "cp-b/converged":
        summary["recovery"] = inputs
    else:
        identities.append(dict(comparison=key, quantity=row["quantity"], model=row["model"], reference=row["reference"],
            absolute_error=error, absolute_criterion=row["atol"], relative_criterion=row["rtol"],
            density_lift_criterion_ratio=inputs.get("lift_criterion_ratio", ""), passed=row["verdict"]))
for name, values in (("derivatives.csv", rows), ("identities.csv", identities),
                     ("recovery.csv", recovery), ("weak-signals.csv", weak)):
    assert values, name
    with (SET / "results" / name).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(values[0]))
        writer.writeheader()
        writer.writerows(values)
summary["comparisons"] = len(comparisons)
summary["all_passed"] = all(row["verdict"] for row in comparisons)
(SET / "results/native-summary.json").write_text(json.dumps(summary, indent=2)+"\n")
assert summary["all_passed"]
print(summary)
