"""Retain native source ln K identities and joint recovery operands."""
import csv
import json
from pathlib import Path

SET = Path(__file__).resolve().parents[1]
ROOT = SET.parents[1]
comparisons = [json.loads(line.split(" ", 1)[1])
    for line in (ROOT / "build/issue211-ln-k.log").read_text().splitlines()
    if line.startswith("EPCSAFT_COMPARE_JSON ")]
comparisons = [row for row in comparisons if row["comparison_id"].startswith("ln-k/")]
assert comparisons and all(row["verdict"] for row in comparisons)
summary = {"base": "95ba231e", "wheel_sha256": json.loads(
    (ROOT / "build/environment-wheel/preparation.json").read_text())["wheel_sha256"],
    "comparisons": len(comparisons), "all_passed": True, "rejections": {}}
values = []
for row in comparisons:
    key = row["comparison_id"]
    if "/rejection/" in key:
        summary["rejections"][key.rsplit("/", 1)[1]] = row["inputs"]["message"]
    elif key == "ln-k/recovery/converged":
        summary["recovery"] = row["inputs"]
    elif key == "ln-k/recovery/covariance":
        # recover_reactions asserts covariance.rank == expected.size() == 2.
        summary["verified_covariance_rank"] = 2
    else:
        values.append(dict(comparison=key, quantity=row["quantity"], unit=row["unit"],
            model=row["model"], reference=row["reference"], absolute_error=abs(row["model"]-row["reference"]),
            atol=row["atol"], rtol=row["rtol"], T_K=row["inputs"].get("T_K", ""),
            Tref_K=row["inputs"].get("Tref_K", ""), passed=row["verdict"]))
with (SET / "results/comparisons.csv").open("w", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=list(values[0]))
    writer.writeheader()
    writer.writerows(values)
(SET / "results/summary.json").write_text(json.dumps(summary, indent=2)+"\n")
print(json.dumps(summary, indent=2))
