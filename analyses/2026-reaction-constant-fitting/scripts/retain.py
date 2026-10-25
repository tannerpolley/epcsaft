"""Retain #187 numerical comparisons from the native and Python evidence logs."""
import argparse
import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
RESULTS = Path(__file__).resolve().parents[1] / "results"


def records(path: Path):
    for line in path.read_text().splitlines():
        if line.startswith("EPCSAFT_COMPARE_JSON "):
            yield json.loads(line.split(" ", 1)[1])


def table(name: str, rows: list[dict]):
    with (RESULTS / name).open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("native", type=Path)
parser.add_argument("exchange", type=Path)
args = parser.parse_args()
derivatives, regularity, identities, recovery, rejections = [], [], [], [], []
summary = {"base": "fe00aa6a", "wheel_sha256": json.loads(
    (ROOT / "build/environment-wheel/preparation.json").read_text())["wheel_sha256"],
    "native_log": str(args.native), "exchange_log": str(args.exchange), "recovery": {}}
for row in records(args.native):
    key, inputs = row["comparison_id"], row["inputs"]
    if not key.startswith("reaction-"):
        continue
    if key.endswith(("/fine", "/coarse")):
        error = abs(row["model"] - row["reference"])
        derivatives.append(dict(comparison=key, exact=row["model"], central_difference=row["reference"],
            absolute_error=error, allowed_error=row["atol"], signal=inputs["fd_signal"],
            prediction=inputs["prediction"], optimizer_step=inputs["u_step"], passed=row["verdict"]))
    elif key.endswith("/regular"):
        regularity.append(dict(case=key.removesuffix("/regular"), **inputs, passed=row["verdict"]))
    elif "/value/" in key:
        recovery.append(dict(case=key.rsplit("/value/", 1)[0], coordinate_row=key.rsplit("/", 1)[1],
            recovered=row["model"], known=row["reference"], absolute_error=abs(row["model"]-row["reference"]),
            atol=row["atol"], rtol=row["rtol"], passed=row["verdict"]))
    elif key.endswith("/converged"):
        summary["recovery"][key.removesuffix("/converged")] = inputs
    elif key.startswith("reaction-preparation/"):
        summary["preparation"] = inputs
        identities.append(dict(comparison=key, quantity=row["quantity"], model=row["model"],
            reference=row["reference"], absolute_error=abs(row["model"]-row["reference"]),
            atol=row["atol"], rtol=row["rtol"], passed=row["verdict"]))
    elif key.startswith("reaction-rejection/"):
        rejections.append(dict(comparison=key, **inputs, passed=row["verdict"]))
    elif row["quantity"] != "verdict":
        identities.append(dict(comparison=key, quantity=row["quantity"], model=row["model"],
            reference=row["reference"], absolute_error=abs(row["model"]-row["reference"]),
            atol=row["atol"], rtol=row["rtol"], passed=row["verdict"]))
for name, rows in (("derivatives.csv", derivatives), ("regularity.csv", regularity),
                   ("identities.csv", identities), ("recovery.csv", recovery)):
    table(name, rows)
exchange, laws = [], []
for row in records(args.exchange):
    if row["comparison_id"].startswith("reaction-c/"):
        if "/law/" in row["comparison_id"]:
            laws.append(dict(comparison=row["comparison_id"], log_activity_ratio=row["model"],
                known_lnK=row["reference"], balance_mol=row["inputs"]["balance_mol"]))
            continue
        criterion = row["criterion"]
        error = abs(row["model"] - row["reference"])
        exchange.append(dict(comparison=row["comparison_id"], engine=row["model"], reference=row["reference"],
            absolute_error=error, atol=criterion["atol"], rtol=criterion["rtol"], passed=error <= criterion["atol"]))
        summary["exchange"] = row["inputs"]
table("exchange.csv", exchange)
table("exchange-laws.csv", laws)
summary["rejections"] = rejections
(RESULTS / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
