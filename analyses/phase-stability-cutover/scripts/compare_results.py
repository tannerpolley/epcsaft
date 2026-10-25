"""Compare retained values only; no Engine evaluation or stability search."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--native-dir", type=Path, required=True)
parser.add_argument("--output-dir", type=Path, required=True)
args = parser.parse_args()
args.native_dir = args.native_dir.resolve()
args.output_dir.mkdir(parents=True, exist_ok=False)
identity = json.loads((args.native_dir / "summary.json").read_text())
tables = HERE / identity["canonical_directory"] if "canonical_directory" in identity else args.native_dir
read = lambda name: list(csv.DictReader((tables / name).open()))
comparisons = read("before-after.csv")
states = {f"{r['case_id']}/{r['declared_topology']}": r for r in read("current-engine.csv")}
rows = []


def retain(group, case, index, old, new, criterion, reason="", one_sided=False):
    difference = "" if old == "" or new == "" else float(new) - float(old)
    error = "" if difference == "" else max(0.0, difference) if one_sided else abs(difference)
    rows.append(dict(group=group, case=case, index=index, before=old, after=new,
                     difference=difference, error=error, criterion=criterion,
                     passes="" if error == "" else error <= criterion, reason=reason))


for row in comparisons:
    case = row["case"]
    state = states[case]
    total = sum(json.loads(state["feed_mol"]))
    old_g = float(row["before_g_over_rt_mol"]) / total
    new_g = float(row["after_g_over_rt_mol"]) / total
    retain("candidate_G_over_RT_per_mol", case, 0, old_g, new_g, 1e-7)
    retain("candidate_objective_increase", case, 0, old_g, new_g,
           1e-8 * max(1.0, old_g), one_sided=True)
    for name in ("dmu_max_abs", "material_max_abs"):
        retain(name, case, 0, row[f"before_{name}"], row[f"after_{name}"], 1e-7)
    for p, (old_x, new_x) in enumerate(zip(json.loads(row["before_phase_mole_fractions"]),
                                         json.loads(row["after_phase_mole_fractions"]), strict=True)):
        for i, (old, new) in enumerate(zip(old_x, new_x, strict=True)):
            retain("candidate_mole_fraction", case, f"{p}/{i}", old, new, 1e-7,
                   "Historical fractions rounded to eight decimal places; different solved coordinates.")
    for p, (old, new) in enumerate(zip(json.loads(row["before_phase_amount_fractions"]),
                                      json.loads(row["after_phase_amounts_mol"]), strict=True)):
        retain("candidate_amount_fraction", case, p, old, new / total, 1e-7)
    retain("fixed_candidate_prediction", case, 0, "", "", 1e-7,
           "No historical raw coordinates, densities or chemical plane retained; not an identical-point replay.")
    retain("candidate_pressure_relative", case, 0, "", state["P_Pa"], 1e-7, "Fixed pressure declaration only; historical central pressure readout absent.")
    retain("candidate_temperature_absolute", case, 0, "", state["T_K"], 1e-6, "Fixed temperature declaration only; historical returned temperature readout absent.")
    for p, density in enumerate(json.loads(state["molar_densities_mol_m3"])):
        retain("candidate_density_relative", case, p, "", density, 1e-7,
               "Historical table has no density column.")
    for branch in ("liquid", "vapor"):
        old, new = row[f"before_tpd_min_{branch}"], row[f"after_tpd_min_{branch}"]
        reason = "Different finite start sets and solved candidates; minimum comparison only."
        if new == "":
            reason = state["stability_verdict"]
        retain("TPD_minimum", case, branch, old, new, 1e-7, reason)
        retain("fixed_trial_prediction", case, branch, old, "", 1e-7,
               reason if new == "" else "No supplied-composition TPD evaluation in stability_test; historical argmin rounded to six decimals.")
        old_x = json.loads(row[f"before_argmin_{branch}"])
        new_x = json.loads(row[f"after_argmin_{branch}"]) if new else [""] * len(old_x)
        for i, (old, new) in enumerate(zip(old_x, new_x, strict=True)):
            retain("TPD_argmin_secondary", case, f"{branch}/{i}", old, new, 1e-7,
                   "Rounded historical coordinates; nonunique minima can select different present phases." if new != "" else reason)
for row in read("declared-residuals.csv"):
    retain("classified_row_closure", row["case"], f"{row['kind']}/{row['phase']}/{row['index']}",
           0.0, row["classified"], 1e-10)
for row in read("quantities.csv"):
    if row["quantity"] in ("charge_residual", "normalization_residual"):
        retain(row["quantity"], row["case"], row["entry"], 0.0, row["value"], 1e-12)
with (args.output_dir / "comparison.csv").open("w", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=rows[0], lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
groups = {}
for row in rows:
    group = groups.setdefault(row["group"], dict(count=0, unavailable=0, failed=0, maximum_error=None, maximum_error_over_criterion=None))
    group["count"] += 1
    if row["error"] == "":
        group["unavailable"] += 1
    else:
        group["failed"] += not row["passes"]
        group["maximum_error"] = max(group["maximum_error"] or 0.0, row["error"])
        group["maximum_error_over_criterion"] = max(group["maximum_error_over_criterion"] or 0.0, row["error"] / row["criterion"])
sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
summary = dict(native_directory=str(args.native_dir.relative_to(HERE)), wheel_sha256=identity["identity"]["wheel_sha256"],
               canonical_directory=str(tables.relative_to(HERE)), inputs={p.name: sha(p) for p in tables.glob("*.csv")},
               script_sha256=sha(Path(__file__)), comparison_sha256=sha(args.output_dir / "comparison.csv"),
               historical_identity=identity["identity"]["historical_identity"], groups=groups,
               completed=False, limits=identity["claim_limits"],
               missing="Fixed-candidate predictions before optimization and common supplied-composition TPD evaluations remain unverified.")
(args.output_dir / "summary.json").write_text(json.dumps(summary, separators=(",", ":")) + "\n")
print(json.dumps(groups, indent=2))
