"""Compare retained Born rows and current native state diagnostics."""

import csv
import hashlib
import json
from pathlib import Path

root = Path(__file__).resolve().parents[1]
new = root / "results/diagnostics/native"
comparisons = []
summary_path = new / "comparison-summary.json"
previous = json.loads(summary_path.read_text()) if summary_path.exists() else {}
reused = {
    r["replay_file"]: r["retained_file"]
    for r in previous.get("reused_identical_files", [])
}


def retained(path):
    return path if path.exists() else root / reused[str(path.relative_to(root))]


def replay_paths(pattern):
    return sorted(
        set((new / "points").rglob(pattern))
        | {
            root / name
            for name in reused
            if Path(name).match(pattern) and "points" in Path(name).parts
        }
    )


def compare(old_path, new_path, fields, units, key):
    old = list(csv.DictReader(old_path.open()))
    after = list(csv.DictReader(retained(new_path).open()))
    before = {tuple(r[k] for k in key): r for r in old}
    for r in after:
        a = before[tuple(r[k] for k in key)]
        for field in fields:
            x, y = float(a[field]), float(r[field])
            comparisons.append(
                {
                    "historical_file": str(old_path.relative_to(root)),
                    "native_file": str(new_path.relative_to(root)),
                    "point": ";".join(r[k] for k in key),
                    "quantity": field,
                    "unit": units[field],
                    "before": x,
                    "after": y,
                    "change": y - x,
                }
            )


compare(
    root / "results/global-scan/row-verification.csv",
    new / "global-scan/row-verification.csv",
    ["original_row_max", "carbon_residual"],
    {"original_row_max": "1", "carbon_residual": "1"},
    ["case", "scan", "attempt"],
)
fields = [
    "row_max",
    "rho_mol_m3",
    "dP_drho",
    "born_over_RT",
    "debye_huckel_over_RT",
    "bulk_relative_permittivity",
]
units = {f: "1" for f in fields}
units.update(rho_mol_m3="mol/m3", dP_drho="Pa m3/mol")
compare(
    root / "results/turning-verification.csv",
    new / "turning-verification.csv",
    fields,
    units,
    ["case", "path", "point", "phase"],
)
for path in sorted(replay_paths("*-point-diagnostics.json")):
    matches = [root / "results" / path.relative_to(new / "points")]
    a, b = (json.loads(retained(p).read_text()) for p in (matches[0], path))
    for field in [
        "row_recheck_max",
        "unscaled_condition",
        "equilibrated_condition",
        "equilibrated_singular_min",
        "equilibrated_singular_max",
        "equilibrated_rank",
    ]:
        comparisons.append(
            {
                "historical_file": str(matches[0].relative_to(root)),
                "native_file": str(path.relative_to(root)),
                "point": path.stem,
                "quantity": field,
                "unit": "1",
                "before": a[field],
                "after": b[field],
                "change": b[field] - a[field],
            }
        )
for path in sorted(replay_paths("*-terms.csv")):
    matches = [root / "results" / path.relative_to(new / "points")]
    fs = [
        "born_over_RT",
        "debye_huckel_over_RT",
        "bulk_relative_permittivity",
        "net_charge_mol",
    ]
    u = {f: "1" for f in fs}
    u["net_charge_mol"] = "mol"
    compare(matches[0], path, fs, u, ["phase"])
with (new / "before-after.csv").open("w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(comparisons[0]), lineterminator="\n")
    w.writeheader()
    w.writerows(comparisons)
summary = {
    "wheel_sha256": hashlib.sha256(
        next(
            (root.parents[1] / "build/environment-wheel").glob("epcsaft-*.whl")
        ).read_bytes()
    ).hexdigest(),
    "historical_source_commit": "e12a5ff84626a2b4bef49e746f155e1df1f9bd73",
    "comparison_count": len(comparisons),
    "changed_values": sum(r["change"] != 0 for r in comparisons),
    "maximum_absolute_change_by_quantity": {
        q: max(abs(r["change"]) for r in comparisons if r["quantity"] == q)
        for q in sorted({r["quantity"] for r in comparisons})
    },
    "global_scan": json.loads(
        retained(new / "global-scan/verification-summary.json").read_text()
    ),
    "turning_phase_states": len(
        list(csv.DictReader(retained(new / "turning-verification.csv").open()))
    ),
    "point_diagnostics": len(replay_paths("*-point-diagnostics.json")),
    "term_diagnostics": "Accepted native states only; failed-iterate diagnostics retired; historical files unchanged",
    "state_meaning": "Turning checks use native pressure-closed states for converged historical roots, not exact compiled-coordinate state readouts",
    "limit": "Numerical verification, not physical validation",
}
# Reuse byte-identical prior values; the comparison table preserves the replay measurements.
original = {}
for path in (root / "results").rglob("*"):
    if path.is_file() and new not in path.parents:
        original.setdefault(hashlib.sha256(path.read_bytes()).hexdigest(), path)
for path in sorted(new.rglob("*")):
    if not path.is_file() or path == summary_path:
        continue
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest in original:
        reused[str(path.relative_to(root))] = str(original[digest].relative_to(root))
        path.unlink()
summary["reused_identical_files"] = [
    {
        "replay_file": name,
        "retained_file": target,
        "sha256": hashlib.sha256((root / target).read_bytes()).hexdigest(),
    }
    for name, target in sorted(reused.items())
]
summary_path.write_text(json.dumps(summary, indent=2) + "\n")
print(
    json.dumps(
        {k: v for k, v in summary.items() if k != "reused_identical_files"}, indent=2
    )
)
print(f"reused {len(reused)} byte-identical files")
