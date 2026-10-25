"""Transport native slice-four rows without recalculating thermodynamics or comparisons."""
import csv
import hashlib
import json
from pathlib import Path

root = Path.cwd()
raw = root / "build/current-calorics/slice4"
out = root / "analyses/2026-reactive-caloric-actions/results/slice4"
rows = []
for path in sorted(raw.glob("qualification-*.log")):
    test_id = ""
    for line in path.read_text().splitlines():
        if line.startswith("EPCSAFT_CASE ") and "test_id=" in line:
            test_id = line.split("test_id=")[1].split()[0]
        if line.startswith("EPCSAFT_COMPARE_JSON "):
            r = json.loads(line.split(" ", 1)[1])
            r["test_id"], r["log"] = test_id, str(path.relative_to(root))
            inputs = r.pop("inputs", {})
            r.update({"input_" + key: value for key, value in inputs.items()})
            rows.append(r)
columns = sorted({key for row in rows for key in row})
with (out / "comparisons.tsv").open("w") as f:
    writer = csv.DictWriter(f, fieldnames=columns, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows({key: json.dumps(value) if isinstance(value, (list, dict)) else value for key, value in row.items()} for row in rows)
preparation = json.loads((root / "build/environment-wheel/preparation.json").read_text())
provenance = {"base": "ec0bafa71b5022821f1774b2b9fffb06622be87c", "preparation": preparation,
              "type": "numerical verification", "review": "delivered review pending",
              "limits": "manufactured branches; no physical calorimetry, parameter adoption, global stability or manuscript promotion",
              "logs": {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(raw.glob("*.log"))},
              "key_solve_work": "build/current-calorics/slice4/key-solves/work.csv",
              "key_solve_timing": "build/current-calorics/slice4/key-solves/summary.json",
              "runtime_diff_sha256": hashlib.sha256((raw / "runtime.patch").read_bytes()).hexdigest()}
(out / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(7, 4), constrained_layout=True)
for index, group in enumerate(("F1", "F2", "F3", "F7")):
    selected = [r for r in rows if r["comparison_id"].startswith(group + "/") and "verdict" in r and "input_error_over_target" in r]
    ax.scatter([index] * len(selected), [r["input_error_over_target"] for r in selected], s=20, label=group)
ax.axhline(1, color="black", linestyle="--", linewidth=1, label="Required limit")
ax.set(ylim=(-1e-11, 1.4), xticks=range(4), xticklabels=["Independent H", "Integrals", "Nonreacting", "Fit columns"], ylabel="Absolute error / native target", title="Apparent Cp: exact calculations versus independent checks")
ax.set_yscale("symlog", linthresh=1e-10)
ax.legend(loc="upper right", ncol=3, fontsize=8)
fig.savefig(root / "analyses/2026-reactive-caloric-actions/figures/slice4-errors.png", dpi=180)
plt.close(fig)
