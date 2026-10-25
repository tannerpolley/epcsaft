"""#572 global-start evidence: read retained multistart fit endpoints (no fitting) into results/start-evidence.csv."""
import collections
import csv
import itertools
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
rows = []


def endpoints(case, path, starts, cost, scale, rule, objective, names, bounds):
    usable = [s for s in starts if s.get("usable") and s.get("status") == "converged"]
    costs = [cost(s) for s in usable]
    spread = max((max(abs(x - y) / scale(j) for j, (x, y) in enumerate(zip(a["optimizer"], b["optimizer"])))
                  for a, b in itertools.combinations(usable, 2)), default=float("nan"))
    declared = "; ".join(f"{n} in [{lo:g}, {hi:g}]" for n, (lo, hi) in zip(names, bounds))
    active = collections.Counter(", ".join(f"{names[j]} at {min(bounds[j], key=lambda b: abs(b - s['optimizer'][j])):g}"
                                           for j in s["covariance"]["active_bounds"]) for s in usable)
    active = "; ".join(f"{k or 'none'} in {n}/{len(usable)}" for k, n in active.items())
    rows.append([case, path, objective, declared, len(starts), len(usable), len(starts) - len(usable), min(costs), max(costs),
                 spread, active, rule])


figiel = "2026-figiel-divalent-chlorides/results/"
stage1 = json.loads((ROOT / figiel / "stage1-native-fit.json").read_text())
endpoints("Figiel stage 1 (7 coordinates)", figiel + "stage1-native-fit.json", stage1["starts"], lambda s: s["cost"],
          lambda j: 6.0 if j == 6 else 2.0, "Ceres f=1e-12 g=1e-10 p=1e-12; 200 iterations/120 s per start; distinct if spread > 0.01",
          "weighted least squares on aqueous LiCl/NaCl/KCl ln gamma and ln phi (sigma 0.02 in log10), hydration holds",
          stage1["coordinate_order"], [(-1, 1)] * 6 + [(3, 9)])  # bounds declared in stage1_fit.py declarations()
for salt in ("CaCl2", "MgCl2"):
    stage2 = json.loads((ROOT / figiel / f"stage2-{salt}-ceres.json").read_text())
    endpoints(f"Figiel stage 2 {salt} (3 coordinates)", figiel + f"stage2-{salt}-ceres.json", stage2["results"],
              lambda s: s["final_cost"], lambda j: 1.0, "Ceres f=1e-12 g=1e-10 p=1e-12; 200 iterations/120 s per start; optimizer units",
              f"weighted least squares on aqueous {salt} ln gamma and ln phi to 4.5 mol/kg (sigma 0.02 in log10), hydration hold",
              [f"sigma {salt[:2]} (angstrom)", f"k {salt[:2]}-water", f"k {salt[:2]}-Cl"],
              [(2.4, 4.4), (-1, 1), (-1, 1)])  # bounds declared in stage2_fit.py declarations()
polar = json.loads((ROOT / "2026-polar-dipole/results/summary.json").read_text())["fits"]
for name, fit in polar.items():
    rows.append([f"Polar {name}", "2026-polar-dipole/results/summary.json", "", "", len(fit["optimizer_starts"]), "", "",
                 fit["final_cost"], fit["final_cost"], "", "", "selected endpoint only; per-start endpoints not retained"])
mea = sorted((ROOT / "2026-mea-parameter-estimation/results").glob("candidate-*-result.json"))
attempted = [json.loads(p.read_text()).get("accounting", {}).get("attempted_starts") for p in mea]
rows.append(["MEA parameter estimation candidates", "2026-mea-parameter-estimation/results/candidate-*-result.json",
             "", "", sum(a or 0 for a in attempted), "", "", "", "", "", "",
             f"{len(mea)} candidate records; attempted starts per record {attempted}; no repeated-start endpoints of one objective"])
with open(ROOT / "2026-correlated-measurement-errors/results/start-evidence.csv", "w", newline="") as out:
    writer = csv.writer(out, lineterminator="\n")
    writer.writerow(["case", "source", "objective", "bounds", "starts", "usable_converged", "unavailable", "cost_min", "cost_max",
                     "max_spread", "active_bounds", "rule"])
    writer.writerows([repr(v) if isinstance(v, float) else v for v in row] for row in rows)
