"""Render retained #374 results and their summary statistics; never run the Engine."""
import csv
import math
import statistics as st
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
CASES = ("aq-literal", "aq-constant-sigma", "aq-cr1-cross", "aq-organic-eps", "both-organic-eps")


def load(name):
    return list(csv.DictReader((RESULTS / name).open()))


comparison = [{k: float(v) for k, v in row.items()} for row in load("comparison.csv")]
summary = []
for case in CASES:
    sweep = load(f"sweep-{case}.csv")
    d = [r[f"pH_{case}_minus_esi"] for r in comparison]
    k = [r[f"K_{case}_rel_esi"] for r in comparison]
    row = dict(
        case=case, solves=len(sweep), success=sum(r["success"] == "True" for r in sweep),
        tolerance_met=sum(r["tolerance_met"] == "True" for r in sweep),
        iterations_median=st.median(float(r["iterations"]) for r in sweep),
        iterations_max=max(float(r["iterations"]) for r in sweep),
        evaluations_max=max(float(r["evaluations"]) for r in sweep),
        hessian_calls_max=max(float(r["hessian_calls"]) for r in sweep),
        solve_ms_median=st.median(float(r["solve_ms"]) for r in sweep),
        classified_residual_max=max(float(r["classified_residual"]) for r in sweep),
        reduced_hessian_min=min(float(r["reduced_hessian_minimum"]) for r in sweep),
        charge_aq_mol_max=max(abs(float(r["charge_aq_mol"])) for r in sweep),
        charge_org_mol_max=max(abs(float(r["charge_org_mol"])) for r in sweep),
        affinity_density_HAc_aq_max=max(abs(float(r["affinity_density_HAc_aq"])) for r in sweep),
        affinity_density_H2O_aq_max=max(abs(float(r["affinity_density_H2O_aq"])) for r in sweep),
        contraction_recomputed_minus_compiled_max=max(abs(float(r[f"contraction_recomputed_minus_compiled_{n}"]))
                                                      for r in sweep for n in ("HAc", "H2O")),
        acid_law_residual_max=max(abs(float(r["acid_law_residual_molality"])) for r in sweep),
        reference_terminal_change_max=max(float(r["reference_terminal_change"]) for r in sweep),
        dpH_esi_mean=st.mean(d), dpH_esi_min=min(d), dpH_esi_max=max(d),
        dpH_within_0p005=sum(abs(v) <= 0.005 for v in d),
        dK_esi_mean=st.mean(k), dK_esi_min=min(k), dK_esi_max=max(k), dK_within_1pct=sum(abs(v) <= 0.01 for v in k),
        AARD_pH_exp_percent=100 * st.mean(abs(r[f"pH_{case}"] - r["pH_exp"]) / r["pH_exp"] for r in comparison),
        bias_pH_exp=st.mean(r[f"pH_{case}"] - r["pH_exp"] for r in comparison),
        coarse_fine_pH_max=max(r[f"pH_{case}_coarse_fine_difference"] for r in comparison),
        coarse_fine_K_max=max(r[f"K_{case}_coarse_fine_difference"] for r in comparison),
    )
    if case.startswith("both"):
        row.update(
            affinity_density_HAc_org_max=max(abs(float(r["affinity_density_HAc_org"])) for r in sweep),
            affinity_density_H2O_org_max=max(abs(float(r["affinity_density_H2O_org"])) for r in sweep),
            mean_ionic_transfer_max=max(abs(float(r["mean_ionic_transfer_residual"])) for r in sweep),
            charge_relative_org_max=max(abs(float(r["charge_relative_org"])) for r in sweep),
            x_ion_org_reconstructed_min=min(float(r["x_ion_org_reconstructed"]) for r in sweep),
            x_ion_org_reconstructed_max=max(float(r["x_ion_org_reconstructed"]) for r in sweep),
            pH_org_reconstructed_min=min(float(r["pH_org_reconstructed"]) for r in sweep),
            pH_org_reconstructed_max=max(float(r["pH_org_reconstructed"]) for r in sweep),
            negative_reduced_hessian=sum(float(r["reduced_hessian_minimum"]) < 0 for r in sweep),
        )
    summary.append(row)
baseline = load("baseline.csv")
neutral = sorted((math.log(float(r["x_acid_aq"])), float(r["K_aq_org"])) for r in baseline if r["case"] == "nonreactive")


def neutral_k(x):
    lx = math.log(x)
    (a, fa), (b, fb) = next(p for p in zip(neutral, neutral[1:]) if p[0][0] <= lx <= p[1][0])
    return fa + (fb - fa) * (lx - a) / (b - a)


k = [neutral_k(r["x_acetic_acid_aq"]) / r["K_esi_model"] - 1 for r in comparison]
binary = next(r for r in baseline if r["case"] == "water-toluene")
summary.append(dict(case="nonreactive-literal", solves=len(neutral) + 1,
                    success=sum(r["success"] == "True" for r in baseline),
                    x_toluene_aq=float(binary["x_toluene_aq"]), x_water_org=float(binary["x_water_org"]),
                    reduced_hessian_min=min(float(r["reduced_hessian_minimum"]) for r in baseline),
                    dK_esi_mean=st.mean(k), dK_esi_min=min(k), dK_esi_max=max(k), dK_within_1pct=sum(abs(v) <= 0.01 for v in k)))
summary.append(dict(case="esi-model", AARD_pH_exp_percent=100 * st.mean(abs(r["pH_esi_model"] - r["pH_exp"]) / r["pH_exp"] for r in comparison),
                    bias_pH_exp=st.mean(r["pH_esi_model"] - r["pH_exp"] for r in comparison)))
with (RESULTS / "summary.csv").open("w", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=list(dict.fromkeys(k for r in summary for k in r)), lineterminator="\n")
    writer.writeheader()
    writer.writerows(summary)

plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
fig, (a, b) = plt.subplots(1, 2, figsize=(9.0, 4.3), layout="constrained")
for case, label, style in (("aq-literal", "Engine, ESI water-only permittivity", "-"),
                           ("both-organic-eps", "Engine, ions in both liquids, CRC organic permittivity", "--")):
    sweep = sorted(load(f"sweep-{case}.csv"), key=lambda r: float(r["x_acid_total_aq"]))
    x = [float(r["x_acid_total_aq"]) for r in sweep]
    a.plot(x, [float(r["pH_aq"]) for r in sweep], style, color="#0072B2", label=label)
    b.plot(x, [float(r["K_aq_org_total"]) for r in sweep], style, color="#0072B2")
xs = [r["x_acetic_acid_aq"] for r in comparison]
a.plot(xs, [r["pH_esi_model"] for r in comparison], "o", mfc="none", color="#D55E00", label="ePC-SAFT advanced, ESI Tables S11–S12")
a.plot(xs, [r["pH_exp"] for r in comparison], "k^", ms=4, label="Experiment (Ziegenfuß and Maurer)")
b.plot(xs, [r["K_esi_model"] for r in comparison], "o", mfc="none", color="#D55E00")
b.plot(xs, [r["K_exp"] for r in comparison], "k^", ms=4)
a.set(xlabel="Acetic acid mole fraction, aqueous liquid (HAc + Ac⁻)", ylabel="Aqueous pH", xlim=(0, 0.04))
b.set(xlabel="Acetic acid mole fraction, aqueous liquid (HAc + Ac⁻)", ylabel="K = x_aq / x_org (total acid)", xlim=(0, 0.04))
fig.legend(*a.get_legend_handles_labels(), loc="outside lower center", ncols=2, frameon=False, fontsize=8)
(ROOT / "figures").mkdir(exist_ok=True)
fig.savefig(ROOT / "figures/ph-and-partition.png", dpi=160)
