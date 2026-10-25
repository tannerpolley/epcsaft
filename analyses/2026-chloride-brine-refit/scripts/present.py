"""Write the notebook's Markdown tables and figure from retained results only; imports no Engine."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parents[1]
RESULTS = HERE / "results"
RUNS = {"stage-a": "binaries + Rard–Miller", "stage-a-calcium": "+ Saad 1980 (calcium branch)"}
CATIONS = ("sodium-cation", "magnesium-cation", "calcium-cation")
SYMBOL = {"sodium-cation": "Na⁺", "magnesium-cation": "Mg²⁺", "calcium-cation": "Ca²⁺"}
NAME = {"H": "Held 2014 form (H)", "F": "Figiel 2025 form (F)"}
HYDRATED = {"sodium-cation": 4.36, "magnesium-cation": 5.98, "calcium-cation": 5.42}  # Marcus 1991 Table 1, 2(r + Δr) / Å
SOURCE = {"robinson-bower-1966": "Robinson–Bower NaCl–CaCl₂ Φ", "dinane-mounir-2003": "Dinane–Mounir NaCl–MgCl₂ a_w",
          "rard-miller-1987": "Rard–Miller NaCl–MgCl₂ Φ (training)"}
LIMIT = {"robinson-bower-1966": (0.01, 0.03), "dinane-mounir-2003": (0.0013, 0.0040), "rard-miller-1987": (0.01, 0.03)}


def read(path: Path) -> list[dict]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def table(header: list[str], rows: list[list]) -> str:
    return "\n".join(["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
                     + ["| " + " | ".join(str(c) for c in row) + " |" for row in rows]) + "\n"


def number(x, digits: int = 4) -> str:
    return "—" if x is None or x == "" else f"{float(x):.{digits}f}"


def verdict(flag) -> str:
    return "pass" if flag in (True, "True") else "fail"


def rms(values: list[float]) -> float:
    return math.sqrt(sum(v * v for v in values) / len(values))


def coordinate_names(labels: list[str]) -> list[str]:
    """Bound flags index coordinates; the packing target shares the σ coordinate."""
    return [label for label in labels if not label.startswith("packing_diameter/")]


def main() -> None:
    selections = {run: json.loads((RESULTS / run / "selection.json").read_text()) for run in RUNS}
    fits = {(run, form): json.loads((RESULTS / run / f"fit-{form}.json").read_text()) for run in RUNS for form in NAME}
    pitzer = {(b["source"], b["role"], b["I_bucket"]): b for b in json.loads((RESULTS / "pitzer.json").read_text())["buckets"]}
    transfer = {(r["probe"], r["source"], r["I_bucket"]): r for r in read(HERE / "inputs/binary-only-transfer.csv")}
    final = "stage-a-calcium"
    declaration = json.loads((HERE / "inputs/fit-declaration.json").read_text())
    owner = declaration["held_out"]["owner_qualification"]

    rows = []
    for form, chosen in selections[final]["forms"].items():
        optimum = next(o for o in fits[final, form]["optima"] if o["starts"][0] == chosen["start"])
        values = dict(zip(chosen["labels"], optimum["physical"]))
        errors = dict(zip(chosen["labels"], optimum["standard_error"] or [None] * len(chosen["labels"])))
        for ion in CATIONS:
            cells = [("—" if key not in values else number(values[key], 3) + ("" if errors[key] is None else f" ± {errors[key]:.2g}"))
                     for key in (f"segment_diameter/{ion}", f"debye_huckel_diameter/{ion}", f"born_diameter/{ion}",
                                 "k_ij/" + "+".join(sorted((ion, "water"))), "k_ij/" + "+".join(sorted((ion, "chloride-anion"))))]
            rows.append([NAME[form], SYMBOL[ion], *cells, number(HYDRATED[ion], 2)])
    (RESULTS / "parameters-table.md").write_text(table(
        ["Form", "Ion", "σ / Å", "d_DH / Å", "d_Born / Å", "k(ion–water)", "k(ion–Cl⁻)", "Marcus hydrated 2(r+Δr) / Å"], rows))

    rows, starts = [], []
    for (run, form), record in fits.items():
        names = coordinate_names(record["starts"][0]["labels"])
        for k, optimum in enumerate(record["optima"], 1):
            rows.append([RUNS[run], NAME[form], k, len(optimum["starts"]), f"{optimum['final_cost']:.5g}",
                         "yes" if optimum["within_band"] else "no", optimum["covariance_status"],
                         ", ".join(names[i] for i in optimum["physical_bound_coordinates"]) or "none",
                         ", ".join(names[i] for i in optimum["numerical_box_coordinates"]) or "none"])
        statuses = [s["status"] for s in record["starts"]]
        starts.append([RUNS[run], NAME[form], len(statuses), statuses.count("converged"), statuses.count("limit_reached"),
                       sum(s["native_failures"] for s in record["starts"]), f"{record['wall_seconds']:.0f}"])
    (RESULTS / "optima-table.md").write_text(table(
        ["Training", "Form", "Optimum", "Starts", "Training cost", "Within 10 %", "Covariance", "At a physical bound",
         "At the numerical box"], rows))
    (RESULTS / "starts-table.md").write_text(table(
        ["Training", "Form", "Starts", "Converged", "Iteration limit", "Native row failures", "Wall / s"], starts))

    rows = []
    for run, selection in selections.items():
        for form, chosen in selection["forms"].items():
            optimum = next(o for o in fits[run, form]["optima"] if o["starts"][0] == chosen["start"])
            for key, value in optimum["training_rms"].items():
                rows.append([RUNS[run], NAME[form], key, number(value), number(optimum["training_max"][key])])
    (RESULTS / "training-table.md").write_text(table(["Training", "Form", "Source / system", "Φ RMS", "Φ maximum"], rows))

    rows, qualification = [], []
    limits = read(RESULTS / "source-agreement-limits.csv")
    realizations = len({r["realization_id"] for r in limits})
    for run, selection in selections.items():
        for form, chosen in selection["forms"].items():
            for optimum in chosen["optima"]:
                folder = RESULTS / run / "held-out" / optimum["label"]
                unfitted = {(b["source"], b["I_bucket"]): b for b in read(folder / "nominal-buckets.csv") if b["declaration"] == "unfitted"}
                for v in read(folder / "verdicts.csv"):
                    if v["declaration"] != "candidate":
                        continue
                    key = (v["source"], v["I_bucket"])
                    probe = transfer.get(("held088_dh" if form == "H" else "figiel_dh", *key))
                    rows.append([RUNS[run], optimum["label"], SOURCE[v["source"]], v["I_bucket"], number(v["RMS"]),
                                 number(v["maximum_absolute_error"]), number(v["limit_RMS"]) + "/" + number(v["limit_maximum"]),
                                 f"{v['passing_realizations']}/{v['frozen_realizations']}", verdict(v["robust_passed"]),
                                 verdict(v["frozen_robust_passed"]), number(unfitted[key]["RMS"]), number(probe["RMS"]) if probe else "—",
                                 number(pitzer[v["source"], "mixture", v["I_bucket"]]["RMS"])])
                qualification.append([RUNS[run], optimum["label"], f"{optimum['realizations_passing_every_bucket']}/{realizations}",
                                      "yes" if optimum["qualified"] else "no",
                                      "yes" if all(b["frozen_robust_passed"] for b in optimum["held_out"]) else "no",
                                      "yes" if (run, form, optimum["label"]) == (owner["run"], owner["form"], chosen["label"]) else "no"])
    failures = [[RUNS[run], optimum["label"], r["I_bucket"], r["realization_id"], number(r["RMS"], 6), number(r["limit_RMS"], 6),
                 number(r["maximum_absolute_error"], 6), number(r["limit_maximum"], 6)]
                for run, selection in selections.items() for chosen in selection["forms"].values() for optimum in chosen["optima"]
                if (RESULTS / run / "held-out" / optimum["label"] / "dinane-mounir-failures.csv").exists()
                for r in read(RESULTS / run / "held-out" / optimum["label"] / "dinane-mounir-failures.csv")]
    (RESULTS / "dinane-mounir-failures-table.md").write_text(table(
        ["Training", "Optimum", "I / (mol/kg)", "Realization", "RMS", "Amended RMS limit", "Maximum", "Amended maximum limit"],
        failures))
    (RESULTS / "held-out-table.md").write_text(table(
        ["Training", "Optimum", "Source", "I / (mol/kg)", "RMS", "Maximum", "Limit RMS/max (nominal)", "Passing realizations",
         "Verdict", "Frozen-limit verdict", "Unfitted /2 RMS", "Binary-only RMS", "Pitzer RMS"], rows))
    (RESULTS / "qualification-table.md").write_text(table(
        ["Training", "Optimum", "Realizations passing every held-out bucket", "Qualified (amended)", "Qualified (frozen limits)",
         f"Qualified (owner, {owner['date']})"], qualification))

    chosen = selections[owner["run"]]["forms"][owner["form"]]
    (miss,) = read(RESULTS / owner["run"] / "held-out" / chosen["label"] / "dinane-mounir-failures.csv")
    training = [float(r["m_NaCl"]) + float(r["m_KCl"]) + 3 * float(r["m_divalent"])
                for r in read(RESULTS / owner["run"] / f"training-residuals-{owner['form']}.csv")]
    held = [float(r["I_mol_kg"]) for r in read(RESULTS / owner["run"] / f"nominal-residuals-{owner['form']}.csv")
            if r["realization_id"] == "r0000" and r["role"] == "mixture" and r["source"] in LIMIT and r["source"] != "rard-miller-1987"]
    (RESULTS / "owner-qualification-table.md").write_text(table(["Quantity", "Value"], [
        ["Packet", f"`{declaration['forms'][owner['form']]['candidate_packet']}/1` ({NAME[owner['form']]}, {chosen['label']})"],
        ["Systems", "aqueous NaCl–MgCl₂–CaCl₂ chloride brines"],
        ["Temperature / K", number(declaration["conditions"]["temperature_K"], 2)],
        ["Training I / (mol/kg)", number(min(training), 3) + "–" + number(max(training), 2)],
        ["Held-out I / (mol/kg)", number(min(held), 2) + "–" + number(max(held), 2)],
        ["Realizations passing every held-out bucket", f"{chosen['realizations_passing_every_bucket']:,} of {realizations:,}"],
        ["Failing realization", f"`{miss['realization_id']}`, every Dinane–Mounir value lowered by half a printed unit"],
        ["Failing bucket", f"Dinane–Mounir NaCl–MgCl₂ a_w, {miss['I_bucket']}"],
        ["RMS / amended RMS limit", number(miss["RMS"], 6) + " / " + number(miss["limit_RMS"], 6)],
        ["RMS above its limit / %", number(100 * (float(miss["RMS"]) / float(miss["limit_RMS"]) - 1), 1)],
        ["Maximum / maximum limit", number(miss["maximum_absolute_error"]) + " / " + number(miss["limit_maximum"])]]))

    rows = []
    for form, chosen in selections[final]["forms"].items():
        endpoints = [float(r["residual"]) for r in read(RESULTS / final / f"nominal-residuals-{form}.csv") if r["role"] == "CaCl2 endpoint"]
        rows.append([NAME[form], "Robinson–Bower CaCl₂ endpoints", len(endpoints), number(rms(endpoints)), number(max(map(abs, endpoints)))])
        checks = read(RESULTS / final / "held-out" / chosen["label"] / "binary-reference-checks.csv")
        for role, label in (("NaCl anchor", "Robinson–Bower NaCl anchors"), ("MgCl2 binary", "Rard–Miller Table III MgCl₂")):
            values = [float(r["residual"]) for r in checks if r["declaration"] == "candidate" and r["role"] == role]
            rows.append([NAME[form], label, len(values), number(rms(values)), number(max(map(abs, values)))])
    (RESULTS / "reference-checks-table.md").write_text(table(["Form", "Check", "Rows", "Φ RMS", "Φ maximum"], rows))

    consistency = read(RESULTS / "source-consistency.csv")
    rows = []
    for bucket in ("I<=3", "I>3"):
        d = [float(r["difference"]) for r in consistency if r["I_bucket"] == bucket]
        rows.append([bucket, len(d), number(sum(d) / len(d)), number(rms(d)), number(max(map(abs, d))), "0.0013/0.0040"])
    (RESULTS / "consistency-table.md").write_text(table(
        ["I / (mol/kg)", "Rows", "Bias", "RMS", "Maximum", "#556 limit RMS/max"], rows))

    rows = []
    for bucket in ("I<=3", "I>3"):
        b = [r for r in limits if r["I_bucket"] == bucket]
        nominal = next(r for r in b if r["realization_id"] == "r0000")
        rows.append([bucket, number(nominal["limit_RMS"]) + "/" + number(nominal["limit_maximum"])]
                    + [number(min(float(r[k]) for r in b)) + "–" + number(max(float(r[k]) for r in b)) for k in ("limit_RMS", "limit_maximum")]
                    + [sum(float(r["disagreement_RMS"]) > LIMIT["dinane-mounir-2003"][0] for r in b), realizations])
    (RESULTS / "amended-limits-table.md").write_text(table(
        ["I / (mol/kg)", "Nominal limit RMS/max", "RMS limit range", "Maximum limit range",
         "Realizations where the disagreement sets the RMS limit", "Realizations"], rows))

    # Claim-led figure: held-out residuals against ionic strength, discrete measured states as points.
    series = [("unfitted held-2014-strategy2/2", read(RESULTS / final / "nominal-residuals-unfitted.csv"), "o", "#9e9e9e"),
              ("H, calcium branch", read(RESULTS / final / "nominal-residuals-H.csv"), "s", "#1f77b4"),
              ("F, calcium branch", read(RESULTS / final / "nominal-residuals-F.csv"), "^", "#d62728"),
              ("PHREEQC Pitzer (reference model)", read(RESULTS / "pitzer-rows.csv"), "x", "#2ca02c")]
    figure, axes = plt.subplots(1, 2, figsize=(10, 4.8), layout="constrained")
    for axis, source, ylabel in ((axes[0], "robinson-bower-1966", "Φ model − Φ measured"),
                                 (axes[1], "dinane-mounir-2003", "a_w model − a_w measured")):
        limit = LIMIT[source][1]
        axis.axhspan(-limit, limit, color="#eeeeee", zorder=0, label="maximum-error limit")
        axis.axhline(0, color="black", linewidth=0.6)
        for label, data, marker, color in series:
            pts = [(float(r.get("I_mol_kg") or float(r["m_NaCl"]) + 3 * float(r["m_divalent"])), float(r["residual"]))
                   for r in data if r["source"] == source and r["role"] == "mixture" and r["residual"] not in ("", "None")]
            axis.scatter(*zip(*pts), marker=marker, s=18, color=color, label=label, linewidths=1)
        if source == "dinane-mounir-2003":
            axis.scatter([float(r["I_mol_kg"]) for r in consistency], [-float(r["difference"]) for r in consistency], marker="D",
                         s=22, facecolors="none", edgecolors="#9467bd", label="Rard–Miller Φ converted to a_w (source, no model)")
        axis.axvline(3, color="#777777", linestyle=":", linewidth=0.8)
        axis.set_xlabel("ionic strength I / (mol/kg)")
        axis.set_ylabel(ylabel)
        axis.set_title(SOURCE[source] + " (held out)")
    handles = {label: handle for axis in axes for handle, label in zip(*axis.get_legend_handles_labels())}
    figure.legend(handles.values(), handles.keys(), fontsize=7, loc="outside lower center", ncol=3)
    (HERE / "figures").mkdir(exist_ok=True)
    figure.savefig(HERE / "figures/held-out-residuals.png", dpi=150)


if __name__ == "__main__":
    main()
