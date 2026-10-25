from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt

CAMPAIGN = Path(__file__).resolve().parents[1]
REPO = CAMPAIGN.parents[1]
COMPONENTS = ("water", "ethanol", "isobutanol", "nacl")
PHASES = ("aqueous", "organic")


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def records(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def ternary(values) -> tuple[float, float]:
    water, ethanol, isobutanol = values
    total = water + ethanol + isobutanol
    return ((ethanol + 0.5 * isobutanol) / total, math.sqrt(3) * isobutanol / (2 * total))


def qualify(results: Path, output: Path) -> None:
    if output.exists() and any(output.iterdir()):
        raise RuntimeError("output directory must be new or empty")
    output.mkdir(parents=True, exist_ok=True)
    refs = json.loads((CAMPAIGN / "inputs/references.json").read_text())
    paths = {name: REPO / item["path"] for name, item in refs.items()}
    if any(sha(paths[name]) != item["sha256"] for name, item in refs.items()):
        raise RuntimeError("qualification input hash mismatch")
    stage_path, stage_summary_path = results / "stage-results.csv", results / "summary.json"
    stage, stage_summary = records(stage_path), json.loads(stage_summary_path.read_text())
    keys = [(int(row["point_id"]), row["phase"]) for row in stage]
    expected_keys = [(point, phase) for point in range(1, 9) for phase in PHASES]
    if sorted(keys) != sorted(expected_keys) or stage_summary["status"] != "accepted" or stage_summary["failures"]:
        raise RuntimeError("current stage evidence is incomplete or unaccepted")
    limits = {"engine_material_balance_max_abs": 1e-10, "engine_phase_charge_max_abs": 1e-10,
              "engine_pressure_max_relative": 1e-8, "engine_classified_residual_max_abs": 1e-10,
              "idaes_component_balance_max_abs": 1e-10, "idaes_charge_residual_abs": 1e-10,
              "idaes_phase_fraction_abs": 1e-12, "idaes_phase_flow_closure_abs": 1e-10,
              "idaes_composition_max_abs": 1e-12}
    for row in stage:
        values = {name: float(row[name]) for name in limits}
        slope, hessian = float(row["pressure_density_derivative"]), float(row["engine_reduced_hessian_minimum"])
        composition = [float(row[f"x_{name}"]) for name in ("water", "ethanol", "isobutanol", "sodium_cation", "chloride_anion")]
        if (not all(math.isfinite(value) for value in (*values.values(), slope, hessian, *composition))
                or any(value < 0 or value > limits[name] for name, value in values.items()) or slope <= 0 or hessian <= 0
                or min(composition) < 0 or abs(math.fsum(composition) - 1) > 1e-10 or abs(composition[3] - composition[4]) > 1e-10
                or int(row["idaes_dof"]) != 0 or row["idaes_termination"] != "optimal" or row["status"] != "accepted"):
            raise RuntimeError(f"current numerical certificate failed for point {row['point_id']} {row['phase']}")
    mapping_rows = records(paths["feed_to_tie_lines"])
    mapping = {int(row["feed_point_id"]): int(row["experimental_tie_line_id"]) for row in mapping_rows}
    if len(mapping_rows) != 8 or mapping != {point: 9 - point for point in range(1, 9)}:
        raise RuntimeError("feed-to-tie-line mapping is not the required bijection")
    experiment_rows = [row for row in records(paths["experimental_tie_lines"])
                       if row["source_table"] == "Table 3" and float(row["salt_wtfrac"]) == 0.05 and float(row["temperature_K"]) == 293.15]
    if sorted(int(row["tie_line"]) for row in experiment_rows) != list(range(1, 9)):
        raise RuntimeError("expected exactly eight 293.15 K, 5 wt% Table 3 rows")
    experiment = {int(row["tie_line"]): row for row in experiment_rows}
    paper = {int(row["point_id"]): row for row in records(paths["paper_epcsaft_curve"])}
    uncertainties = {row["component"]: float(row["standard_uncertainty_mole_fraction"])
                     for row in records(paths["experimental_uncertainties"])}
    table9 = next(row for row in records(paths["paper_table_9"]) if row["T/K"] == "293.15" and row["model"] == "ePC-SAFT")
    paper_grand_aad = float(table9["grand $\\mathrm{AAD}{ }^{b}$"])
    paper_values = [float(paper[point][f"x_{component}_salt_free"]) for point in paper for component in COMPONENTS[:3]]
    if (not math.isfinite(paper_grand_aad) or paper_grand_aad <= 0 or sorted(paper) != list(range(1, 9))
            or not all(math.isfinite(value) and value >= 0 for value in (*paper_values, *uncertainties.values()))):
        raise RuntimeError("paper ePC-SAFT reference is incomplete or nonfinite")
    comparison = []
    for row in stage:
        point, phase, tie = int(row["point_id"]), row["phase"], mapping[int(row["point_id"])]
        ion = (float(row["x_sodium_cation"]) + float(row["x_chloride_anion"])) / 2
        explicit = [float(row[f"x_{name}"]) for name in COMPONENTS[:3]] + [ion]
        formula = [value / math.fsum(explicit) for value in explicit]
        measured = [float(experiment[tie][f"x_{name}_{'aq' if phase == 'aqueous' else 'org'}"]) for name in COMPONENTS]
        feed = [float(row[f"feed_formula_{name}"]) for name in COMPONENTS[:3]]
        if not all(math.isfinite(value) and value >= 0 for value in (*explicit, *formula, *measured, *feed)):
            raise RuntimeError(f"nonfinite comparison value for point {point} {phase}")
        current_sf, measured_sf = ternary(explicit[:3]), ternary(measured[:3])
        feed_sf = ternary(feed)
        for index, component in enumerate(COMPONENTS):
            comparison.append({"point_id": point, "tie_line_id": tie, "phase": phase, "component": component,
                "current_formula": formula[index], "experimental_formula": measured[index], "absolute_error": abs(formula[index] - measured[index]),
                "current_ternary_x": current_sf[0], "current_ternary_y": current_sf[1], "experimental_ternary_x": measured_sf[0], "experimental_ternary_y": measured_sf[1],
                "feed_ternary_x": feed_sf[0], "feed_ternary_y": feed_sf[1], "experimental_standard_uncertainty": uncertainties[component],
                "paper_epcsaft_salt_free": paper[tie].get(f"x_{component}_salt_free", "") if phase == "organic" and component != "nacl" else ""})
    comparison_path = output / "qualification-comparison.csv"
    with comparison_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=tuple(comparison[0]), lineterminator="\n"); writer.writeheader(); writer.writerows(comparison)
    grand_aad = math.fsum(row["absolute_error"] for row in comparison) / 64
    paper_errors = [abs(row["current_formula"] / math.fsum(item["current_formula"] for item in comparison if item["point_id"] == row["point_id"] and item["phase"] == "organic" and item["component"] != "nacl") - float(row["paper_epcsaft_salt_free"])) for row in comparison if row["paper_epcsaft_salt_free"] != ""]
    if not math.isfinite(grand_aad) or not all(math.isfinite(value) for value in paper_errors) or grand_aad > paper_grand_aad:
        raise RuntimeError(f"current grand AAD {grand_aad} exceeds paper ePC-SAFT {paper_grand_aad}")
    fig, ax = plt.subplots(figsize=(9, 7.5)); height = math.sqrt(3) / 2
    ax.plot([0, 1, .5, 0], [0, 0, height, 0], color="black", linewidth=1.2)
    by_key = {(row["point_id"], row["phase"], row["component"]): row for row in comparison}
    for point in range(1, 9):
        for series, xcol, ycol, color, marker in (("experiment", "experimental_ternary_x", "experimental_ternary_y", "#333333", "o"), ("current Engine", "current_ternary_x", "current_ternary_y", "#087f8c", "s")):
            pair = [by_key[(point, phase, "water")] for phase in PHASES]
            ax.plot([float(row[xcol]) for row in pair], [float(row[ycol]) for row in pair], marker=marker, color=color, linewidth=1, markersize=4, alpha=.8, label=series if point == 1 else None)
    feed_rows = [by_key[(point, "aqueous", "water")] for point in range(1, 9)]
    ax.scatter([float(row["feed_ternary_x"]) for row in feed_rows], [float(row["feed_ternary_y"]) for row in feed_rows], marker="^", color="#d17c0f", label="Figure 2 feeds")
    paper_xy = [ternary([float(paper[t][f"x_{c}_salt_free"]) for c in COMPONENTS[:3]]) for t in range(1, 9)]
    ax.scatter([p[0] for p in paper_xy], [p[1] for p in paper_xy], marker="x", color="#7d3c98", label="paper ePC-SAFT organic")
    ax.text(-.02, -.025, "water", ha="right", va="top"); ax.text(1.02, -.025, "ethanol", ha="left", va="top"); ax.text(.5, height + .025, "isobutanol", ha="center", va="bottom")
    ax.set(xlim=(-.08, 1.08), ylim=(-.07, height + .08), aspect="equal", title=f"Current ePC-SAFT / IDAES reference: AAD {grand_aad:.5f} < paper {paper_grand_aad:.4f}")
    ax.axis("off"); ax.legend(loc="upper right"); fig.tight_layout()
    figure_path = output / "ternary-tie-lines.png"; fig.savefig(figure_path, dpi=200); plt.close(fig)
    summary = {"status": "positive_reference", "grand_formula_aad": grand_aad, "paper_epcsaft_grand_aad": paper_grand_aad,
        "paper_organic_salt_free_aad": math.fsum(paper_errors) / len(paper_errors), "paper_organic_salt_free_max_abs": max(paper_errors),
        "counts": {"stage_rows": 16, "formula_comparisons": 64}, "input_sha256": {**{name: item["sha256"] for name, item in refs.items()}, "stage_results": sha(stage_path), "stage_summary": sha(stage_summary_path)},
        "output_sha256": {comparison_path.name: sha(comparison_path), figure_path.name: sha(figure_path)},
        "claim_limits": ["fixed-parameter paper-correlation-scale agreement", "not agreement within experimental uncertainty", "declared local two-liquid branch", "not reactive Lithium qualification"]}
    (output / "qualification.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(); parser.add_argument("--results-dir", type=Path, required=True); parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args(); qualify(args.results_dir.resolve(), args.output_dir.resolve())


if __name__ == "__main__":
    main()
