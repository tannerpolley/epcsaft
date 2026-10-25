"""Bound selected source-coordinate errors with native equilibrium evaluations."""
import os
from concurrent.futures import ProcessPoolExecutor
from time import perf_counter
import json

from calculate import CASE, RESULTS, calculate, figure_tasks, rows, write


def main():
    start = perf_counter()
    tasks = figure_tasks()
    reviewed = rows(RESULTS / "reviewed-state-comparisons.csv")
    final = rows(RESULTS / "state-comparisons.csv")
    def key(row):
        names = "/".join(row["names"]) if isinstance(row["names"], list) else row["names"]
        role = "experiment_paper" if ("experimental" in row["source"] or "markers" in row["source"]) else row["dataset"]
        return (row["figure"], row["series"], row["mode"], row["quantity"], float(row["x"]),
                float(row["T"]) if row["figure"] == "figure-01" else 0., names, role)
    current_index = {key(t): i for i, t in enumerate(tasks)}
    old_by_current = {current_index[key(r)]: r for r in reviewed if key(r) in current_index}
    selected = [i for i, t in enumerate(tasks) if i in old_by_current
                if (t["figure"] == "figure-04" and t["dataset"] == "paper_model" and (t["mode"] == "dew" or t["index"] == 1))
                or (t["figure"] == "figure-03" and t["dataset"] == "paper_model" and t["mode"] == "bubble" and t["index"] in (1, 6))]
    for fig, mode in (("figure-02", "bubble"), ("figure-02", "dew"), ("figure-06", "bubble"), ("figure-07", "bubble")):
        candidates = [i for i, t in enumerate(tasks) if i in old_by_current and t["figure"] == fig and t["mode"] == mode
                      and t["dataset"] == "paper_model" and old_by_current[i]["success"] == "True"]
        selected.append(max(candidates, key=lambda i: abs(float(old_by_current[i]["deviation"]))/tasks[i]["tolerance"]))
    bounds_tasks = []
    for i in selected:
        t = tasks[i]
        # Pure endpoints have exact species identity, not an uncertain binary composition.
        uncertainty = 0. if t["x"] in (0., 1.) else t["input_uncertainty"]
        bounds_tasks.extend({**t, "x": max(0., min(1., t["x"] + sign*uncertainty))} for sign in (-1, 1))
    with ProcessPoolExecutor(max_workers=os.cpu_count()) as pool:
        evaluated = list(pool.map(calculate, bounds_tasks))
    default_statuses = [r["message"] for r in evaluated]
    seeded_calls = 0
    for q, result in enumerate(evaluated):
        if not result["success"]:
            t, center = bounds_tasks[q], final[selected[q//2]]
            fractions, n = json.loads(center["compositions"]), len(t["names"])
            seed = dict(compositions=[fractions[k:k+n] for k in range(0, len(fractions), n)],
                        packing_fractions=json.loads(center["packing_fractions"]), amounts=[1., 0.],
                        pressure=float(center["P_solved_Pa"]))
            evaluated[q] = calculate(t, seed=seed)
            seeded_calls += 1
    bounds = []
    for q, i in enumerate(selected):
        t = tasks[i]
        low, high = evaluated[2*q:2*q+2]
        assert low["success"] and high["success"], (i, low["message"], high["message"])
        values = [low["value"], float(final[i]["value"]), high["value"]]
        lower, upper = min(values), max(values)
        old_reference = float(old_by_current[i]["reference"])
        bounds.append(dict(csv_line=i+2, figure=t["figure"], series=t["series"], mode=t["mode"], quantity=t["quantity"],
                           x=t["x"], input_uncertainty=t["input_uncertainty"], x_low=low["x"], x_high=high["x"],
                           value_low_x=low["value"], value_center=float(final[i]["value"]), value_high_x=high["value"],
                           output_lower=lower, output_upper=upper, tolerance=t["tolerance"], old_reference=old_reference,
                           corrected_reference=t["reference"], old_reference_interval_overlap=lower-t["tolerance"] <= old_reference <= upper+t["tolerance"],
                           maximum_residual=max(low["residual"], high["residual"]),
                           low_default_message=default_statuses[2*q], high_default_message=default_statuses[2*q+1],
                           method="Native solves at composition interval endpoints and center; sampled envelope, not a rigorous interval bound"))
    write(RESULTS / "input-coordinate-bounds.csv", bounds)
    by_line = {b["csv_line"]: b for b in bounds}
    corrections = {}
    for fig in ("figure-01", "figure-02"):
        for correction in rows(CASE / f"figures/{fig}/source/curve-corrections.csv"):
            if correction["disposition"] == "re-extracted":
                corrections[(fig, correction["series"], float(correction["input_coordinate"]))] = correction
    diagnosed = []
    for i, old in enumerate(reviewed):
        if old["within_curve_tolerance"] != "False":
            continue
        current = current_index.get(key(old))
        if current is None:
            diagnosed.append(dict(reviewed_csv_line=i+2, figure=old["figure"], series=old["series"], mode=old["mode"], quantity=old["quantity"],
                                  x=old["x"], old_reference=old["reference"], corrected_reference="", engine_value=old["value"],
                                  old_deviation=old["deviation"], corrected_deviation="", tolerance=old["tolerance"],
                                  classification="obscured model curve", evidence="Experimental symbol hides the solid curve; model row omitted and experimental observation retained. See figure-local curve-corrections.csv and QA overlay."))
            continue
        new, line = final[current], current+2
        classification, evidence = "undiagnosed", "Not included in the bounded source/coordinate diagnosis."
        if new["dataset"] == "experiment_paper":
            classification, evidence = "extraction error", "Source labels identify an experimental marker; it had been counted as a model curve. Reclassified without changing its coordinates."
        elif (old["figure"] == "figure-04" and int(old["index"]) < 11) or (old["figure"] == "figure-03" and int(old["index"]) in (1, 6)):
            classification, evidence = "extraction error", "Original coordinates fail the PDF solid-curve check. Native-raster re-extraction, unchanged composition mapping and unchanged tolerance recover agreement. See figure-local curve-corrections.csv and QA overlay."
        elif (old["figure"], old["series"], float(old["x"])) in corrections:
            classification, evidence = "extraction error", "Marker-overlap proxy re-extracted from the readable solid PC-SAFT stroke; experimental marker retained separately. See figure-local curve-corrections.csv and QA overlay."
        elif line in by_line and by_line[line]["old_reference_interval_overlap"]:
            classification, evidence = "input-coordinate uncertainty", "Reference interval overlaps sampled native output envelope at the retained composition uncertainty. This is consistency with input uncertainty, not proof of extraction accuracy."
        elif line in by_line:
            evidence = "Input-coordinate envelope does not explain the miss; solid-curve extraction has not been rechecked for this row. Calculation disagreement versus extraction remains unresolved."
        diagnosed.append(dict(reviewed_csv_line=i+2, figure=old["figure"], series=old["series"], mode=old["mode"], quantity=old["quantity"],
                              x=old["x"], old_reference=old["reference"], corrected_reference=new["reference"], engine_value=new["value"],
                              old_deviation=old["deviation"], corrected_deviation=new["deviation"], tolerance=old["tolerance"],
                              classification=classification, evidence=evidence))
    write(RESULTS / "model-miss-diagnosis.csv", diagnosed)
    metadata = dict(wall_seconds=perf_counter()-start, workers=os.cpu_count(), selected_rows=len(selected), native_calls=len(evaluated)+seeded_calls, seeded_calls=seeded_calls,
                    input_bounds="Only selected coordinates; no tolerance widened", convention_findings="No convention difference established")
    (RESULTS / "coordinate-diagnosis.json").write_text(json.dumps(metadata, indent=2)+"\n")
    print(json.dumps(metadata))


if __name__ == "__main__":
    main()
