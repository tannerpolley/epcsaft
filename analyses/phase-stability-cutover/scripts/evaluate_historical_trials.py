"""Evaluate retained historical TPD argmins through native supplied trials."""
import argparse
import csv
import importlib.metadata
import json
from pathlib import Path

import epcsaft
from epcsaft import equilibrium

from current_engine import HERE, ROOT, METHANE_ETHANE, TERNARY, CHARGED, KHUDAIDA, KHUDAIDA_C, sha, write_rows

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--output-dir", type=Path, required=True)
output = parser.parse_args().output_dir
output.mkdir(parents=True, exist_ok=False)
read = lambda path: list(csv.DictReader(path.open()))
history_path = HERE / "results/current-engine.csv"
native_path = HERE / "results/native-e8ad39e0/current-engine.csv"
historical = read(history_path)
native = read(native_path)
comparisons, candidates, unavailable, quantities = [], [], [], []
criterion = 1e-7
refusal = "permittivity model has no low-density limit; declare phases"
for old, previous in zip(historical, native, strict=True):
    case, label = old["case_id"], old["declared_topology"]
    assert (case, label) == (previous["case_id"], previous["declared_topology"])
    parameter_path, components = {
        "NEUTRAL_SINGLETON": (METHANE_ETHANE, ("methane", "ethane")),
        "CHARGED_SINGLETON": (CHARGED, None),
        "NEUTRAL_TERNARY_IMMISCIBLE": (TERNARY, None),
        "KHUDAIDA_ROW_7": (KHUDAIDA, KHUDAIDA_C),
    }[case]
    model = epcsaft.Mixture(epcsaft.Parameters.from_json(parameter_path, components=components))
    ids = model.component_ids
    guesses = json.loads(previous["phase_mole_fractions"])
    full_precision = label.startswith("1 ")
    if full_precision:
        # Source-defined binary64 points; rounded historical output is not used.
        guesses = [[0.5, 0.5]] if case == "NEUTRAL_SINGLETON" else (
            [[0.98, 0.01, 0.01]] if case == "CHARGED_SINGLETON" else [[1/3, 1/3, 1/3]])
    kind = "vapor" if label == "1 vapor" else "liquid"
    problem = equilibrium.Problem(T=float(previous["T_K"]), P=float(previous["P_Pa"]),
        feed=equilibrium.Amounts(dict(zip(ids, json.loads(previous["feed_mol"]), strict=True))),
        phases=[equilibrium.Phase(f"phase-{p}", kind=kind, support=ids, composition_guess=x)
                for p, x in enumerate(guesses)], row_tolerance=1e-10)
    central = equilibrium.solve_equilibrium(model, problem)
    if not central.success or not central.requested_tolerance_met:
        raise RuntimeError(f"{case}/{label}: {central.message}")
    compiled = equilibrium.compile_problem(model, problem)
    options = equilibrium.StabilityOptions()
    options.evaluate_supplied_trials = True
    declared_trials = []
    for branch in ("liquid", "vapor"):
        phase = problem._native(model).phases[0]
        phase.kind = getattr(epcsaft._core.PhaseKind, branch.capitalize())
        phase.composition_guess = json.loads(old[f"tpd_argmin_{branch}"])
        # Charged coordinates are historical formula-salt coordinates. Admission
        # refuses the whole charged family before these points can be evaluated.
        phase.packing_guess = None
        declared_trials.append(phase)
    options.trials = declared_trials
    report = None
    try:
        report = equilibrium.stability_test(compiled, central, options)
    except epcsaft.ParameterError as error:
        if str(error) != refusal:
            raise
    key = f"{case}/{label}"
    candidates.append(dict(case=key, T_K=central.temperature, P_Pa=central.pressure,
        input_compositions=json.dumps(guesses), source_defined_full_precision=full_precision,
        candidate_compositions=json.dumps(central.mole_fractions),
        densities_mol_m3=json.dumps(central.molar_densities),
        coordinates=json.dumps(central.coordinates), multipliers=json.dumps(central.multipliers),
        amounts_mol=json.dumps(central.amounts), objective_over_rt_mol=central.objective,
        classified_residual_max=max(map(abs, central.classified_residuals)),
        plane_over_rt="" if report is None else json.dumps(report.plane_over_rt),
        standard_over_rt="" if report is None else json.dumps(report.standard_over_rt)))
    unavailable.append(dict(case=key, comparison="historical candidate snapshot", phase="all",
        reason="No historical density or chemical plane retained; matching c74 wheel not found. Composition alone is not a thermodynamic snapshot."))
    for i, branch in enumerate(("liquid", "vapor")):
        before = float(old[f"tpd_min_{branch}"])
        if report is None:
            unavailable.append(dict(case=key, comparison="historical charged trial", phase=branch, reason=refusal))
            continue
        trial = report.trials[i]
        q = trial.quantities
        if "tpd_over_rt" not in q:
            raise RuntimeError(f"{key}/{branch}: {trial.message}")
        difference = q["tpd_over_rt"] - before
        comparisons.append(dict(case=key, phase=branch, historical_tpd_over_rt=before,
            native_tpd_over_rt=q["tpd_over_rt"], difference=difference, criterion=criterion,
            passes=abs(difference) <= criterion,
            supplied_composition=json.dumps(options.trials[i].composition_guess),
            evaluated_composition=json.dumps(trial.composition),
            displacement_l1=q["input_reconstruction_l1"],
            rounding_linf_bound=0.5e-6, evidence=trial.evidence, message=trial.message))
        quantities.extend(dict(case=key, phase=branch, quantity=k, value=v) for k, v in q.items())
        print(json.dumps(comparisons[-1]), flush=True)
for name, rows in (("comparison", comparisons), ("candidates", candidates),
                   ("unverifiable", unavailable), ("quantities", quantities)):
    write_rows(output / f"{name}.csv", rows)
wheel = Path(json.loads(importlib.metadata.distribution("epcsaft").read_text("direct_url.json"))["url"].removeprefix("file://"))
(output / "summary.json").write_text(json.dumps(dict(
    wheel_sha256=sha(wheel), script_sha256=sha(Path(__file__)), criterion=criterion,
    inputs={str(p.relative_to(ROOT)): sha(p) for p in (history_path, native_path, METHANE_ETHANE, TERNARY, CHARGED, KHUDAIDA)},
    outputs={p.name: sha(p) for p in output.glob("*.csv")},
    method="Native supplied-composition evaluation without packing guesses; branch-selected density roots and newly solved candidate planes. One-phase guesses use historical source definitions; multiphase guesses use retained native candidates.",
    basis="Historical TPD values retain full precision; historical argmins retain six decimals. Stationarity at a minimum makes composition rounding affect TPD at second order, provided the same local branch and plane apply. Unknown exact historical argmins preclude an exact displacement or rigorous rounding-error bound.",
    limits=["not global certificates", "numerical verification, not physical validation",
            "candidate planes are newly solved; historical snapshot equivalence unverifiable",
            "charged whole-family refusal under density-independent permittivity"],
    stage="Stopped before independent delivered review and PR"), indent=2) + "\n")
