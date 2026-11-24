"""Compare independent Full and finite-grid searches on the pinned reactive parents."""
from collections import Counter
from dataclasses import replace
import json
from pathlib import Path
from time import perf_counter

import epcsaft
from epcsaft import equilibrium as eq

HERE = Path(__file__).resolve().parent


def load(name):
    return json.loads((HERE / name).read_text())


def trial_value(trial):
    q = trial.quantities
    return {
        "tpd_over_rt": q["tpd_over_rt"],
        "tpd_error_estimate": q.get("tpd_error_estimate_over_rt", q.get("descent_error_estimate_over_rt")),
    }


def search_record(report, elapsed):
    evaluated = [trial for trial in report.trials if "tpd_over_rt" in trial.quantities]
    descending = [trial for trial in evaluated if trial.resolved_descent]
    witness = min(descending, key=lambda t: t.quantities["split_delta_g_over_rt"], default=None)
    return {
        "completed": report.completed,
        "message": report.message,
        "seconds": elapsed,
        "work": dict(report.work, split_attempts=sum(len(t.split_checks) for t in report.trials)),
        "first_tpd": trial_value(evaluated[0]) if evaluated else None,
        "minimum_tpd": trial_value(min(evaluated, key=lambda t: t.quantities["tpd_over_rt"])) if evaluated else None,
        "refusals": Counter(t.message or t.diagnostic_cause or t.evidence
                            for t in report.trials if "tpd_over_rt" not in t.quantities),
        "witness": None if witness is None else {
            "delta_g_over_rt_mol": witness.quantities["split_delta_g_over_rt"],
            "estimated_error_mol": witness.quantities["split_error_estimate_over_rt"],
            "composition": witness.composition,
        },
    }


def replay(case):
    data = load(f"{case}-parent.json")
    model = epcsaft.Mixture(epcsaft.Parameters.from_mapping(data["parameters"]))
    reactions = [
        eq.Reaction(
            row["stoichiometry"], name=row["name"], correlation=eq.ReactionLogPolynomial(**row["correlation"]),
            reference=eq.ReactionReference(
                getattr(eq.ReferenceSourceBasis, row["reference"]["source_basis"]),
                row["reference"].get("reference_pressure_pa"),
            ),
        )
        for row in data["reactions"]
    ]
    ids = model.component_ids
    problem = eq.Problem(
        T=data["T"], P=data["P"], feed=eq.Amounts(data["feed"]), row_tolerance=data["row_tolerance"],
        phases=[eq.Phase("L", kind="liquid", support=ids)], reactions=reactions,
        neutral_reference=eq.NeutralReference(**data["neutral_reference"]),
    )
    compiled = eq.compile_problem(model, problem)
    parent = compiled.solve()
    assert parent.requested_tolerance_met, parent.message
    liquid, = problem._native(model).phases
    neutral = [i for i, c in enumerate(data["parameters"]["components"])
               if c["fixed"]["charge_number"]["value"]["magnitude"] == 0]
    vapor = type(liquid)()
    vapor.name, vapor.kind, vapor.support, vapor.support_declared = "V", type(liquid.kind).Vapor, neutral, True

    start = perf_counter()
    full = search_record(eq.stability_test(compiled, parent), perf_counter() - start)
    result = dict(T=parent.temperature, P=parent.pressure, parent_x=parent.phase_mole_fractions(0),
                  parent_density=parent.molar_densities[0], full=full, grid={})
    for m in (2, 4):
        grid = eq.stability_grid(compiled, parent, [liquid, vapor], m)
        raw = search_record(grid.raw, grid.elapsed_seconds)
        raw["raw_attempts"] = grid.raw_attempts
        raw["coverage"] = [
            dict(side=row.phase.kind.name, support=[ids[i] for i in row.phase.support],
                 planned=row.planned, visited=row.visited, admitted=row.admitted, h=row.h,
                 rounding_bound=row.rounding_bound, mapping_allowance=row.mapping_allowance,
                 omitted_species=[ids[i] for i in row.omitted_species], cause=row.cause,
                 uninformative=row.uninformative,
                 refusals=Counter(grid.raw.trials[i].message for i in row.nodes
                                  if "tpd_over_rt" not in grid.raw.trials[i].quantities))
            for row in grid.coverage
        ]
        available = {k: r["witness"] for k, r in (("full", full), ("grid", raw)) if r["witness"] is not None}
        raw["union_witness_source"] = min(available, key=lambda k: available[k]["delta_g_over_rt_mol"], default=None)
        assert raw["union_witness_source"] is not None, (case, m, full, raw)
        result["grid"][m] = raw

    # Reference compositions are read only after both blind searches.
    split = load(f"{case}-split.json")
    phases = [eq.Phase(f"L{i}", kind="liquid", support=ids,
                       composition_guess=dict(zip(ids, x, strict=True)), packing_guess=rho * model._model.packing_fraction_coefficient(data["T"], x))
              for i, (x, rho) in enumerate(zip(split["x"], split["rho"], strict=True))]
    state = eq.solve_equilibrium(model, replace(problem, phases=phases))
    assert state.requested_tolerance_met, state.message
    assert all(abs(a - b) <= data["row_tolerance"] * abs(a)
               for a, b in zip(split["rho"], state.molar_densities, strict=True)), state.molar_densities
    result["two_liquid"] = dict(message=state.message, criterion_met=state.requested_tolerance_met,
                                amounts=state.amounts, densities=state.molar_densities,
                                parent_objective=parent.objective, split_objective=state.objective,
                                reference_delta_g_over_rt_mol=split["delta_g_over_rt"], reference_error=split.get("error"))
    return result


if __name__ == "__main__":
    for name in ("mdea", "amp"):
        print(json.dumps({"case": name, **replay(name)}), flush=True)
