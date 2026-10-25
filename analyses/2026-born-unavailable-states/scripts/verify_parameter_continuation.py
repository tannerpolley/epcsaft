"""Compare native paths and verify exact first actions with centered native re-solves."""
import argparse
import csv
import json
from pathlib import Path

import numpy as np
import epcsaft
from epcsaft import equilibrium as eq, _core
from reproduce import HERE, archived, load_case, changed, problem_from_request, saved_state, seed_from_saved, table


def verify_actions(out):
    comparisons, diagnostics = [], []
    for case, value in ((1, 0.5), (2, 3.0)):
        mapping, request, anchor = load_case(case)
        problem = problem_from_request(request, anchor)
        # The retained smooth point supplies guesses; the native driver recomputes it.
        seed_from_saved(problem, case, value)
        family = epcsaft.ParameterFamily.BornCDielectric if case == 1 else epcsaft.ParameterFamily.IonicRegionRelativePermittivity
        active = [epcsaft.ActiveParameter(family, (0,))]
        model = epcsaft.Mixture(epcsaft.Parameters.from_mapping(changed(mapping, case, value)), active_parameters=active)
        central = eq.solve_equilibrium(model, problem)
        assert central.success and central.requested_tolerance_met, central.message
        compiled = eq.compile_problem(model, problem)
        n = len(model.component_ids)
        observables = [_core.SolvedStateObservable(_core.SolvedStateObservableKind.Pressure)]
        for phase in range(central.phase_count):
            observables.append(_core.SolvedStateObservable(_core.SolvedStateObservableKind.PhaseMolarDensity, phase))
            observables.extend(_core.SolvedStateObservable(_core.SolvedStateObservableKind.PhaseMoleFraction, phase, i) for i in range(n))
        batches = _core.solved_state_first_actions(compiled, central, observables,
            [_core.SolvedStateActionDirection(0.0, 0.0, (0.0,) * n, [1.0])])
        batch = batches[0]
        assert all(row.action is not None for row in batch.results), [row.diagnostic_message for row in batch.results]
        exact = np.array([r.action for r in batch.results])
        diag = batch.diagnostics
        diagnostics.append((case, value, diag.rank, diag.condition_estimate, diag.maximum_linear_residual, diag.central_residual, diag.central_floor))
        tangent = np.array(_core.sensitivity(compiled, central.coordinates, (0.0,) * (2+n) + (1.0,)))
        assert tangent.size == len(central.coordinates)
        # Verification of F_y y_theta + F_theta = 0; native owners supply F_y and y_theta.
        row_action = -np.array(compiled.jacobian(central.coordinates)).reshape(len(central.rows), -1) @ tangent
        steps = (1e-3, 1e-4, 1e-5, 1e-6) if case == 1 else (1e-3, 1e-4, 1e-5)
        for h in steps:
            values, rows = [], []
            for sign in (1.0, -1.0):
                shifted = epcsaft.Mixture(epcsaft.Parameters.from_mapping(changed(mapping, case, value+sign*h)), active_parameters=active)
                result = eq.solve_equilibrium(shifted, problem, central)
                assert result.success and result.requested_tolerance_met, result.message
                sample = [result.pressure]
                for phase in range(result.phase_count):
                    sample.append(result.molar_densities[phase])
                    sample.extend(result.mole_fractions[phase*n:(phase+1)*n])
                values.append(np.array(sample))
                rows.append(np.array(eq.compile_problem(shifted, problem).evaluate(central.coordinates).values))
            difference = (values[0]-values[1])/(2*h)
            row_difference = (rows[0]-rows[1])/(2*h)
            for kind, analytic, finite in (("solved", exact, difference), ("row", row_action, row_difference)):
                for i, (a, f) in enumerate(zip(analytic, finite, strict=True)):
                    comparisons.append((case, value, h, kind, i, a, f, abs(a-f)/max(1.0, abs(a), abs(f))))
    table(out, "parameter-actions.csv", ("case", "value", "h", "kind", "index", "exact", "difference", "scaled_error"), comparisons)
    table(out, "parameter-action-diagnostics.csv", ("case", "value", "rank", "condition_estimate", "maximum_linear_residual", "central_residual", "central_floor"), diagnostics)
    for case in (1, 2):
        errors = [max(r[-1] for r in comparisons if r[0] == case and r[2] == h) for h in sorted({r[2] for r in comparisons if r[0] == case}, reverse=True)]
        print("first actions", case, errors, flush=True)
        assert min(errors) <= 1e-5


def compare_paths(out):
    comparison = []
    for run in sorted((HERE / "results/native-parameter").glob("*/case*-states.csv")):
        case = int(run.name[4])
        states = list(csv.DictReader(run.open()))
        values = sorted({float(r["value"]) for r in states}, reverse=True)
        for value in values:
            old = saved_state(case, value)
            if old is None: continue
            historical = list(csv.DictReader(archived()[old].splitlines()))
            current = [r for r in states if float(r["value"]) == value]
            assert len(current) == len(historical)
            for new, previous in zip(current, historical, strict=True):
                assert (new["phase"], new["species"]) == (previous["phase"], previous["species"])
                for column in ("mole_fraction", "rho_mol_m3", "P_Pa"):
                    a, b = float(new[column]), float(previous[column])
                    error = abs(a-b) if column == "mole_fraction" else abs(a-b)/abs(b)
                    comparison.append((run.parent.name, case, value, new["phase"], new["species"], column, a, b, error))
                    assert error <= 1e-7, (case, value, column, error)
            term_old = archived()[old.replace("-state.csv", "-terms.csv")].splitlines()
            terms = [r for r in csv.DictReader(run.with_name(run.name.replace("-states", "-terms")).open()) if float(r["value"]) == value]
            for new, previous in zip(terms, csv.DictReader(term_old), strict=True):
                for column in ("born_over_RT", "debye_huckel_over_RT", "bulk_relative_permittivity"):
                    a, b = float(new[column]), float(previous[column])
                    error = abs(a-b)/abs(b) if b else abs(a-b)
                    comparison.append((run.parent.name, case, value, new["phase"], "", column, a, b, error))
                    assert error <= (1e-7 if b else 1e-10), (case, value, column, error)
        rows = list(csv.DictReader(run.with_name(run.name.replace("-states", "-rows")).open()))
        for column in ("residual", "classified_residual"):
            assert max(abs(float(r[column])) for r in rows if r["accepted"] == "True") <= 1e-10
    print("common point comparisons", len(comparison), "maximum error", max(r[-1] for r in comparison), flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=HERE / "results/current-main")
    parser.add_argument("--actions-only", action="store_true")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    verify_actions(args.out)
    if not args.actions_only: compare_paths(args.out)
    (args.out / "engine-identity.json").write_text(json.dumps({"engine_provenance": list(epcsaft.provenance())}, indent=2)+"\n")


if __name__ == "__main__":
    main()
