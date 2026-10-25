"""Diagnose the copied reactive bubble requests using the installed Engine only."""

import argparse
import copy
import csv
import functools
import hashlib
import json
import os
import math
import re
import tarfile
import tempfile
from pathlib import Path
import time

import epcsaft
from epcsaft import equilibrium as eq

HERE = Path(__file__).resolve().parents[1]
NAMES = ("case1-shell-only-1-0", "case2-eps-ion-2")


def load_case(case):
    name = NAMES[case - 1]
    manifest = json.loads((HERE / "inputs/manifest.json").read_text())
    # The MEA exporter was provenance only (#341 removed it); its hash stays in the manifest.
    for file, record in manifest.items():
        if file == "export_cases.py":
            continue
        assert (
            hashlib.sha256((HERE / "inputs" / file).read_bytes()).hexdigest()
            == record["sha256"]
        )
    mapping = json.loads((HERE / "inputs" / (name + "-parameters.json")).read_text())
    request = json.loads((HERE / "inputs" / (name + "-request.json")).read_text())[
        "reaction_corrected_request"
    ]
    anchor = json.loads((HERE / "inputs" / (name + "-anchor.json")).read_text())[
        "anchor"
    ]
    # The frozen input has no CO2/MEA cross edge; declared to keep that model (#569), not a physical-absence claim.
    mapping["topology"]["no_cross_association"] = [{
        "component_id_a": "carbon-dioxide", "component_id_b": "monoethanolamine",
        "provenance": {"source_id": name, "locator": f"inputs/{name}-parameters.json: no CO2/MEA edge"}}]
    return mapping, request, anchor


def problem_from_request(request, anchor):
    system = request["reaction_system"]
    species = system["species_ids"]
    states = {p["role"]: p for p in request["continuation"]["state"]["phases"]}
    phases = []
    for p in request["phases"]:
        support = (
            species
            if p["support"]["kind"] == "all_components"
            else p["support"]["component_ids"]
        )
        guess = (
            anchor["mole_fractions"]
            if p["fluid_role"] == "liquid"
            else states[p["fluid_role"]]["mole_fractions"]
        )
        phases.append(
            eq.Phase(
                p["identity"],
                support=support,
                kind=p["fluid_role"],
                amount=eq.Pinned(0.0) if p["amount_role"] == "incipient" else eq.Free(),
                composition_guess=guess,
            )
        )
    reactions = []
    for row, r in zip(
        system["reaction_matrix"], system["engine_reactions"], strict=True
    ):
        c = r["engine_correlation"]
        correlation = eq.ReactionLogPolynomial(
            *(
                c[k]
                for k in (
                    "a",
                    "b",
                    "c",
                    "d",
                    "reference_temperature",
                    "temperature_min",
                    "temperature_max",
                    "standard_state_id",
                )
            )
        )
        ref = r["engine_reference"]
        reactions.append(
            eq.Reaction(
                dict(zip(species, row, strict=True)),
                name=r["reaction_id"],
                correlation=correlation,
                reference=eq.ReactionReference(
                    getattr(eq.ReferenceSourceBasis, ref["source_basis"]),
                    reference_pressure_pa=ref["reference_pressure_pa"],
                ),
            )
        )
    source = system["source_standard_state"]
    solvent = species[
        max(range(len(species)), key=source["solvent_composition"].__getitem__)
    ]
    ratios = {
        s: float(z == 0)
        for s, z in zip(species, system["charges"], strict=True)
        if s != solvent
    }
    for pair, coefficients in source["ion_pairs"]:
        for s, c in zip(pair, coefficients, strict=True):
            ratios[s] += c
    ratios = {s: v * 0.125 for s, v in ratios.items()}
    return eq.Problem(
        phases=phases,
        T=request["temperature"]["value"],
        P=eq.Free(anchor["pressure_pa"]),
        feed=eq.Amounts(dict(zip(species, system["feed_amounts_mol"], strict=True))),
        reactions=reactions,
        standard_gibbs_j_per_mol=system.get("standard_gibbs_j_per_mol", ()),
        topology_declared=True,
        neutral_reference=eq.NeutralReference(
            solvent=solvent,
            solute_molality_ratios=ratios,
            standard_molality_mol_per_kg=1.0,
        ),
    )


def changed(mapping, case, value):
    out = copy.deepcopy(mapping)
    if case == 1:
        next(f for f in out["model_families"] if f["kind"] == "electrolyte")[
            "c_dielectric"
        ] = value
    else:
        next(
            c
            for c in out["model_coefficients"]
            if c["family"] == "ionic_region_relative_permittivity"
        )["value"]["magnitude"] = value
    return out


def table(out, name, header, rows):
    rows = list(rows)
    path = out / name
    if not rows:
        path.unlink(missing_ok=True)
        return
    with path.open("w") as file:
        writer = csv.writer(file, lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def record(out, tag, case, value, result, capture, seconds, log, iteration_span=None):
    out.mkdir(parents=True, exist_ok=True)
    convergence = {}
    for key in ("residual", "floor", "amplification", "contraction", "stop"):
        diagnostic = getattr(result.convergence, key)
        convergence[key] = (
            str(diagnostic)
            if key == "stop"
            else (
                None
                if isinstance(diagnostic, float) and not math.isfinite(diagnostic)
                else diagnostic
            )
        )
    summary = {
        "case": case,
        "value": value,
        "success": result.success,
        "solver_status": result.solver_status,
        "message": result.message,
        "iterations": result.iterations,
        "evaluations": result.evaluations,
        "generic_jacobian_rank": result.jacobian_rank,
        "generic_nullity": result.generic_nullity,
        "reduced_hessian_minimum": result.reduced_hessian_minimum
        if result.success
        else None,
        "criterion_met": result.requested_tolerance_met,
        "wall_s": seconds,
        "convergence": convergence,
        "state_file": tag + "-state.csv" if result.success else None,
    }
    (out / (tag + ".json")).write_text(
        json.dumps(summary, indent=2, allow_nan=False) + "\n"
    )
    table(
        out,
        tag + "-rows.csv",
        (
            "row",
            "kind",
            "phase",
            "index",
            "equation",
            "residual",
            "classified_residual",
        ),
        (
            (i, str(row.kind), row.phase, row.index, row.eqid, raw, classified)
            for i, (row, raw, classified) in enumerate(
                zip(result.rows, result.residuals, result.classified_residuals)
            )
        ),
    )
    table(
        out,
        tag + "-iterations.csv",
        (
            "iteration",
            "inf_pr",
            "inf_du",
            "mu",
            "d_norm",
            "regularization",
            "alpha_du",
            "alpha_pr",
            "ls_trials",
        ),
        (
            (e.iteration, *e.values)
            for e in capture.events
            if e.name == "ipopt_iteration" and (e.span == iteration_span if iteration_span is not None
                                               else capture.spans[e.span].parent is None)
        ),
    )
    if not result.success:
        coordinates = (
            []
            if not result.residuals
            else [
                float(v)
                for v in re.findall(r"final x unscaled\[\s*\d+\]=\s*([-+\d.e]+)", log)
            ]
        )
        table(
            out,
            tag + "-failed-coordinates.csv",
            ("coordinate", "value"),
            enumerate(coordinates),
        )
        print(
            tag,
            result.solver_status,
            convergence["residual"],
            result.message,
            flush=True,
        )
        return
    base, request, _ = load_case(case)
    species = request["reaction_system"]["species_ids"]
    model = epcsaft.Mixture(epcsaft.Parameters.from_mapping(changed(base, case, value)))
    states, terms = [], []
    for phase in range(result.phase_count):
        x = result.mole_fractions[phase * len(species) : (phase + 1) * len(species)]
        state = model.state(
            result.temperature, rho=result.molar_densities[phase], x=x
        )
        terms.append(
            (phase, state.born, state.debye_huckel, state.bulk_relative_permittivity)
        )
        states.extend(
            (
                phase,
                name,
                fraction,
                result.amounts[phase],
                result.molar_densities[phase],
                result.packing_fractions[phase],
                result.temperature,
                result.pressure,
            )
            for name, fraction in zip(species, x, strict=True)
        )
    table(
        out,
        tag + "-terms.csv",
        ("phase", "born_over_RT", "debye_huckel_over_RT", "bulk_relative_permittivity"),
        terms,
    )
    table(
        out,
        tag + "-state.csv",
        (
            "phase",
            "species",
            "mole_fraction",
            "amount_mol",
            "rho_mol_m3",
            "packing_fraction",
            "T_K",
            "P_Pa",
        ),
        states,
    )
    table(
        out,
        tag + "-coordinates.csv",
        ("coordinate", "value"),
        enumerate(result.coordinates),
    )
    print(
        tag, result.solver_status, convergence["residual"], result.message, flush=True
    )


def solve(mapping, problem, start=None):
    model = epcsaft.Mixture(epcsaft.Parameters.from_mapping(mapping))
    t = time.monotonic()
    with tempfile.TemporaryDirectory(prefix="ipopt-", dir=HERE) as tmp:
        previous_options = os.environ.get("EPCSAFT_IPOPT_OPTIONS_FILE")
        try:
            options = Path(tmp) / "ipopt.opt"
            logfile = Path(tmp) / "ipopt.log"
            options.write_text(f"output_file {logfile}\nfile_print_level 12\n")
            os.environ["EPCSAFT_IPOPT_OPTIONS_FILE"] = str(options)
            with epcsaft.trace() as capture:
                result = eq.solve_equilibrium(model, problem, start)
            log = logfile.read_text()
        finally:
            if previous_options is None:
                os.environ.pop("EPCSAFT_IPOPT_OPTIONS_FILE", None)
            else:
                os.environ["EPCSAFT_IPOPT_OPTIONS_FILE"] = previous_options
    return result, capture, time.monotonic() - t, log


@functools.cache
def archived():
    """Retired-chart member texts by exact name; #341 B12 retains these inputs once, byte-exact."""
    with tarfile.open(HERE / "results/retired-charts.tar.gz") as tar:
        return {m.name: tar.extractfile(m).read().decode() for m in tar if m.isfile()}


def saved_state(case, value):
    # A restarted stage seeds from the stage that ended at its start value: the coarse (top-level)
    # 0.35/2.25 states start fine/ and the fine/ 0.33/2.04 states start refined/ (README commands),
    # so the coarsest stage holding the value wins over its own re-solve in the next stage.
    return next((name for prefix in ("", "fine/", "refined/")
                 for name in [f"{prefix}case{case}-continuation-{value:.6f}-state.csv"] if name in archived()), None)


def seed_from_saved(problem, case, value):
    name = saved_state(case, value)
    assert name, f"no retained case {case} state at {value}"
    states = list(csv.DictReader(archived()[name].splitlines()))
    problem.P = eq.Free(float(states[0]["P_Pa"]))
    for phase_index, phase in enumerate(problem.phases):
        phase.composition_guess = tuple(
            float(r["mole_fraction"]) for r in states
            if int(r["phase"]) == phase_index and r["species"] in phase.support
        )
        phase.packing_guess = float(next(
            r["packing_fraction"] for r in states if int(r["phase"]) == phase_index
        ))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", type=int, choices=(1, 2), required=True)
    parser.add_argument("--continuation", action="store_true")
    parser.add_argument("--step", type=float, default=None)
    parser.add_argument("--from-value", type=float, default=None)
    parser.add_argument("--out", type=Path, default=HERE / "results")
    args = parser.parse_args()
    mapping, request, anchor = load_case(args.case)
    problem = problem_from_request(request, anchor)
    if not args.continuation:
        value = 0.0 if args.case == 1 else 2.0
        r, cap, sec, log = solve(mapping, problem)
        record(args.out, f"case{args.case}-direct", args.case, value, r, cap, sec, log)
        return
    value, target = (1.0, 0.0) if args.case == 1 else (8.0, 2.0)
    if args.from_value is not None:
        value = args.from_value
        seed_from_saved(problem, args.case, value)
    step = args.step or (0.05 if args.case == 1 else 0.25)
    family = epcsaft.ParameterFamily.BornCDielectric if args.case == 1 else epcsaft.ParameterFamily.IonicRegionRelativePermittivity
    model = epcsaft.Mixture(epcsaft.Parameters.from_mapping(changed(mapping, args.case, value)),
        active_parameters=[epcsaft.ActiveParameter(family, (0,))])
    t = time.monotonic()
    with epcsaft.trace() as capture:
        path = eq.continue_equilibrium(model, problem, problem,
            eq.Continuation(initial_step=step / (value - target), minimum_step=1e-5 / (value - target),
                            max_steps=1000, max_retries_per_step=20), parameter_targets=[target])
    seconds = time.monotonic() - t
    args.out.mkdir(parents=True, exist_ok=True)
    states, rows, coordinates, terms = [], [], [], []
    species = request["reaction_system"]["species_ids"]
    for parameters, result in zip(path.parameter_values, path.states, strict=True):
        physical = parameters[0]
        point_model = epcsaft.Mixture(epcsaft.Parameters.from_mapping(changed(mapping, args.case, physical)))
        for phase in range(result.phase_count):
            x = result.mole_fractions[phase * len(species):(phase + 1) * len(species)]
            state = point_model.state(result.temperature, rho=result.molar_densities[phase], x=x)
            terms.append((physical, phase, state.born, state.debye_huckel, state.bulk_relative_permittivity))
            states.extend((physical, phase, name, fraction, result.amounts[phase], result.molar_densities[phase],
                           result.packing_fractions[phase], result.temperature, result.pressure)
                          for name, fraction in zip(species, x, strict=True))
        rows.extend((physical, True, i, row.eqid, raw, classified) for i, (row, raw, classified) in enumerate(
            zip(result.rows, result.residuals, result.classified_residuals, strict=True)))
        coordinates.extend((physical, i, v) for i, v in enumerate(result.coordinates))
    if path.last_result is not None and path.attempted_path is not None:
        physical = value + path.attempted_path * (target - value)
        result = path.last_result
        rows.extend((physical, False, i, row.eqid, raw, classified) for i, (row, raw, classified) in enumerate(
            zip(result.rows, result.residuals, result.classified_residuals, strict=True)))
    table(args.out, f"case{args.case}-states.csv", ("value", "phase", "species", "mole_fraction", "amount_mol", "rho_mol_m3", "packing_fraction", "T_K", "P_Pa"), states)
    table(args.out, f"case{args.case}-rows.csv", ("value", "accepted", "row", "equation", "residual", "classified_residual"), rows)
    table(args.out, f"case{args.case}-coordinates.csv", ("value", "coordinate", "coordinate_value"), coordinates)
    table(args.out, f"case{args.case}-terms.csv", ("value", "phase", "born_over_RT", "debye_huckel_over_RT", "bulk_relative_permittivity"), terms)
    spans = capture.spans
    table(args.out, f"case{args.case}-trace.csv", ("span", "scope", "parent", "name", "values"),
          ((e.span, spans[e.span].name, spans[e.span].parent, e.name, json.dumps(e.values)) for e in capture.events
           if spans[e.span].name in ("continue_equilibrium", "equilibrium_continuation_attempt", "equilibrium_continuation_initial")
           or e.name == "ipopt_iteration"))
    (args.out / f"case{args.case}-path.json").write_text(json.dumps({
        "success": path.success, "stop": str(path.stop), "message": path.message,
        "reached_s": path.reached_path, "attempted_s": path.attempted_path,
        "cause": str(path.cause_stop), "cause_message": path.cause_message,
        "attempts": path.attempts, "retries": path.retries, "accepted_steps": path.accepted_steps,
        "wall_s": seconds,
        "cause_residual": path.cause_convergence.residual if path.cause_convergence and math.isfinite(path.cause_convergence.residual) else None,
    }, indent=2) + "\n")


def declarations():
    from epcsaft import Mixture, Parameters

    return [
        (Mixture(Parameters.from_mapping(d)), problem_from_request(q, a))
        for d, q, a in map(load_case, (1, 2))
    ]


if __name__ == "__main__":
    main()
