"""Nine-species reactive extraction and same-central-state sensitivity costs (#97)."""

import argparse
from collections import Counter
import csv
import hashlib
from importlib.metadata import distribution
import json
import math
import os
from pathlib import Path
import resource
import subprocess
import time
from urllib.parse import unquote, urlparse

started = time.perf_counter()
import epcsaft  # noqa: E402 - keep import latency separate from numerical work
from epcsaft import equilibrium as eq  # noqa: E402
IMPORT_MS = 1000 * (time.perf_counter() - started)

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parents[1]
CASE = json.loads((ROOT / "inputs/case.json").read_text())
SPECIES = CASE["species"]
N = len(SPECIES)
RESULTS = ROOT / "results"
FAILURES = []
PROVENANCE = {}
OBSERVABLES = [(1, "lithium-carrier"), (1, "sodium-carrier"), (0, "hydrogen-cation")]
DIRECTIONS = {
    "organic_inventory": [CASE["fresh_organic_per_litre_mol"].get(s, 0.0) for s in SPECIES],
    "LiCl": [float(s in ("lithium-cation", "chloride-anion")) for s in SPECIES],
    "NaCl": [float(s in ("sodium-cation", "chloride-anion")) for s in SPECIES],
}


def write_table(name, rows):
    with (RESULTS / name).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(dict.fromkeys(k for row in rows for k in row)), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def retain_run():
    PROVENANCE["failures"] = FAILURES
    (RESULTS / "provenance.json").write_text(json.dumps(PROVENANCE, indent=2) + "\n")


def record_state_failure(state, label):
    FAILURES.append(dict(case=label, message=state.message, solver=state.solver_status,
        requested_tolerance_met=state.requested_tolerance_met, residuals=list(state.residuals),
        classified_residuals=list(state.classified_residuals)))


def shifted(feed, direction, step):
    return {s: feed[s] + step * direction[i] for i, s in enumerate(SPECIES)}


def problem(feed, reactive=True):
    # Disjoint source supports make these declared guesses identical to the source inventory mapping.
    phases = []
    for phase in CASE["phases"]:
        total = math.fsum(feed[s] for s in phase["support"])
        phases.append(eq.Phase(phase["name"], kind=phase["kind"], support=phase["support"],
                               amount=eq.Free(), composition_guess=[feed[s] / total for s in phase["support"]]))
    reactions = [eq.Reaction(r["stoichiometry"], equilibrium_constant=math.exp(r["ln_K_density_basis"]))
                 for r in CASE["reactions"]] if reactive else []
    return eq.Problem(T=CASE["T_K"], P=CASE["P_Pa"], phases=phases,
                      feed=eq.Amounts(feed), reactions=reactions, row_tolerance=1e-10)


def values(state):
    return [state.amounts[p] * state.mole_fractions[p * N + SPECIES.index(s)] for p, s in OBSERVABLES]


def inspect(model, request, state, label, reactive=True, branches=False):
    row = dict(case=label, success=state.success, solver=state.solver_status,
               iterations=state.iterations, classified_residual=state.convergence.residual)
    if not state.success:
        record_state_failure(state, label)
        return row
    compiled = eq.compile_problem(model, request)
    amounts, potentials, pressures, root_errors = [], {}, [], []
    for p, phase in enumerate(request.phases):
        x = state.mole_fractions[p * N:(p + 1) * N]
        amounts.append([state.amounts[p] * v for v in x])
        physical = model.state(request.T, rho=state.molar_densities[p], x=x)
        pressures.append(abs(physical.pressure / request.P - 1))
        if branches:
            independent = model.state(request.T, P=request.P, x=x, phase="liquid")
            root_errors.append(abs(independent.molar_density / physical.molar_density - 1))
        row[f"rho_{phase.name}_mol_m3"] = physical.molar_density
        row[f"volume_{phase.name}_L"] = state.amounts[p] / physical.molar_density * 1000
        for i, s in enumerate(SPECIES):
            row[f"n_{phase.name}_{s}_mol"] = amounts[-1][i]
            if s in phase.support:
                potentials[s] = math.log(x[i] * physical.molar_density) + physical.residual_chemical_potential_over_rt[i]
    delta = [amounts[0][i] + amounts[1][i] - request.feed.values[s] for i, s in enumerate(SPECIES)]
    balance = max(abs(math.fsum(b * d for b, d in zip(basis, delta))) for basis in CASE["conserved_moieties"])
    affinity = max(abs(math.fsum(nu * potentials[s] for s, nu in r["stoichiometry"].items()) - r["ln_K_density_basis"])
                   for r in CASE["reactions"])
    point = state.coordinates[:compiled.primal_size]
    gradient, jacobian = compiled.objective_gradient(point), compiled.jacobian(point)
    rows = [i for i, r in enumerate(compiled.rows) if r.kind.name not in ("Chemical", "Pressure")]
    kkt = max(abs(g + math.fsum(mu * jacobian[i * compiled.primal_size + j]
                  for mu, i in zip(state.constraint_multipliers, rows, strict=True))) for j, g in enumerate(gradient))
    charge = max(abs(math.fsum(z * n for z, n in zip(CASE["charges"], a))) for a in amounts)
    absent = max(abs(amounts[p][i]) for p, phase in enumerate(request.phases)
                 for i, s in enumerate(SPECIES) if s not in phase.support)
    row.update(balance_mol=balance, charge_mol=charge, support_mol=absent,
               original_reaction_log_affinity=affinity, pressure_relative=max(pressures), kkt=kkt,
               curvature=state.reduced_hessian_minimum, rank=state.jacobian_rank,
               dimension=state.unknown_count, density_root_relative=max(root_errors) if root_errors else None)
    accepted = (state.requested_tolerance_met and state.topology_event.name == "None_"
                and state.jacobian_rank == state.unknown_count and state.reduced_hessian_minimum > 0
                and balance <= 1e-9 and charge <= 1e-9 and absent <= 1e-14
                and max(pressures) <= 1e-8 and kkt <= 1e-7
                and max(root_errors, default=0.0) <= 1e-7
                and (affinity <= 1e-6 if reactive else max(abs(v) for v in delta) <= 1e-9))
    row["physical_screens_pass"] = accepted
    if not accepted:
        FAILURES.append(dict(case=label, message="physical screens failed", diagnostics=row))
    return row


def main():
    RESULTS.mkdir(parents=True, exist_ok=True)
    installed = json.loads(distribution("epcsaft").read_text("direct_url.json"))
    installed_wheel = Path(unquote(urlparse(installed["url"]).path))
    wheel_sha256 = hashlib.sha256(installed_wheel.read_bytes()).hexdigest()
    PROVENANCE.update(completed=False,
        engine_commit=subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=REPO, text=True).strip(),
        checkout_dirty=bool(subprocess.check_output(["git", "diff", "HEAD", "--", "engine/native"], cwd=REPO)),
        wheel_sha256=wheel_sha256, import_ms=IMPORT_MS,
        input_sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in (ROOT / "inputs").glob("*.json")},
        cpu_affinity=sorted(os.sched_getaffinity(0)), nice=os.getpriority(os.PRIO_PROCESS, 0),
        warmups=1, samples=5, fd_steps=[0.001, 0.0005],
        threads={s: os.environ.get(s) for s in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")},
        limits="Adopted calibrated model; declared disjoint supports; fixed T/P and lnK; local stability, no global phase search. No external-package runtime comparison.")
    retain_run()
    start = time.perf_counter()
    parameters = epcsaft.Parameters.from_json(ROOT / "inputs/parameters.json", components=SPECIES)
    model = epcsaft.Mixture(parameters)
    setup_ms = 1000 * (time.perf_counter() - start)
    request = problem(CASE["feed_mol"])
    start = time.perf_counter()
    central = eq.solve_equilibrium(model, request)
    first_ms = 1000 * (time.perf_counter() - start)
    if not central.success:
        record_state_failure(central, "central")
        raise RuntimeError(central.message)
    physical_rows = [inspect(model, request, central, "central", branches=True)]
    li_feed = CASE["condition"]["li_mol_l_aqueous"] * CASE["condition"]["aqueous_volume_l"]
    na_feed = CASE["condition"]["na_mol_l_aqueous"] * CASE["condition"]["aqueous_volume_l"]
    sweep = []
    for oa in (0.5, 0.75, 1.0, 1.25, 1.5):
        q = problem(shifted(CASE["feed_mol"], DIRECTIONS["organic_inventory"], oa - 1))
        s = central if oa == 1 else eq.solve_equilibrium(model, q)
        row = inspect(model, q, s, f"OA={oa}", branches=True)
        row["organic_to_aqueous_feed_volume_ratio"] = oa
        if s.success:
            li, na, _ = values(s)
            row.update(li_extraction_percent=100 * li / li_feed, na_extraction_percent=100 * na / na_feed,
                       D_Li_feed_volume=(li / oa) / (li_feed - li), D_Na_feed_volume=(na / oa) / (na_feed - na))
        sweep.append(row)
    for key in ("Li", "Na"):
        reference = CASE["calibration"][f"D_{key}"]
        if abs(sweep[-1].get(f"D_{key}_feed_volume", math.inf) / reference - 1) > 1e-7:
            FAILURES.append(dict(case="calibration", message=f"D_{key} differs from retained fitted endpoint"))
    write_table("inventory-sweep.csv", sweep)

    equilibrated_feed = {s: math.fsum(central.amounts[p] * central.mole_fractions[p * N + i]
                                    for p in (0, 1)) for i, s in enumerate(SPECIES)}
    control = []
    for added in (0.0, 0.0025, 0.005, 0.01, 0.02):
        for reactive in (True, False):
            q = problem(shifted(equilibrated_feed, DIRECTIONS["LiCl"], added), reactive)
            s = eq.solve_equilibrium(model, q, start=central)
            row = inspect(model, q, s, f"LiCl={added}/reactive={reactive}", reactive, branches=True)
            row.update(added_LiCl_mol=added, chemistry="coupled" if reactive else "frozen")
            if added == 0 and s.success and max(abs(a - b) for a, b in zip(values(s), values(central))) > 1e-9:
                FAILURES.append(dict(case=row["case"], message="controls do not share the same baseline"))
            control.append(row)
    write_table("reaction-control.csv", control)

    observables = [eq.SolvedStateObservable(eq.SolvedStateObservableKind.PhaseComponentAmount, p, SPECIES.index(s))
                   for p, s in OBSERVABLES]

    def exact(batched=True):
        compiled = eq.compile_problem(model, request)
        directions = [eq.SolvedStateActionDirection(0, 0, direction, []) for direction in DIRECTIONS.values()]
        batches = eq.solved_state_first_actions(compiled, central, observables, directions) if batched else [
            eq.solved_state_actions(compiled, central, eq.SolvedStateActionRequest(observables, [d])) for d in directions]
        return [[v.action for v in batch.results] for batch in batches], batches

    def finite_difference(h):
        derivatives, states = [], []
        for direction in DIRECTIONS.values():
            pair = []
            for sign in (-1, 1):
                q = problem(shifted(request.feed.values, direction, sign * h))
                s = eq.solve_equilibrium(model, q, start=central)
                states.append((q, s))
                pair.append(values(s) if s.success else [math.nan] * len(observables))
            derivatives.append([(b - a) / (2 * h) for a, b in zip(*pair)])
        return derivatives, states

    derivatives, batches = exact()
    for name, batch in zip(DIRECTIONS, batches):
        if any(v.status != eq.SolvedStateActionStatus.Available for v in batch.results):
            FAILURES.append(dict(case=f"exact/{name}", statuses=[str(v.status) for v in batch.results],
                rank=batch.diagnostics.rank, dimension=batch.diagnostics.state_dimension,
                linear_residual=batch.diagnostics.maximum_linear_residual))
            raise RuntimeError(f"unavailable exact direction {name}: {[str(v.status) for v in batch.results]}")
    derivative_rows = []
    for h in (0.001, 0.0005):
        differences, states = finite_difference(h)
        physical_rows.extend(inspect(model, q, s, f"FD/{h}/{i}") for i, (q, s) in enumerate(states))
        for k, (name, batch) in enumerate(zip(DIRECTIONS, batches)):
            for j, (phase, species) in enumerate(OBSERVABLES):
                a, fd = derivatives[k][j], differences[k][j]
                relative = abs(fd - a) / abs(a) if a else math.nan
                derivative_rows.append(dict(direction=name, observable=species, phase=phase, h=h,
                    step_unit="O/A" if name == "organic_inventory" else "mol", exact=a, finite_difference=fd,
                    absolute_error=abs(fd - a), relative_error=relative,
                    linear_residual=batch.diagnostics.maximum_linear_residual,
                    condition=batch.diagnostics.condition_estimate))
                if h == 0.0005 and (not math.isfinite(relative) or relative > 5e-6):
                    FAILURES.append(dict(case=f"derivative/{name}/{species}", message="A1 comparison exceeds 5e-6"))
    write_table("derivatives.csv", derivative_rows)

    timings = []
    methods = ("cold_seeded_solve", "empty_tape_solve", "previous_state_solve",
               "exact_three_directions", "exact_three_directions_legacy", "finite_difference_three_directions")
    for repetition in range(-1, 5):
        for method in (methods if repetition % 2 else tuple(reversed(methods))):
            # Empty-tape timing excludes model construction, but includes its first recordings.
            timed_model = epcsaft.Mixture(parameters) if method == "empty_tape_solve" else model
            wall, cpu = time.perf_counter(), time.process_time()
            if method.endswith("solve"):
                outcome = eq.solve_equilibrium(timed_model, request, start=central if method.startswith("previous") else None)
            elif method.startswith("exact"):
                outcome = exact(method == "exact_three_directions")
            else:
                outcome = finite_difference(0.0005)
            wall_ms, cpu_ms = 1000 * (time.perf_counter() - wall), 1000 * (time.process_time() - cpu)
            if method.endswith("solve"):
                results = [(request, outcome)]
            else:
                results = outcome[1] if method.startswith("finite") else []
            physical_rows.extend(inspect(model, q, s, f"timing/{method}/{repetition}/{i}")
                                 for i, (q, s) in enumerate(results))
            valid = all(s.success for _, s in results) if results else all(
                v.status == eq.SolvedStateActionStatus.Available for b in outcome[1] for v in b.results)
            if not valid:
                FAILURES.append(dict(case=method, message="timed workload unavailable"))
            if repetition >= 0:
                timings.append(dict(method=method, repetition=repetition, wall_ms=wall_ms, cpu_ms=cpu_ms,
                    nonlinear_solves=len(results), iterations=sum(s.iterations for _, s in results), success=valid))
    write_table("timings.csv", timings)
    write_table("physical-checks.csv", physical_rows)
    with epcsaft.trace() as capture:
        traced = eq.solve_equilibrium(model, request)
    if not traced.success:
        record_state_failure(traced, "traced")
    recordings, phase_signatures, derivative_coordinates = Counter(), Counter(), Counter()
    for event in capture.events:
        span = capture.spans[event.span]
        parent = capture.spans[span.parent].name if span.parent is not None else "none"
        if event.name == "recording":
            recordings[parent] += int(event.values[0])
        if event.name == "predicate_signature" and parent == "equilibrium_prepare":
            phase_signatures[tuple(event.values)] += 1
        if event.name == "derivative_coordinates":
            derivative_coordinates[parent, int(event.values[0])] += 1
    profile = {}
    child_time = [0.0] * len(capture.spans)
    for span in capture.spans:
        if span.parent is not None:
            child_time[span.parent] += span.wall_us
    for i, span in enumerate(capture.spans):
        row = profile.setdefault(span.name, dict(operation=span.name, calls=0, inclusive_ms=0.0, self_ms=0.0))
        row["calls"] += 1
        row["inclusive_ms"] += span.wall_us / 1000
        row["self_ms"] += (span.wall_us - child_time[i]) / 1000
    write_table("profile.csv", sorted(profile.values(), key=lambda r: -r["self_ms"]))
    PROVENANCE.update(completed=not FAILURES, model_setup_ms=setup_ms, first_solve_ms=first_ms,
        peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, cpu_affinity=sorted(os.sched_getaffinity(0)),
        recordings_by_parent=dict(recordings), phase_signature_count=len(phase_signatures),
        phase_evaluations=sum(phase_signatures.values()),
        derivative_coordinate_requests=[dict(parent=parent, coordinates=count, calls=calls)
            for (parent, count), calls in sorted(derivative_coordinates.items())],
        traced_success=traced.success)
    retain_run()
    print(json.dumps(PROVENANCE, indent=2), flush=True)
    assert not FAILURES, "Inspect retained failures; do not relax criteria"


def declarations():
    from epcsaft import Mixture, Parameters

    return (
        Mixture(Parameters.from_json(ROOT / 'inputs/parameters.json', components=SPECIES)),
        problem(CASE['feed_mol']),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, default=RESULTS,
                        help="Separate result directory for a matched installed-wheel replay")
    RESULTS = parser.parse_args().results
    try:
        main()
    except Exception as error:
        FAILURES.append(dict(exception=type(error).__name__, message=str(error)))
        PROVENANCE["completed"] = False
        retain_run()
        raise
