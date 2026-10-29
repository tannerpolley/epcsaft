"""Declare fixed-temperature bubble states and read native solves and exact actions."""
import argparse
import csv
import hashlib
import json
import math
import sys
import time
from pathlib import Path

import epcsaft
from epcsaft import equilibrium as eq

ROOT = Path(__file__).resolve().parents[1]
CONFIG = json.loads((ROOT / "inputs/reactions.json").read_text())
S = tuple(CONFIG["species"])
T = CONFIG["temperature_K"]


def declarations():
    """Construct every chemistry/permittivity variant without solving."""
    xw = float(next(r for r in csv.DictReader((ROOT / "inputs/measurements.csv").open()) if r["selected"] == "True")["x_water"])
    return [(model(v), problem(xw, v)) for v in ("primary", "acid-only", "water-only-permittivity", "no-dissociation")]


def model(variant):
    doc = json.loads((ROOT / "inputs/parameters.json").read_text())
    if variant == "water-only-permittivity":
        for c in doc["components"]:
            if c["component_id"] == "acetic-acid":
                c["coefficients"] = [p for p in c["coefficients"] if p["family"] != "relative_permittivity"]
    species = S[:2] if variant == "no-dissociation" else S[:4] if variant == "acid-only" else S
    return epcsaft.Mixture(epcsaft.Parameters.from_mapping(doc, components=species))


def problem(xw, variant="primary", tolerance=None, reference_log_shifts=(0., 0.)):
    species = S[:2] if variant == "no-dissociation" else S[:4] if variant == "acid-only" else S
    reactions, neutral = [], None
    if variant != "no-dissociation":
        reactions.append(eq.Reaction({S[1]: -1., S[2]: 1., S[3]: 1.}, name="acetic-acid-dissociation",
            reference=eq.ReactionReference(eq.ReferenceSourceBasis.CommonMolalityInfiniteDilution,
                                           CONFIG["acid_reference_pressure_Pa"]),
            correlation=eq.ReactionLogPolynomial(math.log(CONFIG["Ka"]) + reference_log_shifts[0], 0., 0., 0., T, T, T,
                                                "Harned-Ehlers-1933-Table-II-50-C")))
        if len(species) == 5:
            reactions.append(eq.Reaction({S[0]: -1., S[2]: 1., S[4]: 1.}, name="water-ionization",
                reference=eq.ReactionReference(eq.ReferenceSourceBasis.CommonMolalityInfiniteDilution,
                                               CONFIG["water_reference_pressure_Pa"]),
                correlation=eq.ReactionLogPolynomial(-CONFIG["pKw"] * math.log(10) + reference_log_shifts[1], 0., 0., 0., T, T, T,
                                                    "IAPWS-R11-24-Table-3-50-C-formal-H")))
        ratios = {S[1]: 1., S[2]: 2. if len(species) == 5 else 1., S[3]: 1.}
        if len(species) == 5:
            ratios[S[4]] = 1.
        neutral = eq.NeutralReference(solvent=S[0], solute_molality_ratios=ratios,
                                      standard_molality_mol_per_kg=CONFIG["standard_molality_mol_per_kg"])
    guess = [xw, 1 - xw - 2e-5, 1e-5, 1e-5, 1e-16][:len(species)]
    return eq.Problem(T=T, P=eq.Free(11000.), row_tolerance=tolerance or CONFIG["row_tolerance"],
        feed=eq.Amounts({S[0]: xw, S[1]: 1 - xw}),
        phases=[eq.Phase("L", kind="liquid", support=species, composition_guess=guess),
                eq.Phase("V", kind="vapor", support=S[:2], amount=eq.Pinned(0.), composition_guess=[.8, .2])],
        reactions=reactions, neutral_reference=neutral)


def estimate(e):
    return {k: getattr(e, k) for k in ("absolute", "available", "resolved", "cause")}


def solve(xw, variant, tighter=False, actions=False, reference_log_shifts=(0., 0.)):
    mix = model(variant)
    q = problem(xw, variant, CONFIG["tighter_row_tolerance"] if tighter else None, reference_log_shifts)
    compiled = eq.compile_problem(mix, q)
    r = eq.solve_equilibrium(mix, q)
    row = dict(x_water=xw, variant=variant, tighter=tighter, success=r.success,
               requested_tolerance_met=r.requested_tolerance_met, solver=str(r.solver_status),
               row_residual=r.convergence.residual, iterations=r.iterations,
               P_Pa=r.pressure, compositions=list(r.mole_fractions), densities_mol_m3=list(r.molar_densities),
               phase_amounts_mol=list(r.amounts))
    if not r.success:
        return row
    n = len(mix.component_ids)
    x, y = list(r.mole_fractions[:n]), list(r.mole_fractions[n:])
    st = mix.state(T, rho=r.molar_densities[0], x=x)
    vap = mix.state(T, rho=r.molar_densities[1], x=y)
    amounts = [math.fsum(a * r.mole_fractions[p * n + i] for p, a in enumerate(r.amounts)) for i in range(n)]
    charges = [0, 0, 1, -1, -1][:n]
    row.update(charge_error_mol=math.fsum(a * z for a, z in zip(amounts, charges)),
               acid_inventory_error_mol=amounts[1] + (amounts[3] if n > 2 else 0) - (1 - xw),
               water_inventory_error_mol=amounts[0] + (amounts[4] if n == 5 else 0) - xw,
               hydrogen_inventory_error_mol=(amounts[2] - amounts[3] - (amounts[4] if n == 5 else 0)) if n > 2 else 0.,
               pressure_closure_relative=max(abs(st.pressure / r.pressure - 1), abs(vap.pressure / r.pressure - 1)),
               neutral_log_fugacity_difference=[math.log(x[i] / y[i]) + st.log_fugacity_coefficient[i] - vap.log_fugacity_coefficient[i] for i in range(2)],
               reaction_log_residuals=[a.log_quotient - (math.log(CONFIG["Ka"]) if a.reaction_index == 0 else -CONFIG["pKw"] * math.log(10)) - reference_log_shifts[a.reaction_index] for a in eq.reaction_activities(mix, r) if a.phase_index == 0],
               acid_dissociation_fraction=x[3] / (x[1] + x[3]) if n > 2 else 0.)
    row["checks_met"] = all([
        r.requested_tolerance_met, r.pressure > 0, all(d > 0 for d in r.molar_densities),
        abs(row["charge_error_mol"]) <= CONFIG["charge_tolerance_mol"],
        max(abs(row[k]) for k in ("acid_inventory_error_mol", "water_inventory_error_mol", "hydrogen_inventory_error_mol")) <= CONFIG["inventory_tolerance_mol"],
        row["pressure_closure_relative"] <= CONFIG["relative_pressure_tolerance"],
        max(map(abs, row["neutral_log_fugacity_difference"])) <= CONFIG["log_fugacity_tolerance"],
        max(map(abs, row["reaction_log_residuals"]), default=0.) <= CONFIG["log_reaction_tolerance"],
    ])
    row["reaction_references"] = [dict(reaction_index=d.reaction_index,
        source_pressure_Pa=d.reference_pressure_pa, contraction=d.terminal_contraction,
        status=str(d.status), cause=d.message) for d in r.neutral_reference_diagnostics]
    if actions:
        observable = eq.SolvedStateObservable(eq.SolvedStateObservableKind.Pressure)
        direction = eq.SolvedStateActionDirection(0., 0., [1., -1.] + [0.] * (n - 2), [])
        batch = eq.solved_state_actions(compiled, r, eq.SolvedStateActionRequest([observable], [direction]))
        row["pressure_actions"] = [{"order": a.order, "value": a.value, "action": a.action,
                                    "status": str(a.status), "cause": a.diagnostic_message} for a in batch.results]
        numerical = eq.equilibrium_numerical_estimates(mix, r)
        row["numerical_estimate"] = dict(pressure=estimate(numerical.pressure),
            source_reference_method=numerical.source_reference_method, cause=numerical.cause,
            row_inner=list(numerical.row_inner), rank=numerical.rank, dimension=numerical.dimension,
            condition=numerical.condition, backward_error=numerical.backward_error)
    return row


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--declare-only", action="store_true")
    args = parser.parse_args()
    if args.declare_only:
        declarations()
        return
    wheel, digest = epcsaft.provenance()
    if not isinstance(wheel, Path):
        raise RuntimeError("Use the qualified installed-wheel environment, not an editable import")
    build_file = wheel.parent.parent / "build-source.json"
    build = json.loads(build_file.read_text())
    if build["source_revision"] != CONFIG["source_revision"] or build["engine_git_status"]:
        raise RuntimeError("Wheel build record does not identify the accepted clean Engine source")
    identity = dict(source_revision=CONFIG["source_revision"], wheel=str(wheel), wheel_sha256=digest,
                    python=sys.executable, package=epcsaft.__file__, extension=epcsaft._core.__file__,
                    build_source_record=str(build_file), build_source_record_sha256=hashlib.sha256(build_file.read_bytes()).hexdigest(),
                    build_command=build["build_command"], compiler=build["compiler"],
                    script_sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((ROOT / "scripts").glob("*")) if p.suffix in (".py", ".cpp")},
                    input_sha256={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((ROOT / "inputs").glob("*")) if p.is_file()})
    identity["loaded_sha256"] = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in Path(epcsaft.__file__).parent.glob("*.so")}
    identity["loaded_native_maps"] = sorted({line.split()[-1] for line in Path("/proc/self/maps").read_text().splitlines()
                                            if "epcsaft" in line and ".so" in line})
    (ROOT / "results/executable.json").write_text(json.dumps(identity, indent=2) + "\n")
    measurements = [r for r in csv.DictReader((ROOT / "inputs/measurements.csv").open()) if r["selected"] == "True"]
    rows = []
    for measured in measurements:
        center, radius = float(measured["x_water"]), float(measured["U95_x"])
        points = [center - radius + 2 * radius * i / (CONFIG["interval_grid_points"] - 1) for i in range(CONFIG["interval_grid_points"])]
        for xw in points:
            for variant in ("primary", "acid-only", "water-only-permittivity"):
                started = time.perf_counter()
                row = solve(xw, variant, actions=variant == "primary")
                row.update(measurement_row=int(measured["row"]), elapsed_s=time.perf_counter() - started)
                rows.append(row)
                print(f"row {row['measurement_row']} xw={xw:.7f} {variant}: {row['P_Pa']} Pa checks={row.get('checks_met')}", flush=True)
                (ROOT / "results/states.json").write_text(json.dumps(rows, indent=2, allow_nan=False) + "\n")
        for variant, tighter in (("no-dissociation", False), ("primary", True)):
            row = solve(center, variant, tighter=tighter)
            row["measurement_row"] = int(measured["row"])
            rows.append(row)
        (ROOT / "results/states.json").write_text(json.dumps(rows, indent=2, allow_nan=False) + "\n")
    failures = [r for r in rows if not r.get("checks_met", False)]
    (ROOT / "results/failures.json").write_text(json.dumps(dict(attempted=len(rows), failures=failures), indent=2) + "\n")
    if failures:
        raise RuntimeError("Numerical failures retained separately; stop for diagnosis")


if __name__ == "__main__":
    main()
