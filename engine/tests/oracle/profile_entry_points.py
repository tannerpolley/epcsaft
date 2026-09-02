"""Profile every public entry point on the instrumented wheel (perf/instrumentation).

Per call: wall time, tape recordings, recording time, tape ops, mean tape width
(domain) and range, Jacobian/Hessian/Forward calls, sweep time, association
solves and iterations, density-root iterations. Counters reset before each call.
Run from engine/ with OMP_NUM_THREADS=1; writes JSON and a Markdown table.
"""
from __future__ import annotations
import copy, hashlib, json, os, platform, subprocess, sys, time, traceback
from dataclasses import replace
from pathlib import Path

assert os.environ.get("OMP_NUM_THREADS") == "1", "pin OMP_NUM_THREADS=1"
sys.path.insert(0, os.getcwd())  # engine/ so `tests.*` fixtures import
import epcsaft
from epcsaft import equilibrium, regression, eos as eos_module, _core, unit_registry as u
from tests.eos._cases import (NEUTRAL_BINARY, NEUTRAL_TERNARY_IMMISCIBLE, STRONGLY_ASSOCIATING_BINARY,
                              AQUEOUS_ELECTROLYTE, MULTINEUTRAL_REACTIVE, neutral_pure)
from tests.equilibrium._cases import finite_phase_case, solved_pressure_case, ionic_reaction_case, temperature_boundary_case

REPO = Path("~/Workspaces/Engineering/ePC-SAFT-project")
MEA = REPO / "analyses/2026-mea-parameter-estimation/inputs"
OUT = Path(sys.argv[1]); OUT.mkdir(parents=True, exist_ok=True)
WT = Path(os.environ.get("EPCSAFT_WORKTREE", "."))
ROWS: list[dict] = []

def counters(): return dict(_core._instrumentation())
def reset(): _core._instrumentation_reset()

def profile(name, fn, reps=1, warm=1, note=""):
    try:
        for _ in range(warm): fn()
        reset(); t0 = time.perf_counter()
        for _ in range(reps): fn()
        wall = (time.perf_counter() - t0) / reps
        c = {k: v / reps for k, v in counters().items()}
        row = {"entry_point": name, "wall_ms": wall * 1e3, **c,
               "mean_tape_width": (c["tape_domain"] / c["recordings"]) if c["recordings"] else None,
               "mean_tape_ops": (c["tape_ops"] / c["recordings"]) if c["recordings"] else None,
               "note": note}
    except Exception as exc:
        row = {"entry_point": name, "error": f"{type(exc).__name__}: {str(exc)[:160]}", "note": note}
    ROWS.append(row); print(f"{name:60s} " + (f"{row['wall_ms']:10.3f} ms rec={row['recordings']:8.0f} hess={row['hessian_calls']:8.0f} assoc_it={row['association_iterations']:8.0f} dens_it={row['density_iterations']:6.0f} w={row['mean_tape_width'] or 0:5.1f}" if "wall_ms" in row else row["error"]))
    return row

def environment():
    so = next(Path(epcsaft.__file__).parent.glob("_core*.so"))
    wheel = next(iter(sorted((WT / "build/environment-wheel").glob("epcsaft-*.whl"))), None)
    return {"wheel": wheel.name if wheel else None, "wheel_sha256": hashlib.sha256(wheel.read_bytes()).hexdigest() if wheel else None,
            "commit": subprocess.run(["git", "-C", str(WT), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip(),
            "ldd": subprocess.run(["ldd", str(so)], capture_output=True, text=True).stdout,
            "compiler": subprocess.run(["c++", "--version"], capture_output=True, text=True).stdout.splitlines()[0],
            "cpu": [l.split(":", 1)[1].strip() for l in open("/proc/cpuinfo") if l.startswith("model name")][0],
            "python": platform.python_version(), "OMP_NUM_THREADS": os.environ["OMP_NUM_THREADS"], "platform": platform.platform()}

def main():
    env = environment()
    K = u.kelvin; Pa = u.pascal; RHO = u.mole / u.meter**3
    neutral = epcsaft.Mixture(epcsaft.Parameters.from_mapping(NEUTRAL_BINARY))
    assoc = epcsaft.Mixture(epcsaft.Parameters.from_mapping(STRONGLY_ASSOCIATING_BINARY))
    elec_params = epcsaft.Parameters.from_mapping(AQUEOUS_ELECTROLYTE); elec = epcsaft.Mixture(elec_params)
    nine = epcsaft.Mixture(epcsaft.Parameters.from_json(MEA / "neutral-selected-parameters.json"))
    ids = list(nine.component_ids); x9 = [0.0] * 9
    for k, v in {"water": 0.80, "monoethanolamine": 0.10, "carbon-dioxide": 0.05, "protonated-monoethanolamine": 0.02, "carbamate-anion": 0.02, "bicarbonate-anion": 0.005, "hydronium-cation": 0.005}.items(): x9[ids.index(k)] = v
    x9 = tuple(v / sum(x9) for v in x9)

    # --- Mixture.state at fixed density and fixed pressure, both phases ---
    cases = [("neutral", neutral, 300.0, (0.4, 0.6), 1.0e4, 1.0e5),
             ("associating", assoc, 376.88, (0.1385, 0.8615), 1.0e4, 101330.0),
             ("electrolyte3", elec, 298.15, (0.98, 0.01, 0.01), 3.0e4, 1.0e5),
             ("nine_species", nine, 313.15, x9, 3.0e4, 1.0e5)]
    for label, m, T, x, rho, P in cases:
        profile(f"Mixture.state rho | {label}", lambda: m.state(T=T * K, rho=rho * RHO, x=x), reps=20, warm=3)
        for phase in ("liquid", "vapor"):
            profile(f"Mixture.state P {phase} | {label}", lambda: m.state(T=T * K, P=P * Pa, x=x, phase=phase), reps=5, warm=2)
    # --- Mixture.reference ---
    profile("Mixture.reference | electrolyte3", lambda: elec.reference(T=298.15 * K, P=1.0 * u.bar, solvent_x=(1.0, 0.0, 0.0), phase="liquid",
            ion_component_ids=("sodium-cation", "chloride-anion"), formula_stoichiometry=(1, 1), reference_convention="molality-infinite-dilution"), reps=3)
    ref = elec.reference(T=298.15 * K, P=1.0 * u.bar, solvent_x=(1.0, 0.0, 0.0), phase="liquid", ion_component_ids=("sodium-cation", "chloride-anion"), formula_stoichiometry=(1, 1), reference_convention="molality-infinite-dilution")
    profile("ElectrolyteReference.state(molality) | electrolyte3", lambda: ref.state(molality=0.1 * u.mole / u.kilogram), reps=5)
    # --- observable_block / observable_descriptor ---
    active_n = epcsaft.ActiveParameterSet(epcsaft.Parameters.from_mapping(NEUTRAL_BINARY), ("pair/ethane/methane/k_ij",))
    profile("Mixture.observable_descriptor | neutral", lambda: neutral.observable_descriptor(active_n, identity="eos.fixed-state"), reps=20)
    profile("Mixture.observable_block fixed-state | neutral", lambda: neutral.observable_block(active_n, identity="eos.fixed-state", T=300.0 * K, rho=1.0e4 * RHO, x=(0.4, 0.6)), reps=10)
    profile("Mixture.observable_block fixed-pressure-fugacity | neutral", lambda: neutral.observable_block(active_n, identity="eos.fixed-pressure-fugacity", T=300.0 * K, P=1.0e5 * Pa, x=(0.4, 0.6), phase="liquid"), reps=5)
    active_e = epcsaft.ActiveParameterSet(elec_params, ("model/relative_permittivity/ion_fraction_suppression_coefficient",))
    profile("Mixture.observable_block fixed-state | electrolyte3", lambda: elec.observable_block(active_e, identity="eos.fixed-state", T=298.15 * K, rho=3.0e4 * RHO, x=(0.98, 0.01, 0.01)), reps=10)
    # --- Pure.saturation ---
    prop_params = epcsaft.Parameters.from_mapping(neutral_pure("propane")); propane = epcsaft.Pure(prop_params)
    profile("Pure.saturation | propane 300 K", lambda: propane.saturation(T=300.0 * K), reps=3)
    sat = propane.saturation(T=300.0 * K)
    profile("Pure.saturation with start | propane 300 K", lambda: propane.saturation(T=300.0 * K, start=epcsaft.SaturationStart(sat.vapor.molar_density_mol_m3 * RHO, sat.liquid.molar_density_mol_m3 * RHO, sat.saturation_pressure_pa * Pa)), reps=3)
    active_p = epcsaft.ActiveParameterSet(prop_params, ("component/propane/segment_count", "component/propane/segment_diameter", "component/propane/dispersion_energy_over_k"))
    profile("Pure.observable_block pure-saturation (active) | propane", lambda: propane.observable_block(active_p, identity="eos.pure-saturation", T=300.0 * K), reps=3)
    assoc_pure = epcsaft.Pure(epcsaft.Parameters.from_mapping(__import__("tests.eos._cases", fromlist=["ASSOCIATING_PURE"]).ASSOCIATING_PURE))
    profile("Pure.saturation | associating pure 300 K", lambda: assoc_pure.saturation(T=300.0 * K), reps=2)
    # --- equilibrium.solve: declared fixed-T,P (one phase and three phase) ---
    profile("equilibrium.solve declared fixed-T,P one phase | neutral", lambda: equilibrium.solve(neutral, equilibrium.PhaseEquilibriumProblem(
        equilibrium.FixedTPFeed(300.0 * K, 1.0e5 * Pa, (1.0, 1.0), equilibrium.DeclaredTopology((equilibrium.PhaseSpecification("vapor", equilibrium.AllComponents()),))))), reps=3)
    tern = epcsaft.Mixture(epcsaft.Parameters.from_mapping(NEUTRAL_TERNARY_IMMISCIBLE))
    starts = ((0.90, 0.06, 0.04), (0.11, 0.83, 0.06), (0.14, 0.12, 0.74)); fr = (0.26, 0.33, 0.41)
    guesses = tuple(equilibrium.DeclaredPhaseInitialGuess(i, s, float(tern.state(T=300.0 * K, P=1.0e5 * Pa, x=s, phase="liquid").molar_density.to("mole / meter**3").magnitude), fr[i]) for i, s in enumerate(starts))
    spec3 = equilibrium.FixedTPFeed(300.0 * K, 1.0e5 * Pa, (1.0, 1.0, 1.0), equilibrium.DeclaredTopology(tuple(equilibrium.PhaseSpecification("liquid", equilibrium.AllComponents()) for _ in range(3))))
    profile("equilibrium.solve declared fixed-T,P three phase | ternary", lambda: equilibrium.solve(tern, equilibrium.PhaseEquilibriumProblem(spec3, initial_guesses=guesses)), reps=2)
    # --- fixed-T and fixed-P boundaries ---
    two = equilibrium.DeclaredTopology((equilibrium.PhaseSpecification("liquid", equilibrium.AllComponents()), equilibrium.PhaseSpecification("vapor", equilibrium.AllComponents())))
    profile("equilibrium.solve fixed-T coexistence | associating binary", lambda: equilibrium.solve(assoc, equilibrium.PhaseEquilibriumProblem(equilibrium.FixedTCoexistence(300.0 * K, two, 0, (0.5, 0.5)))), reps=2)
    profile("equilibrium.solve fixed-P coexistence (bubble T) | associating binary", lambda: equilibrium.solve(assoc, equilibrium.PhaseEquilibriumProblem(equilibrium.FixedPCoexistence(1200.0 * Pa, (295.0 * K, 305.0 * K), two, 0, (0.5, 0.5)))), reps=1)
    multi = epcsaft.Mixture(epcsaft.Parameters.from_mapping(MULTINEUTRAL_REACTIVE))
    twoN = equilibrium.DeclaredTopology((equilibrium.PhaseSpecification("liquid", equilibrium.AllComponents()), equilibrium.PhaseSpecification("vapor", equilibrium.NeutralComponents())))
    profile("equilibrium.solve fixed-T coexistence | charged 4-comp neutral vapor", lambda: equilibrium.solve(multi, equilibrium.PhaseEquilibriumProblem(equilibrium.FixedTCoexistence(298.15 * K, twoN, 0, (0.70, 0.20, 0.05, 0.05)))), reps=1)
    # --- discovery ---
    profile("equilibrium.solve discovery | neutral binary", lambda: equilibrium.solve(neutral, equilibrium.PhaseEquilibriumProblem(equilibrium.FixedTPFeed(300.0 * K, 1.0e5 * Pa, (0.5, 0.5)))), reps=1)
    profile("equilibrium.solve discovery | charged electrolyte3", lambda: equilibrium.solve(elec, equilibrium.PhaseEquilibriumProblem(equilibrium.FixedTPFeed(298.15 * K, 1.0e5 * Pa, (0.98, 0.01, 0.01)))), reps=1)
    # --- reactive / GREPE ---
    _, fm, fp = finite_phase_case()
    profile("equilibrium.solve GREPE finite two-phase LLE | ternary", lambda: equilibrium.solve(fm, fp), reps=1)
    sp, sm, spr = solved_pressure_case()
    profile("equilibrium.solve GREPE finite+incipient solved P | 4-comp charged", lambda: equilibrium.solve(sm, spr), reps=1)
    sa = epcsaft.ActiveParameterSet(sp, (sp.parameter_specs[0].identity,))
    profile("equilibrium.solve GREPE finite+incipient solved P + active | 4-comp charged", lambda: equilibrium.solve(sm, spr, active_parameters=sa), reps=1)
    im, ip = ionic_reaction_case()
    profile("equilibrium.solve GREPE homogeneous reaction | electrolyte3", lambda: equilibrium.solve(im, ip), reps=1)
    tp, tm, tpr, tid = temperature_boundary_case("liquid")
    profile("equilibrium.solve GREPE solved-T bubble | binary", lambda: equilibrium.solve(tm, tpr), reps=1)
    doc = json.loads((MEA / "candidate-1-mixed-n8-input.json").read_text())
    def fix(o):
        if isinstance(o, dict): return {k: fix(v) for k, v in o.items()}
        if isinstance(o, list): return [fix(v) for v in o]
        if isinstance(o, str): return o.replace("installed-provider-eos", "installed-eos").replace("provider", "eos")
        return o
    req = equilibrium.general_reactive_equilibrium_problem_from_mapping(fix(copy.deepcopy(doc["observations"][0]["request"])))
    profile("equilibrium.solve GREPE reactive bubble point | nine species", lambda: equilibrium.solve(nine, req), reps=1)
    # --- regression ---
    obs_block = propane.observable_block(active_p, identity="eos.pure-saturation", T=300.0 * K)
    observations = tuple(float(v) for v in obs_block.values if v is not None)[:3]
    owner = eos_module.compile_regression_evaluator({"model": propane, "family": "pure-saturation", "temperature_k": 300.0,
        "row_identities": ("saturation_pressure_pa", "vapor_molar_density_mol_per_m3", "liquid_molar_density_mol_per_m3")}, active_p)
    profile("eos.compile_regression_evaluator pure-saturation | propane", lambda: eos_module.compile_regression_evaluator({"model": propane, "family": "pure-saturation", "temperature_k": 300.0,
        "row_identities": ("saturation_pressure_pa", "vapor_molar_density_mol_per_m3", "liquid_molar_density_mol_per_m3")}, active_p), reps=3)
    profile("regression.least_squares pure-saturation 3 params | propane", lambda: regression.least_squares((owner,), observations, (active_p.values[0] + 0.01, active_p.values[1] - 0.01, active_p.values[2]),
        bounds=((1.0, 1.0, 50.0), (5.0, 5.0, 500.0)), residual_scales=observations, maximum_iterations=20, maximum_elapsed_time_seconds=30.0), reps=1)
    profile("regression.curve_fit pure-saturation 3 params | propane", lambda: regression.curve_fit((owner,), observations, active_p.values,
        bounds=((1.0, 1.0, 50.0), (5.0, 5.0, 500.0)), residual_scales=observations, maximum_iterations=4, maximum_elapsed_time_seconds=30.0), reps=1)
    import math
    baseline = equilibrium.solve(sm, spr, active_parameters=sa)
    def anchored_problem():
        fs = baseline.continuation_state.phase(spr.phases[0].identity); total = None
        for row, conserved in zip(spr.reaction_system.balance_matrix, spr.reaction_system.conserved_totals, strict=True):
            bal = math.fsum(c * f for c, f in zip(row, fs.mole_fractions, strict=True))
            if abs(bal) > 1e-12: total = conserved / bal; break
        phases = tuple(replace(ph, start=equilibrium.FinitePhaseStart(tuple(total * f for f in st.mole_fractions), st.molar_volume_m3_per_mol)) if ph.amount_role == "finite" else ph
                       for ph, st in zip(spr.phases, baseline.continuation_state.phases, strict=True))
        return replace(spr, phases=phases, continuation_state=baseline.continuation_state)
    anchored = anchored_problem()
    profile("equilibrium.compile_regression_evaluator | 4-comp charged", lambda: equilibrium.compile_regression_evaluator({"model": sm, "problem": anchored}, sa), reps=2)
    eq_owner = equilibrium.compile_regression_evaluator({"model": sm, "problem": anchored}, sa)
    eq_obs = tuple(float(v) for v in baseline.values)
    profile("regression.least_squares equilibrium owner 1 param | 4-comp charged", lambda: regression.least_squares((eq_owner,), eq_obs, (sa.values[0] + 0.05,),
        bounds=((sa.values[0] - 5.0,), (sa.values[0] + 5.0,)), maximum_iterations=20, maximum_elapsed_time_seconds=120.0), reps=1)

    (OUT / "entry-point-profile.json").write_text(json.dumps({"environment": env, "rows": ROWS}, indent=2, default=str) + "\n")
    cols = ["entry_point", "wall_ms", "recordings", "recording_ns", "mean_tape_ops", "mean_tape_width", "jacobian_calls", "hessian_calls", "forward_calls", "sweep_ns", "association_solves", "association_iterations", "density_iterations"]
    lines = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for r in ROWS:
        if "error" in r: lines.append(f"| {r['entry_point']} | error: {r['error']} |" + " |" * (len(cols) - 2)); continue
        vals = []
        for c in cols:
            v = r.get(c)
            if c == "entry_point": vals.append(str(v))
            elif c in ("recording_ns", "sweep_ns"): vals.append(f"{v/1e6:.2f} ms")
            elif v is None: vals.append("")
            else: vals.append(f"{v:.3f}" if c == "wall_ms" else f"{v:.1f}" if isinstance(v, float) and c in ("mean_tape_ops", "mean_tape_width") else f"{v:.0f}")
        lines.append("| " + " | ".join(vals) + " |")
    header = f"Wheel {env['wheel_sha256']}, commit {env['commit']}, {env['compiler']}, {env['cpu']}, OMP_NUM_THREADS={env['OMP_NUM_THREADS']}\n\n"
    (OUT / "entry-point-profile.md").write_text(header + "\n".join(lines) + "\n")
    print("wrote", OUT / "entry-point-profile.md")

if __name__ == "__main__":
    main()
