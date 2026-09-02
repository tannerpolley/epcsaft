"""Performance-floor measurements on the instrumented wheel (perf/instrumentation).

Records, per case: wall time per call, tape recordings, tape ops, recording time,
Jacobian/Hessian/Forward calls, sweep time, association solves and iterations.
Environment is recorded with every run. OMP_NUM_THREADS must be 1.
"""
from __future__ import annotations
import copy, hashlib, json, os, platform, statistics, subprocess, sys, time
from pathlib import Path

assert os.environ.get("OMP_NUM_THREADS") == "1", "pin OMP_NUM_THREADS=1"
import epcsaft
from epcsaft import equilibrium, _core, unit_registry as u

REPO = Path("~/Workspaces/Engineering/ePC-SAFT-project")
FIX = REPO / "analyses/2026-fixed-pressure-anchor/inputs/fixtures"
MEA = REPO / "analyses/2026-mea-parameter-estimation/inputs"
OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(".")
WT = Path(os.environ.get("EPCSAFT_WORKTREE", "."))

def counters(): return dict(_core._instrumentation())
def reset(): _core._instrumentation_reset()
def delta(a, b): return {k: b[k] - a[k] for k in b}

def timed(fn, reps, warm=3):
    for _ in range(warm): fn()
    reset(); t0 = time.perf_counter();
    for _ in range(reps): fn()
    wall = (time.perf_counter() - t0) / reps
    c = counters()
    return wall, {k: v / reps for k, v in c.items()}

def environment():
    so = next(Path(epcsaft.__file__).parent.glob("_core*.so"))
    wheel = next(iter(sorted((WT / "build/environment-wheel").glob("epcsaft-*.whl"))), None)
    sha = hashlib.sha256(wheel.read_bytes()).hexdigest() if wheel else None
    ldd = subprocess.run(["ldd", str(so)], capture_output=True, text=True).stdout
    cxx = subprocess.run(["c++", "--version"], capture_output=True, text=True).stdout.splitlines()[0]
    cpu = [l.split(":", 1)[1].strip() for l in open("/proc/cpuinfo") if l.startswith("model name")][0]
    commit = subprocess.run(["git", "-C", str(WT), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    return {"wheel": wheel.name if wheel else None, "wheel_sha256": sha, "commit": commit, "ldd": ldd,
            "compiler": cxx, "cpu": cpu, "python": platform.python_version(),
            "OMP_NUM_THREADS": os.environ["OMP_NUM_THREADS"], "platform": platform.platform()}

def state_cases():
    out = {}
    models = {
        "neutral_binary": epcsaft.Mixture(epcsaft.Parameters.from_json(FIX / "neutral_binary.json")),
        "associating_binary": epcsaft.Mixture(epcsaft.Parameters.from_json(FIX / "strongly_associating_binary.json")),
        "aqueous_electrolyte_3": epcsaft.Mixture(epcsaft.Parameters.from_json(FIX / "aqueous_electrolyte.json")),
        "mea_nine_species": epcsaft.Mixture(epcsaft.Parameters.from_json(MEA / "neutral-selected-parameters.json")),
    }
    comps = {"neutral_binary": (0.4, 0.6), "associating_binary": (0.1385, 0.8615),
             "aqueous_electrolyte_3": (0.98, 0.01, 0.01),
             "mea_nine_species": None}
    for name, m in models.items():
        n = len(m.component_ids)
        x = comps[name] or tuple([0.9 - 0.01 * (n - 3)] + [0.01] * (n - 3) + [0.02, 0.02]) if comps[name] else None
        if x is None:
            # nine species: water/MEA/CO2 dominant, ions at 1e-3 with electroneutrality
            ids = list(m.component_ids); x = [0.0] * n
            x[ids.index("water")] = 0.80; x[ids.index("monoethanolamine")] = 0.10; x[ids.index("carbon-dioxide")] = 0.05
            x[ids.index("protonated-monoethanolamine")] = 0.02; x[ids.index("carbamate-anion")] = 0.02
            x[ids.index("bicarbonate-anion")] = 0.005; x[ids.index("hydronium-cation")] = 0.005
            s = sum(x); x = tuple(v / s for v in x)
        T = {"mea_nine_species": 313.15, "associating_binary": 376.88, "aqueous_electrolyte_3": 298.15}.get(name, 300.0)
        rho = 3.0e4 if n > 2 else 1.0e4
        def state_rho(): return m.state(T=T * u.kelvin, rho=rho * u.mole / u.meter**3, x=x)
        wall, c = timed(state_rho, 50)
        out[f"{name}.state_rho"] = {"n": n, "T": T, "wall_us": wall * 1e6, **c}
        # native derivative-rich evaluations
        native = m._native_model(None)
        V = 1.0 / rho; amounts = list(x)
        def deriv2(): return native.evaluate_phase_derivatives(T, amounts, V, [], 2)
        def deriv1(): return native.evaluate_phase_derivatives(T, amounts, V, [], 1)
        def grepe_eval(): return _core.evaluate_parameterized_phase(native, T, tuple(amounts), V, 1.0e-8, 0.74)
        for label, fn in (("phase_derivatives_order2", deriv2), ("phase_derivatives_order1", deriv1), ("grepe_evaluate_parameterized_phase", grepe_eval)):
            try:
                wall, c = timed(fn, 20)
                out[f"{name}.{label}"] = {"n": n, "wall_us": wall * 1e6, **c}
            except Exception as exc:
                out[f"{name}.{label}"] = {"error": f"{type(exc).__name__}: {str(exc)[:120]}"}
    return out

def discovery_case():
    m = epcsaft.Mixture(epcsaft.Parameters.from_json(FIX / "neutral_binary.json"))
    problem = equilibrium.PhaseEquilibriumProblem(equilibrium.FixedTPFeed(300.0 * u.kelvin, 1.0e5 * u.pascal, (0.5, 0.5)))
    reset(); t0 = time.perf_counter(); r = equilibrium.solve(m, problem); wall = time.perf_counter() - t0
    c = counters()
    return {"wall_s": wall, "phases": len(r.phases), "eos_evaluation_count": r.diagnostics.eos_evaluation_count,
            "optimizer_solve_count": r.diagnostics.optimizer_solve_count, "optimizer_iteration_count": r.diagnostics.optimizer_iteration_count, **c}

def reactive_case():
    doc = json.loads((MEA / "candidate-1-mixed-n8-input.json").read_text())
    params = epcsaft.Parameters.from_json(MEA / "neutral-selected-parameters.json")
    m = epcsaft.Mixture(params)
    def fix(o):
        if isinstance(o, dict): return {k: fix(v) for k, v in o.items()}
        if isinstance(o, list): return [fix(v) for v in o]
        if isinstance(o, str): return o.replace("installed-provider-eos", "installed-eos").replace("provider", "eos")
        return o
    results = []
    for obs in doc["observations"][:4]:
        request = fix(copy.deepcopy(obs["request"]))
        problem = equilibrium.general_reactive_equilibrium_problem_from_mapping(request)
        reset(); t0 = time.perf_counter(); r = equilibrium.solve(m, problem); wall = time.perf_counter() - t0
        c = counters(); ev = dict(r.evidence)
        results.append({"identity": obs["identity"], "wall_s": wall, "status": r.status,
                        "eos_evaluation_count": ev.get("eos_evaluation_count"), "ipopt_iterations": ev.get("optimizer_iteration_count") or ev.get("ipopt_iterations"),
                        "evidence_keys": sorted(k for k in ev if "count" in k or "iter" in k)[:12], **c})
    return results

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    result = {"environment": environment()}
    result["state"] = state_cases()
    result["discovery_neutral_binary"] = discovery_case()
    try:
        result["reactive_bubble_point"] = reactive_case()
    except Exception as exc:
        result["reactive_bubble_point"] = {"error": f"{type(exc).__name__}: {str(exc)[:300]}"}
    (OUT / "floor-measurements.json").write_text(json.dumps(result, indent=2, default=str) + "\n")
    for k, v in result["state"].items():
        if "wall_us" in v:
            print(f"{k:48s} {v['wall_us']:9.1f} us  rec={v['recordings']:.1f} ops={v['tape_ops']:.0f} rec_ns={v['recording_ns']:.0f} sweep_ns={v['sweep_ns']:.0f} hess={v['hessian_calls']:.1f} assoc_it={v['association_iterations']:.1f}")
        else:
            print(k, v)
    d = result["discovery_neutral_binary"]; print("discovery", {k: d[k] for k in ("wall_s", "phases", "eos_evaluation_count", "recordings", "hessian_calls")})
    print("reactive", [ (r.get("identity"), r.get("wall_s"), r.get("eos_evaluation_count"), r.get("recordings"), r.get("hessian_calls"), r.get("error")) for r in (result["reactive_bubble_point"] if isinstance(result["reactive_bubble_point"], list) else [result["reactive_bubble_point"]])])

if __name__ == "__main__":
    main()
