"""FeOs 0.10.1 native shared-case values and complete call timings (#222)."""
import json, math, os, subprocess, sys, time
from pathlib import Path

# This file's name must not shadow the installed FeOs package.
sys.path = [p for p in sys.path if Path(p).resolve() != Path(__file__).parent]
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "JULIA_NUM_THREADS"):
    os.environ[name] = "1"
started = time.perf_counter_ns()
import numpy as np
from feos import (BinaryRecord, Contributions, EquationOfState, EquationOfStateAD,
                  Identifier, Parameters, PhaseEquilibrium, Property, PureRecord, State, set_num_threads)
from si_units import JOULE, KELVIN, KB, METER, MOL, NAV, PASCAL, RGAS
set_num_threads(1)
import_ns = time.perf_counter_ns() - started
RHO = MOL / METER**3
CONTROL = dict(tol=1e-12, max_iter=1000)
BUBBLE = dict(tol_inner=1e-12, tol_outer=1e-12, max_iter_inner=1000, max_iter_outer=1000)
TERM = {"Hard Sphere": "hard_chain", "Hard Chain": "hard_chain", "Dispersion": "dispersion",
        "Association": "association", "Ionic": "debye_huckel", "Born": "born", "Dipole": "dipole"}


def dipole_factor(inputs):
    e, f = inputs["constants"]["Engine"], inputs["constants"]["FeOs"]
    return math.sqrt((1e-21 / 299792458.)**2 / (4 * math.pi * e["epsilon0"] * 1e-30 * 1e-19) * f["kB"] / e["kB"])


def prepare(inputs, q):
    spec = inputs["models"][q["case"]]
    ids = [Identifier(name=c["id"]) for c in spec["components"]]
    records = []
    epsilon = 78.09
    if q["case"] == "C3" and q.get("electrostatic", "matched") == "matched":
        e, f = inputs["constants"]["Engine"], inputs["constants"]["FeOs"]
        epsilon *= (f["e"] / e["e"])**2 * e["epsilon0"] / f["epsilon0"] * e["kB"] / f["kB"]
    for ident, c in zip(ids, spec["components"]):
        kw = dict(m=c["m"], sigma=c["sigma"], epsilon_k=c["epsilon"])
        if c.get("kappa_ab"):
            kw["association_sites"] = [dict(kappa_ab=c["kappa_ab"], epsilon_k_ab=c["epsilon_ab"], na=1., nb=1.)]
        if q["case"] == "C3":
            kw.update(z=c.get("z", 0), permittivity_record={"ExperimentalData": {"data": [[298., epsilon], [299., epsilon]]}})
        else:
            kw.update(mu=c.get("mu", 0), nmu=c.get("nmu", 1))
            if q.get("dipole_convention") == "matched": kw["mu"] *= dipole_factor(inputs)
        records.append(PureRecord(ident, c["molar_mass"] * 1000, **kw))
    pairs = spec.get("pairs", [(0, 1, spec["k"])] if "k" in spec else [])
    binaries = [BinaryRecord(ids[i], ids[j], k_ij=q.get("kcl", k) if q["case"] == "C3" and (i, j) == (0, 2) else k)
                for i, j, k in pairs]
    pars = Parameters.from_records(records, binaries)
    eos = EquationOfState.epcsaft(pars, epcsaft_variant="revised") if q["case"] == "C3" else EquationOfState.pcsaft(pars)
    pure, binary = pars.to_json_str(pretty=False)
    return eos, dict(pure_records=json.loads(pure), binary_records=json.loads(binary) if binary else []), spec


def state(eos, q, x):
    kw = dict(temperature=q["T"] * KELVIN, composition=np.array(x))
    if "rho" in q:
        kw["density"] = q["rho"] * RHO
    else:
        kw.update(pressure=q["P"] * PASCAL, density_initialization="liquid")
    return State(eos, **kw)


def phase(s, label, full=False, fixed=False):
    out = dict(phase=label, rho=s.density / RHO, x=s.molefracs.tolist(), pressure=s.pressure() / PASCAL)
    if full or fixed:
        out.update(ar=s.molar_helmholtz_energy(Contributions.Residual) / (RGAS * s.temperature),
                   Z=s.compressibility(), lnphi=s.ln_phi().tolist())
        if out["Z"] <= 0: out["lnphi_reason"] = "Native lnphi contains log(Z); Z <= 0 is outside its real-valued domain"
    if full:
        shares = dict.fromkeys(TERM.values(), 0.)
        for name, value in s.residual_molar_helmholtz_energy_contributions():
            shares[TERM[name]] += value / (RGAS * s.temperature)
        out.update(mu_res=(s.chemical_potential(Contributions.Residual) / (RGAS * s.temperature)).tolist(),
                   pressure_density_derivative=s.dp_drho() / (PASCAL / RHO), contributions=shares)
    return out


def bubble(eos, q):
    return PhaseEquilibrium.bubble_point(eos, q["T"] * KELVIN, np.array(q["x"]), **BUBBLE)


def gradient(eos, spec, q):
    c = spec["components"]
    if q["case"] == "C1":
        p = c[0]
        names = ["m", "sigma", "epsilon_k", "mu"]
        parameters = np.array([p["m"], p["sigma"], p["epsilon"], p.get("mu", 0.)])
        value, g, ok = Property.vapor_pressure_derivatives(EquationOfStateAD.PcSaftNonAssoc, names, parameters, np.array([[q["T"]]], dtype=float))
        extra = dict(mu=float(g[0, 3]))
    elif q["case"] == "C2":
        names = ["k_ij", "epsilon_k_ab1", "kappa_ab1"]
        parameters = np.array([v for p in c for v in (p["m"], p["sigma"], p["epsilon"], p.get("mu", 0.),
                              p.get("kappa_ab", 0.), p.get("epsilon_ab", 0.), float(bool(p.get("kappa_ab"))), float(bool(p.get("kappa_ab"))))] + [spec["k"]])
        guess = bubble(eos, q).vapor.pressure() / PASCAL
        value, g, ok = Property.bubble_point_pressure_derivatives(EquationOfStateAD.PcSaftFull, names, parameters, np.array([[q["T"], q["x"][0], guess]]))
        extra = dict(pressure_guess=guess)
    else:
        raise ValueError("Only C1/C2 exact gradients are selected")
    if not ok[0]:
        return dict(P=float(value[0]), gradients=[None]*3, gradient_names=names[:3], gradient_native_converged=False,
                    gradient_message="Native Property failed; zero-filled gradient storage is unavailable")
    return dict(P=float(value[0]), gradients=g[0, :3].tolist(), gradient_names=names[:3],
                gradient_extra=extra, gradient_native_converged=bool(ok[0]))


def evaluate(q, prepared, helper=None):
    eos, _, spec = prepared
    op = q["operation"]
    if op == "gradient":
        return gradient(eos, spec, q), []
    if op in ("saturation", "bubble", "lle"):
        eq = (PhaseEquilibrium.pure(eos, q["T"] * KELVIN, **CONTROL) if op == "saturation" else
              bubble(eos, q) if op == "bubble" else
              PhaseEquilibrium.tp_flash(eos, q["T"] * KELVIN, q["P"] * PASCAL, np.array(q["x"]) * MOL, **CONTROL))
        states, labels = [eq.liquid, eq.vapor], ["liquid", "vapor"]
        out = dict(P=states[1].pressure() / PASCAL)
        if op == "lle":
            rich = int(states[1].molefracs[0] > states[0].molefracs[0])
            out["beta"] = eq.vapor_phase_fraction if rich else 1. - eq.vapor_phase_fraction
            states, labels = [states[rich], states[1-rich]], ["water-rich", "water-poor"]
    else:
        if op not in ("fixed", "reference", "activity"):
            raise ValueError(f"Unknown operation {op}")
        MW = spec["components"][0]["molar_mass"]
        x = q.get("x")
        if x is None:
            feed = np.array([1. / MW, q["m"], q["m"]])
            x = (feed / feed.sum()).tolist()
        s = state(eos, q, x)
        states, labels, out = [s], [q.get("phase", "liquid")], dict(P=s.pressure() / PASCAL)
        if op == "activity":
            rm = q.get('reference_m', 0.)
            feed = np.array([1./MW, rm, rm]); ref = state(eos, q, (feed/feed.sum()).tolist())
            pure = state(eos,q,[1.,0.,0.]) if rm else ref
            payload = dict(lnphi=s.ln_phi().tolist(), reference_lnphi=ref.ln_phi().tolist(),
                           pure_lnphi=pure.ln_phi().tolist(), m=q["m"], MW=MW, xwater=x[0])
            out.update(helper_payload=payload, reference_inputs=dict(T=q["T"], P=q["P"], m=rm, x=ref.molefracs.tolist()))
            if helper is not None:
                t = time.perf_counter_ns()
                helper.stdin.write(json.dumps(payload) + "\n")
                helper.stdin.flush()
                answer = helper.stdout.readline()
                out["transport_ns"] = time.perf_counter_ns() - t
                if not answer:
                    raise RuntimeError("Julia contraction helper closed stdout")
                out.update(json.loads(answer))
            out["reference_state"] = ref
    out["phases"] = [phase(s, label, fixed=op == "fixed") for s, label in zip(states, labels)]
    return out, list(zip(states, labels))


def finite(value, absent, path):
    if isinstance(value, float) and not math.isfinite(value):
        absent.append(dict(quantity=path, reason="Native nonfinite value; retain state for diagnosis"))
        return None
    if isinstance(value, dict):
        return {k: finite(v, absent, f"{path}.{k}") for k, v in value.items()}
    if isinstance(value, list):
        return [finite(v, absent, f"{path}[{i}]") for i, v in enumerate(value)]
    return value


def main():
    inputs = json.loads(Path(sys.argv[1]).read_text())
    rep = int(sys.argv[3]) if len(sys.argv) > 3 else None
    records, timings, absent, models, serializations = [], [], [], {}, {}
    helper = None; warm_ns=0
    if rep is not None and any(q.get("timing") == "T2" and q["operation"] == "activity" and q.get("diameter") != "sigma" for q in inputs["requests"]):
        here = Path(__file__).resolve().parents[1]
        cache = Path(os.environ.get("XDG_CACHE_HOME", Path.home() / ".cache")) / "epcsaft-comparators/222"
        helper = subprocess.Popen(["systemd-run", "--user", "--quiet", "--scope", "-p", "MemoryMax=6G", "timeout", "180",
                                   str(cache / "julia-1.13.1/bin/julia"), "--startup-file=no", f"--project={here / 'inputs'}",
                                   str(here / "scripts/clapeyron.jl"), "--contract"], stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True)
    def timing(q, kind, ns, sample=0, transport=0):
        timings.append(dict(id=q["id"], tier=q.get("timing", "setup"), selection="warm" if kind == "sample" else "cold",
                            kind=kind, rep=rep, sample=sample, elapsed_ns=ns, transport_ns=transport,
                            work="State+ar/Z/lnphi" if q.get("timing") == "T1" else q["operation"] + ("+Julia contraction/transport" if q["case"] == "C3" else "")))
    try:
        for q in inputs["requests"]:
            if q["case"] == "C3" and q.get("diameter") == "sigma":
                absent.append(dict(id=q["id"], reason="FeOs revised Ionic uses 0.88sigma packing diameters; no sigma screening option"))
                continue
            if rep is not None and not q.get("timing"):
                continue
            key = f"{q['case']}:{q.get('electrostatic', 'matched')}:{q.get('kcl', -.25)}:{q.get('dipole_convention', 'literal')}"
            try:
                if key not in models:
                    t = time.perf_counter_ns()
                    models[key] = prepare(inputs, q)
                    if rep is not None: timing(q, "model", time.perf_counter_ns() - t)
                    serializations[key] = models[key][1]
                t = time.perf_counter_ns()
                out, states = evaluate(q, models[key], helper)
                if rep is not None:
                    timing(q, "first-call", time.perf_counter_ns() - t, transport=out.get("transport_ns", 0))
                    warm_start=time.perf_counter_ns()
                    for _ in range(5): evaluate(q, models[key], helper)
                    for sample in range(25):
                        t = time.perf_counter_ns()
                        sample_out, _ = evaluate(q, models[key], helper)
                        timing(q, "sample", time.perf_counter_ns() - t, sample, sample_out.get("transport_ns", 0))
                    warm_ns+=time.perf_counter_ns()-warm_start
                out["phases"] = [phase(s, label, full=True) for s, label in states]
                if "reference_state" in out: out["reference"] = phase(out.pop("reference_state"), "liquid", full=True)
                records.append(finite(dict(id=q["id"], case=q["case"], T=q["T"], **out), absent, q["id"]))
            except Exception as exc:
                records.append(dict(id=q["id"], case=q["case"], T=q["T"], error=str(exc)))
            print(q["id"], flush=True)
    finally:
        if helper is not None:
            helper.stdin.close()
            try: helper.wait(timeout=30)
            except subprocess.TimeoutExpired: helper.kill(); helper.wait()
    if rep is not None: timing(dict(id="FeOs-import", case="setup", operation="import"), "import", import_ns)
    absent.append(dict(reason="FeOs Python exposes no achieved density/association/equilibrium residual, iterations, site fractions or Delta matrix; AD returns only convergence boolean"))
    result = dict(records=records, timings=timings, warm_elapsed_ns=warm_ns, absent_diagnostics=absent, model_serializations=serializations, dipole_mu_factor=dipole_factor(inputs),
                  quantity_constants=dict(NA=NAV * MOL, kB=KB / (JOULE / KELVIN), R=RGAS / (JOULE / MOL / KELVIN)))
    Path(sys.argv[2]).write_text(json.dumps(result, allow_nan=False, separators=(",", ":")))


if __name__ == "__main__":
    main()
