"""Native envelope trace and supplied liquid TPD evaluations for paper figures F4 and F3."""

import hashlib
import json
import subprocess
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from epcsaft import Mixture, Parameters, equilibrium as eq, _core

HERE = Path(__file__).resolve().parents[1]
ROOT = HERE.parents[1]
S = json.loads((HERE / "inputs/study.json").read_text())
OUT = HERE / "results"


def write(name, rows):
    pd.DataFrame(rows).to_csv(OUT / name, index=False)


def state(r):
    values = [r.temperature, r.pressure, r.molar_densities[0], r.molar_densities[-1], r.mole_fractions[-2]]
    values += [r.success, r.requested_tolerance_met, max(map(abs, r.residuals)), str(r.convergence.stop), r.message]
    keys = "T_K P_Pa rho_feed rho_incipient x0_incipient success tolerance_met residual stop message".split()
    return dict(zip(keys, values), bubble=r.molar_densities[0] > r.molar_densities[-1])


sys.path.insert(0, str((ROOT / S["tpd_source"]).parent))
from feos_compare import engine_parameters


def declarations():
    m = Mixture(Parameters.from_json(ROOT / S["envelope_packet"], components=("methane", "ethane")))
    c = eq.EnvelopeCondition()
    c.kind, c.target = eq.EnvelopeCoordinate.LogPressure, np.log(S["start_pressure_Pa"])
    p = eq.Problem(
        T=eq.Free(S["start_temperature_K"]),
        P=eq.Free(S["start_pressure_Pa"]),
        feed=eq.Amounts(dict(zip(m.component_ids, S["feed_mol"]))),
        row_tolerance=S["row_tolerance"],
        envelope=c,
        phases=[eq.Phase("L", kind="liquid"), eq.Phase("V", kind="vapor", amount=eq.Pinned(0.0))],
    )
    w = Mixture(engine_parameters())
    wp = eq.Problem(
        T=S["tpd_temperature_K"],
        P=S["tpd_pressure_Pa"],
        row_tolerance=1e-10,
        feed=eq.Amounts(dict(zip(w.component_ids, S["feed_mol"]))),
        phases=[eq.Phase("L", kind="liquid")],
    )
    return m, p, w, wp, eq.Continuation(**S["continuation"])


def main():
    m, p, w, wp, o = declarations()
    trace = eq.continue_equilibrium(m, p, p, o)
    rows = [dict(index=i, **state(r)) for i, r in enumerate(trace.states)]
    write("envelope.csv", rows)
    critical = [
        dict(**state(r), eigenvalues=r.critical.eigenvalues, cubic=r.critical.cubic, zero_band=r.critical.zero_band)
        for r in trace.critical_points
    ]
    meta = dict(stop=str(trace.stop), attempts=trace.attempts, retries=trace.retries, critical=critical)
    (OUT / "trace.json").write_text(json.dumps(meta, indent=2) + "\n")
    assert trace.success and len(critical) == 1 and all(r["success"] and r["tolerance_met"] for r in rows)
    checks = []
    quantities = ("P_Pa", "rho_feed", "rho_incipient", "x0_incipient")
    for ref in pd.read_csv(ROOT / S["envelope_reference"]).to_dict("records"):
        candidates = []
        for a, b in zip(rows, rows[1:]):
            t = float(ref["T_K"])
            branch = a["bubble"] == bool(ref["bubble"]) and b["bubble"] == a["bubble"]
            if branch and min(a["T_K"], b["T_K"]) <= t <= max(a["T_K"], b["T_K"]):
                f = (t - a["T_K"]) / (b["T_K"] - a["T_K"])
                errors = {k: abs(a[k] + f * (b[k] - a[k]) - float(ref[k])) / abs(float(ref[k])) for k in quantities}
                candidates.append(dict(T_K=t, bubble=ref["bubble"], segment=a["index"], **errors))
        checks.append(min(candidates, key=lambda x: max(x[k] for k in quantities)))
    write("crosschecks.csv", checks)
    assert max(r[k] for r in checks for k in quantities) <= S["curve_relative_tolerance"], checks
    refs = pd.read_csv(ROOT / S["critical_reference"])
    ref = refs[(refs.source == "feos") & (refs.z0 == 0.5)].iloc[0]
    values = [critical[0][k] for k in ("T_K", "rho_feed", "P_Pa")]
    expected = ref[["T_K", "rho_mol_m3", "P_Pa"]].to_numpy(float)
    assert np.allclose(values, expected, rtol=S["critical_relative_tolerance"], atol=0)
    m, p = w, wp
    r = eq.solve_equilibrium(m, p)
    assert r.success and r.requested_tolerance_met
    compiled = eq.compile_problem(m, p)
    search = eq.stability_test(compiled, r, eq.StabilityOptions())
    trials = [t for t in search.trials if t.phase.kind == _core.PhaseKind.Liquid]
    write(
        "search.csv",
        [
            dict(seed=t.seed, x_water=t.composition[0], resolved_descent=t.resolved_descent, **t.quantities)
            for t in trials
        ],
    )
    v = np.linspace(*S["tpd_logit_limits"], S["tpd_grid_points"])
    x = sorted(set(1 / (1 + np.exp(-v))) | {t.composition[0] for t in trials})
    phases = []
    for value in x:
        phase = _core.PhaseSpec()
        phase.name, phase.kind, phase.support = "liquid-trial", _core.PhaseKind.Liquid, [0, 1]
        phase.composition_guess = [value, 1 - value]
        phases.append(phase)
    options = eq.StabilityOptions()
    options.evaluate_supplied_trials, options.trials = True, phases
    path = eq.stability_test(compiled, r, options)
    write(
        "tpd.csv",
        [dict(x_water=t.composition[0], evidence=t.evidence, message=t.message, **t.quantities) for t in path.trials],
    )
    assert len(path.trials) == len(x) and all("tpd_over_rt" in t.quantities for t in path.trials)
    assert any(t.resolved_descent for t in trials)
    sources = [ROOT / v for k, v in S.items() if k.endswith(("_packet", "_reference", "_source", "_parameters"))]
    sources += [HERE / "inputs/study.json", Path(_core.__file__)]
    meta = dict(
        source_commit=subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(), central=state(r)
    )
    meta["sha256"] = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    (OUT / "provenance.json").write_text(json.dumps(meta, indent=2) + "\n")


if __name__ == "__main__":
    main()
