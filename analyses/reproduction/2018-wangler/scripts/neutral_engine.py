import json
import math
import csv
import sys
import time
import subprocess
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from neutral_models import neutral_mapping, ROOT, LEAF
import epcsaft
from epcsaft import equilibrium as e

OUT = ROOT / "build/issue-857/build-r1"
OUT.mkdir(parents=True, exist_ok=True)
TERMS = ("hard_chain", "dispersion", "association")


def equilibrium(t, x, mode="baseline", pure=False):
    selected = ("mdea",) if pure else ("water", "mdea")
    mapping = neutral_mapping(t, mode, selected)
    model = (epcsaft.Pure if pure else epcsaft.Mixture)(
        epcsaft.Parameters.from_mapping(mapping)
    )
    xx = [1.0] if pure else [1 - x, x]
    feed = dict(zip(selected, xx))
    pguess = (
        math.exp(92.624 - 10367 / t - 9.4699 * math.log(t) + 1.9e-18 * t**6)
        if pure
        else {313: 0.06, 333: 0.15, 353: 0.36, 373: 0.8}.get(t, 0.36) * 1e5
    )
    prob = e.Problem(
        T=t,
        P=e.Free(pguess),
        feed=e.MoleFractions(feed),
        phases=[
            e.Phase("liquid", kind="liquid", support=selected, composition_guess=xx),
            e.Phase(
                "vapor",
                kind="vapor",
                support=selected,
                amount=e.Pinned(0.0),
                composition_guess=[1.0] if pure else [0.9999, 0.0001],
            ),
        ],
        row_tolerance=1e-10,
    )
    tic = time.monotonic()
    r = e.solve_equilibrium(model, prob)
    row = dict(
        code="Engine",
        T=t,
        x_mdea=x,
        mode=mode,
        pure=pure,
        success=r.success,
        message=r.message,
        seconds=time.monotonic() - tic,
    )
    if not r.success:
        raise RuntimeError(row)
    n = len(selected)
    ph = []
    for i, label in enumerate(("liquid", "vapor")):
        rho = r.molar_densities[i]
        xq = list(r.mole_fractions[i * n : (i + 1) * n])
        s = model.state(t, rho=rho, x=xq)
        f = [
            math.exp(v) * rho * epcsaft.GAS_CONSTANT_J_PER_MOL_K * t * xq[j]
            for j, v in enumerate(s.residual_chemical_potential_over_rt)
        ]
        ph.append(
            dict(
                phase=label,
                rho=rho,
                x=xq,
                P=s.pressure,
                mu_res_over_RT=s.residual_chemical_potential_over_rt,
                lnphi=s.log_fugacity_coefficient,
                fugacity_Pa=f,
                terms={k: getattr(s, k) for k in TERMS},
                site_fractions=s.site_fractions,
            )
        )
    row.update(P=r.pressure, phases=ph, residual=max(map(abs, r.residuals)))
    return row, model


def derivative_verification(model, t, rho, x):
    # Verification only: differentiate returned native term energies against FeOs exact derivatives.
    rows = []
    for h in (0.001, 0.0003, 0.0001, 0.00003):
        mu = {k: [] for k in TERMS}
        for i in range(len(x)):
            vals = []
            for sign in (-2, -1, 1, 2):
                ns = list(x)
                ns[i] += sign * h
                nt = sum(ns)
                s = model.state(t, rho=rho * nt, x=[v / nt for v in ns])
                vals.append({k: nt * getattr(s, k) for k in TERMS})
            for k in TERMS:
                mu[k].append(
                    (vals[0][k] - 8 * vals[1][k] + 8 * vals[2][k] - vals[3][k])
                    / (12 * h)
                )
        rows.append(dict(h_mol=h, mu_terms_over_RT=mu))
    return rows


if __name__ == "__main__":
    identity = subprocess.check_output(
        [sys.executable, str(ROOT / "tools/engine_identity.py")], text=True
    )
    OUT.joinpath("engine-identity.json").write_text(identity)
    cross, model = equilibrium(353.0, 0.2)
    cross["term_derivative_verification"] = derivative_verification(
        model, 353.0, cross["phases"][0]["rho"], [0.8, 0.2]
    )
    OUT.joinpath("engine-cross-r2.json").write_text(json.dumps(cross, indent=2) + "\n")
    print("CROSS", cross, flush=True)
    pure = []
    for t in [303.15, 313.0, 333.0, 353.0, 373.0, 393.15, 413.15, 443.15]:
        row, _ = equilibrium(t, 1.0, pure=True)
        row["pahlav_Table1_as_printed_P_Pa"] = math.exp(
            92.624 - 10367 / t - 9.4699 * math.log(t) + 1.9e-18 * t**6
        )
        row["pahlav_Table1_as_printed_rho_mol_m3"] = (
            1000 * 1.0011 / 0.22523 ** (1 + (1 - t / 678.2) ** 0.21515)
        )
        pure.append(row)
        print("PURE", t, row.get("P"), flush=True)
    OUT.joinpath("pure-mdea-r2.json").write_text(json.dumps(pure, indent=2) + "\n")
    sweep = []
    modes = [
        "baseline",
        "eq12-literal",
        "volume-geometric",
        "volume-l-under-root",
        "energy-geometric-mayer",
        "literal-k",
        "literal-l",
        "raw-T-k",
        "k-constant",
        "quadratic-k",
        "cubic-sign-k",
        "water-sigma-constant",
        "mdea-1:1",
        "mdea-4:2",
        "mdea-2:3",
        "induced-retain-cross",
        "induced-zero-self-energy",
        "water-pahlav-2",
        "water-pahlav-4",
        "water-gross-2",
        "mdea-origin-digits",
        "mdea-source-6site",
        "mdea-source-2site",
        "association-diameter-d",
    ]
    for mode in modes:
        for t in (313.0, 333.0, 353.0, 373.0):
            try:
                row, _ = equilibrium(t, 0.2, mode)
            except Exception as ex:
                row = dict(code="Engine", T=t, x_mdea=0.2, mode=mode, error=str(ex))
            sweep.append(row)
            print(
                "SWEEP",
                mode,
                t,
                row.get("P"),
                row.get("error", row.get("message")),
                flush=True,
            )
            OUT.joinpath("convention-sweep-r2.json").write_text(
                json.dumps(sweep, indent=2) + "\n"
            )

    observations = []
    for source in csv.DictReader(
        (LEAF / "inputs/source/kim-2008-neat-mdea.csv").open()
    ):
        t = float(source["T_K"])
        row, _ = equilibrium(t, 1.0, pure=True)
        row["observed_P_Pa"] = float(source["observed_P_Pa"])
        observations.append(row)
    OUT.joinpath("pure-observations-r2.json").write_text(
        json.dumps(observations, indent=2) + "\n"
    )
