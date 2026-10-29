"""Independent FeOs 0.10.1 declaration of Wangler Tables 4, 6–7; one thread."""

import json
import hashlib
from importlib import metadata
import feos
import math
from pathlib import Path
import numpy as np
from feos import (
    BinaryRecord,
    Contributions,
    EquationOfState,
    Identifier,
    Parameters,
    PhaseEquilibrium,
    PureRecord,
    State,
    set_num_threads,
)
from si_units import JOULE, KB, KELVIN, METER, MOL, NAV, PASCAL, RGAS

set_num_threads(1)
OUT = Path(__file__).resolve().parents[4] / "build/issue-857/build-r1"
T = 353.0
sw = 2.7927 + 10.11 * math.exp(-0.01775 * T) - 1.417 * math.exp(-0.01146 * T)
association_correction = 0.600 - 0.00257 * (T - 298.15) - 0.000017 * (T - 298.15) ** 2
k = (
    -0.184
    - 0.00121 * (T - 298.15)
    + 0.0000027 * (T - 298.15) ** 2
    + 0.000000153 * (T - 298.15) ** 3
)
ids = [Identifier(name="water"), Identifier(name="mdea")]
records = [
    PureRecord(
        ids[0],
        18.01528,
        m=1.2046,
        sigma=sw,
        epsilon_k=353.94,
        association_sites=[dict(kappa_ab=0.045, epsilon_k_ab=2425.6, na=1.0, nb=1.0)],
    ),
    PureRecord(
        ids[1],
        119.163,
        m=3.675,
        sigma=3.563,
        epsilon_k=228.71,
        association_sites=[dict(kappa_ab=0.1238, epsilon_k_ab=2046.6, na=2.0, nb=2.0)],
    ),
]
# FeOs association strength multiplies kappa by (sigma_i*sigma_j)^(3/2).
# Engine multiplies kappa by sigma_ij^3: enter the equal physical volume.
kap = math.sqrt(0.045 * 0.1238) * (1 - association_correction)
p = Parameters.from_records(
    records,
    [
        BinaryRecord(
            *ids, k_ij=k, association_sites=[dict(kappa_ab=kap, epsilon_k_ab=2236.1)]
        )
    ],
)
eos = EquationOfState.pcsaft(p, max_iter_cross_assoc=1000, tol_cross_assoc=1e-13)
RHO = MOL / METER**3
terms = {
    "Hard Sphere": "hard_chain",
    "Hard Chain": "hard_chain",
    "Dispersion": "dispersion",
    "Association": "association",
}


def phase(s):
    a = {q: 0.0 for q in set(terms.values())}
    mu = {q: [0.0, 0.0] for q in set(terms.values())}
    for name, v in s.residual_molar_helmholtz_energy_contributions():
        a[terms[name]] += v / (RGAS * s.temperature)
    for j in range(2):
        for name, v in s.chemical_potential_contributions(j, Contributions.Residual):
            mu[terms[name]][j] += v / (RGAS * s.temperature)
    return dict(
        rho=s.density / RHO,
        x=s.molefracs.tolist(),
        P=s.pressure() / PASCAL,
        mu_res_over_RT=(
            s.chemical_potential(Contributions.Residual) / (RGAS * s.temperature)
        ).tolist(),
        lnphi=s.ln_phi().tolist(),
        fugacity_Pa=(
            np.exp(s.ln_phi()) * s.molefracs * (s.pressure() / PASCAL)
        ).tolist(),
        terms=a,
        mu_terms_over_RT=mu,
    )


eng = json.loads(OUT.joinpath("engine-cross-r2.json").read_text())
fixed = phase(
    State(
        eos,
        temperature=T * KELVIN,
        density=eng["phases"][0]["rho"] * RHO,
        composition=np.array([0.8, 0.2]),
    )
)
OUT.joinpath("feos-fixed-r2.json").write_text(json.dumps(fixed, indent=2) + "\n")
print("FIXED", fixed, flush=True)
eq = PhaseEquilibrium.bubble_point(
    eos,
    T * KELVIN,
    np.array([0.8, 0.2]),
    tp_init=36000 * PASCAL,
    vapor_molefracs=np.array([0.9995, 0.0005]),
    tol_inner=1e-9,
    tol_outer=1e-9,
    max_iter_inner=1000,
    max_iter_outer=1000,
)
d = dict(
    code="FeOs",
    T=T,
    k=k,
    association_correction=association_correction,
    sigma_water=sw,
    physical_cross_volume_A3=kap * (sw * 3.563) ** 1.5,
    parameters=[json.loads(v) if v else None for v in p.to_json_str(pretty=False)],
    constants=dict(
        NA=NAV * MOL, kB=KB / (JOULE / KELVIN), R=RGAS / (JOULE / MOL / KELVIN)
    ),
    phases=[phase(eq.liquid), phase(eq.vapor)],
)
e = json.loads(OUT.joinpath("engine-cross-r2.json").read_text())
d["at_engine_liquid"] = phase(
    State(
        eos,
        temperature=T * KELVIN,
        density=e["phases"][0]["rho"] * RHO,
        composition=np.array([0.8, 0.2]),
    )
)
OUT.joinpath("feos-cross-r2.json").write_text(json.dumps(d, indent=2) + "\n")
print(json.dumps(d, indent=2))

ratio = math.sqrt(sw * 3.563) / ((sw + 3.563) / 2)
k_engine = math.sqrt(0.045 * 0.1238 * ratio**3 * (1 - association_correction))
p2 = Parameters.from_records(
    records,
    [
        BinaryRecord(
            *ids,
            k_ij=k,
            association_sites=[dict(kappa_ab=k_engine / ratio**3, epsilon_k_ab=2236.1)],
        )
    ],
)
eos2 = EquationOfState.pcsaft(p2, max_iter_cross_assoc=1000, tol_cross_assoc=1e-13)
r = PhaseEquilibrium.bubble_point(
    eos2,
    353 * KELVIN,
    np.array([0.8, 0.2]),
    tp_init=30000 * PASCAL,
    vapor_molefracs=np.array([0.9995, 0.0005]),
    tol_inner=1e-9,
    tol_outer=1e-9,
    max_iter_inner=1000,
    max_iter_outer=1000,
)
literal = dict(
    parameters=[json.loads(z) if z else None for z in p2.to_json_str(pretty=False)],
    phases=[phase(r.liquid), phase(r.vapor)],
)
OUT.joinpath("feos-literal-eq12-r2.json").write_text(
    json.dumps(literal, indent=2) + "\n"
)
pure = []
eosp = EquationOfState.pcsaft(
    Parameters.from_records([records[1]]),
    max_iter_cross_assoc=1000,
    tol_cross_assoc=1e-13,
)
for row in json.loads(OUT.joinpath("pure-mdea-r2.json").read_text()):
    t = row["T"]
    try:
        eq = PhaseEquilibrium.pure(eosp, t * KELVIN, tol=1e-9, max_iter=1000)
        q = dict(
            T=t,
            P=eq.vapor.pressure() / PASCAL,
            rho=eq.liquid.density / RHO,
            rho_mass=(eq.liquid.density / RHO) * 0.119163,
            mu_liquid=(
                eq.liquid.chemical_potential(Contributions.Residual)
                / (RGAS * eq.liquid.temperature)
            ).tolist(),
            lnphi_liquid=eq.liquid.ln_phi().tolist(),
            lnphi_vapor=eq.vapor.ln_phi().tolist(),
            P_liquid=eq.liquid.pressure() / PASCAL,
        )
    except Exception as ex:
        q = dict(T=t, error=str(ex))
    pure.append(q)
    print("PURE", q, flush=True)
OUT.joinpath("feos-pure-r2.json").write_text(json.dumps(pure, indent=2) + "\n")

# Preserve the independent comparator identity with its native executable hash.
package = Path(feos.__file__).parent
OUT.joinpath("feos-identity.json").write_text(
    json.dumps(
        {
            "version": metadata.version("feos"),
            "native_libraries": {
                str(p): hashlib.sha256(p.read_bytes()).hexdigest()
                for p in package.glob("*.so")
            },
            "threads": 1,
            "constants": d["constants"],
        },
        indent=2,
    )
    + "\n"
)
