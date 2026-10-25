"""Independent teqp reference rows for the Gross-Vrabec dipolar term (#28 falsifier 2).

Run in a separate environment holding teqp 0.23.2 (it is not an Engine dependency):
    python -m venv /tmp/teqp && /tmp/teqp/bin/pip install teqp==0.23.2
    /tmp/teqp/bin/python analyses/2026-polar-dipole/scripts/teqp_reference.py
Each row is alpha_r(polar) - alpha_r(nonpolar) at fixed (T, rho, x), where teqp
receives mu*^2 = mu^2/(4 pi eps0 k m eps sigma^3) computed with the Engine's SI constants.
"""
import csv
import math
from pathlib import Path

import numpy
import teqp

HERE = Path(__file__).resolve().parents[1]
DEBYE, EPS0, K_B = 1.0e-21 / 299792458.0, 8.8541878128e-12, 1.380649e-23
# GV2006 Table 2 (acetone) and Gross-Sadowski 2001 Table 1 (n-pentane).
SPECIES = {
    "acetone": {"m": 2.7447, "sigma": 3.2742, "eps": 232.99, "mu": 2.88},
    "n-pentane": {"m": 2.6896, "sigma": 3.7729, "eps": 231.20, "mu": 0.0},
    "dimethyl-ether": {"m": 2.2634, "sigma": 3.2723, "eps": 210.29, "mu": 1.3},
}
CASES = (  # names, x, T/K, molar densities/(mol/m3): liquid-like, vapor-like and dilute
    (("acetone",), (1.0,), 300.0, (13300.0, 40.0, 1.0e-3)),
    (("acetone",), (1.0,), 450.0, (10000.0, 800.0)),
    (("acetone", "n-pentane"), (0.4, 0.6), 372.7, (9000.0, 300.0)),
    (("dimethyl-ether", "acetone", "n-pentane"), (0.2, 0.3, 0.5), 330.0, (10000.0,)),
)


def model(names: tuple[str, ...], polar: bool) -> object:
    coeffs = []
    for name in names:
        s = SPECIES[name]
        record = {"name": name, "m": s["m"], "sigma_Angstrom": s["sigma"], "epsilon_over_k": s["eps"], "BibTeXKey": ""}
        if polar:
            volume = (s["sigma"] * 1.0e-10) ** 3
            record["(mu^*)^2"] = (s["mu"] * DEBYE) ** 2 / (4.0 * math.pi * EPS0 * K_B * s["eps"] * volume * s["m"])
            record["nmu"] = 1.0 if s["mu"] > 0.0 else 0.0
        coeffs.append(record)
    return teqp.make_model({"kind": "PCSAFT", "model": {"coeffs": coeffs}})


def main() -> None:
    rows = []
    for names, x, temperature, densities in CASES:
        polar, nonpolar = model(names, True), model(names, False)
        z = numpy.array(x)
        for rho in densities:
            delta = polar.get_Ar00(temperature, rho, z) - nonpolar.get_Ar00(temperature, rho, z)
            rows.append({"species": "+".join(names), "x": " ".join(repr(v) for v in x), "T_K": repr(temperature),
                         "rho_mol_m3": repr(rho), "a_dipolar_teqp": repr(delta)})
    out = HERE / "results/teqp-dipolar-rows.csv"
    out.parent.mkdir(exist_ok=True)
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"teqp {teqp.__version__}: {len(rows)} rows -> {out}")


if __name__ == "__main__":
    main()
