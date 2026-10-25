"""FeOs 0.10.1 comparator rows for the Gross (2005) quadrupole term (#574), only at states where every m_i <= 2:
FeOs caps each m_i at 2 (the dipole convention), so benzene is excluded. Its J constants carry more digits than
Table 2, hence the 1e-6 criterion. Moments are rescaled to the Engine's debye bridge as in clapeyron_reference.jl.

    ~/.cache/epcsaft-comparators/222/feos/bin/python analyses/2026-polar-quadrupole/scripts/feos_reference.py
"""
import math
from pathlib import Path

import feos
import numpy as np
import si_units as si

debye = 1e-21 / 299792458.0
scale = math.sqrt(debye**2 / (4 * math.pi * 8.8541878128e-12 * 1e-30) / 1e-19)
CO2 = ("carbon-dioxide", 44.01, 1.5131, 3.1869, 163.33, 4.4)
N2 = ("nitrogen", 28.01, 1.1504, 3.3848, 91.40, 1.43)


def residual(components, x, T, rho, quadrupole):
    records = [feos.PureRecord(feos.Identifier(name=c[0]), c[1], m=c[2], sigma=c[3], epsilon_k=c[4],
                               **({"q": scale * c[5]} if quadrupole else {})) for c in components]
    parameters = feos.Parameters.new_pure(records[0]) if len(records) == 1 else feos.Parameters.new_binary(records)
    state = feos.State(feos.EquationOfState.pcsaft(parameters), T * si.KELVIN, density=rho * si.MOL / si.METER**3,
                       composition=np.array(x))
    return state.molar_helmholtz_energy(feos.Contributions.Residual) / (si.RGAS * T * si.KELVIN)


with (Path(__file__).resolve().parents[1] / "results/feos-term-rows.csv").open("w") as out:
    out.write("species,x,T_K,rho_mol_m3,a_qq_feos\n")
    for components, x, T, rho in (([CO2], [1.0], 250.0, 1000.0), ([CO2], [1.0], 250.0, 22000.0),
                                  ([CO2, N2], [0.7, 0.3], 240.0, 15000.0)):
        value = residual(components, x, T, rho, True) - residual(components, x, T, rho, False)
        out.write(f"{'+'.join(c[0] for c in components)},{' '.join(map(str, x))},{T:.17g},{rho:.17g},{value:.17g}\n")
