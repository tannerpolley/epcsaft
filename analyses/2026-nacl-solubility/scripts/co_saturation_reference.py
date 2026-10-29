"""Convert retained thermochemical tables and explicit salt-content assumptions."""
import json
import math
from pathlib import Path
import numpy as np
from scipy.interpolate import BPoly, CubicSpline
from scipy.integrate import quad

HERE = Path(__file__).resolve().parents[1]
SOURCE = json.loads((HERE / 'inputs/co-saturation.json').read_text())


def thermochemistry(R):
    t = SOURCE['temperature_K']
    na, k, obs = (SOURCE[key] for key in ('NaCl', 'KCl', 'measurement'))
    tn, tr = na['temperature_K'], k['reference_temperature_K']
    thermal = (t - tn) - t * math.log(t / tn)
    gn = na['h_J_mol'] * (1 - t / tn) + na['g_J_mol'] * t / tn + na['cp_J_mol_K'] * thermal
    x, g, ds, cp = (np.array(k[key]) for key in ('aqueous_T_K', 'aqueous_G_increment_J_mol',
        'aqueous_entropy_increment_J_mol_K', 'aqueous_Cp_J_mol_K'))
    sr = k['aqueous_entropy_J_mol_K']
    gaq = BPoly.from_derivatives(x, [[a, -(sr + b), -c / q] for q, a, b, c in zip(x, g, ds, cp)])
    b = k['crystal_b']
    spline = CubicSpline(k['crystal_T_K'], k['crystal_spline_values'], bc_type=((1, -b), (2, 0.)))
    def cpcr(q):
        return (q / (spline(q) + b * q))**3
    gcr = -(t - tr) * k['crystal_entropy_J_mol_K'] + quad(lambda q: cpcr(q) * (1 - t / q), tr, t, epsabs=1e-10)[0]
    gk = k['g_reference_J_mol'] + float(gaq(t)) - gcr
    m, u, rounding, masses = (obs[key] for key in ('molality_mol_kg', 'relative_content_uncertainty',
        'half_last_digit_mol_kg', 'molar_mass_kg_mol'))
    total = 1 + sum(a * b for a, b in zip(m, masses))
    w = [a * b / total for a, b in zip(m, masses)]
    lo, hi = ([a * (1 + sign * b) for a, b in zip(w, u)] for sign in (-1, 1))
    assert abs(cpcr(tr) - 51.44) < .01 and abs(float(gaq(tr))) < 1e-12
    return dict(T_K=t, G_J_mol=[gn, gk], K_source=[math.exp(-v / (R * t)) for v in (gn, gk)],
        G_halfwidth_J_mol=[na['g_se_J_mol'] * t / tn + na['h_se_J_mol'] * abs(1 - t / tn)
            + na['cp_se_J_mol_K'] * abs(thermal) + na['rounding_J_mol'],
            k['g_halfwidth_J_mol'] + k['entropy_halfwidth_J_mol_K'] * abs(t - tr) + k['mapping_allowance_J_mol']],
        aqueous_G_increment_J_mol=float(gaq(t)),
        crystal_G_increment_J_mol=gcr, nominal_mass_fractions=w,
        direct_intervals=[[a * (1 - b) - r, a * (1 + b) + r] for a, b, r in zip(m, u, rounding)],
        mass_fraction_intervals=[[l / (mass * (1 - sum(lo))) - r, h / (mass * (1 - sum(hi))) + r]
            for l, h, mass, r in zip(lo, hi, masses, rounding)])
