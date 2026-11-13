"""Acetic acid + water + toluene LLE with acetic acid and water dissociation (#374).

Every equilibrium, density and reference calculation is native; Python declares feeds and reads states.
pH follows main-paper Eqs. 4-5 (Thomsen): a_H = x_H phi_H / (phi_H^inf,w M_w m0), the infinite-dilution
fugacity coefficient taken in pure water at T, P.
"""
import csv
import json
import math
import sys
import time
from pathlib import Path

import epcsaft
from epcsaft import equilibrium as eq

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
T, P = 298.15, 1.0e5
LN10 = math.log(10.0)
M_W = 0.01801528  # kg/mol, as in the parameter documents
S = ("water", "toluene", "acetic-acid", "hydrogen-cation", "acetate-anion", "hydroxide-anion")
N = len(S)
W, TOL, HAC, H, AC, OH = range(N)
Z = (0, 0, 0, 1, -1, -1)
# pKa: ESI Table S4. pKw substitutes for the ESI's Li-Duan f(T,p) (refs. 30-31, not evaluated): the IAPWS R11-07
# formulation (Bandura and Lvov 2006) at 298.15 K and 997.047 kg/m3 (IAPWS-95, 0.1 MPa) gives 13.9945.
PKA, PKW = 4.76, 13.995
REACTIONS = (("HAc", {"acetic-acid": -1, "hydrogen-cation": 1, "acetate-anion": 1}, PKA),
             ("H2O", {"water": -1, "hydrogen-cation": 1, "hydroxide-anion": 1}, PKW))
RATIOS = {"toluene": 1, "acetic-acid": 1, "hydrogen-cation": 2, "acetate-anion": 1, "hydroxide-anion": 1}
FEED_WATER, FEED_TOLUENE = 1.0, 1.0
TIMING = []




def model(variant, species=S):
    return epcsaft.Mixture(epcsaft.Parameters.from_json(ROOT / f"inputs/parameters-{variant}.json", components=species))


def reactions():
    basis = eq.ReactionReference(eq.ReferenceSourceBasis.CommonMolalityInfiniteDilution)
    return [eq.Reaction(nu, name=name, reference=basis,
                        correlation=eq.ReactionLogPolynomial(-pk * LN10, 0, 0, 0, T, T, T, f"molality-pK-{name}"))
            for name, nu, pk in REACTIONS]


def problem(n_acid, ions_in_organic, reactive=True, species=S):
    feed = dict(zip(species, (FEED_WATER, FEED_TOLUENE, n_acid, 0.0, 0.0, 0.0)))
    ion = 1e-3 * math.sqrt(n_acid / 0.01)
    aq = [1 - n_acid - 1e-4 - 2 * ion, 1e-4, n_acid - ion, ion, ion, 1e-11][:len(species)]
    org = [1e-2, 0.99 - n_acid / 10, n_acid / 10, 1e-28, 1e-28, 1e-40][:len(species)]
    org_support = species if ions_in_organic else species[:3]
    if not reactive:
        aq[0] += sum(aq[3:])
        aq, org = aq[:3], org[:3]
    kwargs = dict(reactions=reactions(), neutral_reference=eq.NeutralReference(
        solvent="water", standard_molality_mol_per_kg=1.0,
        solute_molality_ratios=RATIOS,
    )) if reactive else {}
    return eq.Problem(T=T, P=P, feed=eq.Amounts(feed), row_tolerance=1e-10, phases=[
        eq.Phase("aqueous", kind="liquid", support=species if reactive else species[:3], amount=eq.Free(), composition_guess=aq),
        eq.Phase("organic", kind="liquid", support=org_support if reactive else species[:3], amount=eq.Free(),
                 composition_guess=org[:len(org_support)] if reactive else org)], **kwargs)
