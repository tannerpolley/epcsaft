"""Reacting water/1-butanol/sodium acetate two-liquid state with acetate substantial in both liquids (one tie line)."""
import csv
import json
import math
from collections.abc import Callable
from pathlib import Path
from typing import Any, cast

from epcsaft import Mixture, Parameters, equilibrium as eq

ANALYSIS = Path(__file__).resolve().parents[3] / "engine/tests/cases/sodium-acetate-butanol-reacting-lle"
SPECIES = ("water", "1-butanol", "sodium-cation", "acetate-anion", "acetic-acid", "hydrogen-cation", "hydroxide-anion")
CHARGE = (0, 0, 1, -1, 0, 1, -1)
T, P, LN10, M_W = 298.15, 1.0e5, math.log(10.0), 0.01801528
MASS = {"water": 18.01528, "1-butanol": 74.1216, "salt": 82.03379}


def solve_tie_line_one() -> tuple[Mixture, dict[str, Any], dict[str, float], Any]:
    mix = Mixture(Parameters.from_json(ANALYSIS / "inputs" / "parameters-ternary.json"))
    reactions = json.loads((ANALYSIS / "inputs" / "reactions.json").read_text())["reactions"]
    with (ANALYSIS / "inputs" / "source" / "table-s1-ch3coona.csv").open() as stream:
        row = next(csv.DictReader(line for line in stream if not line.startswith("#")))
    n = [(float(row[f"aq_{c}"]) + float(row[f"org_{c}"])) / 2 / MASS[c] for c in ("water", "1-butanol", "salt")]
    total = sum(n) + n[2]
    feed = dict(zip(SPECIES, (n[0] / total, n[1] / total, n[2] / total, n[2] / total, 0.0, 0.0, 0.0)))
    s = feed["sodium-cation"]
    aq, org = [0.9, 0.01, 1.5 * s, 1.5 * s, 1e-7, 1e-12, 1e-6], [0.3, 0.6, s / 10, s / 10, 1e-7, 1e-12, 1e-8]
    aq[0], org[1] = 1 - sum(aq[1:]), 1 - org[0] - sum(org[2:])
    problem = eq.Problem(
        T=T, P=P, feed=eq.Amounts(feed), row_tolerance=1e-10,
        reactions=[eq.Reaction(r["stoichiometry"], name=name,
                               reference=eq.ReactionReference(eq.ReferenceSourceBasis.CommonMolalityInfiniteDilution),
                               correlation=eq.ReactionLogPolynomial(-r["pK"] * LN10, 0, 0, 0, T, T, T, f"molality-pK-{name}"))
                   for name, r in reactions.items()],
        neutral_reference=eq.NeutralReference(solvent="water", standard_molality_mol_per_kg=1.0,
                                              solute_molality_ratios={c: 1 for c in SPECIES[1:]}),
        phases=[eq.Phase("aqueous", kind="liquid", support=SPECIES, amount=eq.Free(), composition_guess=aq),
                eq.Phase("organic", kind="liquid", support=SPECIES, amount=eq.Free(), composition_guess=org)])
    return mix, reactions, feed, eq.solve_equilibrium(mix, problem)


def test_tie_line_one_conserves_closes_and_holds_acetate_in_both_liquids() -> None:
    """Falsifies wrong reaction stoichiometry or constants, a lost organic charge split, and drift in the pinned point."""
    mix, reactions, feed, r = solve_tie_line_one()
    assert r.success and r.requested_tolerance_met and r.convergence.residual <= 1e-10
    x = [list(r.mole_fractions[:7]), list(r.mole_fractions[7:])]
    n = list(r.amounts)
    state = [mix.state(T, rho=r.molar_densities[p], x=x[p]) for p in range(2)]

    def total(f: Callable[[list[float]], float]) -> float:
        return math.fsum(n[p] * f(x[p]) for p in range(2))
    assert abs(total(lambda v: v[2]) - feed["sodium-cation"]) <= 1e-12
    assert abs(total(lambda v: v[3] + v[4]) - feed["acetate-anion"]) <= 1e-12
    assert abs(total(lambda v: v[1]) - feed["1-butanol"]) <= 1e-12
    assert abs(total(lambda v: v[0] + v[6]) - feed["water"]) <= 1e-12
    assert abs(total(lambda v: v[5] + v[4] - v[6])) <= 1e-12  # proton condition
    for p in range(2):
        assert abs(sum(c * xi for c, xi in zip(CHARGE, x[p]))) / sum(abs(c) * xi for c, xi in zip(CHARGE, x[p])) <= 1e-12
        assert abs(state[p].pressure / P - 1) <= 1e-8
    lna = []
    for p in range(2):
        density, potential = state[p].molar_density, state[p].residual_chemical_potential_over_rt
        assert density is not None and potential is not None
        residual = cast(list[float], potential)
        lna.append([math.log(x[p][i] * density) + residual[i] for i in range(7)])
    d = [lna[1][i] - lna[0][i] for i in range(7)]
    for name, rx in reactions.items():
        nu = rx["stoichiometry"]
        contraction = r.neutral_reference_diagnostics[list(reactions).index(name)].terminal_contraction
        assert contraction is not None
        target = -rx["pK"] * LN10 + (sum(nu.values()) - nu.get("water", 0)) * math.log(M_W) + contraction
        for p in range(2):
            assert abs(math.fsum(q * lna[p][SPECIES.index(c)] for c, q in nu.items()) - target) <= 1e-10
    assert all(abs(d[i]) <= 1e-10 for i in (0, 1, 4))  # neutral phase transfer
    assert abs(d[2] + d[3]) <= 1e-10  # Na+ Ac- mean-ionic transfer
    assert abs(d[5] - d[2]) <= 1e-10 and abs(d[6] - d[3]) <= 1e-10  # shared electric offset
    assert min(x[1][2], x[1][3]) >= 1e-6  # acetate is substantial in the organic liquid
    assert math.isclose(x[1][2], 5.270184253091113e-4, rel_tol=1e-6)  # pinned organic x(Na+), tie line 1


def test_hmf_row_four_without_hmf_acetic_acid_association_meets_row_tolerance() -> None:
    """Table S3 row 4 with HMF/HAc cross association off, feed shares 0.5 and 0.7 (#1094).

    The proton moiety H+ + HAc - OH- sits at the 1e-6 trace size, so a polish trial that crosses it rescaled its row
    by 1e6 and every contracting step read as growth; the solve stalled at 5e-8.
    """
    import importlib.util

    spec = importlib.util.spec_from_file_location("acetate_declaration", ANALYSIS / "scripts" / "run.py")
    assert spec is not None and spec.loader is not None
    run = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(run)
    case = run.CASES["hmf"]
    mix = run.mixture(case, hmf_cross_off=("acetic-acid",))
    for share in (0.5, 0.7):
        r = run.solve(case, run.P.SPECIES["hmf"], mix, 3, share=share)["result"]
        assert r.success and r.requested_tolerance_met and r.convergence.residual <= 1e-10, (share, r.message)
