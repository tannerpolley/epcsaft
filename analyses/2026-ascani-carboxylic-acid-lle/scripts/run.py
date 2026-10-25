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


def write(name, rows):
    with (RESULTS / name).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(dict.fromkeys(k for row in rows for k in row)), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def step(name, function, *args):
    start = time.perf_counter()
    result = function(*args)
    TIMING.append(dict(step=name, wall_s=round(time.perf_counter() - start, 3)))
    return result


def phase_states(mix, result, k):
    out = []
    for p in range(2):
        x = result.mole_fractions[p * k:(p + 1) * k]
        out.append((x, mix.state(T, rho=result.molar_densities[p], x=x)))
    return out


def baseline():
    """Check 1: water + toluene and non-reactive acid distribution with the ESI k_ij."""
    rows = []
    binary = model("literal", S[:2])
    feed = {"water": 1.0, "toluene": 1.0}
    r = eq.solve_equilibrium(binary, eq.Problem(T=T, P=P, feed=eq.Amounts(feed), row_tolerance=1e-10, phases=[
        eq.Phase("aqueous", kind="liquid", amount=eq.Free(), composition_guess=[0.9999, 1e-4]),
        eq.Phase("organic", kind="liquid", amount=eq.Free(), composition_guess=[3e-3, 0.997])]))
    rows.append(dict(case="water-toluene", success=r.success, tolerance_met=r.requested_tolerance_met,
                     iterations=r.iterations, x_toluene_aq=r.mole_fractions[1], x_water_org=r.mole_fractions[2],
                     reduced_hessian_minimum=r.reduced_hessian_minimum))
    ternary = model("literal", S[:3])
    start = None
    for n in grid(41):
        r = eq.solve_equilibrium(ternary, problem(n, False, reactive=False, species=S[:3]), start=start)
        start = r if r.success else None
        rows.append(dict(case="nonreactive", n_acid_feed=n, success=r.success, tolerance_met=r.requested_tolerance_met,
                         iterations=r.iterations, x_acid_aq=r.mole_fractions[2], x_acid_org=r.mole_fractions[5],
                         K_aq_org=r.mole_fractions[2] / r.mole_fractions[5], x_toluene_aq=r.mole_fractions[1],
                         x_water_org=r.mole_fractions[3], reduced_hessian_minimum=r.reduced_hessian_minimum))
    write("baseline.csv", rows)


def grid(points, low=0.002, high=0.05):
    return [low * (high / low) ** (i / (points - 1)) for i in range(points)]


def ln_activity_h(x, state, ln_phi_inf_h):
    # Main-paper Eqs. 4-5: ln a_H = ln(x_H phi_H / (phi_H^inf,w M_w m0)), m0 = 1 mol/kg
    return math.log(x[H]) + state.log_fugacity_coefficient[H] - ln_phi_inf_h - math.log(M_W)


def sweep(case, variant, ions_in_organic, points):
    mix = model(variant)
    reference = mix.reference(T, P, solvent_fractions=[1, 0, 0, 0, 0, 0], ions=("hydrogen-cation", "acetate-anion"),
                              stoichiometry=(1, 1))
    ln_phi_inf_h, ln_phi_inf_ac = reference.log_fugacity
    water_hac = mix.state(T, P=P, x=[1 - 1e-12, 0, 1e-12, 0, 0, 0], phase="liquid")
    ln_phi_inf_hac = water_hac.log_fugacity_coefficient[HAC]
    # Density-basis targets ln K_rho = ln K_m + A_r + C_r (EqID neutral_reaction_reference_transfer). C_r is also
    # recomputed as the Engine's sqrt(m) limit (10 f(1e-14) - f(1e-12))/9 (EqID reference_limit_terminal_pair)
    # from public states at 1e-12 and 1e-14 mol/kg in water.
    def potentials(m):
        c = [1 / M_W] + [m * RATIOS[s] for s in S[1:]]
        st = mix.state(T, P=P, x=[v / sum(c) for v in c], phase="liquid")
        return [math.log(st.molar_density) + mu for mu in st.residual_chemical_potential_over_rt]
    lam = [(10 * fine - middle) / 9 for middle, fine in zip(potentials(1e-12), potentials(1e-14))]
    rows, start = [], None
    for n in grid(points):
        request = problem(n, ions_in_organic)
        t0 = time.perf_counter()
        r = eq.solve_equilibrium(mix, request, start=start)
        solve_ms = 1000 * (time.perf_counter() - t0)
        row = dict(case=case, n_acid_feed=n, success=r.success, tolerance_met=r.requested_tolerance_met,
                   solver=r.solver_status, iterations=r.iterations, evaluations=r.evaluations,
                   hessian_calls=r.hessian_calls, solve_ms=round(solve_ms, 2), classified_residual=r.convergence.residual,
                   reduced_hessian_minimum=r.reduced_hessian_minimum,
                   rank_deficit=r.unknown_count - r.jacobian_rank, topology_event=r.topology_event.name,
                   reference_terminal_change=max((d.observed_terminal_change for d in r.neutral_reference_diagnostics if d.observed_terminal_change is not None), default=None),
                   ln_phi_inf_h_water=ln_phi_inf_h, reference_convergence_error=reference.convergence_error)
        rows.append(row)
        if not r.success:
            continue
        start = r
        (xa, sa), (xo, so) = phase_states(mix, r, N)
        for name, x, s, amount in (("aq", xa, sa, r.amounts[0]), ("org", xo, so, r.amounts[1])):
            row[f"P_rel_{name}"] = s.pressure / P - 1
            for i in (HAC, H, AC, OH):
                row[f"x_{S[i]}_{name}"] = x[i]
                row[f"n_{S[i]}_{name}_mol"] = amount * x[i]
            row[f"charge_{name}_mol"] = amount * math.fsum(z * v for z, v in zip(Z, x))
            ions = math.fsum(abs(z) * v for z, v in zip(Z, x))
            row[f"charge_relative_{name}"] = math.fsum(z * v for z, v in zip(Z, x)) / ions if ions > 0 else 0.0
        row["x_acid_total_aq"] = xa[HAC] + xa[AC]
        row["K_aq_org_molecular"] = xa[HAC] / xo[HAC]
        row["K_aq_org_total"] = (xa[HAC] + xa[AC]) / (xo[HAC] + xo[AC])
        row["dissociation_degree_aq"] = xa[AC] / (xa[HAC] + xa[AC])
        row["molality_h_aq"] = xa[H] / (xa[W] * M_W)
        row["ln_gamma_m_h_aq"] = sa.log_fugacity_coefficient[H] - ln_phi_inf_h + math.log(xa[W])
        row["pH_aq"] = -ln_activity_h(xa, sa, ln_phi_inf_h) / LN10
        for k, (name, nu, pk) in enumerate(REACTIONS):
            d = r.neutral_reference_diagnostics[k]
            assert d.reaction_index == k
            target = -pk * LN10 + (sum(nu.values()) - nu.get("water", 0)) * math.log(M_W) + d.terminal_contraction
            row[f"contraction_recomputed_minus_compiled_{name}"] = (
                math.fsum(v * lam[S.index(s)] for s, v in nu.items()) - d.terminal_contraction)
            for name_p, x, s in (("aq", xa, sa), ("org", xo, so)):
                if all(x[S.index(sp)] > 0 for sp in nu):
                    row[f"affinity_density_{name}_{name_p}"] = math.fsum(
                        v * (math.log(x[S.index(sp)] * s.molar_density) + s.residual_chemical_potential_over_rt[S.index(sp)])
                        for sp, v in nu.items()) - target
        # Reference-consistency diagnostic, not the reaction check: the acid law on molality activities with phi_inf
        # from a separate public reference path (ion pair and HAc each dilute in pure water).
        ln_a = lambda i, ln_inf: math.log(xa[i]) + sa.log_fugacity_coefficient[i] - ln_inf - math.log(M_W)
        row["acid_law_residual_molality"] = (ln_a(H, ln_phi_inf_h) + ln_a(AC, ln_phi_inf_ac) - ln_a(HAC, ln_phi_inf_hac)) + PKA * LN10
        if ions_in_organic:
            row["pH_org_single_ion"] = -ln_activity_h(xo, so, ln_phi_inf_h) / LN10
            d = [math.log(xo[i] / xa[i]) + so.log_fugacity_coefficient[i] - sa.log_fugacity_coefficient[i] for i in (H, AC)]
            row["mean_ionic_transfer_residual"] = d[0] + d[1]  # ln f_H f_Ac equal across liquids
            row["electrochemical_offset"] = (d[0] - d[1]) / 2
            # Reconstructed model estimate, not a returned value: the absolute charge row leaves the trace organic
            # cation/anion split free (#375). Imposing electroneutrality on the returned product gives its square root;
            # pH_org_reconstructed pairs it with phi_H of the returned, charge-unresolved state.
            pair = math.sqrt(xo[H] * xo[AC])
            row["x_ion_org_reconstructed"] = pair
            row["n_ion_org_reconstructed_mol"] = pair * r.amounts[1]
            row["pH_org_reconstructed"] = -(math.log(pair) + so.log_fugacity_coefficient[H] - ln_phi_inf_h - math.log(M_W)) / LN10
    write(f"sweep-{case}.csv", rows)
    return rows


def interpolate(rows, key, x):
    # ESI x_acid is the analytical acid mole fraction, HAc + Ac-.
    pts = sorted((math.log(r["x_acid_total_aq"]), r[key]) for r in rows if r.get("pH_aq") is not None)
    lx = math.log(x)
    for (a, fa), (b, fb) in zip(pts, pts[1:]):
        if a <= lx <= b:
            return fa + (fb - fa) * (lx - a) / (b - a)
    raise ValueError(f"x = {x} outside the sweep")


def compare(sweeps):
    esi_k = list(csv.DictReader((ROOT / "inputs/esi-table-11-partition.csv").open()))
    esi_ph = list(csv.DictReader((ROOT / "inputs/esi-table-12-ph.csv").open()))
    out = []
    for k_row, ph_row in zip(esi_k, esi_ph, strict=True):
        x = float(ph_row["x_acetic_acid_aq"])
        row = dict(x_acetic_acid_aq=x, pH_exp=float(ph_row["pH_exp"]), pH_esi_model=float(ph_row["pH_model"]),
                   K_exp=float(k_row["K_aq_org_exp"]), K_esi_model=float(k_row["K_aq_org_model"]))
        for case, (full, coarse) in sweeps.items():
            ph = interpolate(full, "pH_aq", x)
            k = interpolate(full, "K_aq_org_total", x)
            row[f"pH_{case}"] = ph
            row[f"pH_{case}_minus_esi"] = ph - row["pH_esi_model"]
            row[f"pH_{case}_coarse_fine_difference"] = abs(ph - interpolate(coarse, "pH_aq", x))
            row[f"K_{case}"] = k
            row[f"K_{case}_rel_esi"] = k / row["K_esi_model"] - 1
            row[f"K_{case}_coarse_fine_difference"] = abs(k / interpolate(coarse, "K_aq_org_total", x) - 1)
        out.append(row)
    write("comparison.csv", out)


def main(points):
    RESULTS.mkdir(exist_ok=True)
    started = time.strftime("%Y-%m-%dT%H:%M:%S%z")
    step("baseline", baseline)
    cases = (("aq-literal", "literal", False), ("aq-constant-sigma", "constant-sigma", False),
             ("aq-cr1-cross", "cr1-cross", False),
             ("aq-organic-eps", "organic-eps", False), ("both-organic-eps", "organic-eps", True))
    sweeps = {}
    for case, variant, both in cases:
        full = step(f"sweep {case} ({points} feeds)", sweep, case, variant, both, points)
        coarse = [r for i, r in enumerate(full) if i % 2 == 0]
        sweeps[case] = (full, coarse)
    step("compare", compare, sweeps)
    write("timing.csv", TIMING)
    (RESULTS / "provenance.json").write_text(json.dumps(dict(
        engine=epcsaft.provenance(), started=started, feeds=points, T_K=T, P_Pa=P, pKa=PKA, pKw=PKW,
        feed_water_mol=FEED_WATER, feed_toluene_mol=FEED_TOLUENE, esi_sha256="84c3c8da4592a25f2eef1bddc5c3712096c31463c5df88e4b4fb6f1df6e9ede1",
    ), indent=1) + "\n")


def declarations():
    return [
        (model(v), problem(0.01, both), problem(0.01, False, reactive=False, species=S[:3]))
        for v in ('literal', 'constant-sigma', 'cr1-cross', 'organic-eps')
        for both in (False, True)
    ]


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 81)
