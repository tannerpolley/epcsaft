"""Reacting water + 1-butanol + sodium acetate (+ HMF) two-liquid equilibrium at 298.15 K and 1 bar, all seven (eight) species
permitted in both liquids, with acetic-acid dissociation and water autoprotolysis.

Python declares inputs and reads results; every solve, stability search and certificate is native. Parameters come from
scripts/parameters.py (written to inputs/parameters-*.json); tie lines and standard uncertainties from inputs/source/.

Case 1 ("ternary"): the four Table S1 CH3COONa tie lines, calibration data. Case 2 ("hmf"): the seven Table S3 CH3COONa rows,
a conditional prediction. Feed: the mass midpoint of each measured pair.

Writes results/*.csv, results/summary.json and results/provenance.json (unrounded values).
"""
import csv
import hashlib
import json
import math
import statistics
import sys
import time
from pathlib import Path

import epcsaft
from epcsaft import equilibrium as eq

sys.path.insert(0, str(Path(__file__).resolve().parent))
import parameters as P  # noqa: E402

HERE = Path(__file__).resolve().parents[1]
INPUTS, RESULTS = HERE / "inputs", HERE / "results"
T, PRESSURE = P.T, P.PRESSURE
LN10, M_W = math.log(10.0), 0.01801528
E_B = 1e-10  # numerical uncertainty of the Engine's mass fractions: the row tolerance
BOUNDS = {"row_residual": 1e-10, "cold_vs_continued_x": 1e-10, "conservation_mol": 1e-12, "charge_relative": 1e-12,
          "affinity": 1e-10, "phase_transfer": 1e-10, "pressure_relative": 1e-8, "coupling_w": 1e-5, "pair_sensitivity_w": 1e-5}
ORGANIC_ION_FLOOR = 1e-6  # organic x(Na+), x(Ac-) at least this, 4 orders above the row tolerance
SHARES = (0.15, 0.3, 0.5, 0.7, 0.85)  # aqueous mass share of the feed along a measured tie line
REACTION_PK = json.loads((INPUTS / "reactions.json").read_text())["reactions"]
SOURCE = json.loads((INPUTS / "source" / "source-values.json").read_text())


def table(name):
    with (INPUTS / "source" / name).open() as stream:
        return list(csv.DictReader(line for line in stream if not line.startswith("#")))


def cases():
    s1, s3 = table("table-s1-ch3coona.csv"), table("table-s3-ch3coona.csv")
    u = SOURCE["standard_uncertainty_w"]

    def case(system, names, members, rows, u_table, control):
        ties = [([float(r[f"aq_{n}"]) for n in names], [float(r[f"org_{n}"]) for n in names]) for r in rows]
        return dict(system=system, names=names, members=members, rows=rows, ties=ties, control=control,
                    u=[u_table["aq"][n] for n in names] + [u_table["org"][n] for n in names])
    return {
        "ternary": case("ternary", ("water", "1-butanol", "salt"), ("water", "1-butanol", "sodium-cation"), s1, u["table_s1"],
                        P.SPECIES["ternary"][:4]),
        "hmf": case("hmf", ("water", "1-butanol", "hmf", "salt"), ("water", "1-butanol", "hmf", "sodium-cation"), s3, u["table_s3"],
                    P.SPECIES["hmf"][:5])}


CASES = cases()
M_SALT = P.COMPONENTS["sodium-cation"][1] + P.COMPONENTS["acetate-anion"][1]


def molar_mass(component):
    return M_SALT if component == "sodium-cation" else P.COMPONENTS[component][1]


def mixture(case, species=None, **document_options):
    """Mixture of the retained parameter file (default), or of a sensitivity document built by parameters.py."""
    system = case["system"]
    if not document_options:
        return epcsaft.Mixture(epcsaft.Parameters.from_json(INPUTS / f"parameters-{system}.json", components=species))
    return epcsaft.Mixture(epcsaft.Parameters.from_mapping(P.document(system, species, **document_options)))


def feed(case, tie, share=0.5):
    """Mole feed (sums to 1) from the mass-weighted mix of the aqueous and organic compositions; Na+ and Ac- equal."""
    aq, org = tie
    n = [(share * a + (1 - share) * o) / molar_mass(c) for a, o, c in zip(aq, org, case["members"])]
    total = sum(n) + n[-1]
    f = {c: 0.0 for c in P.SPECIES[case["system"]]}
    for c, v in zip(case["members"], n):
        f[c] = v / total
    f["acetate-anion"] = f["sodium-cation"]
    return f


def problem(species, f, reacting):
    s, h = f["sodium-cation"], f.get("hmf", 0.0)
    aq = {"water": 0.9, "1-butanol": 0.01, "hmf": h, "sodium-cation": 1.5 * s, "acetate-anion": 1.5 * s,
          "acetic-acid": 1e-7, "hydrogen-cation": 1e-12, "hydroxide-anion": 1e-6}
    org = {"water": 0.3, "1-butanol": 0.6, "hmf": 1.5 * h, "sodium-cation": s / 10, "acetate-anion": s / 10,
           "acetic-acid": 1e-7, "hydrogen-cation": 1e-12, "hydroxide-anion": 1e-8}
    ga, go = [aq[c] for c in species], [org[c] for c in species]
    ga[0] = 1 - sum(ga[1:])
    go[1] = 1 - go[0] - sum(go[2:])
    options = {}
    if reacting:
        options = dict(
            reactions=[eq.Reaction(r["stoichiometry"], name=name,
                                   reference=eq.ReactionReference(eq.ReferenceSourceBasis.CommonMolalityInfiniteDilution),
                                   correlation=eq.ReactionLogPolynomial(-r["pK"] * LN10, 0, 0, 0, T, T, T, f"molality-pK-{name}"))
                       for name, r in REACTION_PK.items()],
            neutral_reference=eq.NeutralReference(solvent="water", standard_molality_mol_per_kg=1.0,
                                                  solute_molality_ratios={x: 1 for x in species[1:]}))
    return eq.Problem(T=T, P=PRESSURE, feed=eq.Amounts({c: f[c] for c in species}), row_tolerance=BOUNDS["row_residual"], **options,
                      phases=[eq.Phase("aqueous", kind="liquid", support=species, amount=eq.Free(), composition_guess=ga),
                              eq.Phase("organic", kind="liquid", support=species, amount=eq.Free(), composition_guess=go)])


def declarations():
    """Parameter documents, mixtures and one problem per case; nothing is solved."""
    out = {}
    for name, case in CASES.items():
        species = P.SPECIES[case["system"]]
        out[name] = dict(mixture=mixture(case), problem=problem(species, feed(case, case["ties"][0]), True),
                         control=problem(case["control"], feed(case, case["ties"][0]), False))
    return out


def mass_fractions(case, x, species):
    m = [x[species.index(c)] * molar_mass(c) for c in case["members"]]
    return [q / sum(m) for q in m]


def solve(case, species, mix, tie_index, share=0.5, reacting=True, start=None):
    f = feed(case, case["ties"][tie_index], share)
    f = {c: f[c] for c in species}
    t0 = time.perf_counter()
    r = eq.solve_equilibrium(mix, problem(species, f, reacting), start=start)
    return dict(result=r, feed=f, ms=1000 * (time.perf_counter() - t0), mix=mix, species=species, reacting=reacting, case=case,
                ok=bool(r.success and r.requested_tolerance_met))


def read(s):
    """Falsifier quantities of one solved state, from public state potentials."""
    r, mix, species, case = s["result"], s["mix"], s["species"], s["case"]
    k, ix = len(species), {c: i for i, c in enumerate(species)}
    z = {c: P.COMPONENTS[c][0] for c in species}
    row = dict(success=r.success, tolerance_met=r.requested_tolerance_met, solver_status=r.solver_status, iterations=r.iterations,
               ms=s["ms"], residual=r.convergence.residual, reduced_hessian_minimum=r.reduced_hessian_minimum,
               rank_deficit=r.unknown_count - r.jacobian_rank, topology_event=r.topology_event.name)
    if not r.success:
        row["message"] = str(r.message)[:300]
        return row
    xs = [list(r.mole_fractions[:k]), list(r.mole_fractions[k:])]
    st = [mix.state(T, rho=r.molar_densities[p], x=xs[p]) for p in range(2)]
    amt, f = list(r.amounts), s["feed"]
    g = lambda x, c: x[ix[c]] if c in ix else 0.0  # noqa: E731
    total = lambda fn: math.fsum(amt[p] * fn(xs[p]) for p in range(2))  # noqa: E731
    row["amount_aq"], row["amount_org"] = amt
    for phase, x in zip(("aq", "org"), xs):
        for c in species:
            row[f"x_{phase}_{c}"] = x[ix[c]]
        row.update({f"w_{phase}_{n}": v for n, v in zip(case["names"], mass_fractions(case, x, species))})
        row[f"charge_relative_{phase}"] = math.fsum(z[c] * x[ix[c]] for c in species) / math.fsum(abs(z[c]) * x[ix[c]] for c in species)
    row["balance_sodium"] = total(lambda x: g(x, "sodium-cation")) - f["sodium-cation"]
    row["balance_acetate_moiety"] = total(lambda x: g(x, "acetate-anion") + g(x, "acetic-acid")) - f["acetate-anion"]
    row["balance_1-butanol"] = total(lambda x: g(x, "1-butanol")) - f["1-butanol"]
    row["balance_water_moiety"] = total(lambda x: g(x, "water") + g(x, "hydroxide-anion")) - f["water"]
    if "hmf" in ix:
        row["balance_hmf"] = total(lambda x: g(x, "hmf")) - f["hmf"]
    row["pressure_relative"] = max(abs(st[p].pressure / PRESSURE - 1) for p in range(2))
    lna = [[math.log(xs[p][i] * st[p].molar_density) + st[p].residual_chemical_potential_over_rt[i] for i in range(k)] for p in range(2)]
    d = [lna[1][i] - lna[0][i] for i in range(k)]
    row["transfer_neutral"] = max(abs(d[ix[c]]) for c in species if z[c] == 0)
    row["transfer_mean_ionic_NaAc"] = d[ix["sodium-cation"]] + d[ix["acetate-anion"]]
    row["offset_NaAc"] = (d[ix["sodium-cation"]] - d[ix["acetate-anion"]]) / 2
    if s["reacting"]:
        row["proton_condition"] = total(lambda x: g(x, "hydrogen-cation") + g(x, "acetic-acid") - g(x, "hydroxide-anion"))
        row["offset_H_minus_Na"] = d[ix["hydrogen-cation"]] - d[ix["sodium-cation"]]
        row["offset_OH_minus_Ac"] = d[ix["hydroxide-anion"]] - d[ix["acetate-anion"]]
        for j, (name, rx) in enumerate(REACTION_PK.items()):
            nu = rx["stoichiometry"]
            diag = r.neutral_reference_diagnostics[j]
            target = -rx["pK"] * LN10 + (sum(nu.values()) - nu.get("water", 0)) * math.log(M_W) + diag.terminal_contraction
            for p, phase in enumerate(("aq", "org")):
                row[f"affinity_{name}_{phase}"] = math.fsum(q * lna[p][ix[c]] for c, q in nu.items()) - target
    return row


def with_start(case, species, mix, reacting=True, **options):
    """All tie lines from a cold start; returns solved states in order (None where unsolved)."""
    return [solve(case, species, mix, i, reacting=reacting, **options) for i in range(len(case["ties"]))]


def coordinates(case):
    return [f"{p}_{n}" for p in ("aq", "org") for n in case["names"]]


def w_vector(case, row):
    return [row[f"w_{c}"] for c in coordinates(case)]


def write_csv(name, rows):
    fields = list(dict.fromkeys(k for r in rows for k in r))
    with (RESULTS / name).open("w", newline="") as stream:
        out = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        out.writeheader()
        out.writerows(rows)


def ln_ratio(a, b):
    return abs(math.log(a / b))


def primary(case):
    """Reacting and control states for every tie line, cold; plus each later tie line continued from the previous cold result."""
    species, control = P.SPECIES[case["system"]], case["control"]
    mix = mixture(case)
    out = dict(cold={}, continued={})
    for reacting, sp, m in ((True, species, mix), (False, control, mixture(case, control))):
        out["cold"][reacting] = [solve(case, sp, m, i, reacting=reacting) for i in range(len(case["ties"]))]
        assert all(s["ok"] for s in out["cold"][reacting]), "a primary cold solve failed"
        out["continued"][reacting] = [None] + [solve(case, sp, m, i, reacting=reacting, start=out["cold"][reacting][i - 1]["result"])
                                              for i in range(1, len(case["ties"]))]
    return out, mix


def state_rows(name, case, out):
    return [dict(case=name, row=int(case["rows"][i]["row"]), start=start, reacting=reacting, **read(s))
            for start in ("cold", "continued") for reacting in (True, False) for i, s in enumerate(out[start][reacting]) if s]


def stability(name, case, states):
    rows = []
    for i, s in enumerate(states):
        r = s["result"]
        compiled = eq.compile_problem(s["mix"], problem(s["species"], s["feed"], True))
        report = eq.stability_test(compiled, r, eq.StabilityOptions())
        rows.append(dict(case=name, row=int(case["rows"][i]["row"]), completed=report.completed, no_descent_in_trial_set=report.no_descent_in_trial_set,
                         trials=len(report.trials), resolved_descent=sum(t.resolved_descent for t in report.trials),
                         unresolved=sum(t.evidence != "no_resolved_descent" and not t.resolved_descent for t in report.trials),
                         message=report.message))
    return rows


def comparison(name, case, states, states_sigma):
    rows, coords = [], coordinates(case)
    for i, (s, s_sigma) in enumerate(zip(states, states_sigma)):
        a, a_sigma = w_vector(case, read(s)), w_vector(case, read(s_sigma))
        exp = case["ties"][i][0] + case["ties"][i][1]
        for j, c in enumerate(coords):
            e = case["u"][j] + E_B
            rows.append(dict(case=name, row=int(case["rows"][i]["row"]), coordinate=c, w_engine=a[j], w_experiment=exp[j],
                             deviation=abs(a[j] - exp[j]), e_A=case["u"][j], e_B=E_B, ratio=abs(a[j] - exp[j]) / e,
                             within=abs(a[j] - exp[j]) <= e, w_engine_sigma=a_sigma[j], ratio_sigma=abs(a_sigma[j] - exp[j]) / e,
                             within_sigma=abs(a_sigma[j] - exp[j]) <= e,
                             verdict_changes_with_diameter=(abs(a[j] - exp[j]) <= e) != (abs(a_sigma[j] - exp[j]) <= e)))
    return rows


def by_coordinate(rows):
    out = []
    for name in dict.fromkeys(r["case"] for r in rows):
        sub = [r for r in rows if r["case"] == name]
        stats = []
        for c in dict.fromkeys(r["coordinate"] for r in sub):
            cr = [r for r in sub if r["coordinate"] == c]
            stats.append(dict(case=name, coordinate=c, rows=len(cr), within=sum(r["within"] for r in cr),
                              mean_deviation=statistics.fmean(r["deviation"] for r in cr), e_A=cr[0]["e_A"],
                              mean_ratio=statistics.fmean(r["ratio"] for r in cr), max_ratio=max(r["ratio"] for r in cr),
                              verdict_changes=sum(r["verdict_changes_with_diameter"] for r in cr)))
        for rank, st in enumerate(sorted(stats, key=lambda s: s["mean_ratio"]), 1):
            out.append(dict(st, rank=rank))
    return out


def aad_rows(name, case, states):
    rows = [read(s) for s in states]
    dev = [[abs(w - e) for w, e in zip(w_vector(case, r), case["ties"][i][0] + case["ties"][i][1])] for i, r in enumerate(rows)]
    n = len(case["names"])
    reported = SOURCE["reported_absolute_average_deviation_w"]
    if name == "ternary":
        return [dict(case=name, group="rows 1-4", phase="both", engine=statistics.fmean(d for r in dev for d in r),
                     reported=reported["ternary"], reported_locator="Table 4")]
    out = []
    for group, ids, key in (("hmf-constant rows 1-3", (0, 1, 2), "hmf_constant_rows_1_3"), ("salt-constant rows 4-7", (3, 4, 5, 6), "salt_constant_rows_4_7")):
        for phase, sl in (("aq", slice(0, n)), ("org", slice(n, 2 * n))):
            out.append(dict(case=name, group=group, phase=phase, engine=statistics.fmean(d for i in ids for d in dev[i][sl]),
                            reported=reported[key][phase], reported_locator="Table 5"))
    return out


def aad_of(case, states, ids=None, phase=None):
    n = len(case["names"])
    sl = {None: slice(0, 2 * n), "aq": slice(0, n), "org": slice(n, 2 * n)}[phase]
    ids = range(len(states)) if ids is None else ids
    return statistics.fmean(abs(w - e) for i in ids if states[i] and states[i]["ok"]
                            for w, e in list(zip(w_vector(case, read(states[i])), case["ties"][i][0] + case["ties"][i][1]))[sl])


def retry(case, species, mix, base, **options):
    """Cold start for every tie line; one retry from the base solution where it fails. Returns states and a failure log."""
    out, log = [], []
    for i in range(len(case["ties"])):
        s = solve(case, species, mix, i, **options)
        if not s["ok"]:
            second = solve(case, species, mix, i, start=base[i]["result"], **options)
            log.append(dict(row=int(case["rows"][i]["row"]), cold_residual=s["result"].convergence.residual, cold_message=str(s["result"].message)[:200],
                            retry_from_base_ok=second["ok"]))
            s = second
        out.append(s if s["ok"] else None)
    return out, log


def shift(case, base, other):
    """Largest change in mass fractions and organic ln x of H+, OH-, HAc, salt against the base states."""
    pairs = [(read(b), read(o)) for b, o in zip(base, other) if o]
    dw = max(abs(x - y) for b, o in pairs for x, y in zip(w_vector(case, b), w_vector(case, o)))
    dln = {c: max(ln_ratio(o[f"x_org_{c}"], b[f"x_org_{c}"]) for b, o in pairs) for c in ("hydrogen-cation", "hydroxide-anion", "acetic-acid", "sodium-cation")}
    salt = max(abs(o["w_org_salt"] - b["w_org_salt"]) / b["w_org_salt"] for b, o in pairs)
    return dict(solved=len(pairs), rows=len(base), max_abs_dw=dw, org_dlnx_H=dln["hydrogen-cation"], org_dlnx_OH=dln["hydroxide-anion"],
                org_dlnx_HAc=dln["acetic-acid"], org_dlnx_Na=dln["sodium-cation"], org_salt_relative_shift=salt)


def ph_rows(name, case, states):
    species = P.SPECIES[case["system"]]
    mix = states[0]["mix"]
    ref = mix.reference(T, PRESSURE, solvent_fractions=[1.0] + [0.0] * (len(species) - 1), ions=("hydrogen-cation", "acetate-anion"), stoichiometry=(1, 1))
    out = []
    for i, s in enumerate(states):
        r, ix = s["result"], {c: j for j, c in enumerate(species)}
        xa = list(r.mole_fractions[:len(species)])
        st = mix.state(T, rho=r.molar_densities[0], x=xa)
        ln_a = math.log(xa[ix["hydrogen-cation"]]) + st.log_fugacity_coefficient[ix["hydrogen-cation"]] - ref.log_fugacity[0] - math.log(M_W)
        m_salt = xa[ix["sodium-cation"]] / (xa[ix["water"]] * M_W)
        out.append(dict(case=name, row=int(case["rows"][i]["row"]), molality_NaAc=m_salt, pH=-ln_a / LN10,
                        ideal_weak_base_estimate=7 + 0.5 * REACTION_PK["HAc"]["pK"] + 0.5 * math.log10(m_salt)))
    return out


def worst(rows, key, absolute=True):
    values = [abs(r[key]) if absolute else r[key] for r in rows if key in r and r.get("success", True)]
    return max(values), len(values)


def main():
    RESULTS.mkdir(exist_ok=True)
    states, comp, stab, aad, ph, coupling, sens, assoc, placement, sigma_shift = [], [], [], [], [], [], [], [], [], []
    sigma = {}
    summary = dict(bounds=BOUNDS, organic_ion_floor=ORGANIC_ION_FLOOR)
    for name, case in CASES.items():
        species = P.SPECIES[case["system"]]
        out, mix = primary(case)
        base_states, control_states = out["cold"][True], out["cold"][False]
        states += state_rows(name, case, out)
        sig_mix = mixture(case, dh_factor=1.0)
        sigma[name] = with_start(case, species, sig_mix)
        stab += stability(name, case, base_states)
        comp += comparison(name, case, base_states, sigma[name])
        aad += aad_rows(name, case, base_states)
        ph += ph_rows(name, case, base_states)
        for i, (a, b) in enumerate(zip(base_states, control_states)):
            ra, rb = read(a), read(b)
            coupling.append(dict(case=name, row=int(case["rows"][i]["row"]),
                                 max_abs_dw=max(abs(x - y) for x, y in zip(w_vector(case, ra), w_vector(case, rb))),
                                 max_abs_dx_substantial=max(abs(ra[f"x_{ph}_{c}"] - rb[f"x_{ph}_{c}"]) for ph in ("aq", "org") for c in case["control"])))
        sigma_shift.append(dict(case=name, max_abs_dw_sigma_vs_088=max(abs(x - y) for a, b in zip(base_states, sigma[name])
                                                                       for x, y in zip(w_vector(case, read(a)), w_vector(case, read(b))))))
        # unprinted-pair sensitivity: k_ij = +/-0.25 one pair at a time
        for pair in P.unprinted_pairs(case["system"]):
            for delta in (0.25, -0.25):
                pm = mixture(case, kij={pair: delta})
                solved, log = retry(case, species, pm, base_states)
                sens.append(dict(case=name, pair=f"{pair[0]}/{pair[1]}", k_ij=delta, cold_failures=len(log), **shift(case, base_states, solved)))
        # feed placement along each measured tie line
        for share in SHARES:
            solved, log = retry(case, species, mix, base_states, share=share) if share != 0.5 else (base_states, [])
            groups = {"all": None} if name == "ternary" else {"rows_4_7": (3, 4, 5, 6), "rows_1_3": (0, 1, 2), "all": None}
            for group, ids in groups.items():
                for phase in ("aq", "org"):
                    placement.append(dict(case=name, aqueous_mass_share=share, group=group, phase=phase, solved=sum(s is not None for s in solved),
                                          rows=len(solved), aad=aad_of(case, [s for s in solved], ids, phase)))
        # HMF cross-association variants
        if name == "hmf":
            for off in ((), ("1-butanol",), ("acetic-acid",), ("1-butanol", "acetic-acid")):
                am = mixture(case, hmf_cross_off=off)
                solved, log = retry(case, species, am, base_states)
                row = dict(variant="declared (Wolbach-Sandler)" if not off else "off: " + " + ".join(off), cold_failures=len(log), failures=json.dumps(log),
                           max_residual=max(s["result"].convergence.residual for s in solved if s), **shift(case, base_states, solved))
                for group, ids in (("rows_1_3", (0, 1, 2)), ("rows_4_7", (3, 4, 5, 6))):
                    for phase in ("aq", "org"):
                        row[f"aad_{group}_{phase}"] = aad_of(case, solved, ids, phase)
                assoc.append(row)
    write_csv("states.csv", states)
    write_csv("comparison.csv", comp)
    by = by_coordinate(comp)
    write_csv("comparison-by-coordinate.csv", by)
    write_csv("stability.csv", stab)
    write_csv("deviations.csv", aad)
    write_csv("ph.csv", ph)
    write_csv("coupling.csv", coupling)
    write_csv("unprinted-pair-sensitivity.csv", sens)
    write_csv("feed-placement.csv", placement)
    write_csv("hmf-association.csv", assoc)
    write_csv("diameter-convention.csv", sigma_shift)

    # falsifiers per solved state
    ok_rows = [r for r in states if r.get("success")]
    reacting = [r for r in ok_rows if r["reacting"]]
    f = {}

    def record(key, value, n, bound):
        f[key] = dict(worst=value, bound=bound, ok=value <= bound, n=n)
    record("row_residual", *worst(ok_rows, "residual"), BOUNDS["row_residual"])
    record("tolerance_met_all", sum(not r["tolerance_met"] for r in states), len(states), 0)
    for key in ("balance_sodium", "balance_acetate_moiety", "balance_1-butanol", "balance_water_moiety", "balance_hmf", "proton_condition"):
        subset = [r for r in ok_rows if key in r]
        if subset:
            record(key, *worst(subset, key), BOUNDS["conservation_mol"])
    record("charge_relative", max(abs(r[f"charge_relative_{p}"]) for r in ok_rows for p in ("aq", "org")), len(ok_rows) * 2, BOUNDS["charge_relative"])
    record("affinity", max(abs(r[k]) for r in reacting for k in r if k.startswith("affinity_")), len(reacting) * 4, BOUNDS["affinity"])
    record("transfer_neutral", *worst(ok_rows, "transfer_neutral"), BOUNDS["phase_transfer"])
    record("transfer_mean_ionic_NaAc", *worst(ok_rows, "transfer_mean_ionic_NaAc"), BOUNDS["phase_transfer"])
    record("offset_H_minus_Na", *worst(reacting, "offset_H_minus_Na"), BOUNDS["phase_transfer"])
    record("offset_OH_minus_Ac", *worst(reacting, "offset_OH_minus_Ac"), BOUNDS["phase_transfer"])
    record("pressure_relative", *worst(ok_rows, "pressure_relative"), BOUNDS["pressure_relative"])
    record("organic_ion_floor_violations", sum(min(r["x_org_sodium-cation"], r["x_org_acetate-anion"]) < ORGANIC_ION_FLOOR for r in reacting), len(reacting), 0)
    # cold against continued starts
    gaps = []
    for name in CASES:
        cold = {(r["row"], r["reacting"]): r for r in states if r["case"] == name and r["start"] == "cold"}
        for r in states:
            if r["case"] == name and r["start"] == "continued" and r.get("success"):
                c = cold[(r["row"], r["reacting"])]
                gaps.append(max(abs(r[k] - c[k]) for k in r if k.startswith("x_")))
    record("cold_vs_continued_x", max(gaps), len(gaps), BOUNDS["cold_vs_continued_x"])
    record("stability_failures", sum(not (s["completed"] and s["no_descent_in_trial_set"] and s["unresolved"] == 0 and s["resolved_descent"] == 0) for s in stab), len(stab), 0)
    record("coupling_w_ternary", max(c["max_abs_dw"] for c in coupling if c["case"] == "ternary"), 4, BOUNDS["coupling_w"])
    record("pair_sensitivity_w", max(s["max_abs_dw"] for s in sens), len(sens), BOUNDS["pair_sensitivity_w"])
    summary["falsifiers"] = f
    summary["certificate"] = dict(
        local_minimum="successful, non-singular return of the Engine's reduced-Hessian inertia certificate (ALG-EQUILIBRIUM)",
        states_success=sum(bool(r["success"]) for r in states), states=len(states),
        reduced_hessian_minimum_range=[min(r["reduced_hessian_minimum"] for r in reacting), max(r["reduced_hessian_minimum"] for r in reacting)])
    summary["coupling"] = {n: max(c["max_abs_dw"] for c in coupling if c["case"] == n) for n in CASES}
    summary["comparison"] = {n: dict(within=sum(r["within"] for r in comp if r["case"] == n), coordinates=sum(1 for r in comp if r["case"] == n),
                                     within_sigma=sum(r["within_sigma"] for r in comp if r["case"] == n),
                                     verdict_changes=[(r["row"], r["coordinate"]) for r in comp if r["case"] == n and r["verdict_changes_with_diameter"]])
                             for n in CASES}
    summary["all_falsifiers_ok"] = all(v["ok"] for v in f.values())
    (RESULTS / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    (RESULTS / "provenance.json").write_text(json.dumps(dict(
        engine=list(epcsaft.provenance()), T_K=T, P_Pa=PRESSURE,
        input_sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(INPUTS.glob("*.json")) + sorted((INPUTS / "source").iterdir())}), indent=1) + "\n")
    print(json.dumps({k: (v["worst"], v["bound"], v["ok"]) for k, v in f.items()}, indent=1))
    if not summary["all_falsifiers_ok"]:
        raise SystemExit("falsifier failed; see results/summary.json")


if __name__ == "__main__":
    main()
