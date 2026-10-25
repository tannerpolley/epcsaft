"""Ascani 2023 esterification: simultaneous chemical and liquid-liquid equilibrium, natively in the Engine.

K_a (Table 4) is the activity constant of Eq. 10 with pure-liquid references at the system T and p:
ln K_a = sum nu_i ln(x_i gamma_i), gamma_i = phi_i(x) / phi_i(pure liquid). It is declared as a raw mole-fraction
source correlation with every species on the Engine's PureLiquid reference. Each feed is one native solve with
phase discovery; every K_a value is the Engine's own reaction row. Python declares feeds and reads states.
"""
import csv
import itertools
import json
import math
import os
import sys
import time
from pathlib import Path

import epcsaft
from epcsaft import equilibrium as eq

sys.path.insert(0, str(Path(__file__).resolve().parent))
import parameters  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
P = 1.0e5  # Table 4, 1 bar
NU = (-1, -1, 1, 1)  # acid + alcohol -> ester + water, Eqs. 18-19
RULES = ("printed", "wolbach-sandler")
TOLERANCES = dict(row=1e-10, independent_check=1e-8,  # stated in README before the reacting runs
                  single_liquid_ln_ka=1e-9)  # revision: native reaction row recovers ln K_a at one-liquid states
X = ("x_acetic_acid", "x_1_pentanol", "x_pentyl_acetate", "x_water")
TIMING = []


def model(system, rule):
    path = ROOT / f"inputs/parameters/system-{system}-{rule}.json"
    return epcsaft.Mixture(epcsaft.Parameters.from_json(path, components=parameters.SYSTEMS[system]))


def problem(system, feed, phases=None):
    S, T, ka = parameters.SYSTEMS[system], parameters.temperatures()[system], parameters.ka()[system]
    return eq.Problem(
        T=T, P=P, feed=eq.Amounts(dict(zip(S, feed))), phases=phases or eq.Discover(), row_tolerance=TOLERANCES["row"],
        reactions=[eq.Reaction(dict(zip(S, NU)), name="esterification", correlation=eq.ReactionLogPolynomial(
            math.log(ka), 0, 0, 0, T, T, T, f"ascani-2023-table-4-system-{system}"),
            reference=eq.ReactionReference(eq.ReferenceSourceBasis.RawMoleFractionInfiniteDilution))],
        neutral_reference=eq.NeutralReference(
            solvent="water", standard_molality_mol_per_kg=1.0, solute_molality_ratios={s: 1.0 for s in S[:3]},
            species_reference={s: eq.SpeciesReference.PureLiquid for s in S}))


def native_ln_ka(system, mix):
    """ln K_a(x) from the compiled single-liquid declaration's chemical rows (ALG-EQUILIBRIUM, unscaled).

    Point = (ln x_i, N = 1, ln eta, balance multipliers 0). Row i is ln(x_i/v) + mu_i^r/RT + standard_i - (B^T lambda)_i;
    B nu = 0 removes the multipliers and sum nu = 0 removes ln v and ln Z, so sum nu_i row_i = ln K_a(x) - ln K.
    eta defaults to the liquid root at T and 1 bar.
    """
    T, ln_k = parameters.temperatures()[system], math.log(parameters.ka()[system])
    compiled = eq.compile_problem(mix, problem(system, [0.25] * 4, [eq.Phase("L", kind="liquid")]))

    def ln_ka(x, eta=None):
        eta = eta or mix._model.packing_fraction_coefficient(T, x) * mix.state(T, P=P, x=x, phase="liquid").molar_density
        rows = compiled.evaluate([*map(math.log, x), 1.0, math.log(eta), *[0.0] * compiled.rank]).values
        return math.fsum(n * v for n, v in zip(NU, rows[:4])) + ln_k
    return ln_ka


def solve(system, mix, feed, phases=None):
    T, ka = parameters.temperatures()[system], parameters.ka()[system]
    t0 = time.perf_counter()
    r = eq.solve_equilibrium(mix, problem(system, feed, phases))
    row = dict(success=r.success, tolerance_met=r.requested_tolerance_met, iterations=r.iterations,
               solve_ms=round(1000 * (time.perf_counter() - t0), 1), phases=len(r.amounts) if r.success else 0,
               objective=r.objective)
    if not r.success:
        row["message"] = r.message
        return row, []
    ln_ka = native_ln_ka(system, mix)
    phases = []
    for p, amount in enumerate(r.amounts):
        x = list(r.mole_fractions[4 * p:4 * p + 4])
        s = mix.state(T, rho=r.molar_densities[p], x=x)
        phases.append(dict(x=x, amount=amount, state=s,
                           reaction_residual=ln_ka(x, r.packing_fractions[p]) - math.log(ka), p_rel=s.pressure / P - 1))
    phases.sort(key=lambda ph: -ph["x"][3])  # water-rich first
    ln_f = [[math.log(xi) + lp for xi, lp in zip(ph["x"], ph["state"].log_fugacity_coefficient)] for ph in phases]
    row["reaction_residual_max"] = max(abs(ph["reaction_residual"]) for ph in phases)
    row["pressure_rel_max"] = max(abs(ph["p_rel"]) for ph in phases)
    row["ln_fugacity_difference_max"] = max((abs(a - b) for a, b in zip(ln_f[0], ln_f[-1])), default=0.0)
    row["mass_balance_rel_max"] = max(abs(math.fsum(ph["amount"] * ph["x"][i] for ph in phases)
                                          - (feed[i] + NU[i] * extent(feed, phases))) for i in range(4)) / sum(feed)
    return row, phases


def extent(feed, phases):
    return math.fsum(ph["amount"] * ph["x"][2] for ph in phases) - feed[2]


def read(name):
    with (ROOT / "inputs/data" / name).open() as stream:
        return list(csv.DictReader(stream))


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


def ka_points():
    """Native K_a at each Senina homogeneous point; the band moves every printed x_i by +-0.0005."""
    out = []
    for rule in RULES:
        ln_ka = native_ln_ka(1, model(1, rule))
        for r in read("senina-2022-table-2-homogeneous-ce.csv"):
            x = [float(r[k]) for k in X]
            band = []
            for signs in itertools.product((-1, 1), repeat=4):
                y = [max(xi + 0.0005 * s, 1e-6) for xi, s in zip(x, signs)]
                band.append(math.exp(ln_ka([v / sum(y) for v in y])))
            out.append(dict(rule=rule, point=int(r["point"]), **dict(zip(X, x)), K_a=math.exp(ln_ka(x)),
                            K_a_band_min=min(band), K_a_band_max=max(band)))
    write("ka-at-experimental-points.csv", out)


def homogeneous():
    """Feed = experimental point; the reaction moves it along its stoichiometric line to the model CE.

    A one-liquid answer is rechecked by declared two-liquid solves seeded from each experimental tie-line;
    the lowest converged two-liquid Gibbs objective minus the one-liquid objective is reported.
    """
    seeds = [([max(float(r["aq_" + k]), 1e-4) for k in X], [float(r["org_" + k]) for k in X])
             for r in read("senina-2022-table-3-reactive-tie-lines.csv")]
    out = []
    for rule in RULES:
        mix = model(1, rule)
        for r in read("senina-2022-table-2-homogeneous-ce.csv"):
            x = [float(r[k]) for k in X]
            row, phases = solve(1, mix, x)
            row = dict(rule=rule, point=int(r["point"]), **row)
            if len(phases) == 1:
                two = [solve(1, mix, x, [eq.Phase("aq", kind="liquid", composition_guess=a),
                                         eq.Phase("org", kind="liquid", composition_guess=o)])[0] for a, o in seeds]
                found = [t["objective"] for t in two if t["success"]]
                row["two_liquid_seeds_converged"] = len(found)
                if found:
                    row["two_liquid_minus_one_liquid_objective"] = min(found) - row["objective"]
            if phases:
                row["extent_minus_experiment"] = extent(x, phases)
                for name, ph in zip(("phase1", "phase2"), phases):
                    row |= {f"{name}_{k}": v for k, v in zip(X, ph["x"])} | {f"{name}_amount": ph["amount"]}
            out.append(row)
    write("homogeneous-ce.csv", out)


def tie_lines():
    """Feed = midpoint of each experimental reactive tie-line (Senina Table 3)."""
    out = []
    for rule in RULES:
        mix = model(1, rule)
        for r in read("senina-2022-table-3-reactive-tie-lines.csv"):
            aq = [float(r["aq_" + k]) for k in X]
            org = [float(r["org_" + k]) for k in X]
            feed = [(a + o) / 2 for a, o in zip(aq, org)]
            row, phases = solve(1, mix, feed)
            row = dict(rule=rule, tie_line=int(r["tie_line"]), **row)
            for name, ph in zip(("aq", "org"), phases):
                row |= {f"{name}_{k}": v for k, v in zip(X, ph["x"])}
            out.append(row)
    write("reactive-tie-lines.csv", out)


def binary_feeds():
    """Acid + alcohol feeds without products in both systems: where does the chemical equilibrium split?

    A failed phase discovery is kept in the row; the feed is then solved with declared one- and two-liquid
    topologies (the two-liquid seeded from the previous split) and the lower Gibbs objective is reported.
    """
    out = []
    for system in parameters.SYSTEMS:
        mix = model(system, "printed")
        split = None
        for k in range(1, 20):
            z = k / 20
            feed = [z, 1 - z, 0.0, 0.0]
            row, phases = solve(system, mix, feed)
            route = "discovery"
            if not phases and split:
                message = row["message"]
                one, _ = solve(system, mix, feed, [eq.Phase("L", kind="liquid", composition_guess=split[1]["x"])])
                row, phases = solve(system, mix, feed, [eq.Phase(n, kind="liquid", composition_guess=ph["x"])
                                                        for n, ph in zip(("aq", "org"), split)])
                route = "declared two liquids after discovery failure"
                row |= dict(discovery_message=message, one_liquid_objective=one["objective"])
            if len(phases) == 2:
                split = phases
            row = dict(system=system, z_acid=z, route=route, **row)
            for name, ph in zip(("phase1", "phase2"), phases):
                row |= {f"{name}_x_{s}": v for s, v in zip(("acid", "alcohol", "ester", "water"), ph["x"])}
                row[f"{name}_amount"] = ph["amount"]
            out.append(row)
    write("acid-alcohol-feeds.csv", out)


def discovery_cases():
    """Feeds near the gap edges where phase discovery and a declared two-liquid solve disagree."""
    seeds = {1: ([0.3, 0.003, 0.003, 0.69], [0.4, 0.05, 0.27, 0.28]),
             2: ([0.48, 0.0013, 0.0055, 0.51], [0.58, 0.012, 0.22, 0.19])}
    out = []
    for system, z in ((1, 0.66), (2, 0.75), (2, 0.76), (2, 0.78)):
        mix, feed = model(system, "printed"), [z, 1 - z, 0.0, 0.0]
        aq, org = seeds[system]
        found, _ = solve(system, mix, feed)
        one, _ = solve(system, mix, feed, [eq.Phase("L", kind="liquid", composition_guess=org)])
        two, phases = solve(system, mix, feed, [eq.Phase("aq", kind="liquid", composition_guess=aq),
                                                eq.Phase("org", kind="liquid", composition_guess=org)])
        out.append(dict(system=system, z_acid=z, discovery_success=found["success"], discovery_phases=found["phases"],
                        discovery_message=found.get("message", ""), discovery_objective=found["objective"] if found["success"] else "",
                        one_liquid_objective=one["objective"], two_liquid_success=two["success"],
                        two_liquid_objective=two["objective"], two_liquid_minus_one_liquid=two["objective"] - one["objective"],
                        minor_phase_amount=min(ph["amount"] for ph in phases),
                        x_water_phases=" ".join(f"{ph['x'][3]:.4f}" for ph in phases),
                        reaction_residual_max=two["reaction_residual_max"],
                        ln_fugacity_difference_max=two["ln_fugacity_difference_max"]))
    write("discovery-cases.csv", out)


def main():
    RESULTS.mkdir(exist_ok=True)
    started = time.strftime("%Y-%m-%dT%H:%M:%S%z")
    engine = epcsaft.provenance()  # before results are rewritten
    start = time.perf_counter()
    step("K_a at experimental points", ka_points)
    step("homogeneous chemical equilibrium", homogeneous)
    step("reactive tie-lines", tie_lines)
    step("acid + alcohol feeds", binary_feeds)
    step("discovery cases", discovery_cases)
    TIMING.append(dict(step="total", wall_s=round(time.perf_counter() - start, 3)))
    write("timing.csv", TIMING)
    (RESULTS / "provenance.json").write_text(json.dumps(dict(
        engine=engine, started=started, pressure_pa=P, temperatures_k=parameters.temperatures(),
        k_a=parameters.ka(), tolerances=TOLERANCES, cpu_count=os.cpu_count()), indent=1, default=str) + "\n")


def declarations():
    return [(model(s, r), problem(s, [0.25] * 4)) for s in (1, 2) for r in RULES]


if __name__ == "__main__":
    main()
