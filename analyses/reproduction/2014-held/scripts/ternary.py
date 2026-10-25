"""Held 2014 ternary systems from the installed Engine: Fig. 4a, Fig. 5 / Table 6 and Fig. 6a.

Every density, fugacity and liquid-liquid equilibrium calculation is native; Python declares feeds
and reads phase compositions. Fig. 4a Phi is the exception: a diagnostic Python contraction (see figure_4a).
"""
import csv
import json
import math
import time
from collections import defaultdict

from epcsaft import Mixture, Parameters
from epcsaft import equilibrium as eq

from model import HERE, IONS, NEUTRALS, document

T, P, MW = 298.15, 1.0e5, 0.01801528
ION_MASS = {"Na+": 22.98977, "NH4+": 18.0385, "Cl-": 35.453, "I-": 126.90447}  # g/mol


def write(name, rows):
    with (HERE / "results" / name).open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def mixture(cation, anion, neutral):
    doc, components = document(2, cation, anion, T, neutral=neutral)
    return Mixture(Parameters.from_mapping(doc, components=components)), components


def lle(mix, components, feed, aqueous, organic, organic_support, start):
    problem = eq.Problem(T=T, P=P, feed=eq.Amounts(dict(zip(components, feed))), row_tolerance=1e-10, phases=[
        eq.Phase("aqueous", kind="liquid", amount=eq.Free(), composition_guess=aqueous),
        eq.Phase("organic", kind="liquid", support=organic_support, amount=eq.Free(), composition_guess=organic)])
    result = eq.solve_equilibrium(mix, problem, start=start)
    assert result.success and result.requested_tolerance_met, feed
    return result


def figure_4a(lines, markers):
    """Osmotic coefficient of amino acid + 1 mol/kg KCl + water, Eq. (5) with sum of all solute molalities.

    Not a native result: the Engine's osmotic coefficient (electrolyte_activity) admits only one pure
    solvent and one salt, so a neutral solute plus a salt has no native Phi (capability gap, #387).
    The values below are a Python water-activity contraction of native states, kept as a diagnostic.
    """
    rows = []
    for neutral in ("glycine", "alanine"):
        mix, _ = mixture("K+", "Cl-", neutral)
        pure = mix.state(T, P=P, x=[1, 0, 0, 0], phase="liquid").log_fugacity_coefficient[0]
        paper = dict(lines.get(("figure-4a", f"{neutral} strategy 2"), []))
        # Filled black squares are glycine, light circles alanine (Fig. 4 caption).
        data = {float(r["x"]): float(r["y"]) for r in markers if r["figure"] == "figure-4a"
                and (float(r["core_mean_gray"]) < 100) == (neutral == "glycine")}
        grid = sorted(set(paper) | set(data) | {0.05 * i for i in range(0, 61) if neutral == "glycine" or i <= 40})
        for m in grid:
            n = [1 / MW, m, 1.0, 1.0]
            x = [v / sum(n) for v in n]
            state = mix.state(T, P=P, x=x, phase="liquid")
            assert abs(state.pressure - P) < 0.01 and state.pressure_density_derivative > 0
            phi = -(math.log(x[0]) + state.log_fugacity_coefficient[0] - pure) / (MW * (2.0 + m))
            rows.append(dict(solute=neutral, m_kcl_mol_kg=1.0, m_amino_acid_mol_kg=m, osmotic_python_contraction=phi,
                             osmotic_paper_line=paper.get(m, ""), osmotic_figure_data=data.get(m, "")))
    write("figure-4a.csv", rows)


def figure_5(markers):
    """Benzene solubility in the aqueous phase; ions are confined to the aqueous phase."""
    rows = []
    for anion in ("Cl-", "I-"):
        mix, comp = mixture("Na+", anion, "benzene")
        # Gray squares are NaI, black circles NaCl; the shared salt-free point has neither fill.
        data = {float(r["x"]): float(r["y"]) for r in markers if r["figure"] == "figure-5" and float(r["x"]) > 0.01
                and (float(r["core_mean_gray"]) < 50) == (anion == "Cl-")}
        for m in sorted({max(0.1 * i, 1e-6) for i in range(0, 37)} | set(data)):
            n_w = 1 / MW
            total = n_w + 2 * m
            r = lle(mix, comp, (n_w, 0.5, m, m), [n_w / total * 0.999, 4e-4, m / total, m / total],
                    [3e-3, 0.997], comp[:2], None)  # cold starts: a warm start from 0.1 mol/kg fails at 0.2
            x = r.mole_fractions[:4]
            rows.append(dict(salt="Na" + anion.rstrip("-"), m_salt_mol_kg=x[2] / (x[0] * MW),
                             m_benzene_mmol_kg=1000 * x[1] / (x[0] * MW), x_water_organic=r.mole_fractions[4],
                             m_salt_feed_mol_kg=m, m_benzene_figure_data_mmol_kg=data.get(m, ""),
                             iterations=r.iterations))
    write("figure-5.csv", rows)


def figure_6():
    """Water + 1-butanol + NH4Cl LLE; ions may enter both liquids; mass percent per phase."""
    mix, comp = mixture("NH4+", "Cl-", "1-butanol")
    mass = [MW * 1000, float(NEUTRALS["1-butanol"]["molar_mass_g_mol"]), ION_MASS["NH4+"], ION_MASS["Cl-"]]
    rows, start = [], None
    for i in range(0, 39):
        n = max(0.0025 * i, 1e-6)
        r = lle(mix, comp, (1.0, 0.25, n, n), [0.98 - 2 * n / 1.3, 0.015, n / 1.3, n / 1.3],
                [0.48, 0.52 - 2e-4, 1e-4, 1e-4], comp, start)
        start = r
        phases = []
        for p in range(2):
            x = r.mole_fractions[4 * p:4 * p + 4]
            w = [xi * mi for xi, mi in zip(x, mass)]
            phases.append([100 * v / sum(w) for v in w])
        aqueous, organic = sorted(phases, key=lambda w: -w[0])  # the solver's phase labels are not ordered
        rows.append(dict(w_nh4cl_aqueous=aqueous[2] + aqueous[3], w_butanol_aqueous=aqueous[1],
                         w_butanol_organic=organic[1], w_nh4cl_organic=organic[2] + organic[3],
                         w_water_organic=organic[0], iterations=r.iterations))
    write("figure-6a.csv", rows)


def main():
    start = time.perf_counter()
    lines = defaultdict(list)
    for r in csv.DictReader((HERE / "inputs/paper-lines.csv").open()):
        lines[(r["figure"], r["series"])].append((float(r["x"]), float(r["y"])))
    markers = list(csv.DictReader((HERE / "inputs/paper-markers.csv").open()))
    figure_4a(lines, markers)
    figure_5(markers)
    figure_6()
    identity = json.loads((HERE / "results/run-identity.json").read_text())
    identity["ternary_wall_clock_seconds"] = time.perf_counter() - start
    (HERE / "results/run-identity.json").write_text(json.dumps(identity, indent=2) + "\n")


if __name__ == "__main__":
    assert set(ION_MASS) <= set(IONS)
    main()
