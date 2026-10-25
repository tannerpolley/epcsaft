"""Held 2014 aqueous osmotic coefficients from the installed Engine (Table 5, Figs. 1-3)."""
import csv
import json
import math
import time
from collections import defaultdict
from pathlib import Path

import epcsaft
from epcsaft import Mixture, Parameters

from model import HERE, REPO, document

P, MW = 1.0e5, 0.01801528
CHLORIDES = REPO / "analyses/2026-held-2014-chlorides"
CATIONS = ("H+", "Li+", "Na+", "K+", "Cs+", "NH4+")
ANIONS = ("F-", "Cl-", "Br-", "I-", "ClO4-", "NO3-", "H2PO4-", "OH-", "Ac-", "Fo-")
HW_SALT = {"Ac-": "Ac", "Fo-": "Fo"}  # Hamer-Wu label suffix; others drop the charge sign
FIGURE_SALTS = {"NaCl": ("Na+", "Cl-"), "KBr": ("K+", "Br-"), "KCl": ("K+", "Cl-"),
                "LiAc": ("Li+", "Ac-"), "NaAc": ("Na+", "Ac-"), "KAc": ("K+", "Ac-")}


def water_permittivity(T):
    """Ascani 2022 Table 5 water correlation (Held-group convention); 78.09 at 298.15 K."""
    return -105.2 * math.log(T) + 677.480


def salt_label(cation, anion):
    return cation.rstrip("+") + HW_SALT.get(anion, anion.rstrip("-"))


class Salt:
    """One fully dissociated 1:1 electrolyte; H+ is hydronium formed from one water per proton."""

    def __init__(self, strategy, cation, anion, T, dh_factor=0.88):
        ion = "H3O+" if cation == "H+" else cation
        doc, self.components = document(strategy, ion, anion, T, water_permittivity(T), dh_factor)
        self.model = Mixture(Parameters.from_mapping(doc, components=self.components))
        self.acid, self.T = cation == "H+", T
        self.reference = self.model.reference(
            T, P, solvent_fractions=[1, 0, 0], ions=tuple(self.components[1:]), stoichiometry=(1, 1))
        self.pure = self.model.state(T, P=P, x=[1, 0, 0], phase="liquid").log_fugacity_coefficient[0]

    def osmotic(self, m):
        """Held 2014 Eq. (5), native, with m per kg of initial water; returns (Phi, check diagnostics).

        For an acid, Eq. (10) consumes one water per proton, so the solution equals the native
        formula-unit solution at m' = m / (1 - m Mw) per kg of remaining water; Eq. (5) on the
        initial-water basis is then the exact basis change Phi = Phi' m'/m.
        """
        m_native = m / (1 - m * MW) if self.acid else m
        phi = self.reference.activity(m_native).osmotic_coefficient * m_native / m
        # Independent check only: water-activity contraction of the native state at the same composition.
        n_water = 1 / MW - (m if self.acid else 0.0)
        total = n_water + 2 * m
        x = [n_water / total, m / total, m / total]
        state = self.model.state(self.T, P=P, x=x, phase="liquid")
        phi_check = -(math.log(x[0]) + state.log_fugacity_coefficient[0] - self.pure) / (2 * m * MW)
        return phi, dict(pressure_residual_Pa=state.pressure - P, pressure_slope=state.pressure_density_derivative,
                         contraction_difference=phi - phi_check)


def declarations():
    """Construct, without solving, every Parameters document this analysis uses."""
    printed = {(int(r["strategy"]), r["cation"], r["anion"])
               for r in csv.DictReader((HERE / "inputs/held-2014-table-5-osmotic-ard.csv").open())}
    docs = [document(s, "H3O+" if c == "H+" else c, a, 298.15, water_permittivity(298.15), dh)
            for s, c, a in sorted(printed) for dh in (0.88, 1.0)]
    docs += [document(s, *FIGURE_SALTS[name], T, water_permittivity(T))
             for s in (1, 2) for name in FIGURE_SALTS for T in (273.15, 298.15)]
    docs += [document(2, c, a, 298.15, neutral=n) for c, a, n in (("K+", "Cl-", "glycine"), ("K+", "Cl-", "alanine"),
             ("Na+", "Cl-", "benzene"), ("Na+", "I-", "benzene"), ("NH4+", "Cl-", "1-butanol"))]
    return [Mixture(Parameters.from_mapping(doc, components=components)) for doc, components in docs]


def write(name, rows):
    with (HERE / "results" / name).open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


CHECK = {"acid": 0.0, "salt": 0.0}  # largest |native - Python contraction| in Phi, by electrolyte kind


def checked(salt, m):
    phi, diag = salt.osmotic(m)
    kind = "acid" if salt.acid else "salt"
    CHECK[kind] = max(CHECK[kind], abs(diag["contraction_difference"]))
    assert abs(diag["pressure_residual_Pa"]) < 0.01 and diag["pressure_slope"] > 0
    # Native Phi against the independent Python contraction, in the log domain because Phi divides by 2 m Mw.
    assert abs(diag["contraction_difference"]) * 2 * m * MW < 1e-11 and math.isfinite(phi), (m, diag)
    return phi, diag


def main():
    start = time.perf_counter()
    (HERE / "results").mkdir(exist_ok=True)
    observed = defaultdict(list)
    for r in csv.DictReader((HERE / "inputs/hamer-wu-1972-osmotic.csv").open()):
        if float(r["molality_mol_kg"]) <= 5.0:  # Table 5: concentrations up to 5 m
            observed[r["salt"]].append((float(r["molality_mol_kg"]), float(r["osmotic_coefficient"])))
    # LiCl and NaCl Hamer-Wu rows are reused from the reviewed chloride analysis transcription.
    for r in csv.DictReader((CHLORIDES / "inputs/evaluated.csv").open()):
        if r["salt"] in ("LiCl", "NaCl"):
            observed[r["salt"]].append((float(r["molality_mol_kg"]), float(r["osmotic_evaluated"])))
    observed = {salt: sorted(points) for salt, points in observed.items()}
    printed = {(int(r["strategy"]), r["cation"], r["anion"]): float(r["ard_percent"])
               for r in csv.DictReader((HERE / "inputs/held-2014-table-5-osmotic-ard.csv").open())}

    # Table 5 osmotic-coefficient replay against Hamer-Wu. The sigma closest-approach
    # sensitivity omits LiCl and NaCl: analyses/2026-held-2014-chlorides already holds them.
    rows, summary = [], []
    for dh_factor, strategy in ((0.88, 1), (0.88, 2), (1.0, 1), (1.0, 2)):
        for cation in CATIONS:
            for anion in ANIONS:
                label = salt_label(cation, anion)
                key = (strategy, cation, anion)
                if label not in observed or key not in printed or (dh_factor == 1.0 and label in ("LiCl", "NaCl")):
                    continue
                salt = Salt(strategy, cation, anion, 298.15, dh_factor)
                deviations = []
                for m, phi_hw in observed[label]:
                    phi, diag = checked(salt, m)
                    deviations.append(abs(1 - phi / phi_hw))
                    rows.append(dict(dh_factor=dh_factor, strategy=strategy, salt=label, molality_mol_kg=m, osmotic_engine=phi,
                                     osmotic_hamer_wu=phi_hw, relative_deviation=phi / phi_hw - 1))
                summary.append(dict(dh_factor=dh_factor, strategy=strategy, cation=cation, anion=anion, salt=label,
                                    n_rows=len(deviations), max_molality=observed[label][-1][0],
                                    ard_engine_percent=100 * sum(deviations) / len(deviations),
                                    ard_printed_percent=printed[key]))
    write("table-5-osmotic-rows.csv", rows)
    write("table-5-osmotic-ard.csv", summary)

    # Figures 1-3: Engine against the digitized model lines; figure 2 also at 273.15 K.
    lines = defaultdict(list)
    for r in csv.DictReader((HERE / "inputs/paper-lines.csv").open()):
        if r["figure"] in ("figure-1", "figure-2", "figure-3a", "figure-3b"):  # binaries; ternaries in ternary.py
            lines[(r["figure"], r["series"])].append((float(r["x"]), float(r["y"])))
    curve_rows, line_rows = [], []
    grid = [0.001 * 4000 ** (i / 120) for i in range(121)]
    for (figure, series), points in sorted(lines.items()):
        name = series.split()[0]
        cation, anion = FIGURE_SALTS[name]
        strategy = 1 if "strategy 1" in series else 2
        T = 273.15 if "273.15" in series else 298.15
        salt, text_salt = Salt(strategy, cation, anion, T), Salt(strategy, cation, anion, T, 1.0)
        for m in grid:
            curve_rows.append(dict(figure=figure, series=series, strategy=strategy, temperature_K=T,
                                   molality_mol_kg=m, osmotic_engine=checked(salt, m)[0],
                                   osmotic_engine_sigma_closest_approach=checked(text_salt, m)[0]))
        for m, phi_paper in sorted(points):
            line_rows.append(dict(figure=figure, series=series, strategy=strategy, temperature_K=T,
                                  molality_mol_kg=m, osmotic_paper_line=phi_paper,
                                  osmotic_engine=checked(salt, m)[0],
                                  osmotic_engine_sigma_closest_approach=checked(text_salt, m)[0]))
    write("figure-curves.csv", curve_rows)
    write("figure-line-comparison.csv", line_rows)

    # Figure 2 filled markers are the 273.15 K Lobo-Quaresma points; assign by curve proximity.
    marker_rows = []
    cold = {name: Salt(2, *FIGURE_SALTS[name], 273.15) for name in ("NaCl", "KCl")}
    markers = [(float(r["x"]), float(r["y"])) for r in csv.DictReader((HERE / "inputs/paper-markers.csv").open())
               if r["figure"] == "figure-2"]
    for m, y in markers:
        pair = [q for q in markers if abs(q[0] - m) < 0.02]
        if len(pair) == 2:  # NaCl lies above KCl wherever both are plotted at one molality
            name = "NaCl" if y == max(q[1] for q in pair) else "KCl"
        else:  # otherwise the nearer digitized 273.15 K line
            nearest = {n: min(lines[("figure-2", f"{n} 273.15 K")], key=lambda p: abs(p[0] - m))[1] for n in cold}
            name = min(cold, key=lambda n: abs(nearest[n] - y))
        marker_rows.append(dict(salt=name, temperature_K=273.15, molality_mol_kg=m, osmotic_digitized=y,
                                osmotic_engine=checked(cold[name], m)[0]))
    write("figure-2-cold-markers.csv", marker_rows)

    (HERE / "results/run-identity.json").write_text(json.dumps({
        "engine_provenance": list(epcsaft.provenance()),
        "conditions": "1 bar; 298.15 K except Fig. 2 cold series at 273.15 K; molality per kg (initial) water",
        "water_permittivity": {"298.15": water_permittivity(298.15), "273.15": water_permittivity(273.15),
                               "source": "Ascani 2022 Table 5: -105.2 ln(T/K) + 677.480 (supplement; Held 2014 prints none)"},
        "max_native_minus_python_contraction_osmotic": CHECK,
        "wall_clock_seconds": time.perf_counter() - start,
    }, indent=2) + "\n")
    for row in summary:
        print(row["dh_factor"], row["strategy"], row["salt"], round(row["ard_engine_percent"], 2), row["ard_printed_percent"])


if __name__ == "__main__":
    main()
