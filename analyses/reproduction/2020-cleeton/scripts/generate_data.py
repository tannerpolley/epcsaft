"""Cleeton 2020 native replay: Figure 16 (Dicko 2010 Table 1, 323.15 K, 50 wt% MDEA) and its check ladder.

Printed model only: Table 2 parameters (4C water and MDEA, one-segment non-associating ions, k_ij
with water only), Austgen 1991 Table V constants unconverted on the mole-fraction basis with H3O+
(Cleeton R.1-R.6), water and MDEA at the pure liquid and every other species at infinite dilution
in water (Eqs. 8-9), at the system pressure. Permittivity is water's alone (Eq. 2 without its
pressure term). No parameter is fitted. The main lane uses ion diameters 0.88 sigma (the Held-group
convention); the diagnostic lane uses sigma, and a second diagnostic puts MDEA at infinite dilution.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import runpy
from pathlib import Path

import numpy as np

import epcsaft
from epcsaft import GAS_CONSTANT_J_PER_MOL_K, equilibrium

ROOT = Path(__file__).resolve().parents[4]
CAMPAIGN = Path(__file__).resolve().parents[1]
SOURCE = CAMPAIGN / "inputs" / "source"
DIGITIZED = CAMPAIGN / "inputs" / "digitized"
DATA = CAMPAIGN / "results" / "data"
COMMON = runpy.run_path(
    str(ROOT / "analyses/reproduction/2018-wangler/figures/figure-10/scripts/generate_data.py")
)
PROVENANCE = {"source_id": "cleeton-2020", "locator": "Tables 1-2, Eqs. 2-9", "domain_id": "cleeton-2020"}
for helper in ("_component", "_pair", "_edge"):
    COMMON[helper].__globals__["PROVENANCE"] = PROVENANCE
component, pair, edge = COMMON["_component"], COMMON["_pair"], COMMON["_edge"]
water_permittivity = COMMON["_water_permittivity"]  # Uyan/Wangler transfer of eps_w(T); 69.86 at 323.15 K

# Table 2 (4C water and MDEA). id: (molar mass kg/mol, charge, m, sigma A, u/k K, k_ij with water)
SPECIES = {
    "water": (0.01801528, 0, 2.1945, 2.2290, 141.66, 0.0),
    "mdea": (0.119163, 0, 4.3589, 3.3474, 237.19, 0.0),
    "co2": (0.0440095, 0, 2.6037, 2.5550, 151.04, -0.009),
    "h2s": (0.034081, 0, 1.6686, 3.0349, 229.00, 0.028),
    "h3o-cation": (0.019023, 1, 1.0, 3.4654, 500.00, 0.250),
    "oh-anion": (0.01700734, -1, 1.0, 2.0177, 650.00, -0.250),
    "hco3-anion": (0.06101684, -1, 1.0, 2.9296, 70.00, 0.0),
    "carbonate-anion": (0.0600089, -2, 1.0, 2.4422, 249.26, -0.250),
    "mdeah-cation": (0.12017094, 1, 1.0, 3.3474, 237.19, 0.0),
    "hs-anion": (0.033073, -1, 1.0, 3.0349, 229.00, 0.0),
    "sulfide-anion": (0.032065, -2, 1.0, 3.0349, 229.00, 0.0),
}
# 4C association (2 donor + 2 acceptor sites): eps^AB/k K, kappa^AB.
ASSOCIATION = {"water": (1804.17, 0.20390), "mdea": (1491.60, 0.29836)}
# Cleeton R.1-R.6 with Austgen 1991 Table V (pp. 546-547): ln K_x = A + B/T + C ln T + D T, mole
# fraction, H3O+. K2 A is Austgen's 231.465 (Cleeton Table 1 prints 231.456). Validity in kelvin.
REACTIONS = {
    "K1": ({"water": -2, "oh-anion": 1, "h3o-cation": 1}, (132.899, -13445.9, -22.4773, 0.0), (273.15, 498.15)),
    "K2": ({"co2": -1, "water": -2, "hco3-anion": 1, "h3o-cation": 1}, (231.465, -12092.10, -36.7816, 0.0), (273.15, 498.15)),
    "K3": ({"hco3-anion": -1, "water": -1, "carbonate-anion": 1, "h3o-cation": 1}, (216.049, -12431.70, -35.4819, 0.0), (273.15, 498.15)),
    "K4": ({"mdeah-cation": -1, "water": -1, "mdea": 1, "h3o-cation": 1}, (-9.4165, -4234.98, 0.0, 0.0), (298.15, 333.15)),
    "K5": ({"h2s": -1, "water": -1, "hs-anion": 1, "h3o-cation": 1}, (214.582, -12995.4, -33.5471, 0.0), (273.15, 423.15)),
    "K6": ({"hs-anion": -1, "water": -1, "sulfide-anion": 1, "h3o-cation": 1}, (-32.0, -3338.0, 0.0, 0.0), (287.15, 343.15)),
}
TEMPERATURE = 323.15
REFERENCE_TEMPERATURE = 298.15
MDEA_MASS_FRACTION = 0.5
# Charge-neutral dyadic dilution paths (exact zero charge sum) per species set.
DILUTION_RATIOS = {
    ("co2", "h2s"): {"oh-anion": 0.125, "hco3-anion": 0.125, "carbonate-anion": 0.0625,
                     "hs-anion": 0.0625, "sulfide-anion": 0.03125},
    ("co2",): {"oh-anion": 0.25, "hco3-anion": 0.125, "carbonate-anion": 0.0625},
    ("h2s",): {"oh-anion": 0.25, "hs-anion": 0.125, "sulfide-anion": 0.0625},
}
LANES = {  # lane: (ion diameter / sigma, MDEA reference)
    "main-0.88sigma": (0.88, "pure-liquid"),
    "diagnostic-sigma": (1.0, "pure-liquid"),
    "diagnostic-mdea-infinite-dilution": (0.88, "infinite-dilution"),
}
# Amendment 1 gate, fixed before any native run: per point |ln(p_native/p_figure)| <= 0.05 +
# 0.003 |d ln p_figure / d loading|. 0.034 is 2.5 px of the log axis (170 px per decade), 0.016
# covers the K2 constant (231.465 vs 231.456) and the omitted permittivity pressure term, and
# 0.003 loading is 2.5 px of the loading axis carried through the curve slope.
GATE_BASE, GATE_LOADING = 0.05, 0.003


def species_set(gases: tuple[str, ...]) -> tuple[str, ...]:
    ions = {"co2": ("hco3-anion", "carbonate-anion"), "h2s": ("hs-anion", "sulfide-anion")}
    ids = ["water", "mdea", *gases, "h3o-cation", "oh-anion", "mdeah-cation"]
    for gas in gases:
        ids += ions[gas]
    return tuple(ids)


def build_parameters(ids: tuple[str, ...], ion_factor: float = 0.88) -> epcsaft.Parameters:
    components = []
    for cid in ids:
        mass, charge, m, sigma, u, _ = SPECIES[cid]
        extra = {"packing_diameter": ion_factor * sigma, "debye_huckel_diameter": ion_factor * sigma} if charge else {}
        components.append(component(cid, molar_mass=mass, charge_number=charge, segment_count=m,
                                    segment_diameter=sigma, dispersion_energy_over_k=u, **extra))
    pairs = [pair(a, b, SPECIES[b][5] if a == "water" else SPECIES[a][5] if b == "water" else 0.0)
             for i, a in enumerate(ids) for b in ids[i + 1:]
             if SPECIES[a][1] * SPECIES[b][1] <= 0]  # no same-sign ion pairs (revised ePC-SAFT, l.49)
    associating = [cid for cid in ("water", "mdea") if cid in ids]
    sites = [{"component_id": cid, "site_id": site, "site_role": role, "multiplicity": 2, "provenance": PROVENANCE}
             for cid in associating for site, role in (("a", "donor"), ("b", "acceptor"))]
    edges = [edge((cid, "a"), (cid, "b"), *ASSOCIATION[cid]) for cid in associating]
    if len(associating) == 2:  # Wolbach-Sandler cross association
        (ew, kw), (em, km) = ASSOCIATION["water"], ASSOCIATION["mdea"]
        sw, sm = SPECIES["water"][3], SPECIES["mdea"][3]
        volume = math.sqrt(kw * km) * (math.sqrt(sw * sm) / (0.5 * (sw + sm))) ** 3
        edges += [edge(("water", "a"), ("mdea", "b"), 0.5 * (ew + em), volume),
                  edge(("water", "b"), ("mdea", "a"), 0.5 * (ew + em), volume)]
    families = [
        {"family_id": "model/association", "kind": "association", "choice": "general-site", "provenance": PROVENANCE},
        {"family_id": "model/base", "kind": "base", "choice": "pc-saft", "provenance": PROVENANCE},
    ]
    correlations = []
    ionic = any(SPECIES[cid][1] for cid in ids)
    families += [
        {"family_id": "model/electrolyte", "kind": "electrolyte",
         **({"choice": "fully-dissociated-debye-huckel", "ionic_dispersion": "exclude-same-sign-ion-pairs"}
            if ionic else {"choice": "none"}), "provenance": PROVENANCE},
        {"family_id": "model/relative_permittivity", "kind": "permittivity",
         "choice": "solvent-only" if ionic else "none", "provenance": PROVENANCE},
    ]
    if ionic:
        # Water alone carries a permittivity, so the salt-free solvent rule returns eps_w at every composition.
        correlations.append({
            "correlation_id": "component/water/relative_permittivity/constant", "component_id": "water",
            "family": "relative_permittivity", "form": "constant", "independent_variables": [],
            "constant": {"identity": "component/water/relative_permittivity/constant/constant",
                         "value": {"magnitude": water_permittivity(TEMPERATURE), "unit": "dimensionless"}},
            "provenance": PROVENANCE,
        })
    return epcsaft.Parameters.from_mapping({
        "schema": "epcsaft.parameters", "schema_version": 1, "document_id": "cleeton-2020-table-2",
        "document_version": 1, "purpose": "external-document",
        "sources": [{"source_id": "cleeton-2020",
                     "citation": "Cleeton et al., Fluid Phase Equilibria 511 (2020) 112453.",
                     "use_basis": "Table 2 and Eqs. 2-9; Austgen et al. 1991 Table V constants."}],
        "domains": [{"domain_id": "cleeton-2020", "kind": "reported-conditions",
                     "temperature_min": {"magnitude": 293.0, "unit": "kelvin"},
                     "temperature_max": {"magnitude": 423.0, "unit": "kelvin"},
                     "pressure_min": {"magnitude": 1.0e-3, "unit": "pascal"},
                     "pressure_max": {"magnitude": 2.0e7, "unit": "pascal"},
                     "ion_mole_fraction_max": {"magnitude": 0.60, "unit": "dimensionless"}}],
        "components": components, "pairs": pairs, "model_families": families,
        "model_coefficients": [], "correlations": correlations,
        "topology": {"presets": [], "sites": sites, "edges": edges},
    })


def problem(ids: tuple[str, ...], loadings: dict[str, float], mdea_reference: str,
            pressure_seed_pa: float) -> equilibrium.Problem:
    """CH4-free bubble point of 1 mol MDEA at 50 wt% with the given acid-gas loadings (Cleeton's quaternary)."""
    gases = tuple(g for g in ("co2", "h2s") if g in ids)
    amounts = dict.fromkeys(ids, 0.0)
    amounts["mdea"] = 1.0
    amounts["water"] = (1.0 - MDEA_MASS_FRACTION) / MDEA_MASS_FRACTION * SPECIES["mdea"][0] / SPECIES["water"][0]
    for gas in gases:
        amounts[gas] = loadings[gas]
    # Declared start: half of each gas bound as MDEAH+ with its anion, traces of the other ions.
    extents = {"K1": 1.0e-6, "K2": 0.5 * loadings.get("co2", 0.0), "K3": 1.0e-4 * loadings.get("co2", 0.0),
               "K5": 0.5 * loadings.get("h2s", 0.0), "K6": 1.0e-4 * loadings.get("h2s", 0.0)}
    extents["K4"] = -(extents["K2"] + extents["K3"] + extents["K5"] + extents["K6"])
    reactions = [name for name, (nu, _, _) in REACTIONS.items() if set(nu) <= set(ids)]
    for name in reactions:
        for cid, nu in REACTIONS[name][0].items():
            amounts[cid] += nu * extents[name]
    total = math.fsum(amounts.values())
    neutral = [cid for cid in ("water", "mdea", "co2", "h2s") if cid in ids]
    vapor_guess = {"water": 0.02, "mdea": 1.0e-5}
    for gas in gases:
        vapor_guess[gas] = (1.0 - 0.02 - 1.0e-5) / len(gases)
    ratios = {"mdea": 0.125, "h3o-cation": 0.25, "mdeah-cation": 0.25, **{g: 0.125 for g in gases},
              **DILUTION_RATIOS[gases]}
    return equilibrium.Problem(
        T=TEMPERATURE,
        P=equilibrium.Free(pressure_seed_pa),
        phases=[
            equilibrium.Phase("liquid", kind="liquid", support=ids,
                              composition_guess=tuple(amounts[c] / total for c in ids)),
            equilibrium.Phase("vapor", kind="vapor", support=neutral, amount=equilibrium.Pinned(0.0),
                              composition_guess=tuple(vapor_guess[c] for c in neutral)),
        ],
        feed=equilibrium.Amounts(amounts),
        reactions=[
            equilibrium.Reaction(
                REACTIONS[name][0], name=f"cleeton-2020-{name}",
                correlation=equilibrium.ReactionLogPolynomial(
                    # The native form uses ln(T/Tref); A absorbs C ln(Tref).
                    a + c * math.log(REFERENCE_TEMPERATURE), b, c, d, REFERENCE_TEMPERATURE,
                    *REACTIONS[name][2], "austgen-1991-table-v-mole-fraction"),
                # Eqs. 8-9 evaluate both references at the system pressure.
                reference=equilibrium.ReactionReference(
                    equilibrium.ReferenceSourceBasis.RawMoleFractionInfiniteDilution, reference_pressure_pa=None),
            )
            for name in reactions for (a, b, c, d) in (REACTIONS[name][1],)
        ],
        row_tolerance=1e-10,
        neutral_reference=equilibrium.NeutralReference(
            solvent="water", solute_molality_ratios=ratios, standard_molality_mol_per_kg=1.0,
            species_reference={"mdea": equilibrium.SpeciesReference.PureLiquid}
            if mdea_reference == "pure-liquid" else None),
    )


def stoichiometry(ids):
    return {name: [REACTIONS[name][0].get(c, 0.0) for c in ids]
            for name, (nu, _, _) in REACTIONS.items() if set(nu) <= set(ids)}


def affinity(result, ids) -> float:
    """Largest |sum_i nu_ri r_i| / sum|nu_ri| over the liquid chemical rows (moiety multipliers cancel)."""
    chemical = {row.index: value for row, value in zip(result.rows, result.residuals, strict=True)
                if row.kind.name == "Chemical" and row.phase == 0}
    return max(abs(math.fsum(v * chemical[i] for i, v in enumerate(nu) if v)) / math.fsum(map(abs, nu))
               for nu in stoichiometry(ids).values())


def ln_phi(model, x, pressure=None, density=None):
    state = model.state(TEMPERATURE, P=pressure, x=x, phase="liquid") if density is None else \
        model.state(TEMPERATURE, rho=density, x=x)
    z = state.pressure / (state.molar_density * GAS_CONSTANT_J_PER_MOL_K * TEMPERATURE)
    return state, [mu - math.log(z) for mu in state.residual_chemical_potential_over_rt]


def solve(model, ids, loadings, mdea_reference, seed_pa) -> dict[str, object]:
    """One CH4-free bubble point; returns partial pressures and the per-row diagnostics (amendments 5, 7)."""
    prob = problem(ids, loadings, mdea_reference, seed_pa)
    try:
        result = equilibrium.solve_equilibrium(model, prob)
    except epcsaft.ParameterError as exc:
        return {"success": False, "failure": str(exc)}
    out = {"success": result.success, "solver_status": result.solver_status,
           "failure": "" if result.success else result.message}
    if not result.success:
        return out
    n = len(ids)
    x, y = list(result.mole_fractions[:n]), list(result.mole_fractions[n:])
    neutral = [c for c in ("water", "mdea", "co2", "h2s") if c in ids]
    pressure, density = result.pressure, result.molar_densities[0]
    state, liquid = ln_phi(model, x, density=density)
    ratios = prob.neutral_reference.solute_molality_ratios
    trace = [(1.0 / SPECIES["water"][0] if c == "water" else 0.0) + 1.0e-12 * ratios.get(c, 0.0) for c in ids]
    _, dilute = ln_phi(model, [t / math.fsum(trace) for t in trace], pressure=pressure)
    k = ids.index("mdeah-cation")
    out.update({
        "pressure_pa": pressure,
        **{f"p_{c}_pa": y[neutral.index(c)] * pressure for c in neutral},
        "liquid_density_kg_m3": density * math.fsum(xi * SPECIES[c][0] for xi, c in zip(x, ids)),
        "mdeah_ln_gamma_star": liquid[k] - dilute[k],  # vs infinite dilution in water, same T and P
        "mdeah_fraction_of_mdea": x[k] / (x[k] + x[ids.index("mdea")]),
        "liquid_permittivity": state.bulk_relative_permittivity,
        "reaction_affinity_max": affinity(result, ids),
        "maximum_original_row": max(map(abs, result.residuals)),
        "reference_status": ";".join(sorted({d.status.name for d in result.neutral_reference_diagnostics})),
        "reference_admitted_change_max": max((abs(d.observed_terminal_change)
                                             for d in result.neutral_reference_diagnostics if d.observed_terminal_change is not None), default=None),
        "iterations": result.iterations,
    })
    return out


def ladder() -> list[dict[str, object]]:
    """Amendment 2: pure MDEA (Fig. 2) and the water binaries (Fig. 4 P-x) against the digitized symbols."""
    with (DIGITIZED / "cleeton-2020-figures-2-4-symbols.csv").open(newline="", encoding="utf-8") as stream:
        symbols = list(csv.DictReader(stream))
    rows = []
    mdea = epcsaft.Pure(build_parameters(("mdea",)))
    for s in symbols:
        if not s["panel"].startswith("fig2"):
            continue
        t = float(s["x"])
        r = equilibrium.solve_equilibrium(mdea, mdea.saturation(t))
        observed = float(s["y"])
        native = r.pressure if s["panel"] == "fig2-vapor-pressure" else r.molar_densities[0] * SPECIES["mdea"][0]
        rows.append({"check": s["panel"], "temperature_k": t, "composition": 1.0, "digitized": observed,
                     "native": native if r.success else None, "success": r.success,
                     "relative_deviation": native / observed - 1.0 if r.success else None})
    for gas in ("h2s", "co2"):
        model = epcsaft.Mixture(build_parameters(("water", gas)))
        for s in symbols:
            if s["panel"] != f"fig4-{gas}-water":
                continue
            t, xg, observed = float(s["series"]), float(s["x"]), float(s["y"])
            prob = equilibrium.Problem(
                T=t, P=equilibrium.Free(observed), feed=equilibrium.MoleFractions({"water": 1.0 - xg, gas: xg}),
                phases=[equilibrium.Phase("L", kind="liquid", composition_guess=(1.0 - xg, xg)),
                        equilibrium.Phase("V", kind="vapor", amount=equilibrium.Pinned(0.0),
                                          composition_guess=(0.02, 0.98))])
            r = equilibrium.solve_equilibrium(model, prob)
            rows.append({"check": f"{gas}-water bubble pressure", "temperature_k": t, "composition": xg,
                         "digitized": observed, "native": r.pressure if r.success else None,
                         "success": r.success,
                         "relative_deviation": r.pressure / observed - 1.0 if r.success else None})
    return rows


def mdea_infinite_dilution() -> dict[str, float]:
    """Amendment 3: ln gamma_inf of MDEA in water = ln phi(trace in water) - ln phi(pure MDEA), same T, P."""
    model = epcsaft.Mixture(build_parameters(("water", "mdea")))
    out = {}
    for pressure in (1.0e5, 2.0e6):
        _, dilute = ln_phi(model, [1.0 - 1.0e-12, 1.0e-12], pressure=pressure)
        _, pure = ln_phi(model, [1.0e-12, 1.0 - 1.0e-12], pressure=pressure)
        out[f"ln_gamma_inf_mdea_in_water_{pressure:.0e}_pa"] = dilute[1] - pure[1]
    return out


def load_figure_curves() -> dict[tuple[str, float], tuple[np.ndarray, np.ndarray]]:
    curves: dict[tuple[str, float], list[tuple[float, float]]] = {}
    with (DIGITIZED / "cleeton-2020-figure-16-curves.csv").open(newline="", encoding="utf-8") as stream:
        for r in csv.DictReader(stream):
            curves.setdefault((r["panel"], float(r["series_loading"])), []).append(
                (float(r["loading"]), math.log(float(r["pressure_pa"]))))
    return {k: (np.array([p[0] for p in v]), np.array([p[1] for p in v])) for k, v in curves.items()}


def main() -> None:
    with (SOURCE / "dicko-2010-table-1.csv").open(newline="", encoding="utf-8") as stream:
        dicko = list(csv.DictReader(stream))
    curves = load_figure_curves()
    ladder_rows = ladder()
    write_csv(DATA / "cleeton-2020-ladder.csv", ladder_rows)
    rows, gate = [], []
    for lane, (factor, mdea_reference) in LANES.items():
        models = {}

        def model_for(gases):
            if gases not in models:
                models[gases] = (epcsaft.Mixture(build_parameters(species_set(gases), factor)), species_set(gases))
            return models[gases]

        for d in dicko:
            loadings = {"co2": float(d["co2_loading"]), "h2s": float(d["h2s_loading"])}
            gases = tuple(g for g in ("co2", "h2s") if loadings[g] > 0)
            model, ids = model_for(gases)
            seed = 1.0e6 * (float(d["co2_pressure_mpa"]) + float(d["h2s_pressure_mpa"])) + 1.2e4
            row = {"lane": lane, **d, **solve(model, ids, loadings, mdea_reference, seed)}
            for gas in gases:
                if row["success"]:
                    row[f"{gas}_relative_deviation"] = row[f"p_{gas}_pa"] / (1.0e6 * float(d[f"{gas}_pressure_mpa"])) - 1.0
            rows.append(row)
            print(json.dumps({k: row.get(k) for k in ("lane", "row", "success", "p_co2_pa", "p_h2s_pa")}), flush=True)
        # Amendment 1: native at Figure 16's own states at each symbol. The gate reads the second-gas
        # loading from the caption (0.1/0.3/0.5/0.8); the diagnostic reading uses the Dicko group's own
        # value (CO2: the group's loading; H2S: the mean over the panel's series).
        model, ids = model_for(("co2", "h2s"))
        for d, reading, panel, gas, other, series_column in (
                (d, reading, *spec) for reading in ("caption", "dicko-group") for d in dicko
                for spec in (("left", "co2", "h2s", "figure_16_left_series_h2s_loading"),
                             ("right", "h2s", "co2", "figure_16_right_series_co2_loading"))):
            if not d[series_column]:
                continue
            series = float(d[series_column])
            members = [float(r[f"{other}_loading"]) for r in dicko if r[series_column] == d[series_column]]
            state = series if reading == "caption" else math.fsum(members) / len(members)
            xs, lnp = curves[(panel, series)]
            loading = min(max(float(d[f"{gas}_loading"]), xs.min()), xs.max())  # within the traced span
            figure = float(np.interp(loading, xs, lnp))
            slope = float(np.polyfit(xs[np.abs(xs - loading) <= 0.02], lnp[np.abs(xs - loading) <= 0.02], 1)[0])
            result = solve(model, ids, {gas: loading, other: state}, mdea_reference, 2.0 * math.exp(figure) + 1.2e4)
            entry = {"lane": lane, "reading": reading, "panel": panel, "row": d["row"], "gas": gas,
                     "series_loading": series, "second_gas_loading": state, "loading": loading,
                     "figure_pressure_pa": math.exp(figure), "figure_slope": slope,
                     "tolerance": GATE_BASE + GATE_LOADING * abs(slope), **result}
            if result["success"]:
                entry["ln_ratio"] = math.log(result[f"p_{gas}_pa"]) - figure
                entry["within_tolerance"] = abs(entry["ln_ratio"]) <= entry["tolerance"]
            gate.append(entry)
    write_csv(DATA / "cleeton-2020-dicko-rows.csv", rows)
    write_csv(DATA / "cleeton-2020-figure-16-gate.csv", gate)
    summary = summarize(rows, gate, ladder_rows)
    summary["mdea_infinite_dilution"] = mdea_infinite_dilution()
    (DATA / "cleeton-2020-summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


def aad(values) -> float | None:
    values = [abs(v) for v in values if v is not None]
    return 100.0 * math.fsum(values) / len(values) if values else None


def symmetric(relative) -> float | None:
    """Mean 2|c-e|/(c+e) in percent from relative deviations r = c/e - 1."""
    values = [2.0 * abs(r) / (2.0 + r) for r in relative if r is not None]
    return 100.0 * math.fsum(values) / len(values) if values else None


def curve_deviation(lane_gate, panel, column=None):
    """c/e - 1 of a Figure 16 curve value (native or digitized) against the Dicko point it is drawn beside."""
    with (SOURCE / "dicko-2010-table-1.csv").open(newline="", encoding="utf-8") as stream:
        dicko = {r["row"]: r for r in csv.DictReader(stream)}
    return [(g[column or f"p_{g['gas']}_pa"] / (1.0e6 * float(dicko[g["row"]][f"{g['gas']}_pressure_mpa"])) - 1.0)
            if g["success"] else None for g in lane_gate if g["panel"] == panel]


def summarize(rows, gate, ladder_rows) -> dict[str, object]:
    wheel = next((ROOT / "build/environment-wheel").glob("epcsaft-*.whl"))
    published = {"fig2-vapor-pressure": 6.30, "fig2-density": 0.40, "h2s-water bubble pressure": 7.26,
                 "co2-water bubble pressure": 13.72}
    out = {
        "wheel_sha256": hashlib.sha256(wheel.read_bytes()).hexdigest(),
        "permittivity_expected": water_permittivity(TEMPERATURE),
        "gate_tolerance": f"|ln(p_native/p_figure)| <= {GATE_BASE} + {GATE_LOADING} |d ln p_figure/d loading|",
        "ladder": {check: {"points": sum(r["check"] == check for r in ladder_rows),
                           "solved": sum(r["check"] == check and r["success"] for r in ladder_rows),
                           "native_aad_percent": aad(r["relative_deviation"] for r in ladder_rows if r["check"] == check),
                           "native_symmetric_percent": symmetric(r["relative_deviation"] for r in ladder_rows
                                                                 if r["check"] == check),
                           "published_aad_percent": value}
                   for check, value in published.items()},
        "lanes": {},
    }
    for lane in LANES:
        lane_rows = [r for r in rows if r["lane"] == lane]
        solved = [r for r in lane_rows if r["success"]]
        mixed = [r for r in lane_rows if r["figure_16_right_series_co2_loading"]]
        left = [r for r in lane_rows if r["figure_16_left_series_h2s_loading"]]
        lane_gate = [g for g in gate if g["lane"] == lane and g["reading"] == "caption"]
        group_gate = [g for g in gate if g["lane"] == lane and g["reading"] == "dicko-group"]
        groups = {}
        for r in mixed:
            groups.setdefault(r["co2_loading"], []).append(r)
        out["lanes"][lane] = {
            "rows": len(lane_rows), "solved": len(solved),
            "reaction_affinity_max": max((r["reaction_affinity_max"] for r in solved), default=None),
            "reference_status": sorted({r["reference_status"] for r in solved}),
            "reference_admitted_change_max": max((r["reference_admitted_change_max"] for r in solved), default=None),
            "liquid_density_kg_m3_range": [min(r["liquid_density_kg_m3"] for r in solved), max(r["liquid_density_kg_m3"] for r in solved)],
            "mdeah_ln_gamma_star_range": [min(r["mdeah_ln_gamma_star"] for r in solved), max(r["mdeah_ln_gamma_star"] for r in solved)],
            "permittivity_range": [min(r["liquid_permittivity"] for r in solved), max(r["liquid_permittivity"] for r in solved)] if solved else None,
            "figure_16_gate": {
                "points": len(lane_gate), "solved": sum(g["success"] for g in lane_gate),
                "within_tolerance": sum(bool(g.get("within_tolerance")) for g in lane_gate),
                "max_abs_ln_ratio": max((abs(g["ln_ratio"]) for g in lane_gate if g["success"]), default=None),
                "median_abs_ln_ratio": float(np.median([abs(g["ln_ratio"]) for g in lane_gate if g["success"]])) if lane_gate else None,
            },
            "figure_16_diagnostic_dicko_group_loading": {
                "within_tolerance": sum(bool(g.get("within_tolerance")) for g in group_gate),
                "max_abs_ln_ratio": max((abs(g["ln_ratio"]) for g in group_gate if g["success"]), default=None),
                "median_abs_ln_ratio": float(np.median([abs(g["ln_ratio"]) for g in group_gate if g["success"]])) if group_gate else None,
            },
            # Cleeton's AAD definition is not printed; 2|calc-exp|/(calc+exp) on its own digitized curves gives
            # 55.9% and 79.1% against the published 57.18% and 79.32% (inferred reading).
            "dicko_symmetric_deviation_percent": {
                "co2_figure_16_left_12": symmetric(r.get("co2_relative_deviation") for r in left),
                "h2s_figure_16_right_18": symmetric(r.get("h2s_relative_deviation") for r in mixed),
                "co2_curve_states_left_12": symmetric(curve_deviation(lane_gate, "left")),
                "h2s_curve_states_right_18": symmetric(curve_deviation(lane_gate, "right")),
                "digitized_curve_left_12": symmetric(curve_deviation(lane_gate, "left", "figure_pressure_pa")),
                "digitized_curve_right_18": symmetric(curve_deviation(lane_gate, "right", "figure_pressure_pa")),
            },
            "dicko_aad_percent": {
                "co2_figure_16_left_12": aad(r.get("co2_relative_deviation") for r in left),
                "h2s_figure_16_right_18": aad(r.get("h2s_relative_deviation") for r in mixed),
                "co2_mixed_18": aad(r.get("co2_relative_deviation") for r in mixed),
                "co2_all_23": aad(r.get("co2_relative_deviation") for r in lane_rows if float(r["co2_loading"]) > 0),
                "h2s_all_23": aad(r.get("h2s_relative_deviation") for r in lane_rows if float(r["h2s_loading"]) > 0),
                "published_co2": 57.18, "published_h2s": 79.32,
            },
            "overpredicted_mixed_rows": {gas: sum((r.get(f"{gas}_relative_deviation") or 0.0) > 0 for r in mixed)
                                         for gas in ("co2", "h2s")},
            "co2_pressure_rises_with_h2s_loading": {
                group: all(a["p_co2_pa"] < b["p_co2_pa"] for a, b in zip(g, g[1:]))
                for group, g in groups.items() if all(r["success"] for r in g)},
        }
    return out


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(dict.fromkeys(k for r in rows for k in r)), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def declarations():
    from epcsaft import Mixture, Pure

    return (
        Pure(build_parameters(('mdea',))),
        Mixture(build_parameters(('water', 'co2'))),
        [
            (
                Mixture(build_parameters(species_set(g), f)),
                problem(species_set(g), {'co2': 0.3, 'h2s': 0.3}, ref, 1e5),
            )
            for f, ref in LANES.values()
            for g in (('co2', 'h2s'), ('co2',), ('h2s',))
        ],
    )


if __name__ == "__main__":
    main()
