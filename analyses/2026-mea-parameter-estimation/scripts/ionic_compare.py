"""Rank original Born, SSM and SSM+DS on observables that see the Born form.

The current MEA packet (revision 5) is evaluated with only (c_shell,
c_dielectric) changed, so every difference below is caused by the Born form:

  born-sensitivity.csv  which observables respond at fixed speciation:
    a speciated loaded liquid (density, fugacity coefficients, residual
    enthalpy) and aqueous MEAH+ MEACOO- mean ionic activity/osmotic
    coefficients.
  solvent-fugacity.csv  the comparison: Hilliard (2008) measured water and MEA
    partial pressures over 30 wt% MEA at 313.15/333.15 K against the model
    liquid fugacity x_i phi_i P at the Boettinger (2008) NMR speciation
    measured at the same temperature and MEA strength (interpolated in
    loading). No reaction equilibrium is solved; reaction constants and
    reactive CO2 pressure do not enter.

    .venv/bin/python analyses/2026-mea-parameter-estimation/results/analysis/ionic-electrostatics/compare.py
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path

from epcsaft import Mixture, Parameters

HERE = Path(__file__).resolve().parents[1] / "results/analysis/ionic-electrostatics"
ROOT = HERE.parents[4]
PACKET = ROOT / "data/packets/mea-co2-h2o-nine-species-estimation/5/parameters/parameters.json"
LIBRARY = ROOT / "analyses/2026-mea-parameter-estimation/inputs/source/library/observations"
VLE = LIBRARY / "vapor_liquid_equilibrium/Canonical_VLE_Observations.csv"
SPECIATION = LIBRARY / "liquid_speciation/Canonical_Combined_ChEq.csv"
VARIANTS = {"original-born": (0.0, 0.0), "ssm": (1.0, 0.0), "ssm-plus-ds": (1.0, 1.0)}
MAX_LOADING = 0.50  # above this the NMR bicarbonate/carbonate/CO2 split is unresolved
CARBAMATE_SCALES = (0.95, 1.0, 1.05)  # +-5 % NMR carbamate uncertainty; charge and amine balances re-close
TRACE = 1.0e-12  # mol per kg unloaded solution for species with no measured amount
MOLALITIES = (0.1, 0.5, 1.0, 2.0, 4.0)
IONS = ("protonated-monoethanolamine", "carbamate-anion")
# Assumed measurement scatter in ln(partial pressure) used by the #163 refit (solvent-fugacity-refit*.csv).
MEASUREMENT = {"water": 0.02, "monoethanolamine": 0.10}
# #161: Hilliard (2008) measured states and the stated uncertainties of its UT low-temperature apparatus
# (Zotero item DNQUPRMT, PDF attachment LB7LUETY; printed page = PDF page - 57).
D4 = LIBRARY / "vapor_liquid_equilibrium/Hilliard_2008_AppendixD4_7m_MEA.csv"
APPENDIX_F = LIBRARY / "liquid_speciation/Hilliard_2008_AppendixF_7m_MEA_NMR_speciation.csv"
U_VAPOR = 0.02  # p. 29: relative expanded uncertainty of the FT-IR vapor concentration (Goff 2005); no coverage factor stated, used as standard
U_PRESSURE = 0.0008 * 800.0 / 101.325  # p. 27: transducer +-0.08 % of its 0-8 bar full scale; reactor gas at atmospheric pressure
U_TEMPERATURE = 0.1  # K; p. 27: reactor temperature controlled within +-0.1 C
U_LOADING = 0.02  # relative; p. 29: standard uncertainty of the titrated loading
U_CARBAMATE = 0.05  # project assumption (#163), as CARBAMATE_SCALES; not a Hilliard statement (his pp. 62-63 give the amine-total error)
SLOPE_REFERENCE_T = 323.15  # K, middle of the 313-333 K data; k(T) = k + b (1/T - 1/T_ref), equations.md cross_dispersion_energy
DECISIVE_DELTA_CHI_SQUARE = 10.0  # equal parameter counts, so delta chi2 = delta AIC; >= 10 leaves the worse form essentially no support


def packet_document() -> dict:
    """Packet 5 with its historical CO2/MEA omission declared (#569): no cross edge, not a physical-absence claim."""
    document = json.loads(PACKET.read_text())
    document["topology"]["no_cross_association"] = [{
        "component_id_a": "carbon-dioxide", "component_id_b": "monoethanolamine",
        "provenance": {"source_id": document["document_id"], "locator": "packet 5: no CO2/MEA edge"}}]
    return document


def model(c_shell: float, c_dielectric: float, k_water: float = 0.0, k_mea: float = 0.0,
          families: tuple[str, ...] = (), factor: float = 1.0, water_slope: float = 0.0) -> Mixture:
    """Packet revision 5 with the Born form set; k_* replace the (zero) ion-water and ion-MEA k_ij of MEAH+ and
    MEACOO-, `water_slope` (K) makes the ion-water k_ij temperature dependent, and `factor` scales their listed
    coefficient families."""
    document = packet_document()
    for component in document["components"]:
        if component["component_id"] in IONS:
            for record in component["coefficients"]:
                if record["family"] in families:
                    record["value"]["magnitude"] *= factor
    family = next(f for f in document["model_families"] if f["family_id"] == "model/electrolyte")
    family["c_shell"], family["c_dielectric"] = c_shell, c_dielectric
    for pair in document["pairs"]:
        members = {pair["component_id_a"], pair["component_id_b"]}
        if members & set(IONS):
            for partner, k in (("water", k_water), ("monoethanolamine", k_mea)):
                if partner in members:
                    next(c for c in pair["coefficients"] if c["family"] == "k_ij")["value"]["magnitude"] = k
            if "water" in members and water_slope:
                pair["coefficients"].append(slope_record(pair, water_slope))
    return Mixture(Parameters.from_mapping(document))


def slope_record(pair: dict, slope: float) -> dict:
    return {"identity": f"pair/{pair['component_id_a']}/{pair['component_id_b']}/k_ij/reciprocal_temperature_slope",
            "family": "k_ij_reciprocal_temperature_slope", "value": {"magnitude": slope, "unit": "kelvin"},
            "reference_temperature": {"magnitude": SLOPE_REFERENCE_T, "unit": "kelvin"}}


def molar_masses() -> dict[str, float]:
    document = json.loads(PACKET.read_text())
    return {c["component_id"]: c["fixed"]["molar_mass"]["value"]["magnitude"] for c in document["components"]}


def speciation(temperature: float) -> list[tuple[float, float, float]]:
    """Boettinger 30 wt% rows: (loading, carbamate, bicarbonate aggregate) in mol/kg unloaded solution."""
    rows: dict[float, dict[str, float]] = {}
    with SPECIATION.open() as stream:
        for r in csv.DictReader(stream):
            if (r["source_key"], r["mea_mass_fraction"], float(r["temperature_K"])) == ("Bottinger2008", "0.3", temperature):
                rows.setdefault(float(r["co2_loading_mol_per_mol_mea"]), {})[r["species"]] = float(r["value_mol_per_kg_unloaded_solution"])
    return sorted((a, s["MEACOO-"], s.get("HCO3-", 0.0)) for a, s in rows.items())


def interpolate(table: list[tuple[float, float, float]], loading: float) -> tuple[float, float] | None:
    for (a0, c0, b0), (a1, c1, b1) in zip(table, table[1:]):
        if a0 <= loading <= a1:
            w = (loading - a0) / (a1 - a0)
            return c0 + w * (c1 - c0), b0 + w * (b1 - b0)
    return None


def amounts(loading: float, carbamate: float, bicarbonate: float, masses: dict[str, float]) -> dict[str, float]:
    """True-species amounts per kg unloaded 30 wt% solution; MEAH+ closes charge, MEA closes the amine balance."""
    total_mea = 0.3 / masses["monoethanolamine"]
    protonated = carbamate + bicarbonate + 2 * TRACE  # carbonate trace; hydroxide and hydronium traces cancel
    return {
        "carbon-dioxide": TRACE, "monoethanolamine": total_mea - carbamate - protonated,
        "water": 0.7 / masses["water"] - bicarbonate, "protonated-monoethanolamine": protonated,
        "carbamate-anion": carbamate, "bicarbonate-anion": bicarbonate, "carbonate-anion": TRACE,
        "hydronium-cation": TRACE, "hydroxide-anion": TRACE,
    }


def mole_fractions(ids: tuple[str, ...], n: dict[str, float]) -> list[float]:
    total = sum(n.values())
    return [n[i] / total for i in ids]


def sensitivity(mixtures: dict[str, Mixture], masses: dict[str, float]) -> list[dict]:
    rows = []
    loaded = amounts(0.40, 0.37 * 0.3 / masses["monoethanolamine"], 0.03 * 0.3 / masses["monoethanolamine"], masses)
    for name, mixture in mixtures.items():
        ids = mixture.component_ids
        state = mixture.state(313.15, P=1.0e5, x=mole_fractions(ids, loaded), phase="liquid")
        ln_phi = dict(zip(ids, state.log_fugacity_coefficient))
        for quantity, value in (
            ("molar_density_mol_m3", state.molar_density),
            ("bulk_relative_permittivity", state.bulk_relative_permittivity),
            ("born_residual_helmholtz", state.born),
            ("residual_enthalpy_J_mol", state.residual_enthalpy),
            ("ln_phi_water", ln_phi["water"]),
            ("ln_phi_monoethanolamine", ln_phi["monoethanolamine"]),
            ("ln_phi_carbon_dioxide", ln_phi["carbon-dioxide"]),
        ):
            rows.append({"model": name, "case": "30wt-loading-0.40-313.15K-1bar", "molality_mol_kg": "", "quantity": quantity, "value": value})
        water = ids.index("water")
        reference = mixture.reference(
            298.15, 1.0e5, solvent_fractions=[float(i == water) for i in range(len(ids))],
            ions=("protonated-monoethanolamine", "carbamate-anion"), stoichiometry=(1, 1),
        )
        for m in MOLALITIES:
            activity = reference.activity(m)
            rows.append({"model": name, "case": "aqueous-MEAH+MEACOO--298.15K", "molality_mol_kg": m,
                         "quantity": "mean_ionic_activity_coefficient_molality", "value": activity.mean_ionic_molality_scale})
            rows.append({"model": name, "case": "aqueous-MEAH+MEACOO--298.15K", "molality_mol_kg": m,
                         "quantity": "osmotic_coefficient", "value": activity.osmotic_coefficient})
    return rows


def comparison(mixtures: dict[str, Mixture], masses: dict[str, float], scales=CARBAMATE_SCALES) -> list[dict]:
    tables = {t: speciation(t) for t in (313.15, 333.15)}
    rows = []
    with VLE.open() as stream:
        observations = [r for r in csv.DictReader(stream) if r["source_key"] == "Hilliard2008" and r["MEA_weight_fraction"] == "0.3"]
    for r in observations:
        temperature, loading = float(r["temperature_reported_C"]) + 273.15, float(r["CO2_loading"])
        species = interpolate(tables[temperature], loading) if loading <= MAX_LOADING else None
        if species is None:
            continue  # outside the measured-speciation loading range at this temperature
        pressure = float(r["total_pressure"]) * 1.0e3
        closure = (species[0] + species[1]) / (loading * 0.3 / masses["monoethanolamine"])
        for (name, mixture), scale in ((m, s) for m in mixtures.items() for s in scales):
            n = amounts(loading, species[0] * scale, species[1], masses)
            ids = mixture.component_ids
            x = mole_fractions(ids, n)
            state = mixture.state(temperature, P=pressure, x=x, phase="liquid")
            ln_phi = dict(zip(ids, state.log_fugacity_coefficient))
            for component, column in (("water", "H2O_pressure"), ("monoethanolamine", "MEA_pressure")):
                i = ids.index(component)
                fugacity_kpa = x[i] * math.exp(ln_phi[component]) * pressure / 1.0e3
                observed = float(r[column])
                rows.append({
                    "observation_id": r["observation_id"], "temperature_K": temperature, "loading": loading,
                    "carbamate_scale": scale, "carbamate_mol_kg": species[0] * scale, "bicarbonate_mol_kg": species[1],
                    "nmr_co2_closure": closure,
                    "component": component, "observed_partial_pressure_kPa": observed, "model": name,
                    "model_fugacity_kPa": fugacity_kpa, "ln_model_over_observed": math.log(fugacity_kpa / observed),
                })
    return rows


def trends(rows: list[dict]) -> list[dict]:
    """Least-squares slope of ln(model/observed) against loading; zero means the loading response is reproduced.

    The Born term is inactive at zero loading, so a constant neutral-model offset sits in the
    intercept; the neutral model's loading dependence stays in the slope. The standard error is from row scatter.
    """
    groups: dict[tuple, list[tuple[float, float]]] = {}
    for r in rows:
        key = (r["component"], r["temperature_K"], r["model"], r["carbamate_scale"])
        groups.setdefault(key, []).append((r["loading"], r["ln_model_over_observed"]))
    summary = []
    for (component, temperature, name, scale), points in groups.items():
        n = len(points)
        mean_a = sum(a for a, _ in points) / n
        mean_e = sum(e for _, e in points) / n
        sxx = sum((a - mean_a) ** 2 for a, _ in points)
        slope = sum((a - mean_a) * (e - mean_e) for a, e in points) / sxx
        scatter = sum((e - mean_e - slope * (a - mean_a)) ** 2 for a, e in points) / (n - 2)
        summary.append({"component": component, "temperature_K": temperature, "model": name, "carbamate_scale": scale,
                        "rows": n, "mean_ln_model_over_observed": mean_e,
                        "rms_ln_model_over_observed": math.sqrt(sum(e * e for _, e in points) / n),
                        "slope_per_loading": slope, "slope_standard_error": math.sqrt(scatter / sxx)})
    return summary


PERTURBATIONS = (
    ("dispersion_energy_over_k", ("dispersion_energy_over_k",), (0.9, 1.1)),
    ("diameters", ("segment_diameter", "packing_diameter", "debye_huckel_diameter"), (0.95, 1.05)),
    ("born_diameter", ("born_diameter",), (0.9, 1.1)),
)


def coordinate_sensitivity(masses: dict[str, float]) -> list[dict]:
    """313.15 K loading slopes when one MEAH+/MEACOO- coordinate family or one ion-solvent k_ij is perturbed."""
    cases = [("none", "", {})]
    cases += [(label, f, {"families": families, "factor": f}) for label, families, factors in PERTURBATIONS for f in factors]
    cases += [(f"k_ion_{partner}", k, {f"k_{partner}": k}) for partner in ("water", "mea") for k in (-0.1, 0.1)]
    rows = []
    for name in ("original-born", "ssm-plus-ds"):
        for label, value, change in cases:
            summary = trends(comparison({name: model(*VARIANTS[name], **change)}, masses, scales=(1.0,)))
            slope = {r["component"]: r["slope_per_loading"] for r in summary if r["temperature_K"] == 313.15}
            rows.append({"model": name, "perturbation": label, "value": value,
                         "water_slope_per_loading": slope["water"], "mea_slope_per_loading": slope["monoethanolamine"]})
    return rows


def write(name: str, rows: list[dict]) -> None:
    with (HERE / name).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    masses = molar_masses()
    mixtures = {name: model(*pair) for name, pair in VARIANTS.items()}
    write("born-sensitivity.csv", sensitivity(mixtures, masses))
    rows = comparison(mixtures, masses)
    write("solvent-fugacity.csv", rows)
    summary = trends(rows)
    write("solvent-fugacity-trends.csv", summary)
    for r in summary:
        print("{component:17s} {temperature_K} {model:14s} x{carbamate_scale} n={rows:2d} mean={mean_ln_model_over_observed:+.3f} "
              "slope={slope_per_loading:+.3f}+-{slope_standard_error:.3f}".format(**r))
    write("coordinate-sensitivity.csv", coordinate_sensitivity(masses))
    print("packet sha256", hashlib.sha256(PACKET.read_bytes()).hexdigest())


def declarations():
    return [model(*v) for v in VARIANTS.values()] + [
        model(*VARIANTS[n], families=f, factor=x)
        for n in ('original-born', 'ssm-plus-ds')
        for _, f, xs in PERTURBATIONS
        for x in xs
    ]


if __name__ == "__main__":
    main()
