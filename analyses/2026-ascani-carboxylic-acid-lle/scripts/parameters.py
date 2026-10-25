"""Write the Ascani 2022 ESI acetic acid + water + toluene parameter documents (Tables S1-S5).

Variants differ only in the declared convention under test:
  literal       water is the only dielectric solvent (ESI Table S5 lists water and salts only);
                water sigma = 2.7927 + 10.11 exp(-0.01775 T) - 1.417 exp(-0.01146 T) (Fuchs 2006, the
                ESI's ref. 18; the group's Ascani 2023 Table 1 prints this footnote with the same 2.7927).
  organic-eps   literal + CRC relative permittivities of acetic acid and toluene, needed once ions
                may enter the toluene-rich liquid.
  constant-sigma  literal with the ESI's printed water sigma 2.7927 taken as a constant.
  cr1-cross     literal with the arithmetic-energy, geometric-volume water-acetic acid cross rule.
"""
import itertools
import json
import math
from pathlib import Path

T = 298.15
OUT = Path(__file__).resolve().parents[1] / "inputs"
ESI = "ascani-et-al-2022-esi"
CRC = "crc-handbook-permittivity-liquids"
DOMAIN = "ascani-2022-298-15-k-1-bar"

# component: (charge, molar mass g/mol or None, m, sigma, u/k, ESI locator)
COMPONENTS = {
    "water": (0, 18.01528, 1.2047, None, 353.95, "Table S1 water row"),
    "toluene": (0, 92.13842, 2.8149, 3.7169, 285.69, "Table S1 toluene row"),
    "acetic-acid": (0, 60.05196, 1.3403, 3.8582, 211.59, "Table S1 acetic acid row"),
    "hydrogen-cation": (1, None, 1.0, 3.4656, 500.0, "Table S2 H+ row"),
    "acetate-anion": (-1, None, 1.0, 3.9328, 150.0, "Table S2 Ac- row"),
    "hydroxide-anion": (-1, None, 1.0, 2.0177, 650.0, "Table S2 OH- row"),
}
SELF_ASSOCIATION = {"water": (2425.7, 0.04510), "acetic-acid": (3044.4, 0.07555)}
K_IJ = {("water", "toluene"): 0.02100, ("water", "acetic-acid"): -0.12470, ("toluene", "acetic-acid"): -0.08317}
WATER_EPS = -105.2 * math.log(T) + 677.480  # Table S5
SALT_EPS = 8.0  # Table S5, every ion
# CRC Handbook "Permittivity (Dielectric Constant) of Liquids", eps = A + B T + C T^2, read from the
# table copy shipped in chemicals (Electrolytes/Permittivity (Dielectric Constant) of Liquids.tsv).
ORGANIC_EPS = {"acetic-acid": (-15.731, 0.12662, -0.00017738), "toluene": (3.2584, -0.003441, 0.0000015937)}


def value(magnitude, unit):
    return {"magnitude": magnitude, "unit": unit}


def provenance(locator, source=ESI):
    return {"source_id": source, "locator": locator, "domain_id": DOMAIN}


def coefficient(identity, family, magnitude, unit, locator, source=ESI):
    return {"identity": identity, "family": family, "value": value(magnitude, unit),
            "provenance": provenance(locator, source)}


def document(variant):
    components, correlations, sites, edges, pairs = [], [], [], [], []
    for cid, (charge, mass, m, sigma, u, locator) in COMPONENTS.items():
        fixed = {"charge_number": {"value": value(charge, "elementary-charge"), "provenance": provenance(locator)}}
        if mass is not None:
            fixed["molar_mass"] = {"value": value(mass, "gram / mole"),
                                   "provenance": provenance("IUPAC atomic weights", "iupac-atomic-weights")}
        coefficients = [
            coefficient(f"component/{cid}/segment_count", "segment_count", m, "dimensionless", locator),
            coefficient(f"component/{cid}/dispersion_energy_over_k", "dispersion_energy_over_k", u, "kelvin", locator),
            coefficient(f"component/{cid}/solvation_factor", "solvation_factor", 1.0, "dimensionless",
                        "ESI Eq. 11 uncorrected Born term"),
        ]
        if cid == "water" and variant != "constant-sigma":
            correlations.append({
                "correlation_id": "component/water/segment_diameter/constant-plus-sum-of-exponentials",
                "component_id": "water", "family": "segment_diameter", "form": "constant-plus-sum-of-exponentials",
                "independent_variables": ["temperature"],
                "constant": {"identity": "component/water/segment_diameter/constant-plus-sum-of-exponentials/constant",
                             "value": value(2.7927, "angstrom")},
                "provenance": provenance("Table S1 water sigma 2.7927 with the Fuchs 2006 (ESI ref. 18) temperature terms"),
                "terms": [{"amplitude": {"identity": f"component/water/segment_diameter/constant-plus-sum-of-exponentials/term-{k}/amplitude",
                                         "value": value(a, "angstrom")},
                           "exponent_coefficient": {"identity": f"component/water/segment_diameter/constant-plus-sum-of-exponentials/term-{k}/exponent_coefficient",
                                                    "value": value(b, "1 / kelvin")}}
                          for k, (a, b) in enumerate(((10.11, -0.01775), (-1.417, -0.01146)))],
            })
        else:
            coefficients.append(coefficient(f"component/{cid}/segment_diameter", "segment_diameter",
                                            2.7927 if cid == "water" else sigma, "angstrom", locator))
        if charge:
            for family in ("packing_diameter", "debye_huckel_diameter", "born_diameter"):
                coefficients.append(coefficient(f"component/{cid}/{family}", family, sigma, "angstrom",
                                                f"{locator}; ESI Eq. 11 Born a_i = ion diameter sigma_i"))
            coefficients.append(coefficient(f"component/{cid}/relative_permittivity", "relative_permittivity",
                                            SALT_EPS, "dimensionless", "Table S5 salts"))
        elif cid == "water":
            coefficients.append(coefficient("component/water/relative_permittivity", "relative_permittivity",
                                            WATER_EPS, "dimensionless", "Table S5 water -105.2 ln(T/K) + 677.480 at 298.15 K"))
        elif variant == "organic-eps":
            a, b, c = ORGANIC_EPS[cid]
            coefficients.append(coefficient(f"component/{cid}/relative_permittivity", "relative_permittivity",
                                            a + b * T + c * T * T, "dimensionless",
                                            "A + B T + C T^2 at 298.15 K", CRC))
        components.append({"component_id": cid, "aliases": [], "fixed": fixed, "coefficients": coefficients})
    for cid, (energy, volume) in SELF_ASSOCIATION.items():
        sites += [{"component_id": cid, "site_id": s, "site_role": role, "multiplicity": 1,
                   "provenance": provenance(f"Table S1 {cid} N:N = 1:1")} for s, role in (("a", "donor"), ("b", "acceptor"))]
        edges.append({"endpoint_a": {"component_id": cid, "site_id": "a"}, "endpoint_b": {"component_id": cid, "site_id": "b"},
                      "energy_over_k": {"identity": f"association/{cid}/a/{cid}/b/energy_over_k", "value": value(energy, "kelvin")},
                      "volume": {"identity": f"association/{cid}/a/{cid}/b/volume", "value": value(volume, "dimensionless")},
                      "source": {"kind": "explicit", "provenance": [provenance(f"Table S1 {cid} association")]}})
    rule = "arithmetic-energy-geometric-volume" if variant == "cr1-cross" else "wolbach-sandler"
    for a_site, b_site in (("a", "b"), ("b", "a")):
        edges.append({"endpoint_a": {"component_id": "water", "site_id": a_site},
                      "endpoint_b": {"component_id": "acetic-acid", "site_id": b_site},
                      "source": {"kind": "combining-rule", "rule_id": rule,
                                 "inputs": {"energy_over_k": ["association/water/a/water/b/energy_over_k",
                                                              "association/acetic-acid/a/acetic-acid/b/energy_over_k"],
                                            "volume": ["association/water/a/water/b/volume",
                                                       "association/acetic-acid/a/acetic-acid/b/volume"]},
                                 "provenance": [provenance(f"ESI states no cross-association rule; {rule} under test")]}})
    for a, b in itertools.combinations(COMPONENTS, 2):
        if COMPONENTS[a][0] * COMPONENTS[b][0] > 0:
            continue  # like-charged ions take no dispersion pair
        k =K_IJ.get((a, b), K_IJ.get((b, a), 0.0))
        locator = "Table S3" if k else "Table S2 caption and main text: k_ij = 0 for every pair with an ion"
        pairs.append({"component_id_a": a, "component_id_b": b, "coefficients": [
            coefficient(f"pair/{a}/{b}/k_ij", "k_ij", k, "dimensionless", locator + " (k_ij,T term zero at 298.15 K)"),
            coefficient(f"pair/{a}/{b}/l_ij", "l_ij", 0.0, "dimensionless", "ESI uses no l_ij")]})
    return {
        "schema": "epcsaft.parameters", "schema_version": 1,
        "document_id": f"ascani-2022-acetic-acid-toluene-{variant}", "document_version": 1, "purpose": "external-document",
        "sources": [
            {"source_id": ESI, "citation": "Ascani, M.; Pabsch, D.; Klinksiek, M.; Gajardo-Parra, N.; Sadowski, G.; Held, C. "
             "Prediction of pH in multiphase multicomponent systems with ePC-SAFT advanced. Chem. Commun. 2022, 58, 8436; ESI.",
             "doi": "10.1039/D2CC02943J", "use_basis": "ESI Tables S1-S5 at 298.15 K and 1 bar"},
            {"source_id": CRC, "citation": "CRC Handbook of Chemistry and Physics, Permittivity (Dielectric Constant) of Liquids",
             "use_basis": "acetic acid and toluene relative permittivity, organic-eps variant only"},
            {"source_id": "iupac-atomic-weights", "citation": "IUPAC standard atomic weights", "use_basis": "molar masses"}],
        "domains": [{"domain_id": DOMAIN, "kind": "reported-conditions",
                     "temperature_min": value(T, "kelvin"), "temperature_max": value(T, "kelvin"),
                     "pressure_min": value(1, "bar"), "pressure_max": value(1, "bar")}],
        "components": components, "pairs": pairs,
        "model_families": [
            {"family_id": "model/association", "kind": "association", "choice": "general-site",
             "provenance": provenance("Table S1 association schemes")},
            {"family_id": "model/base", "kind": "base", "choice": "pc-saft", "provenance": provenance("ESI Eq. 9")},
            {"family_id": "model/electrolyte", "kind": "electrolyte", "choice": "born",
             "ionic_dispersion": "exclude-same-sign-ion-pairs", "c_shell": 0.0, "c_dielectric": 0.0,
             "provenance": provenance("ESI Eqs. 10-11; Held 2014 like-ion rule retained by ePC-SAFT advanced")},
            {"family_id": "model/relative_permittivity", "kind": "permittivity", "choice": "component-permittivity-mixing",
             "provenance": provenance("ESI ref. 7 concentration-dependent dielectric constant; Ascani-Held 2021 Eq. 13")}],
        "model_coefficients": [{"identity": "model/relative_permittivity/ionic_region_relative_permittivity",
                                "family": "ionic_region_relative_permittivity", "value": value(1.0, "dimensionless"),
                                "provenance": provenance("uncorrected Born term, ESI Eq. 11")}],
        "correlations": correlations,
        "topology": {"presets": [], "sites": sites, "edges": edges},
    }


if __name__ == "__main__":
    for variant in ("literal", "organic-eps", "constant-sigma", "cr1-cross"):
        (OUT / f"parameters-{variant}.json").write_text(json.dumps(document(variant), indent=1) + "\n")
