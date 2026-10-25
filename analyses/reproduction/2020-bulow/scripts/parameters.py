"""Bülow, Ascani & Held 2020 (ePC-SAFT advanced I) parameter documents, Tables 1-4, verbatim.

One base document holds water, methanol, ethanol and the seven ions; `document(variant)` changes only the
convention under test. Every ion diameter is the unified d_i = 0.88 sigma_i: packing (Held 2014 Eq. 21),
Born (reproduces the Fig. 3 Born bars) and Debye-Hückel closest approach (reproduces the corrigendum's
revised curves and Fig. 6 right-panel contributions; Held 2014 Sec. 2.1 prints sigma, kept as a sensitivity).
"""
import itertools
import json
from pathlib import Path

PAPER = "bulow-2020"
HELD = "held-2014"
DOMAIN = "bulow-2020-298-15-k-1-bar"
# component: (charge, molar mass g/mol or None, m, sigma/A or None (water sigma(T)), u/k, relative permittivity, locator)
COMPONENTS = {
    "water": (0, 18.01528, 1.2047, None, 353.95, 78.09, "Table 1 water row and footnote; Table 4"),
    "methanol": (0, 32.04186, 1.5255, 3.2300, 188.90, 33.05, "Table 1 methanol row; Table 4"),
    "ethanol": (0, 46.06844, 2.3827, 3.1771, 198.24, 24.88, "Table 1 ethanol row; Table 4"),
    "lithium-cation": (1, None, 1.0, 2.8449, 360.00, 8.0, "Table 2 Li+ row; Table 4 salts"),
    "sodium-cation": (1, None, 1.0, 2.8232, 230.00, 8.0, "Table 2 Na+ row; Table 4 salts"),
    "potassium-cation": (1, None, 1.0, 3.3417, 200.00, 8.0, "Table 2 K+ row; Table 4 salts"),
    "fluoride-anion": (-1, None, 1.0, 1.7712, 275.00, 8.0, "Table 2 F- row; Table 4 salts"),
    "chloride-anion": (-1, None, 1.0, 2.7560, 170.00, 8.0, "Table 2 Cl- row; Table 4 salts"),
    "bromide-anion": (-1, None, 1.0, 3.0707, 190.00, 8.0, "Table 2 Br- row; Table 4 salts"),
    "iodide-anion": (-1, None, 1.0, 3.6672, 200.00, 8.0, "Table 2 I- row; Table 4 salts"),
}
ASSOCIATION = {"water": (2425.7, 0.04509), "methanol": (2899.5, 0.03518), "ethanol": (2653.4, 0.03238)}
K_WATER = {"lithium-cation": -0.25, "sodium-cation": 0.0045, "potassium-cation": 0.1997, "fluoride-anion": 0.0,
           "chloride-anion": -0.25, "bromide-anion": -0.25, "iodide-anion": -0.25}  # Table 2
K_SALT = {  # Table 3, anion rows x cation columns
    "fluoride-anion": (0.0, 0.665, 1.000), "chloride-anion": (0.669, 0.317, 0.064),
    "bromide-anion": (0.591, 0.290, -0.102), "iodide-anion": (0.002, 0.018, -0.312)}
CATIONS = ("lithium-cation", "sodium-cation", "potassium-cation")
PACKING = 0.88  # Held 2014 Eq. (21), d_ion = 0.88 sigma_ion

# Variants: (electrolyte choice, permittivity choice, Born diameter factor, DH diameter factor)
VARIANTS = {
    "advanced": ("born", "component-permittivity-mixing", PACKING, PACKING),
    "revised": ("fully-dissociated-debye-huckel", "solvent-only", None, PACKING),
    "born-constant-permittivity": ("born", "solvent-only", PACKING, PACKING),
    "dh-mixed-permittivity": ("fully-dissociated-debye-huckel", "component-permittivity-mixing", None, PACKING),
    "advanced-dh-sigma": ("born", "component-permittivity-mixing", PACKING, 1.0),
    "revised-dh-sigma": ("fully-dissociated-debye-huckel", "solvent-only", None, 1.0),
    "born-diameter-sigma": ("born", "component-permittivity-mixing", 1.0, PACKING),
    "sodium-water-held2014": ("born", "component-permittivity-mixing", PACKING, PACKING),
}
SODIUM_WATER_HELD2014 = 0.00045485000000011766  # data/packets/held-2014-strategy2/1, T expression at 298.15 K


def value(magnitude, unit):
    return {"magnitude": magnitude, "unit": unit}


def provenance(locator, source=PAPER):
    return {"source_id": source, "locator": locator, "domain_id": DOMAIN}


def coefficient(identity, family, magnitude, unit, locator, source=PAPER):
    return {"identity": identity, "family": family, "value": value(magnitude, unit),
            "provenance": provenance(locator, source)}


def pair_k(a, b, variant):
    if "water" in (a, b):
        ion = b if a == "water" else a
        if ion == "sodium-cation" and variant == "sodium-water-held2014":
            return SODIUM_WATER_HELD2014, "Held 2014 Na+-water T expression at 298.15 K (sensitivity)", HELD
        return K_WATER.get(ion, 0.0), "Table 2 k_ion,water", PAPER
    anion = next((c for c in (a, b) if c in K_SALT), None)
    cation = next((c for c in (a, b) if c in CATIONS), None)
    if anion and cation:
        return K_SALT[anion][CATIONS.index(cation)], "Table 3 cation-anion k_ij", PAPER
    return 0.0, "Sec. 3.1: k_ij = 0 for ion-organic solvent and solvent-solvent pairs", PAPER


def document(variant="advanced"):
    electrolyte, permittivity, born_factor, dh_factor = VARIANTS[variant]
    components, correlations, sites, edges, pairs = [], [], [], [], []
    for cid, (charge, mass, m, sigma, u, eps, locator) in COMPONENTS.items():
        fixed = {"charge_number": {"value": value(charge, "elementary-charge"), "provenance": provenance(locator)}}
        if mass is not None:
            fixed["molar_mass"] = {"value": value(mass, "gram / mole"),
                                   "provenance": provenance("IUPAC standard atomic weights", "iupac-atomic-weights")}
        coefficients = [
            coefficient(f"component/{cid}/segment_count", "segment_count", m, "dimensionless", locator),
            coefficient(f"component/{cid}/dispersion_energy_over_k", "dispersion_energy_over_k", u, "kelvin", locator)]
        if electrolyte == "born":
            coefficients.append(coefficient(f"component/{cid}/solvation_factor", "solvation_factor", 1.0, "dimensionless",
                                            "Eq. (2) uncorrected Born term"))
        if sigma is None:
            correlations.append({
                "correlation_id": "component/water/segment_diameter/constant-plus-sum-of-exponentials",
                "component_id": "water", "family": "segment_diameter", "form": "constant-plus-sum-of-exponentials",
                "independent_variables": ["temperature"],
                "constant": {"identity": "component/water/segment_diameter/constant-plus-sum-of-exponentials/constant",
                             "value": value(2.7927, "angstrom")},
                "provenance": provenance("Table 1 footnote *"),
                "terms": [{"amplitude": {"identity": f"component/water/segment_diameter/constant-plus-sum-of-exponentials/term-{k}/amplitude",
                                         "value": value(a, "angstrom")},
                           "exponent_coefficient": {"identity": f"component/water/segment_diameter/constant-plus-sum-of-exponentials/term-{k}/exponent_coefficient",
                                                    "value": value(b, "1 / kelvin")}}
                          for k, (a, b) in enumerate(((10.11, -0.01775), (-1.417, -0.01146)))]})
        else:
            coefficients.append(coefficient(f"component/{cid}/segment_diameter", "segment_diameter", sigma, "angstrom", locator))
        if charge:
            coefficients += [
                coefficient(f"component/{cid}/packing_diameter", "packing_diameter", PACKING * sigma, "angstrom",
                            "Held 2014 Eq. (21) d_ion = 0.88 sigma_ion", HELD),
                coefficient(f"component/{cid}/debye_huckel_diameter", "debye_huckel_diameter", dh_factor * sigma, "angstrom",
                            "unified d_i = 0.88 sigma; reproduces the corrigendum revised curves" if dh_factor == PACKING
                            else "Held 2014 Sec. 2.1 closest approach sigma (sensitivity)", PAPER if dh_factor == PACKING else HELD)]
            if born_factor is not None:
                coefficients.append(coefficient(f"component/{cid}/born_diameter", "born_diameter", born_factor * sigma, "angstrom",
                                                "Sec. 2.1 Born a_i = ePC-SAFT ion diameter; 0.88 sigma reproduces Fig. 3"
                                                if born_factor == PACKING else "Born diameter sigma sensitivity"))
        if not charge or permittivity == "component-permittivity-mixing":
            coefficients.append(coefficient(f"component/{cid}/relative_permittivity", "relative_permittivity", eps,
                                            "dimensionless", locator))
        components.append({"component_id": cid, "aliases": [], "fixed": fixed, "coefficients": coefficients})
    for cid, (energy, volume) in ASSOCIATION.items():
        sites += [{"component_id": cid, "site_id": s, "site_role": role, "multiplicity": 1,
                   "provenance": provenance(f"Table 1 {cid} 2B scheme")} for s, role in (("a", "donor"), ("b", "acceptor"))]
        edges.append({"endpoint_a": {"component_id": cid, "site_id": "a"}, "endpoint_b": {"component_id": cid, "site_id": "b"},
                      "energy_over_k": {"identity": f"association/{cid}/a/{cid}/b/energy_over_k", "value": value(energy, "kelvin")},
                      "volume": {"identity": f"association/{cid}/a/{cid}/b/volume", "value": value(volume, "dimensionless")},
                      "source": {"kind": "explicit", "provenance": [provenance(f"Table 1 {cid} association")]}})
    for a, b in itertools.combinations(COMPONENTS, 2):
        if COMPONENTS[a][0] * COMPONENTS[b][0] > 0:
            continue  # same-sign ions take no dispersion pair
        k, locator, source = pair_k(a, b, variant)
        pairs.append({"component_id_a": a, "component_id_b": b, "coefficients": [
            coefficient(f"pair/{a}/{b}/k_ij", "k_ij", k, "dimensionless", locator, source),
            coefficient(f"pair/{a}/{b}/l_ij", "l_ij", 0.0, "dimensionless", "the paper uses no l_ij")]})
    families = [
        {"family_id": "model/association", "kind": "association", "choice": "general-site",
         "provenance": provenance("Table 1: all solvents 2B")},
        {"family_id": "model/base", "kind": "base", "choice": "pc-saft", "provenance": provenance("Eq. (11)")},
        {"family_id": "model/electrolyte", "kind": "electrolyte", "choice": electrolyte,
         "ionic_dispersion": "exclude-same-sign-ion-pairs",
         "provenance": provenance("Eq. (11); Held 2014 like-ion rule")},
        {"family_id": "model/relative_permittivity", "kind": "permittivity", "choice": permittivity,
         "provenance": provenance("Sec. 3.1 linear mole-fraction rule, Table 4" if permittivity != "solvent-only"
                                  else "salt-independent solvent permittivity (ePC-SAFT revised)")}]
    coefficients = []
    if electrolyte == "born":
        families[2].update(c_shell=0.0, c_dielectric=0.0)
        coefficients.append({"identity": "model/relative_permittivity/ionic_region_relative_permittivity",
                             "family": "ionic_region_relative_permittivity", "value": value(1.0, "dimensionless"),
                             "provenance": provenance("uncorrected Born term, Eq. (2)")})
    return {
        "schema": "epcsaft.parameters", "schema_version": 1,
        "document_id": f"bulow-2020-{variant}", "document_version": 1, "purpose": "external-document",
        "sources": [
            {"source_id": PAPER, "citation": "Bülow, M.; Ascani, M.; Held, C. ePC-SAFT advanced - Part I. Fluid Phase "
             "Equilibria 2021, 535, 112967; corrigendum 548 (2021) 113184.", "doi": "10.1016/j.fluid.2021.112967",
             "use_basis": "Tables 1-4 and Eqs. (2), (7)-(11) at 298.15 K and 1 bar"},
            {"source_id": HELD, "citation": "Held, C.; Reschke, T.; Mohammad, S.; Luza, A.; Sadowski, G. ePC-SAFT revised. "
             "Chem. Eng. Res. Des. 2014, 92, 2884-2897.", "doi": "10.1016/j.cherd.2014.05.017",
             "use_basis": "ion packing and Debye-Hückel diameter conventions of the parameter source"},
            {"source_id": "iupac-atomic-weights", "citation": "IUPAC standard atomic weights", "use_basis": "solvent molar masses"}],
        "domains": [{"domain_id": DOMAIN, "kind": "reported-conditions",
                     "temperature_min": value(298.15, "kelvin"), "temperature_max": value(298.15, "kelvin"),
                     "pressure_min": value(1, "bar"), "pressure_max": value(1, "bar")}],
        "components": components, "pairs": pairs, "model_families": families, "model_coefficients": coefficients,
        "correlations": correlations, "topology": {"presets": [], "sites": sites, "edges": edges, "no_cross_association": [
            {"component_id_a": a, "component_id_b": b, "provenance": provenance("no solvent cross edge: every state holds one solvent")}
            for a, b in itertools.combinations(ASSOCIATION, 2)]},
    }


if __name__ == "__main__":
    out = Path(__file__).resolve().parents[1] / "inputs/parameters/bulow-2020-advanced.json"
    out.write_text(json.dumps(document("advanced"), indent=1) + "\n")
