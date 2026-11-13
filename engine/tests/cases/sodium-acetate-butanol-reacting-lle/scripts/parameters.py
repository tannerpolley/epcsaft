"""Parameter declarations for water + 1-butanol + sodium acetate (+ 5-hydroxymethylfurfural, HMF) at 298.15 K and 1 bar.

The one owner of every value below: each carries its source id and locator, or the source id "assumed" when no
source prints it. Writes inputs/parameters-ternary.json, inputs/parameters-hmf.json and inputs/reactions.json.

Sources (Zotero keys):
  mohammad-2016-fpe   4ZANH7Q4 (SI docx 6FYN5ETB)  Mohammad et al., Fluid Phase Equilib. 428 (2016) 102
  held-2014           3DIW87Q5                      Held et al., ePC-SAFT revised, strategy 2 (Table 2, Table 3, Eq. 21)
  mohammad-2016-jpcb  KW39H7TE (PDF PISHG3WN)       Mohammad et al., J. Phys. Chem. B 120 (2016) 3797
  ascani-2022-esi     AGPHR9MB                      Ascani et al., Chem. Commun. 58 (2022) 8436, ESI Tables S1, S4
  iapws-r11-07        no Zotero item                IAPWS R11-07 (Bandura and Lvov 2006), as in the acetic-acid analyses
"""
import json
import math
from itertools import combinations
from pathlib import Path

T, PRESSURE = 298.15, 1.0e5
OUT = Path(__file__).resolve().parents[1] / "inputs"
DOMAIN = "mohammad-2016-298-15-k-1-bar"
SOURCE_USE = {
    "mohammad-2016-fpe": ("Mohammad, Grundl, Muller, Kunz, Sadowski, Held. Fluid Phase Equilib. 428 (2016) 102; SI Tables S1, S3. Zotero 4ZANH7Q4, SI 6FYN5ETB",
                          "10.1016/j.fluid.2016.05.001", "Tables 2-3, Eqs. 4-9, Sec. 3 permittivity"),
    "held-2014": ("Held, Reschke, Mohammad, Luza, Sadowski. ePC-SAFT revised. Chem. Eng. Res. Des. 92 (2014) 2884. Zotero 3DIW87Q5",
                  "10.1016/j.cherd.2014.05.017", "strategy 2 ion parameters, Tables 2-3, Eq. 21"),
    "mohammad-2016-jpcb": ("Mohammad, Held, Altuntepe, Kose, Sadowski. J. Phys. Chem. B 120 (2016) 3797. Zotero KW39H7TE, PDF PISHG3WN",
                           "10.1021/acs.jpcb.5b11588", "Eqs. 5-6 water/HMF cross association; Table 3 HMF binaries"),
    "ascani-2022-esi": ("Ascani, Pabsch, Klinksiek, Gajardo-Parra, Sadowski, Held. Chem. Commun. 58 (2022) 8436 and ESI. Zotero AGPHR9MB",
                        "10.1039/D2CC02943J", "acetic-acid Table S1 parameters, water pair, pKa Table S4"),
    "iupac-atomic-weights": ("IUPAC standard atomic weights", None, "molar masses"),
    "assumed": ("No source prints this value; declared by the analysis", None,
                "unprinted pair k_ij = l_ij = 0 and unprinted cross-association rules"),
}

# component: charge, molar mass g/mol, m, sigma/A, u/k K, (association eps/k K, kappa) or None, source, locator
FPE, HELD, ESI, JPCB = "mohammad-2016-fpe", "held-2014", "ascani-2022-esi", "mohammad-2016-jpcb"
COMPONENTS = {
    "water": (0, 18.01528, 1.2047, None, 353.95, (2425.67, 0.0451), FPE,
              "Table 2 water row; sigma(T) = 2.7927 + 10.11 exp(-0.01775 T) - 1.417 exp(-0.01146 T), footnote a"),
    "1-butanol": (0, 74.1216, 2.7515, 3.6139, 259.59, (2544.7, 0.0067), FPE, "Table 2 1-butanol row"),
    "hmf": (0, 126.11004, 2.310, 4.038, 320.38, (3167.61, 0.0010), FPE,
            "Table 2 HMF row with sigma 4.038 and m 2.310 (footnote c; Table 2 of the J. Phys. Chem. B paper swaps them)"),
    "sodium-cation": (1, 22.98977, 1.0, 2.8232, 230.0, None, HELD, "Table 2 strategy 2 Na+ row"),
    "acetate-anion": (-1, 59.04402, 1.0, 3.9328, 150.0, None, HELD, "Table 2 strategy 2 Ac- row"),
    "acetic-acid": (0, 60.05196, 1.3403, 3.8582, 211.59, (3044.4, 0.07555), ESI, "Table S1 acetic acid row"),
    "hydrogen-cation": (1, 1.00794, 1.0, 3.4654, 500.0, None, HELD, "Table 2 strategy 2 H3O+ row, used as H+"),
    "hydroxide-anion": (-1, 17.00734, 1.0, 2.0177, 650.0, None, HELD, "Table 2 strategy 2 OH- row"),
}
SPECIES = {"ternary": ("water", "1-butanol", "sodium-cation", "acetate-anion", "acetic-acid", "hydrogen-cation", "hydroxide-anion"),
           "hmf": ("water", "1-butanol", "hmf", "sodium-cation", "acetate-anion", "acetic-acid", "hydrogen-cation", "hydroxide-anion")}
EPS = 78.385  # Mohammad 2016 Sec. 3: one relative permittivity for every mixture
ION_DIAMETER_FACTOR = 0.88  # Held 2014 Eq. 21: packing diameter; the Debye-Hueckel diameter is the sensitivity
K_HB = 0.026  # Mohammad 2016 Eqs. 8-9, water/1-butanol cross association

# (a, b): (k_ij, l_ij, source, locator); a pair not listed is unprinted and takes k_ij = l_ij = 0
PAIRS = {
    ("water", "1-butanol"): (2.94e-4 * (T - 293.15) - 0.01436, -0.0044, FPE, "Table 3 k_ij = 2.94e-4 (T - 293.15 K) - 0.01436 at 298.15 K"),
    ("1-butanol", "sodium-cation"): (0.22, 0.245, FPE, "Table 3"),
    ("1-butanol", "acetate-anion"): (0.45, 0.0, FPE, "Table 3"),
    ("1-butanol", "hmf"): (0.006, 0.003, FPE, "Table 3"),
    ("hmf", "sodium-cation"): (0.92, -0.075, FPE, "Table 3"),
    ("hmf", "acetate-anion"): (-0.19, -0.060, FPE, "Table 3"),
    ("hmf", "water"): (-3.832e-3 * (T - 298.15) - 0.042, 0.0, JPCB, "Table 3 k_ij = -3.832e-3 (T - 298.15 K) - 0.042 at 298.15 K"),
    ("water", "sodium-cation"): (-0.007981 * T + 2.37999, 0.0, HELD, "Table 2 footnote a, at 298.15 K"),
    ("water", "acetate-anion"): (-0.23, 0.0, HELD, "Table 2 strategy 2 k_ion-water"),
    ("water", "hydrogen-cation"): (0.25, 0.0, HELD, "Table 2 strategy 2 k_ion-water (H3O+)"),
    ("water", "hydroxide-anion"): (-0.25, 0.0, HELD, "Table 2 strategy 2 k_ion-water"),
    ("sodium-cation", "acetate-anion"): (0.246, 0.0, HELD, "Table 3"),
    ("sodium-cation", "hydroxide-anion"): (0.649, 0.0, HELD, "Table 3"),
    ("water", "acetic-acid"): (-0.1247, 0.0, ESI, "Table S3 water/acetic acid"),
}
# Wolbach-Sandler cross association without k_hb unless stated (a, b): (source, locator)
CROSS = {
    ("water", "acetic-acid"): ("assumed", "Wolbach-Sandler without k_hb; the Ascani ESI states no rule for this pair"),
    ("1-butanol", "acetic-acid"): ("assumed", "Wolbach-Sandler without k_hb; no source prints this pair"),
    ("water", "hmf"): (JPCB, "Eqs. 5-6 Wolbach-Sandler, no k_hb given for this pair"),
    ("1-butanol", "hmf"): ("assumed", "Wolbach-Sandler without k_hb; Mohammad 2016 prints the rule only for water/1-butanol and water/HMF"),
    ("acetic-acid", "hmf"): ("assumed", "Wolbach-Sandler without k_hb; no source prints this pair"),
}
REACTIONS = {
    "HAc": {"stoichiometry": {"acetic-acid": -1, "hydrogen-cation": 1, "acetate-anion": 1}, "pK": 4.76,
            "source": ESI, "locator": "Table S4 pKa of acetic acid at 298.15 K"},
    "H2O": {"stoichiometry": {"water": -1, "hydrogen-cation": 1, "hydroxide-anion": 1}, "pK": 13.9945,
            "source": "iapws-r11-07", "locator": "IAPWS R11-07 (Bandura and Lvov 2006) pKw at 298.15 K and 0.1 MPa"},
}


def value(magnitude, unit="dimensionless"):
    return {"magnitude": magnitude, "unit": unit}


def provenance(source, locator):
    return {"source_id": source, "locator": locator, "domain_id": DOMAIN}


def coefficient(cid, family, magnitude, unit, source, locator):
    return {"identity": f"component/{cid}/{family}", "family": family, "value": value(magnitude, unit),
            "provenance": provenance(source, locator)}


def sigma_water(t=T):
    return 2.7927 + 10.11 * math.exp(-0.01775 * t) - 1.417 * math.exp(-0.01146 * t)


def unprinted_pairs(system):
    """Opposite-charge or neutral pairs of the system that no source prints."""
    species = SPECIES[system]
    return [(a, b) for a, b in combinations(species, 2)
            if COMPONENTS[a][0] * COMPONENTS[b][0] <= 0 and (a, b) not in PAIRS and (b, a) not in PAIRS]


def combining(other, partner, source, locator):
    ids = [f"association/{other}/a/{other}/b/energy_over_k", f"association/{partner}/a/{partner}/b/energy_over_k"]
    vols = [f"association/{other}/a/{other}/b/volume", f"association/{partner}/a/{partner}/b/volume"]
    return [{"endpoint_a": {"component_id": other, "site_id": a}, "endpoint_b": {"component_id": partner, "site_id": b},
             "source": {"kind": "combining-rule", "rule_id": "wolbach-sandler",
                        "inputs": {"energy_over_k": ids, "volume": vols}, "provenance": [provenance(source, locator)]}}
            for a, b in (("a", "b"), ("b", "a"))]


def document(system, species=None, dh_factor=ION_DIAMETER_FACTOR, kij=None, hmf_cross_off=()):
    """Parameter document of `system` restricted to `species`.

    dh_factor: Debye-Hueckel diameter over sigma (the packing diameter stays 0.88 sigma).
    kij: {(a, b): k_ij} replacing a pair's k_ij (sensitivity); hmf_cross_off: partners whose HMF cross association is removed.
    """
    species = species or SPECIES[system]
    kij = kij or {}
    components, correlations, sites, edges, pairs = [], [], [], [], []
    for cid in species:
        z, mass, m, sigma, u, assoc, source, locator = COMPONENTS[cid]
        fixed = {"charge_number": {"value": value(z, "elementary-charge"), "provenance": provenance(source, locator)},
                 "molar_mass": {"value": value(mass, "gram / mole"), "provenance": provenance("iupac-atomic-weights", "standard atomic weights")}}
        coefficients = [coefficient(cid, "segment_count", m, "dimensionless", source, locator),
                        coefficient(cid, "dispersion_energy_over_k", u, "kelvin", source, locator)]
        if cid == "water":
            terms = ((10.11, -0.01775), (-1.417, -0.01146))
            correlations.append({
                "correlation_id": "component/water/segment_diameter/constant-plus-sum-of-exponentials",
                "component_id": "water", "family": "segment_diameter", "form": "constant-plus-sum-of-exponentials",
                "independent_variables": ["temperature"], "provenance": provenance(source, "Table 2 footnote a"),
                "constant": {"identity": "component/water/segment_diameter/constant-plus-sum-of-exponentials/constant", "value": value(2.7927, "angstrom")},
                "terms": [{"amplitude": {"identity": f"component/water/segment_diameter/constant-plus-sum-of-exponentials/term-{k}/amplitude", "value": value(a, "angstrom")},
                           "exponent_coefficient": {"identity": f"component/water/segment_diameter/constant-plus-sum-of-exponentials/term-{k}/exponent_coefficient", "value": value(b, "1 / kelvin")}}
                          for k, (a, b) in enumerate(terms)]})
        else:
            coefficients.append(coefficient(cid, "segment_diameter", sigma, "angstrom", source, locator))
        if z:
            coefficients.append(coefficient(cid, "packing_diameter", ION_DIAMETER_FACTOR * sigma, "angstrom", HELD, "Eq. 21 packing diameter 0.88 sigma"))
            coefficients.append(coefficient(cid, "debye_huckel_diameter", dh_factor * sigma, "angstrom", HELD,
                                            "Eq. 21 Debye-Hueckel diameter 0.88 sigma" if dh_factor == ION_DIAMETER_FACTOR else f"sensitivity: {dh_factor} sigma"))
        else:
            correlations.append({"correlation_id": f"component/{cid}/relative_permittivity/constant", "component_id": cid,
                                 "family": "relative_permittivity", "form": "constant", "independent_variables": [],
                                 "constant": {"identity": f"component/{cid}/relative_permittivity/constant/constant", "value": value(EPS)},
                                 "provenance": provenance(FPE, "Sec. 3: 78.385 for every mixture")})
        if assoc:
            sites += [{"component_id": cid, "site_id": s, "site_role": role, "multiplicity": 1, "provenance": provenance(source, locator)}
                      for s, role in (("a", "donor"), ("b", "acceptor"))]
            edges.append({"endpoint_a": {"component_id": cid, "site_id": "a"}, "endpoint_b": {"component_id": cid, "site_id": "b"},
                          "energy_over_k": {"identity": f"association/{cid}/a/{cid}/b/energy_over_k", "value": value(assoc[0], "kelvin")},
                          "volume": {"identity": f"association/{cid}/a/{cid}/b/volume", "value": value(assoc[1])},
                          "source": {"kind": "explicit", "provenance": [provenance(source, locator)]}})
        components.append({"component_id": cid, "aliases": [], "fixed": fixed, "coefficients": coefficients})
    # water/1-butanol cross association: Mohammad Eqs. 8-9, Wolbach-Sandler with k_hb, explicit at 298.15 K
    (e_w, k_w), (e_b, k_b) = COMPONENTS["water"][5], COMPONENTS["1-butanol"][5]
    s_w, s_b = sigma_water(), COMPONENTS["1-butanol"][3]
    energy = (e_w + e_b) / 2 * (1 - K_HB)
    volume = math.sqrt(k_w * k_b) * (math.sqrt(s_w * s_b) / (0.5 * (s_w + s_b))) ** 3
    for a, b in (("a", "b"), ("b", "a")):
        edges.append({"endpoint_a": {"component_id": "1-butanol", "site_id": a}, "endpoint_b": {"component_id": "water", "site_id": b},
                      "energy_over_k": {"identity": f"association/1-butanol/{a}/water/{b}/energy_over_k", "value": value(energy, "kelvin")},
                      "volume": {"identity": f"association/1-butanol/{a}/water/{b}/volume", "value": value(volume)},
                      "source": {"kind": "explicit", "provenance": [provenance(FPE, "Eqs. 8-9 Wolbach-Sandler with k_hb = 0.026, evaluated at 298.15 K")]}})
    no_cross = []
    for (a, b), (source, locator) in CROSS.items():
        if a in species and b in species:
            if "hmf" in (a, b) and (a if b == "hmf" else b) in hmf_cross_off:
                no_cross.append({"component_id_a": a, "component_id_b": b, "provenance": [provenance("assumed", "switched off for the sensitivity")]})
            else:
                edges += combining(a, b, source, locator)
    for a, b in combinations(species, 2):
        if COMPONENTS[a][0] * COMPONENTS[b][0] > 0:
            continue  # like-charged ions take no dispersion pair (exclude-same-sign-ion-pairs)
        k, l, source, locator = PAIRS.get((a, b)) or PAIRS.get((b, a)) or (0.0, 0.0, "assumed", "unprinted pair, k_ij = l_ij = 0")
        k = kij.get((a, b), kij.get((b, a), k))
        pairs.append({"component_id_a": a, "component_id_b": b, "coefficients": [
            {"identity": f"pair/{a}/{b}/k_ij", "family": "k_ij", "value": value(k), "provenance": provenance(source, locator)},
            {"identity": f"pair/{a}/{b}/l_ij", "family": "l_ij", "value": value(l), "provenance": provenance(source, locator)}]})
    topology = {"presets": [], "sites": sites, "edges": edges}
    if no_cross:
        topology["no_cross_association"] = no_cross
    return {
        "schema": "epcsaft.parameters", "schema_version": 1, "document_id": f"sodium-acetate-butanol-{system}", "document_version": 1,
        "purpose": "external-document",
        "sources": [{"source_id": sid, "citation": cite, **({"doi": doi} if doi else {}), "use_basis": use}
                    for sid, (cite, doi, use) in SOURCE_USE.items()],
        "domains": [{"domain_id": DOMAIN, "kind": "reported-conditions", "temperature_min": value(T, "kelvin"), "temperature_max": value(T, "kelvin"),
                     "pressure_min": value(1, "bar"), "pressure_max": value(1, "bar")}],
        "components": components, "pairs": pairs,
        "model_families": [
            {"family_id": "model/association", "kind": "association", "choice": "general-site", "provenance": provenance(FPE, "Table 2 2B association")},
            {"family_id": "model/base", "kind": "base", "choice": "pc-saft", "provenance": provenance(FPE, "PC-SAFT hard chain and dispersion")},
            {"family_id": "model/electrolyte", "kind": "electrolyte", "choice": "fully-dissociated-debye-huckel",
             "ionic_dispersion": "exclude-same-sign-ion-pairs", "provenance": provenance(HELD, "Debye-Hueckel term without Born term; like-ion rule")},
            {"family_id": "model/relative_permittivity", "kind": "permittivity", "choice": "solvent-only", "provenance": provenance(FPE, "Sec. 3: constant 78.385")}],
        "model_coefficients": [], "correlations": correlations, "topology": topology}


def reactions_document():
    return {"temperature_K": T, "pressure_Pa": PRESSURE, "reactions": REACTIONS,
            "basis": "molality, infinite dilution in pure water, standard molality 1 mol/kg (CommonMolalityInfiniteDilution, as in the acetic-acid analyses)"}


if __name__ == "__main__":
    for system in SPECIES:
        (OUT / f"parameters-{system}.json").write_text(json.dumps(document(system), indent=1) + "\n")
    (OUT / "reactions.json").write_text(json.dumps(reactions_document(), indent=1) + "\n")
