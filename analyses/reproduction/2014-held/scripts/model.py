"""Declare Held 2014 strategy-1/strategy-2 Parameters documents for water + one salt.

Water, model families, the water sigma(T) correlation and the 2B association topology are
taken from data/packets/held-2014-strategy2/1 (Held 2014 Table 1). Ion rows come from the
retained Table 2/3 transcription. Strategy 1 uses no water-ion k_ij and no cation-anion
dispersion (k_ij = 1); strategy 2 uses Table 2 water-ion and Table 3 cation-anion k_ij.
"""
import copy
import math
import csv
import json
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
REPO = HERE.parents[2]
PACKET = REPO / "data/packets/held-2014-strategy2/1/parameters/parameters.json"
IONS = {r["ion"]: r for r in csv.DictReader((HERE / "inputs/held-2014-ions.csv").open())}
TABLE3 = {r["anion"]: r for r in csv.DictReader((HERE / "inputs/held-2014-table-3.csv").open())}
NEUTRALS = {r["component_id"]: r for r in csv.DictReader((HERE / "inputs/held-2014-table-1.csv").open())}
TABLE4 = {(r["ion"], r["solvent"]): (float(r["k_ij"]), float(r["l_ij"]))
          for r in csv.DictReader((HERE / "inputs/held-2014-table-4.csv").open())}
SOURCE = {"source_id": "held-2014", "locator": "Tables 2-3; analyses/reproduction/2014-held/inputs",
          "domain_id": "held-2014-reproduction"}


def water_ion_kij(ion, T):
    """Table 2 k_ion-water; footnotes a/b are linear in T. Blank (F-) is taken as 0."""
    if ion == "Na+":
        return -0.007981 * T + 2.37999
    if ion == "K+":
        return -0.004012 * T + 1.3959
    return float(IONS[ion]["k_water_strategy2"] or 0.0)


def _value(x, unit):
    return {"magnitude": x, "unit": unit}


def _ion(ion, strategy, dh_factor):
    row = IONS[ion]
    sigma = float(row[f"sigma_strategy{strategy}_A"])
    u = float(row[f"u_strategy{strategy}_K"])
    cid = row["component_id"]
    coefficients = [("segment_count", 1.0, "dimensionless"), ("segment_diameter", sigma, "angstrom"),
                    ("dispersion_energy_over_k", u, "kelvin"), ("debye_huckel_diameter", dh_factor * sigma, "angstrom"),
                    ("packing_diameter", 0.88 * sigma, "angstrom")]  # Eq. (21)
    return {"component_id": cid, "aliases": [ion],
            "fixed": {"charge_number": {"value": _value(int(row["charge"]), "elementary-charge"), "provenance": SOURCE}},
            "coefficients": [{"identity": f"component/{cid}/{family}", "family": family,
                              "value": _value(x, unit), "provenance": SOURCE} for family, x, unit in coefficients]}


def _pair(a, b, k):
    return {"component_id_a": a, "component_id_b": b,
            "coefficients": [{"identity": f"pair/{a}/{b}/k_ij", "family": "k_ij",
                              "value": _value(k, "dimensionless"), "provenance": SOURCE}]}


def water_sigma(T):
    """Water sigma(T) in angstrom, Held 2008 Table 1 / Eq. (7) form (Fuchs 2006).

    Held 2014 Table 1 footnote a prints +0.01146 T in the last exponent, which gives sigma_w = -40.3
    angstrom at 298.15 K; Held 2008 Table 1 and Eq. (7) give -0.01146 T (2.797 angstrom). The Held 2008
    sign is a source-backed correction, not verbatim Held 2014 input.
    """
    return 2.7927 + 10.11 * math.exp(-0.01775 * T) - 1.417 * math.exp(-0.01146 * T)


def _neutral(row):
    cid = row["component_id"]
    coefficients = [("segment_count", float(row["segment_count"]), "dimensionless"),
                    ("segment_diameter", float(row["sigma_A"]), "angstrom"),
                    ("dispersion_energy_over_k", float(row["u_K"]), "kelvin")]
    return {"component_id": cid, "aliases": [],
            "fixed": {"molar_mass": {"value": _value(float(row["molar_mass_g_mol"]), "gram / mole"), "provenance": SOURCE},
                      "charge_number": {"value": _value(0, "elementary-charge"), "provenance": SOURCE}},
            "coefficients": [{"identity": f"component/{cid}/{family}", "family": family,
                              "value": _value(x, unit), "provenance": SOURCE} for family, x, unit in coefficients]}


def _edge(a, site_a, b, site_b, energy, volume):
    ident = f"association/{a}/{site_a}/{b}/{site_b}"
    return {"endpoint_a": {"component_id": a, "site_id": site_a}, "endpoint_b": {"component_id": b, "site_id": site_b},
            "energy_over_k": {"identity": f"{ident}/energy_over_k", "value": _value(energy, "kelvin")},
            "volume": {"identity": f"{ident}/volume", "value": _value(volume, "dimensionless")},
            "source": {"kind": "explicit", "provenance": [SOURCE]}}


def document(strategy, cation, anion, T, water_permittivity=78.09, dh_factor=0.88, neutral=None):
    """Parameters mapping for water [+ one Table 1 neutral] + cation + anion at temperature T.

    Water-ion and water-neutral k_ij are evaluated at T. dh_factor sets the Debye-Huckel closest
    approach as a multiple of sigma_ion: 0.88 is the Eq. (21) diameter d_ion that reproduces the
    published curves; 1.0 is the Section 2.1 wording. Cross association with water follows
    Eqs. (17)-(18) (Wolbach-Sandler) with the Table 1 k_hb, as explicit edges.
    """
    base = json.loads(PACKET.read_text())
    doc = copy.deepcopy(base)
    doc["document_id"] = f"held-2014-strategy{strategy}-{cation}-{anion}"
    doc["domains"].append({"domain_id": "held-2014-reproduction", "kind": "reported-conditions",
                           "temperature_min": _value(273.15, "kelvin"), "temperature_max": _value(298.15, "kelvin"),
                           "pressure_min": _value(1, "bar"), "pressure_max": _value(1, "bar")})
    doc["components"] = [c for c in base["components"] if c["component_id"] == "water"]
    doc["components"] += [_ion(cation, strategy, dh_factor), _ion(anion, strategy, dh_factor)]
    for item in doc["correlations"]:
        if item["family"] == "relative_permittivity":
            item["constant"]["value"]["magnitude"] = water_permittivity
    c, a = IONS[cation]["component_id"], IONS[anion]["component_id"]
    if strategy == 1:
        doc["pairs"] = [_pair(c, "water", 0.0), _pair(a, "water", 0.0), _pair(c, a, 1.0)]
    else:
        doc["pairs"] = [_pair(c, "water", water_ion_kij(cation, T)), _pair(a, "water", water_ion_kij(anion, T)),
                        _pair(c, a, float(TABLE3[anion][cation]))]
    if neutral is None:
        return doc, ["water", c, a]
    row = NEUTRALS[neutral]
    doc["components"].insert(1, _neutral(row))
    k = float(row["kij_water_constant"]) + float(row["kij_water_slope_per_K"]) * (T - float(row["kij_water_reference_K"]))
    water_pair = _pair(neutral, "water", k)
    water_pair["coefficients"].append({"identity": f"pair/{neutral}/water/l_ij", "family": "l_ij",
                                       "value": _value(float(row["lij_water"]), "dimensionless"), "provenance": SOURCE})
    doc["pairs"].append(water_pair)
    for ion, cid in ((cation, c), (anion, a)):  # Table 4; unlisted ion-solvent pairs are 0
        k_ion, l_ion = TABLE4.get((ion, neutral), (0.0, 0.0))
        pair = _pair(cid, neutral, k_ion)
        pair["coefficients"].append({"identity": f"pair/{cid}/{neutral}/l_ij", "family": "l_ij",
                                     "value": _value(l_ion, "dimensionless"), "provenance": SOURCE})
        doc["pairs"].append(pair)
    if row["relative_permittivity"]:
        water_eps = next(item for item in doc["correlations"] if item["family"] == "relative_permittivity")
        eps = copy.deepcopy(water_eps)
        eps["correlation_id"] = eps["constant"]["identity"] = f"component/{neutral}/relative_permittivity/constant"
        eps["component_id"], eps["provenance"] = neutral, SOURCE
        eps["constant"]["value"]["magnitude"] = float(row["relative_permittivity"])
        doc["correlations"].append(eps)
    if row["assoc_energy_K"]:
        energy, volume = float(row["assoc_energy_K"]), float(row["assoc_volume"])
        topology = doc["topology"]
        topology["sites"] += [{"component_id": neutral, "site_id": site, "site_role": role, "multiplicity": 1,
                               "provenance": SOURCE} for site, role in (("a", "donor"), ("b", "acceptor"))]
        water = next(e for e in topology["edges"] if e["endpoint_a"]["component_id"] == "water")
        e_w, v_w = water["energy_over_k"]["value"]["magnitude"], water["volume"]["value"]["magnitude"]
        s_w, s_n = water_sigma(T), float(row["sigma_A"])
        cross_e = 0.5 * (e_w + energy) * (1 - float(row["khb_water"]))
        cross_v = math.sqrt(v_w * volume) * (math.sqrt(s_w * s_n) / (0.5 * (s_w + s_n))) ** 3
        topology["edges"] += [_edge(neutral, "a", neutral, "b", energy, volume),
                              _edge(neutral, "a", "water", "b", cross_e, cross_v),
                              _edge(neutral, "b", "water", "a", cross_e, cross_v)]
    return doc, ["water", neutral, c, a]
