"""Write the Ascani 2023 PC-SAFT parameter documents (Tables 1, 2, 4) for systems 1 and 2.

Each system is evaluated at its one Table 4 temperature, so k_ij(T) of Eq. A8 is declared at that T.
Cross association follows Eq. A6 (arithmetic energy) with two volume rules:
  printed          Eq. A7 as printed: kappa_ij = sqrt(kappa_i kappa_j) [sqrt(s_i s_j) / ((s_i + s_j)/2)]^2,
                   declared as explicit edges with water sigma(T) evaluated at the system temperature;
  wolbach-sandler  the Engine's Gross 2002 rule, the same bracket cubed and evaluated at the state.
Esters carry two sites with zero self-association energy (Table 1), so they associate only by cross edges.
"""
import csv
import itertools
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "results/tables"
SRC = "ascani-2023"
IDS = {"Water": "water", "Acetic acid": "acetic-acid", "1-Pentanol": "1-pentanol", "1-Hexanol": "1-hexanol",
       "Pentyl Acetate": "pentyl-acetate", "Hexyl Acetate": "hexyl-acetate"}
MASS = {"water": 18.01528, "acetic-acid": 60.05196, "1-pentanol": 88.14818, "1-hexanol": 102.17476,
        "pentyl-acetate": 130.18486, "hexyl-acetate": 144.21144}  # IUPAC atomic weights, g/mol
SYSTEMS = {1: ("acetic-acid", "1-pentanol", "pentyl-acetate", "water"),
           2: ("acetic-acid", "1-hexanol", "hexyl-acetate", "water")}
WATER_SIGMA = (2.7927, ((10.11, -0.01775), (-1.417, -0.01146)))  # Table 1 footnote


def rows(table):
    with (TABLES / table / "source" / f"{table}.csv").open() as stream:
        return list(csv.reader(stream))[1:]


def temperatures():
    return {int(r[0]): float(r[1]) for r in rows("table-004")}


def ka():
    return {int(r[0]): float(r[3]) for r in rows("table-004")}


def pure():
    """component_id: (m, sigma or None, u/k, eps_AB/k, kappa_AB) from Table 1."""
    return {IDS[r[0]]: (float(r[1]), None if r[2] == "*" else float(r[2]), float(r[3]), float(r[5]), float(r[6]))
            for r in rows("table-001") if r[0] in IDS}


def k_ij(T):
    """Eq. A8: k_ij(T) = k_ij,298.15 + k_ij,T (T - 298.15)."""
    return {(IDS[r[0]], IDS[r[1]]): float(r[2]) + (0.0 if r[3] == "-" else float(r[3])) * (T - 298.15)
            for r in rows("table-002")}


def water_sigma(T):
    c, terms = WATER_SIGMA
    return c + sum(a * math.exp(b * T) for a, b in terms)


def value(magnitude, unit):
    return {"magnitude": magnitude, "unit": unit}


def provenance(locator, domain):
    return {"source_id": SRC, "locator": locator, "domain_id": domain}


def document(system, rule):
    T = temperatures()[system]
    domain = f"system-{system}-{T}-K-1-bar"
    table1, kij, species = pure(), k_ij(T), SYSTEMS[system]
    components, correlations, sites, edges, pairs = [], [], [], [], []
    for cid in species:
        m, sigma, u, _, _ = table1[cid]
        coefficients = [
            {"identity": f"component/{cid}/segment_count", "family": "segment_count", "value": value(m, "dimensionless"),
             "provenance": provenance(f"Table 1 {cid} row", domain)},
            {"identity": f"component/{cid}/dispersion_energy_over_k", "family": "dispersion_energy_over_k",
             "value": value(u, "kelvin"), "provenance": provenance(f"Table 1 {cid} row", domain)}]
        if sigma is None:
            c, terms = WATER_SIGMA
            stem = "component/water/segment_diameter/constant-plus-sum-of-exponentials"
            correlations.append({
                "correlation_id": stem, "component_id": cid, "family": "segment_diameter",
                "form": "constant-plus-sum-of-exponentials", "independent_variables": ["temperature"],
                "constant": {"identity": f"{stem}/constant", "value": value(c, "angstrom")},
                "provenance": provenance("Table 1 water sigma footnote", domain),
                "terms": [{"amplitude": {"identity": f"{stem}/term-{k}/amplitude", "value": value(a, "angstrom")},
                           "exponent_coefficient": {"identity": f"{stem}/term-{k}/exponent_coefficient",
                                                    "value": value(b, "1 / kelvin")}}
                          for k, (a, b) in enumerate(terms)]})
        else:
            coefficients.append({"identity": f"component/{cid}/segment_diameter", "family": "segment_diameter",
                                 "value": value(sigma, "angstrom"), "provenance": provenance(f"Table 1 {cid} row", domain)})
        components.append({"component_id": cid, "aliases": [], "coefficients": coefficients, "fixed": {
            "charge_number": {"value": value(0, "elementary-charge"), "provenance": provenance("neutral", domain)},
            "molar_mass": {"value": value(MASS[cid], "gram / mole"), "provenance": provenance("IUPAC atomic weights", domain)}}})
    for cid in species:
        _, _, _, energy, volume = table1[cid]
        sites += [{"component_id": cid, "site_id": s, "site_role": role, "multiplicity": 1,
                   "provenance": provenance(f"Table 1 {cid} N_i = 2", domain)} for s, role in (("a", "donor"), ("b", "acceptor"))]
        edges.append({"endpoint_a": {"component_id": cid, "site_id": "a"}, "endpoint_b": {"component_id": cid, "site_id": "b"},
                      "energy_over_k": {"identity": f"association/{cid}/a/{cid}/b/energy_over_k", "value": value(energy, "kelvin")},
                      "volume": {"identity": f"association/{cid}/a/{cid}/b/volume", "value": value(volume, "dimensionless")},
                      "source": {"kind": "explicit", "provenance": [provenance(f"Table 1 {cid} association", domain)]}})
    sigma_at_T = {cid: water_sigma(T) if table1[cid][1] is None else table1[cid][1] for cid in species}
    for i, j in itertools.combinations(species, 2):
        for a_site, b_site in (("a", "b"), ("b", "a")):
            edge = {"endpoint_a": {"component_id": i, "site_id": a_site}, "endpoint_b": {"component_id": j, "site_id": b_site}}
            if rule == "printed":
                si, sj = sigma_at_T[i], sigma_at_T[j]
                kappa = math.sqrt(table1[i][4] * table1[j][4]) * (math.sqrt(si * sj) / (0.5 * (si + sj))) ** 2
                identity = f"association/{i}/{a_site}/{j}/{b_site}"
                edge |= {"energy_over_k": {"identity": f"{identity}/energy_over_k",
                                           "value": value(0.5 * (table1[i][3] + table1[j][3]), "kelvin")},
                         "volume": {"identity": f"{identity}/volume", "value": value(kappa, "dimensionless")},
                         "source": {"kind": "explicit", "provenance": [provenance(
                             f"Eqs. A6-A7 as printed (exponent 2), water sigma at {T} K", domain)]}}
            else:
                edge["source"] = {"kind": "combining-rule", "rule_id": "wolbach-sandler", "inputs": {
                    "energy_over_k": [f"association/{c}/a/{c}/b/energy_over_k" for c in (i, j)],
                    "volume": [f"association/{c}/a/{c}/b/volume" for c in (i, j)]},
                    "provenance": [provenance("Eq. A6; Eq. A7 bracket cubed (Gross 2002 Eqs. 2-3)", domain)]}
            edges.append(edge)
        k = kij.get((i, j), kij.get((j, i)))
        pairs.append({"component_id_a": i, "component_id_b": j, "coefficients": [
            {"identity": f"pair/{i}/{j}/k_ij", "family": "k_ij", "value": value(k, "dimensionless"),
             "provenance": provenance(f"Table 2, Eq. A8 at {T} K", domain)},
            {"identity": f"pair/{i}/{j}/l_ij", "family": "l_ij", "value": value(0.0, "dimensionless"),
             "provenance": provenance("Eq. A4 has no l_ij", domain)}]})
    families = [("association", "general-site", "Table 1 two-site association"), ("base", "pc-saft", "Eq. A3"),
                ("electrolyte", "none", "neutral components"), ("permittivity", "none", "neutral components")]
    return {
        "schema": "epcsaft.parameters", "schema_version": 1, "document_id": f"ascani-2023-system-{system}-{rule}",
        "document_version": 1, "purpose": "external-document",
        "sources": [{"source_id": SRC, "citation": "Ascani, M.; Sadowski, G.; Held, C. Simultaneous Predictions of Chemical "
                     "and Phase Equilibria in Systems with an Esterification Reaction Using PC-SAFT. Molecules 2023, 28, 1768.",
                     "doi": "10.3390/molecules28041768", "use_basis": "Tables 1, 2 and 4; Appendix A Eqs. A4-A8"}],
        "domains": [{"domain_id": domain, "kind": "reported-conditions", "temperature_min": value(T, "kelvin"),
                     "temperature_max": value(T, "kelvin"), "pressure_min": value(1, "bar"), "pressure_max": value(1, "bar")}],
        "components": components, "pairs": pairs,
        "model_families": [{"family_id": f"model/{kind}", "kind": kind, "choice": choice,
                            "provenance": provenance(locator, domain)} for kind, choice, locator in families],
        "model_coefficients": [], "correlations": correlations,
        "topology": {"presets": [], "sites": sites, "edges": edges},
    }


if __name__ == "__main__":
    (ROOT / "inputs/parameters").mkdir(parents=True, exist_ok=True)
    for system in SYSTEMS:
        for rule in ("printed", "wolbach-sandler"):
            (ROOT / f"inputs/parameters/system-{system}-{rule}.json").write_text(json.dumps(document(system, rule), indent=1) + "\n")
