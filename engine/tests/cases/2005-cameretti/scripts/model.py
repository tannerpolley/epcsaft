"""Cameretti, Sadowski and Mollerup (2005) ePC-SAFT as an Engine parameter document.

Python only declares the published parameters; every property is solved natively by the Engine.
"""
import math

from epcsaft import Mixture, Parameters, Pure

ARTICLE, CORRECTION, PERMITTIVITY_SOURCE = "cameretti-2005", "cameretti-2005-correction", "floriano-nascimento-2004"
SOURCES = [
    {"source_id": ARTICLE, "citation": "Cameretti, Sadowski, Mollerup, Ind. Eng. Chem. Res. 2005, 44, 3355-3362", "doi": "10.1021/ie0488142"},
    {"source_id": CORRECTION, "citation": "Cameretti, Sadowski, Mollerup, Additions and Corrections, Ind. Eng. Chem. Res. 2005, 44, 8944", "doi": "10.1021/ie051055i"},
    {"source_id": PERMITTIVITY_SOURCE, "citation": "Floriano, Nascimento, Braz. J. Phys. 2004, 34, 38", "doi": "10.1590/S0103-97332004000100006"},
    {"source_id": "held-2014", "citation": "Held, Reschke, Mohammad, Luza, Sadowski, Chem. Eng. Res. Des. 2014, 92, 2884", "doi": "10.1016/j.cherd.2014.05.017"},
    {"source_id": "nist-srd-69-water", "citation": "NIST Chemistry WebBook, SRD 69, Water, CAS 7732-18-5"},
    {"source_id": "iupac-standard-atomic-weights", "citation": "IUPAC standard atomic weights"},
]
# Article Table 1 (water, 2B); molar masses [g/mol]; ion tables give charge, sigma_j [A], u_j^0 [K], molar mass.
WATER = dict(m=1.09528, sigma=2.88980, u=365.956, eps_ab=2515.6706, kappa_ab=0.0348679836, M=18.01528)
# Corrected Table 2 (Ind. Eng. Chem. Res. 2005, 44, 8944), which replaces the original page-3359 table: a_j [A], u_j^0 [K].
IONS = {
    "Li+": (1, 2.795, 687.096, 6.94), "Na+": (1, 2.911, 157.728, 22.98977),
    "K+": (1, 3.674, 142.057, 39.0983), "Cl-": (-1, 2.986, 348.322, 35.453),
    "Br-": (-1, 3.304, 698.310, 79.904), "I-": (-1, 3.822, 379.295, 126.90447),
    "SO4^2-": (-2, 4.361, 667.230, 96.056),
}
# Original Table 2 of the article (p. 3359), which the correction declares incorrect; it is the table that reproduces the
# article's Figure 2 vapor-pressure curves.
IONS_ORIGINAL = {
    "Li+": (1, 1.8059, 1110.9261, 6.94), "Na+": (1, 1.6262, 119.8060, 22.98977),
    "K+": (1, 2.7602, 8.8773, 39.0983), "Cl-": (-1, 3.5991, 359.6604, 35.453),
    "Br-": (-1, 3.8225, 524.0636, 79.904), "I-": (-1, 4.1766, 413.0494, 126.90447),
    "SO4^2-": (-2, 4.3294, 159.9912, 96.056),
}
SALTS = {"NaCl": ("Na+", "Cl-", 1, 1), "NaBr": ("Na+", "Br-", 1, 1), "NaI": ("Na+", "I-", 1, 1),
         "LiCl": ("Li+", "Cl-", 1, 1), "LiBr": ("Li+", "Br-", 1, 1), "LiI": ("Li+", "I-", 1, 1),
         "KCl": ("K+", "Cl-", 1, 1), "KBr": ("K+", "Br-", 1, 1), "KI": ("K+", "I-", 1, 1),
         "Li2SO4": ("Li+", "SO4^2-", 2, 1), "Na2SO4": ("Na+", "SO4^2-", 2, 1), "K2SO4": ("K+", "SO4^2-", 2, 1)}

# Floriano and Nascimento (2004) Eq. (3), Table I: eps(T, P) = eps(T, 10 MPa) + a0 ln((a1 + P)/(a1 + 10)), P in MPa,
# fitted at four temperatures only. Cameretti evaluate it at 100 kPa; the temperature interpolation is not printed.
FN_TABLE = {273.0: (88.28, 11.9240, 232.5026), 298.0: (78.85, 14.1113, 341.5902),
            323.0: (70.27, 14.3899, 386.4825), 348.0: (62.59, 13.9450, 393.2031)}
FN_POINTS = {T: e + a0 * math.log((a1 + 0.1) / (a1 + 10.0)) for T, (e, a0, a1) in FN_TABLE.items()}
T_REF = 298.15


def _cubic_through_points():
    """Coefficients c0..c3 of the cubic in (T - 298.15) through the four 0.1 MPa values (Lagrange form expanded)."""
    xs = [T - T_REF for T in FN_POINTS]
    ys = list(FN_POINTS.values())
    coefficients = [0.0] * 4
    for i, (xi, yi) in enumerate(zip(xs, ys)):
        poly = [1.0]
        scale = yi
        for j, xj in enumerate(xs):
            if j != i:
                poly = [a - xj * b for a, b in zip([0.0] + poly, poly + [0.0])]
                scale /= xi - xj
        coefficients = [c + scale * p for c, p in zip(coefficients, poly)]
    return coefficients


PERMITTIVITY = _cubic_through_points()


def water_permittivity(T):
    return sum(c * (T - T_REF) ** n for n, c in enumerate(PERMITTIVITY))


PACKING_SOURCE = {"sigma": (ARTICLE, "Section 3, charged hard spheres of diameter sigma_j"),
                  "0.88": ("held-2014", "Eqs. (21)-(22), d_j = 0.88 sigma_j")}


def _q(value, unit):
    return {"magnitude": value, "unit": unit}


def _prov(source, locator):
    return {"source_id": source, "locator": locator, "domain_id": "cameretti-2005-278-to-393-k"}


def _coef(cid, family, value, unit, source, locator):
    return {"identity": f"component/{cid}/{family}", "family": family, "value": _q(value, unit), "provenance": _prov(source, locator)}


def document(ions, packing="sigma", table=None):
    """Parameter document for water plus `ions` from `table` (IONS, the corrected Table 2, or IONS_ORIGINAL).

    packing selects the ionic hard-sphere diameter, which Cameretti do not state as a formula:
      "sigma"  d_j = sigma_j (charged hard spheres of diameter sigma_j, Section 3 text),
      "0.88"   d_j = 0.88 sigma_j (Held et al. 2014, Eqs. 21-22).
    The Debye-Hueckel closest approach is a_j = sigma_j (Section 3) in every case.
    """
    table = table or IONS
    ion_source, ion_table = (ARTICLE, "original Table 2, p. 3359") if table is IONS_ORIGINAL else (CORRECTION, "corrected Table 2")
    components = [{
        "component_id": "water", "fixed": {"molar_mass": {"value": _q(WATER["M"], "gram / mole"), "provenance": _prov("nist-srd-69-water", "molecular weight")},
                                           "charge_number": {"value": _q(0, "elementary-charge"), "provenance": _prov(ARTICLE, "Table 1, water")}},
        "coefficients": [_coef("water", "segment_count", WATER["m"], "dimensionless", ARTICLE, "Table 1"),
                         _coef("water", "dispersion_energy_over_k", WATER["u"], "kelvin", ARTICLE, "Table 1")]}]
    for ion in ions:
        z, sigma, u, M = table[ion]
        d = {"sigma": sigma, "0.88": 0.88 * sigma}[packing]
        components.append({
            "component_id": ion,
            "fixed": {"molar_mass": {"value": _q(M, "gram / mole"), "provenance": _prov("iupac-standard-atomic-weights", "sum of element standard atomic weights")},
                      "charge_number": {"value": _q(z, "elementary-charge"), "provenance": _prov(ion_source, f"{ion_table}, ion label")}},
            "coefficients": [_coef(ion, "segment_count", 1, "dimensionless", ARTICLE, "Section 3, ions as single hard spheres"),
                             _coef(ion, "segment_diameter", sigma, "angstrom", ion_source, f"{ion_table}, sigma_j"),
                             _coef(ion, "dispersion_energy_over_k", u, "kelvin", ion_source, f"{ion_table}, u_j^0"),
                             _coef(ion, "debye_huckel_diameter", sigma, "angstrom", ARTICLE, "Section 3, a_j = sigma_j"),
                             _coef(ion, "packing_diameter", d, "angstrom", *PACKING_SOURCE[packing])]})
    pairs = []
    for i, a in enumerate(ions):
        for b in ions[i + 1:]:
            if IONS[a][0] * IONS[b][0] < 0:  # u_ij = 0 for i != water: cation-anion k_ij = 1; same sign excluded by the family choice
                pairs.append({"component_id_a": a, "component_id_b": b, "coefficients": [
                    {"identity": f"pair/{a}/{b}/k_ij", "family": "k_ij", "value": _q(1.0, "dimensionless"),
                     "provenance": _prov(ARTICLE, "Section 3, u_ij^0 = 0 for i != 1")}]})
        pairs.append({"component_id_a": "water", "component_id_b": a, "coefficients": [
            {"identity": f"pair/water/{a}/k_ij", "family": "k_ij", "value": _q(0.0, "dimensionless"),
             "provenance": _prov(ARTICLE, "Section 3, all k_ij = 0 for water-ion")}]})
    families = [("base", "pc-saft"), ("association", "general-site"), ("permittivity", "solvent-only")]
    model_families = [{"family_id": f"model/{k}", "kind": k, "choice": c, "provenance": _prov(ARTICLE, "Eq. (15)")} for k, c in families]
    model_families.append({"family_id": "model/electrolyte", "kind": "electrolyte", "choice": "fully-dissociated-debye-huckel",
                           "ionic_dispersion": "exclude-same-sign-ion-pairs", "provenance": _prov(CORRECTION, "corrected Eq. (13); Eqs. (14)-(15) and Section 3 of the article")})
    cid = "component/water/relative_permittivity/constant-plus-polynomial-in-temperature"
    correlations = [
        {"correlation_id": cid, "component_id": "water", "family": "relative_permittivity",
         "form": "constant-plus-polynomial-in-temperature", "independent_variables": ["temperature"],
         "reference_temperature": _q(T_REF, "kelvin"), "constant": {"identity": f"{cid}/constant", "value": _q(PERMITTIVITY[0], "dimensionless")},
         "terms": [{"power": n, "coefficient": {"identity": f"{cid}/term-{n}/coefficient", "value": _q(c, "1 / kelvin" if n == 1 else f"1 / kelvin ** {n}")}}
                   for n, c in enumerate(PERMITTIVITY) if n],
         "provenance": _prov(PERMITTIVITY_SOURCE, "Eq. (3) and Table I at 0.1 MPa, cubic through 273/298/323/348 K")},
        {"correlation_id": "component/water/segment_diameter/constant", "component_id": "water", "family": "segment_diameter",
         "form": "constant", "independent_variables": [],
         "constant": {"identity": "component/water/segment_diameter", "value": _q(WATER["sigma"], "angstrom")}, "provenance": _prov(ARTICLE, "Table 1")}]
    edge = {"endpoint_a": {"component_id": "water", "site_id": "a"}, "endpoint_b": {"component_id": "water", "site_id": "b"},
            "energy_over_k": {"identity": "association/water/energy_over_k", "value": _q(WATER["eps_ab"], "kelvin")},
            "volume": {"identity": "association/water/volume", "value": _q(WATER["kappa_ab"], "dimensionless")},
            "source": {"kind": "explicit", "provenance": [_prov(ARTICLE, "Table 1")]}}
    sites = [{"component_id": "water", "site_id": s, "site_role": r, "multiplicity": 1, "provenance": _prov(ARTICLE, "Table 1, two sites")}
             for s, r in (("a", "donor"), ("b", "acceptor"))]
    return {"schema": "epcsaft.parameters", "schema_version": 1, "document_id": "cameretti-2005", "document_version": 1,
            "purpose": "external-document", "sources": SOURCES,
            "domains": [{"domain_id": "cameretti-2005-278-to-393-k", "kind": "reported-conditions"}],
            "model_families": model_families, "components": components, "pairs": pairs, "correlations": correlations,
            "topology": {"presets": [], "sites": sites, "edges": [edge]}}


def mixture(ions, packing="sigma", table=None):
    doc = document(ions, packing, table)
    params = Parameters.from_mapping(doc, components=["water", *ions])
    return Mixture(params)


def water():
    return Pure(Parameters.from_mapping(document([]), components=["water"]))


def declarations():
    return (
        water(),
        [
            mixture(list(i), p, globals()[t])
            for i in [s[:2] for s in SALTS.values()] + [('Na+', 'K+', 'Cl-', 'Br-')]
            for p in ('sigma', '0.88')
            for t in ('IONS', 'IONS_ORIGINAL')
        ],
    )
