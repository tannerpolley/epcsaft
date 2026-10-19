"""#641: a declared charged-molecule ion keeps its self-dispersion and the neutral diameter sigma (1 - 0.12 exp(-3u/kT))
(Held 2014 Eqs. 21-22). Exclude-same-sign-ion-pairs with the declaration equals all-ion-pairs with every other ion's
self-dispersion suppressed by s and the ion's packing diameter set to that d(T) (the Amine #126 construction)."""
import math
from typing import Any

import pytest

from epcsaft import Mixture, ParameterError, Parameters

W, MDEA, MDEAH, HCO3, H = "water", "methyldiethanolamine", "protonated-methyldiethanolamine", "bicarbonate-anion", "hydronium-cation"
# component: molar mass / (kg/mol), charge, m, sigma / angstrom, u/k / K (Uyan 2015 Tables 3-4; Wangler 2018 Tables 4-5)
SPECIES = {W: (0.01801528, 0, 1.2046, 2.7927, 353.94), MDEA: (0.1191622, 0, 3.675, 3.563, 228.711),
           MDEAH: (0.12017014, 1, 3.675, 3.563, 228.71), HCO3: (0.06101684, -1, 1.0, 2.9296, 70.0),
           H: (0.00100794, 1, 1.0, 3.4654, 500.0)}
# k_ij polynomials in T - 298.15 K (Wangler 2018 Table 6); every other pair is zero
KIJ = {(W, MDEA): (-0.184, -1.21e-3, 2.70e-6, 1.53e-7), (W, MDEAH): (-0.184, -1.21e-3, 2.70e-6, 1.53e-7),
       (W, HCO3): (-0.0765, -6.25e-4, -6.25e-5), (W, H): (0.25,)}
SITES = {W: (1, 2425.67, 0.04509), MDEA: (2, 2046.624, 0.123858), MDEAH: (2, 0.0, 0.123858)}  # Uyan 2015 Tables 3-4


def value(magnitude: float, unit: str = "dimensionless") -> dict[str, Any]:
    return {"value": {"magnitude": magnitude, "unit": unit}}


def document(T: float, suppress: float | None = None) -> dict[str, Any]:
    """MDEAH+ declared a charged molecule, or with `suppress` the all-ion-pairs construction at factor s."""
    scaled, factor = ((), 1.0) if suppress is None else ((HCO3, H), suppress)
    components = []
    for cid, (mass, charge, m, sigma, u) in SPECIES.items():
        coefficients = [("segment_count", m, "dimensionless"),
                        ("dispersion_energy_over_k", u * (factor if cid in scaled else 1.0), "kelvin")]
        if cid != W:  # water's sigma(T) is its correlation
            coefficients.append(("segment_diameter", sigma, "angstrom"))
        if cid in (W, MDEA):
            coefficients.append(("relative_permittivity", 78.41 if cid == W else 22.358, "dimensionless"))
        if charge:
            coefficients.append(("debye_huckel_diameter", sigma, "angstrom"))
        if charge and (cid != MDEAH or suppress is not None):
            packing = sigma * (1 - 0.12 * math.exp(-3 * u / T)) if cid == MDEAH else 0.88 * sigma
            coefficients.append(("packing_diameter", packing, "angstrom"))
        fixed = {"molar_mass": value(mass, "kilogram / mole"), "charge_number": value(charge, "elementary-charge")}
        if cid == MDEAH and suppress is None:
            fixed["charged_molecule_ion"] = value(1)
        components.append({"component_id": cid, "fixed": {k: {"identity": f"component/{cid}/{k}", **v} for k, v in fixed.items()},
                           "coefficients": [{"family": f, "identity": f"component/{cid}/{f}", **value(v, unit)}
                                            for f, v, unit in coefficients]})
    pairs = []
    for k, a in enumerate(SPECIES):
        for b in list(SPECIES)[k + 1:]:
            same_sign = SPECIES[a][1] * SPECIES[b][1] > 0
            if same_sign and suppress is None:
                continue
            k_ij = 1.0 if same_sign else sum(c * (T - 298.15) ** p for p, c in enumerate(KIJ.get((a, b), (0.0,))))
            for _ in {a, b} & set(scaled):  # keeps sqrt(u_a u_b) (1 - k_ij) exact
                k_ij = 1 - (1 - k_ij) / math.sqrt(factor)
            pairs.append({"component_id_a": a, "component_id_b": b,
                          "coefficients": [{"family": "k_ij", "identity": f"pair/{a}/{b}/k_ij", **value(k_ij)}]})
    sites = [{"component_id": c, "site_id": s, "site_role": role, "multiplicity": n}
             for c, (n, _, _) in SITES.items() for s, role in (("a", "donor"), ("b", "acceptor"))]
    edges = [{"endpoint_a": {"component_id": c, "site_id": "a"}, "endpoint_b": {"component_id": c, "site_id": "b"},
              "source": {"kind": "explicit"}, "energy_over_k": {"identity": f"association/{c}/energy_over_k", **value(e, "kelvin")},
              "volume": {"identity": f"association/{c}/volume", **value(v)}} for c, (_, e, v) in SITES.items()]
    edges += [{"endpoint_a": {"component_id": d, "site_id": "a"}, "endpoint_b": {"component_id": r, "site_id": "b"},
               "source": {"kind": "combining-rule", "rule_id": "wolbach-sandler", "inputs": {
                   "energy_over_k": [f"association/{i}/energy_over_k", f"association/{j}/energy_over_k"],
                   "volume": [f"association/{i}/volume", f"association/{j}/volume"]}}}
              for i, j in ((W, MDEA), (W, MDEAH), (MDEA, MDEAH)) for d, r in ((i, j), (j, i))]
    water_sigma = "component/water/segment_diameter"
    return {"schema": "epcsaft.parameters", "schema_version": 1, "document_id": "charged-molecule-ion", "document_version": 1,
            "components": components, "pairs": pairs, "topology": {"sites": sites, "edges": edges},
            "correlations": [{"component_id": W, "correlation_id": water_sigma, "family": "segment_diameter",
                              "form": "constant-plus-sum-of-exponentials", "constant": {"identity": f"{water_sigma}/constant", **value(2.7927, "angstrom")},
                              "terms": [{"amplitude": {"identity": f"{water_sigma}/{k}/amplitude", **value(a, "angstrom")},
                                         "exponent_coefficient": {"identity": f"{water_sigma}/{k}/exponent", **value(b, "1 / kelvin")}}
                                        for k, (a, b) in enumerate(((10.11, -0.01775), (-1.417, -0.01146)))]}],
            "model_families": [{"kind": "base", "choice": "pc-saft"}, {"kind": "association", "choice": "general-site"},
                               {"kind": "electrolyte", "choice": "fully-dissociated-debye-huckel",
                                "ionic_dispersion": "all-ion-pairs" if suppress else "exclude-same-sign-ion-pairs"},
                               {"kind": "permittivity", "choice": "solvent-only"}]}


def composition(loading: float, m_mdea: float = 2.0) -> list[float]:
    n = [1 / 0.01801528, m_mdea * (1 - loading) + 1e-6, m_mdea * loading, m_mdea * loading + 1e-6, 1e-6]
    return [v / sum(n) for v in n]


@pytest.mark.parametrize("T", [313.15, 373.15, 413.15])
def test_declaration_equals_the_suppressed_workaround(T: float) -> None:
    declared, workaround = (Mixture(Parameters.from_mapping(document(T, s))) for s in (None, 1e-12))
    for loading in (0.1, 0.5, 1.0):
        x = composition(loading)
        a = declared.state(T, P=1e5, x=x, phase="liquid")
        b = workaround.state(T, rho=a.molar_density, x=x)
        assert a.log_fugacity_coefficient is not None and b.log_fugacity_coefficient is not None
        pairs = [(a.residual_helmholtz, b.residual_helmholtz), (a.compressibility_factor, b.compressibility_factor),
                 *zip(a.log_fugacity_coefficient, b.log_fugacity_coefficient)]
        for left, right in pairs:
            assert left is not None and right is not None
            assert abs(left - right) <= 1e-10 * max(1.0, abs(right)), (loading, left, right)


def test_declaration_is_refused_where_it_does_not_apply() -> None:
    def refused(match: str, edit: Any) -> None:
        mapping = document(313.15)
        edit({c["component_id"]: c for c in mapping["components"]})
        with pytest.raises(ParameterError, match=match):
            Parameters.from_mapping(mapping)
    refused("unconsumed component/water/charged_molecule_ion",
            lambda c: c[W]["fixed"].update(charged_molecule_ion={"identity": "w", **value(1)}))
    refused("unconsumed component/protonated-methyldiethanolamine/packing_diameter",
            lambda c: c[MDEAH]["coefficients"].append({"family": "packing_diameter", "identity": "p", **value(3.1, "angstrom")}))
    refused("charged_molecule_ion is not 1", lambda c: c[MDEAH]["fixed"]["charged_molecule_ion"]["value"].update(magnitude=2))
