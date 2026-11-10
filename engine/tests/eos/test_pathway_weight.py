"""Hong et al. 2027 pathway-weighted reactive cross-association declared through ``Parameters.from_mapping``.

Tables 1-3 values; each reactive CO2 site bonds MEA site b only, a test declaration and not a chemistry choice. The
native falsifiers (properties_checks PathwayWeight) own the numerical verification; this test covers the document
mapping: block membership, the reactive-site correction, the weight and the complement k_ij law.
"""

from __future__ import annotations

import math
from typing import Any

import pytest

from epcsaft.eos import Mixture
from epcsaft.parameters import ParameterError, Parameters

T_REF = 298.15
PURE = {"mea": (3.1200, 3.0229, 287.94, 0.06108352), "water": (2.1945, 2.2290, 141.66, 0.01801528),
        "co2": (2.0729, 2.7852, 169.21, 0.0440095)}
PATHWAY = {"sym": (3.675e-3, -2.3981, 1.555e-3, 4.1160, 3546.55, 2), "single": (5.939e-3, -1.9778, 6.5e-5, 12.5398, 1197.42, 1)}
EPS_MW = 1.0374 * math.sqrt(2355.03 * 1804.17)
KAPPA_MW = ((0.0313 ** (1 / 3) + 0.2039 ** (1 / 3)) / 2) ** 3
L_STATE = (333.15, 42424.79, [0.10404183733215595, 0.8231288765353348, 0.07282928613250916])


def value(identity: str, magnitude: float, unit: str = "dimensionless") -> dict[str, Any]:
    return {"identity": identity, "value": {"magnitude": magnitude, "unit": unit}}


def law(identity: str, slope: float, intercept: float, units: tuple[str, str] = ("dimensionless", "1 / kelvin")) -> dict[str, Any]:
    terms = [{"power": 0, "coefficient": value(f"{identity}/a0", intercept + slope * T_REF, units[0])},
             {"power": 1, "coefficient": value(f"{identity}/a1", slope, units[1])}]
    return {"identity": identity, "temperature_law": {"form": "polynomial-plus-reciprocal-temperature",
                                                      "reference_temperature": {"magnitude": T_REF, "unit": "kelvin"}, "terms": terms}}


def document(pathways: tuple[str, ...], intercept: float = 1.7) -> dict[str, Any]:
    """One block per pathway; two pathways weigh the corrected symmetric block by w and the single-site block by 1 - w."""
    blend = len(pathways) == 2
    sites, edges = [], []
    cross = {"identity": "association/mea-water/energy_over_k", "temperature_law": {
        "form": "polynomial-plus-reciprocal-temperature", "reference_temperature": {"magnitude": T_REF, "unit": "kelvin"},
        "terms": [{"power": 0, "coefficient": value("association/mea-water/energy_over_k/a0", EPS_MW, "kelvin")}]}}
    for block, pathway in enumerate(pathways):
        _, _, kappa, c, d, nu = PATHWAY[pathway]
        tag = f"/{block}" if blend else ""
        for component, site, multiplicity in (("mea", "a", 1), ("mea", "b", 1), ("water", "a", 2), ("water", "b", 2), ("co2", "r", nu)):
            sites.append({"component_id": component, "site_id": site + tag, "multiplicity": multiplicity}
                         | ({"pathway": ("weighted", "complement")[block]} if blend else {})
                         | ({"symmetric_reactive_site_correction": True} if pathway == "sym" and site == "r" and blend else {}))

        def edge(a: tuple[str, str], b: tuple[str, str], energy: dict[str, Any], volume: float, tag: str = tag) -> dict[str, Any]:
            name = f"association/{a[0]}/{a[1]}{tag}/{b[0]}/{b[1]}{tag}"
            return {"endpoint_a": {"component_id": a[0], "site_id": a[1] + tag}, "endpoint_b": {"component_id": b[0], "site_id": b[1] + tag},
                    "source": {"kind": "explicit"}, "energy_over_k": energy if "temperature_law" in energy else energy | {"identity": f"{name}/energy"},
                    "volume": value(f"{name}/volume", volume)}

        edges += [edge(("mea", "a"), ("mea", "b"), value("", 2355.03, "kelvin"), 0.0313),
                  edge(("water", "a"), ("water", "b"), value("", 1804.17, "kelvin"), 0.2039),
                  edge(("mea", "a"), ("water", "b"), dict(cross), KAPPA_MW), edge(("mea", "b"), ("water", "a"), dict(cross), KAPPA_MW),
                  edge(("co2", "r"), ("mea", "b"), law(f"association/r{tag}/energy", c, d, ("kelvin", "dimensionless")), kappa)]
    first = PATHWAY[pathways[0]]
    mea_co2 = [law("pair/mea/co2/k_ij", first[0], first[1]) | {"family": "k_ij"}]
    if blend:
        mea_co2.append(law("pair/mea/co2/k_ij_complement", *PATHWAY["single"][:2]) | {"family": "k_ij_complement"})
    result: dict[str, Any] = {
        "schema": "epcsaft.parameters", "schema_version": 1,
        "model_families": [{"kind": "base", "choice": "pc-saft"}, {"kind": "association", "choice": "general-site"},
                           {"kind": "electrolyte", "choice": "none"}, {"kind": "permittivity", "choice": "none"}],
        "components": [{"component_id": cid, "fixed": {"charge_number": {"value": {"magnitude": 0, "unit": "elementary-charge"}},
                                                       "molar_mass": {"value": {"magnitude": mass, "unit": "kilogram / mole"}}},
                        "coefficients": [value(f"{cid}/m", m) | {"family": "segment_count"},
                                         value(f"{cid}/sigma", sigma, "angstrom") | {"family": "segment_diameter"},
                                         value(f"{cid}/eps", eps, "kelvin") | {"family": "dispersion_energy_over_k"}]}
                       for cid, (m, sigma, eps, mass) in PURE.items()],
        "pairs": [{"component_id_a": "mea", "component_id_b": "water", "coefficients": [law("pair/mea/water/k_ij", 9.0e-5, -0.0603) | {"family": "k_ij"}]},
                  {"component_id_a": "water", "component_id_b": "co2", "coefficients": [law("pair/water/co2/k_ij", 2.46e-4, -0.0633) | {"family": "k_ij"}]},
                  {"component_id_a": "mea", "component_id_b": "co2", "coefficients": mea_co2}],
        "topology": {"sites": sites, "edges": edges,
                     "no_cross_association": [{"component_id_a": "water", "component_id_b": "co2", "provenance": {"locator": "Hong 2027 2.2"}}]}}
    if blend:
        result["pathway_weight"] = {"loading": {"numerator": "co2", "denominator": "mea"}, "coefficients": [
            value("w/b0", intercept) | {"family": "intercept"}, value("w/b_alpha", -10.0) | {"family": "loading_slope"},
            value("w/b_T", -1.5) | {"family": "temperature_slope"}]}
    return result


def test_declared_blend_reports_its_weight_and_reaches_the_single_site_limit() -> None:
    temperature, density, x = L_STATE
    blend = Mixture(Parameters.from_mapping(document(("sym", "single")))).state(temperature, rho=density, x=x)
    z = 1.7 - 1.5 * (temperature - 353.15) / 40.0 - 10.0 * (x[2] / x[0] - 0.5)
    assert blend.pathway_weight == pytest.approx(1.0 / (1.0 + math.exp(-z)), rel=1e-15, abs=0.0)
    assert len(blend.site_fractions) == 10
    limit = Mixture(Parameters.from_mapping(document(("sym", "single"), intercept=-800.0))).state(temperature, rho=density, x=x)
    single = Mixture(Parameters.from_mapping(document(("single",)))).state(temperature, rho=density, x=x)
    assert limit.pathway_weight == 0.0
    assert abs(limit.residual_helmholtz - single.residual_helmholtz) <= 1e-10 * abs(single.residual_helmholtz)


def test_declared_blend_matches_the_native_declaration_at_l() -> None:
    """The same blend declared natively (properties_checks PathwayWeight, exact rows of F3/L/caught) at state L."""
    temperature, density, x = L_STATE
    state = Mixture(Parameters.from_mapping(document(("sym", "single")))).state(temperature, rho=density, x=x)
    native_mu = (-24.70357855175673, -8.357368873159238, 9.244225781433842)
    assert state.residual_chemical_potential_over_rt == pytest.approx(native_mu, rel=1e-10, abs=0.0)
    assert state.residual_internal_energy == pytest.approx(-47321.65391952424, rel=1e-10, abs=0.0)


def test_an_edge_across_pathway_blocks_is_refused() -> None:
    crossing = document(("sym", "single"))
    crossing["topology"]["edges"][-1]["endpoint_b"]["site_id"] = "b/0"  # single-site r bonds the weighted MEA b
    with pytest.raises(ParameterError, match="crosses pathway blocks"):
        Parameters.from_mapping(crossing)
