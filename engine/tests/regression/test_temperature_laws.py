"""#567: temperature-law declarations travel to fit workers.

The pinned four-site case of the native temperature-laws check, declared as a parameter document: k01, l01,
E00, v00, E11 and the shared Wolbach--Sandler volume correction c01 are q = sum_p a_p (T - 300 K)^p
+ s (1/T - 1/300 K). One process and two workers return the same residuals and 3 x 30 Jacobian at u = 0;
the 320 K state group is the one the round-robin assignment sends to the worker. The native check verifies
these columns against rebuilt-model differences.
"""

import math

from epcsaft import ActiveParameter, Mixture, ParameterFamily, Parameters, regression

T_R = {"magnitude": 300.0, "unit": "kelvin"}
DIMENSIONLESS = ("dimensionless", "1 / kelvin", "1 / kelvin ** 2", "1 / kelvin ** 3", "kelvin")
ENERGY = ("kelvin", "dimensionless", "1 / kelvin", "1 / kelvin ** 2", "kelvin ** 2")
LAWS = {  # [a0, a1, a2, a3, s], coordinate scales
    "k01": ([-0.052, 1e-4, 2e-6, 1e-8, 10.0], DIMENSIONLESS), "l01": ([0.01, 2e-5, 1e-7, -1e-9, 1.0], DIMENSIONLESS),
    "E00": ([2383.4744, 0.8, 0.002, 1e-5, 3000.0], ENERGY), "v00": ([0.118488, 2e-4, -5e-7, 2e-9, 0.05], DIMENSIONLESS),
    "E11": ([1804.17, -0.4, 0.001, -3e-6, 1500.0], ENERGY), "c01": ([0.15, 1e-3, -2e-6, 1e-8, 3.0], DIMENSIONLESS),
}
SCALES = {DIMENSIONLESS: (0.01, 1e-3, 1e-5, 1e-7, 10.0), ENERGY: (100.0, 1.0, 0.01, 1e-4, 3000.0)}
X = {"a": 0.4, "b": 0.6}


def _law(name: str) -> dict[str, object]:
    values, units = LAWS[name]
    terms = [{"power": p, "coefficient": {"identity": f"manufactured/{name}/a{p}", "value": {"magnitude": values[p], "unit": units[p]}}}
             for p in range(4)]
    return {"identity": f"manufactured/{name}", "temperature_law": {
        "form": "polynomial-plus-reciprocal-temperature", "reference_temperature": T_R, "terms": terms,
        "reciprocal_temperature_slope": {"identity": f"manufactured/{name}/s", "value": {"magnitude": values[4], "unit": units[4]}}}}


def _document() -> dict[str, object]:
    def component(cid: str, m: float, sigma: float, eps: float, mass: float) -> dict[str, object]:
        return {"component_id": cid, "fixed": {"molar_mass": {"value": {"magnitude": mass, "unit": "kilogram / mole"}},
                                                "charge_number": {"value": {"magnitude": 0, "unit": "elementary-charge"}}},
                "coefficients": [{"identity": f"{cid}/{f}", "family": f, "value": {"magnitude": v, "unit": u}} for f, v, u in
                                 (("segment_count", m, "dimensionless"), ("segment_diameter", sigma, "angstrom"),
                                  ("dispersion_energy_over_k", eps, "kelvin"))]}

    def end(cid: str, site: str) -> dict[str, str]:
        return {"component_id": cid, "site_id": site}

    rule_inputs = {"energy_over_k": ["manufactured/E00", "manufactured/E11"], "volume": ["manufactured/v00", "v11"]}
    rule = {"kind": "combining-rule", "rule_id": "wolbach-sandler", "inputs": rule_inputs}
    return {
        "schema": "epcsaft.parameters", "schema_version": 1,
        "components": [component("a", 4.5354, 2.6019, 204.0438, 0.06108), component("b", 2.1945, 2.229, 141.66, 0.01801528)],
        "pairs": [{"component_id_a": "a", "component_id_b": "b",
                   "coefficients": [{"family": "k_ij", **_law("k01")}, {"family": "l_ij", **_law("l01")}]}],
        "model_families": [{"kind": "base", "choice": "pc-saft"}, {"kind": "association", "choice": "general-site"},
                           {"kind": "electrolyte", "choice": "none"}, {"kind": "permittivity", "choice": "none"}],
        "topology": {
            "sites": [{"component_id": c, "site_id": s, "multiplicity": n} for c, s, n in (("a", "A", 2), ("a", "B", 1), ("b", "A", 2), ("b", "B", 2))],
            "edges": [
                {"endpoint_a": end("a", "A"), "endpoint_b": end("a", "B"), "source": {"kind": "explicit"},
                 "energy_over_k": _law("E00"), "volume": _law("v00")},
                {"endpoint_a": end("b", "A"), "endpoint_b": end("b", "B"), "source": {"kind": "explicit"},
                 "energy_over_k": _law("E11"), "volume": {"identity": "v11", "value": {"magnitude": 0.2039, "unit": "dimensionless"}}},
                {"endpoint_a": end("a", "A"), "endpoint_b": end("b", "B"), "source": rule, "volume_correction": _law("c01")},
                {"endpoint_a": end("a", "B"), "endpoint_b": end("b", "A"), "source": rule, "volume_correction": _law("c01")},
            ]},
    }


def test_temperature_law_coefficients_reach_fit_workers() -> None:
    parameters = Parameters.from_mapping(_document())
    identities = parameters._native.temperature_law_coefficient_identities
    assert len(identities) == 30  # c01 is declared on both cross edges and owns one set of five slots
    coordinates = []
    for name, (values, units) in LAWS.items():
        for p, (value, scale) in enumerate(zip(values, SCALES[units], strict=True)):
            identity = f"manufactured/{name}/" + ("s" if p == 4 else f"a{p}")
            coordinates.append(regression.coordinate(parameters, [("temperature_law_coefficient", identity, 1.0)], origin=value,
                                                     scale=scale, bounds=(value - 0.1 * scale, value + 0.1 * scale)))
    model = Mixture(parameters)
    observations = [regression.density(parameters, T, model.state(T, rho=100.0, x=X).pressure, X, phase="vapor", observed=101.0,
                                       form="difference", scale=100.0, mass_basis=False) for T in (300.0, 320.0, 340.0)]
    serial, parallel = (regression.fit(parameters, coordinates, observations, controls=regression.FitControls(
        maximum_iterations=1, maximum_elapsed_time_seconds=30.0, gradient_tolerance=1e6, workers=workers)) for workers in (1, 2))
    for result in (serial, parallel):
        assert list(result.optimizer) == [0.0] * 30, result.message
        assert all(abs(r + 0.01) <= 1e-10 for r in result.residuals)
    for one, two in ((serial.residuals, parallel.residuals), (serial.optimizer_jacobian, parallel.optimizer_jacobian)):
        assert len(one) == len(two) and all(abs(a - b) <= 1e-10 * max(1.0, abs(a), abs(b)) for a, b in zip(one, two, strict=True))
    # Every law reaches the 320 and 340 K rows (T != T_r), so detached fields cannot pass by equality of zeros.
    jacobian = serial.optimizer_jacobian
    for law in range(6):
        assert any(jacobian[row * 30 + 5 * law + p] != 0.0 for row in (1, 2) for p in range(5))
    # The fit's columns are the observation owner's exact gradient in u = (coefficient - origin) / scale.
    active = [ActiveParameter(ParameterFamily.TemperatureLawCoefficient, [slot]) for slot in range(30)]
    gradient_model = Mixture(parameters, active_parameters=active)
    scales = [coordinate.scale for coordinate in coordinates]
    for row, item in enumerate(observations):
        gradient = regression.differentiate(gradient_model, item).value
        for column in range(30):
            exact = scales[column] * gradient[column] / 100.0
            assert math.isclose(jacobian[row * 30 + column], exact, rel_tol=1e-10, abs_tol=1e-12)
