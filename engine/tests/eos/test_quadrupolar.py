"""#574: optional Gross (2005) quadrupole-quadrupole term (ADR 0014, full EOS 4.2).

a_qq is isolated as a^r(polar family) - a^r(same parameters without it) at fixed
(T, rho, x). The checks falsify the equations against Clapeyron 0.6.29 QPCPSAFT rows
(moments rescaled to the Engine's debye constant), FeOs 0.10.1 where every m_i <= 2,
the second-virial limit written
independently here from Gross (2005) Eq. 9 and Table 2 without a chain-length cap,
CO2 and benzene saturation pressure, the CO2/N2 bubble point (code identity of the
mixing rule, not validation), the derivative chains, and the absent-quadrupole limit.
"""

import copy
import csv
import json
import math
from pathlib import Path
from typing import Any

from epcsaft import DensityActionDirection, Mixture, Parameters, PropertyActionStatus, PropertyObservable, Pure
from epcsaft.equilibrium import Amounts, Free, Phase, Pinned, Problem, solve_equilibrium
from tests.native_inputs import mixed_difference
from tests.packets import parameters_json

RESULTS = Path(__file__).resolve().parents[3] / "analyses/2026-polar-quadrupole/results"
AVOGADRO, BOLTZMANN, EPS0, DEBYE = 6.02214076e23, 1.380649e-23, 8.8541878128e-12, 1.0e-21 / 299792458.0


def _document(names: tuple[str, ...], polar: bool) -> dict[str, Any]:
    """Gross (2005) Table A1 records plus Gross-Sadowski 2001 methane (Q = 0), all k_ij = 0."""
    gross = json.loads(parameters_json("gross-2005", 1).read_text())
    gs = json.loads(parameters_json("gross-sadowski-2001", 3).read_text())
    methane = copy.deepcopy(next(c for c in gs["components"] if c["component_id"] == "methane"))
    methane["coefficients"] += [
        {"identity": f"methane/{family}", "family": family, "value": {"magnitude": 0.0, "unit": unit}}
        for family, unit in (("quadrupole_moment", "debye * angstrom"), ("quadrupole_count", "dimensionless"))]
    document: dict[str, Any] = copy.deepcopy(gross)
    document["components"] = [c for c in gross["components"] + [methane] if c["component_id"] in names]
    document["pairs"] = [
        {"component_id_a": a, "component_id_b": b,
         "coefficients": [{"identity": f"{a}/{b}/k_ij", "family": "k_ij", "value": {"magnitude": 0.0, "unit": "dimensionless"}}]}
        for i, a in enumerate(names) for b in names[i + 1:]]
    if not polar:
        document["model_families"] = [f for f in document["model_families"] if f["kind"] != "polar"]
        for component in document["components"]:
            component["coefficients"] = [
                c for c in component["coefficients"] if c["family"] not in ("quadrupole_moment", "quadrupole_count")]
    return document


def _model(names: tuple[str, ...], polar: bool = True) -> Mixture | Pure:
    parameters = Parameters.from_mapping(_document(names, polar), components=names)
    return Pure(parameters) if len(names) == 1 else Mixture(parameters)


def _residual(model: Mixture | Pure, temperature: float, density: float, x: tuple[float, ...]) -> float:
    if isinstance(model, Pure):
        return model.state(temperature, rho=density).residual_helmholtz
    return model.state(temperature, rho=density, x=x).residual_helmholtz


def _quadrupolar(names: tuple[str, ...], temperature: float, density: float, x: tuple[float, ...]) -> float:
    return (_residual(_model(names, True), temperature, density, x)
            - _residual(_model(names, False), temperature, density, x))


def _rows(name: str) -> list[dict[str, str]]:
    with (RESULTS / name).open() as handle:
        return list(csv.DictReader(handle))


def test_quadrupolar_helmholtz_matches_clapeyron_and_feos() -> None:
    # Benzene (m = 2.2463 > 2) fails by 21.5 % if the dipole's chain-length cap reaches the quadrupole. FeOs caps
    # each m_i and carries longer J constants, so it is compared only where every m_i <= 2, at 1e-6.
    for name, column, count, tolerance in (("clapeyron-term-rows.csv", "a_qq_clapeyron", 4, 1e-10),
                                            ("feos-term-rows.csv", "a_qq_feos", 3, 1e-6)):
        rows = _rows(name)
        assert len(rows) == count
        for row in rows:
            names = tuple(row["species"].split("+"))
            x = tuple(float(v) for v in row["x"].split())
            value = _quadrupolar(names, float(row["T_K"]), float(row["rho_mol_m3"]), x)
            reference = float(row[column])
            assert abs(value - reference) <= tolerance * abs(reference), row


def test_quadrupolar_second_virial_limit() -> None:
    # Eq. 9 at eta = 0: a_qq / rho_N -> -(9 pi/16) sigma^3 q^2 J2(0), q = Q^2/(4 pi eps0 k T m sigma^5), m uncapped.
    record = {c["family"]: c["value"]["magnitude"] for c in next(
        c for c in json.loads(parameters_json("gross-2005", 1).read_text())["components"]
        if c["component_id"] == "benzene")["coefficients"]}
    m, sigma, epsilon, moment = (record[k] for k in ("segment_count", "segment_diameter", "dispersion_energy_over_k", "quadrupole_moment"))
    temperature = 400.0
    first, second = (m - 1.0) / m, (m - 1.0) * (m - 2.0) / m**2
    j2 = (1.2378308 + first * 1.2854109 + second * 1.7942954
          + (0.4542718 + first * -0.8137340 + second * 6.8682675) * epsilon / temperature)
    q = moment**2 * DEBYE**2 * 1e-20 / (4.0 * math.pi * EPS0 * BOLTZMANN * temperature * m * (sigma * 1e-10) ** 5)
    virial = -9.0 * math.pi / 16.0 * sigma**3 * q * q * j2
    errors = []
    for density in (1.0e-1, 1.0, 10.0):
        number_density = density * AVOGADRO * 1.0e-30
        errors.append(_quadrupolar(("benzene",), temperature, density, (1.0,)) / number_density / virial - 1.0)
    # The limit holds and the deviation is O(rho).
    intercept = errors[0] - (errors[1] - errors[0]) / 9.0
    assert abs(intercept) < 1e-8 and abs(errors[0]) > 1e-6
    for coarse, fine in zip(errors[1:], errors[:-1]):
        assert math.isclose(coarse / fine, 10.0, rel_tol=2e-2)


def _saturation_pressure(model: Mixture | Pure, component: str, temperature: float, guess: float) -> float:
    problem = Problem(phases=[Phase("L", kind="liquid", composition_guess=[1.0]),
                              Phase("V", kind="vapor", amount=Pinned(0.0), composition_guess=[1.0])],
                      T=temperature, P=Free(guess), feed=Amounts({component: 1.0}), row_tolerance=1.0e-12)
    result = solve_equilibrium(model, problem)
    assert result.success and result.pressure is not None, (component, temperature, result.message)
    return result.pressure


def test_saturation_pressure_reproduces_gross_2005() -> None:
    # Code identity with Clapeyron (1e-9), then the physical reproduction: CO2 vs Span-Wagner (paper 0.34 %, nonpolar
    # Gross-Sadowski 2001 3.01 % here) and benzene vs DIPPR-101 (paper 0.29 %; the m <= 2 cap gives 3.00 %).
    rows = _rows("clapeyron-saturation-rows.csv")
    deviations: dict[str, list[float]] = {"carbon-dioxide": [], "benzene": []}
    for row in rows:
        component, temperature, reference = row["species"], float(row["T_K"]), float(row["p_sat_clapeyron_Pa"])
        value = _saturation_pressure(_model((component,)), component, temperature, reference)
        assert abs(value / reference - 1.0) <= 1e-9, row
        if component == "carbon-dioxide":
            deviations[component].append(abs(value / float(row["p_sat_span_wagner_Pa"]) - 1.0))
        else:
            dippr = math.exp(83.107 - 6486.2 / temperature - 9.2194 * math.log(temperature) + 6.9844e-6 * temperature**2)
            deviations[component].append(abs(value / dippr - 1.0))
    assert len(deviations["carbon-dioxide"]) == 17 and len(deviations["benzene"]) == 26
    assert 100.0 * sum(deviations["carbon-dioxide"]) / 17 <= 0.5
    assert 100.0 * sum(deviations["benzene"]) / 26 <= 0.35


def test_co2_nitrogen_bubble_point_matches_clapeyron() -> None:
    # Code identity of the Gross (2005) mixing rule; the source finds it unsuitable for real mixtures with more
    # than one quadrupolar component, so this is not validation.
    model = _model(("carbon-dioxide", "nitrogen"))
    for row in _rows("clapeyron-bubble-rows.csv"):
        x = float(row["x_nitrogen"])
        problem = Problem(phases=[Phase("L", kind="liquid"), Phase("V", kind="vapor", amount=Pinned(0.0))],
                          T=float(row["T_K"]), P=Free(float(row["p_bubble_clapeyron_Pa"])),
                          feed=Amounts({"carbon-dioxide": 1.0 - x, "nitrogen": x}), row_tolerance=1.0e-12)
        result = solve_equilibrium(model, problem)
        assert result.success and result.pressure is not None, result.message
        assert abs(result.pressure / float(row["p_bubble_clapeyron_Pa"]) - 1.0) <= 1e-9, row
        assert abs(result.phase_mole_fractions(1)[1] / float(row["y_nitrogen_clapeyron"]) - 1.0) <= 1e-9, row


def test_quadrupolar_derivative_actions_and_identities() -> None:
    names, temperature, density, x = ("carbon-dioxide", "methane"), 240.0, 15000.0, (0.7, 0.3)
    model, nonpolar = _model(names, True), _model(names, False)
    assert isinstance(model, Mixture) and isinstance(nonpolar, Mixture)
    state = model.state(temperature, rho=density, x=x)
    raw = ((1.0, 0.0, (0.0, 0.0)), (0.0, 20.0, (0.0, 0.0)), (0.0, 0.0, (0.02, -0.02)))
    directions = tuple(DensityActionDirection(*item) for item in raw)

    def enthalpy(t: float, rho: float, z: tuple[float, ...]) -> float:
        return model.state(t, rho=rho, x=z).residual_enthalpy

    # Richardson (h, h/2) mixed stencils; predeclared relative criteria 1e-8, 1e-6, 1e-4.
    for order, tolerance in ((1, 1e-8), (2, 1e-6), (3, 1e-4)):
        batch = model.density_actions(temperature, density, x=x,
                                      observables=(PropertyObservable.ResidualEnthalpy,), directions=directions[:order])
        result = batch.results[0]
        assert result.status == PropertyActionStatus.Available and result.action is not None
        step = 0.2
        coarse = mixed_difference(enthalpy, temperature, density, x, raw[:order], step)
        fine = mixed_difference(enthalpy, temperature, density, x, raw[:order], step / 2)
        reference = (4.0 * fine - coarse) / 3.0
        assert abs(result.action - reference) <= tolerance * abs(reference), (order, result.action, reference)
        without = nonpolar.density_actions(temperature, density, x=x,
                                           observables=(PropertyObservable.ResidualEnthalpy,), directions=directions[:order])
        # The quadrupolar share of each action exceeds its criterion by many orders (non-vacuous).
        assert without.results[0].action is not None
        assert abs(result.action - without.results[0].action) > 1e2 * tolerance * abs(reference)

    # Z - 1 = rho a_rho and u_res = -R T^2 a_T through the quadrupolar chain.
    def slope(f: Any, center: float, h: float) -> float:
        return float(4.0 * (f(center + h / 2) - f(center - h / 2)) / h - (f(center + h) - f(center - h)) / (2 * h)) / 3.0

    a_rho = slope(lambda r: model.state(temperature, rho=r, x=x).residual_helmholtz, density, 1e-3 * density)
    assert abs(state.compressibility_factor - 1.0 - density * a_rho) < 1e-9
    a_t = slope(lambda t: model.state(t, rho=density, x=x).residual_helmholtz, temperature, 1e-3 * temperature)
    u = -8.31446261815324 * temperature**2 * a_t
    assert abs(state.residual_internal_energy - u) < 1e-8 * abs(u)


def test_absent_quadrupolar_species_limit() -> None:
    # With CO2 absent the declared family adds nothing, a_qq ~ A2 is quadratic in x_CO2, and an insertion of CO2
    # reports the dipole's support-change status (the A2-only branch omits A3).
    names, temperature, density = ("carbon-dioxide", "methane"), 240.0, 15000.0
    polar, nonpolar = _model(names, True), _model(names, False)
    assert isinstance(polar, Mixture) and isinstance(nonpolar, Mixture)
    with_family = polar.state(temperature, rho=density, x=(0.0, 1.0))
    without = nonpolar.state(temperature, rho=density, x=(0.0, 1.0))
    assert with_family.residual_helmholtz == without.residual_helmholtz
    assert with_family.pressure == without.pressure
    assert with_family.residual_chemical_potential_over_rt == without.residual_chemical_potential_over_rt
    small = [_quadrupolar(names, temperature, density, (s, 1.0 - s)) for s in (1e-4, 2e-4)]
    assert math.isclose(small[1] / small[0], 4.0, rel_tol=1e-3)
    insertion = (DensityActionDirection(0.0, 0.0, (1.0, -1.0)),) * 2
    for model, status in ((polar, PropertyActionStatus.SupportChange), (nonpolar, PropertyActionStatus.Available)):
        batch = model.density_actions(temperature, density, x=(0.0, 1.0),
                                      observables=(PropertyObservable.ResidualEnthalpy,), directions=insertion)
        assert batch.results[0].status == status
