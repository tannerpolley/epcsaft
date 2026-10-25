"""#27: admitted density-dependent permittivity forms (full EOS 6.2, 6.4).

The selections are the retained Figiel-base hybrids of
analyses/2026-dielectric-alternatives; each exchanges only the permittivity
formulation. The checks falsify packet-to-native unit conversion, the ion-free
limit, the single epsilon_r consumed by DH and Born, and the epsilon_r(T, rho)
derivative chain through pressure and internal energy.
"""

import copy
import json
import math
import runpy
from collections.abc import Callable, Sequence
from pathlib import Path
from typing import Any

import pytest

from epcsaft import GAS_CONSTANT_J_PER_MOL_K, ActiveParameter, Mixture, ParameterFamily, Parameters, _core, regression
from epcsaft.equilibrium import (
    Amounts, NeutralReference, Phase, Problem, Reaction, ReactionLogPolynomial, ReactionReference,
    ReferenceSourceBasis, compile_problem, solve_equilibrium,
)

SELECTIONS = Path(__file__).resolve().parents[3] / "analyses/2026-dielectric-alternatives/inputs/selections"
COMPONENTS = ("water", "sodium-cation", "chloride-anion")
T, RHO, X = 298.15, 55000.0, (0.96, 0.02, 0.02)
AVOGADRO, BOLTZMANN = 6.02214076e23, 1.380649e-23
# Debye^2 in J A^3 from the SI Debye (1e-21/c C m), CODATA 2018 epsilon_0 and the angstrom.
DEBYE_SQUARED_J_A3 = (1.0e-21 / 299792458.0) ** 2 / (4.0 * math.pi * 8.8541878128e-12 * 1.0e-30)


def test_component_permittivity_zero_water() -> None:
    """Bülow 2019 Eq. (7), including support transitions on the same native model."""
    declarations = runpy.run_path(str(SELECTIONS.parents[2] / "reproduction/2019-bulow/scripts/reproduce.py"))["declarations"]
    model, _ = declarations()
    for path in ((0.5, 0.0, 0.5), (0.0, 0.5, 0.0), (1.0, 0.9, 1.0e-8, 1.0e-12)):
        for t in path:
            x = (t, (1.0 - t) / 2.0, (1.0 - t) / 2.0)
            state = model.state(298.15, rho=1000.0, x=x)
            reference = {0.0: 11.0, 1.0: 80.0, 0.9: 73.1, 0.5: 45.5}.get(t)
            assert math.isfinite(state.bulk_relative_permittivity)
            if reference is not None:
                assert abs(state.bulk_relative_permittivity - reference) <= 1.0e-12
            fresh, _ = declarations()
            expected = fresh.state(298.15, rho=1000.0, x=x)
            for quantity in ("bulk_relative_permittivity", "pressure", "residual_internal_energy", "compressibility_factor"):
                assert getattr(state, quantity).hex() == getattr(expected, quantity).hex()


def _document(name: str) -> dict[str, Any]:
    document: dict[str, Any] = json.loads((SELECTIONS / f"{name}.json").read_text())
    return document


def _record(document: dict[str, Any], component: str, family: str) -> float:
    item = next(c for c in document["components"] if c["component_id"] == component)
    return float(next(r["value"]["magnitude"] for r in item["coefficients"] if r["family"] == family))


def _pure_water_permittivity(name: str, document: dict[str, Any], temperature: float, density: float) -> float:
    """Full EOS 6.2 / 6.4 written independently from the packet records (SI inputs)."""
    if name == "density-temperature":
        volume = _record(document, "water", "dielectric_volume") * 1.0e-3  # liter/mole -> m3/mol
        return 1.0 + density * volume * (_record(document, "water", "dielectric_temperature") / temperature - 1.0)
    number_density = density * AVOGADRO * 1.0e-30  # 1/A^3
    y = 4.0 * math.pi / 9.0 * number_density * (
        _record(document, "water", "dipole_scaling") * DEBYE_SQUARED_J_A3 / (BOLTZMANN * temperature)
        + 3.0 * _record(document, "water", "polarizability_scaling"))
    integral = 1.0 + _record(document, "water", "correlation_integral_parameter") * (math.exp(-y) - 1.0)
    return 1.0 + 3.0 * y * (1.0 + y + (17.0 / 16.0 * integral - 1.0) * y * y)


def _constant_permittivity_model(epsilon: float) -> Mixture:
    """Figiel base with a composition-independent water epsilon_r (solvent-only form)."""
    document = copy.deepcopy(_document("ion-fraction-suppression"))
    next(f for f in document["model_families"] if f["kind"] == "permittivity")["choice"] = "solvent-only"
    document["model_coefficients"] = [
        m for m in document["model_coefficients"] if m["family"] != "ion_fraction_suppression_coefficient"]
    water = next(c for c in document["correlations"]
                 if c["component_id"] == "water" and c["family"] == "relative_permittivity")
    water["constant"]["value"]["magnitude"] = epsilon
    return Mixture(Parameters.from_mapping(document, components=COMPONENTS))


def _richardson(value: Callable[[float], float], center: float, step: float) -> float:
    coarse = (value(center + step) - value(center - step)) / (2.0 * step)
    fine = (value(center + step / 2) - value(center - step / 2)) / step
    return float(4.0 * fine - coarse) / 3.0


@pytest.mark.parametrize("name", ["density-temperature", "dipolar-perturbation"])
def test_density_dependent_permittivity_limits_and_chains(name: str) -> None:
    document = _document(name)
    model = Mixture(Parameters.from_mapping(document, components=COMPONENTS))

    # Ion-free limit reproduces the pure-solvent relation on the admitted units.
    for density in (RHO, 20.0):
        limit = model.state(T, rho=density, x=(1.0, 0.0, 0.0)).bulk_relative_permittivity
        assert math.isclose(limit, _pure_water_permittivity(name, document, T, density), rel_tol=1e-13)

    state = model.state(T, rho=RHO, x=X)
    epsilon = state.bulk_relative_permittivity
    assert 1.0 < epsilon < _pure_water_permittivity(name, document, T, RHO)

    # DH and Born consume the same epsilon_r: a constant-epsilon model at that value agrees.
    constant = _constant_permittivity_model(epsilon).state(T, rho=RHO, x=X)
    assert math.isclose(constant.bulk_relative_permittivity, epsilon, rel_tol=1e-14)
    assert math.isclose(state.debye_huckel, constant.debye_huckel, rel_tol=1e-12)
    assert math.isclose(state.born, constant.born, rel_tol=1e-12)
    assert math.isclose(state.residual_helmholtz, constant.residual_helmholtz, rel_tol=1e-12)

    # epsilon_r(rho) and epsilon_r(T) chains: Z - 1 = rho a_rho, u_res = -R T^2 a_T.
    # Richardson (h, h/2) differences resolve both to <1e-11 relative here; the
    # epsilon_r chain moves Z by ~0.03 and u_res by ~75 J/mol (constant-epsilon model),
    # >1e5 times each criterion, so a dropped chain cannot pass.
    density_slope = _richardson(lambda r: model.state(T, rho=r, x=X).residual_helmholtz, RHO, 1.0e-3 * RHO)
    assert abs(state.compressibility_factor - 1.0 - RHO * density_slope) < 1e-9
    assert abs(state.compressibility_factor - constant.compressibility_factor) > 1e-2
    temperature_slope = _richardson(lambda t: model.state(t, rho=RHO, x=X).residual_helmholtz, T, 1.0e-3 * T)
    internal_energy = -GAS_CONSTANT_J_PER_MOL_K * T * T * temperature_slope
    assert abs(state.residual_internal_energy - internal_energy) < 1e-9 * abs(internal_energy)
    assert abs(state.residual_internal_energy - constant.residual_internal_energy) > 1e-3 * abs(internal_energy)


DIELECTRIC = [ActiveParameter(ParameterFamily.DielectricVolume, (0,)), ActiveParameter(ParameterFamily.DielectricTemperature, (0,))]
DIELECTRIC_STEPS = (3.777e-4, 1.403)  # 1e-3 of water's d_V (L/mol) and d_T (K)


def _dielectric(volume: float = 0.0, temperature: float = 0.0) -> Parameters:
    """The density-temperature selection with water's d_V (L/mol) and d_T (K) moved by the given amounts."""
    document = _document("density-temperature")
    water = next(c for c in document["components"] if c["component_id"] == "water")
    for item in water["coefficients"]:
        item["value"]["magnitude"] += {"dielectric_volume": volume, "dielectric_temperature": temperature}.get(item["family"], 0.0)
    return Parameters.from_mapping(document, components=COMPONENTS)


def _rebuilt_column(value: Callable[[Parameters], list[float]], column: int) -> list[float]:
    """Richardson (h, h/2) differences of rebuilt-model values, per m3/mol for d_V and per K for d_T."""
    h = DIELECTRIC_STEPS[column]
    plus, minus, half_plus, half_minus = (value(_dielectric(*((s * h, 0.0) if column == 0 else (0.0, s * h))))
                                          for s in (1.0, -1.0, 0.5, -0.5))
    unit = 1.0e3 if column == 0 else 1.0
    return [unit * (4.0 * (c - d) / h - (a - b) / (2.0 * h)) / 3.0 for a, b, c, d in zip(plus, minus, half_plus, half_minus)]


def _numbers(values: Sequence[float | None]) -> list[float]:
    assert all(value is not None for value in values)
    return [float(value) for value in values if value is not None]


def test_dielectric_coefficient_columns_and_admission() -> None:
    """#441: at fixed rho, Yang 2024 Eqs. 2a/2b give d eps/d d_V,w = rho x_w (d_T/T - 1) and d eps/d d_T,w = rho x_w d_V/T.
    The pressure-closed permittivity column less its density chain (eps - 1)/rho drho/dtheta equals them to 1e-12; a
    missing 1e-3 L/mol factor or a swapped slot misses by orders of magnitude. The ln(phi) columns, which carry the
    pressure closure, equal rebuilt-model differences (resolved to ~2e-9). An ion slot or another form is refused."""
    parameters = _dielectric()
    model = Mixture(parameters, active_parameters=DIELECTRIC)
    x = dict(zip(COMPONENTS, X))
    permittivity = regression.electrolyte(parameters, "bulk_relative_permittivity", T, 1.0e5, x, ions=(COMPONENTS[1], COMPONENTS[2]),
                                          stoichiometry=(1, 1), molality=0.0, observed=1.0, form="difference", scale=1.0)
    density = regression.density(parameters, T, 1.0e5, x, phase="liquid", observed=1.0, form="difference", mass_basis=False, scale=1.0)
    epsilon, rho = _numbers([regression.predict(model, item).value for item in (permittivity, density)])
    d_epsilon, d_rho = (regression.differentiate(model, item).value for item in (permittivity, density))
    volume, temperature = 1.0e-3 * _record(_document("density-temperature"), "water", "dielectric_volume"), _record(
        _document("density-temperature"), "water", "dielectric_temperature")
    for column, analytic in enumerate((rho * X[0] * (temperature / T - 1.0), rho * X[0] * volume / T)):
        assert abs(d_epsilon[column] - (epsilon - 1.0) / rho * d_rho[column] - analytic) <= 1e-12 * abs(analytic)

    exact = _numbers(model.state(T, P=1.0e5, x=X, phase="liquid").log_fugacity_parameter_derivative)
    for column in range(2):
        rebuilt = _rebuilt_column(lambda p: _numbers(Mixture(p).state(T, P=1.0e5, x=X, phase="liquid").log_fugacity_coefficient or []), column)
        for species, difference in enumerate(rebuilt):
            assert abs(exact[species * 2 + column] - difference) <= 1e-8 * abs(difference)

    suppressed = Parameters.from_mapping(_document("ion-fraction-suppression"), components=COMPONENTS)
    for source, family, slot in ((parameters, ParameterFamily.DielectricTemperature, 1), (suppressed, ParameterFamily.DielectricVolume, 0)):
        with pytest.raises(ValueError, match="dielectric coordinate"):
            Mixture(source, active_parameters=[ActiveParameter(family, (slot,))])


def test_dielectric_coefficients_recover_in_one_process_and_with_workers() -> None:
    """#441: water d_V (an L/mol coordinate, target coefficient 1e-3) and d_T recover to 1e-6 from manufactured bulk
    permittivity at two temperatures, brine densities at both (two non-electrolyte state groups, the second sent to a
    worker process) and ln gamma/phi at 0.1, 1 and 4.5 mol/kg. One process and two agree to 1e-10; every Jacobian column
    equals differences of rebuilt models and references in optimizer coordinates (resolved to ~1e-8)."""
    truth, ions, brine = _dielectric(), (COMPONENTS[1], COMPONENTS[2]), dict(zip(COMPONENTS, X))
    templates = [regression.electrolyte(truth, "bulk_relative_permittivity", t, 1.0e5, brine, ions=ions, stoichiometry=(1, 1),
                                        molality=0.0, observed=1.0, form="log_ratio") for t in (298.15, 323.15)]
    templates += [regression.density(truth, t, 1.0e5, brine, phase="liquid", observed=1.0, form="log_ratio", mass_basis=False)
                  for t in (298.15, 323.15)]
    templates += [regression.electrolyte(truth, kind, 298.15, 1.0e5, {"water": 1.0}, ions=ions, stoichiometry=(1, 1), molality=m,
                                         observed=1.0, form="log_ratio")
                  for m in (0.1, 1.0, 4.5) for kind in ("mean_ionic_activity_coefficient", "osmotic_coefficient")]
    observations = _core.observations_at(truth._native, [], [], [], templates, [0.0] * len(templates))
    start, scales = _dielectric(-0.0177, -53.0), (0.05, 100.0)
    coordinates = [
        regression.coordinate(start, [("dielectric_volume", "water", 1.0e-3)], origin=0.35, scale=scales[0], bounds=(0.2, 0.6)),
        regression.coordinate(start, [("dielectric_temperature", "water", 1.0)], origin=1400.0, scale=scales[1], bounds=(900.0, 2200.0))]
    serial, parallel = (regression.fit(start, coordinates, observations, controls=regression.FitControls(
        maximum_iterations=100, maximum_elapsed_time_seconds=600.0, function_tolerance=1e-15, gradient_tolerance=1e-14,
        parameter_tolerance=1e-14, workers=workers)) for workers in (1, 2))
    assert serial.status == "converged", serial.message
    assert abs(serial.physical[0] / 3.777e-4 - 1.0) <= 1e-6 and abs(serial.physical[1] / 1403.0 - 1.0) <= 1e-6
    for one, two in ((serial.physical, parallel.physical), (serial.optimizer_jacobian, parallel.optimizer_jacobian)):
        assert all(math.isclose(a, b, rel_tol=1e-10) for a, b in zip(one, two, strict=True))

    for column, step in enumerate(h / scale for h, scale in zip(DIELECTRIC_STEPS, scales)):
        def residuals(shift: float) -> list[float]:
            point = list(serial.optimizer)
            point[column] += shift
            chart = regression.evaluate_coordinates(start, coordinates, point, [])
            values = _core.observations_at(chart.parameters, [], [], [], templates, [0.0] * len(templates))
            return [math.log(value.observed / item.observed) for value, item in zip(values, observations)]
        plus, minus, half_plus, half_minus = (residuals(s * step) for s in (1.0, -1.0, 0.5, -0.5))
        for row in range(len(templates)):
            difference = (4.0 * (half_plus[row] - half_minus[row]) / step - (plus[row] - minus[row]) / (2.0 * step)) / 3.0
            assert abs(serial.optimizer_jacobian[row * 2 + column] - difference) <= 1e-7 * abs(difference)


def test_dielectric_tangent_at_a_solved_ionic_state() -> None:
    """#441: at the referenced water = Na+ + Cl- state (K = 1e-8), whose infinite-dilution reference carries its own
    d_V and d_T columns, R_z z_theta + R_theta = 0 with R_theta from rebuilt models and references at the solved
    point, and z_theta equals centred re-solves, each to 1e-8 of its scale (differences resolve both to ~2e-10)."""
    law = ReactionLogPolynomial(math.log(1e-8), 0.0, 0.0, 0.0, 298.15, 250.0, 350.0, "analytic-source")
    problem = Problem(
        T=298.15, P=1.0e5, feed=Amounts({"water": 1.0, "sodium-cation": 1e-6, "chloride-anion": 1e-6}), phases=[Phase("L", kind="liquid")],
        neutral_reference=NeutralReference(solvent="water", solute_molality_ratios={"sodium-cation": 1.0, "chloride-anion": 1.0},
                                           standard_molality_mol_per_kg=1.0),
        reactions=[Reaction({"water": -1.0, "sodium-cation": 1.0, "chloride-anion": 1.0}, correlation=law,
                            reference=ReactionReference(ReferenceSourceBasis.RawMoleFractionInfiniteDilution, None))])

    def solved(parameters: Parameters) -> list[float]:
        result = solve_equilibrium(Mixture(parameters, active_parameters=DIELECTRIC), problem)
        assert result.success, result.message
        return list(result.coordinates)

    model = Mixture(_dielectric(), active_parameters=DIELECTRIC)
    point = solved(_dielectric())
    compiled = compile_problem(model, problem)
    jacobian, width = compiled.jacobian(point), len(point)
    for column in range(2):
        tangent = _core.sensitivity(compiled, point, [0.0] * 5 + [float(column == k) for k in range(2)])
        rows = _rebuilt_column(lambda p: list(compile_problem(Mixture(p, active_parameters=DIELECTRIC), problem).evaluate(point).values), column)
        for row, r_theta in enumerate(rows):
            assert abs(sum(jacobian[row * width + j] * tangent[j] for j in range(width)) + r_theta) <= 1e-8 * max(map(abs, rows))
        resolved = _rebuilt_column(solved, column)
        assert max(abs(a - b) for a, b in zip(tangent, resolved)) <= 1e-8 * max(map(abs, tangent))
