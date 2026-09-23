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
from collections.abc import Callable
from pathlib import Path
from typing import Any

import pytest

from epcsaft import GAS_CONSTANT_J_PER_MOL_K, Mixture, Parameters

SELECTIONS = Path(__file__).resolve().parents[3] / "analyses/2026-dielectric-alternatives/inputs/selections"
COMPONENTS = ("water", "sodium-cation", "chloride-anion")
T, RHO, X = 298.15, 55000.0, (0.96, 0.02, 0.02)
AVOGADRO, BOLTZMANN = 6.02214076e23, 1.380649e-23


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
        _record(document, "water", "dipole_scaling") * 1.0e-19 / (BOLTZMANN * temperature)
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
