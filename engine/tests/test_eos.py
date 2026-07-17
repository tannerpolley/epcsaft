from dataclasses import FrozenInstanceError
import math
from pathlib import Path

import pytest
from pint import UnitRegistry

import epcsaft
from epcsaft import EosError, ParameterBundle, ParameterError
from epcsaft import unit_registry
from epcsaft.eos import EPCSAFT


GOLDENS = (
    # T [K], rho_m [mol/m3], x, a_hc, a_disp, a_res, Z, P [Pa]
    (
        200.0,
        1000.0,
        (1.0, 0.0),
        0.06290111429923406,
        -0.16906790398822233,
        -0.10616678968898827,
        0.8963961463906834,
        1490610.4500443912,
    ),
    (
        250.0,
        1000.0,
        (0.0, 1.0),
        0.10804736768399885,
        -0.3566531475792243,
        -0.24860577989522545,
        0.7598730786270726,
        1579484.0766964532,
    ),
    (
        250.0,
        8000.0,
        (0.4, 0.6),
        0.8611345413201740,
        -1.9282055648091894,
        -1.0670710234890153,
        0.19067336346364548,
        3170693.1055920523,
    ),
)

AVOGADRO = 6.02214076e23
GAS_CONSTANT = 8.31446261815324


def _model(component_ids: tuple[str, ...] = ("methane", "ethane")) -> EPCSAFT:
    return EPCSAFT(
        ParameterBundle.from_catalog(
            "gross-2001-methane-ethane", version=1
        ).select(component_ids)
    )


def test_native_model_boundary_is_ordered_strict_and_immutable() -> None:
    bundle = ParameterBundle.from_catalog("gross-2001-methane-ethane", version=1)
    parameters = bundle.select(("methane", "ethane"))
    model = EPCSAFT(parameters)

    assert model.component_ids == ("methane", "ethane")
    assert model.parameter_fingerprint == parameters.fingerprint
    assert EPCSAFT(bundle.select(("ethane", "methane"))).component_ids == (
        "ethane",
        "methane",
    )
    with pytest.raises(FrozenInstanceError):
        model._component_ids = ("ethane", "methane")

    fixture = ParameterBundle.from_path(
        Path(__file__).parent / "data/bundles/figiel-2025-contract/1"
    )
    with pytest.raises(ParameterError) as exc_info:
        EPCSAFT(fixture.select(("water",)))
    assert exc_info.value.code == "missing-record"
    assert dict(exc_info.value.context) == {
        "component_id": "water",
        "family": "segment_count",
    }

    assert epcsaft.__all__ == (
        "ParameterBundle",
        "ParameterSet",
        "ParameterError",
        "unit_registry",
        "EPCSAFT",
        "EosResult",
        "EosError",
    )


@pytest.mark.parametrize(
    "temperature,density,fractions,hard_chain,dispersion,residual,z,pressure",
    GOLDENS,
)
def test_gross_hard_chain_dispersion_and_pressure_goldens(
    temperature: float,
    density: float,
    fractions: tuple[float, float],
    hard_chain: float,
    dispersion: float,
    residual: float,
    z: float,
    pressure: float,
) -> None:
    model = _model()
    foreign = UnitRegistry()
    result = model.evaluate(
        temperature=temperature * 1.8 * foreign.degree_Rankine,
        molar_density=density / 1000 * foreign.kilomole / foreign.meter**3,
        mole_fractions=fractions,
    )

    assert isinstance(result, epcsaft.EosResult)
    assert result.hard_chain == pytest.approx(hard_chain, rel=5e-12, abs=5e-13)
    assert result.dispersion == pytest.approx(dispersion, rel=5e-12, abs=5e-13)
    assert result.residual_helmholtz == pytest.approx(
        residual, rel=5e-12, abs=5e-13
    )
    assert result.compressibility_factor == pytest.approx(
        z, rel=5e-12, abs=5e-13
    )
    assert result.pressure._REGISTRY is unit_registry
    assert result.pressure.to("pascal").magnitude == pytest.approx(
        pressure, rel=5e-12, abs=1e-6
    )

    step = density * 1.0e-3

    def residual_at(value: float) -> float:
        return model.evaluate(
            temperature=temperature * unit_registry.kelvin,
            molar_density=value * unit_registry.mole / unit_registry.meter**3,
            mole_fractions=fractions,
        ).residual_helmholtz

    derivative = (
        -residual_at(density + 2 * step)
        + 8 * residual_at(density + step)
        - 8 * residual_at(density - step)
        + residual_at(density - 2 * step)
    ) / (12 * step)
    reconstructed = (
        (1 + density * derivative) * density * GAS_CONSTANT * temperature
    )
    assert result.pressure.to("pascal").magnitude == pytest.approx(
        reconstructed, rel=2e-7
    )


def test_evaluate_is_permutation_invariant() -> None:
    state = {
        "temperature": 250 * unit_registry.kelvin,
        "molar_density": 8000 * unit_registry.mole / unit_registry.meter**3,
    }
    forward = _model().evaluate(**state, mole_fractions=(0.4, 0.6))
    reverse = _model(("ethane", "methane")).evaluate(
        **state, mole_fractions=(0.6, 0.4)
    )

    for field in (
        "hard_chain",
        "dispersion",
        "residual_helmholtz",
        "compressibility_factor",
    ):
        assert getattr(reverse, field) == pytest.approx(getattr(forward, field))
    assert reverse.pressure.to("pascal").magnitude == pytest.approx(
        forward.pressure.to("pascal").magnitude
    )


def test_evaluate_rejects_invalid_states_and_has_physical_limits() -> None:
    model = _model()
    valid = {
        "temperature": 250 * unit_registry.kelvin,
        "molar_density": 8000 * unit_registry.mole / unit_registry.meter**3,
        "mole_fractions": (0.4, 0.6),
    }
    invalid_states = (
        {"temperature": 250.0},
        {"molar_density": 8000.0},
        {"temperature": 250 * unit_registry.meter},
        {"molar_density": 8000 * unit_registry.kelvin},
        {"temperature": 0 * unit_registry.kelvin},
        {"temperature": -1 * unit_registry.kelvin},
        {"temperature": math.nan * unit_registry.kelvin},
        {"temperature": math.inf * unit_registry.kelvin},
        {"molar_density": 0 * unit_registry.mole / unit_registry.meter**3},
        {"molar_density": -1 * unit_registry.mole / unit_registry.meter**3},
        {
            "molar_density": math.nan
            * unit_registry.mole
            / unit_registry.meter**3
        },
        {
            "molar_density": math.inf
            * unit_registry.mole
            / unit_registry.meter**3
        },
        {"mole_fractions": (1.0,)},
        {"mole_fractions": (-0.1, 1.1)},
        {"mole_fractions": (math.nan, math.nan)},
        {"mole_fractions": (0.4, 0.600000000002)},
    )
    for override in invalid_states:
        state = {**valid, **override}
        with pytest.raises(EosError):
            model.evaluate(**state)

    low_density = 1.0e-6
    low = model.evaluate(
        temperature=250 * unit_registry.kelvin,
        molar_density=low_density * unit_registry.mole / unit_registry.meter**3,
        mole_fractions=(0.4, 0.6),
    )
    assert abs(low.residual_helmholtz) < 1.0e-8
    assert abs(low.compressibility_factor - 1.0) < 1.0e-8
    assert low.pressure.to("pascal").magnitude == pytest.approx(
        low_density * GAS_CONSTANT * 250, rel=1.0e-8
    )

    temperature = 250.0
    diameter = 3.7039 * (1 - 0.12 * math.exp(-3 * 150.03 / temperature))
    close_packing = math.pi / (3 * math.sqrt(2))
    close_density = close_packing * 6 / (
        math.pi * (AVOGADRO / 1.0e30) * diameter**3
    )
    with pytest.raises(EosError) as exc_info:
        model.evaluate(
            temperature=temperature * unit_registry.kelvin,
            molar_density=(close_density * (1 + 1.0e-12))
            * unit_registry.mole
            / unit_registry.meter**3,
            mole_fractions=(1.0, 0.0),
        )
    assert exc_info.value.code == "evaluation-failed"
    assert exc_info.value.context


def test_eos_error_contract_is_stable_and_immutable() -> None:
    error = EosError("invalid-state", (("field", "temperature"),))
    assert isinstance(error, ValueError)
    assert error.code == "invalid-state"
    assert error.context == (("field", "temperature"),)
    assert str(error) == "invalid-state (field=temperature)"
