from dataclasses import FrozenInstanceError
from pathlib import Path

import pytest

import epcsaft
from epcsaft import ParameterBundle, ParameterError
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


def test_native_model_boundary_is_ordered_strict_and_immutable() -> None:
    bundle = ParameterBundle.from_catalog(
        "gross-2001-methane-ethane", version=1
    )
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
    assert exc_info.value.code == "invalid-value"
    assert dict(exc_info.value.context) == {
        "record_id": "water-association-energy",
        "family": "association_energy_over_k",
    }

    assert epcsaft.__all__ == (
        "ParameterBundle",
        "ParameterSet",
        "ParameterError",
        "unit_registry",
        "EPCSAFT",
        "EosResult",
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
    parameters = ParameterBundle.from_catalog(
        "gross-2001-methane-ethane", version=1
    ).select(("methane", "ethane"))
    result = EPCSAFT(parameters).evaluate(
        temperature=temperature * unit_registry.kelvin,
        molar_density=density * unit_registry.mole / unit_registry.meter**3,
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
