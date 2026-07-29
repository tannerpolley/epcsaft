from dataclasses import FrozenInstanceError
import math
from pathlib import Path

import pytest
from pint import UnitRegistry

import epcsaft
from epcsaft import EosError, _ParameterBundle, ParameterError
from epcsaft import unit_registry
from epcsaft.eos import Mixture


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


def _model(component_ids: tuple[str, ...] = ("methane", "ethane")) -> Mixture:
    return Mixture(
        _ParameterBundle.from_catalog("gross-2001-methane-ethane", version=1).select(component_ids)
    )


def test_native_model_boundary_is_ordered_strict_and_immutable() -> None:
    bundle = _ParameterBundle.from_catalog("gross-2001-methane-ethane", version=1)
    parameters = bundle.select(("methane", "ethane"))
    model = Mixture(parameters)

    assert model.component_ids == ("methane", "ethane")
    assert model.parameter_fingerprint == parameters.fingerprint
    assert Mixture(bundle.select(("ethane", "methane"))).component_ids == (
        "ethane",
        "methane",
    )
    with pytest.raises(FrozenInstanceError):
        model._component_ids = ("ethane", "methane")

    fixture = _ParameterBundle.from_path(
        Path(__file__).parent / "data/bundles/figiel-2025-contract/1"
    )
    with pytest.raises(ParameterError) as exc_info:
        Mixture(fixture.select(("water",)))
    assert exc_info.value.code == "missing-record"
    assert dict(exc_info.value.context) == {
        "component_id": "water",
        "family": "segment_count",
    }

    assert epcsaft.__all__ == (
        "Parameters",
        "BundlePurpose",
        "ParameterError",
        "unit_registry",
        "Mixture",
        "State",
        "Fugacity",
        "DensityDiagnostics",
        "EosError",
        "Phase",
        "PropertyBasis",
        "native_sdk",
        "export_native_model",
        "ElectrolyteReference",
        "Activity",
        "TransferFreeEnergy",
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
    result = model.state(
        T=temperature * 1.8 * foreign.degree_Rankine,
        rho=density / 1000 * foreign.kilomole / foreign.meter**3,
        x=fractions,
    )

    assert isinstance(result, epcsaft.State)
    assert result.component_ids == model.component_ids
    assert result.mole_fractions == fractions
    assert result.temperature.to("kelvin").magnitude == pytest.approx(temperature)
    assert result.parameter_fingerprint == model.parameter_fingerprint
    assert result.hard_chain == pytest.approx(hard_chain, rel=5e-12, abs=5e-13)
    assert result.dispersion == pytest.approx(dispersion, rel=5e-12, abs=5e-13)
    assert result.association == 0.0
    assert result.born == 0.0
    assert not hasattr(result, "born_ssm_ds")
    assert result._residual_helmholtz == pytest.approx(residual, rel=5e-12, abs=5e-13)
    assert result.compressibility_factor == pytest.approx(z, rel=5e-12, abs=5e-13)
    assert result.pressure._REGISTRY is unit_registry
    assert result.pressure.to("pascal").magnitude == pytest.approx(pressure, rel=5e-12, abs=1e-6)

    step = density * 1.0e-3

    def residual_at(value: float) -> float:
        return model.state(
            T=temperature * unit_registry.kelvin,
            rho=value * unit_registry.mole / unit_registry.meter**3,
            x=fractions,
        )._residual_helmholtz

    derivative = (
        -residual_at(density + 2 * step)
        + 8 * residual_at(density + step)
        - 8 * residual_at(density - step)
        + residual_at(density - 2 * step)
    ) / (12 * step)
    reconstructed = (1 + density * derivative) * density * GAS_CONSTANT * temperature
    assert result.pressure.to("pascal").magnitude == pytest.approx(reconstructed, rel=2e-7)


def test_figiel_result_exposes_bulk_relative_permittivity() -> None:
    model = Mixture(
        _ParameterBundle.from_catalog("figiel-2025-reference-electrolytes", version=1).select(
            ("water", "sodium-cation", "chloride-anion")
        )
    )
    result = model.state(
        T=298.15 * unit_registry.kelvin,
        rho=40.0 * unit_registry.mole / unit_registry.meter**3,
        x=(0.9, 0.05, 0.05),
    )

    assert model.relative_permittivity_formulation == "ion-fraction-suppression"
    assert result.bulk_relative_permittivity == pytest.approx(78.09 / (1.0 + 7.01 * 0.1))


def test_evaluate_is_permutation_invariant() -> None:
    state = {
        "T": 250 * unit_registry.kelvin,
        "rho": 8000 * unit_registry.mole / unit_registry.meter**3,
    }
    forward = _model().state(**state, x=(0.4, 0.6))
    reverse = _model(("ethane", "methane")).state(**state, x=(0.6, 0.4))

    for field in (
        "hard_chain",
        "dispersion",
        "_residual_helmholtz",
        "compressibility_factor",
    ):
        assert getattr(reverse, field) == pytest.approx(getattr(forward, field))
    assert reverse.pressure.to("pascal").magnitude == pytest.approx(
        forward.pressure.to("pascal").magnitude
    )


def test_evaluate_rejects_invalid_states_and_has_physical_limits() -> None:
    model = _model()
    valid = {
        "T": 250 * unit_registry.kelvin,
        "rho": 8000 * unit_registry.mole / unit_registry.meter**3,
        "x": (0.4, 0.6),
    }
    invalid_states = (
        {"T": 250.0},
        {"rho": 8000.0},
        {"T": 250 * unit_registry.meter},
        {"rho": 8000 * unit_registry.kelvin},
        {"T": 0 * unit_registry.kelvin},
        {"T": -1 * unit_registry.kelvin},
        {"T": math.nan * unit_registry.kelvin},
        {"T": math.inf * unit_registry.kelvin},
        {"T": 50 * unit_registry.kelvin},
        {"T": 301 * unit_registry.kelvin},
        {"rho": 0 * unit_registry.mole / unit_registry.meter**3},
        {"rho": -1 * unit_registry.mole / unit_registry.meter**3},
        {"rho": math.nan * unit_registry.mole / unit_registry.meter**3},
        {"rho": math.inf * unit_registry.mole / unit_registry.meter**3},
        {"x": (1.0,)},
        {"x": (-0.1, 1.1)},
        {"x": (math.nan, math.nan)},
        {"x": (0.4, 0.600000000002)},
    )
    for override in invalid_states:
        state = {**valid, **override}
        with pytest.raises(EosError):
            model.state(**state)

    low_density = 1.0e-6
    low = model.state(
        T=250 * unit_registry.kelvin,
        rho=low_density * unit_registry.mole / unit_registry.meter**3,
        x=(0.4, 0.6),
    )
    assert abs(low._residual_helmholtz) < 1.0e-8
    assert abs(low.compressibility_factor - 1.0) < 1.0e-8
    assert low.pressure.to("pascal").magnitude == pytest.approx(
        low_density * GAS_CONSTANT * 250, rel=1.0e-8
    )

    temperature = 250.0
    diameter = 3.7039 * (1 - 0.12 * math.exp(-3 * 150.03 / temperature))
    close_packing = math.pi / (3 * math.sqrt(2))
    close_density = close_packing * 6 / (math.pi * (AVOGADRO / 1.0e30) * diameter**3)
    with pytest.raises(EosError) as exc_info:
        model.state(
            T=temperature * unit_registry.kelvin,
            rho=(close_density * (1 + 1.0e-12))
            * unit_registry.mole
            / unit_registry.meter**3,
            x=(1.0, 0.0),
        )
    assert exc_info.value.code == "evaluation-failed"
    assert exc_info.value.context


@pytest.mark.parametrize(
    "override,expected_context",
    (
        (
            {"T": 250.0},
            {
                    "field": "T",
                "reason": "invalid-unit",
                "expected_unit": "kelvin",
            },
        ),
        (
            {"T": 50 * unit_registry.kelvin},
            {
                    "field": "T",
                "reason": "outside-domain",
                "value": "50.0",
                "minimum": "97.0",
                "maximum": "300.0",
                "unit": "kelvin",
            },
        ),
        (
            {"x": (0.4, 0.61)},
            {
                    "field": "x",
                "reason": "not-normalized",
                "value": "1.01",
                "expected": "1.0",
                "tolerance": "1e-12",
            },
        ),
        (
            {"rho": None, "P": 1 * unit_registry.meter},
            {
                    "field": "P",
                "reason": "invalid-unit",
                "expected_unit": "pascal",
            },
        ),
        (
            {
                "rho": -1 * unit_registry.mole / unit_registry.meter**3,
            },
            {
                    "field": "rho",
                "reason": "nonpositive",
                "value": "-1.0",
            },
        ),
        (
            {"T": float("nan") * unit_registry.kelvin},
            {
                    "field": "T",
                "reason": "nonfinite",
                "value": "nan",
            },
        ),
    ),
)
def test_invalid_state_errors_explain_the_rejected_input(
    override: dict[str, object], expected_context: dict[str, str]
) -> None:
    state = {
        "T": 250 * unit_registry.kelvin,
        "rho": 8000 * unit_registry.mole / unit_registry.meter**3,
        "x": (0.4, 0.6),
        **override,
    }

    with pytest.raises(EosError) as exc_info:
        _model().state(**state)

    assert exc_info.value.code == "invalid-state"
    assert dict(exc_info.value.context) == expected_context


def test_eos_error_contract_is_stable_and_immutable() -> None:
    error = EosError("invalid-state", (("field", "temperature"),))
    assert isinstance(error, ValueError)
    assert error.code == "invalid-state"
    assert error.context == (("field", "temperature"),)
    assert str(error) == "invalid-state (field=temperature)"
