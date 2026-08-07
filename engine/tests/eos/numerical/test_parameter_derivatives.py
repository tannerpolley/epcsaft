from __future__ import annotations

from collections.abc import Callable
from copy import deepcopy

import pytest

from epcsaft import ActiveParameterSet, Mixture, Parameters, Pure
from epcsaft import unit_registry as u
from tests.eos._raw_cases import (
    AQUEOUS_ELECTROLYTE,
    ASSOCIATING_PURE,
    POLAR_BINARY,
    POLAR_PURE,
    _component,
    _pair,
    neutral_pure,
)
from tests.equilibrium.parameter_dictionaries import PURE_SATURATION_PARAMETERS

Values = tuple[float, ...]


def _set_parameter(mapping: dict[str, object], identity: str, value: float) -> None:
    pending: list[object] = [mapping]
    while pending:
        item = pending.pop()
        if isinstance(item, dict):
            if item.get("identity") == identity:
                quantity = item["value"]
                assert isinstance(quantity, dict)
                quantity["magnitude"] = value
                return
            pending.extend(item.values())
        elif isinstance(item, list):
            pending.extend(item)
    raise AssertionError(f"parameter identity not found: {identity}")


def _parameter_value(mapping: dict[str, object], identity: str) -> float:
    parameters = Parameters.from_mapping(mapping)
    value = next(spec.value for spec in parameters.parameter_specs if spec.identity == identity)
    return float(getattr(value, "magnitude", value))


def _perturbed(
    mapping: dict[str, object], identity: str, value: float
) -> dict[str, object]:
    changed = deepcopy(mapping)
    _set_parameter(changed, identity, value)
    return changed


def _central_difference(
    mapping: dict[str, object],
    identity: str,
    step: float,
    evaluate: Callable[[dict[str, object]], Values],
) -> Values:
    value = _parameter_value(mapping, identity)
    lower = evaluate(_perturbed(mapping, identity, value - step))
    upper = evaluate(_perturbed(mapping, identity, value + step))
    return tuple((high - low) / (2.0 * step) for low, high in zip(lower, upper, strict=True))


def _fixed_state_values(
    mapping: dict[str, object],
    temperature_k: float,
    density: float,
    composition: tuple[float, ...],
) -> Values:
    parameters = Parameters.from_mapping(mapping)
    model = Pure(parameters) if len(composition) == 1 else Mixture(parameters)
    state = (
        model.state(T=temperature_k * u.kelvin, rho=density * u.mole / u.meter**3)
        if isinstance(model, Pure)
        else model.state(
            T=temperature_k * u.kelvin,
            rho=density * u.mole / u.meter**3,
            x=composition,
        )
    )
    return (
        float(state.ares(basis="reduced")),
        float(state.pressure.to("pascal").magnitude),
        *state.residual_chemical_potential_over_rt,
    )


def _neutral_associating_mixture() -> dict[str, object]:
    mapping = deepcopy(ASSOCIATING_PURE)
    mapping["document_id"] = "neutral-associating-mixture-fixture"
    components = mapping["components"]
    pairs = mapping["pairs"]
    assert isinstance(components, list)
    assert isinstance(pairs, list)
    components.append(
        _component(
            "inert-fluid",
            segment_count=1.4,
            segment_diameter=3.4,
            dispersion_energy_over_k=180.0,
        )
    )
    pairs.append(_pair("associating-fluid", "inert-fluid", k_ij=0.04))
    return mapping


@pytest.mark.parametrize(
    ("mapping", "temperature_k", "density", "composition", "identities", "steps"),
    (
        (
            neutral_pure("propane"),
            300.0,
            4000.0,
            (1.0,),
            (
                "component/propane/segment_count",
                "component/propane/segment_diameter",
                "component/propane/dispersion_energy_over_k",
            ),
            (1.0e-5, 1.0e-5, 1.0e-3),
        ),
        (
            POLAR_PURE,
            315.0,
            4200.0,
            (1.0,),
            (
                "component/polar-fluid/segment_count",
                "component/polar-fluid/segment_diameter",
                "component/polar-fluid/dispersion_energy_over_k",
            ),
            (1.0e-5, 1.0e-5, 1.0e-3),
        ),
        (
            _neutral_associating_mixture(),
            300.0,
            3500.0,
            (0.6, 0.4),
            (
                "association/associating-fluid/a/associating-fluid/b/energy_over_k",
                "association/associating-fluid/a/associating-fluid/b/volume",
            ),
            (1.0e-3, 1.0e-6),
        ),
        (
            AQUEOUS_ELECTROLYTE,
            298.15,
            2000.0,
            (0.98, 0.01, 0.01),
            (
                "association/water/a/water/b/energy_over_k",
                "association/water/a/water/b/volume",
            ),
            (1.0e-3, 1.0e-6),
        ),
    ),
    ids=(
        "base",
        "polar-fixed-moment",
        "neutral-mixture-association",
        "electrolyte-temperature-dependent-association",
    ),
)
def test_fixed_state_selected_parameter_derivatives_match_independent_differences(
    mapping: dict[str, object],
    temperature_k: float,
    density: float,
    composition: tuple[float, ...],
    identities: tuple[str, ...],
    steps: tuple[float, ...],
) -> None:
    parameters = Parameters.from_mapping(mapping)
    active = ActiveParameterSet(parameters, identities)
    model = Pure(parameters) if len(composition) == 1 else Mixture(parameters)
    exact = (
        model.phase_parameter_derivatives(
            active,
            T=temperature_k * u.kelvin,
            rho=density * u.mole / u.meter**3,
        )
        if isinstance(model, Pure)
        else model.phase_parameter_derivatives(
            active,
            T=temperature_k * u.kelvin,
            rho=density * u.mole / u.meter**3,
            x=composition,
        )
    )

    for column, (identity, step) in enumerate(zip(identities, steps, strict=True)):
        numerical = _central_difference(
            mapping,
            identity,
            step,
            lambda candidate: _fixed_state_values(
                candidate, temperature_k, density, composition
            ),
        )
        for row, expected in enumerate(numerical):
            assert exact.derivatives[row][column] == pytest.approx(
                expected, rel=3.0e-5, abs=2.0e-7
            )


@pytest.mark.parametrize(
    ("identity", "step"),
    (
        ("pair/dipolar/quadrupolar/k_ij", 1.0e-5),
        ("pair/dipolar/quadrupolar/l_ij", 1.0e-5),
    ),
)
def test_binary_interaction_derivatives_match_independent_differences(
    identity: str, step: float
) -> None:
    parameters = Parameters.from_mapping(POLAR_BINARY)
    exact = Mixture(parameters).phase_parameter_derivatives(
        ActiveParameterSet(parameters, (identity,)),
        T=300.0 * u.kelvin,
        rho=5000.0 * u.mole / u.meter**3,
        x=(0.4, 0.6),
    )
    numerical = _central_difference(
        POLAR_BINARY,
        identity,
        step,
        lambda candidate: _fixed_state_values(candidate, 300.0, 5000.0, (0.4, 0.6)),
    )
    for row, expected in enumerate(numerical):
        assert exact.derivatives[row][0] == pytest.approx(
            expected, rel=3.0e-5, abs=2.0e-7
        )


@pytest.mark.parametrize(
    ("mapping", "temperature_k", "identities", "steps"),
    (
        (
            PURE_SATURATION_PARAMETERS["propane"],
            300.0,
            (
                "component/propane/segment_count",
                "component/propane/segment_diameter",
                "component/propane/dispersion_energy_over_k",
            ),
            (1.0e-5, 1.0e-5, 1.0e-3),
        ),
        (
            ASSOCIATING_PURE,
            300.0,
            (
                "association/associating-fluid/a/associating-fluid/b/energy_over_k",
                "association/associating-fluid/a/associating-fluid/b/volume",
            ),
            (1.0e-3, 1.0e-6),
        ),
    ),
    ids=("base", "association"),
)
def test_saturation_implicit_parameter_derivatives_match_resolved_boundaries(
    mapping: dict[str, object],
    temperature_k: float,
    identities: tuple[str, ...],
    steps: tuple[float, ...],
) -> None:
    parameters = Parameters.from_mapping(mapping)
    exact = Pure(parameters).saturation_parameter_derivatives(
        ActiveParameterSet(parameters, identities), T=temperature_k * u.kelvin
    )

    def values(candidate: dict[str, object]) -> Values:
        result = Pure(Parameters.from_mapping(candidate)).saturation(
            T=temperature_k * u.kelvin
        )
        return (
            result.saturation_pressure_pa,
            result.vapor.molar_density_mol_m3,
            result.liquid.molar_density_mol_m3,
        )

    for column, (identity, step) in enumerate(zip(identities, steps, strict=True)):
        numerical = _central_difference(mapping, identity, step, values)
        for row, expected in enumerate(numerical):
            assert exact.derivatives[row][column] == pytest.approx(
                expected, rel=3.0e-5, abs=2.0e-6
            )


def test_fixed_pressure_fugacity_derivatives_include_density_closure() -> None:
    mapping = neutral_pure("propane")
    identities = (
        "component/propane/segment_diameter",
        "component/propane/dispersion_energy_over_k",
    )
    parameters = Parameters.from_mapping(mapping)
    exact = Pure(parameters).fugacity_parameter_derivatives(
        ActiveParameterSet(parameters, identities),
        T=300.0 * u.kelvin,
        P=1.0e5 * u.pascal,
        phase="vapor",
    )

    def value(candidate: dict[str, object]) -> Values:
        state = Pure(Parameters.from_mapping(candidate)).state(
            T=300.0 * u.kelvin, P=1.0e5 * u.pascal, phase="vapor"
        )
        assert state.fugacity is not None
        return state.fugacity.ln_coefficient

    for column, (identity, step) in enumerate(
        zip(identities, (1.0e-5, 1.0e-3), strict=True)
    ):
        numerical = _central_difference(mapping, identity, step, value)
        assert exact.derivatives[0][column] == pytest.approx(
            numerical[0], rel=3.0e-5, abs=2.0e-8
        )


@pytest.mark.parametrize(
    ("identity", "step", "evaluate"),
    (
        (
            "component/sodium-cation/born_diameter",
            1.0e-5,
            lambda model, active: model.ion_solvation_parameter_derivatives(
                active,
                reference=model.reference(
                    T=298.15 * u.kelvin,
                    P=1.0e5 * u.pascal,
                    solvent_x=(1.0, 0.0, 0.0),
                    phase="liquid",
                ),
            ),
        ),
        (
            "model/relative_permittivity/ionic_region_relative_permittivity",
            1.0e-5,
            lambda model, active: model.ion_solvation_parameter_derivatives(
                active,
                reference=model.reference(
                    T=298.15 * u.kelvin,
                    P=1.0e5 * u.pascal,
                    solvent_x=(1.0, 0.0, 0.0),
                    phase="liquid",
                ),
            ),
        ),
        (
            "component/water/relative_permittivity/constant/constant",
            1.0e-4,
            lambda model, active: model.ion_solvation_parameter_derivatives(
                active,
                reference=model.reference(
                    T=298.15 * u.kelvin,
                    P=1.0e5 * u.pascal,
                    solvent_x=(1.0, 0.0, 0.0),
                    phase="liquid",
                ),
            ),
        ),
        (
            "component/water/solvation_factor",
            1.0e-5,
            lambda model, active: model.mean_ionic_activity_parameter_derivatives(
                active,
                T=298.15 * u.kelvin,
                P=1.0e5 * u.pascal,
                molality=1.0 * u.mole / u.kilogram,
            ),
        ),
        (
            "model/relative_permittivity/ion_fraction_suppression_coefficient",
            1.0e-5,
            lambda model, active: model.relative_permittivity_parameter_derivatives(
                active, ion_mole_fraction=0.02
            ),
        ),
    ),
    ids=("born", "ionic-permittivity", "correlated-permittivity", "solvation", "suppression"),
)
def test_electrolyte_auxiliary_derivatives_match_independent_differences(
    identity: str,
    step: float,
    evaluate: Callable[[Mixture, ActiveParameterSet], object],
) -> None:
    def result(mapping: dict[str, object]):
        parameters = Parameters.from_mapping(mapping)
        return evaluate(Mixture(parameters), ActiveParameterSet(parameters, (identity,)))

    exact = result(AQUEOUS_ELECTROLYTE)
    numerical = _central_difference(
        AQUEOUS_ELECTROLYTE,
        identity,
        step,
        lambda candidate: tuple(result(candidate).values),
    )
    for row, expected in enumerate(numerical):
        assert exact.derivatives[row][0] == pytest.approx(
            expected, rel=3.0e-5, abs=2.0e-7
        )


def test_ion_solvation_derivatives_are_component_order_independent() -> None:
    identities = (
        "component/sodium-cation/born_diameter",
        "component/water/solvation_factor",
        "component/water/relative_permittivity/constant/constant",
        "component/water/dispersion_energy_over_k",
    )
    canonical_parameters = Parameters.from_mapping(AQUEOUS_ELECTROLYTE)
    canonical = Mixture(canonical_parameters)
    canonical_result = canonical.ion_solvation_parameter_derivatives(
        ActiveParameterSet(canonical_parameters, identities),
        reference=canonical.reference(
            T=298.15 * u.kelvin,
            P=1.0e5 * u.pascal,
            solvent_x=(1.0, 0.0, 0.0),
            phase="liquid",
        ),
    )

    reordered_parameters = Parameters.from_mapping(
        AQUEOUS_ELECTROLYTE,
        components=("chloride-anion", "water", "sodium-cation"),
    )
    reordered = Mixture(reordered_parameters)
    reordered_result = reordered.ion_solvation_parameter_derivatives(
        ActiveParameterSet(reordered_parameters, identities),
        reference=reordered.reference(
            T=298.15 * u.kelvin,
            P=1.0e5 * u.pascal,
            solvent_x=(0.0, 1.0, 0.0),
            phase="liquid",
            ion_component_ids=("sodium-cation", "chloride-anion"),
        ),
    )

    assert reordered_result.values == pytest.approx(canonical_result.values, rel=2.0e-12)
    assert reordered_result.derivatives[0] == pytest.approx(
        canonical_result.derivatives[0], rel=2.0e-10
    )
    assert reordered_result.derivatives[1] == pytest.approx(
        canonical_result.derivatives[1], rel=2.0e-10
    )
