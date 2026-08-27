from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from typing import Literal, cast

import pytest

from epcsaft import Mixture, NonEvaluableTrial, Parameters, State
from epcsaft import unit_registry as u
from tests.eos._raw_cases import (
    AQUEOUS_ELECTROLYTE,
    NEUTRAL_BINARY,
    STRONGLY_ASSOCIATING_BINARY,
)

TOLERANCE = 3.0e-6
FLOOR = 1.0e-8
type Phase = Literal["vapor", "liquid"]


def _state(
    mapping: Mapping[str, object],
    temperature_k: float,
    pressure_pa: float,
    composition: Sequence[float],
    phase: Phase,
    components: Sequence[str] | None = None,
) -> State:
    parameters = (
        Parameters.from_mapping(mapping)
        if components is None
        else Parameters.from_mapping(mapping, components=components)
    )
    return Mixture(parameters).state(
        T=temperature_k * u.kelvin,
        P=pressure_pa * u.pascal,
        x=composition,
        phase=phase,
    )


def _available(state: State) -> Mapping[str, object]:
    block = state.fixed_pressure_composition_derivatives
    assert block is not None and block["status"] == "available", block
    assert block["failure"] is None
    return block


def _matrix(block: Mapping[str, object], key: str) -> tuple[tuple[float, ...], ...]:
    return cast(tuple[tuple[float, ...], ...], block[key])


def _perturb(x: Sequence[float], basis: Sequence[float], step: float) -> tuple[float, ...]:
    values = tuple(xi * math.exp(step * bi) for xi, bi in zip(x, basis, strict=True))
    total = math.fsum(values)
    return tuple(value / total for value in values)


def _log_activities(state: State) -> tuple[float, ...]:
    assert state.fugacity is not None
    return tuple(
        math.log(xi) + ln_phi
        for xi, ln_phi in zip(state.mole_fractions, state.fugacity.ln_coefficient, strict=True)
    )


def _disagreement(actual: float, expected: float) -> float:
    return abs(actual - expected) / max(abs(actual), abs(expected), FLOOR)


def _centered_check(
    mapping: Mapping[str, object],
    temperature_k: float,
    pressure_pa: float,
    composition: tuple[float, ...],
    phase: Phase,
    step: float,
) -> tuple[float, float, float]:
    center = _state(mapping, temperature_k, pressure_pa, composition, phase)
    block = _available(center)
    basis = tuple(row[0] for row in _matrix(block, "log_composition_basis"))
    minus, plus = (
        _state(
            mapping,
            temperature_k,
            pressure_pa,
            _perturb(composition, basis, direction * step),
            phase,
        )
        for direction in (-1.0, 1.0)
    )
    finite_difference = tuple(
        (right - left) / (2.0 * step)
        for left, right in zip(_log_activities(minus), _log_activities(plus), strict=True)
    )
    gamma = _matrix(block, "chemical_potential_derivatives_over_rt")
    density = cast(tuple[float, ...], block["molar_density_derivatives_mol_per_m3"])
    density_difference = (
        plus.molar_density.to("mole / meter ** 3").magnitude
        - minus.molar_density.to("mole / meter ** 3").magnitude
    ) / (2.0 * step)
    return (
        max(
            _disagreement(row[0], expected)
            for row, expected in zip(gamma, finite_difference, strict=True)
        ),
        _disagreement(density[0], density_difference),
        max(
            abs(state.pressure.to("pascal").magnitude - pressure_pa) / pressure_pa
            for state in (minus, plus)
        ),
    )


def test_neutral_ideal_limit_centered_steps_gibbs_duhem_and_density_closure() -> None:
    ideal_x = (0.4, 0.6)
    ideal_block = _available(_state(NEUTRAL_BINARY, 300.0, 1.0, ideal_x, "vapor"))
    ideal_basis = _matrix(ideal_block, "log_composition_basis")
    ideal_gamma = _matrix(ideal_block, "chemical_potential_derivatives_over_rt")
    assert (
        max(
            _disagreement(actual, expected)
            for actual_row, expected_row in zip(ideal_gamma, ideal_basis, strict=True)
            for actual, expected in zip(actual_row, expected_row, strict=True)
        )
        <= TOLERANCE
    )
    liquid = _state(
        STRONGLY_ASSOCIATING_BINARY,
        376.88,
        101330.0,
        (0.1385, 0.8615),
        "liquid",
    )
    liquid_block = _available(liquid)
    liquid_gamma = _matrix(liquid_block, "chemical_potential_derivatives_over_rt")
    assert (
        abs(
            math.fsum(
                xi * row[0] for xi, row in zip(liquid.mole_fractions, liquid_gamma, strict=True)
            )
        )
        <= 1.0e-12
    )
    assert 1.0 <= cast(float, liquid_block["condition_measure"]) <= 1.0e12
    for phase, temperature_k, pressure_pa, composition in (
        ("liquid", 376.88, 101330.0, (0.1385, 0.8615)),
        ("vapor", 350.0, 500000.0, (0.3, 0.7)),
    ):
        for step in (1.0e-4, 5.0e-5):
            gamma_error, density_error, pressure_error = _centered_check(
                STRONGLY_ASSOCIATING_BINARY,
                temperature_k,
                pressure_pa,
                composition,
                phase,
                step,
            )
            assert gamma_error <= TOLERANCE
            assert density_error <= TOLERANCE
            assert pressure_error <= 1.0e-11


def test_charged_tangent_and_component_permutation_are_invariant() -> None:
    composition = (0.98, 0.01, 0.01)
    state = _state(AQUEOUS_ELECTROLYTE, 298.15, 100000.0, composition, "liquid")
    block = _available(state)
    assert block["constraint_identities"] == ("normalization", "electroneutrality")
    basis = tuple(row[0] for row in _matrix(block, "log_composition_basis"))
    perturbed = tuple(_perturb(composition, basis, sign * 2.0e-4) for sign in (-1.0, 1.0))
    assert all(abs(math.fsum(values) - 1.0) <= 2.0e-16 for values in perturbed)
    assert all(abs(values[1] - values[2]) <= 2.0e-16 for values in perturbed)
    assert (
        _centered_check(AQUEOUS_ELECTROLYTE, 298.15, 100000.0, composition, "liquid", 2.0e-4)[0]
        <= TOLERANCE
    )

    reordered = _state(
        AQUEOUS_ELECTROLYTE,
        298.15,
        100000.0,
        (0.01, 0.98, 0.01),
        "liquid",
        ("chloride-anion", "water", "sodium-cation"),
    )
    original_by_id = dict(
        zip(
            state.component_ids,
            (row[0] for row in _matrix(block, "chemical_potential_derivatives_over_rt")),
            strict=True,
        )
    )
    reordered_by_id = dict(
        zip(
            reordered.component_ids,
            (
                row[0]
                for row in _matrix(_available(reordered), "chemical_potential_derivatives_over_rt")
            ),
            strict=True,
        )
    )
    assert reordered_by_id == pytest.approx(original_by_id, rel=2.0e-12, abs=1.0e-13)
    assert (
        abs(
            math.fsum(
                xi * row[0]
                for xi, row in zip(
                    state.mole_fractions,
                    _matrix(block, "chemical_potential_derivatives_over_rt"),
                    strict=True,
                )
            )
        )
        <= 1.0e-12
    )


def test_boundary_returns_typed_failure_without_matrices() -> None:
    block = _state(
        NEUTRAL_BINARY, 300.0, 100000.0, (0.0, 1.0), "liquid"
    ).fixed_pressure_composition_derivatives
    assert block is not None and block["status"] == "non_evaluable"
    assert isinstance(block["failure"], NonEvaluableTrial)
    assert block["failure"].code == "provider_domain_rejection"
    assert "log_composition_basis" not in block
    assert "chemical_potential_derivatives_over_rt" not in block
