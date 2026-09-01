from __future__ import annotations

import math
from collections.abc import Mapping
from dataclasses import replace
from typing import cast

import pytest

from epcsaft import EosError, Mixture, Parameters, State
from epcsaft import unit_registry as u
from tests.eos._cases import AQUEOUS_ELECTROLYTE, STRONGLY_ASSOCIATING_BINARY


def _state(model: Mixture, composition: tuple[float, ...], anchor: State | None = None) -> State:
    return model.state(
        T=376.88 * u.kelvin,
        P=101330.0 * u.pascal,
        x=composition,
        phase="liquid",
        density_anchor=anchor,
    )


def _targets() -> tuple[tuple[float, float], ...]:
    return tuple(
        (
            scaled / (scaled + 0.8615),
            0.8615 / (scaled + 0.8615),
        )
        for step in range(21)
        for scaled in (0.1385 * math.exp(-0.05 + 0.1 * step / 20.0),)
    )


def _derivatives(state: State) -> Mapping[str, object]:
    block = state.fixed_pressure_composition_derivatives
    assert block is not None and block["status"] == "available"
    return block


@pytest.mark.parametrize("targets", (_targets(), tuple(reversed(_targets()))))
def test_density_anchor_preserves_local_branch_values_and_exact_derivatives(
    targets: tuple[tuple[float, float], ...],
) -> None:
    model = Mixture(Parameters.from_mapping(STRONGLY_ASSOCIATING_BINARY))
    ordinary = tuple(_state(model, composition) for composition in targets)
    anchored = [ordinary[0]]
    for composition in targets[1:]:
        anchored.append(_state(model, composition, anchored[-1]))

    for cold, warm in zip(ordinary, anchored, strict=True):
        assert warm.certified_branch == anchored[0].certified_branch
        assert warm.requested_phase == cold.requested_phase == "liquid"
        assert warm.density_diagnostics is not None
        assert cold.density_diagnostics is not None
        assert warm.certified_branch["stable_root_index"] == (
            cold.density_diagnostics.selected_stable_root_index
        )
        assert warm.density_diagnostics.anchor_status in {"not_supplied", "accepted"}
        assert warm.pressure.to("pascal").magnitude == pytest.approx(
            cold.pressure.to("pascal").magnitude, rel=1.0e-10, abs=1.0e-12
        )
        assert warm.molar_density.to("mole / meter ** 3").magnitude == pytest.approx(
            cold.molar_density.to("mole / meter ** 3").magnitude, rel=1.0e-10, abs=1.0e-12
        )
        assert warm.residual_chemical_potential_over_rt == pytest.approx(
            cold.residual_chemical_potential_over_rt, rel=1.0e-10, abs=1.0e-12
        )
        warm_derivatives = cast(
            tuple[tuple[float, ...], ...],
            _derivatives(warm)["chemical_potential_derivatives_over_rt"],
        )
        cold_derivatives = cast(
            tuple[tuple[float, ...], ...],
            _derivatives(cold)["chemical_potential_derivatives_over_rt"],
        )
        assert tuple(value for row in warm_derivatives for value in row) == pytest.approx(
            tuple(value for row in cold_derivatives for value in row), rel=3.0e-6, abs=1.0e-8
        )


def test_density_anchor_rejects_mismatched_or_unreachable_branches() -> None:
    model = Mixture(Parameters.from_mapping(STRONGLY_ASSOCIATING_BINARY))
    anchor = _state(model, (0.1385, 0.8615))

    for request in (
        {
            "T": 376.88 * u.kelvin,
            "rho": anchor.molar_density,
            "x": (0.1385, 0.8615),
            "density_anchor": anchor,
        },
        {
            "T": 376.88 * u.kelvin,
            "P": 101330.0 * u.pascal,
            "x": (0.1385, 0.8615),
            "density_anchor": anchor,
        },
    ):
        with pytest.raises(EosError) as captured:
            model.state(**request)
        assert captured.value.code == "invalid-density-anchor"

    for density in (math.nan, 0.0):
        malformed = replace(
            anchor,
            molar_density=density * u.mole / u.meter**3,
        )
        with pytest.raises(EosError) as captured:
            _state(model, (0.1385, 0.8615), malformed)
        assert captured.value.code == "invalid-density-anchor"

    for code, kwargs in (
        ("invalid-density-anchor", {"density_anchor": object()}),
        ("density-anchor-phase-mismatch", {"phase": "vapor", "density_anchor": anchor}),
        ("density-anchor-state-mismatch", {"T": 376.89 * u.kelvin, "density_anchor": anchor}),
    ):
        request = {
            "T": 376.88 * u.kelvin,
            "P": 101330.0 * u.pascal,
            "x": (0.1385, 0.8615),
            "phase": "liquid",
            **kwargs,
        }
        with pytest.raises(EosError) as captured:
            model.state(**request)
        assert captured.value.code == code

    electrolyte = Mixture(Parameters.from_mapping(AQUEOUS_ELECTROLYTE)).state(
        T=298.15 * u.kelvin,
        P=100000.0 * u.pascal,
        x=(0.98, 0.01, 0.01),
        phase="liquid",
    )
    with pytest.raises(EosError) as model_mismatch:
        _state(model, (0.1385, 0.8615), electrolyte)
    assert model_mismatch.value.code == "density-anchor-model-mismatch"

    with pytest.raises(EosError) as unreachable:
        _state(model, (0.99, 0.01), anchor)
    assert unreachable.value.code == "density-anchor-branch-unreachable"
