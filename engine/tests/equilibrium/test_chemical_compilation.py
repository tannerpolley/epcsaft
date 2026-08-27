from __future__ import annotations

import math

import pytest

import epcsaft
from epcsaft import _core
from tests.eos._raw_cases import NEUTRAL_TERNARY_IMMISCIBLE

from .chemical_equilibrium_cases import base_system as _base_system
from .chemical_equilibrium_cases import bind_records as _bind_record
from .chemical_equilibrium_cases import typed_problem as _typed_problem


def _solve(
    spec: dict[str, object],
) -> object:
    model = epcsaft.Mixture(epcsaft.Parameters.from_mapping(NEUTRAL_TERNARY_IMMISCIBLE))
    feed = tuple(float(value) for value in spec["feed_amounts"])
    if "conserved_totals" not in spec:
        spec["conserved_totals"] = tuple(
            math.fsum(coefficient * amount for coefficient, amount in zip(row, feed, strict=True))
            for row in spec["balance_matrix"]
        )
    phase = epcsaft.equilibrium.ReactivePhase(
        "ideal",
        "vapor",
        "finite",
        epcsaft.equilibrium.AllComponents(),
        epcsaft.equilibrium.IdealGasModel("ideal-gas"),
        epcsaft.equilibrium.FinitePhaseStart(feed, 1.0),
    )
    problem = epcsaft.equilibrium.GeneralReactiveEquilibriumProblem(
        "chemical-compilation",
        epcsaft.equilibrium.Fixed(float(spec["temperature_k"]) * epcsaft.unit_registry.kelvin),
        epcsaft.equilibrium.Fixed(float(spec["pressure_pa"]) * epcsaft.unit_registry.pascal),
        (phase,),
        _typed_problem(spec, minimum_amount_mol=1.0e-30),
        (phase.identity,),
        tuple(
            epcsaft.equilibrium.EquilibriumOutput(
                f"amount-{species}",
                "phase.amount",
                "mole",
                "true-species-amount",
                phase.identity,
                tuple(float(index == species) for index in range(len(feed))),
            )
            for species in range(len(feed))
        ),
    )
    return _core._solve_grepe(model._native, model._native, problem, None, "").equilibrium


@pytest.mark.parametrize("inconsistent", (False, True))
def test_redundant_reactions_preserve_equilibrium_or_reject_bad_cycles(
    inconsistent: bool,
) -> None:
    spec = {
        **_base_system(),
        "species_ids": ("A", "B", "C"),
        "charges": (0, 0, 0),
        "molar_masses_kg_per_mol": (1.0, 1.0, 1.0),
        "balance_matrix": ((1.0, 1.0, 1.0),),
        "reaction_matrix": (
            (-1.0, 1.0, 0.0),
            (0.0, -1.0, 1.0),
            (-1.0, 0.0, 1.0),
        ),
        "feed_amounts": (1.0, 0.0, 0.0),
        "ln_k": (
            math.log(2.0),
            math.log(3.0),
            math.log(7.0) if inconsistent else math.log(6.0),
        ),
    }
    _bind_record(spec)

    if inconsistent:
        result = _solve(spec)
        assert not result.accepted
        assert "inconsistent" in result.failure_reason
        return

    result = _solve(spec)
    assert result.accepted is True
    assert result.output_values == pytest.approx((1.0 / 9.0, 2.0 / 9.0, 2.0 / 3.0))
    assert result.reaction_affinity_inf_norm <= 1.0e-7


def _nearly_dependent_system(
    reaction_scale: float, delta: float, *, reverse_order: bool = False
) -> dict[str, object]:
    reactions: tuple[tuple[float, ...], ...] = (
        (-1.0, 1.0, 0.0, 0.0),
        (
            -reaction_scale,
            reaction_scale,
            reaction_scale * delta,
            -reaction_scale * delta,
        ),
    )
    ln_k: tuple[float, ...] = (2.0, reaction_scale * (2.0 - delta))
    if reverse_order:
        reactions = tuple(reversed(reactions))
        ln_k = tuple(reversed(ln_k))
    spec = {
        **_base_system(),
        "species_ids": ("A", "B", "C", "D"),
        "charges": (0, 0, 0, 0),
        "molar_masses_kg_per_mol": (1.0, 1.0, 1.0, 1.0),
        "balance_matrix": (
            (1.0, 1.0, 1.0, 1.0),
            (1.0, 1.0, 0.0, 0.0),
        ),
        "reaction_matrix": reactions,
        "feed_amounts": (1.0, 1.0, 1.0, 1.0),
        "ln_k": ln_k,
    }
    _bind_record(spec)
    return spec


def test_scaled_nearly_dependent_reactions_are_stable_within_condition_limit() -> None:
    expected = (
        2.0 / (1.0 + math.exp(2.0)),
        2.0 * math.exp(2.0) / (1.0 + math.exp(2.0)),
        2.0 / (1.0 + math.e),
        2.0 * math.e / (1.0 + math.e),
    )
    for scale, reverse_order in (
        (1.0, False),
        (1.0, True),
        (1.0e8, False),
        (1.0e8, True),
    ):
        result = _solve(_nearly_dependent_system(scale, 1.5e-6, reverse_order=reverse_order))
        assert result.output_values == pytest.approx(expected, rel=2.0e-8)

    with pytest.raises(ValueError):
        _solve(_nearly_dependent_system(1.0, 1.0e-14))


def test_reaction_system_rejects_charge_inconsistency() -> None:
    spec = _base_system()
    spec.update({"charges": (1, -1), "feed_amounts": (1.0, 0.0)})
    with pytest.raises(ValueError):
        _solve(spec)

    spec["feed_amounts"] = (1.0, 1.0)
    with pytest.raises(ValueError):
        _solve(spec)


def test_accessible_face_preserves_reaction_combinations_and_exact_zeros() -> None:
    spec = {
        **_base_system(),
        "species_ids": ("A", "B", "C", "X"),
        "charges": (0, 0, 0, 0),
        "molar_masses_kg_per_mol": (1.0, 2.0, 1.0, 1.0),
        "balance_matrix": (
            (1.0, 1.0, 1.0, 0.0),
            (0.0, 1.0, 0.0, 1.0),
        ),
        "reaction_matrix": (
            (-1.0, 1.0, 0.0, -1.0),
            (0.0, -1.0, 1.0, 1.0),
        ),
        "feed_amounts": (1.0, 0.0, 0.0, 0.0),
        "ln_k": (math.log(2.0), math.log(3.0)),
    }
    _bind_record(spec)

    result = _solve(spec)

    assert result.output_values == pytest.approx(
        (1.0 / 7.0, 0.0, 6.0 / 7.0, 0.0), rel=3.0e-8, abs=0.0
    )
    assert result.phase_amounts_mol[0] == pytest.approx(
        (1.0 / 7.0, 0.0, 6.0 / 7.0, 0.0), rel=3.0e-8, abs=0.0
    )
