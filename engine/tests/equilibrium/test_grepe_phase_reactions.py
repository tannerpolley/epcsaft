from __future__ import annotations

import math

import pytest

from epcsaft.equilibrium import (
    PhaseReactionRow,
    PhaseReactionStandardState,
    PhaseReactionSystem,
    PhaseSpeciesTerm,
    canonical_stacked_directions,
    ideal_exchange_dxi_dln_k,
    ideal_exchange_extent,
    ideal_exchange_ln_q,
    preflight_phase_reaction_system,
)


def _exchange_system() -> PhaseReactionSystem:
    return PhaseReactionSystem(
        phase_ids=("aqueous", "organic"),
        species_ids=("Li+", "Na+", "Cl-", "NaBTA·TOPO", "LiBTA·TOPO"),
        reactions=(
            PhaseReactionRow(
                "li-na-carrier-exchange",
                (
                    PhaseSpeciesTerm("aqueous", "Li+", -1.0),
                    PhaseSpeciesTerm("aqueous", "Na+", 1.0),
                    PhaseSpeciesTerm("organic", "NaBTA·TOPO", -1.0),
                    PhaseSpeciesTerm("organic", "LiBTA·TOPO", 1.0),
                ),
                ln_equilibrium_constant=0.0,
                standard_state=PhaseReactionStandardState("dimensionless-activity"),
            ),
        ),
        charges=(1.0, 1.0, -1.0, 0.0, 0.0),
        balance_matrix=(
            (1.0, 0.0, 0.0, 0.0, 1.0),
            (0.0, 1.0, 0.0, 1.0, 0.0),
            (0.0, 0.0, 0.0, 1.0, 1.0),
        ),
    )


def test_phase_local_row_is_a_special_case_of_stacked_direction() -> None:
    system = PhaseReactionSystem(
        phase_ids=("liquid", "vapor"),
        species_ids=("A", "B", "C"),
        reactions=(
            PhaseReactionRow(
                "liquid-local",
                (
                    PhaseSpeciesTerm("liquid", "A", -1.0),
                    PhaseSpeciesTerm("liquid", "B", 1.0),
                ),
            ),
        ),
        charges=(0.0, 0.0, 0.0),
        balance_matrix=((1.0, 1.0, 0.0), (0.0, 0.0, 1.0)),
    )
    preflight = preflight_phase_reaction_system(system)
    assert preflight.rank == 1
    assert preflight.coordinate_order == (
        ("liquid", "A"),
        ("liquid", "B"),
        ("liquid", "C"),
        ("vapor", "A"),
        ("vapor", "B"),
        ("vapor", "C"),
    )
    assert preflight.direction_matrix == ((-1.0, 1.0, 0.0, 0.0, 0.0, 0.0),)


def test_cross_phase_li_na_exchange_conserves_global_material_and_phase_charge() -> None:
    system = _exchange_system()
    preflight = preflight_phase_reaction_system(system)
    assert preflight.rank == 1
    assert canonical_stacked_directions(system) == (
        (-1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 1.0),
    )


def test_ideal_exchange_extent_and_exact_log_k_sensitivity_are_independent_oracle() -> None:
    initial = (0.8, 0.2, 0.7, 0.3)
    target_extent = 0.1
    ln_k = ideal_exchange_ln_q(initial, target_extent)
    assert ideal_exchange_extent(initial, ln_k) == pytest.approx(target_extent, abs=1.0e-13)
    expected = 1.0 / (1.0 / 0.3 + 1.0 / 0.7 + 1.0 / 0.4 + 1.0 / 0.6)
    assert ideal_exchange_dxi_dln_k(initial, target_extent) == pytest.approx(expected)
    step = 1.0e-7
    plus = ideal_exchange_extent(initial, ln_k + step)
    minus = ideal_exchange_extent(initial, ln_k - step)
    assert (plus - minus) / (2.0 * step) == pytest.approx(expected, rel=1.0e-7)


@pytest.mark.parametrize(
    ("mutator", "message"),
    [
        (
            lambda system: PhaseReactionSystem(
                system.phase_ids,
                system.species_ids,
                (
                    PhaseReactionRow(
                        "li-na-carrier-exchange",
                        system.reactions[0].terms + (PhaseSpeciesTerm("aqueous", "Li+", 0.5),),
                    ),
                ),
                system.charges,
                system.balance_matrix,
            ),
            "repeats coordinate",
        ),
        (
            lambda system: PhaseReactionSystem(
                system.phase_ids,
                system.species_ids,
                (
                    PhaseReactionRow(
                        "unknown-phase",
                        system.reactions[0].terms[:1] + (PhaseSpeciesTerm("missing", "Na+", 1.0),),
                    ),
                ),
                system.charges,
                system.balance_matrix,
            ),
            "unknown phase",
        ),
        (
            lambda system: PhaseReactionSystem(
                system.phase_ids,
                system.species_ids,
                (
                    PhaseReactionRow(
                        "nonfinite",
                        (PhaseSpeciesTerm("aqueous", "Li+", math.inf),),
                    ),
                ),
                system.charges,
                system.balance_matrix,
            ),
            "must be finite",
        ),
    ],
)
def test_preflight_rejects_invalid_sparse_rows(mutator, message: str) -> None:
    with pytest.raises((TypeError, ValueError), match=message):
        preflight_phase_reaction_system(mutator(_exchange_system()))


def test_preflight_rejects_dependent_rows_and_conservation_or_charge_violation() -> None:
    base = _exchange_system()
    dependent = PhaseReactionSystem(
        base.phase_ids,
        base.species_ids,
        base.reactions
        + (
            PhaseReactionRow(
                "duplicate-direction",
                tuple(
                    PhaseSpeciesTerm(term.phase_id, term.species_id, 2.0 * term.coefficient)
                    for term in base.reactions[0].terms
                ),
            ),
        ),
        base.charges,
        base.balance_matrix,
    )
    with pytest.raises(ValueError, match="independent"):
        preflight_phase_reaction_system(dependent)

    bad_balance = PhaseReactionSystem(
        base.phase_ids,
        base.species_ids,
        base.reactions,
        base.charges,
        base.balance_matrix[:-1] + ((1.0, 0.0, 0.0, 0.0, 0.0),),
    )
    with pytest.raises(ValueError, match="global balance"):
        preflight_phase_reaction_system(bad_balance)

    bad_charge = PhaseReactionSystem(
        base.phase_ids,
        base.species_ids,
        base.reactions,
        (1.0, 0.0, -1.0, 0.0, 0.0),
        base.balance_matrix,
    )
    with pytest.raises(ValueError, match="charge"):
        preflight_phase_reaction_system(bad_charge)
