from __future__ import annotations

import math

import pint
import pytest

from epcsaft import equilibrium


def test_feed_bases_resolve_to_one_amount_vector_and_fail_closed() -> None:
    units = pint.UnitRegistry()
    component_ids = ("water", "cosolvent", "salt")
    molar_masses = {
        "water": 0.01801528 * units.kilogram / units.mole,
        "cosolvent": 0.04606844 * units.kilogram / units.mole,
        "salt": 0.05844 * units.kilogram / units.mole,
    }
    solvent = equilibrium.solvent_composition_molalities(
        component_ids,
        {"water": 0.7, "cosolvent": 0.3},
        1.0 * units.kilogram,
        {"salt": 1.5 * units.mole / units.kilogram},
        {key: value for key, value in molar_masses.items() if key != "salt"},
    )
    total_moles = math.fsum(solvent.amounts_mol)
    total_mass = math.fsum(
        amount * float(molar_masses[identity].to("kilogram / mole").magnitude)
        for amount, identity in zip(solvent.amounts_mol, component_ids, strict=True)
    )
    forms = (
        equilibrium.direct_amounts(
            component_ids, tuple(value * units.mole for value in solvent.amounts_mol)
        ),
        equilibrium.mole_fractions(
            component_ids,
            tuple(value / total_moles for value in solvent.amounts_mol),
            total_moles * units.mole,
        ),
        equilibrium.mass_fractions(
            component_ids,
            tuple(
                amount * float(molar_masses[identity].to("kilogram / mole").magnitude) / total_mass
                for amount, identity in zip(solvent.amounts_mol, component_ids, strict=True)
            ),
            total_mass * units.kilogram,
            molar_masses,
        ),
        solvent,
    )
    for candidate in forms:
        assert candidate.amounts_mol == pytest.approx(solvent.amounts_mol, rel=0.0, abs=2.0e-14)
    with pytest.raises(ValueError):
        equilibrium.mole_fractions(component_ids, (0.5, 0.5, 0.1), 1.0)
    with pytest.raises(ValueError):
        equilibrium.direct_amounts(component_ids, (1.0 * units.kelvin, 1.0, 1.0))


def test_affine_feed_preserves_exact_jacobian_bounds_and_permutation() -> None:
    component_ids = ("water", "cosolvent", "salt")
    feed = equilibrium.AffineFeed(
        component_ids,
        (0.7, 0.3, 0.0),
        ((0.0, 0.2), (0.0, 0.1), (1.0, -0.3)),
        (
            equilibrium.FeedCoordinate("coordinate-a", (0.0, 1.0)),
            equilibrium.FeedCoordinate("coordinate-b", (0.0, 1.0)),
        ),
    )
    resolved = feed.resolve({"coordinate-b": 0.5, "coordinate-a": 0.25})
    assert resolved.amounts_mol == pytest.approx((0.8, 0.35, 0.1), abs=1.0e-14)
    assert resolved.dn_dq == feed.basis_matrix_mol_per_coordinate
    permuted = feed.reorder(("salt", "water", "cosolvent")).resolve((0.25, 0.5))
    assert permuted.amounts_mol == pytest.approx((0.1, 0.8, 0.35), abs=1.0e-14)
    assert permuted.dn_dq == ((1.0, -0.3), (0.0, 0.2), (0.0, 0.1))
    with pytest.raises(ValueError):
        feed.resolve((1.1, 0.0))


def test_output_expressions_preserve_chain_rule_and_reject_invalid_algebra() -> None:
    rows = (
        equilibrium.EvaluatedRow("x", 2.0, (1.0, 0.5)),
        equilibrium.EvaluatedRow("y", 3.0, (-2.0, 4.0)),
    )
    result = equilibrium.evaluate_expression(
        equilibrium.LinearSum(
            (
                (2.0, equilibrium.RowReference("x")),
                (
                    -3.0,
                    equilibrium.Product(
                        equilibrium.RowReference("x"), equilibrium.RowReference("y")
                    ),
                ),
            ),
            constant=4.0,
        ),
        rows,
    )
    assert result.value == pytest.approx(-10.0)
    assert result.jacobian == pytest.approx((5.0, -27.5))

    invalid = (
        (
            equilibrium.Ratio(
                equilibrium.RowReference("x"),
                equilibrium.RowReference("small"),
                denominator_floor=1.0e-8,
            ),
            rows + (equilibrium.EvaluatedRow("small", 1.0e-9, (0.0, 0.0)),),
        ),
        (
            equilibrium.RowReference("missing"),
            rows,
        ),
        (
            equilibrium.RowReference("x"),
            (equilibrium.EvaluatedRow("x", float("nan"), (1.0,)),),
        ),
    )
    for expression, candidate_rows in invalid:
        with pytest.raises(equilibrium.ExpressionEvaluationError):
            equilibrium.evaluate_expression(expression, candidate_rows)
