from __future__ import annotations

import math

import pint
import pytest

from epcsaft.equilibrium import (
    AffineFeed,
    FeedCoordinate,
    direct_amounts,
    mass_fractions,
    mole_fractions,
    solvent_composition_molalities,
)

_UNIT_REGISTRY = pint.UnitRegistry()
_COMPONENTS = ("water", "cosolvent", "salt")
_MOLAR_MASSES = {
    "water": 0.01801528 * _UNIT_REGISTRY.kilogram / _UNIT_REGISTRY.mole,
    "cosolvent": 0.04606844 * _UNIT_REGISTRY.kilogram / _UNIT_REGISTRY.mole,
    "salt": 0.05844 * _UNIT_REGISTRY.kilogram / _UNIT_REGISTRY.mole,
}


def test_feed_basis_forms_resolve_to_one_canonical_amount_vector() -> None:
    """Different source bases must describe the same conserved component amounts."""

    solvent = solvent_composition_molalities(
        _COMPONENTS,
        {"water": 0.7, "cosolvent": 0.3},
        1.0 * _UNIT_REGISTRY.kilogram,
        {"salt": 1.5 * _UNIT_REGISTRY.mole / _UNIT_REGISTRY.kilogram},
        {key: value for key, value in _MOLAR_MASSES.items() if key != "salt"},
    )
    total_moles = math.fsum(solvent.amounts_mol)
    fractions = tuple(value / total_moles for value in solvent.amounts_mol)
    total_mass = math.fsum(
        amount * float(_MOLAR_MASSES[identity].to("kilogram / mole").magnitude)
        for amount, identity in zip(solvent.amounts_mol, _COMPONENTS, strict=True)
    )
    mass_fraction_values = tuple(
        amount * float(_MOLAR_MASSES[identity].to("kilogram / mole").magnitude) / total_mass
        for amount, identity in zip(solvent.amounts_mol, _COMPONENTS, strict=True)
    )
    forms = (
        direct_amounts(
            _COMPONENTS,
            tuple(value * _UNIT_REGISTRY.mole for value in solvent.amounts_mol),
        ),
        mole_fractions(
            _COMPONENTS,
            fractions,
            total_moles * _UNIT_REGISTRY.mole,
        ),
        mass_fractions(
            _COMPONENTS,
            mass_fraction_values,
            total_mass * _UNIT_REGISTRY.kilogram,
            _MOLAR_MASSES,
        ),
        solvent,
    )
    for candidate in forms:
        assert candidate.amounts_mol == pytest.approx(solvent.amounts_mol, rel=0.0, abs=2.0e-14)
    assert solvent.reorder(("salt", "water", "cosolvent")).amounts_mol == pytest.approx(
        (solvent.amounts_mol[2], solvent.amounts_mol[0], solvent.amounts_mol[1]),
        rel=0.0,
        abs=2.0e-14,
    )
    with pytest.raises(ValueError, match="normalized"):
        mole_fractions(_COMPONENTS, (0.5, 0.5, 0.1), 1.0)
    with pytest.raises(ValueError, match="convertible"):
        direct_amounts(_COMPONENTS, (1.0 * _UNIT_REGISTRY.kelvin, 1.0, 1.0))


def test_affine_feed_enforces_bounds_and_preserves_exact_jacobian_under_permutation() -> None:
    """Affine conserved coordinates expose B exactly and remain identity-based."""

    feed = AffineFeed(
        _COMPONENTS,
        (0.7, 0.3, 0.0),
        ((0.0, 0.2), (0.0, 0.1), (1.0, -0.3)),
        (
            FeedCoordinate("feed-coordinate-a", (0.0, 1.0)),
            FeedCoordinate("feed-coordinate-b", (0.0, 1.0)),
        ),
    )
    resolved = feed.resolve({"feed-coordinate-b": 0.5, "feed-coordinate-a": 0.25})
    assert resolved.amounts_mol == pytest.approx((0.8, 0.35, 0.1), rel=0.0, abs=1.0e-14)
    assert resolved.dn_dq == feed.basis_matrix_mol_per_coordinate
    permuted = feed.reorder(("salt", "water", "cosolvent")).resolve((0.25, 0.5))
    assert permuted.amounts_mol == pytest.approx((0.1, 0.8, 0.35), rel=0.0, abs=1.0e-14)
    assert permuted.dn_dq == (
        (1.0, -0.3),
        (0.0, 0.2),
        (0.0, 0.1),
    )
    with pytest.raises(ValueError, match="outside"):
        feed.resolve((1.1, 0.0))
