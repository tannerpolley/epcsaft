"""Source-neutral conserved-feed normalization for GREPE.

The functions in this module convert common experimental composition bases to
one canonical amount vector.  They deliberately do not know about a species,
paper, phase, or solver formulation.  ``AffineFeed`` is the corresponding
generic conserved-feed parameterization, ``n = base + Bq``.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from numbers import Real

from pint import Quantity


def _component_ids(value: Sequence[str], *, field: str = "component_ids") -> tuple[str, ...]:
    if isinstance(value, (str, bytes)):
        raise TypeError(f"{field} must be a sequence of names")
    result = tuple(value)
    if (
        not result
        or any(not isinstance(item, str) or not item or item != item.strip() for item in result)
        or len(set(result)) != len(result)
    ):
        raise ValueError(f"{field} must contain unique nonempty names")
    return result


def _scalar(value: object, unit: str, *, field: str, nonnegative: bool = False) -> float:
    if isinstance(value, Quantity):
        try:
            result = float(value.to(unit).magnitude)
        except (AttributeError, TypeError, ValueError) as error:
            raise ValueError(f"{field} must be convertible to {unit}") from error
    elif isinstance(value, Real) and not isinstance(value, bool):
        result = float(value)
    else:
        raise TypeError(f"{field} must be a number or quantity convertible to {unit}")
    if not math.isfinite(result):
        raise ValueError(f"{field} must be finite")
    if nonnegative and result < 0.0:
        raise ValueError(f"{field} must be nonnegative")
    return result


def _vector(
    value: Sequence[object] | Mapping[str, object],
    component_ids: tuple[str, ...],
    unit: str,
    *,
    field: str,
    nonnegative: bool = False,
) -> tuple[float, ...]:
    if isinstance(value, Mapping):
        if set(value) != set(component_ids):
            raise ValueError(f"{field} mapping must name every component exactly once")
        raw = tuple(value[identity] for identity in component_ids)
    elif isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
        if len(value) != len(component_ids):
            raise ValueError(f"{field} must match component_ids length")
        raw = tuple(value)
    else:
        raise TypeError(f"{field} must be a sequence or component mapping")
    return tuple(
        _scalar(item, unit, field=f"{field}[{index}]", nonnegative=nonnegative)
        for index, item in enumerate(raw)
    )


def _fractions(
    value: Sequence[object] | Mapping[str, object],
    component_ids: tuple[str, ...],
    *,
    field: str,
) -> tuple[float, ...]:
    result = _vector(value, component_ids, "dimensionless", field=field, nonnegative=True)
    if not math.isclose(math.fsum(result), 1.0, rel_tol=0.0, abs_tol=1.0e-12):
        raise ValueError(f"{field} must be normalized to one")
    return result


@dataclass(frozen=True, slots=True)
class ResolvedFeed:
    """Canonical component amounts and an optional exact ``dn/dq`` matrix."""

    component_ids: tuple[str, ...]
    amounts_mol: tuple[float, ...]
    coordinate_ids: tuple[str, ...] = ()
    amount_jacobian_mol_per_coordinate: tuple[tuple[float, ...], ...] = ()

    def __post_init__(self) -> None:
        component_ids = _component_ids(self.component_ids)
        amounts = tuple(
            _scalar(value, "mole", field=f"amounts_mol[{index}]", nonnegative=True)
            for index, value in enumerate(self.amounts_mol)
        )
        if len(amounts) != len(component_ids) or not math.fsum(amounts) > 0.0:
            raise ValueError("amounts_mol must match components and have positive total")
        coordinate_ids = tuple(self.coordinate_ids)
        if any(not isinstance(item, str) or not item for item in coordinate_ids):
            raise ValueError("coordinate_ids must contain nonempty names")
        if len(set(coordinate_ids)) != len(coordinate_ids):
            raise ValueError("coordinate_ids must be unique")
        jacobian = tuple(tuple(float(value) for value in row) for row in self.amount_jacobian_mol_per_coordinate)
        if jacobian:
            if len(jacobian) != len(component_ids) or any(
                len(row) != len(coordinate_ids) for row in jacobian
            ):
                raise ValueError("amount_jacobian must have component-by-coordinate shape")
            if any(not math.isfinite(value) for row in jacobian for value in row):
                raise ValueError("amount_jacobian must be finite")
        elif coordinate_ids:
            raise ValueError("coordinate_ids require an amount_jacobian")
        object.__setattr__(self, "component_ids", component_ids)
        object.__setattr__(self, "amounts_mol", amounts)
        object.__setattr__(self, "coordinate_ids", coordinate_ids)
        object.__setattr__(self, "amount_jacobian_mol_per_coordinate", jacobian)

    @property
    def dn_dq(self) -> tuple[tuple[float, ...], ...]:
        """Exact component-by-coordinate derivative of the conserved feed."""

        return self.amount_jacobian_mol_per_coordinate

    def reorder(self, component_ids: Sequence[str]) -> ResolvedFeed:
        """Return the same feed in a validated permutation of component order."""

        requested = _component_ids(component_ids)
        if set(requested) != set(self.component_ids):
            raise ValueError("component permutation must contain the same component identities")
        index = {identity: position for position, identity in enumerate(self.component_ids)}
        order = tuple(index[identity] for identity in requested)
        return ResolvedFeed(
            requested,
            tuple(self.amounts_mol[position] for position in order),
            self.coordinate_ids,
            tuple(self.amount_jacobian_mol_per_coordinate[position] for position in order)
            if self.amount_jacobian_mol_per_coordinate
            else (),
        )


def direct_amounts(
    component_ids: Sequence[str], amounts_mol: Sequence[object] | Mapping[str, object]
) -> ResolvedFeed:
    """Normalize canonical component amounts, optionally supplied as quantities."""

    ids = _component_ids(component_ids)
    return ResolvedFeed(ids, _vector(amounts_mol, ids, "mole", field="amounts_mol", nonnegative=True))


def mole_fractions(
    component_ids: Sequence[str], fractions: Sequence[object] | Mapping[str, object], total_amount_mol: object
) -> ResolvedFeed:
    """Convert normalized mole fractions and a total amount to component moles."""

    ids = _component_ids(component_ids)
    total = _scalar(total_amount_mol, "mole", field="total_amount_mol")
    if total <= 0.0:
        raise ValueError("total_amount_mol must be positive")
    values = _fractions(fractions, ids, field="mole_fractions")
    return ResolvedFeed(ids, tuple(total * value for value in values))


def mass_fractions(
    component_ids: Sequence[str],
    fractions: Sequence[object] | Mapping[str, object],
    total_mass_kg: object,
    molar_masses_kg_per_mol: Sequence[object] | Mapping[str, object],
) -> ResolvedFeed:
    """Convert normalized mass fractions and total mass to component moles."""

    ids = _component_ids(component_ids)
    total_mass = _scalar(total_mass_kg, "kilogram", field="total_mass_kg")
    if total_mass <= 0.0:
        raise ValueError("total_mass_kg must be positive")
    values = _fractions(fractions, ids, field="mass_fractions")
    molar_masses = _vector(
        molar_masses_kg_per_mol,
        ids,
        "kilogram / mole",
        field="molar_masses_kg_per_mol",
    )
    if any(value <= 0.0 for value in molar_masses):
        raise ValueError("molar masses must be positive")
    return ResolvedFeed(
        ids,
        tuple(total_mass * fraction / molar_mass for fraction, molar_mass in zip(values, molar_masses, strict=True)),
    )


def solvent_composition_molalities(
    component_ids: Sequence[str],
    solvent_composition: Mapping[str, object],
    solvent_mass_kg: object,
    solute_molalities_mol_per_kg: Mapping[str, object],
    solvent_molar_masses_kg_per_mol: Mapping[str, object],
) -> ResolvedFeed:
    """Convert a solvent mole basis plus named solute molalities to moles.

    ``solvent_composition`` is normalized over solvent components.  The
    supplied solvent mass is used to determine total solvent moles from the
    composition-weighted solvent molar mass; each named solute amount is then
    ``molality * solvent_mass``.
    """

    ids = _component_ids(component_ids)
    if not isinstance(solvent_composition, Mapping) or not solvent_composition:
        raise TypeError("solvent_composition must be a nonempty component mapping")
    if not isinstance(solute_molalities_mol_per_kg, Mapping):
        raise TypeError("solute_molalities_mol_per_kg must be a component mapping")
    solvent_ids = _component_ids(tuple(solvent_composition), field="solvent component_ids")
    solute_ids = _component_ids(tuple(solute_molalities_mol_per_kg), field="solute component_ids")
    if any(identity not in ids for identity in (*solvent_ids, *solute_ids)):
        raise ValueError("solvent and solute identities must be declared component_ids")
    if set(solvent_ids) & set(solute_ids):
        raise ValueError("solvent and solute identities must be disjoint")
    solvent_fraction = _fractions(solvent_composition, solvent_ids, field="solvent_composition")
    solvent_mass = _scalar(solvent_mass_kg, "kilogram", field="solvent_mass_kg")
    if solvent_mass <= 0.0:
        raise ValueError("solvent_mass_kg must be positive")
    if set(solvent_molar_masses_kg_per_mol) != set(solvent_ids):
        raise ValueError("solvent molar masses must name solvent components exactly once")
    solvent_masses = tuple(
        _scalar(solvent_molar_masses_kg_per_mol[identity], "kilogram / mole", field=f"molar_mass[{identity}]")
        for identity in solvent_ids
    )
    if any(value <= 0.0 for value in solvent_masses):
        raise ValueError("solvent molar masses must be positive")
    average_molar_mass = math.fsum(
        fraction * molar_mass
        for fraction, molar_mass in zip(solvent_fraction, solvent_masses, strict=True)
    )
    total_solvent_moles = solvent_mass / average_molar_mass
    amounts = {identity: 0.0 for identity in ids}
    amounts.update(
        {
            identity: total_solvent_moles * fraction
            for identity, fraction in zip(solvent_ids, solvent_fraction, strict=True)
        }
    )
    for identity, raw_molality in solute_molalities_mol_per_kg.items():
        molality = _scalar(raw_molality, "mole / kilogram", field=f"molality[{identity}]", nonnegative=True)
        amounts[identity] = molality * solvent_mass
    return ResolvedFeed(ids, tuple(amounts[identity] for identity in ids))


@dataclass(frozen=True, slots=True)
class FeedCoordinate:
    """One dimensionless bounded coordinate in an affine conserved feed."""

    identity: str
    bounds: tuple[float, float]

    def __post_init__(self) -> None:
        if not isinstance(self.identity, str) or not self.identity or self.identity != self.identity.strip():
            raise ValueError("feed coordinate identity must be a nonempty name")
        if len(self.bounds) != 2 or any(not math.isfinite(float(value)) for value in self.bounds):
            raise ValueError("feed coordinate bounds must be two finite values")
        lower, upper = (float(value) for value in self.bounds)
        if not lower < upper:
            raise ValueError("feed coordinate lower bound must be below upper bound")
        object.__setattr__(self, "bounds", (lower, upper))


@dataclass(frozen=True, slots=True)
class AffineFeed:
    """Generic conserved feed ``n(q) = base + Bq`` with exact ``dn/dq = B``."""

    component_ids: tuple[str, ...]
    base_amounts_mol: tuple[float, ...]
    basis_matrix_mol_per_coordinate: tuple[tuple[float, ...], ...]
    coordinates: tuple[FeedCoordinate, ...]

    def __post_init__(self) -> None:
        ids = _component_ids(self.component_ids)
        base = tuple(
            _scalar(value, "mole", field=f"base_amounts_mol[{index}]", nonnegative=True)
            for index, value in enumerate(self.base_amounts_mol)
        )
        if len(base) != len(ids):
            raise ValueError("base_amounts_mol must match component_ids")
        coordinates = tuple(self.coordinates)
        if not coordinates or len({item.identity for item in coordinates}) != len(coordinates):
            raise ValueError("affine feed requires unique coordinates")
        matrix = tuple(tuple(float(value) for value in row) for row in self.basis_matrix_mol_per_coordinate)
        if len(matrix) != len(ids) or any(len(row) != len(coordinates) for row in matrix):
            raise ValueError("basis_matrix must have component-by-coordinate shape")
        if any(not math.isfinite(value) for row in matrix for value in row):
            raise ValueError("basis_matrix must be finite")
        object.__setattr__(self, "component_ids", ids)
        object.__setattr__(self, "base_amounts_mol", base)
        object.__setattr__(self, "basis_matrix_mol_per_coordinate", matrix)
        object.__setattr__(self, "coordinates", coordinates)

    @property
    def dn_dq(self) -> tuple[tuple[float, ...], ...]:
        """Exact component-by-coordinate derivative matrix ``B``."""

        return self.basis_matrix_mol_per_coordinate

    def resolve(self, values: Sequence[object] | Mapping[str, object]) -> ResolvedFeed:
        """Evaluate the affine feed at bounded coordinate values."""

        if isinstance(values, Mapping):
            expected = {coordinate.identity for coordinate in self.coordinates}
            if set(values) != expected:
                raise ValueError("affine values must name every coordinate exactly once")
            q = tuple(
                _scalar(values[coordinate.identity], "dimensionless", field=f"q[{coordinate.identity}]")
                for coordinate in self.coordinates
            )
        elif isinstance(values, Sequence) and not isinstance(values, (str, bytes)):
            if len(values) != len(self.coordinates):
                raise ValueError("affine values must match coordinate count")
            q = tuple(
                _scalar(value, "dimensionless", field=f"q[{index}]")
                for index, value in enumerate(values)
            )
        else:
            raise TypeError("affine values must be a sequence or coordinate mapping")
        for value, coordinate in zip(q, self.coordinates, strict=True):
            lower, upper = coordinate.bounds
            if not lower <= value <= upper:
                raise ValueError(f"q[{coordinate.identity}] is outside its bounds")
        amounts = tuple(
            base + math.fsum(entry * value for entry, value in zip(row, q, strict=True))
            for base, row in zip(self.base_amounts_mol, self.basis_matrix_mol_per_coordinate, strict=True)
        )
        return ResolvedFeed(
            self.component_ids,
            amounts,
            tuple(coordinate.identity for coordinate in self.coordinates),
            self.basis_matrix_mol_per_coordinate,
        )

    def reorder(self, component_ids: Sequence[str]) -> AffineFeed:
        """Return the same affine feed in a validated permutation of components."""

        requested = _component_ids(component_ids)
        if set(requested) != set(self.component_ids):
            raise ValueError("component permutation must contain the same component identities")
        index = {identity: position for position, identity in enumerate(self.component_ids)}
        order = tuple(index[identity] for identity in requested)
        return AffineFeed(
            requested,
            tuple(self.base_amounts_mol[position] for position in order),
            tuple(self.basis_matrix_mol_per_coordinate[position] for position in order),
            self.coordinates,
        )


__all__ = (
    "AffineFeed",
    "FeedCoordinate",
    "ResolvedFeed",
    "direct_amounts",
    "mass_fractions",
    "mole_fractions",
    "solvent_composition_molalities",
)
