"""Source-neutral declarations for phase-aware reaction directions.

The runtime equilibrium owner can use this module to validate the structural
part of a reaction declaration before it builds a Gibbs/KKT problem.  A
reaction row is sparse in ``(phase_id, species_id)`` coordinates; its canonical
dense representation is phase-major and species-minor.  Nothing in this
module chooses a solver route, phase topology, or standard-state convention.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass
from numbers import Real

_DEFAULT_RANK_TOLERANCE = 1.0e-12
_DEFAULT_CONSERVATION_TOLERANCE = 1.0e-12


def _identity(value: object, *, field: str) -> str:
    if not isinstance(value, str) or not value or value != value.strip():
        raise ValueError(f"{field} must be a nonempty stripped string")
    return value


def _identities(values: Sequence[str], *, field: str) -> tuple[str, ...]:
    if isinstance(values, (str, bytes)):
        raise TypeError(f"{field} must be a sequence of identities")
    result = tuple(values)
    for index, value in enumerate(result):
        _identity(value, field=f"{field}[{index}]")
    if len(set(result)) != len(result):
        raise ValueError(f"{field} must contain unique identities")
    if not result:
        raise ValueError(f"{field} must not be empty")
    return result


def _finite(value: object, *, field: str) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{field} must be a real number")
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{field} must be finite")
    return result


@dataclass(frozen=True, slots=True)
class PhaseSpeciesTerm:
    """One sparse stoichiometric coefficient at a phase/species coordinate."""

    phase_id: str
    species_id: str
    coefficient: float


@dataclass(frozen=True, slots=True)
class PhaseReactionStandardState:
    """Minimal source-neutral standard-state metadata for one reaction row.

    The activity-scale identity is deliberately opaque to this structural
    module.  EOS code owns its interpretation and any reference
    chemical potentials.
    """

    activity_scale_id: str
    reference_pressure_pa: float | None = None


@dataclass(frozen=True, slots=True)
class PhaseReactionRow:
    """One reaction direction, sparse over phase/species coordinates."""

    identity: str
    terms: tuple[PhaseSpeciesTerm, ...]
    ln_equilibrium_constant: float | None = None
    standard_state: PhaseReactionStandardState | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "terms", tuple(self.terms))


@dataclass(frozen=True, slots=True)
class PhaseReactionSystem:
    """Phase/species declarations and conserved quantities for preflight."""

    phase_ids: tuple[str, ...]
    species_ids: tuple[str, ...]
    reactions: tuple[PhaseReactionRow, ...]
    charges: tuple[float, ...] = ()
    balance_matrix: tuple[tuple[float, ...], ...] = ()

    def __post_init__(self) -> None:
        object.__setattr__(self, "phase_ids", tuple(self.phase_ids))
        object.__setattr__(self, "species_ids", tuple(self.species_ids))
        object.__setattr__(self, "reactions", tuple(self.reactions))
        object.__setattr__(self, "charges", tuple(self.charges))
        object.__setattr__(self, "balance_matrix", tuple(tuple(row) for row in self.balance_matrix))


@dataclass(frozen=True, slots=True)
class PhaseReactionPreflight:
    """Validated canonical direction data consumed by a native owner later."""

    phase_ids: tuple[str, ...]
    species_ids: tuple[str, ...]
    row_identities: tuple[str, ...]
    coordinate_order: tuple[tuple[str, str], ...]
    direction_matrix: tuple[tuple[float, ...], ...]
    charges: tuple[float, ...]
    balance_matrix: tuple[tuple[float, ...], ...]
    rank: int

    @property
    def dimension(self) -> int:
        return len(self.coordinate_order)

    @property
    def reaction_count(self) -> int:
        return len(self.direction_matrix)


def phase_species_coordinate_order(
    phase_ids: Sequence[str], species_ids: Sequence[str]
) -> tuple[tuple[str, str], ...]:
    """Return the canonical phase-major/species-minor coordinate order."""

    phases = _identities(phase_ids, field="phase_ids")
    species = _identities(species_ids, field="species_ids")
    return tuple((phase_id, species_id) for phase_id in phases for species_id in species)


def _matrix_rank(rows: tuple[tuple[float, ...], ...], tolerance: float) -> int:
    if not rows:
        return 0
    work = [list(row) for row in rows]
    row_count = len(work)
    column_count = len(work[0])
    scale = max(1.0, max(abs(value) for row in work for value in row))
    threshold = tolerance * scale
    rank = 0
    for column in range(column_count):
        pivot = max(range(rank, row_count), key=lambda index: abs(work[index][column]))
        pivot_value = work[pivot][column]
        if abs(pivot_value) <= threshold:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        pivot_value = work[rank][column]
        for index in range(rank + 1, row_count):
            factor = work[index][column] / pivot_value
            if factor == 0.0:
                continue
            for entry in range(column, column_count):
                work[index][entry] -= factor * work[rank][entry]
        rank += 1
        if rank == row_count:
            break
    return rank


def _validate_tolerance(value: float, *, field: str) -> float:
    result = _finite(value, field=field)
    if result <= 0.0:
        raise ValueError(f"{field} must be positive")
    return result


def _validated_direction_matrix(
    system: PhaseReactionSystem,
) -> tuple[tuple[float, ...], ...]:
    phase_ids = _identities(system.phase_ids, field="phase_ids")
    species_ids = _identities(system.species_ids, field="species_ids")
    phase_index = {identity: index for index, identity in enumerate(phase_ids)}
    species_index = {identity: index for index, identity in enumerate(species_ids)}
    species_count = len(species_ids)
    rows: list[tuple[float, ...]] = []
    identities: set[str] = set()
    for row_index, reaction in enumerate(system.reactions):
        if not isinstance(reaction, PhaseReactionRow):
            raise TypeError(f"reactions[{row_index}] must be a PhaseReactionRow")
        identity = _identity(reaction.identity, field=f"reactions[{row_index}].identity")
        if identity in identities:
            raise ValueError(f"reaction identities must be unique: {identity!r}")
        identities.add(identity)
        if not reaction.terms:
            raise ValueError(f"reactions[{row_index}] must contain at least one term")
        values = [0.0] * (len(phase_ids) * species_count)
        coordinates: set[tuple[str, str]] = set()
        for term_index, term in enumerate(reaction.terms):
            if not isinstance(term, PhaseSpeciesTerm):
                raise TypeError(f"reactions[{row_index}].terms[{term_index}] has invalid type")
            phase_id = _identity(
                term.phase_id,
                field=f"reactions[{row_index}].terms[{term_index}].phase_id",
            )
            species_id = _identity(
                term.species_id,
                field=f"reactions[{row_index}].terms[{term_index}].species_id",
            )
            if phase_id not in phase_index:
                raise ValueError(f"reaction term names unknown phase {phase_id!r}")
            if species_id not in species_index:
                raise ValueError(f"reaction term names unknown species {species_id!r}")
            coordinate = (phase_id, species_id)
            if coordinate in coordinates:
                raise ValueError(f"reaction row {identity!r} repeats coordinate {coordinate!r}")
            coordinates.add(coordinate)
            coefficient = _finite(
                term.coefficient,
                field=f"reactions[{row_index}].terms[{term_index}].coefficient",
            )
            if coefficient == 0.0:
                raise ValueError(f"reaction row {identity!r} contains a zero coefficient")
            offset = phase_index[phase_id] * species_count + species_index[species_id]
            values[offset] = coefficient
        if reaction.ln_equilibrium_constant is not None:
            _finite(
                reaction.ln_equilibrium_constant,
                field=f"reactions[{row_index}].ln_equilibrium_constant",
            )
        if reaction.standard_state is not None:
            if not isinstance(reaction.standard_state, PhaseReactionStandardState):
                raise TypeError(f"reaction row {identity!r} has invalid standard_state")
            _identity(
                reaction.standard_state.activity_scale_id,
                field=f"reactions[{row_index}].standard_state.activity_scale_id",
            )
            pressure = reaction.standard_state.reference_pressure_pa
            if pressure is not None and _finite(pressure, field="reference_pressure_pa") <= 0.0:
                raise ValueError("reference_pressure_pa must be positive")
        rows.append(tuple(values))
    return tuple(rows)


def canonical_stacked_directions(system: PhaseReactionSystem) -> tuple[tuple[float, ...], ...]:
    """Return validated dense reaction rows in canonical stacked coordinates."""

    return _validated_direction_matrix(system)


def preflight_phase_reaction_system(
    system: PhaseReactionSystem,
    *,
    rank_tolerance: float = _DEFAULT_RANK_TOLERANCE,
    conservation_tolerance: float = _DEFAULT_CONSERVATION_TOLERANCE,
) -> PhaseReactionPreflight:
    """Validate sparse phase reactions and return their canonical directions.

    ``balance_matrix`` is oriented as conserved-quantity-by-species.  Each
    balance is global across phases, while charge is checked independently in
    every phase.  Consequently a cross-phase transfer must preserve charge in
    each participating phase; global charge cancellation alone is insufficient.
    """

    if not isinstance(system, PhaseReactionSystem):
        raise TypeError("system must be a PhaseReactionSystem")
    rank_tolerance = _validate_tolerance(rank_tolerance, field="rank_tolerance")
    conservation_tolerance = _validate_tolerance(
        conservation_tolerance, field="conservation_tolerance"
    )
    phase_ids = _identities(system.phase_ids, field="phase_ids")
    species_ids = _identities(system.species_ids, field="species_ids")
    species_count = len(species_ids)
    charges = tuple(system.charges)
    if not charges:
        charges = (0.0,) * species_count
    if len(charges) != species_count:
        raise ValueError("charges must match species_ids")
    charges = tuple(_finite(value, field=f"charges[{index}]") for index, value in enumerate(charges))
    balances = tuple(tuple(row) for row in system.balance_matrix)
    for row_index, row in enumerate(balances):
        if len(row) != species_count:
            raise ValueError(f"balance_matrix[{row_index}] must match species_ids")
        for value_index, value in enumerate(row):
            _finite(value, field=f"balance_matrix[{row_index}][{value_index}]")
    directions = _validated_direction_matrix(system)
    rank = _matrix_rank(directions, rank_tolerance)
    if rank != len(directions):
        raise ValueError(
            f"reaction directions must be independent: rank {rank} != {len(directions)}"
        )
    for row_index, direction in enumerate(directions):
        for phase_index in range(len(phase_ids)):
            start = phase_index * species_count
            phase_direction = direction[start : start + species_count]
            charge_change = math.fsum(
                coefficient * charge
                for coefficient, charge in zip(phase_direction, charges, strict=True)
            )
            if abs(charge_change) > conservation_tolerance:
                raise ValueError(
                    f"reaction {system.reactions[row_index].identity!r} violates charge "
                    f"in phase {phase_ids[phase_index]!r}"
                )
        for balance_index, balance in enumerate(balances):
            change = math.fsum(
                direction[phase_index * species_count + species_index] * balance[species_index]
                for phase_index in range(len(phase_ids))
                for species_index in range(species_count)
            )
            if abs(change) > conservation_tolerance:
                raise ValueError(
                    f"reaction {system.reactions[row_index].identity!r} violates "
                    f"global balance {balance_index}"
                )
    return PhaseReactionPreflight(
        phase_ids,
        species_ids,
        tuple(reaction.identity for reaction in system.reactions),
        phase_species_coordinate_order(phase_ids, species_ids),
        directions,
        charges,
        balances,
        rank,
    )


__all__ = [
    "PhaseReactionPreflight",
    "PhaseReactionRow",
    "PhaseReactionStandardState",
    "PhaseReactionSystem",
    "PhaseSpeciesTerm",
    "canonical_stacked_directions",
    "phase_species_coordinate_order",
    "preflight_phase_reaction_system",
]
