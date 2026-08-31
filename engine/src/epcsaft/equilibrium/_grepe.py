from __future__ import annotations

import hashlib
import json
import math
from dataclasses import asdict, dataclass, replace
from typing import Any, Literal, cast

from pint import Quantity

import epcsaft

from .. import _core
from .._canonical import _canonical_json
from .._trials import NON_EVALUABLE_TRIAL_CODES, NonEvaluableTrial
from ._api import (
    ActiveParameterTrialReceipt,
    AllComponents,
    ChemicalEquilibriumProblem,
    ChemicalStandardState,
    DeclaredComponents,
    DeclaredNeutralComponents,
    EquilibriumContinuationState,
    EquilibriumPhaseState,
    NeutralComponents,
    PhaseSupport,
    _support_active_parameters,
)
from ._feed import AffineFeed, ResolvedFeed
from ._phase_reactions import PhaseReactionSystem, preflight_phase_reaction_system

type AmountRole = Literal["finite", "incipient"]
type ReactivePhaseSupport = PhaseSupport | DeclaredComponents
type OutputSupport = Literal[
    "real",
    "positive",
    "nonnegative",
    "left_censored",
    "right_censored",
]
type EquilibriumSelector = Literal[
    "system.pressure",
    "system.temperature",
    "phase.amount",
    "phase.mole_fraction",
    "phase.molality",
    "phase.molar_density",
    "phase.compressibility_factor",
    "phase.pressure",
    "phase.fugacity",
    "phase.log_activity",
    "phase.log_fugacity",
    "phase.chemical_potential_over_rt",
    "phase.partial_pressure",
]
type IntensiveBoundarySelector = Literal[
    "system.pressure",
    "phase.log_activity",
    "phase.log_fugacity",
    "phase.chemical_potential_over_rt",
]

_EQUILIBRIUM_SELECTORS = frozenset(
    {
        "system.pressure",
        "system.temperature",
        "phase.amount",
        "phase.mole_fraction",
        "phase.molality",
        "phase.molar_density",
        "phase.compressibility_factor",
        "phase.pressure",
        "phase.fugacity",
        "phase.log_activity",
        "phase.log_fugacity",
        "phase.chemical_potential_over_rt",
        "phase.partial_pressure",
    }
)
_SELECTOR_UNITS = {
    "system.pressure": "pascal",
    "system.temperature": "kelvin",
    "phase.amount": "mole",
    "phase.mole_fraction": "dimensionless",
    "phase.molality": "mole / kilogram",
    "phase.molar_density": "mole / meter**3",
    "phase.compressibility_factor": "dimensionless",
    "phase.pressure": "pascal",
    "phase.fugacity": "pascal",
    "phase.log_activity": "dimensionless",
    "phase.log_fugacity": "dimensionless",
    "phase.chemical_potential_over_rt": "dimensionless",
    "phase.partial_pressure": "pascal",
}


@dataclass(frozen=True, slots=True)
class Fixed:
    """A fixed intensive-state value with explicit units."""

    value: Quantity[Any]


@dataclass(frozen=True, slots=True)
class Solved:
    """A solved intensive-state value with explicit initial point, bounds, and starts."""

    initial: Quantity[Any]
    bounds: tuple[Quantity[Any], Quantity[Any]]
    starts: tuple[Quantity[Any], ...]

    def __post_init__(self) -> None:
        object.__setattr__(self, "bounds", tuple(self.bounds))
        object.__setattr__(self, "starts", tuple(self.starts))
        if len(self.bounds) != 2 or not self.starts:
            raise ValueError("a solved state requires two bounds and at least one start")


@dataclass(frozen=True, slots=True)
class ProviderModel:
    """Installed EOS model declaration for one phase role."""

    admissible_packing_fraction_interval: tuple[float, float]
    reference_id: str

    def __post_init__(self) -> None:
        interval = tuple(self.admissible_packing_fraction_interval)
        object.__setattr__(self, "admissible_packing_fraction_interval", interval)
        if (
            len(interval) != 2
            or not 0.0 <= interval[0] < interval[1] < 1.0
            or not self.reference_id
        ):
            raise ValueError("provider phase model requires a physical packing interval")


@dataclass(frozen=True, slots=True)
class IdealGasModel:
    """Explicit comparison-only ideal-gas phase model."""

    reference_id: str

    def __post_init__(self) -> None:
        if not self.reference_id:
            raise ValueError("ideal-gas phase model requires a reference identity")


type ReactivePhaseModel = ProviderModel | IdealGasModel


@dataclass(frozen=True, slots=True)
class FinitePhaseStart:
    """Strictly positive local start for one finite GREPE phase."""

    amounts_mol: tuple[float, ...]
    molar_volume_m3_per_mol: float

    def __post_init__(self) -> None:
        object.__setattr__(self, "amounts_mol", tuple(self.amounts_mol))
        if (
            not self.amounts_mol
            or any(not math.isfinite(value) or value < 0.0 for value in self.amounts_mol)
            or not any(value > 0.0 for value in self.amounts_mol)
            or not math.isfinite(self.molar_volume_m3_per_mol)
            or self.molar_volume_m3_per_mol <= 0.0
        ):
            raise ValueError(
                "finite phase starts require nonnegative amounts, positive total, and molar volume"
            )


@dataclass(frozen=True, slots=True)
class ReactivePhase:
    """One declared fluid phase, amount role, support, model, and local anchor."""

    identity: str
    fluid_role: Literal["liquid", "vapor"]
    amount_role: AmountRole
    support: ReactivePhaseSupport
    model: ReactivePhaseModel
    start: FinitePhaseStart | None = None

    def __post_init__(self) -> None:
        if (
            not self.identity
            or self.fluid_role not in {"liquid", "vapor"}
            or self.amount_role not in {"finite", "incipient"}
            or not isinstance(
                self.support,
                (AllComponents, NeutralComponents, DeclaredNeutralComponents, DeclaredComponents),
            )
            or not isinstance(self.model, (ProviderModel, IdealGasModel))
        ):
            raise ValueError("reactive phase declaration is invalid")
        if self.amount_role == "incipient" and self.start is not None:
            raise ValueError("incipient phases have no finite amount start")
        if isinstance(self.model, IdealGasModel) and self.fluid_role != "vapor":
            raise ValueError("an ideal-gas phase must have the vapor role")
        if self.start is not None and self.amount_role != "finite":
            raise ValueError("only finite phases accept finite phase starts")


@dataclass(frozen=True, slots=True)
class EquilibriumOutput:
    """One source-neutral output selector and contraction."""

    identity: str
    selector: EquilibriumSelector
    unit: str
    basis: str
    phase_identity: str | None = None
    coefficients: tuple[float, ...] = ()
    solvent_mass_coefficients_kg_per_mol: tuple[float, ...] = ()
    support: OutputSupport = "real"
    censor_limit: float | None = None
    aggregate_identity: str | None = None
    covariance_identity: str | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "coefficients", tuple(self.coefficients))
        object.__setattr__(
            self,
            "solvent_mass_coefficients_kg_per_mol",
            tuple(self.solvent_mass_coefficients_kg_per_mol),
        )
        if (
            not self.identity
            or self.selector not in _EQUILIBRIUM_SELECTORS
            or self.unit != _SELECTOR_UNITS.get(self.selector)
            or not self.basis
        ):
            raise ValueError("equilibrium output identity or selector is invalid")
        if self.selector in {"system.pressure", "system.temperature"}:
            if self.phase_identity is not None or self.coefficients:
                raise ValueError("system intensives have no phase or contraction coefficients")
        elif not self.phase_identity:
            raise ValueError("phase output selectors require a phase identity")
        if self.support not in {
            "real",
            "positive",
            "nonnegative",
            "left_censored",
            "right_censored",
        }:
            raise ValueError("equilibrium output support is invalid")
        if self.support.endswith("censored") != (self.censor_limit is not None):
            raise ValueError("censored outputs require exactly one censor limit")
        if bool(self.aggregate_identity) != bool(self.covariance_identity):
            raise ValueError("aggregate and covariance output identities must be paired")


@dataclass(frozen=True, slots=True)
class IntensiveBoundary:
    """One imposed intensive value conjugate to an affine feed coordinate."""

    identity: str
    coordinate_identity: str
    phase_identity: str | None
    selector: IntensiveBoundarySelector
    coefficients: tuple[float, ...]
    target: float
    initial: float

    def __post_init__(self) -> None:
        object.__setattr__(self, "coefficients", tuple(float(value) for value in self.coefficients))
        phase_selector = self.selector in {
            "phase.log_activity",
            "phase.log_fugacity",
            "phase.chemical_potential_over_rt",
        }
        pressure_selector = self.selector == "system.pressure"
        if (
            not self.identity
            or not self.coordinate_identity
            or not math.isfinite(float(self.target))
            or not math.isfinite(float(self.initial))
            or (
                pressure_selector
                and (
                    self.phase_identity is not None
                    or self.coefficients
                    or float(self.target) <= 0.0
                )
            )
            or (
                phase_selector
                and (
                    not self.phase_identity
                    or not self.coefficients
                    or not any(value != 0.0 for value in self.coefficients)
                    or not all(math.isfinite(value) for value in self.coefficients)
                )
            )
            or not (phase_selector or pressure_selector)
        ):
            raise ValueError("intensive boundary declaration is invalid")
        object.__setattr__(self, "target", float(self.target))
        object.__setattr__(self, "initial", float(self.initial))


@dataclass(frozen=True, slots=True)
class _FinitePhaseOwnerDescriptor:
    rows: tuple[EquilibriumOutput, ...]
    row_identities: tuple[str, ...]
    active_parameter_identities: tuple[str, ...]
    active_parameter_units: tuple[str, ...]
    parameter_fingerprint: str
    active_parameter_selection_fingerprint: str
    problem_fingerprint: str
    reaction_fingerprint: str
    topology_fingerprint: str
    reference_fingerprint: str
    domain_fingerprint: str
    jacobian_shape: tuple[int, int]
    fingerprint: str


@dataclass(frozen=True, slots=True)
class GeneralReactiveEquilibriumProblem:
    """One role-declared GREPE condition over the implemented formulations."""

    identity: str
    temperature: Fixed | Solved
    pressure: Fixed | Solved
    phases: tuple[ReactivePhase, ...]
    reaction_system: ChemicalEquilibriumProblem
    reaction_phase_ids: tuple[str, ...]
    outputs: tuple[EquilibriumOutput, ...]
    continuation_identity: str | None = None
    continuation_state: EquilibriumContinuationState | None = None
    feed: ResolvedFeed | AffineFeed | None = None
    phase_reactions: PhaseReactionSystem | None = None
    intensive_boundaries: tuple[IntensiveBoundary, ...] = ()

    def __post_init__(self) -> None:
        object.__setattr__(self, "phases", tuple(self.phases))
        object.__setattr__(self, "reaction_phase_ids", tuple(self.reaction_phase_ids))
        object.__setattr__(self, "outputs", tuple(self.outputs))
        object.__setattr__(self, "intensive_boundaries", tuple(self.intensive_boundaries))
        phase_ids = tuple(phase.identity for phase in self.phases)
        output_ids = tuple(output.identity for output in self.outputs)
        if not isinstance(self.reaction_system, ChemicalEquilibriumProblem):
            raise TypeError("general reactive equilibrium request is incomplete")
        if self.feed is not None:
            if not isinstance(self.feed, (ResolvedFeed, AffineFeed)):
                raise TypeError("GREPE feed must be resolved or affine")
            if set(self.feed.component_ids) != set(self.reaction_system.species_ids):
                raise ValueError("GREPE feed identities must equal the reaction species identities")
        if self.phase_reactions is not None:
            if not isinstance(self.phase_reactions, PhaseReactionSystem):
                raise TypeError("GREPE phase reactions must be a phase/species direction system")
            phase_reaction_preflight = preflight_phase_reaction_system(self.phase_reactions)
            if phase_reaction_preflight.species_ids != self.reaction_system.species_ids:
                raise ValueError(
                    "GREPE phase-reaction species order must equal the reaction system"
                )
            if phase_reaction_preflight.phase_ids != phase_ids:
                raise ValueError("GREPE phase-reaction order must equal the declared phase order")
            if phase_reaction_preflight.charges != tuple(self.reaction_system.charges):
                raise ValueError("GREPE phase-reaction charges differ from the reaction system")
            if phase_reaction_preflight.balance_matrix != self.reaction_system.balance_matrix:
                raise ValueError("GREPE phase-reaction balances differ from the reaction system")
        has_reactions = bool(self.reaction_system.reaction_matrix)
        if (
            not self.identity
            or not self.phases
            or len(set(phase_ids)) != len(phase_ids)
            or not self.outputs
            or len(set(output_ids)) != len(output_ids)
            or has_reactions != bool(self.reaction_phase_ids)
            or len(set(self.reaction_phase_ids)) != len(self.reaction_phase_ids)
            or any(identity not in phase_ids for identity in self.reaction_phase_ids)
        ):
            raise ValueError("general reactive equilibrium request is incomplete")
        finite_ids = {phase.identity for phase in self.phases if phase.amount_role == "finite"}
        if any(identity not in finite_ids for identity in self.reaction_phase_ids):
            raise ValueError("reactions may be declared only in finite phases")
        for output in self.outputs:
            if output.phase_identity is not None and output.phase_identity not in phase_ids:
                raise ValueError("equilibrium output names an undeclared phase")
        if self.intensive_boundaries:
            if not isinstance(self.feed, AffineFeed):
                raise ValueError("intensive boundaries require an affine feed")
            if len(self.intensive_boundaries) != 1:
                raise ValueError("GREPE supports exactly one intensive boundary per request")
            boundary_ids = tuple(boundary.identity for boundary in self.intensive_boundaries)
            coordinate_ids = tuple(
                boundary.coordinate_identity for boundary in self.intensive_boundaries
            )
            feed_coordinate_ids = tuple(coordinate.identity for coordinate in self.feed.coordinates)
            if (
                len(set(boundary_ids)) != len(boundary_ids)
                or len(set(coordinate_ids)) != len(coordinate_ids)
                or set(coordinate_ids) != set(feed_coordinate_ids)
            ):
                raise ValueError(
                    "intensive boundaries must name each affine coordinate exactly once"
                )
            coordinate_by_id = {
                coordinate.identity: coordinate for coordinate in self.feed.coordinates
            }
            for boundary in self.intensive_boundaries:
                if boundary.selector == "system.pressure":
                    valid_boundary = boundary.phase_identity is None and not boundary.coefficients
                else:
                    valid_boundary = boundary.phase_identity in finite_ids and len(
                        boundary.coefficients
                    ) == len(self.reaction_system.species_ids)
                if not valid_boundary:
                    raise ValueError("intensive boundary phase or component contraction is invalid")
                lower, upper = coordinate_by_id[boundary.coordinate_identity].bounds
                if not lower <= boundary.initial <= upper:
                    raise ValueError("intensive boundary initial lies outside feed bounds")
        elif isinstance(self.feed, AffineFeed):
            raise ValueError("an affine GREPE feed requires intensive boundaries")
        incipient_count = sum(phase.amount_role == "incipient" for phase in self.phases)
        if incipient_count and not self.continuation_identity:
            raise ValueError("an incipient formulation requires a continuation identity")
        if self.continuation_state is not None and not self.continuation_identity:
            raise ValueError("a continuation state requires a continuation identity")


def _effective_reaction_system(
    problem: GeneralReactiveEquilibriumProblem,
) -> ChemicalEquilibriumProblem:
    """Return the chemical declaration with the request's normalized fixed feed.

    Composition-basis conversion is deliberately outside the thermodynamic
    solve.  Once resolved, the feed is just a canonical amount vector and its
    conserved totals are obtained from the declared balance matrix.
    """

    reaction = problem.reaction_system
    if problem.feed is None:
        return reaction
    if isinstance(problem.feed, AffineFeed):
        raise TypeError("affine feed must be resolved by the intensive-boundary solver")
    feed = problem.feed.reorder(reaction.species_ids)
    conserved_totals = tuple(
        math.fsum(
            coefficient * amount for coefficient, amount in zip(row, feed.amounts_mol, strict=True)
        )
        for row in reaction.balance_matrix
    )
    return replace(
        reaction,
        feed_amounts_mol=feed.amounts_mol,
        conserved_totals=conserved_totals,
    )


@dataclass(frozen=True, slots=True)
class GeneralEquilibriumRowResult:
    identity: str
    selector: EquilibriumSelector
    unit: str
    basis: str
    status: Literal["evaluated", "non_evaluable"]
    value: float | None
    jacobian: tuple[float, ...] | None
    reason: str


@dataclass(frozen=True, slots=True)
class GeneralReactiveEquilibriumDescriptor:
    """Closed GREPE request/output/derivative identity consumed by Regression."""

    identity: Literal["equilibrium.general-reactive"]
    version: int
    problem_identity: str
    component_ids: tuple[str, ...]
    phase_identities: tuple[str, ...]
    outputs: tuple[EquilibriumOutput, ...]
    active_parameter_identities: tuple[str, ...]
    active_parameter_units: tuple[str, ...]
    parameter_fingerprint: str
    active_parameter_selection_fingerprint: str
    problem_fingerprint: str
    reaction_fingerprint: str
    topology_fingerprint: str
    reference_fingerprint: str
    domain_fingerprint: str
    certificate_scope: str
    jacobian_shape: tuple[int, int]
    fingerprint: str
    formulation: Literal["finite", "finite_incipient", "intensive_boundary", "temperature_boundary"]


@dataclass(frozen=True, slots=True)
class GeneralEquilibriumOwnerObservation:
    target_identity: str
    prediction_identity: str
    aggregate_identity: str | None = None
    covariance_identity: str | None = None

    def __post_init__(self) -> None:
        if not self.target_identity or not self.prediction_identity:
            raise ValueError("GREPE owner observation identities are required")
        if bool(self.aggregate_identity) != bool(self.covariance_identity):
            raise ValueError("GREPE aggregate and covariance identities must be paired")


@dataclass(frozen=True, slots=True)
class GeneralReactiveEquilibriumResult:
    """One generic result from the configured GREPE formulation."""

    problem: GeneralReactiveEquilibriumProblem
    descriptor: GeneralReactiveEquilibriumDescriptor | None
    status: Literal["evaluated", "non_evaluable"]
    solver_status: str
    numerical_status: str
    physical_status: str
    provider_domain_status: str
    rows: tuple[GeneralEquilibriumRowResult, ...]
    phases: tuple[EquilibriumPhaseState, ...]
    failure: NonEvaluableTrial | None
    certificate_scope: str
    globality_claim: Literal["not_established"]
    trial: ActiveParameterTrialReceipt | None
    branch_signature: tuple[object, ...] | None
    evidence: tuple[tuple[str, float | int | str], ...]
    resolved_policy: tuple[tuple[str, float], ...]
    continuation_state: EquilibriumContinuationState | None = None

    @property
    def values(self) -> tuple[float | None, ...]:
        return tuple(row.value for row in self.rows)

    @property
    def jacobian(self) -> tuple[tuple[float, ...] | None, ...]:
        return tuple(row.jacobian for row in self.rows)


def _with_continuation_state(
    result: GeneralReactiveEquilibriumResult,
) -> GeneralReactiveEquilibriumResult:
    if (
        result.status != "evaluated"
        or result.descriptor is None
        or result.trial is None
        or not result.phases
    ):
        return result
    return replace(
        result,
        continuation_state=EquilibriumContinuationState(
            tuple(phase.identity for phase in result.problem.phases),
            result.phases,
            next(
                (
                    row.value
                    for row in result.rows
                    if row.selector == "system.pressure" and row.value is not None
                ),
                result.phases[0].pressure_pa,
            ),
            result.descriptor.problem_fingerprint,
            result.descriptor.reaction_fingerprint,
            result.descriptor.topology_fingerprint,
            result.trial.parameter_fingerprint,
        ),
    )


def _result_evidence(owner: object | None) -> tuple[tuple[str, float | int | str], ...]:
    if isinstance(owner, _FinitePhaseOwnerResult):
        names: tuple[str, ...] = (
            "objective_over_rt",
            "balance_inf_norm",
            "pressure_relative_inf_norm",
            "reaction_affinity_inf_norm",
            "transfer_chemical_potential_inf_norm",
            "kkt_stationarity_inf_norm",
            "source_reference_representation_residual_inf_norm",
            "phase_distinctness",
            "reduced_hessian_min_eigenvalue",
            "declared_direction_rank",
            "admissible_nullity",
            "optimizer_iterations",
            "provider_evaluation_count",
        )
        return tuple((name, cast(float | int, getattr(owner, name))) for name in names)
    return ()


def _native_field(native: object, name: str) -> object:
    return getattr(native, name)


def _native_resolved_policy(native: object) -> tuple[tuple[str, float], ...]:
    records = tuple(
        (str(record[0]), float(cast(float, record[1])))
        for record in cast(
            tuple[tuple[object, object], ...], _native_field(native, "resolved_policy")
        )
    )
    if (
        not records
        or len({name for name, _ in records}) != len(records)
        or any(not name or not math.isfinite(value) or value < 0.0 for name, value in records)
    ):
        raise ValueError("native GREPE numerical policy is invalid")
    return records


def _quantity(value: Quantity[Any], unit: str, name: str) -> float:
    try:
        result = float(value.to(unit).magnitude)
    except (AttributeError, TypeError, ValueError) as error:
        raise ValueError(f"GREPE {name} must have units convertible to {unit}") from error
    if not math.isfinite(result):
        raise ValueError(f"GREPE {name} must be finite")
    return result


def _matrix_rank(rows: tuple[tuple[float, ...], ...]) -> int:
    if not rows:
        return 0
    work = [list(row) for row in rows]
    row_count = len(work)
    column_count = len(work[0])
    scale = max(1.0, max(abs(value) for row in work for value in row))
    rank = 0
    for column in range(column_count):
        pivot = max(range(rank, row_count), key=lambda index: abs(work[index][column]))
        if abs(work[pivot][column]) <= 1.0e-12 * scale:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        pivot_value = work[rank][column]
        for index in range(rank + 1, row_count):
            factor = work[index][column] / pivot_value
            for entry in range(column, column_count):
                work[index][entry] -= factor * work[rank][entry]
        rank += 1
        if rank == row_count:
            break
    return rank


def _preflight(model: epcsaft.Mixture, problem: GeneralReactiveEquilibriumProblem) -> None:
    reaction = _effective_reaction_system(problem)
    component_ids = model.component_ids
    if reaction.species_ids != component_ids:
        raise ValueError("GREPE species order must equal the Provider component order")
    component_count = len(component_ids)
    phase_supports = {
        phase.identity: frozenset(
            _phase_supported_component_ids(phase, component_ids, tuple(reaction.charges))
        )
        for phase in problem.phases
    }
    if any(not supported for supported in phase_supports.values()):
        raise ValueError("GREPE finite phase support must not be empty")
    if problem.phase_reactions is not None:
        for row in problem.phase_reactions.reactions:
            for term in row.terms:
                if term.species_id not in phase_supports[term.phase_id]:
                    raise ValueError(
                        f"phase direction {row.identity!r} names species {term.species_id!r} "
                        f"outside phase {term.phase_id!r} support"
                    )
    if (
        len(reaction.charges) != component_count
        or len(reaction.molar_masses_kg_per_mol) != component_count
        or len(reaction.feed_amounts_mol) != component_count
        or any(len(row) != component_count for row in reaction.balance_matrix)
        or any(len(row) != component_count for row in reaction.reaction_matrix)
        or len(reaction.equilibrium_constants) != len(reaction.reaction_matrix)
        or _matrix_rank(reaction.reaction_matrix) != len(reaction.reaction_matrix)
    ):
        raise ValueError("GREPE reaction system dimensions or independent rank are invalid")
    for stoichiometry in reaction.reaction_matrix:
        charge_change = math.fsum(
            coefficient * charge
            for coefficient, charge in zip(stoichiometry, reaction.charges, strict=True)
        )
        conservation = tuple(
            math.fsum(
                coefficient * conserved
                for coefficient, conserved in zip(stoichiometry, row, strict=True)
            )
            for row in reaction.balance_matrix
        )
        if abs(charge_change) > 1.0e-12 or any(abs(value) > 1.0e-12 for value in conservation):
            raise ValueError("GREPE reaction stoichiometry violates charge or conservation")
    standard = reaction.source_standard_state
    if not reaction.reaction_matrix and standard is not None:
        raise ValueError("a nonreactive GREPE request cannot declare a reaction standard state")
    if standard is not None:
        eligible = {
            phase.identity: phase
            for phase in problem.phases
            if phase.identity in problem.reaction_phase_ids
        }
        if standard.source_reference_phase != "liquid" or any(
            phase.fluid_role != "liquid" for phase in eligible.values()
        ):
            raise ValueError("reaction standard state is inconsistent with its eligible phase")


def _standard_state_payload(standard: ChemicalStandardState | None) -> object:
    if standard is None:
        return None
    return {
        "id": standard.id,
        "activity_scale_id": standard.activity_scale_id,
        "log_activity_scale_factors": standard.log_activity_scale_factors,
        "reference_pressure_pa": standard.reference_pressure_pa,
        "source_reference_component_ids": standard.source_reference_component_ids,
        "source_reference_solvent_composition": (standard.source_reference_solvent_composition),
        "source_reference_ion_pairs": tuple(
            (pair.component_ids, pair.stoichiometry) for pair in standard.source_reference_ion_pairs
        ),
        "source_reference_pure_components": tuple(
            (component.component_id, component.phase)
            for component in standard.source_reference_pure_components
        ),
        "source_reference_phase": standard.source_reference_phase,
        "source_reference_convention": standard.source_reference_convention,
        "source_reference_activity_convention_id": (
            standard.source_reference_activity_convention_id
        ),
        "source_reference_standard_molality_mol_per_kg": (
            standard.source_reference_standard_molality_mol_per_kg
        ),
    }


def _source_reference_fingerprint(
    standard: ChemicalStandardState | None,
    *,
    active_parameters_requested: bool,
) -> str:
    if standard is None:
        return ""
    required_derivatives = (
        ("pressure", "parameters") if active_parameters_requested else ("pressure",)
    )
    payload = {
        "schema": "epcsaft-source-reference-declaration-v1",
        "reference_state_id": standard.id,
        "component_ids": standard.source_reference_component_ids,
        "solvent_mole_fractions": standard.source_reference_solvent_composition,
        "ion_pairs": tuple(
            (pair.component_ids, pair.stoichiometry) for pair in standard.source_reference_ion_pairs
        ),
        "phase": standard.source_reference_phase,
        "reference_convention": standard.source_reference_convention,
        "activity_convention_id": standard.source_reference_activity_convention_id,
        "standard_molality_mol_per_kg": (standard.source_reference_standard_molality_mol_per_kg),
        "reference_pressure_pa": standard.reference_pressure_pa,
        "required_derivatives": required_derivatives,
        "pure_components": tuple(
            (component.component_id, component.phase)
            for component in standard.source_reference_pure_components
        ),
    }
    return "sha256:" + hashlib.sha256(_canonical_json(payload).encode()).hexdigest()


def _bubble_quantity(
    output: EquilibriumOutput, finite: ReactivePhase, incipient: ReactivePhase
) -> str:
    if output.selector == "system.pressure":
        return "bubble_pressure_pa"
    if output.phase_identity == finite.identity:
        quantities = {
            "phase.amount": "liquid_amount_mol",
            "phase.mole_fraction": "liquid_mole_fraction",
            "phase.molality": "liquid_molality_mol_per_kg",
            "phase.molar_density": "liquid_molar_density_mol_per_m3",
        }
    elif output.phase_identity == incipient.identity:
        quantities = {
            "phase.mole_fraction": "vapor_mole_fraction",
            "phase.partial_pressure": "vapor_partial_pressure_pa",
        }
    else:
        quantities = {}
    try:
        return quantities[output.selector]
    except KeyError as error:
        raise ValueError("output selector is unavailable for the declared phase role") from error


def _generic_rows(
    outputs: tuple[EquilibriumOutput, ...],
    rows: tuple[object, ...],
) -> tuple[GeneralEquilibriumRowResult, ...]:
    return tuple(
        GeneralEquilibriumRowResult(
            identity=output.identity,
            selector=output.selector,
            unit=output.unit,
            basis=output.basis,
            status=cast(Any, row).status,
            value=cast(Any, row).value,
            jacobian=cast(Any, row).jacobian,
            reason=(
                cast(Any, row).reason
                if hasattr(row, "reason")
                else ""
                if cast(Any, row).failure is None
                else cast(Any, row).failure.diagnostic
            ),
        )
        for output, row in zip(outputs, rows, strict=True)
    )


def _fingerprint(schema: str, value: object) -> str:
    payload = _canonical_json({"schema": schema, "value": value}).encode("utf-8")
    return f"sha256:{hashlib.sha256(payload).hexdigest()}"


def _state_role_payload(role: Fixed | Solved, unit: str, name: str) -> dict[str, object]:
    if isinstance(role, Fixed):
        return {"role": "fixed", "value": _quantity(role.value, unit, name), "unit": unit}
    return {
        "role": "solved",
        "initial": _quantity(role.initial, unit, f"{name} initial"),
        "bounds": tuple(_quantity(item, unit, f"{name} bound") for item in role.bounds),
        "starts": tuple(_quantity(item, unit, f"{name} start") for item in role.starts),
        "unit": unit,
    }


def _phase_payload(phase: ReactivePhase) -> dict[str, object]:
    if isinstance(phase.support, DeclaredComponents):
        support = {
            "kind": "declared_components",
            "component_ids": phase.support.component_ids,
        }
    elif isinstance(phase.support, DeclaredNeutralComponents):
        support = {
            "kind": "declared_neutral_components",
            "component_ids": phase.support.component_ids,
        }
    elif isinstance(phase.support, NeutralComponents):
        support = {"kind": "neutral_components", "component_ids": ()}
    else:
        support = {"kind": "all_components", "component_ids": ()}
    if isinstance(phase.model, ProviderModel):
        model: dict[str, object] = {
            "kind": "provider",
            "reference_id": phase.model.reference_id,
            "admissible_packing_fraction_interval": (
                phase.model.admissible_packing_fraction_interval
            ),
        }
    else:
        model = {
            "kind": "ideal_gas",
            "reference_id": phase.model.reference_id,
            "admissible_packing_fraction_interval": None,
        }
    return {
        "identity": phase.identity,
        "fluid_role": phase.fluid_role,
        "amount_role": phase.amount_role,
        "support": support,
        "model": model,
        "start": (
            None
            if phase.start is None
            else {
                "amounts_mol": phase.start.amounts_mol,
                "molar_volume_m3_per_mol": phase.start.molar_volume_m3_per_mol,
            }
        ),
    }


def _phase_supported_component_ids(
    phase: ReactivePhase,
    component_ids: tuple[str, ...],
    charges: tuple[int, ...],
) -> tuple[str, ...]:
    if isinstance(phase.support, AllComponents):
        return component_ids
    if isinstance(phase.support, NeutralComponents):
        return tuple(
            identity for identity, charge in zip(component_ids, charges, strict=True) if charge == 0
        )
    declared = phase.support.component_ids
    unknown = tuple(identity for identity in declared if identity not in component_ids)
    if unknown:
        raise ValueError(f"GREPE phase support names unknown components: {unknown!r}")
    if isinstance(phase.support, DeclaredNeutralComponents):
        charge_by_id = dict(zip(component_ids, charges, strict=True))
        if any(charge_by_id[identity] != 0 for identity in declared):
            raise ValueError("declared neutral support contains a charged component")
    declared_set = frozenset(declared)
    return tuple(identity for identity in component_ids if identity in declared_set)


def _validate_finite_phase_electrolyte_supports(
    model: epcsaft.Mixture,
    phases: tuple[ReactivePhase, ...],
    supported_component_ids: tuple[tuple[str, ...], ...],
) -> None:
    resolved = model._model
    if (
        resolved.topology != "ionic-fully-dissociated-debye-huckel"
        or resolved.relative_permittivity_formulation != "solvent-only"
    ):
        return
    charges = dict(zip(resolved.component_ids, resolved.charge_number, strict=True))
    dielectric_solvents = frozenset(
        component_id
        for component_id, charge, permittivity in zip(
            resolved.component_ids,
            resolved.charge_number,
            resolved.relative_permittivity,
            strict=True,
        )
        if charge == 0 and permittivity is not None
    )
    for phase, supported in zip(phases, supported_component_ids, strict=True):
        if not any(charges[component_id] != 0 for component_id in supported):
            continue
        if not dielectric_solvents.intersection(supported):
            raise epcsaft.ParameterError(
                "invalid-value",
                (
                    ("phase_identity", phase.identity),
                    ("reason", "dielectric-solvent-required"),
                ),
            )


def _feed_payload(feed: ResolvedFeed | AffineFeed | None) -> object:
    if feed is None:
        return None
    if isinstance(feed, AffineFeed):
        return {
            "kind": "affine",
            "component_ids": feed.component_ids,
            "base_amounts_mol": feed.base_amounts_mol,
            "basis_matrix_mol_per_coordinate": feed.basis_matrix_mol_per_coordinate,
            "coordinates": tuple(
                (coordinate.identity, coordinate.bounds) for coordinate in feed.coordinates
            ),
        }
    return {
        "kind": "resolved",
        "component_ids": feed.component_ids,
        "amounts_mol": feed.amounts_mol,
        "coordinate_ids": feed.coordinate_ids,
        "amount_jacobian_mol_per_coordinate": feed.amount_jacobian_mol_per_coordinate,
    }


def _descriptor(
    model: epcsaft.Mixture,
    problem: GeneralReactiveEquilibriumProblem,
    formulation: Literal["finite", "finite_incipient", "temperature_boundary"],
    owner: _FinitePhaseOwnerDescriptor,
    certificate_scope: str,
) -> GeneralReactiveEquilibriumDescriptor:
    outputs = tuple(problem.outputs)
    units_match = all(
        epcsaft.unit_registry.Quantity(1.0, output.unit).to(cast(Any, row).unit).magnitude == 1.0
        for output, row in zip(outputs, owner.rows, strict=True)
    )
    if tuple(output.identity for output in outputs) != owner.row_identities or not units_match:
        raise ValueError("GREPE outputs do not match the compiled owner descriptor")
    payload = {
        "problem_identity": problem.identity,
        "component_ids": model.component_ids,
        "temperature": _state_role_payload(problem.temperature, "kelvin", "temperature"),
        "pressure": _state_role_payload(problem.pressure, "pascal", "pressure"),
        "phases": tuple(_phase_payload(phase) for phase in problem.phases),
        "feed": _feed_payload(problem.feed),
        "phase_reactions": _phase_reaction_payload(problem.phase_reactions),
        "reaction_phase_ids": problem.reaction_phase_ids,
        "coupled_continuation_identity": problem.continuation_identity,
        "continuation_state_fingerprint": (
            None if problem.continuation_state is None else problem.continuation_state.fingerprint
        ),
        "outputs": tuple(
            {
                "identity": output.identity,
                "selector": output.selector,
                "unit": output.unit,
                "basis": output.basis,
                "phase_identity": output.phase_identity,
                "coefficients": output.coefficients,
                "solvent_mass_coefficients_kg_per_mol": (
                    output.solvent_mass_coefficients_kg_per_mol
                ),
                "support": output.support,
                "censor_limit": output.censor_limit,
                "aggregate_identity": output.aggregate_identity,
                "covariance_identity": output.covariance_identity,
            }
            for output in outputs
        ),
        "formulation": formulation,
        "owner_descriptor_fingerprint": owner.fingerprint,
        "certificate_scope": certificate_scope,
    }
    return GeneralReactiveEquilibriumDescriptor(
        "equilibrium.general-reactive",
        2,
        problem.identity,
        model.component_ids,
        tuple(phase.identity for phase in problem.phases),
        outputs,
        owner.active_parameter_identities,
        owner.active_parameter_units,
        owner.parameter_fingerprint,
        owner.active_parameter_selection_fingerprint,
        owner.problem_fingerprint,
        owner.reaction_fingerprint,
        owner.topology_fingerprint,
        owner.reference_fingerprint,
        owner.domain_fingerprint,
        certificate_scope,
        owner.jacobian_shape,
        _fingerprint("epcsaft.grepe-descriptor.v2", payload),
        formulation,
    )


def _unavailable_formulation(
    problem: GeneralReactiveEquilibriumProblem,
    diagnostic: str,
) -> GeneralReactiveEquilibriumResult:
    failure = NonEvaluableTrial("unsupported_topology", diagnostic)
    return GeneralReactiveEquilibriumResult(
        problem=problem,
        descriptor=None,
        status="non_evaluable",
        solver_status="not_started",
        numerical_status="unavailable",
        physical_status="not_adjudicated",
        provider_domain_status="not_adjudicated",
        rows=tuple(
            GeneralEquilibriumRowResult(
                output.identity,
                output.selector,
                output.unit,
                output.basis,
                "non_evaluable",
                None,
                None,
                diagnostic,
            )
            for output in problem.outputs
        ),
        phases=(),
        failure=failure,
        certificate_scope="unavailable_formulation",
        globality_claim="not_established",
        trial=None,
        branch_signature=None,
        evidence=(),
        resolved_policy=(),
    )


@dataclass(frozen=True, slots=True)
class _FinitePhaseOwnerResult:
    descriptor: _FinitePhaseOwnerDescriptor
    rows: tuple[GeneralEquilibriumRowResult, ...]
    phases: tuple[EquilibriumPhaseState, ...]
    trial: ActiveParameterTrialReceipt
    temperature_k: float
    pressure_pa: float
    objective_over_rt: float
    balance_inf_norm: float
    pressure_relative_inf_norm: float
    reaction_affinity_inf_norm: float
    transfer_chemical_potential_inf_norm: float
    kkt_stationarity_inf_norm: float
    source_reference_representation_residual_inf_norm: float
    phase_distinctness: float
    reduced_hessian_min_eigenvalue: float
    declared_direction_rank: int
    admissible_nullity: int
    optimizer_iterations: int
    provider_evaluation_count: int
    branch_signature: tuple[object, ...]


def _reaction_payload(problem: ChemicalEquilibriumProblem) -> dict[str, object]:
    return {
        "species_ids": problem.species_ids,
        "charges": problem.charges,
        "molar_masses_kg_per_mol": problem.molar_masses_kg_per_mol,
        "balance_matrix": problem.balance_matrix,
        "conserved_totals": problem.conserved_totals,
        "reaction_matrix": problem.reaction_matrix,
        "feed_amounts_mol": problem.feed_amounts_mol,
        "equilibrium_constants": tuple(
            {
                "ln_value": record.ln_value,
                "source_id": record.source_id,
                "reference_id": record.reference_id,
                "reaction_orientation": record.reaction_orientation,
                "conversion_id": record.conversion_id,
                "dimensionless": record.dimensionless,
            }
            for record in problem.equilibrium_constants
        ),
        "strict_interior_amount_floor_mol": problem.strict_interior_amount_floor_mol,
        "source_standard_state": _standard_state_payload(problem.source_standard_state),
    }


def _phase_reaction_payload(system: PhaseReactionSystem | None) -> object:
    if system is None:
        return None
    return {
        "phase_ids": system.phase_ids,
        "species_ids": system.species_ids,
        "charges": system.charges,
        "balance_matrix": system.balance_matrix,
        "reactions": tuple(
            {
                "identity": reaction.identity,
                "terms": tuple(
                    (term.phase_id, term.species_id, term.coefficient) for term in reaction.terms
                ),
                "ln_equilibrium_constant": reaction.ln_equilibrium_constant,
                "standard_state": (
                    None
                    if reaction.standard_state is None
                    else {
                        "activity_scale_id": reaction.standard_state.activity_scale_id,
                        "reference_pressure_pa": (reaction.standard_state.reference_pressure_pa),
                    }
                ),
            }
            for reaction in system.reactions
        ),
    }


def _validate_finite_phase_outputs(
    problem: GeneralReactiveEquilibriumProblem,
) -> None:
    reaction = _effective_reaction_system(problem)
    component_count = len(reaction.species_ids)
    phase_by_id = {phase.identity: phase for phase in problem.phases}
    phase_ids = {phase.identity for phase in problem.phases}
    vector_selectors = {
        "phase.amount",
        "phase.fugacity",
        "phase.log_fugacity",
        "phase.molality",
        "phase.mole_fraction",
        "phase.chemical_potential_over_rt",
    }
    scalar_selectors = {
        "phase.compressibility_factor",
        "phase.molar_density",
        "phase.pressure",
    }
    for output in problem.outputs:
        if output.selector == "system.pressure":
            continue
        if output.phase_identity not in phase_ids:
            raise ValueError("finite-phase GREPE output names an undeclared phase")
        if output.selector in vector_selectors:
            if len(output.coefficients) != component_count:
                raise ValueError("finite-phase GREPE vector output requires full coefficients")
            if (
                output.selector in {"phase.amount", "phase.molality", "phase.mole_fraction"}
                and sum(value != 0.0 for value in output.coefficients) > 1
                and (not output.aggregate_identity or not output.covariance_identity)
            ):
                raise ValueError(
                    "finite-phase GREPE aggregate output requires aggregate and covariance identities"
                )
            if output.selector == "phase.molality" and (
                len(output.solvent_mass_coefficients_kg_per_mol) != component_count
                or not any(output.solvent_mass_coefficients_kg_per_mol)
            ):
                raise ValueError("finite-phase molality requires a declared solvent-mass basis")
            if output.selector == "phase.chemical_potential_over_rt":
                phase = phase_by_id[output.phase_identity]
                supported = frozenset(
                    _phase_supported_component_ids(
                        phase, reaction.species_ids, tuple(reaction.charges)
                    )
                )
                if any(
                    coefficient != 0.0 and species not in supported
                    for species, coefficient in zip(
                        reaction.species_ids, output.coefficients, strict=True
                    )
                ):
                    raise ValueError(
                        "finite-phase chemical-potential output names species outside support"
                    )
        elif output.selector in scalar_selectors:
            if output.coefficients:
                raise ValueError("finite-phase GREPE scalar output does not accept coefficients")
        else:
            raise ValueError("output selector is unavailable for finite-phase GREPE")


def _finite_phase_owner_descriptor(
    model: epcsaft.Mixture,
    problem: GeneralReactiveEquilibriumProblem,
    temperature_k: float,
    pressure_pa: float,
    active_parameters: epcsaft.ActiveParameterSet | None,
) -> _FinitePhaseOwnerDescriptor:
    _validate_finite_phase_outputs(problem)
    active_ids = () if active_parameters is None else active_parameters.identities
    active_units = (
        ()
        if active_parameters is None
        else tuple(parameter.unit for parameter in active_parameters.coordinates)
    )
    selection_fingerprint = (
        _fingerprint("epcsaft.empty-active-selection", ())
        if active_parameters is None
        else _fingerprint(
            "epcsaft.active-parameter-selection",
            {
                "parameter_fingerprint": model.parameter_fingerprint,
                "identities": active_ids,
                "units": active_units,
            },
        )
    )
    reaction = _effective_reaction_system(problem)
    reaction_payload = {
        "chemical": _reaction_payload(reaction),
        "phase_directions": _phase_reaction_payload(problem.phase_reactions),
    }
    reaction_fingerprint = _fingerprint("epcsaft.finite-phases-reactions", reaction_payload)
    topology_fingerprint = _fingerprint(
        "epcsaft.finite-phases-topology",
        {
            "component_ids": model.component_ids,
            "phases": tuple(_phase_payload(phase) for phase in problem.phases),
            "reaction_phase_ids": problem.reaction_phase_ids,
            "reaction_matrix": reaction.reaction_matrix,
            "balance_matrix": reaction.balance_matrix,
            "phase_reactions": _phase_reaction_payload(problem.phase_reactions),
        },
    )
    problem_fingerprint = _fingerprint(
        "epcsaft.finite-phases-problem",
        {
            "identity": problem.identity,
            "temperature_k": temperature_k,
            "pressure_pa": pressure_pa,
            "topology_fingerprint": topology_fingerprint,
            "reaction_fingerprint": reaction_fingerprint,
            "outputs": tuple(
                {
                    "identity": output.identity,
                    "selector": output.selector,
                    "unit": output.unit,
                    "basis": output.basis,
                    "phase_identity": output.phase_identity,
                    "coefficients": output.coefficients,
                    "support": output.support,
                    "censor_limit": output.censor_limit,
                    "aggregate_identity": output.aggregate_identity,
                    "covariance_identity": output.covariance_identity,
                }
                for output in problem.outputs
            ),
        },
    )
    domain_fingerprint = _fingerprint(
        "epcsaft.finite-phases-domain",
        {
            "temperature_k": temperature_k,
            "pressure_pa": pressure_pa,
            "packing_intervals": tuple(
                (
                    phase.model.admissible_packing_fraction_interval
                    if isinstance(phase.model, ProviderModel)
                    else None
                )
                for phase in problem.phases
            ),
            "strict_interior_amount_floor_mol": (reaction.strict_interior_amount_floor_mol),
            "local_scope": "fixed_TP_configured_finite_phases_local_minimum",
        },
    )
    reference_fingerprint = _fingerprint(
        "epcsaft.finite-phases-starts",
        {
            "mode": (
                "generated_reaction_feed_partitions_v1"
                if all(phase.start is None for phase in problem.phases)
                else "declared"
            ),
            "phases": tuple(
                {
                    "phase_identity": phase.identity,
                    "start": (
                        None
                        if phase.start is None
                        else {
                            "amounts_mol": phase.start.amounts_mol,
                            "molar_volume_m3_per_mol": (phase.start.molar_volume_m3_per_mol),
                        }
                    ),
                }
                for phase in problem.phases
            ),
        },
    )
    payload = {
        "problem_fingerprint": problem_fingerprint,
        "reaction_fingerprint": reaction_fingerprint,
        "topology_fingerprint": topology_fingerprint,
        "reference_fingerprint": reference_fingerprint,
        "domain_fingerprint": domain_fingerprint,
        "parameter_fingerprint": model.parameter_fingerprint,
        "active_parameter_selection_fingerprint": selection_fingerprint,
        "active_parameter_identities": active_ids,
        "active_parameter_units": active_units,
    }
    return _FinitePhaseOwnerDescriptor(
        problem.outputs,
        tuple(output.identity for output in problem.outputs),
        active_ids,
        active_units,
        model.parameter_fingerprint,
        selection_fingerprint,
        problem_fingerprint,
        reaction_fingerprint,
        topology_fingerprint,
        reference_fingerprint,
        domain_fingerprint,
        (len(problem.outputs), len(active_ids)),
        _fingerprint("epcsaft.finite-phases-descriptor", payload),
    )


def _finite_incipient_owner_descriptor(
    model: epcsaft.Mixture,
    problem: GeneralReactiveEquilibriumProblem,
    temperature_k: float,
    active_parameters: epcsaft.ActiveParameterSet | None,
) -> _FinitePhaseOwnerDescriptor:
    finite = next(phase for phase in problem.phases if phase.amount_role == "finite")
    incipient = next(phase for phase in problem.phases if phase.amount_role == "incipient")
    component_count = len(model.component_ids)
    vapor_count = len(cast(DeclaredNeutralComponents, incipient.support).component_ids)
    for output in problem.outputs:
        quantity = _bubble_quantity(output, finite, incipient)
        count = vapor_count if quantity.startswith("vapor_") else component_count
        if quantity in {"bubble_pressure_pa", "liquid_molar_density_mol_per_m3"}:
            if output.coefficients:
                raise ValueError("scalar finite-incipient output has coefficients")
        elif len(output.coefficients) != count:
            raise ValueError("finite-incipient output coefficients differ from phase support")
        if (
            quantity == "liquid_molality_mol_per_kg"
            and len(output.solvent_mass_coefficients_kg_per_mol) != component_count
        ):
            raise ValueError("finite-incipient molality requires a solvent-mass basis")
    active_ids = () if active_parameters is None else active_parameters.identities
    active_units = (
        ()
        if active_parameters is None
        else tuple(parameter.unit for parameter in active_parameters.coordinates)
    )
    reaction = _effective_reaction_system(problem)
    reaction_fingerprint = _fingerprint("epcsaft.grepe-reactions", _reaction_payload(reaction))
    topology_fingerprint = _fingerprint(
        "epcsaft.grepe-topology",
        {
            "component_ids": model.component_ids,
            "phases": tuple(_phase_payload(phase) for phase in problem.phases),
            "reaction_phase_ids": problem.reaction_phase_ids,
        },
    )
    problem_fingerprint = _fingerprint(
        "epcsaft.grepe-finite-incipient-problem",
        {
            "identity": problem.identity,
            "temperature_k": temperature_k,
            "pressure": _state_role_payload(problem.pressure, "pascal", "pressure"),
            "topology_fingerprint": topology_fingerprint,
            "reaction_fingerprint": reaction_fingerprint,
            "outputs": tuple(asdict(output) for output in problem.outputs),
        },
    )
    selection_fingerprint = _fingerprint(
        "epcsaft.active-parameter-selection",
        (model.parameter_fingerprint, active_ids, active_units),
    )
    reference_fingerprint = (
        _fingerprint("epcsaft.no-continuation", ())
        if problem.continuation_state is None
        else problem.continuation_state.fingerprint
    )
    domain_fingerprint = _fingerprint(
        "epcsaft.grepe-finite-incipient-domain",
        {
            "pressure": _state_role_payload(problem.pressure, "pascal", "pressure"),
            "packing": tuple(
                phase.model.admissible_packing_fraction_interval
                for phase in problem.phases
                if isinstance(phase.model, ProviderModel)
            ),
        },
    )
    payload = (
        problem_fingerprint,
        reaction_fingerprint,
        topology_fingerprint,
        reference_fingerprint,
        domain_fingerprint,
        model.parameter_fingerprint,
        selection_fingerprint,
        active_ids,
        active_units,
    )
    return _FinitePhaseOwnerDescriptor(
        problem.outputs,
        tuple(output.identity for output in problem.outputs),
        active_ids,
        active_units,
        model.parameter_fingerprint,
        selection_fingerprint,
        problem_fingerprint,
        reaction_fingerprint,
        topology_fingerprint,
        reference_fingerprint,
        domain_fingerprint,
        (len(problem.outputs), len(active_ids)),
        _fingerprint("epcsaft.grepe-finite-incipient-descriptor", payload),
    )


def _temperature_boundary_owner_descriptor(
    model: epcsaft.Mixture,
    problem: GeneralReactiveEquilibriumProblem,
    reference_phase_index: int,
    reference_mole_fractions: tuple[float, ...],
    active_parameters: epcsaft.ActiveParameterSet | None,
) -> _FinitePhaseOwnerDescriptor:
    component_count = len(model.component_ids)
    phase_ids = {phase.identity for phase in problem.phases}
    for output in problem.outputs:
        if output.selector in {"system.temperature", "system.pressure"}:
            if output.coefficients:
                raise ValueError("temperature-boundary system outputs have no coefficients")
        elif output.selector == "phase.mole_fraction":
            if (
                output.phase_identity not in phase_ids
                or len(output.coefficients) != component_count
            ):
                raise ValueError(
                    "temperature-boundary composition outputs require full coefficients"
                )
            if sum(value != 0.0 for value in output.coefficients) > 1 and (
                not output.aggregate_identity or not output.covariance_identity
            ):
                raise ValueError("temperature-boundary aggregate composition requires identities")
        elif output.selector in {"phase.molar_density", "phase.pressure"}:
            if output.phase_identity not in phase_ids or output.coefficients:
                raise ValueError("temperature-boundary scalar phase output is invalid")
        else:
            raise ValueError("output selector is unavailable for fixed-pressure coexistence")

    active_ids = () if active_parameters is None else active_parameters.identities
    active_units = (
        ()
        if active_parameters is None
        else tuple(parameter.unit for parameter in active_parameters.coordinates)
    )
    selection_fingerprint = (
        _fingerprint("epcsaft.empty-active-selection", ())
        if active_parameters is None
        else _fingerprint(
            "epcsaft.active-parameter-selection",
            {
                "parameter_fingerprint": model.parameter_fingerprint,
                "identities": active_ids,
                "units": active_units,
            },
        )
    )
    reaction_fingerprint = _fingerprint(
        "epcsaft.temperature-boundary-reactions",
        _reaction_payload(_effective_reaction_system(problem)),
    )
    topology_fingerprint = _fingerprint(
        "epcsaft.temperature-boundary-topology",
        {
            "component_ids": model.component_ids,
            "phases": tuple(_phase_payload(phase) for phase in problem.phases),
            "reference_phase_index": reference_phase_index,
        },
    )
    problem_fingerprint = _fingerprint(
        "epcsaft.temperature-boundary-problem",
        {
            "identity": problem.identity,
            "temperature": _state_role_payload(problem.temperature, "kelvin", "temperature"),
            "pressure": _state_role_payload(problem.pressure, "pascal", "pressure"),
            "topology_fingerprint": topology_fingerprint,
            "reaction_fingerprint": reaction_fingerprint,
            "outputs": tuple(
                (
                    output.identity,
                    output.selector,
                    output.unit,
                    output.basis,
                    output.phase_identity,
                    output.coefficients,
                    output.support,
                    output.censor_limit,
                    output.aggregate_identity,
                    output.covariance_identity,
                )
                for output in problem.outputs
            ),
        },
    )
    reference_fingerprint = _fingerprint(
        "epcsaft.temperature-boundary-reference",
        {
            "feed": _feed_payload(problem.feed),
            "reference_phase_index": reference_phase_index,
            "reference_mole_fractions": reference_mole_fractions,
            "temperature": _state_role_payload(problem.temperature, "kelvin", "temperature"),
        },
    )
    domain_fingerprint = _fingerprint(
        "epcsaft.temperature-boundary-domain",
        {
            "temperature_bounds": _state_role_payload(problem.temperature, "kelvin", "temperature")[
                "bounds"
            ],
            "pressure_pa": _quantity(cast(Fixed, problem.pressure).value, "pascal", "pressure"),
            "packing_intervals": tuple(
                cast(ProviderModel, phase.model).admissible_packing_fraction_interval
                for phase in problem.phases
            ),
        },
    )
    payload = {
        "problem_fingerprint": problem_fingerprint,
        "reaction_fingerprint": reaction_fingerprint,
        "topology_fingerprint": topology_fingerprint,
        "reference_fingerprint": reference_fingerprint,
        "domain_fingerprint": domain_fingerprint,
        "parameter_fingerprint": model.parameter_fingerprint,
        "active_parameter_selection_fingerprint": selection_fingerprint,
    }
    return _FinitePhaseOwnerDescriptor(
        problem.outputs,
        tuple(output.identity for output in problem.outputs),
        active_ids,
        active_units,
        model.parameter_fingerprint,
        selection_fingerprint,
        problem_fingerprint,
        reaction_fingerprint,
        topology_fingerprint,
        reference_fingerprint,
        domain_fingerprint,
        (len(problem.outputs), len(active_ids)),
        _fingerprint("epcsaft.temperature-boundary-descriptor", payload),
    )


def _intensive_boundary_output(boundary: IntensiveBoundary) -> EquilibriumOutput:
    return EquilibriumOutput(
        identity=f"__intensive_boundary__/{boundary.identity}",
        selector=boundary.selector,
        unit="pascal" if boundary.selector == "system.pressure" else "dimensionless",
        basis=(
            "imposed-system-pressure"
            if boundary.selector == "system.pressure"
            else "imposed-dimensionless-intensive-contraction"
        ),
        phase_identity=boundary.phase_identity,
        coefficients=boundary.coefficients,
    )


def _resolved_intensive_problem(
    problem: GeneralReactiveEquilibriumProblem,
    coordinate_value: float,
    outputs: tuple[EquilibriumOutput, ...],
) -> GeneralReactiveEquilibriumProblem:
    feed = cast(AffineFeed, problem.feed)
    boundary = problem.intensive_boundaries[0]
    return replace(
        problem,
        feed=feed.resolve({boundary.coordinate_identity: coordinate_value}),
        outputs=outputs,
        intensive_boundaries=(),
    )


def _describe_general_reactive_equilibrium(
    model: epcsaft.Mixture,
    problem: GeneralReactiveEquilibriumProblem,
    active_parameters: epcsaft.ActiveParameterSet | None,
) -> tuple[GeneralReactiveEquilibriumDescriptor, _FinitePhaseOwnerDescriptor] | str:
    """Validate the public request and preserve its result/persistence identity."""

    if not isinstance(model, epcsaft.Mixture):
        raise TypeError("GREPE currently requires an epcsaft.Mixture")
    if not isinstance(problem, GeneralReactiveEquilibriumProblem):
        raise TypeError("problem must be a GeneralReactiveEquilibriumProblem")

    if problem.intensive_boundaries:
        if not isinstance(problem.temperature, Fixed):
            return "intensive boundaries currently require fixed temperature"
        boundary = problem.intensive_boundaries[0]
        augmented = problem.outputs + (_intensive_boundary_output(boundary),)
        internal = _resolved_intensive_problem(problem, boundary.initial, augmented)
        described = _describe_general_reactive_equilibrium(model, internal, active_parameters)
        if isinstance(described, str):
            return described
        base, owner = described
        if base.formulation not in {"finite", "finite_incipient"}:
            return "intensive boundaries require a differentiable finite or incipient formulation"
        if boundary.selector == "system.pressure" and base.formulation != "finite_incipient":
            return "system pressure boundaries require an incipient formulation"
        if base.formulation == "finite_incipient" and boundary.selector != "system.pressure":
            return "finite-incipient boundaries require exactly one imposed system pressure"
        if base.formulation == "finite_incipient":
            pressure = cast(Solved, problem.pressure)
            lower, upper = (
                _quantity(value, "pascal", "pressure bound") for value in pressure.bounds
            )
            if not lower <= boundary.target <= upper:
                return "imposed system pressure lies outside the declared pressure bounds"
        feed = cast(AffineFeed, problem.feed)
        payload = {
            "schema": "epcsaft.grepe-intensive-boundary.v2",
            "problem_identity": problem.identity,
            "base_descriptor": base.fingerprint,
            "feed": {
                "component_ids": feed.component_ids,
                "base_amounts_mol": feed.base_amounts_mol,
                "basis_matrix_mol_per_coordinate": feed.basis_matrix_mol_per_coordinate,
                "coordinates": tuple(
                    (coordinate.identity, coordinate.bounds) for coordinate in feed.coordinates
                ),
            },
            "boundaries": tuple(
                {
                    "identity": item.identity,
                    "coordinate_identity": item.coordinate_identity,
                    "phase_identity": item.phase_identity,
                    "selector": item.selector,
                    "coefficients": item.coefficients,
                    "target": item.target,
                    "initial": item.initial,
                }
                for item in problem.intensive_boundaries
            ),
            "outputs": tuple(output.identity for output in problem.outputs),
        }
        descriptor = replace(
            base,
            outputs=problem.outputs,
            certificate_scope=(
                "fixed_TP_affine_feed_imposed_intensive_boundaries"
                if base.formulation == "finite"
                else "fixed_T_one_finite_liquid_one_incipient_vapor_affine_feed_imposed_boundaries"
            ),
            jacobian_shape=(len(problem.outputs), len(base.active_parameter_identities)),
            fingerprint=_fingerprint("epcsaft.grepe-intensive-boundary-descriptor.v2", payload),
            formulation="intensive_boundary",
        )
        return descriptor, owner

    _preflight(model, problem)
    finite = tuple(phase for phase in problem.phases if phase.amount_role == "finite")
    incipient = tuple(phase for phase in problem.phases if phase.amount_role == "incipient")

    if (
        isinstance(problem.temperature, Solved)
        and isinstance(problem.pressure, Fixed)
        and len(finite) == len(incipient) == 1
        and not _effective_reaction_system(problem).reaction_matrix
        and not problem.reaction_phase_ids
        and problem.phase_reactions is None
        and all(isinstance(phase.model, ProviderModel) for phase in problem.phases)
        and all(isinstance(phase.support, AllComponents) for phase in problem.phases)
        and {phase.fluid_role for phase in problem.phases} == {"liquid", "vapor"}
    ):
        lower, upper = (
            _quantity(value, "kelvin", "temperature bound") for value in problem.temperature.bounds
        )
        initial = _quantity(problem.temperature.initial, "kelvin", "temperature initial")
        starts = tuple(
            _quantity(value, "kelvin", "temperature start") for value in problem.temperature.starts
        )
        if not lower < initial < upper or any(not lower < value < upper for value in starts):
            raise ValueError("GREPE temperature starts must lie inside their bounds")
        reaction = _effective_reaction_system(problem)
        total = math.fsum(reaction.feed_amounts_mol)
        if total <= 0.0 or any(amount < 0.0 for amount in reaction.feed_amounts_mol):
            raise ValueError("fixed-pressure coexistence requires a nonnegative reference feed")
        reference = tuple(amount / total for amount in reaction.feed_amounts_mol)
        reference_index = problem.phases.index(finite[0])
        owner = _temperature_boundary_owner_descriptor(
            model, problem, reference_index, reference, active_parameters
        )
        return (
            _descriptor(
                model,
                problem,
                "temperature_boundary",
                owner,
                "fixed_P_one_finite_one_incipient_declared_coexistence",
            ),
            owner,
        )

    if not isinstance(problem.temperature, Fixed):
        return "solved temperature has no implemented GREPE formulation"
    temperature_k = _quantity(problem.temperature.value, "kelvin", "temperature")
    if problem.phase_reactions is not None and incipient:
        return "explicit phase/species directions currently require finite configured phases"

    if finite and not incipient and isinstance(problem.pressure, Fixed):
        pressure_pa = _quantity(problem.pressure.value, "pascal", "pressure")
        if all(isinstance(phase.model, ProviderModel) for phase in finite):
            supported = tuple(
                _phase_supported_component_ids(
                    phase,
                    model.component_ids,
                    tuple(_effective_reaction_system(problem).charges),
                )
                for phase in finite
            )
            _validate_finite_phase_electrolyte_supports(model, finite, supported)
        elif not (
            len(finite) == 1
            and isinstance(finite[0].model, IdealGasModel)
            and active_parameters is None
        ):
            return "declared reactive phase/state topology has no implemented GREPE formulation"
        phase_starts = tuple(phase.start for phase in finite)
        if any(start is None for start in phase_starts) and any(
            start is not None for start in phase_starts
        ):
            raise ValueError(
                "configured finite phases must either all supply starts or all omit them"
            )
        if problem.phase_reactions is not None and not phase_starts[0]:
            raise ValueError(
                "explicit phase/species directions require declared reference phase amounts"
            )
        owner = _finite_phase_owner_descriptor(
            model, problem, temperature_k, pressure_pa, active_parameters
        )
        scope = (
            "fixed_TP_single_finite_reacting_phase"
            if len(finite) == 1
            else "fixed_TP_configured_finite_phases_local_minimum"
        )
        return _descriptor(model, problem, "finite", owner, scope), owner

    if (
        len(finite) == len(incipient) == 1
        and isinstance(problem.pressure, Solved)
        and finite[0].fluid_role == "liquid"
        and incipient[0].fluid_role == "vapor"
        and (
            not _effective_reaction_system(problem).reaction_matrix
            or problem.reaction_phase_ids == (finite[0].identity,)
        )
        and isinstance(finite[0].support, AllComponents)
        and isinstance(incipient[0].support, DeclaredNeutralComponents)
    ):
        if not isinstance(finite[0].model, ProviderModel):
            return "the configured finite reacting liquid requires the Provider EOS"
        if finite[0].model.reference_id != "installed-provider-eos":
            return "GREPE provider model reference is not registered"
        lower, upper = (
            _quantity(value, "pascal", "pressure bound") for value in problem.pressure.bounds
        )
        initial = _quantity(problem.pressure.initial, "pascal", "pressure initial")
        if not lower <= initial <= upper:
            raise ValueError("GREPE pressure initial lies outside its bounds")
        owner = _finite_incipient_owner_descriptor(model, problem, temperature_k, active_parameters)
        return (
            _descriptor(
                model,
                problem,
                "finite_incipient",
                owner,
                "fixed_T_one_finite_liquid_one_incipient_vapor",
            ),
            owner,
        )

    return "declared reactive phase/state topology has no implemented GREPE formulation"


def general_reactive_equilibrium_descriptor(
    model: epcsaft.Mixture,
    problem: GeneralReactiveEquilibriumProblem,
    *,
    active_parameters: epcsaft.ActiveParameterSet | None = None,
) -> GeneralReactiveEquilibriumDescriptor:
    """Describe one configured GREPE request without evaluating thermodynamics."""

    described = _describe_general_reactive_equilibrium(model, problem, active_parameters)
    if isinstance(described, str):
        raise ValueError(described)  # noqa: TRY004 - invalid configured formulation
    return described[0]


def _finite_phase_trial(
    descriptor: _FinitePhaseOwnerDescriptor,
    active_parameters: epcsaft.ActiveParameterSet | None,
) -> ActiveParameterTrialReceipt:
    identities = () if active_parameters is None else active_parameters.identities
    values = () if active_parameters is None else active_parameters.values
    units = (
        ()
        if active_parameters is None
        else tuple(parameter.unit for parameter in active_parameters.coordinates)
    )
    active_fingerprint = (
        descriptor.active_parameter_selection_fingerprint
        if active_parameters is None
        else active_parameters.fingerprint
    )
    payload = {
        "parameter_fingerprint": descriptor.parameter_fingerprint,
        "identities": identities,
        "values": values,
        "units": units,
        "active_parameter_fingerprint": active_fingerprint,
    }
    return ActiveParameterTrialReceipt(
        descriptor.parameter_fingerprint,
        identities,
        values,
        units,
        active_fingerprint,
        _fingerprint("epcsaft.active-parameter-trial", payload),
    )


def _finite_phase_failed(
    problem: GeneralReactiveEquilibriumProblem,
    descriptor: GeneralReactiveEquilibriumDescriptor,
    trial: ActiveParameterTrialReceipt,
    code: Literal[
        "physical_rejection",
        "numerical_convergence_failure",
        "first_order_failure",
        "singular_sensitivity",
        "provider_domain_rejection",
        "owner_contract_failure",
    ],
    diagnostic: str,
    *,
    solver_status: str,
    numerical_status: str,
    physical_status: str,
    provider_domain_status: str,
    phases: tuple[EquilibriumPhaseState, ...] = (),
    rows: tuple[GeneralEquilibriumRowResult, ...] | None = None,
    evidence: tuple[tuple[str, float | int | str], ...] = (),
    resolved_policy: tuple[tuple[str, float], ...] = (),
) -> GeneralReactiveEquilibriumResult:
    return GeneralReactiveEquilibriumResult(
        problem=problem,
        descriptor=descriptor,
        status="non_evaluable",
        solver_status=solver_status,
        numerical_status=numerical_status,
        physical_status=physical_status,
        provider_domain_status=provider_domain_status,
        rows=(
            tuple(
                replace(row, status="non_evaluable", jacobian=None, reason=diagnostic)
                for row in rows
            )
            if rows is not None
            else tuple(
            GeneralEquilibriumRowResult(
                output.identity,
                output.selector,
                output.unit,
                output.basis,
                "non_evaluable",
                None,
                None,
                diagnostic,
            )
            for output in problem.outputs
            )
        ),
        phases=phases,
        failure=NonEvaluableTrial(code, diagnostic),
        certificate_scope=descriptor.certificate_scope,
        globality_claim="not_established",
        trial=None,
        branch_signature=None,
        evidence=evidence,
        resolved_policy=resolved_policy,
    )


def _native_grepe_failure(
    problem: GeneralReactiveEquilibriumProblem,
    descriptor: GeneralReactiveEquilibriumDescriptor,
    trial: ActiveParameterTrialReceipt,
    native_result: object,
    native: Any,
    phases: tuple[EquilibriumPhaseState, ...],
) -> GeneralReactiveEquilibriumResult:
    raw_code = str(
        getattr(native_result, "raw_failure_code", "")
        or getattr(native_result, "failure_code", "")
    )
    code = raw_code if raw_code in NON_EVALUABLE_TRIAL_CODES else "owner_contract_failure"
    raw_message = str(
        getattr(native_result, "failure_message", "")
        or native.failure_reason
        or native.solver_status
    )
    diagnostic = (
        f"raw_native_code={raw_code}; {raw_message}"
        if raw_code and code != raw_code
        else ": ".join(part for part in (raw_message, str(native.callback_error)) if part)
    )
    evidence = (
        tuple((str(name), float(value)) for name, value in getattr(native_result, "evidence", ()))
        + tuple((str(name), str(value)) for name, value in getattr(native_result, "text_evidence", ()))
        + (("raw_native_failure_code", raw_code),)
        + (("raw_native_failure_message", raw_message),)
        + ((
            "partial_selection_reason",
            str(getattr(native_result, "partial_selection_reason", "none")),
        ),)
    )
    return _finite_phase_failed(
        problem,
        descriptor,
        trial,
        cast(Any, code),
        diagnostic,
        solver_status=native.solver_status,
        numerical_status=native.numerical_status,
        physical_status=native.physical_status,
        provider_domain_status=native.provider_domain_status,
        phases=phases,
        rows=_finite_phase_rows(
            problem,
            descriptor,
            native,
            aggregate_evaluable=False,
            output_value_available=tuple(getattr(native_result, "output_value_available", ())),
            output_derivative_available=tuple(
                getattr(native_result, "output_derivative_available", ())
            ),
        ),
        evidence=evidence,
        resolved_policy=_native_resolved_policy(native),
    )


def _native_evidence(native_result: object) -> tuple[tuple[str, float | str], ...]:
    return (
        tuple((str(name), float(value)) for name, value in getattr(native_result, "evidence", ()))
        + tuple((str(name), str(value)) for name, value in getattr(native_result, "text_evidence", ()))
    )


def _failure_axes(failure: NonEvaluableTrial) -> tuple[str, str, str]:
    if failure.code == "provider_domain_rejection":
        return "not_adjudicated", "not_adjudicated", "failed"
    if failure.code in {"physical_rejection", "no_incipient_phase", "branch_loss"}:
        return "passed", "failed", "passed"
    return "failed", "not_adjudicated", "not_adjudicated"


def _native_float_vector(value: object) -> tuple[float, ...]:
    return tuple(float(cast(Any, item)) for item in cast(tuple[object, ...], value))


def _native_float_matrix(value: object) -> tuple[tuple[float, ...], ...]:
    return tuple(_native_float_vector(row) for row in cast(tuple[object, ...], value))


def _finite_phase_states(
    model: epcsaft.Mixture,
    problem: GeneralReactiveEquilibriumProblem,
    native: object,
    *,
    mechanical_class: str | tuple[str, ...],
    retained_phase_request_indices: tuple[int, ...] | None = None,
) -> tuple[EquilibriumPhaseState, ...]:
    native_amounts = _native_float_matrix(_native_field(native, "phase_amounts_mol"))
    native_volumes = _native_float_vector(_native_field(native, "phase_volumes_m3"))
    native_pressures = _native_float_vector(_native_field(native, "phase_pressures_pa"))
    native_packing = _native_float_vector(_native_field(native, "phase_packing_fractions"))
    native_mu = _native_float_matrix(_native_field(native, "phase_chemical_potentials_over_rt"))
    phase_count = len(native_amounts)
    if not all(
        len(values) == phase_count
        for values in (
            native_amounts,
            native_volumes,
            native_pressures,
            native_packing,
            native_mu,
        )
    ):
        return ()
    component_count = len(model.component_ids)
    if not (
        all(len(amounts) == component_count and all(math.isfinite(value) for value in amounts)
            and math.fsum(amounts) > 0.0 for amounts in native_amounts)
        and all(len(values) == component_count and all(math.isfinite(value) for value in values)
            for values in native_mu)
        and all(math.isfinite(value) and value > 0.0 for value in native_volumes)
        and all(math.isfinite(value) and value > 0.0 for value in native_pressures)
        and all(math.isfinite(value) for value in native_packing)
    ):
        return ()
    retained = (
        tuple(int(value) for value in retained_phase_request_indices)
        if retained_phase_request_indices is not None
        else tuple(range(phase_count))
    )
    if len(retained) != phase_count or len(set(retained)) != len(retained) or any(
        index < 0 or index >= len(problem.phases) for index in retained
    ):
        return ()
    phases = tuple(problem.phases[index] for index in retained)
    supported_component_ids = tuple(
        _phase_supported_component_ids(
            phase,
            problem.reaction_system.species_ids,
            tuple(problem.reaction_system.charges),
        )
        for phase in phases
    )
    component_index = {identity: index for index, identity in enumerate(model.component_ids)}
    native_states = tuple(
        EquilibriumPhaseState(
            role=phase.fluid_role,
            support=(
                "all_components"
                if isinstance(phase.support, AllComponents)
                else "neutral_components"
                if isinstance(phase.support, NeutralComponents)
                else "declared_neutral_components"
                if isinstance(phase.support, DeclaredNeutralComponents)
                else "declared_components"
            ),
            supported_component_ids=supported,
            mole_fractions=(
                tuple(amount / math.fsum(amounts) for amount in amounts)
                if phase.amount_role == "finite"
                else tuple(
                    amounts[component_index[identity]] / math.fsum(amounts)
                    for identity in supported
                )
            ),
            molar_volume_m3_per_mol=volume / math.fsum(amounts),
            molar_density_mol_m3=math.fsum(amounts) / volume,
            packing_fraction=packing,
            pressure_pa=pressure,
            chemical_potential_component_ids=supported,
            chemical_potential_over_rt=tuple(
                chemical_potentials[component_index[identity]] for identity in supported
            ),
            mechanical_class=(
                mechanical_class[phase_index]
                if isinstance(mechanical_class, tuple)
                else mechanical_class
            ),
        )
        for phase_index, (
            phase,
            supported,
            amounts,
            volume,
            packing,
            pressure,
            chemical_potentials,
        ) in enumerate(
            zip(
                phases,
                supported_component_ids,
                native_amounts,
                native_volumes,
                native_packing,
                native_pressures,
                native_mu,
                strict=True,
            )
        )
    )
    return tuple(
        state for _, state in sorted(zip(retained, native_states, strict=True), key=lambda item: item[0])
    )


def _map_native_incipient_result(
    model: epcsaft.Mixture,
    problem: GeneralReactiveEquilibriumProblem,
    owner_descriptor: _FinitePhaseOwnerDescriptor,
    descriptor: GeneralReactiveEquilibriumDescriptor,
    native_result: object,
    active_parameters: epcsaft.ActiveParameterSet | None,
) -> GeneralReactiveEquilibriumResult:
    native = cast(Any, native_result).equilibrium
    trial = _finite_phase_trial(owner_descriptor, active_parameters)
    mechanical_classes = tuple(cast(Any, native_result).phase_mechanical_classes)
    phases = _finite_phase_states(
        model,
        problem,
        native,
        mechanical_class=mechanical_classes,
        retained_phase_request_indices=tuple(
            int(value) for value in getattr(native_result, "retained_phase_request_indices", ())
        ),
    )
    if not native.accepted:
        return _native_grepe_failure(
            problem, descriptor, trial, native_result, native, phases
        )
    if tuple(native.parameter_order) != descriptor.active_parameter_identities:
        raise ValueError("finite-incipient GREPE parameter derivative order changed")
    rows = _finite_phase_rows(
        problem,
        descriptor,
        native,
        output_value_available=tuple(getattr(native_result, "output_value_available", ())),
        output_derivative_available=tuple(
            getattr(native_result, "output_derivative_available", ())
        ),
    )
    vapor_index = next(
        index for index, phase in enumerate(problem.phases) if phase.amount_role == "incipient"
    )
    vapor = phases[vapor_index]
    liquid_state_identity = _fingerprint(
        "epcsaft.incipient-liquid-state",
        {
            "phase_identity": next(
                phase.identity for phase in problem.phases if phase.amount_role == "finite"
            ),
            "temperature_k": float(cast(Any, native_result).temperature_k),
            "pressure_pa": float(cast(Any, native_result).pressure_pa),
            "amounts_mol": _native_float_matrix(native.phase_amounts_mol)[0],
            "volume_m3": _native_float_vector(native.phase_volumes_m3)[0],
            "trial": trial.fingerprint,
        },
    )
    branch_signature: tuple[object, ...] = (
        "finite_incipient_local_branch",
        liquid_state_identity,
        float(cast(Any, native_result).pressure_pa),
        vapor.mole_fractions,
        vapor.molar_volume_m3_per_mol,
    )
    return GeneralReactiveEquilibriumResult(
        problem,
        descriptor,
        "evaluated",
        native.solver_status,
        native.numerical_status,
        native.physical_status,
        native.provider_domain_status,
        rows,
        phases,
        None,
        descriptor.certificate_scope,
        "not_established",
        trial,
        branch_signature,
        _native_evidence(native_result),
        _native_resolved_policy(native),
    )


def _map_native_boundary_result(
    model: epcsaft.Mixture,
    problem: GeneralReactiveEquilibriumProblem,
    owner_descriptor: _FinitePhaseOwnerDescriptor,
    descriptor: GeneralReactiveEquilibriumDescriptor,
    native_result: object,
    active_parameters: epcsaft.ActiveParameterSet | None,
) -> GeneralReactiveEquilibriumResult:
    native = cast(Any, native_result).equilibrium
    trial = _finite_phase_trial(owner_descriptor, active_parameters)
    phases = _finite_phase_states(
        model,
        problem,
        native,
        mechanical_class=tuple(cast(Any, native_result).phase_mechanical_classes),
        retained_phase_request_indices=tuple(
            int(value) for value in getattr(native_result, "retained_phase_request_indices", ())
        ),
    )
    if not native.accepted:
        return _native_grepe_failure(
            problem, descriptor, trial, native_result, native, phases
        )
    rows = _finite_phase_rows(problem, descriptor, native)
    coordinates = tuple(float(value) for value in cast(Any, native_result).branch_coordinates)
    nested: tuple[object, ...]
    if any(phase.amount_role == "incipient" for phase in problem.phases):
        vapor = phases[next(
            index for index, phase in enumerate(problem.phases) if phase.amount_role == "incipient"
        )]
        nested = (
            "finite_incipient_local_branch",
            descriptor.problem_fingerprint,
            float(cast(Any, native_result).pressure_pa),
            vapor.mole_fractions,
            vapor.molar_volume_m3_per_mol,
        )
    else:
        nested = (
            "finite_phases_local_minimum",
            descriptor.problem_fingerprint,
            tuple(phase.identity for phase in problem.phases),
        )
    native_evidence = _native_evidence(native_result)
    evidence = native_evidence or (
        ("objective_over_rt", native.objective_over_rt),
        ("balance_inf_norm", native.balance_inf_norm),
        ("pressure_relative_inf_norm", native.pressure_relative_inf_norm),
        ("reaction_affinity_inf_norm", native.reaction_affinity_inf_norm),
        ("transfer_chemical_potential_inf_norm", native.transfer_chemical_potential_inf_norm),
        ("kkt_stationarity_inf_norm", native.kkt_stationarity_inf_norm),
        ("reduced_hessian_min_eigenvalue", native.reduced_hessian_min_eigenvalue),
        ("optimizer_iterations", native.optimizer_iterations),
        ("provider_evaluation_count", native.provider_evaluation_count),
    )
    return GeneralReactiveEquilibriumResult(
        problem,
        descriptor,
        "evaluated",
        native.solver_status,
        native.numerical_status,
        native.physical_status,
        native.provider_domain_status,
        rows,
        phases,
        None,
        descriptor.certificate_scope,
        "not_established",
        trial,
        (
            "intensive_boundary",
            ((problem.intensive_boundaries[0].coordinate_identity, coordinates[0]),),
            nested,
        ),
        evidence,
        _native_resolved_policy(native),
    )


def _map_native_temperature_result(
    model: epcsaft.Mixture,
    problem: GeneralReactiveEquilibriumProblem,
    owner_descriptor: _FinitePhaseOwnerDescriptor,
    descriptor: GeneralReactiveEquilibriumDescriptor,
    native_result: object,
    active_parameters: epcsaft.ActiveParameterSet | None,
) -> GeneralReactiveEquilibriumResult:
    native = cast(Any, native_result).equilibrium
    trial = _finite_phase_trial(owner_descriptor, active_parameters)
    phases = _finite_phase_states(
        model,
        problem,
        native,
        mechanical_class=tuple(cast(Any, native_result).phase_mechanical_classes),
    )
    if not native.accepted:
        diagnostic = str(native.failure_reason or native.solver_status)
        code = str(cast(Any, native_result).failure_code or "physical_rejection")
        failure = NonEvaluableTrial(cast(Any, code), diagnostic)
        numerical, physical, provider = _failure_axes(failure)
        return GeneralReactiveEquilibriumResult(
            problem,
            descriptor,
            "non_evaluable",
            native.solver_status,
            numerical,
            physical,
            provider,
            tuple(
                GeneralEquilibriumRowResult(
                    output.identity,
                    output.selector,
                    output.unit,
                    output.basis,
                    "non_evaluable",
                    None,
                    None,
                    diagnostic,
                )
                for output in problem.outputs
            ),
            phases,
            failure,
            descriptor.certificate_scope,
            "not_established",
            trial,
            None,
            tuple((str(name), float(value)) for name, value in cast(Any, native_result).evidence),
            (),
        )
    return GeneralReactiveEquilibriumResult(
        problem,
        descriptor,
        "evaluated",
        native.solver_status,
        native.numerical_status,
        native.physical_status,
        native.provider_domain_status,
        _finite_phase_rows(problem, descriptor, native),
        phases,
        None,
        descriptor.certificate_scope,
        "not_established",
        trial,
        (
            "fixed_p_declared_coexistence",
            int(cast(Any, native_result).branch_coordinates[1]),
            tuple((phase.fluid_role, type(phase.support).__name__) for phase in problem.phases),
        ),
        tuple((str(name), float(value)) for name, value in cast(Any, native_result).evidence),
        (),
    )


def _finite_phase_rows(
    problem: GeneralReactiveEquilibriumProblem,
    descriptor: GeneralReactiveEquilibriumDescriptor,
    native: object,
    input_identities: tuple[str, ...] | None = None,
    *,
    aggregate_evaluable: bool = True,
    output_value_available: tuple[bool, ...] = (),
    output_derivative_available: tuple[bool, ...] = (),
) -> tuple[GeneralEquilibriumRowResult, ...]:
    values = _native_float_vector(_native_field(native, "output_values"))
    input_order = tuple(
        str(value) for value in cast(tuple[object, ...], _native_field(native, "input_order"))
    )
    selected = (
        descriptor.active_parameter_identities if input_identities is None else input_identities
    )
    indices = tuple(input_order.index(identity) for identity in selected if identity in input_order)
    derivatives = _native_float_vector(_native_field(native, "output_input_derivatives"))
    width = len(input_order)
    value_mask_valid = not output_value_available or (
        len(output_value_available) == len(problem.outputs)
        and all(value in (0, 1, False, True) for value in output_value_available)
    )
    derivative_mask_valid = not output_derivative_available or (
        len(output_derivative_available) == len(problem.outputs)
        and all(value in (0, 1, False, True) for value in output_derivative_available)
    )
    value_mask = tuple(bool(value) for value in output_value_available)
    derivative_mask = tuple(bool(value) for value in output_derivative_available)
    values_complete = value_mask_valid and len(values) == len(problem.outputs)
    derivatives_complete = (
        derivative_mask_valid and len(derivatives) == len(problem.outputs) * width
    )
    return tuple(
        _finite_phase_row(
            output, row, values, derivatives, width, indices, selected,
            aggregate_evaluable, values_complete, derivatives_complete,
            value_mask, derivative_mask,
        )
        for row, output in enumerate(problem.outputs)
    )


def _finite_phase_row(
    output: EquilibriumOutput,
    row: int,
    values: tuple[float, ...],
    derivatives: tuple[float, ...],
    width: int,
    indices: tuple[int, ...],
    selected: tuple[str, ...],
    aggregate_evaluable: bool,
    values_complete: bool,
    derivatives_complete: bool,
    value_mask: tuple[bool, ...],
    derivative_mask: tuple[bool, ...],
) -> GeneralEquilibriumRowResult:
    value_available = values_complete and math.isfinite(values[row]) and (
        not value_mask or value_mask[row]
    )
    derivative_available = (
        value_available
        and derivatives_complete
        and (not selected or width)
        and len(indices) == len(selected)
        and (not derivative_mask or derivative_mask[row])
        and all(math.isfinite(derivatives[row * width + index]) for index in indices)
    )
    return GeneralEquilibriumRowResult(
        output.identity,
        output.selector,
        output.unit,
        output.basis,
        "evaluated" if aggregate_evaluable and value_available else "non_evaluable",
        values[row] if value_available else None,
        tuple(derivatives[row * width + index] for index in indices)
        if derivative_available
        else None,
        "" if aggregate_evaluable and value_available else "aggregate_certificate_unavailable",
    )


def _map_native_finite_result(
    model: epcsaft.Mixture,
    problem: GeneralReactiveEquilibriumProblem,
    owner_descriptor: _FinitePhaseOwnerDescriptor,
    descriptor: GeneralReactiveEquilibriumDescriptor,
    native_result: object,
    active_parameters: epcsaft.ActiveParameterSet | None,
    *,
    mechanical_class: str | tuple[str, ...],
) -> GeneralReactiveEquilibriumResult:
    native = cast(Any, native_result).equilibrium
    trial = _finite_phase_trial(owner_descriptor, active_parameters)
    if not native.accepted:
        return _native_grepe_failure(
            problem,
            descriptor,
            trial,
            native_result,
            native,
            _finite_phase_states(
                model,
                problem,
                native,
                mechanical_class="rejected_configured_candidate",
                retained_phase_request_indices=tuple(
                    int(value) for value in getattr(native_result, "retained_phase_request_indices", ())
                ),
            ),
        )
    if tuple(native.parameter_order) != descriptor.active_parameter_identities:
        return _finite_phase_failed(
            problem,
            descriptor,
            trial,
            "first_order_failure",
            "finite-phase GREPE parameter derivative order differs from its descriptor",
            solver_status=native.solver_status,
            numerical_status=native.numerical_status,
            physical_status=native.physical_status,
            provider_domain_status=native.provider_domain_status,
        )
    rows = _finite_phase_rows(
        problem,
        descriptor,
        native,
        output_value_available=tuple(getattr(native_result, "output_value_available", ())),
        output_derivative_available=tuple(
            getattr(native_result, "output_derivative_available", ())
        ),
    )
    phases = _finite_phase_states(
        model,
        problem,
        native,
        mechanical_class=mechanical_class,
        retained_phase_request_indices=tuple(
            int(value) for value in getattr(native_result, "retained_phase_request_indices", ())
        ),
    )
    branch_signature: tuple[object, ...] = (
        "finite_phases_local_minimum",
        descriptor.problem_fingerprint,
        tuple(
            (
                phase.identity,
                problem.continuation_identity,
                None
                if problem.continuation_state is None
                else problem.continuation_state.fingerprint,
            )
            for phase in problem.phases
        ),
    )
    owner = _FinitePhaseOwnerResult(
        owner_descriptor,
        rows,
        phases,
        trial,
        float(cast(Any, native_result).temperature_k),
        float(cast(Any, native_result).pressure_pa),
        native.objective_over_rt,
        native.balance_inf_norm,
        native.pressure_relative_inf_norm,
        native.reaction_affinity_inf_norm,
        native.transfer_chemical_potential_inf_norm,
        native.kkt_stationarity_inf_norm,
        float(cast(Any, native_result).source_reference_representation_residual_inf_norm),
        native.phase_distinctness,
        native.reduced_hessian_min_eigenvalue,
        native.declared_direction_rank,
        native.admissible_nullity,
        native.optimizer_iterations,
        native.provider_evaluation_count,
        branch_signature,
    )
    return GeneralReactiveEquilibriumResult(
        problem=problem,
        descriptor=descriptor,
        status="evaluated",
        solver_status=native.solver_status,
        numerical_status=native.numerical_status,
        physical_status=native.physical_status,
        provider_domain_status=native.provider_domain_status,
        rows=rows,
        phases=phases,
        failure=None,
        certificate_scope=descriptor.certificate_scope,
        globality_claim="not_established",
        trial=trial,
        branch_signature=branch_signature,
        evidence=_result_evidence(owner) + (
            ("local_evaluation_valid", int(native.local_evaluation_valid)),
            ("equilibrium_value_valid", int(native.equilibrium_value_valid)),
            ("local_stability_status", native.local_stability_status),
        ),
        resolved_policy=_native_resolved_policy(native),
    )


def _native_general_reactive_equilibrium_inputs(
    model: epcsaft.Mixture,
    problem: GeneralReactiveEquilibriumProblem,
    *,
    active_parameters: epcsaft.ActiveParameterSet | None = None,
    continuation_start: tuple[tuple[float, ...], float] | None = None,
) -> tuple[_core.NativeModel, GeneralReactiveEquilibriumProblem, str]:
    incipient = next(
        (
            phase
            for phase in problem.phases
            if phase.amount_role == "incipient" and isinstance(phase.model, ProviderModel)
        ),
        None,
    )
    incipient_model = model._native_model(active_parameters)
    if incipient is not None and isinstance(incipient.support, DeclaredNeutralComponents):
        declared = frozenset(incipient.support.component_ids)
        supported = tuple(identity for identity in model.component_ids if identity in declared)
        parameters = model._parameters.select(supported)
        partner = epcsaft.Pure(parameters) if len(supported) == 1 else epcsaft.Mixture(parameters)
        incipient_model = partner._native_model(
            _support_active_parameters(parameters, active_parameters)
        )
    native_problem = problem
    if (
        continuation_start is not None
        and len(problem.phases) == 1
        and isinstance(problem.phases[0].model, IdealGasModel)
    ):
        amounts, volume = continuation_start
        native_problem = replace(
            problem,
            phases=(
                replace(
                    problem.phases[0],
                    start=FinitePhaseStart(amounts, volume / math.fsum(amounts)),
                ),
            ),
        )
    return (
        incipient_model,
        native_problem,
        _source_reference_fingerprint(
            problem.reaction_system.source_standard_state,
            active_parameters_requested=active_parameters is not None,
        ),
    )


def _map_native_general_reactive_equilibrium(
    model: epcsaft.Mixture,
    problem: GeneralReactiveEquilibriumProblem,
    active_parameters: epcsaft.ActiveParameterSet | None,
    native: object,
    described: tuple[GeneralReactiveEquilibriumDescriptor, _FinitePhaseOwnerDescriptor]
    | None = None,
) -> GeneralReactiveEquilibriumResult:
    resolved = (
        _describe_general_reactive_equilibrium(model, problem, active_parameters)
        if described is None
        else described
    )
    if isinstance(resolved, str):
        return _unavailable_formulation(problem, resolved)
    descriptor, owner_descriptor = resolved
    formulation = str(cast(Any, native).formulation)
    if formulation == "finite":
        return _with_continuation_state(
            _map_native_finite_result(
                model,
                problem,
                owner_descriptor,
                descriptor,
                native,
                active_parameters,
                mechanical_class=tuple(cast(Any, native).phase_mechanical_classes),
            )
        )
    if formulation == "finite_incipient":
        return _with_continuation_state(
            _map_native_incipient_result(
                model,
                problem,
                owner_descriptor,
                descriptor,
                native,
                active_parameters,
            )
        )
    if formulation == "intensive_boundary":
        return _with_continuation_state(
            _map_native_boundary_result(
                model,
                problem,
                owner_descriptor,
                descriptor,
                native,
                active_parameters,
            )
        )
    if formulation == "temperature_boundary":
        return _with_continuation_state(
            _map_native_temperature_result(
                model,
                problem,
                owner_descriptor,
                descriptor,
                native,
                active_parameters,
            )
        )
    raise RuntimeError("native GREPE formulation differs from the validated request")


def _solve_general_reactive_equilibrium(
    model: epcsaft.Mixture,
    problem: GeneralReactiveEquilibriumProblem,
    *,
    active_parameters: epcsaft.ActiveParameterSet | None = None,
    continuation_start: tuple[tuple[float, ...], float] | None = None,
) -> GeneralReactiveEquilibriumResult:
    """Validate, enter the native GREPE owner once, and map its final result."""

    described = _describe_general_reactive_equilibrium(model, problem, active_parameters)
    if isinstance(described, str):
        return _unavailable_formulation(problem, described)
    incipient_model, native_problem, reference_fingerprint = (
        _native_general_reactive_equilibrium_inputs(
            model,
            problem,
            active_parameters=active_parameters,
            continuation_start=continuation_start,
        )
    )
    native = _core._solve_grepe(
        model._native,
        incipient_model,
        native_problem,
        active_parameters,
        reference_fingerprint,
    )
    return _map_native_general_reactive_equilibrium(
        model, problem, active_parameters, native, described
    )


def _owner_observation_payload(
    observation: GeneralEquilibriumOwnerObservation,
) -> dict[str, object]:
    return {
        "target_identity": observation.target_identity,
        "prediction_identity": observation.prediction_identity,
        "aggregate_identity": observation.aggregate_identity,
        "covariance_identity": observation.covariance_identity,
    }


def _validated_observation_selection(
    descriptor: GeneralReactiveEquilibriumDescriptor,
    observations: tuple[GeneralEquilibriumOwnerObservation, ...],
) -> tuple[GeneralEquilibriumOwnerObservation, ...]:
    selected = tuple(observations)
    outputs = {output.identity: output for output in descriptor.outputs}
    if not selected or len({item.target_identity for item in selected}) != len(selected):
        raise ValueError("GREPE owner observations must be nonempty and target-unique")
    for observation in selected:
        output = outputs.get(observation.prediction_identity)
        if (
            output is None
            or output.aggregate_identity != observation.aggregate_identity
            or output.covariance_identity != observation.covariance_identity
        ):
            raise ValueError("GREPE owner observation selection is inconsistent")
    return selected


def _general_receipt_payload(
    result: GeneralReactiveEquilibriumResult,
    batch_identity: str,
    observations: tuple[GeneralEquilibriumOwnerObservation, ...],
) -> dict[str, object]:
    descriptor = cast(GeneralReactiveEquilibriumDescriptor, result.descriptor)
    trial = cast(ActiveParameterTrialReceipt, result.trial)
    selected_rows = {row.identity: row for row in result.rows}
    return {
        "schema": "epcsaft.grepe-owner-receipt.v3",
        "batch_identity": batch_identity,
        "descriptor_fingerprint": descriptor.fingerprint,
        "problem_fingerprint": descriptor.problem_fingerprint,
        "reaction_fingerprint": descriptor.reaction_fingerprint,
        "topology_fingerprint": descriptor.topology_fingerprint,
        "reference_fingerprint": descriptor.reference_fingerprint,
        "domain_fingerprint": descriptor.domain_fingerprint,
        "formulation": descriptor.formulation,
        "certificate_scope": descriptor.certificate_scope,
        "globality_claim": result.globality_claim,
        "outcomes": {
            "solver_convergence": {"status": "passed", "reason": result.solver_status},
            "numerical_convergence": {
                "status": "passed",
                "reason": result.numerical_status,
            },
            "physical_certification": {
                "status": "passed",
                "reason": result.physical_status,
            },
            "scientific_fit_acceptance": {
                "status": "not_applicable",
                "reason": "GREPE owner receipt is not a fit acceptance decision",
            },
        },
        "provider_domain_status": result.provider_domain_status,
        "resolved_policy": result.resolved_policy,
        "trial": {
            "parameter_fingerprint": trial.parameter_fingerprint,
            "identities": trial.identities,
            "values": trial.values,
            "units": trial.units,
            "active_parameter_fingerprint": trial.active_parameter_fingerprint,
            "fingerprint": trial.fingerprint,
        },
        "branch_signature": result.branch_signature,
        "evidence": result.evidence,
        "phases": tuple(
            {
                "identity": phase.identity,
                "role": state.role,
                "support": state.support,
                "supported_component_ids": state.supported_component_ids,
                "mole_fractions": state.mole_fractions,
                "molar_volume_m3_per_mol": state.molar_volume_m3_per_mol,
                "molar_density_mol_m3": state.molar_density_mol_m3,
                "packing_fraction": state.packing_fraction,
                "pressure_pa": state.pressure_pa,
                "chemical_potential_component_ids": state.chemical_potential_component_ids,
                "chemical_potential_over_rt": state.chemical_potential_over_rt,
                "mechanical_class": state.mechanical_class,
            }
            for phase, state in zip(result.problem.phases, result.phases)
        ),
        "observations": tuple(
            {
                **_owner_observation_payload(observation),
                "selector": selected_rows[observation.prediction_identity].selector,
                "unit": selected_rows[observation.prediction_identity].unit,
                "value": selected_rows[observation.prediction_identity].value,
                "jacobian": selected_rows[observation.prediction_identity].jacobian,
            }
            for observation in observations
        ),
    }


def serialize_general_reactive_equilibrium_receipt(
    result: GeneralReactiveEquilibriumResult,
    *,
    batch_identity: str,
    observations: tuple[GeneralEquilibriumOwnerObservation, ...],
) -> str:
    """Serialize one canonical opaque GREPE owner receipt."""

    descriptor = result.descriptor
    if result.status != "evaluated" or descriptor is None or result.trial is None:
        raise ValueError("only evaluated GREPE results have owner receipts")
    selected = _validated_observation_selection(descriptor, observations)
    return _canonical_json(_general_receipt_payload(result, batch_identity, selected))


def validate_general_reactive_equilibrium_receipt(
    receipt_json: str,
    *,
    descriptor: GeneralReactiveEquilibriumDescriptor,
    problem: GeneralReactiveEquilibriumProblem,
    batch_identity: str,
    active_parameters: epcsaft.ActiveParameterSet,
    observations: tuple[GeneralEquilibriumOwnerObservation, ...],
) -> str:
    """Validate canonical GREPE owner bytes and return their SHA-256 identity."""

    payload = json.loads(receipt_json)
    if not isinstance(payload, dict) or _canonical_json(payload) != receipt_json:
        raise ValueError("GREPE owner receipt is not canonical JSON")
    selected = _validated_observation_selection(descriptor, observations)
    trial = payload.get("trial")
    rows = payload.get("observations")
    phases = payload.get("phases")
    evidence = payload.get("evidence")
    outcomes = payload.get("outcomes")
    resolved_policy = payload.get("resolved_policy")
    expected_trial = {
        "parameter_fingerprint": descriptor.parameter_fingerprint,
        "identities": active_parameters.identities,
        "values": active_parameters.values,
        "units": tuple(parameter.unit for parameter in active_parameters.coordinates),
        "active_parameter_fingerprint": active_parameters.fingerprint,
    }
    if (
        payload.get("schema") != "epcsaft.grepe-owner-receipt.v3"
        or payload.get("batch_identity") != batch_identity
        or payload.get("descriptor_fingerprint") != descriptor.fingerprint
        or payload.get("problem_fingerprint") != descriptor.problem_fingerprint
        or payload.get("reaction_fingerprint") != descriptor.reaction_fingerprint
        or payload.get("topology_fingerprint") != descriptor.topology_fingerprint
        or payload.get("reference_fingerprint") != descriptor.reference_fingerprint
        or payload.get("domain_fingerprint") != descriptor.domain_fingerprint
        or payload.get("formulation") != descriptor.formulation
        or payload.get("certificate_scope") != descriptor.certificate_scope
        or payload.get("globality_claim") != "not_established"
        or payload.get("provider_domain_status") not in {"passed", "not_adjudicated"}
        or problem.identity != descriptor.problem_identity
        or tuple(phase.identity for phase in problem.phases) != descriptor.phase_identities
        or not isinstance(trial, dict)
        or not isinstance(rows, list)
        or not isinstance(phases, list)
        or not isinstance(evidence, list)
        or not isinstance(outcomes, dict)
        or not isinstance(resolved_policy, list)
        or len(rows) != len(selected)
    ):
        raise ValueError("GREPE owner receipt identity or observations do not match")
    if (
        trial.get("parameter_fingerprint") != expected_trial["parameter_fingerprint"]
        or tuple(trial.get("identities", ())) != active_parameters.identities
        or tuple(float(value) for value in trial.get("values", ())) != active_parameters.values
        or tuple(trial.get("units", ())) != expected_trial["units"]
        or trial.get("active_parameter_fingerprint") != active_parameters.fingerprint
        or trial.get("fingerprint")
        != _fingerprint("epcsaft.active-parameter-trial", expected_trial)
    ):
        raise ValueError("GREPE owner receipt trial does not match")
    expected_outcome_statuses = {
        "solver_convergence": "passed",
        "numerical_convergence": "passed",
        "physical_certification": "passed",
        "scientific_fit_acceptance": "not_applicable",
    }
    if set(outcomes) != set(expected_outcome_statuses) or any(
        not isinstance(outcome, dict)
        or outcome.get("status") != expected_outcome_statuses[name]
        or not isinstance(outcome.get("reason"), str)
        or not outcome["reason"]
        for name, outcome in outcomes.items()
    ):
        raise ValueError("GREPE owner receipt outcomes are invalid")
    policy_records = tuple(
        (str(record[0]), float(record[1]))
        for record in resolved_policy
        if isinstance(record, list) and len(record) == 2
    )
    if len(policy_records) != len(resolved_policy) or (
        policy_records
        and (
            len({name for name, _ in policy_records}) != len(policy_records)
            or any(
                not name or not math.isfinite(value) or value < 0.0
                for name, value in policy_records
            )
        )
    ):
        raise ValueError("GREPE owner receipt numerical policy is invalid")
    output_by_identity = {output.identity: output for output in descriptor.outputs}
    for record, observation in zip(rows, selected, strict=True):
        output = output_by_identity[observation.prediction_identity]
        if (
            not isinstance(record, dict)
            or record.get("target_identity") != observation.target_identity
            or record.get("prediction_identity") != observation.prediction_identity
            or record.get("aggregate_identity") != observation.aggregate_identity
            or record.get("covariance_identity") != observation.covariance_identity
            or record.get("selector") != output.selector
            or record.get("unit") != output.unit
            or not isinstance(record.get("value"), (int, float))
            or not math.isfinite(float(record["value"]))
            or not isinstance(record.get("jacobian"), list)
            or len(cast(list[object], record["jacobian"])) != len(active_parameters.identities)
            or any(
                not isinstance(value, (int, float)) or not math.isfinite(float(value))
                for value in cast(list[object], record["jacobian"])
            )
        ):
            raise ValueError("GREPE owner receipt observation does not match")
    if phases:
        for record, phase in zip(phases, problem.phases, strict=True):
            if (
                not isinstance(record, dict)
                or record.get("identity") != phase.identity
                or record.get("role") != phase.fluid_role
                or not isinstance(record.get("mole_fractions"), list)
                or not math.isclose(
                    math.fsum(float(value) for value in record["mole_fractions"]),
                    1.0,
                    rel_tol=0.0,
                    abs_tol=2.0e-12,
                )
            ):
                raise ValueError("GREPE owner receipt phase does not match")
    if len(phases) not in {0, len(problem.phases)} or any(
        not isinstance(item, list)
        or len(item) != 2
        or not isinstance(item[0], str)
        or not isinstance(item[1], (int, float, str))
        for item in evidence
    ):
        raise ValueError("GREPE owner receipt evidence is invalid")
    return f"sha256:{hashlib.sha256(receipt_json.encode('utf-8')).hexdigest()}"


__all__ = (
    "EquilibriumOutput",
    "FinitePhaseStart",
    "Fixed",
    "GeneralEquilibriumOwnerObservation",
    "GeneralEquilibriumRowResult",
    "GeneralReactiveEquilibriumDescriptor",
    "GeneralReactiveEquilibriumProblem",
    "GeneralReactiveEquilibriumResult",
    "IdealGasModel",
    "IntensiveBoundary",
    "ProviderModel",
    "ReactivePhase",
    "Solved",
    "general_reactive_equilibrium_descriptor",
    "serialize_general_reactive_equilibrium_receipt",
    "validate_general_reactive_equilibrium_receipt",
)
