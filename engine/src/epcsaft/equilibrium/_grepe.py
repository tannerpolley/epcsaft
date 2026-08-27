from __future__ import annotations

import hashlib
import json
import math
from collections.abc import Iterator
from dataclasses import asdict, dataclass, replace
from typing import Any, Literal, cast

from pint import Quantity

import epcsaft
from epcsaft.eos import _native_transport

from .._trials import NonEvaluableTrial
from . import _equilibrium
from ._api import (
    ActiveParameterTrialReceipt,
    AllComponents,
    ChemicalEquilibriumError,
    ChemicalEquilibriumProblem,
    ChemicalEquilibriumSensitivityRequest,
    ChemicalStandardState,
    DeclaredComponents,
    DeclaredNeutralComponents,
    DeclaredTopology,
    EquilibriumContinuationState,
    EquilibriumPhaseState,
    FixedPCoexistence,
    IdealGasPhase,
    IncipientTemperatureBoundary,
    NeutralComponents,
    PhaseEquilibriumError,
    PhaseEquilibriumProblem,
    PhaseSpecification,
    PhaseSupport,
    ProviderPhase,
    _chemical_equilibrium_impl,
    _continuation_start,
    _contract_phase_output,
    _ScalarBoundarySolveError,
    _solve_bracketed_scalar,
    _source_reference_payload,
    _source_reference_transfer,
    solve,
)
from ._feed import AffineFeed, ResolvedFeed
from ._incipient import (
    _closure,
    _evaluate_incipient,
    _FiniteIncipientKernel,
    _fixed_pressure_coordinate_derivatives,
    _fixed_pressure_coordinate_result,
    _IdealVapor,
    _IncipientFailure,
    _ProviderVapor,
)
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
_FINITE_PHASE_FAILURE_RESIDUALS = (
    "balance_inf_norm",
    "pressure_relative_inf_norm",
    "reaction_affinity_inf_norm",
    "transfer_chemical_potential_inf_norm",
    "kkt_stationarity_inf_norm",
)


def _finite_phase_failure_score(
    item: tuple[int, dict[str, object]],
) -> tuple[float, int]:
    score = 0.0
    for name in _FINITE_PHASE_FAILURE_RESIDUALS:
        raw_value = item[1][name]
        if raw_value is None:
            score += 1.0e100
            continue
        value = float(cast(float, raw_value))
        score += abs(value) if math.isfinite(value) else 1.0e100
    return score, item[0]


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


def _native_resolved_policy(native: dict[str, object]) -> tuple[tuple[str, float], ...]:
    records = tuple(
        (str(record[0]), float(cast(float, record[1])))
        for record in cast(tuple[tuple[object, object], ...], native["resolved_policy"])
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


def _phase_model(
    model: epcsaft.Mixture,
    phase: ReactivePhase,
) -> ProviderPhase | IdealGasPhase:
    if isinstance(phase.model, ProviderModel):
        if phase.model.reference_id != "installed-provider-eos":
            raise ValueError("GREPE provider model reference is not registered")
        return ProviderPhase(model, phase.model.admissible_packing_fraction_interval)
    return IdealGasPhase(model.parameter_fingerprint, phase.model.reference_id)


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


def _canonical_json(value: object) -> str:
    return json.dumps(
        value,
        ensure_ascii=False,
        allow_nan=False,
        separators=(",", ":"),
        sort_keys=True,
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


def _incipient_branch_signature(closure: Any) -> tuple[object, ...]:
    return (
        "finite_incipient_local_branch",
        closure.liquid_state_identity,
        closure.pressure_pa,
        closure.phases[1].mole_fractions,
        closure.phases[1].molar_volume_m3_per_mol,
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
class _IncipientCompilation:
    temperature: Quantity[Any]
    owner_problem: _FiniteIncipientKernel
    owner_descriptor: _FinitePhaseOwnerDescriptor
    descriptor: GeneralReactiveEquilibriumDescriptor


@dataclass(frozen=True, slots=True)
class _FinitePhaseCompilation:
    temperature_k: float
    pressure_pa: float
    phases: tuple[ReactivePhase, ...]
    supported_component_ids: tuple[tuple[str, ...], ...]
    declared_start_set: tuple[FinitePhaseStart, ...] | None
    generated_allocations: tuple[tuple[tuple[float, ...], ...], ...]
    owner_descriptor: _FinitePhaseOwnerDescriptor
    descriptor: GeneralReactiveEquilibriumDescriptor


@dataclass(frozen=True, slots=True)
class _IdealPhaseCompilation:
    temperature_k: float
    pressure_pa: float
    phase: ReactivePhase
    owner_descriptor: _FinitePhaseOwnerDescriptor
    descriptor: GeneralReactiveEquilibriumDescriptor


@dataclass(frozen=True, slots=True)
class _TemperatureBoundaryCompilation:
    owner_problem: PhaseEquilibriumProblem
    owner_descriptor: _FinitePhaseOwnerDescriptor
    descriptor: GeneralReactiveEquilibriumDescriptor


@dataclass(frozen=True, slots=True)
class _IntensiveBoundaryCompilation:
    descriptor: GeneralReactiveEquilibriumDescriptor
    augmented_outputs: tuple[EquilibriumOutput, ...]
    base_formulation: Literal["finite", "finite_incipient"]


type _IntensiveOwnerResult = GeneralReactiveEquilibriumResult


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


def _generated_finite_phase_allocations(
    problem: GeneralReactiveEquilibriumProblem,
    phases: tuple[ReactivePhase, ...],
) -> tuple[tuple[tuple[float, ...], ...], ...]:
    """Build deterministic feasible amount partitions for declared phases."""

    reaction_system = _effective_reaction_system(problem)
    feed = reaction_system.feed_amounts_mol
    phase_count = len(phases)
    component_count = len(feed)
    floor = reaction_system.strict_interior_amount_floor_mol
    if any(amount <= phase_count * 0.1 * floor for amount in feed):
        raise ValueError("automatic finite-phase starts require a strict-interior species feed")

    allocations: list[tuple[tuple[float, ...], ...]] = []
    if phase_count == 2:
        neutral_indices = tuple(
            index for index, charge in enumerate(reaction_system.charges) if charge == 0
        )
        for reaction in reaction_system.reaction_matrix:
            first = tuple(
                (
                    0.5
                    if reaction_system.charges[component] != 0 or coefficient == 0.0
                    else 0.8
                    if coefficient > 0.0
                    else 0.2
                )
                for component, coefficient in enumerate(reaction)
            )
            if any(value == 0.8 for value in first) and any(value == 0.2 for value in first):
                allocations.append((first, tuple(1.0 - value for value in first)))
        for anchor in neutral_indices:
            first = tuple(
                0.8 if component == anchor else 0.2 for component in range(component_count)
            )
            if len(neutral_indices) != component_count:
                first = tuple(
                    value if reaction_system.charges[component] == 0 else 0.5
                    for component, value in enumerate(first)
                )
            anchor_allocation = (first, tuple(1.0 - value for value in first))
            if anchor_allocation not in allocations:
                allocations.append(anchor_allocation)
    equal_allocation = tuple((1.0 / phase_count,) * component_count for _ in range(phase_count))
    if equal_allocation not in allocations:
        allocations.append(equal_allocation)

    component_index = {
        identity: index for index, identity in enumerate(reaction_system.species_ids)
    }
    phase_supports = tuple(
        frozenset(
            component_index[identity]
            for identity in _phase_supported_component_ids(
                phase, reaction_system.species_ids, reaction_system.charges
            )
        )
        for phase in phases
    )
    projected: list[tuple[tuple[float, ...], ...]] = []
    for allocation in allocations:
        rows = [[0.0] * component_count for _ in phases]
        for component in range(component_count):
            eligible = tuple(
                phase for phase, support in enumerate(phase_supports) if component in support
            )
            if not eligible:
                raise ValueError(
                    "every fed species must belong to at least one declared finite-phase support"
                )
            scale = math.fsum(allocation[phase][component] for phase in eligible)
            for phase in eligible:
                rows[phase][component] = allocation[phase][component] / scale
        candidate = tuple(tuple(row) for row in rows)
        if candidate not in projected:
            projected.append(candidate)
    return tuple(projected)


def _finite_phase_start_set_from_allocation(
    model: epcsaft.Mixture,
    phases: tuple[ReactivePhase, ...],
    feed: tuple[float, ...],
    allocation: tuple[tuple[float, ...], ...],
    temperature_k: float,
    pressure_pa: float,
) -> tuple[FinitePhaseStart, ...]:
    temperature = temperature_k * epcsaft.unit_registry.kelvin
    pressure = pressure_pa * epcsaft.unit_registry.pascal
    starts: list[FinitePhaseStart] = []
    for phase, fractions in zip(phases, allocation, strict=True):
        amounts = tuple(amount * fraction for amount, fraction in zip(feed, fractions, strict=True))
        total = math.fsum(amounts)
        composition = tuple(amount / total for amount in amounts)
        state = model.state(
            T=temperature,
            P=pressure,
            x=composition,
            phase=phase.fluid_role,
        )
        density = float(state.molar_density.to("mole / meter**3").magnitude)
        starts.append(FinitePhaseStart(amounts, 1.0 / density))
    return tuple(starts)


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


def _compile_intensive_boundary(
    model: epcsaft.Mixture,
    problem: GeneralReactiveEquilibriumProblem,
    active_parameters: epcsaft.ActiveParameterSet | None,
) -> _IntensiveBoundaryCompilation | str:
    if not isinstance(problem.temperature, Fixed):
        return "intensive boundaries currently require fixed temperature"
    boundary = problem.intensive_boundaries[0]
    augmented_outputs = problem.outputs + (_intensive_boundary_output(boundary),)
    internal_problem = _resolved_intensive_problem(problem, boundary.initial, augmented_outputs)
    compiled = _compile_intensive_owner(
        model,
        internal_problem,
        active_parameters=active_parameters,
    )
    if isinstance(compiled, str):
        return compiled
    if not isinstance(
        compiled, (_FinitePhaseCompilation, _IdealPhaseCompilation, _IncipientCompilation)
    ):
        return "intensive boundaries require a differentiable finite or incipient formulation"
    pressure_boundary = boundary.selector == "system.pressure"
    if (
        isinstance(compiled, (_FinitePhaseCompilation, _IdealPhaseCompilation))
        and pressure_boundary
    ):
        return "system pressure boundaries require an incipient formulation"
    if isinstance(compiled, _IncipientCompilation) and not pressure_boundary:
        return "finite-incipient boundaries require exactly one imposed system pressure"
    if isinstance(compiled, _IncipientCompilation):
        pressure = cast(Solved, problem.pressure)
        lower, upper = (_quantity(value, "pascal", "pressure bound") for value in pressure.bounds)
        if not lower <= boundary.target <= upper:
            return "imposed system pressure lies outside the declared pressure bounds"
    feed = cast(AffineFeed, problem.feed)
    payload = {
        "schema": "epcsaft.grepe-intensive-boundary.v2",
        "problem_identity": problem.identity,
        "base_descriptor": compiled.descriptor.fingerprint,
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
                "identity": boundary.identity,
                "coordinate_identity": boundary.coordinate_identity,
                "phase_identity": boundary.phase_identity,
                "selector": boundary.selector,
                "coefficients": boundary.coefficients,
                "target": boundary.target,
                "initial": boundary.initial,
            }
            for boundary in problem.intensive_boundaries
        ),
        "outputs": tuple(output.identity for output in problem.outputs),
    }
    descriptor = replace(
        compiled.descriptor,
        outputs=problem.outputs,
        certificate_scope=(
            "fixed_TP_affine_feed_imposed_intensive_boundaries"
            if isinstance(compiled, (_FinitePhaseCompilation, _IdealPhaseCompilation))
            else "fixed_T_one_finite_liquid_one_incipient_vapor_affine_feed_imposed_boundaries"
        ),
        jacobian_shape=(
            len(problem.outputs),
            len(compiled.descriptor.active_parameter_identities),
        ),
        fingerprint=_fingerprint("epcsaft.grepe-intensive-boundary-descriptor.v2", payload),
        formulation="intensive_boundary",
    )
    return _IntensiveBoundaryCompilation(
        descriptor,
        augmented_outputs,
        cast(Literal["finite", "finite_incipient"], compiled.descriptor.formulation),
    )


def _compile_intensive_owner(
    model: epcsaft.Mixture,
    problem: GeneralReactiveEquilibriumProblem,
    *,
    active_parameters: epcsaft.ActiveParameterSet | None,
) -> (
    _IncipientCompilation
    | _FinitePhaseCompilation
    | _IdealPhaseCompilation
    | _TemperatureBoundaryCompilation
    | _IntensiveBoundaryCompilation
    | str
):
    compiled = _compile_general_reactive_equilibrium(
        model, problem, active_parameters=active_parameters
    )
    return compiled


def _compile_general_reactive_equilibrium(
    model: epcsaft.Mixture,
    problem: GeneralReactiveEquilibriumProblem,
    *,
    active_parameters: epcsaft.ActiveParameterSet | None = None,
) -> (
    _IncipientCompilation
    | _FinitePhaseCompilation
    | _IdealPhaseCompilation
    | _TemperatureBoundaryCompilation
    | _IntensiveBoundaryCompilation
    | str
):
    if not isinstance(model, epcsaft.Mixture):
        raise TypeError("GREPE currently requires an epcsaft.Mixture")
    if not isinstance(problem, GeneralReactiveEquilibriumProblem):
        raise TypeError("problem must be a GeneralReactiveEquilibriumProblem")
    if problem.intensive_boundaries:
        return _compile_intensive_boundary(model, problem, active_parameters)
    _preflight(model, problem)
    finite = tuple(phase for phase in problem.phases if phase.amount_role == "finite")
    incipient = tuple(phase for phase in problem.phases if phase.amount_role == "incipient")

    if (
        isinstance(problem.temperature, Solved)
        and isinstance(problem.pressure, Fixed)
        and len(finite) == 1
        and len(incipient) == 1
        and not _effective_reaction_system(problem).reaction_matrix
        and not problem.reaction_phase_ids
        and problem.phase_reactions is None
        and all(isinstance(phase.model, ProviderModel) for phase in problem.phases)
        and all(isinstance(phase.support, AllComponents) for phase in problem.phases)
        and {phase.fluid_role for phase in problem.phases} == {"liquid", "vapor"}
    ):
        bounds = tuple(value.to("kelvin") for value in problem.temperature.bounds)
        initial = problem.temperature.initial.to("kelvin")
        starts = tuple(value.to("kelvin") for value in problem.temperature.starts)
        lower, upper = (float(value.magnitude) for value in bounds)
        initial_k = float(initial.magnitude)
        if not lower < initial_k < upper or any(
            not lower < float(value.magnitude) < upper for value in starts
        ):
            raise ValueError("GREPE temperature starts must lie inside their bounds")
        reaction = _effective_reaction_system(problem)
        total = math.fsum(reaction.feed_amounts_mol)
        if total <= 0.0 or any(amount < 0.0 for amount in reaction.feed_amounts_mol):
            raise ValueError("fixed-pressure coexistence requires a nonnegative reference feed")
        reference = tuple(amount / total for amount in reaction.feed_amounts_mol)
        reference_index = problem.phases.index(finite[0])
        owner_problem = PhaseEquilibriumProblem(
            FixedPCoexistence(
                problem.pressure.value,
                cast(tuple[Quantity[Any], Quantity[Any]], bounds),
                DeclaredTopology(
                    tuple(
                        PhaseSpecification(phase.fluid_role, cast(AllComponents, phase.support))
                        for phase in problem.phases
                    )
                ),
                reference_index,
                reference,
                initial,
                tuple(
                    cast(ProviderModel, phase.model).admissible_packing_fraction_interval
                    for phase in problem.phases
                ),
            )
        )
        temperature_owner_descriptor = _temperature_boundary_owner_descriptor(
            model,
            problem,
            reference_index,
            reference,
            active_parameters,
        )
        return _TemperatureBoundaryCompilation(
            owner_problem,
            temperature_owner_descriptor,
            _descriptor(
                model,
                problem,
                "temperature_boundary",
                temperature_owner_descriptor,
                "fixed_P_one_finite_one_incipient_declared_coexistence",
            ),
        )

    if not isinstance(problem.temperature, Fixed):
        return "solved temperature has no implemented GREPE formulation"
    temperature_k = _quantity(problem.temperature.value, "kelvin", "temperature")
    temperature = temperature_k * epcsaft.unit_registry.kelvin

    if problem.phase_reactions is not None and incipient:
        return "explicit phase/species directions currently require finite configured phases"

    if (
        finite
        and not incipient
        and isinstance(problem.pressure, Fixed)
        and all(isinstance(phase.model, ProviderModel) for phase in finite)
    ):
        pressure_pa = _quantity(problem.pressure.value, "pascal", "pressure")
        supported_component_ids = tuple(
            _phase_supported_component_ids(
                phase,
                model.component_ids,
                tuple(_effective_reaction_system(problem).charges),
            )
            for phase in finite
        )
        _validate_finite_phase_electrolyte_supports(model, finite, supported_component_ids)
        supplied_starts = tuple(phase.start for phase in finite)
        if any(start is None for start in supplied_starts) and any(
            start is not None for start in supplied_starts
        ):
            raise ValueError(
                "configured finite phases must either all supply starts or all omit them"
            )
        component_count = len(model.component_ids)
        if supplied_starts[0] is not None and any(
            len(cast(FinitePhaseStart, start).amounts_mol) != component_count
            for start in supplied_starts
        ):
            raise ValueError(
                "configured finite-phase starts must match the Provider component order"
            )
        component_index = {identity: index for index, identity in enumerate(model.component_ids)}
        for phase, supported, start in zip(
            finite, supported_component_ids, supplied_starts, strict=True
        ):
            if start is None:
                continue
            support_indices = {component_index[identity] for identity in supported}
            for index, amount in enumerate(start.amounts_mol):
                if index in support_indices and amount <= 0.0:
                    raise ValueError(
                        f"finite phase {phase.identity!r} start must be positive on its support"
                    )
                if index not in support_indices and amount != 0.0:
                    raise ValueError(
                        f"finite phase {phase.identity!r} start must be zero outside its support"
                    )
        declared_start_set = (
            tuple(cast(FinitePhaseStart, start) for start in supplied_starts)
            if supplied_starts[0] is not None
            else None
        )
        if problem.continuation_state is not None:
            if declared_start_set is not None or len(finite) != 1:
                raise ValueError(
                    "a continuation state replaces the explicit start of one finite phase"
                )
            continued = _continuation_start(
                _effective_reaction_system(problem),
                problem.continuation_state.phase(finite[0].identity),
            )
            if continued is None:
                raise ValueError("the finite continuation state is incomplete")
            amounts, volume = continued
            declared_start_set = (FinitePhaseStart(amounts, volume / math.fsum(amounts)),)
        if problem.phase_reactions is not None and declared_start_set is None:
            raise ValueError(
                "explicit phase/species directions require declared reference phase amounts"
            )
        reaction = _effective_reaction_system(problem)
        generated_allocations = (
            ()
            if declared_start_set is not None
            and any(
                amount <= len(finite) * 0.1 * reaction.strict_interior_amount_floor_mol
                for amount in reaction.feed_amounts_mol
            )
            else _generated_finite_phase_allocations(problem, finite)
        )
        owner_descriptor = _finite_phase_owner_descriptor(
            model,
            problem,
            temperature_k,
            pressure_pa,
            active_parameters,
        )
        certificate_scope = (
            "fixed_TP_single_finite_reacting_phase"
            if len(finite) == 1
            else "fixed_TP_configured_finite_phases_local_minimum"
        )
        descriptor = _descriptor(
            model,
            problem,
            "finite",
            owner_descriptor,
            certificate_scope,
        )
        return _FinitePhaseCompilation(
            temperature_k,
            pressure_pa,
            finite,
            supported_component_ids,
            declared_start_set,
            generated_allocations,
            owner_descriptor,
            descriptor,
        )

    if (
        len(finite) == 1
        and not incipient
        and isinstance(problem.pressure, Fixed)
        and isinstance(finite[0].model, IdealGasModel)
    ):
        if active_parameters is not None:
            return "ideal-gas GREPE has no Provider parameter coordinates"
        pressure_pa = _quantity(problem.pressure.value, "pascal", "pressure")
        owner_descriptor = _finite_phase_owner_descriptor(
            model, problem, temperature_k, pressure_pa, None
        )
        return _IdealPhaseCompilation(
            temperature_k,
            pressure_pa,
            finite[0],
            owner_descriptor,
            _descriptor(
                model,
                problem,
                "finite",
                owner_descriptor,
                "fixed_TP_single_finite_reacting_phase",
            ),
        )

    if (
        len(finite) == 1
        and len(incipient) == 1
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
        liquid, vapor = finite[0], incipient[0]
        liquid_model = _phase_model(model, liquid)
        if not isinstance(liquid_model, ProviderPhase):
            return "the configured finite reacting liquid requires the Provider EOS"
        incipient_problem = _FiniteIncipientKernel(
            liquid.identity,
            liquid_model,
            _effective_reaction_system(problem),
            vapor.identity,
            cast(DeclaredNeutralComponents, vapor.support).component_ids,
            _ProviderVapor(vapor.model.reference_id)
            if isinstance(vapor.model, ProviderModel)
            else _IdealVapor(vapor.model.reference_id),
            cast(
                tuple[float, float],
                tuple(
                    _quantity(value, "pascal", "pressure bound")
                    for value in problem.pressure.bounds
                ),
            ),
            tuple(
                _quantity(value, "pascal", "pressure start") for value in problem.pressure.starts
            ),
            problem.continuation_state,
        )
        initial_pressure = _quantity(problem.pressure.initial, "pascal", "pressure initial")
        if (
            not incipient_problem.pressure_interval_pa[0]
            <= initial_pressure
            <= incipient_problem.pressure_interval_pa[1]
        ):
            raise ValueError("GREPE pressure initial lies outside its bounds")
        bubble_descriptor = _finite_incipient_owner_descriptor(
            model,
            problem,
            temperature_k,
            active_parameters,
        )
        descriptor = _descriptor(
            model,
            problem,
            "finite_incipient",
            bubble_descriptor,
            "fixed_T_one_finite_liquid_one_incipient_vapor",
        )
        return _IncipientCompilation(
            temperature,
            incipient_problem,
            bubble_descriptor,
            descriptor,
        )

    return "declared reactive phase/state topology has no implemented GREPE formulation"


def general_reactive_equilibrium_descriptor(
    model: epcsaft.Mixture,
    problem: GeneralReactiveEquilibriumProblem,
    *,
    active_parameters: epcsaft.ActiveParameterSet | None = None,
) -> GeneralReactiveEquilibriumDescriptor:
    """Describe one configured GREPE request without evaluating thermodynamics."""

    compiled = _compile_general_reactive_equilibrium(
        model, problem, active_parameters=active_parameters
    )
    if isinstance(compiled, str):
        raise ValueError(compiled)  # noqa: TRY004 - invalid configured formulation
    return compiled.descriptor


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


def _finite_phase_spec(
    problem: GeneralReactiveEquilibriumProblem,
    temperature_k: float,
    pressure_pa: float,
) -> dict[str, object]:
    reaction = _effective_reaction_system(problem)
    standard = reaction.source_standard_state
    return {
        "species_ids": reaction.species_ids,
        "charges": reaction.charges,
        "molar_masses_kg_per_mol": reaction.molar_masses_kg_per_mol,
        "balance_matrix": reaction.balance_matrix,
        "conserved_totals": reaction.conserved_totals,
        "reaction_matrix": reaction.reaction_matrix,
        "feed_amounts": reaction.feed_amounts_mol,
        "ln_k": tuple(record.ln_value for record in reaction.equilibrium_constants),
        "equilibrium_constant_records": tuple(
            {
                "source_id": record.source_id,
                "reference_id": record.reference_id,
                "reaction_orientation": record.reaction_orientation,
                "conversion_id": record.conversion_id,
                "dimensionless": record.dimensionless,
                "temperature_k": temperature_k,
                "pressure_pa": (
                    standard.reference_pressure_pa if standard is not None else pressure_pa
                ),
            }
            for record in reaction.equilibrium_constants
        ),
        "temperature_k": temperature_k,
        "pressure_pa": pressure_pa,
    }


def _finite_phase_direction_spec(
    problem: GeneralReactiveEquilibriumProblem,
) -> dict[str, object] | None:
    system = problem.phase_reactions
    if system is None:
        return None
    preflight = preflight_phase_reaction_system(system)
    return {
        "phase_ids": preflight.phase_ids,
        "species_ids": preflight.species_ids,
        "charges": preflight.charges,
        "balance_matrix": preflight.balance_matrix,
        "direction_rows": preflight.direction_matrix,
        "driving_constants": tuple(
            0.0 if reaction.ln_equilibrium_constant is None else reaction.ln_equilibrium_constant
            for reaction in system.reactions
        ),
    }


def _finite_phase_support_indices(
    model: epcsaft.Mixture,
    compilation: _FinitePhaseCompilation,
) -> tuple[tuple[int, ...], ...]:
    component_index = {identity: index for index, identity in enumerate(model.component_ids)}
    return tuple(
        tuple(component_index[identity] for identity in supported)
        for supported in compilation.supported_component_ids
    )


def _reactive_finite_source_standard_state(
    model: epcsaft.Mixture,
    problem: GeneralReactiveEquilibriumProblem,
    phase: ReactivePhase,
    temperature_k: float,
    pressure_pa: float,
    active_parameters: epcsaft.ActiveParameterSet | None,
) -> dict[str, object] | None:
    reaction = _effective_reaction_system(problem)
    standard = reaction.source_standard_state
    if standard is None:
        return None
    provider_model = cast(ProviderModel, phase.model)
    transfer = _source_reference_transfer(
        ProviderPhase(model, provider_model.admissible_packing_fraction_interval),
        standard,
        species_ids=reaction.species_ids,
        temperature_k=temperature_k,
        pressure_pa=pressure_pa,
        sensitivities_requested=active_parameters is not None,
        active_parameters=active_parameters,
    )
    return {
        "id": standard.id,
        "activity_scale_id": standard.activity_scale_id,
        "log_activity_scale_factors": standard.log_activity_scale_factors,
        "reference_pressure_pa": standard.reference_pressure_pa,
        "source_reference_activity_convention_id": (
            standard.source_reference_activity_convention_id
        ),
        "source_reference_standard_molality_mol_per_kg": (
            standard.source_reference_standard_molality_mol_per_kg
        ),
        "provider_transfer": _source_reference_payload(transfer),
    }


def _finite_phase_failed(
    problem: GeneralReactiveEquilibriumProblem,
    compilation: _FinitePhaseCompilation,
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
) -> GeneralReactiveEquilibriumResult:
    return GeneralReactiveEquilibriumResult(
        problem=problem,
        descriptor=compilation.descriptor,
        status="non_evaluable",
        solver_status=solver_status,
        numerical_status=numerical_status,
        physical_status=physical_status,
        provider_domain_status=provider_domain_status,
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
        phases=phases,
        failure=NonEvaluableTrial(code, diagnostic),
        certificate_scope=compilation.descriptor.certificate_scope,
        globality_claim="not_established",
        trial=trial,
        branch_signature=None,
        evidence=(),
        resolved_policy=(),
    )


def _contract_parameter_derivatives(
    coefficients: tuple[float, ...],
    derivatives: tuple[float, ...],
    *,
    phase_index: int,
    component_count: int,
    parameter_count: int,
) -> tuple[float, ...]:
    return tuple(
        math.fsum(
            coefficient
            * derivatives[
                ((phase_index * component_count + component) * parameter_count) + parameter
            ]
            for component, coefficient in enumerate(coefficients)
        )
        for parameter in range(parameter_count)
    )


def _native_float_vector(value: object) -> tuple[float, ...]:
    return tuple(float(cast(Any, item)) for item in cast(tuple[object, ...], value))


def _native_float_matrix(value: object) -> tuple[tuple[float, ...], ...]:
    return tuple(_native_float_vector(row) for row in cast(tuple[object, ...], value))


def _finite_phase_states(
    model: epcsaft.Mixture,
    compilation: _FinitePhaseCompilation,
    native: dict[str, object],
    *,
    mechanical_class: str,
) -> tuple[EquilibriumPhaseState, ...]:
    native_amounts = _native_float_matrix(native["phase_amounts_mol"])
    native_volumes = _native_float_vector(native["phase_volumes_m3"])
    native_pressures = _native_float_vector(native["phase_pressures_pa"])
    native_packing = _native_float_vector(native["phase_packing_fractions"])
    native_mu = _native_float_matrix(native["phase_chemical_potentials_over_rt"])
    if not all(
        len(values) == len(compilation.phases)
        for values in (
            native_amounts,
            native_volumes,
            native_pressures,
            native_packing,
            native_mu,
        )
    ):
        return ()
    component_index = {identity: index for index, identity in enumerate(model.component_ids)}
    return tuple(
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
            mole_fractions=tuple(amount / math.fsum(amounts) for amount in amounts),
            molar_volume_m3_per_mol=volume / math.fsum(amounts),
            molar_density_mol_m3=math.fsum(amounts) / volume,
            packing_fraction=packing,
            pressure_pa=pressure,
            chemical_potential_component_ids=supported,
            chemical_potential_over_rt=tuple(
                chemical_potentials[component_index[identity]] for identity in supported
            ),
            mechanical_class=mechanical_class,
        )
        for phase, supported, amounts, volume, packing, pressure, chemical_potentials in zip(
            compilation.phases,
            compilation.supported_component_ids,
            native_amounts,
            native_volumes,
            native_packing,
            native_pressures,
            native_mu,
            strict=True,
        )
    )


def _finite_phase_rows(
    problem: GeneralReactiveEquilibriumProblem,
    compilation: _FinitePhaseCompilation,
    native: dict[str, object],
) -> tuple[GeneralEquilibriumRowResult, ...]:
    component_count = len(_effective_reaction_system(problem).species_ids)
    parameter_count = len(compilation.descriptor.active_parameter_identities)
    phase_index_by_id = {phase.identity: index for index, phase in enumerate(compilation.phases)}
    amounts = _native_float_matrix(native["phase_amounts_mol"])
    volumes = _native_float_vector(native["phase_volumes_m3"])
    chemical_potentials = _native_float_matrix(native["phase_chemical_potentials_over_rt"])
    log_fugacities = _native_float_matrix(native["phase_log_fugacities"])
    amount_derivatives = _native_float_vector(native["amount_derivatives"])
    volume_derivatives = _native_float_vector(native["volume_derivatives"])
    chemical_potential_derivatives = _native_float_vector(native["chemical_potential_derivatives"])
    log_fugacity_derivatives = _native_float_vector(native["log_fugacity_derivatives"])
    rows = []
    for output in problem.outputs:
        if output.selector == "system.pressure":
            value = compilation.pressure_pa
            jacobian = (0.0,) * parameter_count
        else:
            phase_index = phase_index_by_id[cast(str, output.phase_identity)]
            phase_amounts = amounts[phase_index]
            total_amount = math.fsum(phase_amounts)
            volume = volumes[phase_index]
            if output.selector in {
                "phase.amount",
                "phase.mole_fraction",
                "phase.molality",
                "phase.molar_density",
            }:
                value, jacobian = _contract_phase_output(
                    output.selector,
                    output.coefficients,
                    output.solvent_mass_coefficients_kg_per_mol,
                    phase_amounts,
                    volume,
                    tuple(
                        tuple(
                            amount_derivatives[
                                ((phase_index * component_count + component) * parameter_count)
                                + parameter
                            ]
                            for component in range(component_count)
                        )
                        for parameter in range(parameter_count)
                    ),
                    tuple(
                        volume_derivatives[phase_index * parameter_count + parameter]
                        for parameter in range(parameter_count)
                    ),
                )
            elif output.selector == "phase.compressibility_factor":
                value = (
                    compilation.pressure_pa
                    * volume
                    / (total_amount * 8.31446261815324 * compilation.temperature_k)
                )
                jacobian = tuple(
                    value
                    * (
                        volume_derivatives[phase_index * parameter_count + parameter] / volume
                        - math.fsum(
                            amount_derivatives[
                                ((phase_index * component_count + component) * parameter_count)
                                + parameter
                            ]
                            for component in range(component_count)
                        )
                        / total_amount
                    )
                    for parameter in range(parameter_count)
                )
            elif output.selector == "phase.pressure":
                value = compilation.pressure_pa
                jacobian = (0.0,) * parameter_count
            else:
                assert output.selector in {
                    "phase.chemical_potential_over_rt",
                    "phase.fugacity",
                    "phase.log_fugacity",
                }
                potentials = (
                    log_fugacities[phase_index]
                    if output.selector in {"phase.fugacity", "phase.log_fugacity"}
                    else chemical_potentials[phase_index]
                )
                derivatives = (
                    log_fugacity_derivatives
                    if output.selector in {"phase.fugacity", "phase.log_fugacity"}
                    else chemical_potential_derivatives
                )
                value = math.fsum(
                    coefficient * chemical_potential
                    for coefficient, chemical_potential in zip(
                        output.coefficients, potentials, strict=True
                    )
                )
                jacobian = _contract_parameter_derivatives(
                    output.coefficients,
                    derivatives,
                    phase_index=phase_index,
                    component_count=component_count,
                    parameter_count=parameter_count,
                )
                if output.selector == "phase.fugacity":
                    value = math.exp(value)
                    jacobian = tuple(value * derivative for derivative in jacobian)
        rows.append(
            GeneralEquilibriumRowResult(
                output.identity,
                output.selector,
                output.unit,
                output.basis,
                "evaluated",
                value,
                jacobian,
                "",
            )
        )
    return tuple(rows)


def _evaluate_finite_phase(
    model: epcsaft.Mixture,
    problem: GeneralReactiveEquilibriumProblem,
    compilation: _FinitePhaseCompilation,
    active_parameters: epcsaft.ActiveParameterSet | None,
) -> GeneralReactiveEquilibriumResult:
    trial = _finite_phase_trial(compilation.owner_descriptor, active_parameters)
    try:
        source_standard_state = _reactive_finite_source_standard_state(
            model,
            problem,
            compilation.phases[0],
            compilation.temperature_k,
            compilation.pressure_pa,
            active_parameters,
        )
        native_attempt_list: list[dict[str, object]] = []

        def candidate_start_sets() -> Iterator[tuple[FinitePhaseStart, ...]]:
            if compilation.declared_start_set is not None:
                yield compilation.declared_start_set
            for allocation in compilation.generated_allocations:
                try:
                    yield _finite_phase_start_set_from_allocation(
                        model,
                        compilation.phases,
                        _effective_reaction_system(problem).feed_amounts_mol,
                        allocation,
                        compilation.temperature_k,
                        compilation.pressure_pa,
                    )
                except (ArithmeticError, RuntimeError, ValueError):
                    continue

        for start_set in candidate_start_sets():
            attempt = cast(
                dict[str, object],
                _equilibrium._reactive_finite_phases(
                    _native_transport(model, active_parameters),
                    _finite_phase_spec(problem, compilation.temperature_k, compilation.pressure_pa),
                    source_standard_state,
                    tuple(
                        cast(ProviderModel, phase.model).admissible_packing_fraction_interval
                        for phase in compilation.phases
                    ),
                    tuple(
                        {
                            "amounts_mol": start.amounts_mol,
                            "molar_volume_m3_per_mol": start.molar_volume_m3_per_mol,
                        }
                        for start in start_set
                    ),
                    _effective_reaction_system(problem).strict_interior_amount_floor_mol,
                    active_parameters is not None,
                    False,
                    _finite_phase_direction_spec(problem),
                    _finite_phase_support_indices(model, compilation),
                ),
            )
            native_attempt_list.append(attempt)
            if bool(attempt["accepted"]):
                break
        if not native_attempt_list:
            raise ValueError("Provider could not construct a finite-phase candidate start")
        native_attempts = tuple(native_attempt_list)
    except (ArithmeticError, RuntimeError, ValueError) as error:
        diagnostic = str(error)
        if (
            "Provider" not in diagnostic
            and "source domain" not in diagnostic
            and "outside the declared packing interval" not in diagnostic
        ):
            raise
        return _finite_phase_failed(
            problem,
            compilation,
            trial,
            "provider_domain_rejection",
            diagnostic,
            solver_status="not_started",
            numerical_status="not_adjudicated",
            physical_status="not_adjudicated",
            provider_domain_status="failed",
        )
    accepted_attempts = tuple(
        (ordinal, attempt)
        for ordinal, attempt in enumerate(native_attempts)
        if bool(attempt["accepted"])
    )
    if accepted_attempts:
        _, selected = accepted_attempts[0]
    else:
        _, selected = min(enumerate(native_attempts), key=_finite_phase_failure_score)
    native = dict(selected)
    native["optimizer_iterations"] = sum(
        int(cast(int, attempt["optimizer_iterations"])) for attempt in native_attempts
    )
    native["provider_evaluation_count"] = sum(
        int(cast(int, attempt["provider_evaluation_count"])) for attempt in native_attempts
    )
    if not bool(native["accepted"]):
        diagnostic = ": ".join(
            part
            for part in (
                str(native["failure_reason"] or native["solver_status"]),
                str(native["callback_error"]),
            )
            if part
        )
        if str(native["physical_status"]) == "failed":
            diagnostic += (
                "; balance={balance}, pressure={pressure}, affinity={affinity}, "
                "transfer={transfer}, stationarity={stationarity}, "
                "distinctness={distinctness}, reduced_hessian={reduced_hessian}, "
                "phase_totals={phase_totals}"
            ).format(
                balance=native["balance_inf_norm"],
                pressure=native["pressure_relative_inf_norm"],
                affinity=native["reaction_affinity_inf_norm"],
                transfer=native["transfer_chemical_potential_inf_norm"],
                stationarity=native["kkt_stationarity_inf_norm"],
                distinctness=native["phase_distinctness"],
                reduced_hessian=native["reduced_hessian_min_eigenvalue"],
                phase_totals=tuple(
                    math.fsum(row) for row in _native_float_matrix(native["phase_amounts_mol"])
                ),
            )
        code = (
            "provider_domain_rejection"
            if str(native["provider_domain_status"]) == "failed"
            else "owner_contract_failure"
            if str(native["callback_failure_kind"]) == "internal_contract"
            else "first_order_failure"
            if "sensitivity" in diagnostic
            else "physical_rejection"
            if str(native["physical_status"]) == "failed"
            else "numerical_convergence_failure"
        )
        return _finite_phase_failed(
            problem,
            compilation,
            trial,
            cast(Any, code),
            diagnostic,
            solver_status=str(native["solver_status"]),
            numerical_status=str(native["numerical_status"]),
            physical_status=str(native["physical_status"]),
            provider_domain_status=str(native["provider_domain_status"]),
            phases=_finite_phase_states(
                model,
                compilation,
                native,
                mechanical_class="rejected_configured_candidate",
            ),
        )
    if tuple(str(value) for value in cast(tuple[object, ...], native["parameter_order"])) != (
        compilation.descriptor.active_parameter_identities
    ):
        return _finite_phase_failed(
            problem,
            compilation,
            trial,
            "first_order_failure",
            "finite-phase GREPE parameter derivative order differs from its descriptor",
            solver_status=str(native["solver_status"]),
            numerical_status=str(native["numerical_status"]),
            physical_status=str(native["physical_status"]),
            provider_domain_status=str(native["provider_domain_status"]),
        )
    rows = _finite_phase_rows(problem, compilation, native)
    phases = _finite_phase_states(
        model,
        compilation,
        native,
        mechanical_class="strict_stable_local_minimum",
    )
    branch_signature: tuple[object, ...] = (
        "finite_phases_local_minimum",
        compilation.descriptor.problem_fingerprint,
        tuple(
            (
                phase.identity,
                problem.continuation_identity,
                None
                if problem.continuation_state is None
                else problem.continuation_state.fingerprint,
            )
            for phase in compilation.phases
        ),
    )
    owner = _FinitePhaseOwnerResult(
        compilation.owner_descriptor,
        rows,
        phases,
        trial,
        compilation.temperature_k,
        compilation.pressure_pa,
        float(cast(float, native["objective_over_rt"])),
        float(cast(float, native["balance_inf_norm"])),
        float(cast(float, native["pressure_relative_inf_norm"])),
        float(cast(float, native["reaction_affinity_inf_norm"])),
        float(cast(float, native["transfer_chemical_potential_inf_norm"])),
        float(cast(float, native["kkt_stationarity_inf_norm"])),
        float(cast(float, native["source_reference_representation_residual_inf_norm"])),
        float(cast(float, native["phase_distinctness"])),
        float(cast(float, native["reduced_hessian_min_eigenvalue"])),
        int(cast(int, native["declared_direction_rank"])),
        int(cast(int, native["admissible_nullity"])),
        int(cast(int, native["optimizer_iterations"])),
        int(cast(int, native["provider_evaluation_count"])),
        branch_signature,
    )
    return GeneralReactiveEquilibriumResult(
        problem=problem,
        descriptor=compilation.descriptor,
        status="evaluated",
        solver_status=str(native["solver_status"]),
        numerical_status=str(native["numerical_status"]),
        physical_status=str(native["physical_status"]),
        provider_domain_status=str(native["provider_domain_status"]),
        rows=rows,
        phases=phases,
        failure=None,
        certificate_scope=compilation.descriptor.certificate_scope,
        globality_claim="not_established",
        trial=trial,
        branch_signature=branch_signature,
        evidence=_result_evidence(owner),
        resolved_policy=_native_resolved_policy(native),
    )


def _evaluate_single_finite_phase_with_balance_derivatives(
    model: epcsaft.Mixture,
    problem: GeneralReactiveEquilibriumProblem,
    compilation: _FinitePhaseCompilation,
    active_parameters: epcsaft.ActiveParameterSet | None,
) -> tuple[GeneralReactiveEquilibriumResult, tuple[tuple[float, ...], ...]]:
    phase = compilation.phases[0]
    trial = _finite_phase_trial(compilation.owner_descriptor, active_parameters)
    start = compilation.declared_start_set
    continuation_start = (
        None
        if problem.continuation_state is None or start is None
        else (
            start[0].amounts_mol,
            math.fsum(start[0].amounts_mol) * start[0].molar_volume_m3_per_mol,
        )
    )
    try:
        solved = _chemical_equilibrium_impl(
            ProviderPhase(
                model,
                cast(ProviderModel, phase.model).admissible_packing_fraction_interval,
            ),
            compilation.temperature_k * epcsaft.unit_registry.kelvin,
            compilation.pressure_pa * epcsaft.unit_registry.pascal,
            _effective_reaction_system(problem),
            sensitivity_request=ChemicalEquilibriumSensitivityRequest(active_parameters),
            continuation_start=continuation_start,
        )
    except ChemicalEquilibriumError as error:
        diagnostics = error.diagnostics
        code: Literal[
            "provider_domain_rejection",
            "singular_sensitivity",
            "physical_rejection",
            "numerical_convergence_failure",
        ]
        code = (
            "provider_domain_rejection"
            if diagnostics.provider_domain_status == "failed"
            else "singular_sensitivity"
            if diagnostics.first_failed_numerical_criterion
            == "kkt_condition_number_inf"
            else "physical_rejection"
            if diagnostics.physical_status == "failed"
            else "numerical_convergence_failure"
        )
        return _finite_phase_failed(
            problem,
            compilation,
            trial,
            code,
            str(error),
            solver_status=diagnostics.solver_status,
            numerical_status=diagnostics.numerical_status,
            physical_status=diagnostics.physical_status,
            provider_domain_status=diagnostics.provider_domain_status,
        ), ()
    sensitivity = solved.sensitivity
    if sensitivity is None or sensitivity.status != "available":
        return _finite_phase_failed(
            problem,
            compilation,
            trial,
            "singular_sensitivity",
            "finite-phase GREPE requires exact conditioned sensitivities",
            solver_status=solved.diagnostics.solver_status,
            numerical_status="failed",
            physical_status=solved.diagnostics.physical_status,
            provider_domain_status=solved.diagnostics.provider_domain_status,
        ), ()
    amounts = solved.amounts_mol
    volume = solved.volume_m3
    active = {
        parameter.identity: (dn, dv)
        for parameter, dn, dv in zip(
            sensitivity.parameters,
            sensitivity.amount_derivatives,
            sensitivity.volume_derivatives,
            strict=True,
        )
        if parameter.identity is not None
    }
    balance = tuple(
        (dn, dv)
        for parameter, dn, dv in zip(
            sensitivity.parameters,
            sensitivity.amount_derivatives,
            sensitivity.volume_derivatives,
            strict=True,
        )
        if parameter.kind == "compiled_balance_total"
    )
    component_count = len(amounts)
    active_count = 0 if active_parameters is None else len(active_parameters.identities)
    width = component_count + 1 + active_count
    if active_parameters is None:
        phase_block = _equilibrium.evaluate_parameterized_phase(
            _native_transport(model, None),
            compilation.temperature_k,
            amounts,
            volume,
            *cast(ProviderModel, phase.model).admissible_packing_fraction_interval,
        )
        hessian = _native_float_vector(phase_block["state_hessian"])
        potentials = _native_float_vector(phase_block["state_gradient"])[:component_count]
    else:
        derivative_block = model._native.evaluate_phase_derivatives(
            compilation.temperature_k,
            amounts,
            volume,
            model._active_coordinates(active_parameters),
            2,
        )
        hessian = tuple(float(value) for value in derivative_block.hessian)
        potentials = tuple(float(value) for value in derivative_block.gradient[:component_count])

    def native_rows(
        columns: tuple[tuple[tuple[float, ...], float], ...],
        *,
        direct: bool,
    ) -> dict[str, object]:
        parameter_count = len(columns)
        potential_derivatives = tuple(
            math.fsum(
                hessian[component * width + species] * dn[species]
                for species in range(component_count)
            )
            + hessian[component * width + component_count] * dv
            + (hessian[component * width + component_count + 1 + parameter] if direct else 0.0)
            for component in range(component_count)
            for parameter, (dn, dv) in enumerate(columns)
        )
        return {
            "phase_amounts_mol": (amounts,),
            "phase_volumes_m3": (volume,),
            "phase_pressures_pa": (solved.pressure_pa,),
            "phase_packing_fractions": (cast(float, solved.diagnostics.packing_fraction),),
            "phase_chemical_potentials_over_rt": (potentials,),
            "phase_log_fugacities": (
                tuple(
                    value + math.log(8.31446261815324 * compilation.temperature_k)
                    for value in potentials
                ),
            ),
            "amount_derivatives": tuple(
                dn[component] for component in range(component_count) for dn, _ in columns
            ),
            "volume_derivatives": tuple(dv for _, dv in columns),
            "chemical_potential_derivatives": potential_derivatives,
            "log_fugacity_derivatives": potential_derivatives,
            "parameter_order": tuple(range(parameter_count)),
        }

    active_columns = tuple(
        active[identity] for identity in compilation.descriptor.active_parameter_identities
    )
    active_native = native_rows(active_columns, direct=True)
    rows = _finite_phase_rows(problem, compilation, active_native)
    phases = _finite_phase_states(
        model,
        compilation,
        active_native,
        mechanical_class="strict_stable_local_minimum",
    )
    balance_compilation = replace(
        compilation,
        descriptor=replace(
            compilation.descriptor,
            active_parameter_identities=tuple(f"balance[{index}]" for index in range(len(balance))),
        ),
    )
    balance_rows = _finite_phase_rows(
        problem,
        balance_compilation,
        native_rows(balance, direct=False),
    )
    diagnostics = solved.diagnostics
    branch_signature = (
        "finite_phases_local_minimum",
        compilation.descriptor.problem_fingerprint,
        ((phase.identity, problem.continuation_identity, None),),
    )
    result = GeneralReactiveEquilibriumResult(
        problem,
        compilation.descriptor,
        "evaluated",
        diagnostics.solver_status,
        diagnostics.numerical_status,
        diagnostics.physical_status,
        diagnostics.provider_domain_status,
        rows,
        phases,
        None,
        compilation.descriptor.certificate_scope,
        "not_established",
        trial,
        branch_signature,
        tuple(
            (name, value)
            for name, value in (
                ("balance_inf_norm", diagnostics.balance_inf_norm),
                ("electroneutrality_inf", diagnostics.charge_inf_norm),
                ("reaction_affinity_inf_norm", diagnostics.reaction_affinity_inf_norm),
                ("pressure_relative_inf_norm", diagnostics.pressure_relative_residual),
                ("kkt_stationarity_inf_norm", diagnostics.kkt_stationarity_inf_norm),
                (
                    "source_reference_representation_residual_inf_norm",
                    diagnostics.reference_representation_residual_inf_norm,
                ),
                ("provider_evaluation_count", diagnostics.provider_evaluation_count),
            )
            if value is not None
        ),
        solved.resolved_policy,
    )
    return result, tuple(cast(tuple[float, ...], row.jacobian) for row in balance_rows)


def _evaluate_ideal_phase(
    model: epcsaft.Mixture,
    problem: GeneralReactiveEquilibriumProblem,
    compilation: _IdealPhaseCompilation,
    continuation_start: tuple[tuple[float, ...], float] | None,
) -> tuple[GeneralReactiveEquilibriumResult, tuple[tuple[float, ...], ...]]:
    trial = _finite_phase_trial(compilation.owner_descriptor, None)
    phase = compilation.phase
    if continuation_start is None and phase.start is not None:
        continuation_start = (
            phase.start.amounts_mol,
            math.fsum(phase.start.amounts_mol) * phase.start.molar_volume_m3_per_mol,
        )
    try:
        solved = _chemical_equilibrium_impl(
            IdealGasPhase(
                model.parameter_fingerprint, cast(IdealGasModel, phase.model).reference_id
            ),
            compilation.temperature_k * epcsaft.unit_registry.kelvin,
            compilation.pressure_pa * epcsaft.unit_registry.pascal,
            _effective_reaction_system(problem),
            sensitivity_request=ChemicalEquilibriumSensitivityRequest(),
            continuation_start=continuation_start,
        )
    except ChemicalEquilibriumError as error:
        diagnostics = error.diagnostics
        failure = NonEvaluableTrial(
            (
                "singular_sensitivity"
                if diagnostics.first_failed_numerical_criterion
                == "kkt_condition_number_inf"
                else "numerical_convergence_failure"
            ),
            str(error),
        )
        return GeneralReactiveEquilibriumResult(
            problem,
            compilation.descriptor,
            "non_evaluable",
            diagnostics.solver_status,
            diagnostics.numerical_status,
            diagnostics.physical_status,
            diagnostics.provider_domain_status,
            tuple(
                GeneralEquilibriumRowResult(
                    output.identity,
                    output.selector,
                    output.unit,
                    output.basis,
                    "non_evaluable",
                    None,
                    None,
                    str(error),
                )
                for output in problem.outputs
            ),
            (),
            failure,
            compilation.descriptor.certificate_scope,
            "not_established",
            trial,
            None,
            (),
            (),
        ), ()
    amounts = solved.amounts_mol
    volume = solved.volume_m3
    potentials = tuple(math.log(amount / volume) for amount in amounts)
    pseudo = _FinitePhaseCompilation(
        compilation.temperature_k,
        compilation.pressure_pa,
        (phase,),
        (model.component_ids,),
        None,
        (),
        compilation.owner_descriptor,
        compilation.descriptor,
    )

    def native(
        columns: tuple[tuple[tuple[float, ...], float], ...],
    ) -> dict[str, object]:
        derivatives = tuple(
            dn[component] / amounts[component] - dv / volume
            for component in range(len(amounts))
            for dn, dv in columns
        )
        return {
            "phase_amounts_mol": (amounts,),
            "phase_volumes_m3": (volume,),
            "phase_pressures_pa": (solved.pressure_pa,),
            "phase_packing_fractions": (0.0,),
            "phase_chemical_potentials_over_rt": (potentials,),
            "phase_log_fugacities": (
                tuple(
                    value + math.log(8.31446261815324 * solved.temperature_k)
                    for value in potentials
                ),
            ),
            "amount_derivatives": tuple(
                dn[component] for component in range(len(amounts)) for dn, _ in columns
            ),
            "volume_derivatives": tuple(dv for _, dv in columns),
            "chemical_potential_derivatives": derivatives,
            "log_fugacity_derivatives": derivatives,
        }

    rows = _finite_phase_rows(problem, pseudo, native(()))
    phases = _finite_phase_states(
        model, pseudo, native(()), mechanical_class="analytical_ideal_gas"
    )
    sensitivity = cast(Any, solved.sensitivity)
    balance = tuple(
        (dn, dv)
        for parameter, dn, dv in zip(
            sensitivity.parameters,
            sensitivity.amount_derivatives,
            sensitivity.volume_derivatives,
            strict=True,
        )
        if parameter.kind == "compiled_balance_total"
    )
    balance_pseudo = replace(
        pseudo,
        descriptor=replace(
            pseudo.descriptor,
            active_parameter_identities=tuple(f"balance[{index}]" for index in range(len(balance))),
        ),
    )
    balance_rows = _finite_phase_rows(problem, balance_pseudo, native(balance))
    diagnostics = solved.diagnostics
    return GeneralReactiveEquilibriumResult(
        problem,
        compilation.descriptor,
        "evaluated",
        diagnostics.solver_status,
        diagnostics.numerical_status,
        diagnostics.physical_status,
        diagnostics.provider_domain_status,
        rows,
        phases,
        None,
        compilation.descriptor.certificate_scope,
        "not_established",
        trial,
        (
            "finite_phases_local_minimum",
            compilation.descriptor.problem_fingerprint,
            ((phase.identity, problem.continuation_identity, None),),
        ),
        tuple(
            (name, value)
            for name, value in (
                ("balance_inf_norm", diagnostics.balance_inf_norm),
                ("electroneutrality_inf", diagnostics.charge_inf_norm),
                ("reaction_affinity_inf_norm", diagnostics.reaction_affinity_inf_norm),
                ("pressure_relative_inf_norm", diagnostics.pressure_relative_residual),
                ("kkt_stationarity_inf_norm", diagnostics.kkt_stationarity_inf_norm),
                ("provider_evaluation_count", diagnostics.provider_evaluation_count),
            )
            if value is not None
        ),
        solved.resolved_policy,
    ), tuple(cast(tuple[float, ...], row.jacobian) for row in balance_rows)


def _failure_axes(failure: NonEvaluableTrial | None) -> tuple[str, str, str]:
    if failure is None:
        return "passed", "passed", "passed"
    if failure.code == "provider_domain_rejection":
        return "not_adjudicated", "not_adjudicated", "failed"
    if failure.code in {
        "numerical_convergence_failure",
        "singular_sensitivity",
        "ill_conditioned_sensitivity",
    }:
        return "failed", "not_adjudicated", "passed"
    if failure.code == "first_order_failure":
        return "passed", "passed", "passed"
    return "passed", "failed", "passed"


def _compiled_balance_rows(
    reaction: ChemicalEquilibriumProblem,
) -> tuple[tuple[float, ...], ...]:
    """Mirror the native deterministic independent-balance row selection."""

    selected: list[tuple[float, ...]] = [tuple(reaction.molar_masses_kg_per_mol)]
    basis: list[tuple[float, ...]] = list(selected)
    charge = tuple(float(value) for value in reaction.charges)
    if any(charge) and _matrix_rank(tuple(basis + [charge])) > len(basis):
        basis.append(charge)
    for row in reaction.balance_matrix:
        candidate = tuple(row)
        if _matrix_rank(tuple(basis + [candidate])) > len(basis):
            basis.append(candidate)
            selected.append(candidate)
    return tuple(selected)


def _intensive_failure(
    problem: GeneralReactiveEquilibriumProblem,
    compilation: _IntensiveBoundaryCompilation,
    failure: NonEvaluableTrial,
    owner: _IntensiveOwnerResult | None = None,
) -> GeneralReactiveEquilibriumResult:
    numerical_status, physical_status, provider_status = _failure_axes(failure)
    return GeneralReactiveEquilibriumResult(
        problem=problem,
        descriptor=compilation.descriptor,
        status="non_evaluable",
        solver_status=failure.code,
        numerical_status=numerical_status,
        physical_status=physical_status,
        provider_domain_status=provider_status,
        rows=tuple(
            GeneralEquilibriumRowResult(
                output.identity,
                output.selector,
                output.unit,
                output.basis,
                "non_evaluable",
                None,
                None,
                failure.diagnostic,
            )
            for output in problem.outputs
        ),
        phases=(() if owner is None else owner.phases),
        failure=failure,
        certificate_scope=compilation.descriptor.certificate_scope,
        globality_claim="not_established",
        trial=None if owner is None else owner.trial,
        branch_signature=None if owner is None else owner.branch_signature,
        evidence=() if owner is None else owner.evidence,
        resolved_policy=(),
    )


class _IntensiveOwnerFailure(RuntimeError):
    """Preserve a wrapped owner's typed failure across the outer solve."""

    def __init__(self, owner: _IntensiveOwnerResult) -> None:
        self.owner = owner
        diagnostic = "wrapped equilibrium state is non-evaluable"
        if owner.failure is not None:
            diagnostic = owner.failure.diagnostic
        super().__init__(diagnostic)


def _evaluate_fixed_pressure_boundary(
    model: epcsaft.Mixture,
    problem: GeneralReactiveEquilibriumProblem,
    compilation: _IntensiveBoundaryCompilation,
    active_parameters: epcsaft.ActiveParameterSet | None,
) -> GeneralReactiveEquilibriumResult:
    """Close one declared affine feed against an imposed bubble pressure."""

    feed = cast(AffineFeed, problem.feed).reorder(problem.reaction_system.species_ids)
    boundary = problem.intensive_boundaries[0]
    coordinate = feed.coordinates[0]
    balance_direction = tuple(
        math.fsum(
            balance[species] * feed.basis_matrix_mol_per_coordinate[species][0]
            for species in range(len(feed.component_ids))
        )
        for balance in _compiled_balance_rows(problem.reaction_system)
    )

    def evaluate(
        value: float,
    ) -> tuple[
        _FiniteIncipientKernel,
        _FinitePhaseOwnerDescriptor,
        Any,
        float,
    ]:
        internal_problem = _resolved_intensive_problem(
            problem, value, compilation.augmented_outputs
        )
        internal = _compile_intensive_owner(
            model, internal_problem, active_parameters=active_parameters
        )
        if not isinstance(internal, _IncipientCompilation):
            raise _IncipientFailure(
                "owner_contract_failure",
                "fixed-pressure affine feed did not compile to finite-incipient equations",
            )
        owner_problem = replace(
            internal.owner_problem,
            continuation_state=None,
        )
        owner_descriptor = internal.owner_descriptor
        closure = _closure(
            owner_problem,
            internal.temperature,
            owner_descriptor,
            active_parameters,
            boundary.target,
        )
        slope, _ = _fixed_pressure_coordinate_derivatives(closure, balance_direction)
        return owner_problem, owner_descriptor, closure, slope

    lower, upper = coordinate.bounds
    width = upper - lower

    def evaluate_scalar(
        value: float,
    ) -> tuple[
        tuple[_FiniteIncipientKernel, _FinitePhaseOwnerDescriptor, Any],
        float,
        float,
    ]:
        owner_problem, owner_descriptor, closure, slope = evaluate(value)
        return (owner_problem, owner_descriptor, closure), closure.residual, slope

    try:
        solved_coordinate, selected_state, _ = _solve_bracketed_scalar(
            evaluate_scalar,
            lower,
            upper,
            initial=boundary.initial,
            bracket_points=(
                boundary.initial,
                lower + 0.25 * width,
                lower + 0.5 * width,
                lower + 0.75 * width,
            ),
            tolerance=2.0e-10,
            maximum_iterations=64,
        )
    except _IncipientFailure as error:
        return _intensive_failure(
            problem, compilation, NonEvaluableTrial(error.code, str(error)), None
        )
    except _ScalarBoundarySolveError as error:
        code: Literal["no_incipient_phase", "numerical_convergence_failure"] = (
            "no_incipient_phase"
            if error.reason == "not_bracketed"
            else "numerical_convergence_failure"
        )
        return _intensive_failure(problem, compilation, NonEvaluableTrial(code, str(error)), None)
    owner_problem, owner_descriptor, closure = selected_state
    try:
        row_values, evidence = _fixed_pressure_coordinate_result(
            owner_problem, owner_descriptor, closure, balance_direction
        )
    except _IncipientFailure as error:
        return _intensive_failure(
            problem, compilation, NonEvaluableTrial(error.code, str(error)), None
        )
    rows = tuple(
        GeneralEquilibriumRowResult(
            output.identity,
            output.selector,
            output.unit,
            output.basis,
            "evaluated",
            row_values[index][0],
            row_values[index][1],
            "",
        )
        for index, output in enumerate(problem.outputs)
    )
    diagnostics = closure.equilibrium.diagnostics
    return GeneralReactiveEquilibriumResult(
        problem=problem,
        descriptor=compilation.descriptor,
        status="evaluated",
        solver_status="converged",
        numerical_status="passed",
        physical_status="passed",
        provider_domain_status=(
            "passed" if diagnostics is None else diagnostics.provider_domain_status
        ),
        rows=rows,
        phases=closure.phases,
        failure=None,
        certificate_scope=compilation.descriptor.certificate_scope,
        globality_claim="not_established",
        trial=closure.trial,
        branch_signature=(
            "intensive_boundary",
            ((boundary.coordinate_identity, solved_coordinate),),
            _incipient_branch_signature(closure),
        ),
        evidence=evidence,
        resolved_policy=closure.equilibrium.resolved_policy,
    )


def _evaluate_intensive_boundary(
    model: epcsaft.Mixture,
    problem: GeneralReactiveEquilibriumProblem,
    compilation: _IntensiveBoundaryCompilation,
    active_parameters: epcsaft.ActiveParameterSet | None,
) -> GeneralReactiveEquilibriumResult:
    if compilation.base_formulation == "finite_incipient":
        return _evaluate_fixed_pressure_boundary(model, problem, compilation, active_parameters)
    feed = cast(AffineFeed, problem.feed).reorder(problem.reaction_system.species_ids)
    boundary = problem.intensive_boundaries[0]
    coordinate_index = {
        coordinate.identity: index for index, coordinate in enumerate(feed.coordinates)
    }
    feed_column = coordinate_index[boundary.coordinate_identity]
    balance_rows = _compiled_balance_rows(problem.reaction_system)
    balance_direction = tuple(
        math.fsum(
            balance[species] * feed.basis_matrix_mol_per_coordinate[species][feed_column]
            for species in range(len(feed.component_ids))
        )
        for balance in balance_rows
    )
    bounds_by_id = {coordinate.identity: coordinate.bounds for coordinate in feed.coordinates}
    lower, upper = bounds_by_id[boundary.coordinate_identity]

    def evaluate(
        coordinate: float,
    ) -> tuple[
        _IntensiveOwnerResult,
        float,
        float,
        tuple[float, ...],
    ]:
        internal_problem = _resolved_intensive_problem(
            problem, coordinate, compilation.augmented_outputs
        )
        internal = _compile_intensive_owner(
            model, internal_problem, active_parameters=active_parameters
        )
        if isinstance(internal, _FinitePhaseCompilation):
            owner, balance_derivatives = _evaluate_single_finite_phase_with_balance_derivatives(
                model, internal_problem, internal, active_parameters
            )
        elif isinstance(internal, _IdealPhaseCompilation):
            owner, balance_derivatives = _evaluate_ideal_phase(
                model, internal_problem, internal, None
            )
        else:
            raise TypeError("intensive boundary compiled to a nondifferentiable formulation")
        if owner.status != "evaluated":
            raise _IntensiveOwnerFailure(owner)
        values = tuple(float(cast(float, row.value)) for row in owner.rows)
        if any(len(row) != len(balance_rows) for row in balance_derivatives):
            raise ValueError("finite formulation balance basis does not match GREPE")
        coordinate_derivatives = tuple(
            math.fsum(
                derivative * balance_direction[balance] for balance, derivative in enumerate(row)
            )
            for row in balance_derivatives
        )
        boundary_value = values[len(problem.outputs)]
        residual = (
            math.log(boundary_value / boundary.target)
            if boundary.selector == "system.pressure"
            else boundary_value - boundary.target
        )
        boundary_derivative = coordinate_derivatives[len(problem.outputs)]
        if boundary.selector == "system.pressure":
            boundary_derivative /= boundary_value
        return owner, residual, boundary_derivative, coordinate_derivatives

    coordinate = boundary.initial
    owner: _IntensiveOwnerResult | None = None
    coordinate_derivatives: tuple[float, ...] = ()
    boundary_derivative = 0.0
    try:
        for _ in range(32):
            owner, residual, boundary_derivative, coordinate_derivatives = evaluate(coordinate)
            if abs(residual) <= 2.0e-11:
                break
            if abs(boundary_derivative) <= 1.0e-13:
                raise ValueError("intensive-boundary Jacobian is singular")
            step = -residual / boundary_derivative
            maximum_step = 1.0
            if step > 0.0:
                maximum_step = min(maximum_step, 0.95 * (upper - coordinate) / step)
            elif step < 0.0:
                maximum_step = min(maximum_step, 0.95 * (lower - coordinate) / step)
            accepted = False
            scale = max(0.0, maximum_step)
            for _ in range(20):
                candidate = coordinate + scale * step
                candidate_owner, candidate_residual, candidate_derivative, candidate_outputs = (
                    evaluate(candidate)
                )
                if abs(candidate_residual) < abs(residual):
                    coordinate = candidate
                    owner = candidate_owner
                    boundary_derivative = candidate_derivative
                    coordinate_derivatives = candidate_outputs
                    accepted = True
                    break
                scale *= 0.5
            if not accepted:
                raise ValueError("intensive-boundary Newton step did not reduce the residual")
        else:
            raise ValueError("intensive-boundary solve exceeded its iteration budget")
        if owner is None:
            raise ValueError("intensive-boundary solve returned no state")
        owner, residual, boundary_derivative, coordinate_derivatives = evaluate(coordinate)
        if abs(residual) > 2.0e-11:
            raise ValueError("intensive-boundary residual did not converge")
        if abs(boundary_derivative) <= 1.0e-13:
            raise ValueError("intensive-boundary Jacobian is singular")
        parameter_count = len(compilation.descriptor.active_parameter_identities)
        boundary_parameter_jacobian = cast(
            tuple[float, ...], owner.rows[len(problem.outputs)].jacobian
        )
        if boundary.selector == "system.pressure":
            boundary_parameter_jacobian = tuple(
                derivative / cast(float, owner.rows[len(problem.outputs)].value)
                for derivative in boundary_parameter_jacobian
            )
        coordinate_parameter_jacobian = tuple(
            -boundary_parameter_jacobian[parameter] / boundary_derivative
            for parameter in range(parameter_count)
        )
        generic_rows = []
        for index, output in enumerate(problem.outputs):
            row = owner.rows[index]
            fixed_coordinate_jacobian = cast(tuple[float, ...], row.jacobian)
            total_jacobian = tuple(
                fixed_coordinate_jacobian[parameter]
                + coordinate_derivatives[index] * coordinate_parameter_jacobian[parameter]
                for parameter in range(parameter_count)
            )
            generic_rows.append(
                GeneralEquilibriumRowResult(
                    output.identity,
                    output.selector,
                    output.unit,
                    output.basis,
                    "evaluated",
                    row.value,
                    total_jacobian,
                    "",
                )
            )
    except _IntensiveOwnerFailure as error:
        failure = error.owner.failure or NonEvaluableTrial("owner_contract_failure", str(error))
        return _intensive_failure(problem, compilation, failure, error.owner)
    except (RuntimeError, TypeError, ValueError) as error:
        failure = (
            NonEvaluableTrial("singular_sensitivity", str(error))
            if "Jacobian is singular" in str(error)
            else NonEvaluableTrial("numerical_convergence_failure", str(error))
        )
        return _intensive_failure(
            problem,
            compilation,
            failure,
            owner,
        )
    return GeneralReactiveEquilibriumResult(
        problem=problem,
        descriptor=compilation.descriptor,
        status="evaluated",
        solver_status="converged",
        numerical_status="passed",
        physical_status="passed",
        provider_domain_status=owner.provider_domain_status,
        rows=tuple(generic_rows),
        phases=owner.phases,
        failure=None,
        certificate_scope=compilation.descriptor.certificate_scope,
        globality_claim="not_established",
        trial=owner.trial,
        branch_signature=(
            "intensive_boundary",
            ((boundary.coordinate_identity, coordinate),),
            owner.branch_signature,
        ),
        evidence=owner.evidence,
        resolved_policy=owner.resolved_policy,
    )


def _evaluate_temperature_boundary(
    model: epcsaft.Mixture,
    problem: GeneralReactiveEquilibriumProblem,
    compilation: _TemperatureBoundaryCompilation,
    active_parameters: epcsaft.ActiveParameterSet | None,
) -> GeneralReactiveEquilibriumResult:
    trial = _finite_phase_trial(compilation.owner_descriptor, active_parameters)
    try:
        solved = solve(
            model,
            compilation.owner_problem,
            active_parameters=active_parameters,
        )
    except PhaseEquilibriumError as error:
        reason = error.diagnostics.failure_reason or str(error)
        if "domain" in reason or "provider" in reason:
            code = "provider_domain_rejection"
        elif "singular" in reason:
            code = "singular_sensitivity"
        elif "condition" in reason:
            code = "ill_conditioned_sensitivity"
        elif "iteration" in reason:
            code = "numerical_convergence_failure"
        else:
            code = "physical_rejection"
        failure = NonEvaluableTrial(cast(Any, code), reason)
        numerical, physical, provider = _failure_axes(failure)
        return GeneralReactiveEquilibriumResult(
            problem,
            compilation.descriptor,
            "non_evaluable",
            error.diagnostics.solver_status,
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
                    reason,
                )
                for output in problem.outputs
            ),
            error.best_candidate,
            failure,
            compilation.descriptor.certificate_scope,
            "not_established",
            trial,
            None,
            (("provider_evaluation_count", error.diagnostics.provider_evaluation_count),),
            (),
        )

    extent = solved.extent
    if not isinstance(extent, IncipientTemperatureBoundary):
        raise TypeError("fixed-pressure coexistence returned the wrong extent")
    parameter_count = len(compilation.owner_descriptor.active_parameter_identities)
    sensitivity = solved.sensitivity
    if active_parameters is not None and (sensitivity is None or sensitivity.status != "available"):
        reason = (
            "fixed-pressure coexistence did not provide exact total derivatives"
            if sensitivity is None
            else sensitivity.failure_reason
            or "fixed-pressure coexistence sensitivity is unavailable"
        )
        failure = NonEvaluableTrial("owner_contract_failure", reason)
        return GeneralReactiveEquilibriumResult(
            problem,
            compilation.descriptor,
            "non_evaluable",
            solved.diagnostics.solver_status,
            "failed",
            "passed",
            "passed",
            tuple(
                GeneralEquilibriumRowResult(
                    output.identity,
                    output.selector,
                    output.unit,
                    output.basis,
                    "non_evaluable",
                    None,
                    None,
                    reason,
                )
                for output in problem.outputs
            ),
            solved.phases,
            failure,
            compilation.descriptor.certificate_scope,
            "not_established",
            trial,
            None,
            (),
            (),
        )

    derivative_by_identity = (
        {}
        if sensitivity is None
        else dict(
            zip(
                sensitivity.output_identities,
                sensitivity.output_parameter_jacobian,
                strict=True,
            )
        )
    )
    phase_index = {phase.identity: index for index, phase in enumerate(problem.phases)}
    rows: list[GeneralEquilibriumRowResult] = []
    for output in problem.outputs:
        if output.selector == "system.temperature":
            value = extent.temperature_k
            jacobian = derivative_by_identity.get(
                "boundary.temperature_k", (0.0,) * parameter_count
            )
        elif output.selector == "system.pressure":
            value = extent.pressure_pa
            jacobian = (0.0,) * parameter_count
        else:
            index = phase_index[cast(str, output.phase_identity)]
            phase = solved.phases[index]
            if output.selector == "phase.mole_fraction":
                value = math.fsum(
                    coefficient * fraction
                    for coefficient, fraction in zip(
                        output.coefficients, phase.mole_fractions, strict=True
                    )
                )
                jacobian = tuple(
                    math.fsum(
                        coefficient
                        * derivative_by_identity[f"phase[{index}].mole_fraction[{component_id}]"][
                            parameter
                        ]
                        for coefficient, component_id in zip(
                            output.coefficients, phase.supported_component_ids, strict=True
                        )
                    )
                    for parameter in range(parameter_count)
                )
            elif output.selector == "phase.molar_density":
                value = phase.molar_density_mol_m3
                volume_derivative = derivative_by_identity.get(
                    f"phase[{index}].molar_volume_m3_per_mol",
                    (0.0,) * parameter_count,
                )
                jacobian = tuple(-derivative * value * value for derivative in volume_derivative)
            else:
                value = phase.pressure_pa
                jacobian = (0.0,) * parameter_count
        rows.append(
            GeneralEquilibriumRowResult(
                output.identity,
                output.selector,
                output.unit,
                output.basis,
                "evaluated",
                value,
                tuple(jacobian),
                "",
            )
        )

    diagnostics = solved.diagnostics
    return GeneralReactiveEquilibriumResult(
        problem,
        compilation.descriptor,
        "evaluated",
        diagnostics.solver_status,
        diagnostics.numerical_status,
        diagnostics.physical_status,
        "passed",
        tuple(rows),
        solved.phases,
        None,
        compilation.descriptor.certificate_scope,
        "not_established",
        trial,
        (
            "fixed_p_declared_coexistence",
            extent.reference_phase_index,
            tuple((phase.fluid_role, type(phase.support).__name__) for phase in problem.phases),
        ),
        (
            ("temperature_k", extent.temperature_k),
            ("pressure_residual_inf", diagnostics.pressure_residual_inf),
            ("optimizer_solve_count", diagnostics.optimizer_solve_count),
            ("provider_evaluation_count", diagnostics.provider_evaluation_count),
        ),
        (),
    )


def _solve_general_reactive_equilibrium(
    model: epcsaft.Mixture,
    problem: GeneralReactiveEquilibriumProblem,
    *,
    active_parameters: epcsaft.ActiveParameterSet | None = None,
    continuation_start: tuple[tuple[float, ...], float] | None = None,
) -> GeneralReactiveEquilibriumResult:
    """Compile and solve one role-declared GREPE condition."""

    compiled = _compile_general_reactive_equilibrium(
        model, problem, active_parameters=active_parameters
    )
    if isinstance(compiled, str):
        return _unavailable_formulation(problem, compiled)
    if isinstance(compiled, _IntensiveBoundaryCompilation):
        return _with_continuation_state(
            _evaluate_intensive_boundary(model, problem, compiled, active_parameters)
        )
    if isinstance(compiled, _TemperatureBoundaryCompilation):
        return _with_continuation_state(
            _evaluate_temperature_boundary(model, problem, compiled, active_parameters)
        )
    if isinstance(compiled, _FinitePhaseCompilation):
        return _with_continuation_state(
            _evaluate_single_finite_phase_with_balance_derivatives(
                model, problem, compiled, active_parameters
            )[0]
            if len(compiled.phases) == 1
            else _evaluate_finite_phase(model, problem, compiled, active_parameters)
        )
    if isinstance(compiled, _IdealPhaseCompilation):
        return _with_continuation_state(
            _evaluate_ideal_phase(model, problem, compiled, continuation_start)[0]
        )
    closure, failure, row_values, _, evidence = _evaluate_incipient(
        compiled.owner_problem,
        compiled.temperature,
        compiled.owner_descriptor,
        active_parameters=active_parameters,
    )
    numerical_status, physical_status, provider_status = _failure_axes(failure)
    diagnostics = None if closure is None else closure.equilibrium.diagnostics
    rows = (
        tuple(
            GeneralEquilibriumRowResult(
                output.identity,
                output.selector,
                output.unit,
                output.basis,
                "non_evaluable",
                None,
                None,
                "incipient state is non-evaluable" if failure is None else failure.diagnostic,
            )
            for output in problem.outputs
        )
        if closure is None
        else tuple(
            GeneralEquilibriumRowResult(
                output.identity,
                output.selector,
                output.unit,
                output.basis,
                "evaluated",
                value,
                jacobian,
                "",
            )
            for output, (value, jacobian) in zip(problem.outputs, row_values, strict=True)
        )
    )
    return _with_continuation_state(
        GeneralReactiveEquilibriumResult(
            problem=problem,
            descriptor=compiled.descriptor,
            status="non_evaluable" if closure is None else "evaluated",
            solver_status=(
                str(diagnostics.solver_status)
                if diagnostics is not None
                else failure.code
                if failure is not None
                else "not_started"
            ),
            numerical_status="passed" if closure is not None else numerical_status,
            physical_status="passed" if closure is not None else physical_status,
            provider_domain_status=(
                str(diagnostics.provider_domain_status)
                if diagnostics is not None
                else provider_status
            ),
            rows=rows,
            phases=() if closure is None else closure.phases,
            failure=failure,
            certificate_scope=compiled.descriptor.certificate_scope,
            globality_claim="not_established",
            trial=None if closure is None else closure.trial,
            branch_signature=None if closure is None else _incipient_branch_signature(closure),
            evidence=evidence,
            resolved_policy=() if closure is None else closure.equilibrium.resolved_policy,
        )
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
