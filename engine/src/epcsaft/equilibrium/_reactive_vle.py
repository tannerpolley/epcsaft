from __future__ import annotations

import hashlib
import json
import math
from collections.abc import Sequence
from dataclasses import dataclass, replace
from typing import Any, Literal, cast

import epcsaft

from .._trials import NonEvaluableTrial, NonEvaluableTrialCode
from ..eos import _native_transport
from . import _equilibrium
from ._api import (
    AllComponents,
    ChemicalEquilibriumSensitivity,
    DeclaredNeutralComponents,
    DeclaredTopology,
    EquilibriumPhaseState,
    FixedTCoexistence,
    IncipientBoundary,
    PhaseEquilibriumDiagnostics,
    PhaseEquilibriumError,
    PhaseEquilibriumProblem,
    PhaseEquilibriumResult,
    PhaseSpecification,
    ProviderPhase,
    _fixed_t_coexistence_with_reference_directions,
    _ReferenceDirectionSensitivity,
)
from ._observations import (
    ActiveParameterTrialReceipt,
    CertifiedHomogeneousState,
    HomogeneousReactiveObservationProblem,
    HomogeneousReactiveObservationRow,
    PhysicalBranchReceipt,
    _canonical_fingerprint,
    evaluate_homogeneous_reactive_observations,
    homogeneous_reactive_observation_descriptor,
)

type ReactiveBubbleQuantity = Literal[
    "bubble_pressure_pa",
    "vapor_mole_fraction",
    "vapor_partial_pressure_pa",
    "liquid_amount_mol",
    "liquid_mole_fraction",
    "liquid_molality_mol_per_kg",
    "liquid_molar_density_mol_per_m3",
]
type ReactiveVaporModelKind = Literal["provider_nonideal", "ideal_gas_comparison"]
type ObservationSupport = Literal["real", "positive", "left_censored", "right_censored"]

_ROW_UNITS: dict[ReactiveBubbleQuantity, str] = {
    "bubble_pressure_pa": "pascal",
    "vapor_mole_fraction": "dimensionless",
    "vapor_partial_pressure_pa": "pascal",
    "liquid_amount_mol": "mole",
    "liquid_mole_fraction": "dimensionless",
    "liquid_molality_mol_per_kg": "mole / kilogram",
    "liquid_molar_density_mol_per_m3": "mole / meter ** 3",
}


def _canonical_json(value: object) -> str:
    return json.dumps(
        value,
        ensure_ascii=False,
        allow_nan=False,
        separators=(",", ":"),
        sort_keys=True,
    )


@dataclass(frozen=True, slots=True)
class ProviderNonidealVapor:
    """The declared neutral subset evaluated by the installed nonideal EOS."""

    reference_id: str

    def __post_init__(self) -> None:
        if not self.reference_id:
            raise ValueError("nonideal vapor model requires a reference identity")


@dataclass(frozen=True, slots=True)
class IdealGasVaporComparison:
    """Explicit ideal-vapor comparison model; never a failure fallback."""

    reference_id: str

    def __post_init__(self) -> None:
        if not self.reference_id:
            raise ValueError("ideal-vapor comparison requires a reference identity")


ReactiveVaporModel = ProviderNonidealVapor | IdealGasVaporComparison


@dataclass(frozen=True, slots=True)
class ReactiveBubbleContinuationReference:
    """Immutable local coupled-branch anchor used by pure trial evaluations."""

    pressure_pa: float
    vapor_mole_fractions: tuple[float, ...]
    liquid_molar_volume_m3_per_mol: float
    vapor_molar_volume_m3_per_mol: float
    maximum_log_pressure_distance: float
    maximum_log_composition_distance: float
    maximum_log_volume_distance: float
    anchor_state_identity: str
    anchor_problem_fingerprint: str
    anchor_topology_fingerprint: str
    anchor_parameter_fingerprint: str

    def __post_init__(self) -> None:
        object.__setattr__(self, "vapor_mole_fractions", tuple(self.vapor_mole_fractions))
        if (
            not math.isfinite(self.pressure_pa)
            or self.pressure_pa <= 0.0
            or not self.vapor_mole_fractions
            or any(not math.isfinite(value) or value <= 0.0 for value in self.vapor_mole_fractions)
            or not math.isclose(
                math.fsum(self.vapor_mole_fractions), 1.0, rel_tol=0.0, abs_tol=1.0e-12
            )
            or any(
                not math.isfinite(value) or value <= 0.0
                for value in (
                    self.liquid_molar_volume_m3_per_mol,
                    self.vapor_molar_volume_m3_per_mol,
                    self.maximum_log_pressure_distance,
                    self.maximum_log_composition_distance,
                    self.maximum_log_volume_distance,
                )
            )
            or any(
                not value.startswith("sha256:")
                for value in (
                    self.anchor_state_identity,
                    self.anchor_problem_fingerprint,
                    self.anchor_topology_fingerprint,
                    self.anchor_parameter_fingerprint,
                )
            )
        ):
            raise ValueError("reactive bubble continuation reference is incomplete")

    def to_mapping(self) -> dict[str, object]:
        return {
            "pressure_pa": self.pressure_pa,
            "vapor_mole_fractions": self.vapor_mole_fractions,
            "liquid_molar_volume_m3_per_mol": self.liquid_molar_volume_m3_per_mol,
            "vapor_molar_volume_m3_per_mol": self.vapor_molar_volume_m3_per_mol,
            "maximum_log_pressure_distance": self.maximum_log_pressure_distance,
            "maximum_log_composition_distance": self.maximum_log_composition_distance,
            "maximum_log_volume_distance": self.maximum_log_volume_distance,
            "anchor_state_identity": self.anchor_state_identity,
            "anchor_problem_fingerprint": self.anchor_problem_fingerprint,
            "anchor_topology_fingerprint": self.anchor_topology_fingerprint,
            "anchor_parameter_fingerprint": self.anchor_parameter_fingerprint,
        }

    @property
    def fingerprint(self) -> str:
        return _canonical_fingerprint(
            "epcsaft.reactive-bubble-continuation-reference", self.to_mapping()
        )


@dataclass(frozen=True, slots=True)
class ReactiveBubbleVLEProblem:
    """One source-neutral fixed-topology reactive incipient-vapor problem."""

    identity: str
    liquid_problem: HomogeneousReactiveObservationProblem
    vapor_phase_identity: str
    vapor_component_ids: tuple[str, ...]
    vapor_model: ReactiveVaporModel
    pressure_interval_pa: tuple[float, float]
    pressure_starts_pa: tuple[float, ...]
    continuation_identity: str
    continuation_reference: ReactiveBubbleContinuationReference | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "vapor_component_ids", tuple(self.vapor_component_ids))
        object.__setattr__(self, "pressure_interval_pa", tuple(self.pressure_interval_pa))
        object.__setattr__(self, "pressure_starts_pa", tuple(self.pressure_starts_pa))
        species = self.liquid_problem.reaction_system.species_ids
        charges = self.liquid_problem.reaction_system.charges
        charge_by_id = dict(zip(species, charges, strict=True))
        if (
            not self.identity
            or not self.vapor_phase_identity
            or not self.continuation_identity
            or not isinstance(self.liquid_problem.phase, ProviderPhase)
            or self.liquid_problem.phase_role != "liquid"
            or not self.vapor_component_ids
            or len(set(self.vapor_component_ids)) != len(self.vapor_component_ids)
            or any(identity not in charge_by_id for identity in self.vapor_component_ids)
            or any(charge_by_id[identity] != 0 for identity in self.vapor_component_ids)
            or len(self.pressure_interval_pa) != 2
            or not 0.0 < self.pressure_interval_pa[0] < self.pressure_interval_pa[1]
            or not self.pressure_starts_pa
            or len(set(self.pressure_starts_pa)) != len(self.pressure_starts_pa)
            or any(
                not math.isfinite(value)
                or not self.pressure_interval_pa[0] <= value <= self.pressure_interval_pa[1]
                for value in self.pressure_starts_pa
            )
        ):
            raise ValueError("reactive bubble problem topology or pressure starts are invalid")


@dataclass(frozen=True, slots=True)
class ReactiveBubbleObservationRow:
    """One coupled reactive-bubble observable in its declared measurement basis."""

    identity: str
    quantity: ReactiveBubbleQuantity
    coefficients: tuple[float, ...] = ()
    solvent_mass_coefficients_kg_per_mol: tuple[float, ...] = ()
    support: ObservationSupport = "real"
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
        if not self.identity or self.quantity not in _ROW_UNITS:
            raise ValueError("reactive bubble row identity or quantity is invalid")
        if self.support not in {"real", "positive", "left_censored", "right_censored"}:
            raise ValueError("reactive bubble observation support is invalid")
        if (self.support.endswith("censored")) != (self.censor_limit is not None):
            raise ValueError("censored reactive bubble rows require exactly one censor limit")

    @property
    def unit(self) -> str:
        return _ROW_UNITS[self.quantity]


@dataclass(frozen=True, slots=True)
class ReactiveBubbleObservationDescriptor:
    """Closed coupled observation schema and immutable topology identity."""

    identity: Literal["equilibrium.reactive-bubble-vle"]
    version: int
    component_ids: tuple[str, ...]
    charges: tuple[int, ...]
    vapor_component_ids: tuple[str, ...]
    rows: tuple[ReactiveBubbleObservationRow, ...]
    row_identities: tuple[str, ...]
    row_units: tuple[str, ...]
    active_parameter_identities: tuple[str, ...]
    active_parameter_units: tuple[str, ...]
    problem_fingerprint: str
    liquid_problem_fingerprint: str
    reaction_fingerprint: str
    topology_fingerprint: str
    reference_fingerprint: str
    domain_fingerprint: str
    parameter_fingerprint: str
    partner_parameter_fingerprint: str
    partner_component_ids: tuple[str, ...]
    active_parameter_selection_fingerprint: str
    jacobian_shape: tuple[int, int]
    local_scope: Literal["fixed_T_reactive_incipient_vapor_local_branch"]
    fingerprint: str


@dataclass(frozen=True, slots=True)
class ReactiveBubbleStartAttempt:
    ordinal: int
    start_identity: str
    start_pressure_pa: float
    status: Literal["accepted", "rejected", "failed"]
    iteration_count: int
    final_pressure_pa: float | None
    final_residual: float | None
    failure_code: NonEvaluableTrialCode | None
    diagnostic: str


@dataclass(frozen=True, slots=True)
class ReactiveBubbleSearch:
    attempts: tuple[ReactiveBubbleStartAttempt, ...]
    closure_evaluation_count: int
    provider_evaluation_count: int
    accepted_candidate_count: int
    distinct_branch_count: int
    selected_branch_ordinal: int | None
    completeness_claim: Literal["finite_declared_starts_only"]


@dataclass(frozen=True, slots=True)
class ReactiveBubblePhysicalBranch:
    branch_family_identity: str
    state_identity: str
    problem_fingerprint: str
    continuation_identity: str
    liquid_branch_family_identity: str
    liquid_state_binding_identity: str
    partner_parameter_fingerprint: str
    partner_component_ids: tuple[str, ...]
    active_trial_fingerprint: str
    pressure_pa: float
    liquid_pressure_pa: float
    liquid_mole_fractions: tuple[float, ...]
    vapor_mole_fractions: tuple[float, ...]
    liquid_molar_volume_m3_per_mol: float
    vapor_molar_volume_m3_per_mol: float
    globality_claim: Literal["not_established"]
    fingerprint: str


@dataclass(frozen=True, slots=True)
class ReactiveBubbleObservationRowResult:
    identity: str
    unit: str
    status: Literal["evaluated", "non_evaluable"]
    value: float | None
    jacobian: tuple[float, ...] | None
    failure: NonEvaluableTrial | None


@dataclass(frozen=True, slots=True)
class ReactiveBubbleCertification:
    status: Literal["certified_local_equilibrium"]
    globality_claim: Literal["not_established"]
    numerical_status: Literal["passed"]
    physical_status: Literal["passed"]
    first_order_status: Literal["passed"]
    second_order_status: Literal["passed"]
    bubble_log_closure_abs: float
    bubble_log_pressure_derivative: float
    liquid_balance_inf: float
    liquid_electroneutrality_inf: float
    liquid_reaction_affinity_inf: float
    liquid_pressure_relative_residual: float
    phase_pressure_residual_inf_pa: float
    phase_chemical_potential_residual_inf: float


@dataclass(frozen=True, slots=True)
class ReactiveBubbleObservationResult:
    descriptor: ReactiveBubbleObservationDescriptor
    status: Literal["evaluated", "non_evaluable"]
    failure: NonEvaluableTrial | None
    rows: tuple[ReactiveBubbleObservationRowResult, ...]
    trial: ActiveParameterTrialReceipt | None
    liquid_state: CertifiedHomogeneousState | None
    liquid_branch: PhysicalBranchReceipt | None
    phases: tuple[EquilibriumPhaseState, ...]
    pressure_pa: float | None
    physical_branch: ReactiveBubblePhysicalBranch | None
    certification: ReactiveBubbleCertification | None
    search: ReactiveBubbleSearch
    chemical_diagnostics: object | None
    phase_diagnostics: PhaseEquilibriumDiagnostics | None

    @property
    def values(self) -> tuple[float | None, ...]:
        return tuple(row.value for row in self.rows)

    @property
    def jacobian(self) -> tuple[tuple[float, ...] | None, ...]:
        return tuple(row.jacobian for row in self.rows)


@dataclass(frozen=True)
class _Closure:
    pressure_pa: float
    residual: float
    residual_log_pressure_derivative: float
    liquid: Any
    phases: tuple[EquilibriumPhaseState, EquilibriumPhaseState]
    boundary_pressure_pa: float
    boundary_directional_derivatives: dict[str, float]
    boundary_direct_derivatives: dict[str, float]
    vapor_fraction_directional_derivatives: dict[str, tuple[float, ...]]
    vapor_fraction_direct_derivatives: dict[str, tuple[float, ...]]
    coexistence: PhaseEquilibriumResult | None
    directions: _ReferenceDirectionSensitivity | None
    phase_diagnostics: PhaseEquilibriumDiagnostics | None
    partner_parameter_fingerprint: str
    partner_component_ids: tuple[str, ...]
    chemical_amount_directions: dict[str, tuple[float, ...]]
    chemical_volume_directions: dict[str, float]
    chemical_fraction_directions: dict[str, tuple[float, ...]]


class _BubbleFailure(RuntimeError):
    def __init__(self, code: NonEvaluableTrialCode, diagnostic: str) -> None:
        super().__init__(diagnostic)
        self.code = code


def _problem_payload(problem: ReactiveBubbleVLEProblem) -> dict[str, object]:
    liquid = problem.liquid_problem
    return {
        "identity": problem.identity,
        "liquid_problem_identity": liquid.identity,
        "liquid_phase_identity": liquid.phase_identity,
        "vapor_phase_identity": problem.vapor_phase_identity,
        "component_ids": liquid.reaction_system.species_ids,
        "charges": liquid.reaction_system.charges,
        "vapor_component_ids": problem.vapor_component_ids,
        "vapor_model": type(problem.vapor_model).__name__,
        "vapor_model_reference": problem.vapor_model.reference_id,
        "pressure_interval_pa": problem.pressure_interval_pa,
        "continuation_identity": problem.continuation_identity,
        "liquid_continuation_identity": liquid.continuation_identity,
        "liquid_continuation_reference": (
            None
            if liquid.continuation_reference is None
            else liquid.continuation_reference.fingerprint
        ),
    }


def reactive_bubble_observation_descriptor(
    problem: ReactiveBubbleVLEProblem,
    rows: tuple[ReactiveBubbleObservationRow, ...],
    *,
    active_parameters: epcsaft.ActiveParameterSet | None = None,
) -> ReactiveBubbleObservationDescriptor:
    """Describe one coupled batch without evaluating its thermodynamic state."""

    rows = tuple(rows)
    if not rows or len({row.identity for row in rows}) != len(rows):
        raise ValueError("reactive bubble rows must be nonempty with unique identities")
    species = problem.liquid_problem.reaction_system.species_ids
    vapor_count = len(problem.vapor_component_ids)
    for row in rows:
        count = vapor_count if row.quantity.startswith("vapor_") else len(species)
        if row.quantity == "bubble_pressure_pa":
            if row.coefficients:
                raise ValueError("bubble pressure row does not accept coefficients")
        elif len(row.coefficients) != count:
            raise ValueError("reactive bubble row coefficients do not match their phase support")
        if row.quantity == "liquid_molality_mol_per_kg" and len(
            row.solvent_mass_coefficients_kg_per_mol
        ) != len(species):
            raise ValueError("liquid molality row requires full solvent-mass coefficients")
        nonzero = sum(value != 0.0 for value in row.coefficients)
        if nonzero > 1 and (not row.aggregate_identity or not row.covariance_identity):
            raise ValueError("linear aggregate rows require aggregate and covariance identities")
    phase = cast(ProviderPhase, problem.liquid_problem.phase)
    liquid_contract = homogeneous_reactive_observation_descriptor(
        problem.liquid_problem,
        (
            HomogeneousReactiveObservationRow(
                "reactive-bubble-liquid-contract",
                "mole_fraction",
                (1.0,) + (0.0,) * (len(species) - 1),
            ),
        ),
        active_parameters=active_parameters,
    )
    partner_parameter_fingerprint = (
        phase.model._parameters.select(problem.vapor_component_ids).fingerprint
        if isinstance(problem.vapor_model, ProviderNonidealVapor)
        else _canonical_fingerprint(
            "epcsaft.ideal-vapor-comparison",
            {
                "reference_id": problem.vapor_model.reference_id,
                "component_ids": problem.vapor_component_ids,
            },
        )
    )
    active_ids = () if active_parameters is None else active_parameters.identities
    active_units = (
        ()
        if active_parameters is None
        else tuple(parameter.unit for parameter in active_parameters.coordinates)
    )
    topology_fingerprint = _canonical_fingerprint(
        "epcsaft.reactive-bubble-topology",
        {
            "liquid": species,
            "charges": problem.liquid_problem.reaction_system.charges,
            "vapor": problem.vapor_component_ids,
            "roles": ("liquid", "vapor"),
            "vapor_model": type(problem.vapor_model).__name__,
        },
    )
    problem_fingerprint = _canonical_fingerprint(
        "epcsaft.reactive-bubble-problem",
        {
            "declaration": _problem_payload(problem),
            "liquid_problem_fingerprint": liquid_contract.problem_fingerprint,
            "reaction_fingerprint": liquid_contract.reaction_fingerprint,
        },
    )
    reference_fingerprint = (
        "sha256:" + "0" * 64
        if problem.continuation_reference is None
        else problem.continuation_reference.fingerprint
    )
    domain_fingerprint = _canonical_fingerprint(
        "epcsaft.reactive-bubble-domain",
        {
            "pressure_interval_pa": problem.pressure_interval_pa,
            "liquid_packing": phase.admissible_packing_fraction_interval,
            "local_scope": "fixed_T_reactive_incipient_vapor_local_branch",
        },
    )
    selection_fingerprint = (
        _canonical_fingerprint("epcsaft.empty-active-selection", ())
        if active_parameters is None
        else _canonical_fingerprint(
            "epcsaft.active-parameter-selection",
            {
                "parameter_fingerprint": phase.model.parameter_fingerprint,
                "identities": active_ids,
                "units": active_units,
            },
        )
    )
    payload = {
        "problem_fingerprint": problem_fingerprint,
        "liquid_problem_fingerprint": liquid_contract.problem_fingerprint,
        "reaction_fingerprint": liquid_contract.reaction_fingerprint,
        "topology_fingerprint": topology_fingerprint,
        "reference_fingerprint": reference_fingerprint,
        "domain_fingerprint": domain_fingerprint,
        "rows": tuple(
            {
                "identity": row.identity,
                "quantity": row.quantity,
                "coefficients": row.coefficients,
                "solvent_mass_coefficients_kg_per_mol": (row.solvent_mass_coefficients_kg_per_mol),
                "support": row.support,
                "censor_limit": row.censor_limit,
                "aggregate_identity": row.aggregate_identity,
                "covariance_identity": row.covariance_identity,
            }
            for row in rows
        ),
        "active_parameter_identities": active_ids,
        "active_parameter_units": active_units,
        "partner_parameter_fingerprint": partner_parameter_fingerprint,
        "partner_component_ids": problem.vapor_component_ids,
    }
    return ReactiveBubbleObservationDescriptor(
        identity="equilibrium.reactive-bubble-vle",
        version=1,
        component_ids=species,
        charges=problem.liquid_problem.reaction_system.charges,
        vapor_component_ids=problem.vapor_component_ids,
        rows=rows,
        row_identities=tuple(row.identity for row in rows),
        row_units=tuple(row.unit for row in rows),
        active_parameter_identities=active_ids,
        active_parameter_units=active_units,
        problem_fingerprint=problem_fingerprint,
        liquid_problem_fingerprint=liquid_contract.problem_fingerprint,
        reaction_fingerprint=liquid_contract.reaction_fingerprint,
        topology_fingerprint=topology_fingerprint,
        reference_fingerprint=reference_fingerprint,
        domain_fingerprint=domain_fingerprint,
        parameter_fingerprint=phase.model.parameter_fingerprint,
        partner_parameter_fingerprint=partner_parameter_fingerprint,
        partner_component_ids=problem.vapor_component_ids,
        active_parameter_selection_fingerprint=selection_fingerprint,
        jacobian_shape=(len(rows), len(active_ids)),
        local_scope="fixed_T_reactive_incipient_vapor_local_branch",
        fingerprint=_canonical_fingerprint("epcsaft.reactive-bubble-descriptor", payload),
    )


def _failed_result(
    descriptor: ReactiveBubbleObservationDescriptor,
    code: NonEvaluableTrialCode,
    diagnostic: str,
    *,
    search: ReactiveBubbleSearch | None = None,
) -> ReactiveBubbleObservationResult:
    failure = NonEvaluableTrial(code, diagnostic)
    return ReactiveBubbleObservationResult(
        descriptor=descriptor,
        status="non_evaluable",
        failure=failure,
        rows=tuple(
            ReactiveBubbleObservationRowResult(
                row.identity, row.unit, "non_evaluable", None, None, failure
            )
            for row in descriptor.rows
        ),
        trial=None,
        liquid_state=None,
        liquid_branch=None,
        phases=(),
        pressure_pa=None,
        physical_branch=None,
        certification=None,
        search=(
            ReactiveBubbleSearch((), 0, 0, 0, 0, None, "finite_declared_starts_only")
            if search is None
            else search
        ),
        chemical_diagnostics=None,
        phase_diagnostics=None,
    )


def _chemical_directions(
    sensitivity: ChemicalEquilibriumSensitivity,
    active_ids: tuple[str, ...],
    amounts: tuple[float, ...],
) -> tuple[
    dict[str, tuple[float, ...]],
    dict[str, float],
    dict[str, tuple[float, ...]],
]:
    if sensitivity.status != "available":
        raise _BubbleFailure(
            "first_order_failure",
            sensitivity.failure_reason or "chemical sensitivity is unavailable",
        )
    requested = ("pressure_pa", *active_ids)
    amount_rows: dict[str, tuple[float, ...]] = {}
    volume_rows: dict[str, float] = {}
    for parameter, amount_row, volume in zip(
        sensitivity.parameters,
        sensitivity.amount_derivatives,
        sensitivity.volume_derivatives,
        strict=True,
    ):
        identity = "pressure_pa" if parameter.name == "pressure_pa" else parameter.identity
        if identity in requested:
            amount_rows[identity] = amount_row
            volume_rows[identity] = volume
    if tuple(amount_rows) != requested:
        raise _BubbleFailure(
            "parameter_topology_mismatch",
            "chemical sensitivity omitted or reordered a coupled input direction",
        )
    total = math.fsum(amounts)
    fraction_rows: dict[str, tuple[float, ...]] = {}
    for identity in requested:
        derivative = amount_rows[identity]
        dtotal = math.fsum(derivative)
        fraction_rows[identity] = tuple(
            (dn * total - amount * dtotal) / (total * total)
            for amount, dn in zip(amounts, derivative, strict=True)
        )
    return amount_rows, volume_rows, fraction_rows


def _liquid_trial_state(
    problem: ReactiveBubbleVLEProblem,
    temperature: Any,
    descriptor: ReactiveBubbleObservationDescriptor,
    active_parameters: epcsaft.ActiveParameterSet | None,
    pressure_pa: float,
) -> tuple[
    Any,
    Any,
    dict[str, tuple[float, ...]],
    dict[str, float],
    dict[str, tuple[float, ...]],
]:
    reaction_system = problem.liquid_problem.reaction_system
    if any(
        amount <= reaction_system.strict_interior_amount_floor_mol
        for amount in reaction_system.feed_amounts_mol
    ):
        raise _BubbleFailure(
            "active_bound_ambiguity",
            "the declared liquid feed is outside the strict positive derivative chart",
        )
    internal_row = HomogeneousReactiveObservationRow(
        identity="coupled-liquid-state-sentinel",
        quantity="mole_fraction",
        coefficients=(1.0,) + (0.0,) * (len(descriptor.component_ids) - 1),
    )
    liquid_descriptor = homogeneous_reactive_observation_descriptor(
        problem.liquid_problem,
        (internal_row,),
        active_parameters=active_parameters,
    )
    liquid = evaluate_homogeneous_reactive_observations(
        problem.liquid_problem,
        temperature,
        pressure_pa * epcsaft.unit_registry.pascal,
        liquid_descriptor,
        active_parameters=active_parameters,
    )
    if liquid.status != "evaluated" or liquid.equilibrium is None:
        if liquid.equilibrium is not None and liquid.equilibrium.diagnostics.active_trace_species:
            raise _BubbleFailure(
                "active_bound_ambiguity",
                "the reactive liquid derivative chart contains an active trace-bound species",
            )
        failure = liquid.failure
        raise _BubbleFailure(
            "numerical_convergence_failure" if failure is None else failure.code,
            "liquid reactive state is non-evaluable" if failure is None else failure.diagnostic,
        )
    equilibrium = liquid.equilibrium
    if equilibrium.sensitivity is None:
        raise _BubbleFailure("first_order_failure", "chemical pressure sensitivity is missing")
    active_ids = descriptor.active_parameter_identities
    amount_rows, volume_rows, fraction_rows = _chemical_directions(
        equilibrium.sensitivity, active_ids, equilibrium.amounts_mol
    )
    if equilibrium.diagnostics.active_trace_species:
        raise _BubbleFailure(
            "active_bound_ambiguity",
            "the reactive liquid contains a trace-bound species at the coupled state",
        )
    return liquid, equilibrium, amount_rows, volume_rows, fraction_rows


def _nonideal_closure(
    problem: ReactiveBubbleVLEProblem,
    temperature: Any,
    descriptor: ReactiveBubbleObservationDescriptor,
    active_parameters: epcsaft.ActiveParameterSet | None,
    pressure_pa: float,
) -> _Closure:
    liquid, equilibrium, amount_rows, volume_rows, fraction_rows = _liquid_trial_state(
        problem, temperature, descriptor, active_parameters, pressure_pa
    )
    active_ids = descriptor.active_parameter_identities
    topology = DeclaredTopology(
        (
            PhaseSpecification("liquid", AllComponents()),
            PhaseSpecification("vapor", DeclaredNeutralComponents(problem.vapor_component_ids)),
        )
    )
    coexistence_problem = PhaseEquilibriumProblem(
        FixedTCoexistence(temperature, topology, 0, equilibrium.mole_fractions)
    )
    direction_ids = ("pressure_pa", *active_ids)
    try:
        coexistence, directions, partner_fingerprint, partner_ids = (
            _fixed_t_coexistence_with_reference_directions(
                cast(ProviderPhase, problem.liquid_problem.phase).model,
                coexistence_problem,
                active_parameters=active_parameters,
                direction_identities=direction_ids,
                reference_composition_directions=tuple(
                    fraction_rows[identity] for identity in direction_ids
                ),
            )
        )
    except PhaseEquilibriumError as error:
        if error.diagnostics.failure_stage == "branch_selection":
            code: NonEvaluableTrialCode = "multiple_unresolved_branches"
        elif error.diagnostics.failure_stage in {
            "local_equilibrium",
            "phase_role_certificate",
        }:
            code = "no_incipient_phase"
        else:
            code = "physical_rejection"
        raise _BubbleFailure(code, str(error)) from error
    except (RuntimeError, TypeError, ValueError) as error:
        raise _BubbleFailure("parameter_topology_mismatch", str(error)) from error
    if directions.status != "available":
        reason = directions.failure_reason or "coexistence direction sensitivity unavailable"
        code = "ill_conditioned_sensitivity" if "condition" in reason else "singular_sensitivity"
        raise _BubbleFailure(cast(NonEvaluableTrialCode, code), reason)
    if (
        partner_fingerprint != descriptor.partner_parameter_fingerprint
        or partner_ids != descriptor.partner_component_ids
    ):
        raise _BubbleFailure(
            "parameter_topology_mismatch",
            "coexistence partner Provider does not match the declared vapor support",
        )
    extent = cast(IncipientBoundary, coexistence.extent)
    boundary_row = directions.output_identities.index("boundary.pressure_pa")
    direction_outputs = {
        identity: tuple(row[column] for row in directions.output_direction_jacobian)
        for column, identity in enumerate(directions.identities)
    }
    direct_outputs: dict[str, tuple[float, ...]] = {
        identity: (0.0,) * len(directions.output_identities) for identity in active_ids
    }
    if active_ids:
        sensitivity = coexistence.sensitivity
        if sensitivity is None or sensitivity.status != "available":
            raise _BubbleFailure(
                "first_order_failure", "coexistence active-parameter sensitivity is unavailable"
            )
        if tuple(parameter.identity for parameter in sensitivity.parameters) != active_ids:
            raise _BubbleFailure(
                "parameter_topology_mismatch", "coexistence active parameter order changed"
            )
        direct_outputs = {
            identity: tuple(row[column] for row in sensitivity.output_parameter_jacobian)
            for column, identity in enumerate(active_ids)
        }
    boundary_directional = {
        identity: values[boundary_row] for identity, values in direction_outputs.items()
    }
    boundary_direct = {
        identity: values[boundary_row] for identity, values in direct_outputs.items()
    }
    vapor_directional = {
        identity: tuple(
            values[directions.output_identities.index(f"phase[1].mole_fraction[{component}]")]
            for component in descriptor.vapor_component_ids
        )
        for identity, values in direction_outputs.items()
    }
    vapor_direct = {
        identity: tuple(
            values[directions.output_identities.index(f"phase[1].mole_fraction[{component}]")]
            for component in descriptor.vapor_component_ids
        )
        for identity, values in direct_outputs.items()
    }
    pressure_direction = boundary_directional["pressure_pa"]
    residual = math.log(extent.pressure_pa / pressure_pa)
    residual_derivative = pressure_pa / extent.pressure_pa * pressure_direction - 1.0
    if not math.isfinite(residual_derivative):
        raise _BubbleFailure("first_order_failure", "bubble closure derivative is non-finite")
    return _Closure(
        pressure_pa,
        residual,
        residual_derivative,
        liquid,
        cast(tuple[EquilibriumPhaseState, EquilibriumPhaseState], coexistence.phases),
        extent.pressure_pa,
        boundary_directional,
        boundary_direct,
        vapor_directional,
        vapor_direct,
        coexistence,
        directions,
        coexistence.diagnostics,
        partner_fingerprint,
        partner_ids,
        amount_rows,
        volume_rows,
        fraction_rows,
    )


def _ideal_closure(
    problem: ReactiveBubbleVLEProblem,
    temperature: Any,
    descriptor: ReactiveBubbleObservationDescriptor,
    active_parameters: epcsaft.ActiveParameterSet | None,
    pressure_pa: float,
) -> _Closure:
    liquid, equilibrium, amount_rows, volume_rows, fraction_rows = _liquid_trial_state(
        problem, temperature, descriptor, active_parameters, pressure_pa
    )
    phase = cast(ProviderPhase, problem.liquid_problem.phase)
    try:
        block = _equilibrium.evaluate_parameterized_phase(
            _native_transport(phase.model, active_parameters),
            equilibrium.temperature_k,
            equilibrium.amounts_mol,
            equilibrium.volume_m3,
            *phase.admissible_packing_fraction_interval,
        )
    except (RuntimeError, TypeError, ValueError) as error:
        raise _BubbleFailure("provider_domain_rejection", str(error)) from error
    active_ids = descriptor.active_parameter_identities
    returned_ids = tuple(
        str(value) for value in cast(Sequence[object], block["active_parameter_identities"])
    )
    if returned_ids != active_ids:
        raise _BubbleFailure(
            "parameter_topology_mismatch", "ideal comparison Provider parameter order changed"
        )
    amounts = equilibrium.amounts_mol
    component_count = len(amounts)
    state_count = component_count + 1
    state_gradient = tuple(
        float(cast(Any, value)) for value in cast(Sequence[object], block["state_gradient"])
    )
    state_hessian = tuple(
        float(cast(Any, value)) for value in cast(Sequence[object], block["state_hessian"])
    )
    direct_mu = tuple(
        float(cast(Any, value))
        for value in cast(
            Sequence[object], block["chemical_potential_parameter_derivatives_over_rt"]
        )
    )
    if (
        len(state_gradient) != state_count
        or len(state_hessian) != state_count * state_count
        or len(direct_mu) != component_count * len(active_ids)
    ):
        raise _BubbleFailure("owner_contract_failure", "ideal comparison derivative layout changed")
    component_index = {identity: index for index, identity in enumerate(descriptor.component_ids)}
    vapor_indices = tuple(component_index[identity] for identity in problem.vapor_component_ids)
    rt = 8.31446261815324 * equilibrium.temperature_k
    fugacities = tuple(rt * math.exp(state_gradient[index]) for index in vapor_indices)
    boundary_pressure = math.fsum(fugacities)
    if not math.isfinite(boundary_pressure) or boundary_pressure <= 0.0:
        raise _BubbleFailure(
            "no_incipient_phase", "ideal comparison has no positive incipient pressure"
        )
    direction_ids = ("pressure_pa", *active_ids)
    boundary_directional: dict[str, float] = {}
    vapor_directional: dict[str, tuple[float, ...]] = {}
    for identity in direction_ids:
        mu_direction = tuple(
            math.fsum(
                state_hessian[component * state_count + species] * derivative
                for species, derivative in enumerate(amount_rows[identity])
            )
            + state_hessian[component * state_count + component_count] * volume_rows[identity]
            for component in vapor_indices
        )
        boundary_derivative = math.fsum(
            fugacity * derivative
            for fugacity, derivative in zip(fugacities, mu_direction, strict=True)
        )
        boundary_directional[identity] = boundary_derivative
        vapor_directional[identity] = tuple(
            fugacity
            * (derivative * boundary_pressure - boundary_derivative)
            / (boundary_pressure * boundary_pressure)
            for fugacity, derivative in zip(fugacities, mu_direction, strict=True)
        )
    boundary_direct: dict[str, float] = {}
    vapor_direct: dict[str, tuple[float, ...]] = {}
    for parameter, identity in enumerate(active_ids):
        mu_partial = tuple(
            direct_mu[component * len(active_ids) + parameter] for component in vapor_indices
        )
        boundary_derivative = math.fsum(
            fugacity * derivative
            for fugacity, derivative in zip(fugacities, mu_partial, strict=True)
        )
        boundary_direct[identity] = boundary_derivative
        vapor_direct[identity] = tuple(
            fugacity
            * (derivative * boundary_pressure - boundary_derivative)
            / (boundary_pressure * boundary_pressure)
            for fugacity, derivative in zip(fugacities, mu_partial, strict=True)
        )
    total = math.fsum(amounts)
    molar_volume = equilibrium.volume_m3 / total
    liquid_phase = EquilibriumPhaseState(
        role="liquid",
        support="all_components",
        supported_component_ids=descriptor.component_ids,
        mole_fractions=equilibrium.mole_fractions,
        molar_volume_m3_per_mol=molar_volume,
        molar_density_mol_m3=1.0 / molar_volume,
        packing_fraction=float(cast(Any, block["packing_fraction"])),
        pressure_pa=float(cast(Any, block["pressure_pa"])),
        chemical_potential_component_ids=descriptor.component_ids,
        chemical_potential_over_rt=state_gradient[:component_count],
        mechanical_class=liquid.state.phase_role_receipt.mechanical_class,
    )
    vapor_fractions = tuple(value / boundary_pressure for value in fugacities)
    vapor_volume = rt / boundary_pressure
    vapor_phase = EquilibriumPhaseState(
        role="vapor",
        support="declared_neutral_components",
        supported_component_ids=descriptor.vapor_component_ids,
        mole_fractions=vapor_fractions,
        molar_volume_m3_per_mol=vapor_volume,
        molar_density_mol_m3=1.0 / vapor_volume,
        packing_fraction=0.0,
        pressure_pa=boundary_pressure,
        chemical_potential_component_ids=descriptor.vapor_component_ids,
        chemical_potential_over_rt=tuple(state_gradient[index] for index in vapor_indices),
        mechanical_class="strict_stable_ideal_gas",
    )
    residual = math.log(boundary_pressure / pressure_pa)
    residual_derivative = (
        pressure_pa / boundary_pressure * boundary_directional["pressure_pa"] - 1.0
    )
    if not math.isfinite(residual_derivative):
        raise _BubbleFailure("first_order_failure", "ideal bubble closure derivative is non-finite")
    return _Closure(
        pressure_pa,
        residual,
        residual_derivative,
        liquid,
        (liquid_phase, vapor_phase),
        boundary_pressure,
        boundary_directional,
        boundary_direct,
        vapor_directional,
        vapor_direct,
        None,
        None,
        None,
        descriptor.partner_parameter_fingerprint,
        descriptor.partner_component_ids,
        amount_rows,
        volume_rows,
        fraction_rows,
    )


def _closure(
    problem: ReactiveBubbleVLEProblem,
    temperature: Any,
    descriptor: ReactiveBubbleObservationDescriptor,
    active_parameters: epcsaft.ActiveParameterSet | None,
    pressure_pa: float,
) -> _Closure:
    if isinstance(problem.vapor_model, IdealGasVaporComparison):
        return _ideal_closure(problem, temperature, descriptor, active_parameters, pressure_pa)
    return _nonideal_closure(problem, temperature, descriptor, active_parameters, pressure_pa)


def _candidate_distance(
    left: _Closure,
    right: _Closure,
) -> float:
    left_vapor = left.phases[1]
    right_vapor = right.phases[1]
    return max(
        abs(math.log(left.pressure_pa / right.pressure_pa)),
        max(
            abs(math.log(a / b))
            for a, b in zip(left_vapor.mole_fractions, right_vapor.mole_fractions, strict=True)
        ),
    )


def _validate_same_liquid_state(closure: _Closure) -> None:
    """Bind the converged coexistence liquid root to the exact chemical T,n,V state."""

    equilibrium = closure.liquid.equilibrium
    liquid_state = closure.liquid.state
    liquid_branch = closure.liquid.physical_branch
    liquid_trial = closure.liquid.trial
    liquid_phase = closure.phases[0]
    if equilibrium is None:
        raise _BubbleFailure("owner_contract_failure", "chemical result is missing")
    chemical_molar_volume = equilibrium.volume_m3 / math.fsum(equilibrium.amounts_mol)
    if (
        liquid_state is None
        or liquid_branch is None
        or liquid_trial is None
        or liquid_phase.role != "liquid"
        or liquid_phase.support != "all_components"
        or liquid_phase.supported_component_ids != equilibrium.species_ids
        or any(
            abs(actual - expected) > 2.0e-13
            for actual, expected in zip(
                liquid_phase.mole_fractions, equilibrium.mole_fractions, strict=True
            )
        )
        or not math.isclose(
            liquid_phase.molar_volume_m3_per_mol,
            chemical_molar_volume,
            rel_tol=2.0e-9,
            abs_tol=2.0e-14,
        )
        or not math.isclose(
            liquid_phase.pressure_pa,
            equilibrium.pressure_pa,
            rel_tol=2.0e-9,
            abs_tol=2.0e-5,
        )
        or liquid_state.state_identity != liquid_branch.state_identity
        or liquid_trial.parameter_fingerprint != closure.liquid.descriptor.parameter_fingerprint
        or (
            closure.coexistence is not None
            and liquid_trial.parameter_fingerprint != closure.coexistence.parameter_fingerprint
        )
    ):
        raise _BubbleFailure(
            "branch_loss",
            "converged coexistence liquid root does not bind to the certified chemical branch",
        )


def _solve_starts(
    problem: ReactiveBubbleVLEProblem,
    temperature: Any,
    descriptor: ReactiveBubbleObservationDescriptor,
    active_parameters: epcsaft.ActiveParameterSet | None,
) -> tuple[tuple[_Closure, ...], ReactiveBubbleSearch]:
    lower, upper = problem.pressure_interval_pa
    starts = (
        (
            _canonical_fingerprint(
                "epcsaft.reactive-bubble-pressure-start",
                {"problem_identity": problem.identity, "kind": "interval_lower", "value": lower},
            ),
            lower,
        ),
        *(
            (
                _canonical_fingerprint(
                    "epcsaft.reactive-bubble-pressure-start",
                    {
                        "problem_identity": problem.identity,
                        "kind": "declared",
                        "value": value,
                    },
                ),
                value,
            )
            for value in problem.pressure_starts_pa
        ),
        (
            _canonical_fingerprint(
                "epcsaft.reactive-bubble-pressure-start",
                {"problem_identity": problem.identity, "kind": "interval_upper", "value": upper},
            ),
            upper,
        ),
    )
    cache: dict[float, _Closure] = {}

    def evaluate(pressure: float) -> _Closure:
        if pressure not in cache:
            cache[pressure] = _closure(
                problem, temperature, descriptor, active_parameters, pressure
            )
        return cache[pressure]

    attempts: list[ReactiveBubbleStartAttempt] = []
    accepted: list[_Closure] = []
    for ordinal, (start_identity, start) in enumerate(starts):
        pressure = start
        final: _Closure | None = None
        failure: _BubbleFailure | None = None
        iteration_count = 0
        for iteration in range(18):
            iteration_count = iteration + 1
            try:
                current = evaluate(pressure)
            except _BubbleFailure as error:
                failure = error
                break
            if abs(current.residual) <= 2.0e-10:
                final = current
                break
            slope = current.residual_log_pressure_derivative
            if abs(slope) < 1.0e-10:
                failure = _BubbleFailure("singular_sensitivity", "bubble closure slope is singular")
                break
            step = max(-1.0, min(1.0, -current.residual / slope))
            candidate_log = min(
                math.log(upper),
                max(math.log(lower), math.log(pressure) + step),
            )
            candidate = math.exp(candidate_log)
            if candidate == pressure:
                failure = _BubbleFailure(
                    "no_incipient_phase", "declared pressure interval excludes the Newton step"
                )
                break
            pressure = candidate
        if final is not None:
            attempts.append(
                ReactiveBubbleStartAttempt(
                    ordinal,
                    start_identity,
                    start,
                    "accepted",
                    iteration_count,
                    final.pressure_pa,
                    final.residual,
                    None,
                    "",
                )
            )
            if not any(_candidate_distance(final, existing) <= 2.0e-8 for existing in accepted):
                accepted.append(final)
        else:
            code = "numerical_convergence_failure" if failure is None else failure.code
            diagnostic = (
                "bubble Newton iteration limit reached" if failure is None else str(failure)
            )
            attempts.append(
                ReactiveBubbleStartAttempt(
                    ordinal,
                    start_identity,
                    start,
                    "rejected" if failure is not None else "failed",
                    iteration_count,
                    None,
                    None,
                    code,
                    diagnostic,
                )
            )
    search = ReactiveBubbleSearch(
        attempts=tuple(attempts),
        closure_evaluation_count=len(cache),
        provider_evaluation_count=sum(
            candidate.liquid.state.phase_role_receipt.provider_evaluation_count
            + (
                1
                if candidate.phase_diagnostics is None
                else candidate.phase_diagnostics.provider_evaluation_count
            )
            for candidate in cache.values()
        ),
        accepted_candidate_count=sum(attempt.status == "accepted" for attempt in attempts),
        distinct_branch_count=len(accepted),
        selected_branch_ordinal=None,
        completeness_claim="finite_declared_starts_only",
    )
    return tuple(accepted), search


def _select_branch(
    problem: ReactiveBubbleVLEProblem,
    descriptor: ReactiveBubbleObservationDescriptor,
    candidates: tuple[_Closure, ...],
    search: ReactiveBubbleSearch,
) -> tuple[int, _Closure]:
    if not candidates:
        priority: tuple[NonEvaluableTrialCode, ...] = (
            "multiple_unresolved_branches",
            "branch_loss",
            "active_bound_ambiguity",
            "component_order_mismatch",
            "parameter_topology_mismatch",
            "provider_domain_rejection",
            "source_reference_derivative_unavailable",
            "source_reference_limit_failure",
            "source_reference_failure",
            "ill_conditioned_sensitivity",
            "singular_sensitivity",
            "first_order_failure",
            "physical_rejection",
            "no_incipient_phase",
            "numerical_convergence_failure",
        )
        failures = {
            attempt.failure_code: attempt.diagnostic
            for attempt in search.attempts
            if attempt.failure_code is not None
        }
        for code in priority:
            if code in failures:
                raise _BubbleFailure(code, failures[code])
        raise _BubbleFailure("no_incipient_phase", "no declared pressure start found a branch")
    reference = problem.continuation_reference
    if reference is None:
        if len(candidates) != 1:
            raise _BubbleFailure(
                "multiple_unresolved_branches", "declared starts found multiple bubble branches"
            )
        return 0, candidates[0]
    first_trial = candidates[0].liquid.trial
    if (
        first_trial is None
        or reference.anchor_problem_fingerprint != descriptor.problem_fingerprint
        or reference.anchor_topology_fingerprint != descriptor.topology_fingerprint
        or reference.anchor_parameter_fingerprint != first_trial.parameter_fingerprint
    ):
        raise _BubbleFailure(
            "parameter_topology_mismatch",
            "reactive bubble continuation anchor does not match this trial topology",
        )
    admitted: list[tuple[int, _Closure]] = []
    for ordinal, candidate in enumerate(candidates):
        liquid = candidate.phases[0]
        vapor = candidate.phases[1]
        pressure_distance = abs(math.log(candidate.pressure_pa / reference.pressure_pa))
        composition_distance = max(
            abs(math.log(value / anchor))
            for value, anchor in zip(
                vapor.mole_fractions, reference.vapor_mole_fractions, strict=True
            )
        )
        volume_distance = max(
            abs(
                math.log(liquid.molar_volume_m3_per_mol / reference.liquid_molar_volume_m3_per_mol)
            ),
            abs(math.log(vapor.molar_volume_m3_per_mol / reference.vapor_molar_volume_m3_per_mol)),
        )
        if (
            pressure_distance <= reference.maximum_log_pressure_distance
            and composition_distance <= reference.maximum_log_composition_distance
            and volume_distance <= reference.maximum_log_volume_distance
        ):
            admitted.append((ordinal, candidate))
    if not admitted:
        raise _BubbleFailure("branch_loss", "no bubble branch remains inside the reference region")
    if len(admitted) != 1:
        raise _BubbleFailure(
            "multiple_unresolved_branches",
            "multiple bubble branches remain inside the reference region",
        )
    return admitted[0]


def _row_values(
    descriptor: ReactiveBubbleObservationDescriptor,
    closure: _Closure,
) -> tuple[ReactiveBubbleObservationRowResult, ...]:
    active_ids = descriptor.active_parameter_identities
    pressure_partial = closure.boundary_directional_derivatives["pressure_pa"]
    schur_pressure = pressure_partial - 1.0
    schur_scale = max(1.0, abs(pressure_partial))
    if not math.isfinite(schur_pressure) or abs(schur_pressure) < 1.0e-12 * schur_scale:
        raise _BubbleFailure("singular_sensitivity", "coupled bubble Schur complement is singular")
    if abs(schur_pressure) < 1.0e-8 * schur_scale:
        raise _BubbleFailure(
            "ill_conditioned_sensitivity",
            "coupled bubble Schur complement is ill-conditioned",
        )
    pressure_derivatives: dict[str, float] = {}
    for identity in active_ids:
        fixed_pressure_parameter_derivative = (
            closure.boundary_direct_derivatives[identity]
            + closure.boundary_directional_derivatives[identity]
        )
        pressure_derivatives[identity] = -fixed_pressure_parameter_derivative / schur_pressure
    equilibrium = closure.liquid.equilibrium
    assert equilibrium is not None
    amounts = equilibrium.amounts_mol
    volume = equilibrium.volume_m3
    total = math.fsum(amounts)
    liquid_amount_total: dict[str, tuple[float, ...]] = {}
    liquid_volume_total: dict[str, float] = {}
    liquid_fraction_total: dict[str, tuple[float, ...]] = {}
    for identity in active_ids:
        pressure_derivative = pressure_derivatives[identity]
        amount_row = tuple(
            direct + pressure * pressure_derivative
            for direct, pressure in zip(
                closure.chemical_amount_directions[identity],
                closure.chemical_amount_directions["pressure_pa"],
                strict=True,
            )
        )
        liquid_amount_total[identity] = amount_row
        liquid_volume_total[identity] = (
            closure.chemical_volume_directions[identity]
            + closure.chemical_volume_directions["pressure_pa"] * pressure_derivative
        )
        dtotal = math.fsum(amount_row)
        liquid_fraction_total[identity] = tuple(
            (dn * total - amount * dtotal) / (total * total)
            for amount, dn in zip(amounts, amount_row, strict=True)
        )
    vapor = closure.phases[1]
    vapor_derivatives: dict[str, tuple[float, ...]] = {}
    for identity in active_ids:
        vapor_derivatives[identity] = tuple(
            direct + directional + pressure * pressure_derivatives[identity]
            for direct, directional, pressure in zip(
                closure.vapor_fraction_direct_derivatives[identity],
                closure.vapor_fraction_directional_derivatives[identity],
                closure.vapor_fraction_directional_derivatives["pressure_pa"],
                strict=True,
            )
        )
    results: list[ReactiveBubbleObservationRowResult] = []
    for row in descriptor.rows:
        if row.quantity == "bubble_pressure_pa":
            value = closure.pressure_pa
            jacobian = tuple(pressure_derivatives[identity] for identity in active_ids)
        elif row.quantity == "vapor_mole_fraction":
            value = math.fsum(
                coefficient * fraction
                for coefficient, fraction in zip(
                    row.coefficients, vapor.mole_fractions, strict=True
                )
            )
            jacobian = tuple(
                math.fsum(
                    coefficient * derivative
                    for coefficient, derivative in zip(
                        row.coefficients, vapor_derivatives[identity], strict=True
                    )
                )
                for identity in active_ids
            )
        elif row.quantity == "vapor_partial_pressure_pa":
            vapor_fraction = math.fsum(
                coefficient * fraction
                for coefficient, fraction in zip(
                    row.coefficients, vapor.mole_fractions, strict=True
                )
            )
            value = closure.pressure_pa * vapor_fraction
            jacobian = tuple(
                pressure_derivatives[identity] * vapor_fraction
                + closure.pressure_pa
                * math.fsum(
                    coefficient * derivative
                    for coefficient, derivative in zip(
                        row.coefficients, vapor_derivatives[identity], strict=True
                    )
                )
                for identity in active_ids
            )
        elif row.quantity == "liquid_amount_mol":
            value = math.fsum(
                coefficient * amount
                for coefficient, amount in zip(row.coefficients, amounts, strict=True)
            )
            jacobian = tuple(
                math.fsum(
                    coefficient * derivative
                    for coefficient, derivative in zip(
                        row.coefficients, liquid_amount_total[identity], strict=True
                    )
                )
                for identity in active_ids
            )
        elif row.quantity == "liquid_mole_fraction":
            value = math.fsum(
                coefficient * fraction
                for coefficient, fraction in zip(
                    row.coefficients, equilibrium.mole_fractions, strict=True
                )
            )
            jacobian = tuple(
                math.fsum(
                    coefficient * derivative
                    for coefficient, derivative in zip(
                        row.coefficients, liquid_fraction_total[identity], strict=True
                    )
                )
                for identity in active_ids
            )
        elif row.quantity == "liquid_molality_mol_per_kg":
            numerator = math.fsum(
                coefficient * amount
                for coefficient, amount in zip(row.coefficients, amounts, strict=True)
            )
            denominator_mass = math.fsum(
                coefficient * amount
                for coefficient, amount in zip(
                    row.solvent_mass_coefficients_kg_per_mol, amounts, strict=True
                )
            )
            if denominator_mass <= 0.0:
                raise _BubbleFailure(
                    "observation_out_of_support", "liquid molality denominator is not positive"
                )
            value = numerator / denominator_mass
            jacobian_values = []
            for identity in active_ids:
                dn = liquid_amount_total[identity]
                dnumerator = math.fsum(
                    coefficient * derivative
                    for coefficient, derivative in zip(row.coefficients, dn, strict=True)
                )
                ddenominator = math.fsum(
                    coefficient * derivative
                    for coefficient, derivative in zip(
                        row.solvent_mass_coefficients_kg_per_mol, dn, strict=True
                    )
                )
                jacobian_values.append(
                    (dnumerator * denominator_mass - numerator * ddenominator)
                    / (denominator_mass * denominator_mass)
                )
            jacobian = tuple(jacobian_values)
        else:
            value = total / volume
            jacobian = tuple(
                math.fsum(liquid_amount_total[identity]) / volume
                - total * liquid_volume_total[identity] / (volume * volume)
                for identity in active_ids
            )
        if not math.isfinite(value) or any(not math.isfinite(item) for item in jacobian):
            raise _BubbleFailure(
                "owner_contract_failure", "reactive bubble row returned non-finite data"
            )
        if row.support == "positive" and value <= 0.0:
            raise _BubbleFailure(
                "observation_out_of_support", "positive reactive bubble row is not positive"
            )
        results.append(
            ReactiveBubbleObservationRowResult(
                row.identity, row.unit, "evaluated", value, jacobian, None
            )
        )
    return tuple(results)


def _physical_branch(
    problem: ReactiveBubbleVLEProblem,
    descriptor: ReactiveBubbleObservationDescriptor,
    closure: _Closure,
) -> ReactiveBubblePhysicalBranch:
    liquid_state = closure.liquid.state
    liquid_branch = closure.liquid.physical_branch
    liquid_trial = closure.liquid.trial
    if liquid_state is None or liquid_branch is None or liquid_trial is None:
        raise _BubbleFailure("owner_contract_failure", "liquid branch receipt is missing")
    liquid_phase, vapor_phase = closure.phases
    branch_family = _canonical_fingerprint(
        "epcsaft.reactive-bubble-branch-family",
        {
            "problem_fingerprint": descriptor.problem_fingerprint,
            "topology_fingerprint": descriptor.topology_fingerprint,
            "reference_fingerprint": descriptor.reference_fingerprint,
            "continuation_identity": problem.continuation_identity,
            "liquid_branch_family_identity": liquid_branch.branch_family_identity,
        },
    )
    liquid_state_binding_identity = _canonical_fingerprint(
        "epcsaft.reactive-bubble-liquid-state-binding",
        {
            "component_ids": liquid_phase.supported_component_ids,
            "mole_fractions": liquid_phase.mole_fractions,
            "molar_volume_m3_per_mol": liquid_phase.molar_volume_m3_per_mol,
            "pressure_pa": liquid_phase.pressure_pa,
            "liquid_state_identity": liquid_state.state_identity,
            "liquid_branch_family_identity": liquid_branch.branch_family_identity,
            "parameter_fingerprint": liquid_trial.parameter_fingerprint,
        },
    )
    state_identity = _canonical_fingerprint(
        "epcsaft.reactive-bubble-state",
        {
            "pressure_pa": closure.pressure_pa,
            "liquid_state_identity": liquid_state.state_identity,
            "liquid_mole_fractions": liquid_phase.mole_fractions,
            "vapor_mole_fractions": vapor_phase.mole_fractions,
            "liquid_molar_volume_m3_per_mol": liquid_phase.molar_volume_m3_per_mol,
            "vapor_molar_volume_m3_per_mol": vapor_phase.molar_volume_m3_per_mol,
            "topology_fingerprint": descriptor.topology_fingerprint,
        },
    )
    payload = {
        "branch_family_identity": branch_family,
        "state_identity": state_identity,
        "problem_fingerprint": descriptor.problem_fingerprint,
        "continuation_identity": problem.continuation_identity,
        "liquid_branch_family_identity": liquid_branch.branch_family_identity,
        "liquid_state_binding_identity": liquid_state_binding_identity,
        "partner_parameter_fingerprint": closure.partner_parameter_fingerprint,
        "partner_component_ids": closure.partner_component_ids,
        "active_trial_fingerprint": liquid_trial.fingerprint,
        "pressure_pa": closure.pressure_pa,
        "liquid_pressure_pa": liquid_phase.pressure_pa,
        "liquid_mole_fractions": liquid_phase.mole_fractions,
        "vapor_mole_fractions": vapor_phase.mole_fractions,
        "liquid_molar_volume_m3_per_mol": liquid_phase.molar_volume_m3_per_mol,
        "vapor_molar_volume_m3_per_mol": vapor_phase.molar_volume_m3_per_mol,
        "globality_claim": "not_established",
    }
    return ReactiveBubblePhysicalBranch(
        **payload,
        fingerprint=_canonical_fingerprint("epcsaft.reactive-bubble-physical-branch", payload),
    )


def _certification(closure: _Closure) -> ReactiveBubbleCertification:
    diagnostics = closure.liquid.equilibrium.diagnostics
    liquid_phase, vapor_phase = closure.phases
    potential_residual = max(
        abs(left - right)
        for left, right in zip(
            (
                liquid_phase.chemical_potential_over_rt[
                    closure.liquid.equilibrium.species_ids.index(identity)
                ]
                for identity in vapor_phase.chemical_potential_component_ids
            ),
            vapor_phase.chemical_potential_over_rt,
            strict=True,
        )
    )
    phase_pressure_residual = (
        abs(liquid_phase.pressure_pa - vapor_phase.pressure_pa)
        if closure.phase_diagnostics is None
        else closure.phase_diagnostics.pressure_residual_inf
    )
    if diagnostics.local_minimum_status != "passed" or any(
        "strict_stable" not in phase.mechanical_class for phase in closure.phases
    ):
        raise _BubbleFailure(
            "second_order_failure",
            "coupled bubble state lacks strict local chemical/mechanical stability",
        )
    required = (
        diagnostics.balance_inf_norm,
        diagnostics.charge_inf_norm,
        diagnostics.reaction_affinity_inf_norm,
        diagnostics.pressure_relative_residual,
    )
    if any(value is None or not math.isfinite(value) for value in required):
        raise _BubbleFailure("physical_rejection", "coupled physical residuals are incomplete")
    if (
        abs(closure.residual) > 2.0e-10
        or abs(closure.residual_log_pressure_derivative) < 1.0e-10
        or not math.isfinite(potential_residual)
        or not math.isfinite(phase_pressure_residual)
    ):
        raise _BubbleFailure("first_order_failure", "coupled bubble closure is uncertified")
    return ReactiveBubbleCertification(
        "certified_local_equilibrium",
        "not_established",
        "passed",
        "passed",
        "passed",
        "passed",
        abs(closure.residual),
        closure.residual_log_pressure_derivative,
        cast(float, diagnostics.balance_inf_norm),
        cast(float, diagnostics.charge_inf_norm),
        cast(float, diagnostics.reaction_affinity_inf_norm),
        cast(float, diagnostics.pressure_relative_residual),
        phase_pressure_residual,
        potential_residual,
    )


def evaluate_reactive_bubble_observations(
    problem: ReactiveBubbleVLEProblem,
    temperature: Any,
    descriptor: ReactiveBubbleObservationDescriptor,
    *,
    active_parameters: epcsaft.ActiveParameterSet | None = None,
) -> ReactiveBubbleObservationResult:
    """Solve and evaluate one pure, locally certified coupled bubble batch."""

    try:
        expected = reactive_bubble_observation_descriptor(
            problem, descriptor.rows, active_parameters=active_parameters
        )
    except (TypeError, ValueError) as error:
        return _failed_result(descriptor, "parameter_topology_mismatch", str(error))
    if expected.fingerprint != descriptor.fingerprint:
        return _failed_result(
            descriptor, "parameter_topology_mismatch", "reactive bubble descriptor changed"
        )
    candidates, search = _solve_starts(problem, temperature, descriptor, active_parameters)
    try:
        selected_ordinal, closure = _select_branch(problem, descriptor, candidates, search)
        search = replace(search, selected_branch_ordinal=selected_ordinal)
        _validate_same_liquid_state(closure)
        rows = _row_values(descriptor, closure)
        branch = _physical_branch(problem, descriptor, closure)
        certification = _certification(closure)
    except _BubbleFailure as error:
        return _failed_result(descriptor, error.code, str(error), search=search)
    return ReactiveBubbleObservationResult(
        descriptor=descriptor,
        status="evaluated",
        failure=None,
        rows=rows,
        trial=closure.liquid.trial,
        liquid_state=closure.liquid.state,
        liquid_branch=closure.liquid.physical_branch,
        phases=closure.phases,
        pressure_pa=closure.pressure_pa,
        physical_branch=branch,
        certification=certification,
        search=search,
        chemical_diagnostics=closure.liquid.diagnostics,
        phase_diagnostics=closure.phase_diagnostics,
    )


def certify_reactive_bubble_continuation_reference(
    problem: ReactiveBubbleVLEProblem,
    temperature: Any,
    *,
    active_parameters: epcsaft.ActiveParameterSet | None = None,
    maximum_log_pressure_distance: float,
    maximum_log_composition_distance: float,
    maximum_log_volume_distance: float,
) -> ReactiveBubbleContinuationReference:
    """Create an immutable local coupled-branch anchor from one unambiguous solve."""

    if problem.continuation_reference is not None:
        raise ValueError(
            "reactive bubble continuation certification requires an unanchored problem"
        )
    row = ReactiveBubbleObservationRow("bubble-pressure-anchor", "bubble_pressure_pa")
    descriptor = reactive_bubble_observation_descriptor(
        problem, (row,), active_parameters=active_parameters
    )
    result = evaluate_reactive_bubble_observations(
        problem, temperature, descriptor, active_parameters=active_parameters
    )
    if (
        result.status != "evaluated"
        or result.pressure_pa is None
        or result.physical_branch is None
        or result.certification is None
        or result.trial is None
        or len(result.phases) != 2
    ):
        detail = "reactive bubble anchor is non-evaluable"
        if result.failure is not None:
            detail = result.failure.diagnostic
        raise ValueError(detail)
    liquid, vapor = result.phases
    return ReactiveBubbleContinuationReference(
        pressure_pa=result.pressure_pa,
        vapor_mole_fractions=vapor.mole_fractions,
        liquid_molar_volume_m3_per_mol=liquid.molar_volume_m3_per_mol,
        vapor_molar_volume_m3_per_mol=vapor.molar_volume_m3_per_mol,
        maximum_log_pressure_distance=maximum_log_pressure_distance,
        maximum_log_composition_distance=maximum_log_composition_distance,
        maximum_log_volume_distance=maximum_log_volume_distance,
        anchor_state_identity=result.physical_branch.state_identity,
        anchor_problem_fingerprint=descriptor.problem_fingerprint,
        anchor_topology_fingerprint=descriptor.topology_fingerprint,
        anchor_parameter_fingerprint=result.trial.parameter_fingerprint,
    )


def serialize_reactive_bubble_owner_receipt(
    result: ReactiveBubbleObservationResult,
    *,
    batch_identity: str,
) -> str:
    """Serialize one evaluated CAP-13 owner result as canonical opaque bytes."""

    if (
        result.status != "evaluated"
        or result.trial is None
        or result.liquid_state is None
        or result.liquid_branch is None
        or result.physical_branch is None
        or result.certification is None
        or result.pressure_pa is None
    ):
        raise ValueError("only evaluated reactive bubble results have owner receipts")
    certification = result.certification
    assert certification is not None
    payload = {
        "schema": "epcsaft.reactive-bubble-owner-receipt.v1",
        "batch_identity": batch_identity,
        "descriptor_fingerprint": result.descriptor.fingerprint,
        "problem_fingerprint": result.descriptor.problem_fingerprint,
        "topology_fingerprint": result.descriptor.topology_fingerprint,
        "reference_fingerprint": result.descriptor.reference_fingerprint,
        "trial": {
            "parameter_fingerprint": result.trial.parameter_fingerprint,
            "identities": result.trial.identities,
            "values": result.trial.values,
            "units": result.trial.units,
            "active_parameter_fingerprint": result.trial.active_parameter_fingerprint,
            "fingerprint": result.trial.fingerprint,
        },
        "pressure_pa": result.pressure_pa,
        "temperature_k": result.liquid_state.temperature_k,
        "certification": {
            "status": certification.status,
            "globality_claim": certification.globality_claim,
            "numerical_status": certification.numerical_status,
            "physical_status": certification.physical_status,
            "first_order_status": certification.first_order_status,
            "second_order_status": certification.second_order_status,
            "bubble_log_closure_abs": certification.bubble_log_closure_abs,
            "bubble_log_pressure_derivative": certification.bubble_log_pressure_derivative,
            "liquid_balance_inf": certification.liquid_balance_inf,
            "liquid_electroneutrality_inf": certification.liquid_electroneutrality_inf,
            "liquid_reaction_affinity_inf": certification.liquid_reaction_affinity_inf,
            "liquid_pressure_relative_residual": (certification.liquid_pressure_relative_residual),
            "phase_pressure_residual_inf_pa": certification.phase_pressure_residual_inf_pa,
            "phase_chemical_potential_residual_inf": (
                certification.phase_chemical_potential_residual_inf
            ),
        },
        "liquid_state_identity": result.liquid_state.state_identity,
        "liquid_branch_family_identity": result.liquid_branch.branch_family_identity,
        "physical_branch": {
            "branch_family_identity": result.physical_branch.branch_family_identity,
            "state_identity": result.physical_branch.state_identity,
            "problem_fingerprint": result.physical_branch.problem_fingerprint,
            "continuation_identity": result.physical_branch.continuation_identity,
            "liquid_branch_family_identity": (result.physical_branch.liquid_branch_family_identity),
            "fingerprint": result.physical_branch.fingerprint,
            "globality_claim": result.physical_branch.globality_claim,
            "liquid_state_binding_identity": (result.physical_branch.liquid_state_binding_identity),
            "partner_parameter_fingerprint": (result.physical_branch.partner_parameter_fingerprint),
            "partner_component_ids": result.physical_branch.partner_component_ids,
            "active_trial_fingerprint": result.physical_branch.active_trial_fingerprint,
            "pressure_pa": result.physical_branch.pressure_pa,
            "liquid_pressure_pa": result.physical_branch.liquid_pressure_pa,
            "liquid_mole_fractions": result.physical_branch.liquid_mole_fractions,
            "vapor_mole_fractions": result.physical_branch.vapor_mole_fractions,
            "liquid_molar_volume_m3_per_mol": (
                result.physical_branch.liquid_molar_volume_m3_per_mol
            ),
            "vapor_molar_volume_m3_per_mol": (result.physical_branch.vapor_molar_volume_m3_per_mol),
        },
        "search": {
            "closure_evaluation_count": result.search.closure_evaluation_count,
            "provider_evaluation_count": result.search.provider_evaluation_count,
            "accepted_candidate_count": result.search.accepted_candidate_count,
            "distinct_branch_count": result.search.distinct_branch_count,
            "selected_branch_ordinal": result.search.selected_branch_ordinal,
            "completeness_claim": result.search.completeness_claim,
            "attempts": tuple(
                {
                    "ordinal": attempt.ordinal,
                    "start_identity": attempt.start_identity,
                    "start_pressure_pa": attempt.start_pressure_pa,
                    "status": attempt.status,
                    "iteration_count": attempt.iteration_count,
                    "final_pressure_pa": attempt.final_pressure_pa,
                    "final_residual": attempt.final_residual,
                    "failure_code": attempt.failure_code,
                    "diagnostic": attempt.diagnostic,
                }
                for attempt in result.search.attempts
            ),
        },
        "observations": tuple(
            {
                "identity": row.identity,
                "unit": row.unit,
                "value": row.value,
                "jacobian": row.jacobian,
            }
            for row in result.rows
        ),
    }
    return _canonical_json(payload)


def validate_reactive_bubble_owner_receipt(
    receipt_json: str,
    *,
    descriptor: ReactiveBubbleObservationDescriptor,
    batch_identity: str,
    requested_temperature: Any,
    active_parameters: epcsaft.ActiveParameterSet,
) -> str:
    """Validate canonical CAP-13 owner bytes and return their digest."""

    payload = json.loads(receipt_json)
    if not isinstance(payload, dict) or _canonical_json(payload) != receipt_json:
        raise ValueError("reactive bubble owner receipt is not canonical JSON")
    trial = payload.get("trial")
    branch = payload.get("physical_branch")
    certification = payload.get("certification")
    if (
        not isinstance(trial, dict)
        or not isinstance(branch, dict)
        or not isinstance(certification, dict)
    ):
        raise TypeError("reactive bubble owner receipt is incomplete")
    try:
        temperature_k = float(requested_temperature.to("kelvin").magnitude)
    except (AttributeError, TypeError, ValueError) as error:
        raise TypeError("reactive bubble receipt validation requires a temperature") from error
    expected_trial = _canonical_fingerprint(
        "epcsaft.active-parameter-trial",
        {
            "parameter_fingerprint": descriptor.parameter_fingerprint,
            "identities": active_parameters.identities,
            "values": active_parameters.values,
            "units": tuple(parameter.unit for parameter in active_parameters.coordinates),
            "active_parameter_fingerprint": active_parameters.fingerprint,
        },
    )
    observations = payload.get("observations")
    search = payload.get("search")
    if not isinstance(observations, list) or len(observations) != len(descriptor.rows):
        raise ValueError("reactive bubble owner receipt observation rows are incomplete")
    if not isinstance(search, dict):
        raise TypeError("reactive bubble owner receipt search accounting is incomplete")
    certification_values = tuple(
        certification.get(field)
        for field in (
            "bubble_log_closure_abs",
            "bubble_log_pressure_derivative",
            "liquid_balance_inf",
            "liquid_electroneutrality_inf",
            "liquid_reaction_affinity_inf",
            "liquid_pressure_relative_residual",
            "phase_pressure_residual_inf_pa",
            "phase_chemical_potential_residual_inf",
        )
    )
    if (
        payload.get("schema") != "epcsaft.reactive-bubble-owner-receipt.v1"
        or payload.get("batch_identity") != batch_identity
        or payload.get("descriptor_fingerprint") != descriptor.fingerprint
        or payload.get("problem_fingerprint") != descriptor.problem_fingerprint
        or payload.get("topology_fingerprint") != descriptor.topology_fingerprint
        or payload.get("reference_fingerprint") != descriptor.reference_fingerprint
        or payload.get("temperature_k") != temperature_k
        or certification.get("status") != "certified_local_equilibrium"
        or certification.get("globality_claim") != "not_established"
        or any(
            certification.get(field) != "passed"
            for field in (
                "numerical_status",
                "physical_status",
                "first_order_status",
                "second_order_status",
            )
        )
        or any(
            isinstance(value, bool)
            or not isinstance(value, (int, float))
            or not math.isfinite(float(value))
            for value in certification_values
        )
        or float(certification["bubble_log_closure_abs"]) > 2.0e-10
        or abs(float(certification["bubble_log_pressure_derivative"])) < 1.0e-10
        or search.get("distinct_branch_count") != 1
        or search.get("selected_branch_ordinal") is None
        or search.get("completeness_claim") != "finite_declared_starts_only"
        or trial.get("parameter_fingerprint") != descriptor.parameter_fingerprint
        or tuple(trial.get("identities", ())) != active_parameters.identities
        or tuple(float(value) for value in trial.get("values", ())) != active_parameters.values
        or tuple(trial.get("units", ()))
        != tuple(parameter.unit for parameter in active_parameters.coordinates)
        or trial.get("active_parameter_fingerprint") != active_parameters.fingerprint
        or trial.get("fingerprint") != expected_trial
        or branch.get("problem_fingerprint") != descriptor.problem_fingerprint
        or branch.get("globality_claim") != "not_established"
        or branch.get("continuation_identity") is None
        or tuple(branch.get("partner_component_ids", ())) != descriptor.partner_component_ids
        or branch.get("partner_parameter_fingerprint") != descriptor.partner_parameter_fingerprint
        or branch.get("active_trial_fingerprint") != trial.get("fingerprint")
        or not str(branch.get("liquid_state_binding_identity", "")).startswith("sha256:")
        or not str(branch.get("branch_family_identity", "")).startswith("sha256:")
        or not str(branch.get("state_identity", "")).startswith("sha256:")
        or not str(branch.get("fingerprint", "")).startswith("sha256:")
        or tuple(row.get("identity") for row in observations) != descriptor.row_identities
        or tuple(row.get("unit") for row in observations) != descriptor.row_units
        or any(
            not isinstance(row.get("jacobian"), list)
            or len(row["jacobian"]) != len(active_parameters.identities)
            or not isinstance(row.get("value"), (int, float))
            or isinstance(row.get("value"), bool)
            or not math.isfinite(float(row["value"]))
            or any(
                isinstance(value, bool)
                or not isinstance(value, (int, float))
                or not math.isfinite(float(value))
                for value in row["jacobian"]
            )
            for row in observations
        )
    ):
        raise ValueError("reactive bubble owner receipt identity mismatch")
    branch_payload = {
        key: branch[key]
        for key in (
            "branch_family_identity",
            "state_identity",
            "problem_fingerprint",
            "continuation_identity",
            "liquid_branch_family_identity",
            "liquid_state_binding_identity",
            "partner_parameter_fingerprint",
            "partner_component_ids",
            "active_trial_fingerprint",
            "pressure_pa",
            "liquid_pressure_pa",
            "liquid_mole_fractions",
            "vapor_mole_fractions",
            "liquid_molar_volume_m3_per_mol",
            "vapor_molar_volume_m3_per_mol",
            "globality_claim",
        )
    }
    expected_branch = _canonical_fingerprint(
        "epcsaft.reactive-bubble-physical-branch", branch_payload
    )
    expected_branch_family = _canonical_fingerprint(
        "epcsaft.reactive-bubble-branch-family",
        {
            "problem_fingerprint": descriptor.problem_fingerprint,
            "topology_fingerprint": descriptor.topology_fingerprint,
            "reference_fingerprint": descriptor.reference_fingerprint,
            "continuation_identity": branch["continuation_identity"],
            "liquid_branch_family_identity": branch["liquid_branch_family_identity"],
        },
    )
    expected_binding = _canonical_fingerprint(
        "epcsaft.reactive-bubble-liquid-state-binding",
        {
            "component_ids": descriptor.component_ids,
            "mole_fractions": branch["liquid_mole_fractions"],
            "molar_volume_m3_per_mol": branch["liquid_molar_volume_m3_per_mol"],
            "pressure_pa": branch["liquid_pressure_pa"],
            "liquid_state_identity": payload["liquid_state_identity"],
            "liquid_branch_family_identity": payload["liquid_branch_family_identity"],
            "parameter_fingerprint": trial["parameter_fingerprint"],
        },
    )
    expected_state = _canonical_fingerprint(
        "epcsaft.reactive-bubble-state",
        {
            "pressure_pa": branch["pressure_pa"],
            "liquid_state_identity": payload["liquid_state_identity"],
            "liquid_mole_fractions": branch["liquid_mole_fractions"],
            "vapor_mole_fractions": branch["vapor_mole_fractions"],
            "liquid_molar_volume_m3_per_mol": branch["liquid_molar_volume_m3_per_mol"],
            "vapor_molar_volume_m3_per_mol": branch["vapor_molar_volume_m3_per_mol"],
            "topology_fingerprint": descriptor.topology_fingerprint,
        },
    )
    if (
        branch.get("fingerprint") != expected_branch
        or payload.get("pressure_pa") != branch.get("pressure_pa")
        or branch.get("branch_family_identity") != expected_branch_family
        or branch.get("liquid_state_binding_identity") != expected_binding
        or branch.get("state_identity") != expected_state
    ):
        raise ValueError("reactive bubble owner receipt physical branch is inconsistent")
    return "sha256:" + hashlib.sha256(receipt_json.encode()).hexdigest()


__all__ = [
    "IdealGasVaporComparison",
    "ProviderNonidealVapor",
    "ReactiveBubbleCertification",
    "ReactiveBubbleContinuationReference",
    "ReactiveBubbleObservationDescriptor",
    "ReactiveBubbleObservationResult",
    "ReactiveBubbleObservationRow",
    "ReactiveBubbleObservationRowResult",
    "ReactiveBubblePhysicalBranch",
    "ReactiveBubbleSearch",
    "ReactiveBubbleStartAttempt",
    "ReactiveBubbleVLEProblem",
    "certify_reactive_bubble_continuation_reference",
    "evaluate_reactive_bubble_observations",
    "reactive_bubble_observation_descriptor",
    "serialize_reactive_bubble_owner_receipt",
    "validate_reactive_bubble_owner_receipt",
]
