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
    ChemicalEquilibriumSensitivity,
    EquilibriumPhaseState,
    ProviderPhase,
    _support_active_parameters,
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
    boundary_directional_derivatives: dict[str, float]
    boundary_direct_derivatives: dict[str, float]
    vapor_fraction_directional_derivatives: dict[str, tuple[float, ...]]
    vapor_fraction_direct_derivatives: dict[str, tuple[float, ...]]
    partner_parameter_fingerprint: str
    partner_component_ids: tuple[str, ...]
    chemical_amount_directions: dict[str, tuple[float, ...]]
    chemical_volume_directions: dict[str, float]
    balance_direction_ids: tuple[str, ...]


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
        if row.quantity in {
            "bubble_pressure_pa",
            "liquid_molar_density_mol_per_m3",
        }:
            if row.coefficients:
                raise ValueError("scalar reactive bubble row does not accept coefficients")
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
    )


def _chemical_directions(
    sensitivity: ChemicalEquilibriumSensitivity,
    active_ids: tuple[str, ...],
) -> tuple[
    dict[str, tuple[float, ...]],
    dict[str, float],
    tuple[str, ...],
]:
    if sensitivity.status != "available":
        raise _BubbleFailure(
            "first_order_failure",
            sensitivity.failure_reason or "chemical sensitivity is unavailable",
        )
    requested = ("pressure_pa", *active_ids)
    found_amount_rows: dict[str, tuple[float, ...]] = {}
    found_volume_rows: dict[str, float] = {}
    balance_ids: list[str] = []
    for parameter, amount_row, volume in zip(
        sensitivity.parameters,
        sensitivity.amount_derivatives,
        sensitivity.volume_derivatives,
        strict=True,
    ):
        identity: str | None
        if parameter.kind == "compiled_balance_total":
            if parameter.source_index != len(balance_ids):
                raise _BubbleFailure(
                    "parameter_topology_mismatch",
                    "compiled balance sensitivity order changed",
                )
            identity = f"compiled-balance-total/{parameter.source_index}"
            balance_ids.append(identity)
        else:
            identity = "pressure_pa" if parameter.name == "pressure_pa" else parameter.identity
            if identity is None:
                continue
        if identity in requested or identity in balance_ids:
            found_amount_rows[identity] = amount_row
            found_volume_rows[identity] = volume
    if any(identity not in found_amount_rows for identity in requested):
        raise _BubbleFailure(
            "parameter_topology_mismatch",
            "chemical sensitivity omitted or reordered a coupled input direction",
        )
    ordered = (*requested, *balance_ids)
    return (
        {identity: found_amount_rows[identity] for identity in ordered},
        {identity: found_volume_rows[identity] for identity in ordered},
        tuple(balance_ids),
    )


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
    tuple[str, ...],
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
    amount_rows, volume_rows, balance_ids = _chemical_directions(
        equilibrium.sensitivity, active_ids
    )
    if equilibrium.diagnostics.active_trace_species:
        raise _BubbleFailure(
            "active_bound_ambiguity",
            "the reactive liquid contains a trace-bound species at the coupled state",
        )
    return liquid, equilibrium, amount_rows, volume_rows, balance_ids


@dataclass(frozen=True)
class _IncipientVaporPoint:
    mole_fractions: tuple[float, ...]
    volume_m3: float
    pressure_pa: float
    packing_fraction: float
    chemical_potential_over_rt: tuple[float, ...]
    mu_log_ratio_jacobian: tuple[tuple[float, ...], ...]
    mu_pressure_derivative: tuple[float, ...]
    mu_parameter_derivatives: dict[str, tuple[float, ...]]


def _softmax_reference_chart(log_ratios: tuple[float, ...]) -> tuple[float, ...]:
    logits = (*log_ratios, 0.0)
    maximum = max(logits)
    weights = tuple(math.exp(value - maximum) for value in logits)
    total = math.fsum(weights)
    fractions = tuple(value / total for value in weights)
    if any(not math.isfinite(value) or value <= 0.0 for value in fractions):
        raise _BubbleFailure(
            "active_bound_ambiguity",
            "declared vapor support left its strict positive composition chart",
        )
    return fractions


def _solve_incipient_linear_system(
    matrix: tuple[tuple[float, ...], ...],
    right_hand_side: tuple[float, ...],
) -> tuple[float, ...]:
    """Solve one small dense tangent system with fail-closed pivot checks."""

    dimension = len(right_hand_side)
    if dimension == 0:
        return ()
    work = [list(row) + [value] for row, value in zip(matrix, right_hand_side, strict=True)]
    scale = max(1.0, max(abs(value) for row in matrix for value in row))
    minimum_pivot = math.inf
    maximum_pivot = 0.0
    for column in range(dimension):
        pivot = max(range(column, dimension), key=lambda row: abs(work[row][column]))
        pivot_value = abs(work[pivot][column])
        if not math.isfinite(pivot_value) or pivot_value <= 1.0e-13 * scale:
            raise _BubbleFailure(
                "singular_sensitivity",
                "incipient-vapor composition Jacobian is singular",
            )
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
        minimum_pivot = min(minimum_pivot, pivot_value)
        maximum_pivot = max(maximum_pivot, pivot_value)
        for row in range(column + 1, dimension):
            multiplier = work[row][column] / work[column][column]
            for entry in range(column, dimension + 1):
                work[row][entry] -= multiplier * work[column][entry]
    if maximum_pivot / minimum_pivot > 1.0e10:
        raise _BubbleFailure(
            "ill_conditioned_sensitivity",
            "incipient-vapor composition Jacobian is ill-conditioned",
        )
    solution = [0.0] * dimension
    for row in range(dimension - 1, -1, -1):
        solution[row] = (
            work[row][dimension]
            - math.fsum(
                work[row][column] * solution[column] for column in range(row + 1, dimension)
            )
        ) / work[row][row]
    if any(not math.isfinite(value) for value in solution):
        raise _BubbleFailure(
            "singular_sensitivity",
            "incipient-vapor linear solve returned nonfinite values",
        )
    return tuple(solution)


def _evaluate_incipient_vapor_point(
    *,
    model: epcsaft.Pure | epcsaft.Mixture,
    active_parameters: epcsaft.ActiveParameterSet | None,
    all_active_identities: tuple[str, ...],
    temperature: Any,
    temperature_k: float,
    pressure_pa: float,
    log_ratios: tuple[float, ...],
    packing_fraction_interval: tuple[float, float],
) -> _IncipientVaporPoint:
    """Evaluate one declared vapor root and its exact local tangent blocks."""

    fractions = _softmax_reference_chart(log_ratios)
    active_ids = () if active_parameters is None else active_parameters.identities
    active_position = {identity: index for index, identity in enumerate(active_ids)}
    active_native = (
        model._native
        if active_parameters is None
        else model._native.with_active_coordinates(model._active_coordinates(active_parameters))
    )
    component_count = len(fractions)
    state_count = component_count + 1
    volume_coordinate = component_count
    try:
        state = model._state(
            T=temperature,
            P=pressure_pa * epcsaft.unit_registry.pascal,
            x=fractions,
            phase="vapor",
            native_model=active_native,
        )
        diagnostics = state.density_diagnostics
        if diagnostics is None or not diagnostics.stable or diagnostics.branch != "vapor":
            raise _BubbleFailure(
                "no_incipient_phase",
                "declared vapor Provider did not return a strict-stable vapor root",
            )
        volume = 1.0 / float(state.molar_density.to("mole / meter ** 3").magnitude)
        block = _equilibrium.evaluate_parameterized_phase(
            _native_transport(model, active_parameters),
            temperature_k,
            fractions,
            volume,
            *packing_fraction_interval,
        )
        for _ in range(4):
            pressure_residual = float(cast(Any, block["pressure_pa"])) - pressure_pa
            if abs(pressure_residual) <= 1.0e-6:
                break
            curvature = tuple(
                float(cast(Any, value))
                for value in cast(Sequence[object], block["state_hessian"])
            )[volume_coordinate * state_count + volume_coordinate]
            if not math.isfinite(curvature) or curvature <= 0.0:
                raise _BubbleFailure(
                    "second_order_failure",
                    "declared vapor mechanical polish lacks positive curvature",
                )
            candidate = volume + pressure_residual / (
                8.31446261815324 * temperature_k * curvature
            )
            if not math.isfinite(candidate) or candidate <= 0.0:
                raise _BubbleFailure(
                    "numerical_convergence_failure",
                    "declared vapor mechanical polish left the positive-volume domain",
                )
            volume = candidate
            block = _equilibrium.evaluate_parameterized_phase(
                _native_transport(model, active_parameters),
                temperature_k,
                fractions,
                volume,
                *packing_fraction_interval,
            )
        if abs(float(cast(Any, block["pressure_pa"])) - pressure_pa) > 1.0e-6:
            raise _BubbleFailure(
                "numerical_convergence_failure",
                "declared vapor mechanical root did not reach the requested pressure",
            )
    except _BubbleFailure:
        raise
    except (RuntimeError, TypeError, ValueError) as error:
        raise _BubbleFailure("provider_domain_rejection", str(error)) from error

    gradient = tuple(
        float(cast(Any, value)) for value in cast(Sequence[object], block["state_gradient"])
    )
    hessian = tuple(
        float(cast(Any, value)) for value in cast(Sequence[object], block["state_hessian"])
    )
    direct_mu = tuple(
        float(cast(Any, value))
        for value in cast(
            Sequence[object], block["chemical_potential_parameter_derivatives_over_rt"]
        )
    )
    pressure_parameter = tuple(
        float(cast(Any, value))
        for value in cast(Sequence[object], block["pressure_parameter_derivatives_pa"])
    )
    if (
        len(gradient) != state_count
        or len(hessian) != state_count * state_count
        or len(direct_mu) != component_count * len(active_ids)
        or len(pressure_parameter) != len(active_ids)
    ):
        raise _BubbleFailure("owner_contract_failure", "vapor Provider derivative layout changed")

    h_vv = hessian[volume_coordinate * state_count + volume_coordinate]
    h_scale = max(
        1.0,
        max(abs(hessian[volume_coordinate * state_count + index]) for index in range(state_count)),
    )
    if not math.isfinite(h_vv) or h_vv <= 1.0e-14 * h_scale:
        raise _BubbleFailure(
            "second_order_failure",
            "declared vapor root lacks strict mechanical stability",
        )

    reference = component_count - 1
    chart = tuple(
        tuple(
            fractions[component]
            * ((1.0 if component == coordinate else 0.0) - fractions[coordinate])
            for coordinate in range(reference)
        )
        for component in range(component_count)
    )
    volume_log_ratio = tuple(
        -math.fsum(
            hessian[volume_coordinate * state_count + component] * chart[component][coordinate]
            for component in range(component_count)
        )
        / h_vv
        for coordinate in range(reference)
    )
    mu_log_ratio = tuple(
        tuple(
            math.fsum(
                hessian[component * state_count + species] * chart[species][coordinate]
                for species in range(component_count)
            )
            + hessian[component * state_count + volume_coordinate] * volume_log_ratio[coordinate]
            for coordinate in range(reference)
        )
        for component in range(component_count)
    )
    rt = 8.31446261815324 * temperature_k
    volume_pressure = -1.0 / (rt * h_vv)
    mu_pressure = tuple(
        hessian[component * state_count + volume_coordinate] * volume_pressure
        for component in range(component_count)
    )
    mu_parameters: dict[str, tuple[float, ...]] = {}
    for identity in all_active_identities:
        parameter = active_position.get(identity)
        if parameter is None:
            mu_parameters[identity] = (0.0,) * component_count
            continue
        volume_parameter = pressure_parameter[parameter] / (rt * h_vv)
        mu_parameters[identity] = tuple(
            direct_mu[component * len(active_ids) + parameter]
            + hessian[component * state_count + volume_coordinate] * volume_parameter
            for component in range(component_count)
        )
    return _IncipientVaporPoint(
        mole_fractions=fractions,
        volume_m3=volume,
        pressure_pa=float(cast(Any, block["pressure_pa"])),
        packing_fraction=float(cast(Any, block["packing_fraction"])),
        chemical_potential_over_rt=gradient[:component_count],
        mu_log_ratio_jacobian=mu_log_ratio,
        mu_pressure_derivative=mu_pressure,
        mu_parameter_derivatives=mu_parameters,
    )


def _solve_declared_incipient_vapor(
    *,
    model: epcsaft.Pure | epcsaft.Mixture,
    active_parameters: epcsaft.ActiveParameterSet | None,
    all_active_identities: tuple[str, ...],
    temperature: Any,
    temperature_k: float,
    pressure_pa: float,
    liquid_potentials: tuple[float, ...],
    continuation_fractions: tuple[float, ...] | None,
    packing_fraction_interval: tuple[float, float],
) -> _IncipientVaporPoint:
    """Solve composition inside the declared vapor support.

    The vapor role is already fixed by the problem topology. An unanchored
    boundary therefore retains the lowest-density stationary composition
    reached inside that declared support. A supplied continuation composition
    instead owns local branch selection. Neither case discovers a phase or
    changes the configured phase topology.
    """

    component_count = len(liquid_potentials)
    reference = component_count - 1
    shifted = tuple(value - max(liquid_potentials) for value in liquid_potentials)
    weights = tuple(math.exp(max(-700.0, value)) for value in shifted)
    weight_total = math.fsum(weights)
    starts = [tuple(value / weight_total for value in weights)]
    starts.append(
        (1.0 / component_count,) * component_count
        if continuation_fractions is None
        else continuation_fractions
    )
    # Fugacity equality also admits a dense same-phase stationary solution
    # when the vapor support contains the liquid's neutral components. Generic
    # component-dominant starts expose dilute vapor basins without assuming a
    # volatile species or changing the declared topology.
    minor_fraction = 1.0e-4 / max(1, component_count - 1)
    for dominant in range(component_count):
        starts.append(
            tuple(
                1.0 - minor_fraction * (component_count - 1)
                if component == dominant
                else minor_fraction
                for component in range(component_count)
            )
        )
    cache: dict[tuple[float, ...], _IncipientVaporPoint] = {}

    def evaluate(log_ratios: tuple[float, ...]) -> _IncipientVaporPoint:
        if log_ratios not in cache:
            cache[log_ratios] = _evaluate_incipient_vapor_point(
                model=model,
                active_parameters=active_parameters,
                all_active_identities=all_active_identities,
                temperature=temperature,
                temperature_k=temperature_k,
                pressure_pa=pressure_pa,
                log_ratios=log_ratios,
                packing_fraction_interval=packing_fraction_interval,
            )
        return cache[log_ratios]

    accepted: list[_IncipientVaporPoint] = []
    for fractions in starts:
        log_ratios = tuple(
            math.log(fractions[index] / fractions[reference]) for index in range(reference)
        )
        point: _IncipientVaporPoint | None = None
        for _iteration in range(20):
            try:
                point = evaluate(log_ratios)
            except _BubbleFailure:
                point = None
                break
            residual = tuple(
                (
                    point.chemical_potential_over_rt[index]
                    - point.chemical_potential_over_rt[reference]
                )
                - (liquid_potentials[index] - liquid_potentials[reference])
                for index in range(reference)
            )
            residual_norm = max((abs(value) for value in residual), default=0.0)
            if residual_norm <= 2.0e-10:
                break
            jacobian = tuple(
                tuple(
                    point.mu_log_ratio_jacobian[row][column]
                    - point.mu_log_ratio_jacobian[reference][column]
                    for column in range(reference)
                )
                for row in range(reference)
            )
            step = _solve_incipient_linear_system(jacobian, tuple(-value for value in residual))
            maximum_step = max((abs(value) for value in step), default=0.0)
            step_scale = 1.0 if maximum_step <= 2.0 else 2.0 / maximum_step
            accepted_step = False
            for reduction in range(8):
                factor = step_scale * 0.5**reduction
                candidate = tuple(
                    value + factor * delta for value, delta in zip(log_ratios, step, strict=True)
                )
                try:
                    trial = evaluate(candidate)
                except _BubbleFailure:
                    continue
                trial_residual = tuple(
                    (
                        trial.chemical_potential_over_rt[index]
                        - trial.chemical_potential_over_rt[reference]
                    )
                    - (liquid_potentials[index] - liquid_potentials[reference])
                    for index in range(reference)
                )
                if max((abs(value) for value in trial_residual), default=0.0) < residual_norm:
                    log_ratios = candidate
                    accepted_step = True
                    break
            if not accepted_step:
                point = None
                break
        if point is None:
            continue
        final_residual = tuple(
            (point.chemical_potential_over_rt[index] - point.chemical_potential_over_rt[reference])
            - (liquid_potentials[index] - liquid_potentials[reference])
            for index in range(reference)
        )
        if max((abs(value) for value in final_residual), default=0.0) > 2.0e-10:
            continue
        if not any(
            max(
                abs(math.log(left / right))
                for left, right in zip(point.mole_fractions, existing.mole_fractions, strict=True)
            )
            <= 2.0e-8
            for existing in accepted
        ):
            accepted.append(point)
    if not accepted:
        raise _BubbleFailure(
            "no_incipient_phase",
            "no declared strict-stable vapor composition satisfies fugacity equality",
        )
    if continuation_fractions is not None:
        ranked_by_continuation = sorted(
            accepted,
            key=lambda point: (
                max(
                    abs(math.log(value / anchor))
                    for value, anchor in zip(
                        point.mole_fractions, continuation_fractions, strict=True
                    )
                ),
                point.mole_fractions,
            ),
        )
        if len(ranked_by_continuation) > 1:
            best_distance = max(
                abs(math.log(value / anchor))
                for value, anchor in zip(
                    ranked_by_continuation[0].mole_fractions,
                    continuation_fractions,
                    strict=True,
                )
            )
            next_distance = max(
                abs(math.log(value / anchor))
                for value, anchor in zip(
                    ranked_by_continuation[1].mole_fractions,
                    continuation_fractions,
                    strict=True,
                )
            )
            if abs(next_distance - best_distance) <= 2.0e-8:
                raise _BubbleFailure(
                    "multiple_unresolved_branches",
                    "multiple vapor compositions are equally close to the continuation anchor",
                )
        return ranked_by_continuation[0]
    # The phase role is an input, so branch selection follows that role rather
    # than a tangent-plane-distance ranking. At fixed T and P the largest
    # molar volume is the vapor branch among the stationary compositions found
    # inside the declared support.
    ranked = sorted(
        accepted,
        key=lambda point: (
            -point.volume_m3,
            point.mole_fractions,
        ),
    )
    best = ranked[0]
    if len(ranked) > 1:
        next_volume = ranked[1].volume_m3
    else:
        next_volume = 0.0
    if next_volume > 0.0 and abs(math.log(best.volume_m3 / next_volume)) <= 2.0e-10:
        raise _BubbleFailure(
            "multiple_unresolved_branches",
            "multiple incipient-vapor compositions have indistinguishable vapor volumes",
        )
    return best


def _nonideal_closure(
    problem: ReactiveBubbleVLEProblem,
    temperature: Any,
    descriptor: ReactiveBubbleObservationDescriptor,
    active_parameters: epcsaft.ActiveParameterSet | None,
    pressure_pa: float,
) -> _Closure:
    """Close the declared vapor against the exact certified reactive liquid.

    Phase count and roles are immutable inputs. Chemical equilibrium owns the
    liquid ``(T,n,V)`` state; this operation solves only the declared neutral
    vapor composition and its strict-stable Provider density branch.
    """

    liquid, equilibrium, amount_rows, volume_rows, balance_ids = _liquid_trial_state(
        problem, temperature, descriptor, active_parameters, pressure_pa
    )
    active_ids = descriptor.active_parameter_identities
    phase = cast(ProviderPhase, problem.liquid_problem.phase)
    partner_parameters = phase.model._parameters.select(problem.vapor_component_ids)
    partner_model: epcsaft.Pure | epcsaft.Mixture = (
        epcsaft.Pure(partner_parameters)
        if len(problem.vapor_component_ids) == 1
        else epcsaft.Mixture(partner_parameters)
    )
    partner_active = _support_active_parameters(partner_parameters, active_parameters)
    if (
        partner_model.parameter_fingerprint != descriptor.partner_parameter_fingerprint
        or partner_model.component_ids != descriptor.partner_component_ids
    ):
        raise _BubbleFailure(
            "parameter_topology_mismatch",
            "declared vapor Provider does not match the ordered neutral support",
        )
    try:
        liquid_block = _equilibrium.evaluate_parameterized_phase(
            _native_transport(phase.model, active_parameters),
            equilibrium.temperature_k,
            equilibrium.amounts_mol,
            equilibrium.volume_m3,
            *phase.admissible_packing_fraction_interval,
        )
    except (RuntimeError, TypeError, ValueError) as error:
        raise _BubbleFailure("provider_domain_rejection", str(error)) from error

    component_count = len(descriptor.component_ids)
    state_count = component_count + 1
    liquid_gradient = tuple(
        float(cast(Any, value)) for value in cast(Sequence[object], liquid_block["state_gradient"])
    )
    liquid_hessian = tuple(
        float(cast(Any, value)) for value in cast(Sequence[object], liquid_block["state_hessian"])
    )
    liquid_direct_mu = tuple(
        float(cast(Any, value))
        for value in cast(
            Sequence[object], liquid_block["chemical_potential_parameter_derivatives_over_rt"]
        )
    )
    if (
        len(liquid_gradient) != state_count
        or len(liquid_hessian) != state_count * state_count
        or len(liquid_direct_mu) != component_count * len(active_ids)
    ):
        raise _BubbleFailure("owner_contract_failure", "reactive liquid derivative layout changed")

    component_position = {
        identity: index for index, identity in enumerate(descriptor.component_ids)
    }
    vapor_positions = tuple(
        component_position[identity] for identity in problem.vapor_component_ids
    )
    liquid_potentials = tuple(liquid_gradient[index] for index in vapor_positions)
    direction_ids = tuple(amount_rows)
    liquid_mu_directions: dict[str, tuple[float, ...]] = {}
    for identity in direction_ids:
        parameter = active_ids.index(identity) if identity in active_ids else None
        liquid_mu_directions[identity] = tuple(
            math.fsum(
                liquid_hessian[component * state_count + species] * derivative
                for species, derivative in enumerate(amount_rows[identity])
            )
            + liquid_hessian[component * state_count + component_count] * volume_rows[identity]
            + (
                0.0
                if parameter is None
                else liquid_direct_mu[component * len(active_ids) + parameter]
            )
            for component in vapor_positions
        )

    vapor = _solve_declared_incipient_vapor(
        model=partner_model,
        active_parameters=partner_active,
        all_active_identities=active_ids,
        temperature=temperature,
        temperature_k=equilibrium.temperature_k,
        pressure_pa=pressure_pa,
        liquid_potentials=liquid_potentials,
        continuation_fractions=(
            None
            if problem.continuation_reference is None
            else problem.continuation_reference.vapor_mole_fractions
        ),
        packing_fraction_interval=phase.admissible_packing_fraction_interval,
    )

    vapor_count = len(problem.vapor_component_ids)
    reference = vapor_count - 1
    g_jacobian = tuple(
        tuple(
            vapor.mu_log_ratio_jacobian[row][column]
            - vapor.mu_log_ratio_jacobian[reference][column]
            for column in range(reference)
        )
        for row in range(reference)
    )
    h_log_ratio = vapor.mu_log_ratio_jacobian[reference]
    chart = tuple(
        tuple(
            vapor.mole_fractions[component]
            * ((1.0 if component == coordinate else 0.0) - vapor.mole_fractions[coordinate])
            for coordinate in range(reference)
        )
        for component in range(vapor_count)
    )
    closure_directions: dict[str, float] = {}
    vapor_fraction_directions: dict[str, tuple[float, ...]] = {}
    for identity in direction_ids:
        vapor_base = (
            vapor.mu_pressure_derivative
            if identity == "pressure_pa"
            else vapor.mu_parameter_derivatives[identity]
            if identity in active_ids
            else (0.0,) * vapor_count
        )
        potential_direction = tuple(
            vapor_value - liquid_value
            for vapor_value, liquid_value in zip(
                vapor_base, liquid_mu_directions[identity], strict=True
            )
        )
        g_direction = tuple(
            potential_direction[index] - potential_direction[reference]
            for index in range(reference)
        )
        log_ratio_direction = _solve_incipient_linear_system(
            g_jacobian, tuple(-value for value in g_direction)
        )
        closure_directions[identity] = potential_direction[reference] + math.fsum(
            derivative * direction
            for derivative, direction in zip(h_log_ratio, log_ratio_direction, strict=True)
        )
        vapor_fraction_directions[identity] = tuple(
            math.fsum(
                chart[component][coordinate] * log_ratio_direction[coordinate]
                for coordinate in range(reference)
            )
            for component in range(vapor_count)
        )

    residual = vapor.chemical_potential_over_rt[reference] - liquid_potentials[reference]
    residual_log_pressure_derivative = pressure_pa * closure_directions["pressure_pa"]
    if not math.isfinite(residual_log_pressure_derivative):
        raise _BubbleFailure("first_order_failure", "bubble closure derivative is non-finite")

    # _row_values consumes the generic B(P,x(theta,P))-P Schur layout.  For
    # the reduced chemical-potential closure h=0, encode B_P-1=P*h_P and
    # B_theta=P*h_theta; the resulting total dP/dtheta is identical.
    boundary_directional = {
        "pressure_pa": 1.0 + residual_log_pressure_derivative,
        **{
            identity: pressure_pa * closure_directions[identity]
            for identity in direction_ids
            if identity != "pressure_pa"
        },
    }
    boundary_direct = {identity: 0.0 for identity in active_ids}
    vapor_direct = {identity: (0.0,) * vapor_count for identity in active_ids}

    total = math.fsum(equilibrium.amounts_mol)
    liquid_molar_volume = equilibrium.volume_m3 / total
    liquid_phase = EquilibriumPhaseState(
        role="liquid",
        support="all_components",
        supported_component_ids=descriptor.component_ids,
        mole_fractions=equilibrium.mole_fractions,
        molar_volume_m3_per_mol=liquid_molar_volume,
        molar_density_mol_m3=1.0 / liquid_molar_volume,
        packing_fraction=float(cast(Any, liquid_block["packing_fraction"])),
        pressure_pa=float(cast(Any, liquid_block["pressure_pa"])),
        chemical_potential_component_ids=descriptor.component_ids,
        chemical_potential_over_rt=liquid_gradient[:component_count],
        mechanical_class=liquid.state.phase_role_receipt.mechanical_class,
    )
    vapor_phase = EquilibriumPhaseState(
        role="vapor",
        support="declared_neutral_components",
        supported_component_ids=descriptor.vapor_component_ids,
        mole_fractions=vapor.mole_fractions,
        molar_volume_m3_per_mol=vapor.volume_m3,
        molar_density_mol_m3=1.0 / vapor.volume_m3,
        packing_fraction=vapor.packing_fraction,
        pressure_pa=vapor.pressure_pa,
        chemical_potential_component_ids=descriptor.vapor_component_ids,
        chemical_potential_over_rt=vapor.chemical_potential_over_rt,
        mechanical_class="strict_stable_vapor",
    )
    return _Closure(
        pressure_pa,
        residual,
        residual_log_pressure_derivative,
        liquid,
        (liquid_phase, vapor_phase),
        boundary_directional,
        boundary_direct,
        vapor_fraction_directions,
        vapor_direct,
        partner_model.parameter_fingerprint,
        partner_model.component_ids,
        amount_rows,
        volume_rows,
        balance_ids,
    )


def _ideal_closure(
    problem: ReactiveBubbleVLEProblem,
    temperature: Any,
    descriptor: ReactiveBubbleObservationDescriptor,
    active_parameters: epcsaft.ActiveParameterSet | None,
    pressure_pa: float,
) -> _Closure:
    liquid, equilibrium, amount_rows, volume_rows, balance_ids = _liquid_trial_state(
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
    direction_ids = tuple(amount_rows)
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
        boundary_directional,
        boundary_direct,
        vapor_directional,
        vapor_direct,
        descriptor.partner_parameter_fingerprint,
        descriptor.partner_component_ids,
        amount_rows,
        volume_rows,
        balance_ids,
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


def _fixed_pressure_coordinate_derivatives(
    closure: _Closure,
    balance_direction: tuple[float, ...],
) -> tuple[float, dict[str, float]]:
    balance_ids = closure.balance_direction_ids
    if len(balance_ids) != len(balance_direction) or not balance_ids:
        raise _BubbleFailure(
            "parameter_topology_mismatch",
            "fixed-pressure incipient closure lacks the conserved-balance directions",
        )
    vapor = closure.phases[1]
    residual_scale = (
        vapor.pressure_pa
        if vapor.mechanical_class == "strict_stable_ideal_gas"
        else closure.pressure_pa
    )
    if not math.isfinite(residual_scale) or residual_scale <= 0.0:
        raise _BubbleFailure("first_order_failure", "incipient residual scale is invalid")

    def residual_direction(identity: str) -> float:
        return (
            closure.boundary_directional_derivatives[identity]
            + closure.boundary_direct_derivatives.get(identity, 0.0)
        ) / residual_scale

    coordinate_slope = math.fsum(
        residual_direction(identity) * direction
        for identity, direction in zip(balance_ids, balance_direction, strict=True)
    )
    slope_scale = max(
        1.0,
        *(abs(residual_direction(identity)) for identity in balance_ids),
    )
    if abs(coordinate_slope) <= 1.0e-12 * slope_scale:
        raise _BubbleFailure(
            "singular_sensitivity", "fixed-pressure affine-feed closure is singular"
        )
    if abs(coordinate_slope) <= 1.0e-8 * slope_scale:
        raise _BubbleFailure(
            "ill_conditioned_sensitivity",
            "fixed-pressure affine-feed closure is ill-conditioned",
        )
    return coordinate_slope, {
        identity: residual_direction(identity)
        for identity in closure.boundary_directional_derivatives
    }


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
    """Verify the phase receipt carries the exact chemical ``(T,n,V)`` liquid."""

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
    ):
        raise _BubbleFailure(
            "branch_loss",
            "coupled result does not bind to the certified chemical liquid branch",
        )


def _solve_starts(
    problem: ReactiveBubbleVLEProblem,
    temperature: Any,
    descriptor: ReactiveBubbleObservationDescriptor,
    active_parameters: epcsaft.ActiveParameterSet | None,
) -> tuple[tuple[_Closure, ...], ReactiveBubbleSearch]:
    lower, upper = problem.pressure_interval_pa
    starts = tuple(
        (
            _canonical_fingerprint(
                "epcsaft.reactive-bubble-pressure-start",
                {
                    "problem_identity": problem.identity,
                    "kind": "declared",
                    "ordinal": ordinal,
                    "value": value,
                },
            ),
            value,
        )
        for ordinal, value in enumerate(problem.pressure_starts_pa)
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
            duplicate = next(
                (
                    index
                    for index, existing in enumerate(accepted)
                    if _candidate_distance(final, existing) <= 2.0e-8
                ),
                None,
            )
            if duplicate is None:
                accepted.append(final)
            else:
                existing = accepted[duplicate]
                final_key = (
                    final.pressure_pa,
                    final.phases[1].mole_fractions,
                    final.phases[0].molar_volume_m3_per_mol,
                    final.phases[1].molar_volume_m3_per_mol,
                )
                existing_key = (
                    existing.pressure_pa,
                    existing.phases[1].mole_fractions,
                    existing.phases[0].molar_volume_m3_per_mol,
                    existing.phases[1].molar_volume_m3_per_mol,
                )
                if final_key < existing_key:
                    accepted[duplicate] = final
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
            candidate.liquid.state.phase_role_receipt.provider_evaluation_count + 1
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
    derivative_ids: tuple[str, ...] | None = None,
    fixed_balance_direction: tuple[float, ...] | None = None,
) -> tuple[ReactiveBubbleObservationRowResult, ...]:
    direction_ids = (
        descriptor.active_parameter_identities if derivative_ids is None else derivative_ids
    )
    equilibrium = closure.liquid.equilibrium
    assert equilibrium is not None
    amounts = equilibrium.amounts_mol
    volume = equilibrium.volume_m3
    total = math.fsum(amounts)
    liquid_amount_total: dict[str, tuple[float, ...]] = {}
    liquid_volume_total: dict[str, float] = {}
    vapor = closure.phases[1]
    vapor_derivatives: dict[str, tuple[float, ...]] = {}
    if fixed_balance_direction is None:
        pressure_partial = closure.boundary_directional_derivatives["pressure_pa"]
        schur_pressure = pressure_partial - 1.0
        schur_scale = max(1.0, abs(pressure_partial))
        if not math.isfinite(schur_pressure) or abs(schur_pressure) < 1.0e-12 * schur_scale:
            raise _BubbleFailure(
                "singular_sensitivity", "coupled bubble Schur complement is singular"
            )
        if abs(schur_pressure) < 1.0e-8 * schur_scale:
            raise _BubbleFailure(
                "ill_conditioned_sensitivity",
                "coupled bubble Schur complement is ill-conditioned",
            )
        pressure_derivatives = {
            identity: -(
                closure.boundary_direct_derivatives.get(identity, 0.0)
                + closure.boundary_directional_derivatives[identity]
            )
            / schur_pressure
            for identity in direction_ids
        }
        for identity in direction_ids:
            pressure_derivative = pressure_derivatives[identity]
            liquid_amount_total[identity] = tuple(
                direct + pressure * pressure_derivative
                for direct, pressure in zip(
                    closure.chemical_amount_directions[identity],
                    closure.chemical_amount_directions["pressure_pa"],
                    strict=True,
                )
            )
            liquid_volume_total[identity] = (
                closure.chemical_volume_directions[identity]
                + closure.chemical_volume_directions["pressure_pa"] * pressure_derivative
            )
            vapor_derivatives[identity] = tuple(
                direct + directional + pressure * pressure_derivative
                for direct, directional, pressure in zip(
                    closure.vapor_fraction_direct_derivatives.get(
                        identity, (0.0,) * len(vapor.mole_fractions)
                    ),
                    closure.vapor_fraction_directional_derivatives[identity],
                    closure.vapor_fraction_directional_derivatives["pressure_pa"],
                    strict=True,
                )
            )
    else:
        if derivative_ids is not None:
            raise _BubbleFailure(
                "parameter_topology_mismatch",
                "fixed-pressure rows require active parameter directions",
            )
        slope, residual_directions = _fixed_pressure_coordinate_derivatives(
            closure, fixed_balance_direction
        )
        coordinate_derivatives = {
            identity: -residual_directions[identity] / slope for identity in direction_ids
        }
        pressure_derivatives = {identity: 0.0 for identity in direction_ids}
        balance_ids = closure.balance_direction_ids
        coordinate_amount = tuple(
            math.fsum(
                closure.chemical_amount_directions[identity][component] * direction
                for identity, direction in zip(
                    balance_ids, fixed_balance_direction, strict=True
                )
            )
            for component in range(len(descriptor.component_ids))
        )
        coordinate_volume = math.fsum(
            closure.chemical_volume_directions[identity] * direction
            for identity, direction in zip(balance_ids, fixed_balance_direction, strict=True)
        )
        coordinate_vapor = tuple(
            math.fsum(
                closure.vapor_fraction_directional_derivatives[identity][component] * direction
                for identity, direction in zip(
                    balance_ids, fixed_balance_direction, strict=True
                )
            )
            for component in range(len(vapor.mole_fractions))
        )
        for identity in direction_ids:
            dq = coordinate_derivatives[identity]
            liquid_amount_total[identity] = tuple(
                direct + coordinate * dq
                for direct, coordinate in zip(
                    closure.chemical_amount_directions[identity],
                    coordinate_amount,
                    strict=True,
                )
            )
            liquid_volume_total[identity] = (
                closure.chemical_volume_directions[identity] + coordinate_volume * dq
            )
            vapor_derivatives[identity] = tuple(
                direct + directional + coordinate * dq
                for direct, directional, coordinate in zip(
                    closure.vapor_fraction_direct_derivatives.get(
                        identity, (0.0,) * len(vapor.mole_fractions)
                    ),
                    closure.vapor_fraction_directional_derivatives[identity],
                    coordinate_vapor,
                    strict=True,
                )
            )
    liquid_fraction_total: dict[str, tuple[float, ...]] = {}
    for identity, amount_row in liquid_amount_total.items():
        dtotal = math.fsum(amount_row)
        liquid_fraction_total[identity] = tuple(
            (dn * total - amount * dtotal) / (total * total)
            for amount, dn in zip(amounts, amount_row, strict=True)
        )
    results: list[ReactiveBubbleObservationRowResult] = []
    for row in descriptor.rows:
        if row.quantity == "bubble_pressure_pa":
            value = closure.pressure_pa
            jacobian = tuple(pressure_derivatives[identity] for identity in direction_ids)
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
                for identity in direction_ids
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
                for identity in direction_ids
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
                for identity in direction_ids
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
                for identity in direction_ids
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
            for identity in direction_ids:
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
                for identity in direction_ids
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
    phase_pressure_residual = abs(liquid_phase.pressure_pa - vapor_phase.pressure_pa)
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


def _fixed_pressure_coordinate_result(
    problem: ReactiveBubbleVLEProblem,
    descriptor: ReactiveBubbleObservationDescriptor,
    closure: _Closure,
    balance_direction: tuple[float, ...],
) -> ReactiveBubbleObservationResult:
    """Certify one fixed-pressure state owned by the reactive-bubble formulation."""

    _validate_same_liquid_state(closure)
    rows = _row_values(
        descriptor, closure, fixed_balance_direction=balance_direction
    )
    branch = _physical_branch(problem, descriptor, closure)
    certification = _certification(closure)
    liquid_state = closure.liquid.state
    assert liquid_state is not None
    search = ReactiveBubbleSearch(
        attempts=(
            ReactiveBubbleStartAttempt(
                0,
                _canonical_fingerprint(
                    "epcsaft.reactive-bubble-fixed-pressure",
                    {
                        "problem_identity": problem.identity,
                        "pressure_pa": closure.pressure_pa,
                    },
                ),
                closure.pressure_pa,
                "accepted",
                1,
                closure.pressure_pa,
                closure.residual,
                None,
                "",
            ),
        ),
        closure_evaluation_count=1,
        provider_evaluation_count=(
            liquid_state.phase_role_receipt.provider_evaluation_count + 1
        ),
        accepted_candidate_count=1,
        distinct_branch_count=1,
        selected_branch_ordinal=0,
        completeness_claim="finite_declared_starts_only",
    )
    return ReactiveBubbleObservationResult(
        descriptor=descriptor,
        status="evaluated",
        failure=None,
        rows=rows,
        trial=closure.liquid.trial,
        liquid_state=liquid_state,
        liquid_branch=closure.liquid.physical_branch,
        phases=closure.phases,
        pressure_pa=closure.pressure_pa,
        physical_branch=branch,
        certification=certification,
        search=search,
        chemical_diagnostics=closure.liquid.diagnostics,
    )


def _evaluate_reactive_bubble_observations(
    problem: ReactiveBubbleVLEProblem,
    temperature: Any,
    descriptor: ReactiveBubbleObservationDescriptor,
    *,
    active_parameters: epcsaft.ActiveParameterSet | None = None,
) -> tuple[ReactiveBubbleObservationResult, tuple[tuple[float, ...], ...]]:
    """Solve and evaluate one pure, locally certified coupled bubble batch."""

    try:
        expected = reactive_bubble_observation_descriptor(
            problem, descriptor.rows, active_parameters=active_parameters
        )
    except (TypeError, ValueError) as error:
        return _failed_result(descriptor, "parameter_topology_mismatch", str(error)), ()
    if expected.fingerprint != descriptor.fingerprint:
        return (
            _failed_result(
                descriptor, "parameter_topology_mismatch", "reactive bubble descriptor changed"
            ),
            (),
        )
    candidates, search = _solve_starts(problem, temperature, descriptor, active_parameters)
    try:
        selected_ordinal, closure = _select_branch(problem, descriptor, candidates, search)
        search = replace(search, selected_branch_ordinal=selected_ordinal)
        _validate_same_liquid_state(closure)
        rows = _row_values(descriptor, closure)
        balance_rows = _row_values(descriptor, closure, closure.balance_direction_ids)
        branch = _physical_branch(problem, descriptor, closure)
        certification = _certification(closure)
    except _BubbleFailure as error:
        return _failed_result(descriptor, error.code, str(error), search=search), ()
    return (
        ReactiveBubbleObservationResult(
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
        ),
        tuple(cast(tuple[float, ...], row.jacobian) for row in balance_rows),
    )


def evaluate_reactive_bubble_observations(
    problem: ReactiveBubbleVLEProblem,
    temperature: Any,
    descriptor: ReactiveBubbleObservationDescriptor,
    *,
    active_parameters: epcsaft.ActiveParameterSet | None = None,
) -> ReactiveBubbleObservationResult:
    """Solve and evaluate one pure, locally certified coupled bubble batch."""

    result, _ = _evaluate_reactive_bubble_observations(
        problem,
        temperature,
        descriptor,
        active_parameters=active_parameters,
    )
    return result


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
