from __future__ import annotations

import hashlib
import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, replace
from functools import cache
from importlib import metadata
from typing import Any, Literal, cast

from pint import Quantity
from pint.errors import DimensionalityError

import epcsaft
from epcsaft.eos import _native_transport

from . import _equilibrium


@dataclass(frozen=True)
class _HeldPhaseState:
    """One explicit phase state returned by the retained HELD route."""

    amount_mol: float
    mole_fractions: tuple[float, ...]
    volume_m3: float
    molar_density_mol_m3: float
    pressure_pa: float
    chemical_potential_over_rt: tuple[float, ...]


@dataclass(frozen=True)
class HeldDiagnostics:
    """Compact evidence from the bounded HELD search and local refinement."""

    outcome: str
    search_status: str
    solver_status: str
    numerical_status: str
    physical_status: str
    attempts: int
    major_iterations: int
    material_balance_max_abs: float | None
    pressure_stationarity_max_relative: float | None
    kkt_stationarity_max_abs: float | None
    chemical_potential_max_relative: float | None
    search_profiles: tuple[str, ...]
    globality_certificate: str
    provider_evaluation_count: int
    failure_reason: str


@dataclass(frozen=True)
class Held2Diagnostics:
    """Compact evidence from the retained HELD2 controller."""

    outcome: str
    failure_stage: str | None
    solver_status: str
    numerical_status: str
    physical_status: str
    optimizer_solves: int
    optimizer_iterations: int
    major_iterations: int
    material_balance_max_abs: float | None
    electroneutrality_inf: float | None
    pressure_stationarity_max_relative: float | None
    kkt_stationarity_max_abs: float | None
    globality_certificate: str
    provider_evaluation_count: int
    failure_reason: str


PhaseRole = Literal["liquid", "vapor"]


@dataclass(frozen=True, slots=True)
class AllComponents:
    """Every model component is admitted as a positive phase coordinate."""


@dataclass(frozen=True, slots=True)
class NeutralComponents:
    """Only the model's neutral components are present in this phase."""


PhaseSupport = AllComponents | NeutralComponents


@dataclass(frozen=True, slots=True)
class PhaseSpecification:
    """One ordered phase identity, role, and structural component support."""

    role: PhaseRole
    support: PhaseSupport


@dataclass(frozen=True, slots=True)
class DeclaredTopology:
    """Ordered phase specifications; tuple length is the phase count."""

    phases: tuple[PhaseSpecification, ...]

    def __post_init__(self) -> None:
        object.__setattr__(self, "phases", tuple(self.phases))


@dataclass(frozen=True, slots=True)
class FixedTPFeed:
    """Fixed-T,P feed; omitted topology selects HELD/HELD2 discovery."""

    temperature: Quantity[Any]
    pressure: Quantity[Any]
    feed_amounts_mol: tuple[float, ...]
    topology: DeclaredTopology | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "feed_amounts_mol", tuple(self.feed_amounts_mol))


@dataclass(frozen=True, slots=True)
class FixedTCoexistence:
    """Fixed-temperature two-phase boundary with one fixed reference composition."""

    temperature: Quantity[Any]
    topology: DeclaredTopology
    reference_phase_index: int
    reference_mole_fractions: tuple[float, ...]

    def __post_init__(self) -> None:
        object.__setattr__(self, "reference_mole_fractions", tuple(self.reference_mole_fractions))


PhaseEquilibriumSpecification = FixedTPFeed | FixedTCoexistence


@dataclass(frozen=True, slots=True)
class PhaseHint:
    """Non-authoritative start-order hint bound to one phase index."""

    phase_index: int
    mole_fractions: tuple[float, ...]
    molar_density_mol_m3: float

    def __post_init__(self) -> None:
        object.__setattr__(self, "mole_fractions", tuple(self.mole_fractions))


@dataclass(frozen=True, slots=True)
class PhaseEquilibriumProblem:
    """One typed HELD-lineage feed or coexistence specification."""

    specification: PhaseEquilibriumSpecification
    hints: tuple[PhaseHint, ...] = ()

    def __post_init__(self) -> None:
        object.__setattr__(self, "hints", tuple(self.hints))


@dataclass(frozen=True, slots=True)
class EquilibriumPhaseState:
    """One intensive phase state on its declared component support."""

    role: PhaseRole
    support: Literal["all_components", "neutral_components"]
    supported_component_ids: tuple[str, ...]
    mole_fractions: tuple[float, ...]
    molar_volume_m3_per_mol: float
    molar_density_mol_m3: float
    packing_fraction: float
    pressure_pa: float
    chemical_potential_component_ids: tuple[str, ...]
    chemical_potential_over_rt: tuple[float, ...]
    mechanical_class: str


@dataclass(frozen=True)
class PhaseEquilibriumDiagnostics:
    """Route-native numerical, physical, support, and selection evidence."""

    route: str
    solver_status: str
    numerical_status: str
    physical_status: str
    failure_stage: str | None
    failure_reason: str | None
    certificate_scope: str
    phase_selection_claim: str
    globality_claim: str
    search_scope: str
    search_budget: str
    search_exclusions: tuple[str, ...]
    completeness_criterion: str
    role_evidence: tuple[str, ...]
    optimizer_solve_count: int
    optimizer_iteration_count: int
    provider_evaluation_count: int
    hint_count: int
    hint_reordered_start_count: int
    material_balance_inf: float
    electroneutrality_inf: float
    pressure_residual_inf: float
    chemical_potential_residual_inf: float
    kkt_stationarity_inf: float
    dual_sign_violation_inf: float
    complementarity_inf: float
    phase_distinctness: float


@dataclass(frozen=True)
class PhaseEquilibriumSensitivity:
    """Exact first-order derivative of one accepted fixed-topology phase result."""

    status: Literal["available", "unavailable"]
    failure_reason: str | None
    derivative_order: Literal[1]
    higher_order_status: Literal["unavailable"]
    higher_order_reason: str
    tensor_layout: Literal["row-major"]
    parameters: tuple[epcsaft.ActiveParameter, ...]
    residual_identities: tuple[str, ...]
    state_identities: tuple[str, ...]
    output_identities: tuple[str, ...]
    residual_state_jacobian: tuple[tuple[float, ...], ...]
    residual_parameter_jacobian: tuple[tuple[float, ...], ...]
    state_parameter_jacobian: tuple[tuple[float, ...], ...]
    output_parameter_jacobian: tuple[tuple[float, ...], ...]
    rank: int
    dimension: int
    scaled_condition_number_2: float | None
    linear_solve_backward_error: float | None
    reference_phase_index: int


@dataclass(frozen=True)
class _HeldDiscoveryResult:
    """Internal certified result from the HELD or HELD2 controller."""

    phases: tuple[_HeldPhaseState, ...]
    phase_fractions: tuple[float, ...]
    total_free_energy_over_rt: float
    parameter_fingerprint: str
    diagnostics: HeldDiagnostics | Held2Diagnostics
    sensitivity: PhaseEquilibriumSensitivity | None


class _HeldDiscoveryError(RuntimeError):
    """Raised when internal HELD discovery cannot return an accepted result."""

    def __init__(
        self,
        message: str,
        diagnostics: HeldDiagnostics | Held2Diagnostics,
    ) -> None:
        super().__init__(message)
        self.diagnostics = diagnostics


@dataclass(frozen=True, slots=True)
class FeedSplit:
    """Extensive fixed-feed split; amounts and fractions preserve phase order."""

    feed_amounts_mol: tuple[float, ...]
    phase_amounts_mol: tuple[float, ...]
    phase_fractions: tuple[float, ...]
    total_free_energy_over_rt: float


@dataclass(frozen=True, slots=True)
class IncipientBoundary:
    """Two-phase boundary without a thermodynamically determined split fraction."""

    reference_phase_index: int
    pressure_pa: float


EquilibriumExtent = FeedSplit | IncipientBoundary


@dataclass(frozen=True, slots=True)
class PhaseEquilibriumResult:
    """Certified phase states plus the typed extensive or boundary meaning."""

    problem: PhaseEquilibriumProblem
    phases: tuple[EquilibriumPhaseState, ...]
    extent: EquilibriumExtent
    parameter_fingerprint: str
    diagnostics: PhaseEquilibriumDiagnostics
    sensitivity: PhaseEquilibriumSensitivity | None


class PhaseEquilibriumError(RuntimeError):
    """Raised when HELD fails or its retained state violates the role scope."""

    def __init__(
        self,
        message: str,
        diagnostics: PhaseEquilibriumDiagnostics,
        best_candidate: tuple[EquilibriumPhaseState, ...],
    ) -> None:
        super().__init__(message)
        self.diagnostics = diagnostics
        self.best_candidate = best_candidate


@dataclass(frozen=True)
class ChemicalEquilibriumConstant:
    """One dimensionless reaction constant and its immutable provenance."""

    ln_value: float
    source_id: str
    reference_id: str
    reaction_orientation: str
    conversion_id: str
    dimensionless: bool


@dataclass(frozen=True)
class ChemicalStandardState:
    """Source activity convention and p° provenance transformed at each trial P."""

    id: str
    activity_scale_id: str
    log_activity_scale_factors: tuple[float, ...]
    reference_pressure_pa: float
    source_reference_component_ids: tuple[str, ...]
    source_reference_solvent_composition: tuple[float, ...]
    source_reference_ion_pairs: tuple[epcsaft.SourceReferenceIonPair, ...]
    source_reference_phase: Literal["liquid"]
    source_reference_convention: Literal["molality-infinite-dilution"]
    source_reference_activity_convention_id: str
    source_reference_standard_molality_mol_per_kg: float

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "log_activity_scale_factors",
            tuple(self.log_activity_scale_factors),
        )
        object.__setattr__(
            self,
            "source_reference_component_ids",
            tuple(self.source_reference_component_ids),
        )
        object.__setattr__(
            self,
            "source_reference_solvent_composition",
            tuple(self.source_reference_solvent_composition),
        )
        object.__setattr__(
            self, "source_reference_ion_pairs", tuple(self.source_reference_ion_pairs)
        )


@dataclass(frozen=True)
class IdealGasPhase:
    """Source-bound ideal-gas restriction used for nonpredictive benchmarks."""

    model_fingerprint: str
    reference_id: str


@dataclass(frozen=True)
class ProviderPhase:
    """Installed Provider model and application-admitted physical packing domain."""

    model: epcsaft.Mixture
    admissible_packing_fraction_interval: tuple[float, float]


@dataclass(frozen=True)
class ChemicalEquilibriumProblem:
    """One reaction system with an explicit strict-interior admission boundary."""

    species_ids: tuple[str, ...]
    charges: tuple[int, ...]
    molar_masses_kg_per_mol: tuple[float, ...]
    balance_matrix: tuple[tuple[float, ...], ...]
    conserved_totals: tuple[float, ...]
    reaction_matrix: tuple[tuple[float, ...], ...]
    feed_amounts_mol: tuple[float, ...]
    equilibrium_constants: tuple[ChemicalEquilibriumConstant, ...]
    strict_interior_amount_floor_mol: float
    source_standard_state: ChemicalStandardState | None = None


@dataclass(frozen=True)
class ChemicalEquilibriumAttempt:
    """One deterministic primary or recovery terminal in the finite search."""

    ordinal: int
    primary_ordinal: int
    kind: str
    parent_ordinal: int | None
    start_identity: str
    start_construction_status: str
    retraction_status: str
    provider_domain_status: str
    solver_status: str
    callback_error: str
    terminal_status: str
    amounts_mol: tuple[float, ...]
    volume_m3: float | None
    objective: float | None
    balance_inf_norm: float | None
    charge_inf_norm: float | None
    pressure_relative_residual: float | None
    reaction_affinity_inf_norm: float | None
    kkt_stationarity_inf_norm: float | None
    complementarity_inf_norm: float | None
    kkt_dimension: int
    kkt_rank: int
    condition_number_inf: float | None
    local_minimum_status: str
    trace_status: str
    basin_ordinal: int | None
    recovery_seed_count: int
    recovery_solve_count: int


@dataclass(frozen=True)
class ChemicalEquilibriumBasin:
    """One materially distinct certified strict local minimum."""

    ordinal: int
    representative_attempt_ordinal: int
    amounts_mol: tuple[float, ...]
    volume_m3: float
    objective: float


@dataclass(frozen=True)
class ChemicalEquilibriumBudgetPrefix:
    """A projection of the fixed run at one nested primary-start budget."""

    primary_budget: int
    attempted_primary_ordinals: tuple[int, ...]
    basin_ordinals: tuple[int, ...]
    selected_basin_ordinal: int | None
    selection_changed: bool


@dataclass(frozen=True)
class ChemicalEquilibriumSearch:
    """Complete deterministic finite-basin search receipt."""

    status: str
    primary_budget: int
    primary_attempt_count: int
    generated_start_count: int
    budget_truncated_start_count: int
    duplicate_start_count: int
    infeasible_start_count: int
    evaluated_start_count: int
    domain_rejected_start_count: int
    construction_rejected_start_count: int
    attempts: tuple[ChemicalEquilibriumAttempt, ...]
    basins: tuple[ChemicalEquilibriumBasin, ...]
    budget_prefixes: tuple[ChemicalEquilibriumBudgetPrefix, ...]
    selected_basin_ordinal: int | None
    selected_objective: float | None
    selection_label: str


@dataclass(frozen=True)
class ChemicalEquilibriumCriterion:
    """One named scalar certification check with its declared threshold."""

    name: str
    status: str
    value: float
    limit: float
    comparison: str


@dataclass(frozen=True)
class ChemicalEquilibriumDiagnostics:
    """Independent local chemical-equilibrium admission axes and residuals."""

    solver_status: str
    callback_error: str
    failure_kind: str
    failure_reason: str
    chemical_certification_level: str
    boundary_status: str
    structural_zero_species_indices: tuple[int, ...]
    numerical_status: str
    physical_status: str
    predictive_status: str
    provider_domain_status: str
    local_minimum_status: str
    negative_curvature_recovery_status: str
    negative_curvature_recovery_attempts: int
    negative_curvature_recovery_selected_sign: int
    trace_status: str
    globality_status: str
    search: ChemicalEquilibriumSearch
    amounts_mol: tuple[float, ...]
    volume_m3: float | None
    balance_inf_norm: float | None
    charge_inf_norm: float | None
    pressure_relative_residual: float | None
    reaction_affinity_inf_norm: float | None
    reaction_affinity_residuals: tuple[float, ...]
    packing_fraction: float | None
    packing_fraction_min: float | None
    packing_fraction_max: float | None
    total_ion_fraction: float | None
    total_ion_fraction_max: float | None
    minimum_amount_mol: float | None
    trace_floor_mol: float | None
    kkt_stationarity_inf_norm: float | None
    physical_stationarity_residuals: tuple[float, ...]
    physical_equality_multipliers: tuple[float, ...]
    physical_constraint_jacobian: tuple[float, ...]
    physical_constraint_shape: tuple[int, int]
    physical_lagrangian_gradient: tuple[float, ...]
    physical_objective_gradient: tuple[float, ...]
    physical_to_chart_jacobian: tuple[float, ...]
    physical_to_chart_shape: tuple[int, int]
    chart_physical_pullback_residual_inf_norm: float | None
    complementarity_inf_norm: float | None
    numerical_criteria: tuple[ChemicalEquilibriumCriterion, ...]
    physical_criteria: tuple[ChemicalEquilibriumCriterion, ...]
    first_failed_numerical_criterion: str | None
    first_failed_physical_criterion: str | None
    kkt_dimension: int
    kkt_rank: int
    condition_number_inf: float | None
    active_lower_bounds: tuple[int, ...]
    active_upper_bounds: tuple[int, ...]
    active_constraint_bounds: tuple[int, ...]
    active_trace_species: tuple[int, ...]
    chart_topology: str
    reduced_hessian_status: str
    reduced_hessian: tuple[float, ...]
    reduced_hessian_nullspace_basis: tuple[float, ...]
    reduced_hessian_nullspace_shape: tuple[int, int]
    reduced_hessian_eigenvalues: tuple[float, ...]
    reduced_hessian_spectrum_status: str
    reduced_hessian_raw_inertia: tuple[int, int, int]
    reduced_hessian_inertia: tuple[int, int, int]
    reduced_hessian_scale: float | None
    reduced_hessian_eigenvalue_tolerance: float | None
    objective_gradient: tuple[float, ...]
    constraint_values: tuple[float, ...]
    constraint_jacobian: tuple[float, ...]
    lagrangian_gradient: tuple[float, ...]
    equality_multipliers: tuple[float, ...]
    chart_stationarity_inf_norm: float | None
    lagrangian_hessian: tuple[float, ...]
    covariant_lagrangian_hessian: tuple[float, ...]
    derivative_check_step: float | None
    objective_gradient_check_relative_error: float | None
    physical_objective_gradient_check_relative_error: float | None
    constraint_jacobian_check_relative_error: float | None
    lagrangian_hessian_check_relative_error: float | None
    reduced_hessian_check_relative_error: float | None
    derivative_check_worst_entry: str
    derivative_check_worst_relative_error: float | None
    derivative_check_worst_analytic_value: float | None
    derivative_check_worst_finite_difference_value: float | None
    derivative_check_worst_step: float | None
    derivative_coordinate_order: tuple[str, ...]
    derivative_objective_basis: str
    derivative_constraint_basis: str
    derivative_constraint_order: tuple[str, ...]
    objective_gradient_shape: tuple[int, ...]
    constraint_jacobian_shape: tuple[int, int]
    lagrangian_hessian_shape: tuple[int, int]
    kkt_root_jacobian: tuple[float, ...]
    kkt_root_shape: tuple[int, int]
    kkt_root_status: str
    kkt_root_jacobian_check_relative_error: float | None
    reference_representation_residual_inf_norm: float | None
    reference_convergence_error: float | None
    reference_derivative_availability: int | None


@dataclass(frozen=True)
class ChemicalEquilibriumSensitivityRequest:
    """Request exact operation columns and ordered installed-Provider coordinates."""

    active_parameters: epcsaft.ActiveParameterSet | None = None


@dataclass(frozen=True)
class ChemicalEquilibriumSensitivityParameter:
    """One ordered scalar input and the units of its returned state derivatives."""

    name: str
    kind: str
    source_index: int | None
    input_unit: str
    amount_derivative_unit: str
    volume_derivative_unit: str
    identity: str | None = None
    component_ids: tuple[str, ...] = ()
    value: float | None = None


@dataclass(frozen=True)
class ChemicalEquilibriumSensitivity:
    """Exact conditioned state-input Jacobian or fail-closed unavailability evidence."""

    status: str
    failure_reason: str
    parameters: tuple[ChemicalEquilibriumSensitivityParameter, ...]
    amount_state_order: tuple[str, ...]
    amount_derivatives: tuple[tuple[float, ...], ...]
    volume_derivatives: tuple[float, ...]
    kkt_dimension: int
    kkt_rank: int
    condition_number_inf: float
    active_lower_bounds: tuple[int, ...]
    active_upper_bounds: tuple[int, ...]
    active_constraint_bounds: tuple[int, ...]
    active_trace_species: tuple[int, ...]
    chart_topology: str
    parameter_fingerprint: str
    provider_parameter_status: str
    provider_parameter_failure_reason: str
    reference_parameter_status: str
    reference_parameter_failure_reason: str


@dataclass(frozen=True)
class ChemicalArtifactIdentity:
    """Installed distributions consumed by one chemical result."""

    equilibrium_distribution: str
    equilibrium_version: str
    equilibrium_record_sha256: str
    provider_distribution: str | None
    provider_version: str | None
    provider_record_sha256: str | None
    provider_transport_capsule_name: str | None


@dataclass(frozen=True)
class ChemicalEquilibriumResult:
    """Certified local fixed-T,P homogeneous chemical-equilibrium value."""

    temperature_k: float
    pressure_pa: float
    species_ids: tuple[str, ...]
    charges: tuple[int, ...]
    amounts_mol: tuple[float, ...]
    mole_fractions: tuple[float, ...]
    volume_m3: float
    thermodynamic_model: str
    model_fingerprint: str
    provider_component_ids: tuple[str, ...] | None
    provider_parameter_fingerprint: str | None
    equilibrium_constants: tuple[ChemicalEquilibriumConstant, ...]
    source_standard_state: ChemicalStandardState | None
    source_reference_transfer: epcsaft.SourceReferenceTransfer | None
    provider_reference_id: str | None
    standard_offsets: tuple[float, ...] | None
    ln_k_provider_basis: tuple[float, ...] | None
    local_scope: str
    diagnostics: ChemicalEquilibriumDiagnostics
    response_kind: str
    artifact_identity: ChemicalArtifactIdentity
    sensitivity: ChemicalEquilibriumSensitivity | None


class ChemicalEquilibriumError(RuntimeError):
    """Raised when no typed, certified local chemical state can be admitted."""

    def __init__(self, message: str, diagnostics: ChemicalEquilibriumDiagnostics) -> None:
        super().__init__(message)
        self.diagnostics = diagnostics


_GAS_CONSTANT_J_MOL_K = 8.31446261815324
_LIQUID_COMPRESSIBILITY_FACTOR_MAXIMUM = 0.25
_VAPOR_COMPRESSIBILITY_FACTOR_MINIMUM = 0.45
_HELD2_PACKING_FRACTION_MINIMUM = 1.0e-6
_HELD2_PACKING_FRACTION_MAXIMUM = 0.74


def _float_tuple(payload: object) -> tuple[float, ...]:
    if isinstance(payload, (int, float)):
        return (float(payload),)
    values = cast(list[float] | tuple[float, ...], payload)
    return tuple(float(value) for value in values)


def _vector(payload: object, size: int, name: str) -> tuple[float, ...]:
    values = _float_tuple(payload)
    if len(values) != size:
        raise RuntimeError(f"native {name} has the wrong size")
    return values


def _phase(payload: Mapping[str, object]) -> _HeldPhaseState:
    mole_fractions = _float_tuple(payload.get("mole_fractions", (1.0,)))
    chemical_potential = _vector(
        payload["chemical_potential_over_rt"],
        len(mole_fractions),
        "phase chemical-potential vector",
    )
    return _HeldPhaseState(
        amount_mol=float(cast(float, payload["amount_mol"])),
        mole_fractions=mole_fractions,
        volume_m3=float(cast(float, payload["volume_m3"])),
        molar_density_mol_m3=float(cast(float, payload["molar_density_mol_m3"])),
        pressure_pa=float(cast(float, payload["pressure_pa"])),
        chemical_potential_over_rt=chemical_potential,
    )


def _quantity(quantity: object, units: str, name: str, operation: str) -> float:
    if not isinstance(quantity, Quantity):
        raise TypeError(f"{operation} requires a Pint {name} quantity")
    try:
        value = float(quantity.to(units).magnitude)
    except DimensionalityError as error:
        raise ValueError(f"{name} units must be convertible to {units}") from error
    except (TypeError, ValueError) as error:
        raise TypeError(f"{operation} requires one scalar Pint {name} quantity") from error
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be positive and finite")
    return value


def _phase_feed(
    feed_amounts_mol: Sequence[float],
    component_count: int,
) -> tuple[float, ...]:
    if isinstance(feed_amounts_mol, (str, bytes)):
        raise TypeError("feed amounts must be a numeric sequence")
    try:
        values = tuple(float(value) for value in feed_amounts_mol)
    except (TypeError, ValueError) as error:
        raise TypeError("feed amounts must be a numeric sequence") from error
    if len(values) != component_count:
        raise ValueError(f"phase_equilibrium requires exactly {component_count} components")
    if not all(math.isfinite(value) and value > 0.0 for value in values):
        raise ValueError("feed amounts must be positive and finite")
    return values


def _optional_float(payload: Mapping[str, object], name: str) -> float | None:
    value = payload[name]
    return None if value is None else float(cast(float, value))


def _held_diagnostics(payload: Mapping[str, object]) -> HeldDiagnostics:
    profiles = cast(Sequence[object], payload["search_profiles"])
    statuses = tuple(
        payload[name] for name in ("solver_status", "numerical_status", "physical_status")
    )
    if not all(
        isinstance(status, str) and status in {"passed", "failed", "not_adjudicated"}
        for status in statuses
    ):
        raise ValueError("native HELD payload has an invalid status vocabulary")
    typed_statuses = cast(tuple[str, str, str], statuses)
    return HeldDiagnostics(
        outcome=str(payload["outcome"]),
        search_status=str(payload["search_status"]),
        solver_status=typed_statuses[0],
        numerical_status=typed_statuses[1],
        physical_status=typed_statuses[2],
        attempts=int(cast(int, payload["attempts"])),
        major_iterations=int(cast(int, payload["major_iterations"])),
        material_balance_max_abs=_optional_float(payload, "material_balance_max_abs"),
        pressure_stationarity_max_relative=_optional_float(
            payload, "pressure_stationarity_max_relative"
        ),
        kkt_stationarity_max_abs=_optional_float(payload, "kkt_stationarity_max_abs"),
        chemical_potential_max_relative=_optional_float(payload, "chemical_potential_max_relative"),
        search_profiles=tuple(str(profile) for profile in profiles),
        globality_certificate=str(payload["globality_certificate"]),
        provider_evaluation_count=int(cast(int, payload["provider_evaluation_count"])),
        failure_reason=str(payload["failure_reason"]),
    )


def _failed_held_diagnostics(outcome: str, search_status: str, reason: str) -> HeldDiagnostics:
    return HeldDiagnostics(
        outcome=outcome,
        search_status=search_status,
        solver_status="not_adjudicated",
        numerical_status="not_adjudicated",
        physical_status="not_adjudicated",
        attempts=0,
        major_iterations=0,
        material_balance_max_abs=None,
        pressure_stationarity_max_relative=None,
        kkt_stationarity_max_abs=None,
        chemical_potential_max_relative=None,
        search_profiles=(),
        globality_certificate="not_guaranteed",
        provider_evaluation_count=0,
        failure_reason=reason,
    )


def _held2_diagnostics(payload: Mapping[str, object]) -> Held2Diagnostics:
    outcome = str(payload["outcome"])
    accepted = outcome == "physical_equilibrium_accepted"
    one_phase = outcome == "one_phase_no_negative_witness_detected"
    failure = payload.get("failure_reason")
    completed = accepted or one_phase
    adjudication_status = "passed" if completed else "not_adjudicated"
    physical_status = (
        "passed"
        if completed
        else "failed"
        if payload["material_balance_max_abs"] is not None
        else "not_adjudicated"
    )
    return Held2Diagnostics(
        outcome=outcome,
        failure_stage=(
            None if payload.get("failure_stage") is None else str(payload["failure_stage"])
        ),
        solver_status=adjudication_status,
        numerical_status=adjudication_status,
        physical_status=physical_status,
        optimizer_solves=int(cast(int, payload["optimizer_solves"])),
        optimizer_iterations=int(cast(int, payload["optimizer_iterations"])),
        major_iterations=int(cast(int, payload["major_iterations"])),
        material_balance_max_abs=_optional_float(payload, "material_balance_max_abs"),
        electroneutrality_inf=_optional_float(payload, "electroneutrality_inf"),
        pressure_stationarity_max_relative=_optional_float(
            payload, "pressure_stationarity_max_relative"
        ),
        kkt_stationarity_max_abs=_optional_float(payload, "kkt_stationarity_max_abs"),
        globality_certificate=str(payload["globality_certificate"]),
        provider_evaluation_count=int(cast(int, payload["provider_evaluation_count"])),
        failure_reason="" if failure is None else str(failure),
    )


def _held2_phase_state(
    component_count: int,
    payload: Mapping[str, object],
) -> _HeldPhaseState:
    fractions = _vector(payload["mole_fractions"], component_count, "HELD2 phase composition")
    phase_fraction = float(cast(float, payload["phase_fraction"]))
    molar_volume = float(cast(float, payload["molar_volume_m3_mol"]))
    if not math.isfinite(phase_fraction) or phase_fraction <= 0.0:
        raise ValueError("native HELD2 phase fraction must be positive and finite")
    if not math.isfinite(molar_volume) or molar_volume <= 0.0:
        raise ValueError("native HELD2 phase volume must be positive and finite")
    pressure_pa = float(cast(float, payload["pressure_pa"]))
    if not math.isfinite(pressure_pa):
        raise ValueError("native HELD2 phase pressure must be finite")
    return _HeldPhaseState(
        amount_mol=phase_fraction,
        mole_fractions=fractions,
        volume_m3=phase_fraction * molar_volume,
        molar_density_mol_m3=1.0 / molar_volume,
        pressure_pa=pressure_pa,
        chemical_potential_over_rt=_vector(
            payload["chemical_potential_over_rt"],
            component_count,
            "HELD2 chemical-potential vector",
        ),
    )


def _held2_result(
    feed: tuple[float, ...],
    fingerprint: str,
    native: Mapping[str, object],
    active_parameters: epcsaft.ActiveParameterSet | None,
) -> _HeldDiscoveryResult:
    diagnostics = _held2_diagnostics(native)
    if diagnostics.globality_certificate != "not_guaranteed":
        raise ValueError("native HELD2 result has an invalid globality certificate")
    if diagnostics.outcome not in {
        "one_phase_no_negative_witness_detected",
        "physical_equilibrium_accepted",
    }:
        reason = diagnostics.failure_reason or "HELD2 did not return a certified equilibrium"
        raise _HeldDiscoveryError(reason, diagnostics)
    if str(native["parameter_fingerprint"]) != fingerprint:
        raise ValueError("native HELD2 result has the wrong Provider fingerprint")
    sensitivity = _phase_equilibrium_sensitivity(native.get("sensitivity"), active_parameters)

    phase_payloads = cast(Sequence[Mapping[str, object]], native["phases"])
    phases = tuple(_held2_phase_state(len(feed), payload) for payload in phase_payloads)
    phase_fractions = tuple(phase.amount_mol for phase in phases)
    expected_phase_count = 1 if diagnostics.outcome.startswith("one_phase") else 2
    if len(phases) != expected_phase_count or not math.isclose(
        math.fsum(phase_fractions), 1.0, rel_tol=0.0, abs_tol=1.0e-10
    ):
        raise ValueError("native HELD2 phase result is incomplete")
    total_gibbs = float(cast(float, native["total_free_energy_over_rt"]))
    if not math.isfinite(total_gibbs):
        raise ValueError("native HELD2 objective must be finite")
    return _HeldDiscoveryResult(
        phases=phases,
        phase_fractions=phase_fractions,
        total_free_energy_over_rt=total_gibbs,
        parameter_fingerprint=fingerprint,
        diagnostics=diagnostics,
        sensitivity=sensitivity,
    )


def _row_major_matrix(
    value: object,
    rows: int,
    columns: int,
    field: str,
    *,
    unavailable: bool,
) -> tuple[tuple[float, ...], ...]:
    flat = _float_tuple(value)
    if unavailable and not flat:
        return ()
    if len(flat) != rows * columns:
        raise ValueError(f"native {field} has the wrong row-major shape")
    if any(not math.isfinite(entry) for entry in flat):
        raise ValueError(f"native {field} must be finite")
    return tuple(flat[row * columns : (row + 1) * columns] for row in range(rows))


def _phase_equilibrium_sensitivity(
    value: object,
    active_parameters: epcsaft.ActiveParameterSet | None,
) -> PhaseEquilibriumSensitivity | None:
    if active_parameters is None:
        if value is not None:
            raise ValueError("native phase sensitivity was not requested")
        return None
    if not isinstance(value, Mapping):
        raise TypeError("native phase sensitivity is missing")
    status = str(value["status"])
    if status not in {"available", "unavailable"}:
        raise ValueError("native phase sensitivity status is invalid")
    unavailable = status == "unavailable"
    failure = value["failure_reason"]
    if (unavailable and not isinstance(failure, str)) or (not unavailable and failure is not None):
        raise ValueError("native phase sensitivity failure state is invalid")
    parameter_identities = tuple(
        str(item) for item in cast(Sequence[object], value["parameter_identities"])
    )
    if parameter_identities != active_parameters.identities:
        raise ValueError("native phase sensitivity parameter identity or order changed")
    derivative_order = int(cast(int, value["derivative_order"]))
    higher_order_status = str(value["higher_order_status"])
    tensor_layout = str(value["tensor_layout"])
    if derivative_order != 1 or higher_order_status != "unavailable":
        raise ValueError("native phase sensitivity derivative-order declaration is invalid")
    if tensor_layout != "row-major":
        raise ValueError("native phase sensitivity tensor layout is invalid")
    dimension = int(cast(int, value["dimension"]))
    parameter_count = int(cast(int, value["parameter_count"]))
    output_count = int(cast(int, value["output_count"]))
    rank = int(cast(int, value["rank"]))
    if (
        parameter_count != len(active_parameters.coordinates)
        or min(dimension, output_count, rank) < 0
    ):
        raise ValueError("native phase sensitivity dimensions are invalid")
    residual_identities = tuple(
        str(item) for item in cast(Sequence[object], value["residual_identities"])
    )
    state_identities = tuple(
        str(item) for item in cast(Sequence[object], value["state_identities"])
    )
    output_identities = tuple(
        str(item) for item in cast(Sequence[object], value["output_identities"])
    )
    if len(residual_identities) != dimension or len(state_identities) != dimension:
        raise ValueError("native phase sensitivity residual/state identities are incomplete")
    if len(output_identities) != output_count:
        raise ValueError("native phase sensitivity output identities are incomplete")
    condition_value = value["scaled_condition_number_2"]
    backward_value = value["linear_solve_backward_error"]
    condition = None if condition_value is None else float(cast(float, condition_value))
    backward = None if backward_value is None else float(cast(float, backward_value))
    if not unavailable and (
        rank != dimension
        or condition is None
        or not math.isfinite(condition)
        or condition > 1.0e8
        or backward is None
        or not math.isfinite(backward)
    ):
        raise ValueError("native phase sensitivity linear solve is uncertified")
    return PhaseEquilibriumSensitivity(
        status=cast(Literal["available", "unavailable"], status),
        failure_reason=None if failure is None else str(failure),
        derivative_order=cast(Literal[1], derivative_order),
        higher_order_status=cast(Literal["unavailable"], higher_order_status),
        higher_order_reason=str(value["higher_order_reason"]),
        tensor_layout=cast(Literal["row-major"], tensor_layout),
        parameters=active_parameters.coordinates,
        residual_identities=residual_identities,
        state_identities=state_identities,
        output_identities=output_identities,
        residual_state_jacobian=_row_major_matrix(
            value["residual_state_jacobian"],
            dimension,
            dimension,
            "residual-state Jacobian",
            unavailable=unavailable,
        ),
        residual_parameter_jacobian=_row_major_matrix(
            value["residual_parameter_jacobian"],
            dimension,
            parameter_count,
            "residual-parameter Jacobian",
            unavailable=unavailable,
        ),
        state_parameter_jacobian=_row_major_matrix(
            value["state_parameter_jacobian"],
            dimension,
            parameter_count,
            "state-parameter Jacobian",
            unavailable=unavailable,
        ),
        output_parameter_jacobian=_row_major_matrix(
            value["output_parameter_jacobian"],
            output_count,
            parameter_count,
            "output-parameter Jacobian",
            unavailable=unavailable,
        ),
        rank=rank,
        dimension=dimension,
        scaled_condition_number_2=condition,
        linear_solve_backward_error=backward,
        reference_phase_index=int(cast(int, value["reference_phase_index"])),
    )


def _held_discovery(
    model: epcsaft.Mixture,
    temperature_k: float,
    pressure_pa: float,
    feed: tuple[float, ...],
    active_parameters: epcsaft.ActiveParameterSet | None = None,
) -> _HeldDiscoveryResult:
    """Run the internal bounded HELD or strong-electrolyte HELD2 controller."""

    try:
        capsule = _native_transport(model, active_parameters)
    except (TypeError, ValueError) as error:
        diagnostics = _failed_held_diagnostics("invalid_input", "input_rejected", str(error))
        raise _HeldDiscoveryError(str(error), diagnostics) from error

    try:
        native = cast(
            Mapping[str, object],
            _equilibrium._solve_held_discovery(
                capsule,
                temperature_k,
                pressure_pa,
                feed,
            ),
        )
    except (RuntimeError, TypeError, ValueError) as error:
        diagnostics = _failed_held_diagnostics("error", "native_exception", str(error))
        raise _HeldDiscoveryError(str(error), diagnostics) from error

    controller = native.get("controller")
    if controller == "held2_steps_1_to_10_v1":
        try:
            return _held2_result(
                feed,
                model.parameter_fingerprint,
                native,
                active_parameters,
            )
        except _HeldDiscoveryError:
            raise
        except (KeyError, TypeError, ValueError) as error:
            diagnostics = _failed_held_diagnostics("error", "payload_error", str(error))
            raise _HeldDiscoveryError(
                "native HELD2 payload does not match the typed contract", diagnostics
            ) from error

    if controller != "held_binary_v1":
        diagnostics = _failed_held_diagnostics(
            "error",
            "payload_error",
            f"unknown native HELD controller: {controller!r}",
        )
        raise _HeldDiscoveryError(
            "native HELD payload has an unknown controller",
            diagnostics,
        )

    try:
        diagnostics = _held_diagnostics(native)
        if diagnostics.globality_certificate != "not_guaranteed":
            raise ValueError("native HELD result has an invalid globality certificate")
        if diagnostics.outcome not in {"one_phase", "accepted"}:
            reason = diagnostics.failure_reason or "HELD search did not return an accepted result"
            raise _HeldDiscoveryError(reason, diagnostics)
        phase_payloads = cast(Sequence[Mapping[str, object]], native["phases"])
        phases = tuple(_phase(payload) for payload in phase_payloads)
        phase_fractions = _float_tuple(native["phase_fractions"])
        if len(phases) != len(phase_fractions):
            raise ValueError("native HELD phase count does not match its outcome")
        if not math.isclose(sum(phase_fractions), 1.0, rel_tol=0.0, abs_tol=1.0e-10):
            raise ValueError("native HELD phase fractions do not sum to one")
        if str(native["parameter_fingerprint"]) != model.parameter_fingerprint:
            raise ValueError("native HELD result has the wrong provider fingerprint")
        sensitivity = _phase_equilibrium_sensitivity(native.get("sensitivity"), active_parameters)
        return _HeldDiscoveryResult(
            phases=phases,
            phase_fractions=phase_fractions,
            total_free_energy_over_rt=float(cast(float, native["total_free_energy_over_rt"])),
            parameter_fingerprint=str(native["parameter_fingerprint"]),
            diagnostics=diagnostics,
            sensitivity=sensitivity,
        )
    except _HeldDiscoveryError:
        raise
    except (KeyError, TypeError, ValueError) as error:
        diagnostics = _failed_held_diagnostics("error", "payload_error", str(error))
        raise _HeldDiscoveryError(
            "native HELD payload does not match the typed contract", diagnostics
        ) from error


def _native_phase_diagnostics(
    payload: Mapping[str, object],
) -> PhaseEquilibriumDiagnostics:
    return PhaseEquilibriumDiagnostics(
        route=str(payload["route"]),
        solver_status=str(payload["solver_status"]),
        numerical_status=str(payload["numerical_status"]),
        physical_status=str(payload["physical_status"]),
        failure_stage=None if payload["failure_stage"] is None else str(payload["failure_stage"]),
        failure_reason=None
        if payload["failure_reason"] is None
        else str(payload["failure_reason"]),
        certificate_scope=str(payload["certificate_scope"]),
        phase_selection_claim=str(payload["phase_selection_claim"]),
        globality_claim=str(payload["globality_claim"]),
        search_scope=str(payload["search_scope"]),
        search_budget=str(payload["search_budget"]),
        search_exclusions=tuple(
            str(value) for value in cast(Sequence[object], payload["search_exclusions"])
        ),
        completeness_criterion=str(payload["completeness_criterion"]),
        role_evidence=tuple(
            str(value) for value in cast(Sequence[object], payload["role_evidence"])
        ),
        optimizer_solve_count=int(cast(int, payload["optimizer_solve_count"])),
        optimizer_iteration_count=int(cast(int, payload["optimizer_iteration_count"])),
        provider_evaluation_count=int(cast(int, payload["provider_evaluation_count"])),
        hint_count=int(cast(int, payload["hint_count"])),
        hint_reordered_start_count=int(cast(int, payload["hint_reordered_start_count"])),
        material_balance_inf=float(cast(float, payload["material_balance_inf"])),
        electroneutrality_inf=float(cast(float, payload["electroneutrality_inf"])),
        pressure_residual_inf=float(cast(float, payload["pressure_residual_inf"])),
        chemical_potential_residual_inf=float(
            cast(float, payload["chemical_potential_residual_inf"])
        ),
        kkt_stationarity_inf=float(cast(float, payload["kkt_stationarity_inf"])),
        dual_sign_violation_inf=float(cast(float, payload["dual_sign_violation_inf"])),
        complementarity_inf=float(cast(float, payload["complementarity_inf"])),
        phase_distinctness=float(cast(float, payload["phase_distinctness"])),
    )


def _native_equilibrium_phase(payload: Mapping[str, object]) -> EquilibriumPhaseState:
    role = str(payload["role"])
    support = str(payload["support"])
    if role not in {"liquid", "vapor"} or support not in {
        "all_components",
        "neutral_components",
    }:
        raise ValueError("native phase identity is invalid")
    component_ids = tuple(
        str(value) for value in cast(Sequence[object], payload["supported_component_ids"])
    )
    chemical_ids = tuple(
        str(value) for value in cast(Sequence[object], payload["chemical_potential_component_ids"])
    )
    return EquilibriumPhaseState(
        role=cast(PhaseRole, role),
        support=cast(Literal["all_components", "neutral_components"], support),
        supported_component_ids=component_ids,
        mole_fractions=_vector(payload["mole_fractions"], len(component_ids), "phase composition"),
        molar_volume_m3_per_mol=float(cast(float, payload["molar_volume_m3_per_mol"])),
        molar_density_mol_m3=float(cast(float, payload["molar_density_mol_m3"])),
        packing_fraction=float(cast(float, payload["packing_fraction"])),
        pressure_pa=float(cast(float, payload["pressure_pa"])),
        chemical_potential_component_ids=chemical_ids,
        chemical_potential_over_rt=_vector(
            payload["chemical_potential_over_rt"],
            len(chemical_ids),
            "phase chemical potentials",
        ),
        mechanical_class=str(payload["mechanical_class"]),
    )


def _phase_support_name(support: PhaseSupport) -> str:
    if type(support) is AllComponents:
        return "all_components"
    if type(support) is NeutralComponents:
        return "neutral_components"
    raise TypeError("phase support must be AllComponents or NeutralComponents")


def _support_active_parameters(
    parameters: epcsaft.Parameters,
    active_parameters: epcsaft.ActiveParameterSet | None,
) -> epcsaft.ActiveParameterSet | None:
    if active_parameters is None:
        return None
    selected_components = set(parameters.component_ids)
    coordinates = tuple(
        coordinate
        for coordinate in active_parameters.coordinates
        if coordinate.applicability is not None
        and (
            coordinate.applicability.scope == "model"
            or (
                bool(coordinate.applicability.component_ids)
                and all(
                    component in selected_components
                    for component in coordinate.applicability.component_ids
                )
            )
        )
    )
    if not coordinates:
        return None
    return epcsaft.ActiveParameterSet(
        parameters,
        tuple(coordinate.identity for coordinate in coordinates),
        values=tuple(coordinate.value for coordinate in coordinates),
    )


def _validated_phase_hints(
    hints: tuple[PhaseHint, ...],
    phase_count: int,
    component_counts: tuple[int, ...],
) -> tuple[tuple[int, ...], tuple[tuple[float, ...], ...], tuple[float, ...]]:
    if any(type(hint) is not PhaseHint for hint in hints):
        raise TypeError("hints must contain PhaseHint values")
    indices = tuple(hint.phase_index for hint in hints)
    if len(indices) != len(set(indices)) or any(
        index < 0 or index >= phase_count for index in indices
    ):
        raise ValueError("phase hint indices must be unique and in range")
    for hint in hints:
        if len(hint.mole_fractions) != component_counts[hint.phase_index]:
            raise ValueError("phase hint composition has the wrong component count")
        if not all(math.isfinite(value) and value > 0.0 for value in hint.mole_fractions):
            raise ValueError("phase hint composition must be positive and finite")
        if not math.isclose(math.fsum(hint.mole_fractions), 1.0, rel_tol=0.0, abs_tol=1.0e-12):
            raise ValueError("phase hint composition must sum to one within 1e-12")
        if not math.isfinite(hint.molar_density_mol_m3) or hint.molar_density_mol_m3 <= 0.0:
            raise ValueError("phase hint density must be positive and finite")
    return (
        indices,
        tuple(hint.mole_fractions for hint in hints),
        tuple(hint.molar_density_mol_m3 for hint in hints),
    )


def _protected_discovery_diagnostics(
    retained: HeldDiagnostics | Held2Diagnostics,
    hint_count: int,
    *,
    numerical_status: str | None = None,
    physical_status: str | None = None,
    failure_stage: str | None = None,
    failure_reason: str | None = None,
    phase_selection_claim: str,
    role_evidence: tuple[str, ...] = (),
) -> PhaseEquilibriumDiagnostics:
    search_exclusions: tuple[str, ...]
    if isinstance(retained, Held2Diagnostics):
        route = "held2_discovery"
        native_failure_stage = retained.failure_stage
        search_scope = "retained_held2_steps_1_to_10_finite_domain"
        search_budget = (
            f"optimizer_solves={retained.optimizer_solves};"
            f"optimizer_iterations={retained.optimizer_iterations}"
        )
        search_exclusions = (
            "three_or_more_phase_topologies",
            "reactive_phases",
            "solid_phases",
        )
        optimizer_solve_count = retained.optimizer_solves
        optimizer_iteration_count = retained.optimizer_iterations
    else:
        route = "held_discovery"
        native_failure_stage = None
        search_scope = "retained_held_profiles:" + ",".join(retained.search_profiles)
        search_budget = f"attempts={retained.attempts};major_iterations={retained.major_iterations}"
        search_exclusions = (
            "charged_topologies",
            "three_or_more_phase_topologies",
        )
        optimizer_solve_count = retained.attempts
        optimizer_iteration_count = retained.major_iterations
    return PhaseEquilibriumDiagnostics(
        route=route,
        solver_status=retained.solver_status,
        numerical_status=(
            retained.numerical_status if numerical_status is None else numerical_status
        ),
        physical_status=(retained.physical_status if physical_status is None else physical_status),
        failure_stage=(native_failure_stage if failure_stage is None else failure_stage),
        failure_reason=failure_reason,
        certificate_scope="protected_held_finite_discovery",
        phase_selection_claim=phase_selection_claim,
        globality_claim="not_established",
        search_scope=search_scope,
        search_budget=search_budget,
        search_exclusions=search_exclusions,
        completeness_criterion=(
            "retained_controller_terminal_certificate;global_completeness_not_established"
        ),
        role_evidence=role_evidence,
        optimizer_solve_count=optimizer_solve_count,
        optimizer_iteration_count=optimizer_iteration_count,
        provider_evaluation_count=retained.provider_evaluation_count,
        hint_count=hint_count,
        hint_reordered_start_count=0,
        material_balance_inf=retained.material_balance_max_abs or 0.0,
        electroneutrality_inf=(
            retained.electroneutrality_inf or 0.0 if isinstance(retained, Held2Diagnostics) else 0.0
        ),
        pressure_residual_inf=retained.pressure_stationarity_max_relative or 0.0,
        chemical_potential_residual_inf=(
            retained.chemical_potential_max_relative or 0.0
            if isinstance(retained, HeldDiagnostics)
            else 0.0
        ),
        kkt_stationarity_inf=retained.kkt_stationarity_max_abs or 0.0,
        dual_sign_violation_inf=0.0,
        complementarity_inf=0.0,
        phase_distinctness=0.0,
    )


def _discovery_result(
    model: epcsaft.Mixture,
    problem: PhaseEquilibriumProblem,
    specification: FixedTPFeed,
    active_parameters: epcsaft.ActiveParameterSet | None,
) -> PhaseEquilibriumResult:
    temperature_k = _quantity(
        specification.temperature, "kelvin", "temperature", "phase_equilibrium"
    )
    pressure_pa = _quantity(specification.pressure, "pascal", "pressure", "phase_equilibrium")
    feed = _phase_feed(specification.feed_amounts_mol, len(model.component_ids))
    total = math.fsum(feed)
    normalized = tuple(value / total for value in feed)
    try:
        retained = _held_discovery(
            model,
            temperature_k,
            pressure_pa,
            normalized,
            active_parameters,
        )
    except _HeldDiscoveryError as error:
        retained_diagnostics = error.diagnostics
        diagnostics = _protected_discovery_diagnostics(
            retained_diagnostics,
            len(problem.hints),
            physical_status="rejected",
            failure_stage=(
                retained_diagnostics.failure_stage
                if isinstance(retained_diagnostics, Held2Diagnostics)
                else "protected_discovery"
            ),
            failure_reason=str(error),
            phase_selection_claim="no_phase_selection_claim",
        )
        raise PhaseEquilibriumError(str(error), diagnostics, ()) from error
    sensitivity = retained.sensitivity
    if sensitivity is not None and total != 1.0:
        sensitivity = replace(
            sensitivity,
            residual_state_jacobian=tuple(
                tuple(value / total for value in row) for row in sensitivity.residual_state_jacobian
            ),
            state_parameter_jacobian=tuple(
                tuple(total * value for value in row)
                for row in sensitivity.state_parameter_jacobian
            ),
            output_parameter_jacobian=tuple(
                tuple(total * value for value in row)
                for row in sensitivity.output_parameter_jacobian
            ),
        )
    capsule = _native_transport(model, active_parameters)
    roles: list[PhaseRole] = []
    for phase in retained.phases:
        molar_volume = phase.volume_m3 / phase.amount_mol
        compressibility = pressure_pa * molar_volume / (_GAS_CONSTANT_J_MOL_K * temperature_k)
        role: PhaseRole = (
            "liquid"
            if compressibility <= _LIQUID_COMPRESSIBILITY_FACTOR_MAXIMUM
            else "vapor"
            if compressibility >= _VAPOR_COMPRESSIBILITY_FACTOR_MINIMUM
            else cast(PhaseRole, "")
        )
        if not role:
            diagnostics = _protected_discovery_diagnostics(
                retained.diagnostics,
                len(problem.hints),
                numerical_status="uncertified",
                physical_status="rejected",
                failure_stage="phase_role_certificate",
                failure_reason="retained HELD phase role is thermodynamically inconclusive",
                phase_selection_claim="no_phase_role_claim",
            )
            raise PhaseEquilibriumError(
                diagnostics.failure_reason or "phase role was not certified",
                diagnostics,
                (),
            )
        roles.append(role)
    try:
        phases = tuple(
            _native_equilibrium_phase(
                cast(
                    Mapping[str, object],
                    _equilibrium._audit_phase_state(
                        capsule,
                        role,
                        temperature_k,
                        retained_phase.mole_fractions,
                        retained_phase.volume_m3 / retained_phase.amount_mol,
                    ),
                )
            )
            for role, retained_phase in zip(roles, retained.phases, strict=True)
        )
    except (RuntimeError, TypeError, ValueError) as error:
        diagnostics = _protected_discovery_diagnostics(
            retained.diagnostics,
            len(problem.hints),
            numerical_status="uncertified",
            physical_status="rejected",
            failure_stage="phase_role_certificate",
            failure_reason=f"native phase-role audit failed: {error}",
            phase_selection_claim="no_phase_role_claim",
        )
        raise PhaseEquilibriumError(
            diagnostics.failure_reason or "phase role audit failed",
            diagnostics,
            (),
        ) from error
    if any(phase.mechanical_class != "strict_stable" for phase in phases):
        diagnostics = _protected_discovery_diagnostics(
            retained.diagnostics,
            len(problem.hints),
            numerical_status="uncertified",
            physical_status="rejected",
            failure_stage="phase_role_certificate",
            failure_reason="retained HELD phase is not mechanically strict-stable",
            phase_selection_claim="no_phase_role_claim",
        )
        raise PhaseEquilibriumError(
            diagnostics.failure_reason or "phase role was not certified",
            diagnostics,
            phases,
        )
    diagnostics = _protected_discovery_diagnostics(
        retained.diagnostics,
        len(problem.hints),
        phase_selection_claim="protected_held_finite_discovery_result",
        role_evidence=(
            "native_pressure_derivative_strict_stability",
            "compressibility_branch_with_fail_closed_gap",
        ),
    )
    return PhaseEquilibriumResult(
        problem=problem,
        phases=phases,
        extent=FeedSplit(
            feed_amounts_mol=feed,
            phase_amounts_mol=tuple(total * value for value in retained.phase_fractions),
            phase_fractions=retained.phase_fractions,
            total_free_energy_over_rt=total * retained.total_free_energy_over_rt,
        ),
        parameter_fingerprint=retained.parameter_fingerprint,
        diagnostics=diagnostics,
        sensitivity=sensitivity,
    )


def phase_equilibrium(
    model: epcsaft.Mixture,
    problem: PhaseEquilibriumProblem,
    *,
    active_parameters: epcsaft.ActiveParameterSet | None = None,
) -> PhaseEquilibriumResult:
    """Solve one typed HELD-lineage feed or two-phase coexistence specification."""

    if not isinstance(model, epcsaft.Mixture):
        raise TypeError("phase_equilibrium requires an epcsaft.Mixture")
    if type(problem) is not PhaseEquilibriumProblem:
        raise TypeError("problem must be a PhaseEquilibriumProblem")
    if active_parameters is not None and type(active_parameters) is not epcsaft.ActiveParameterSet:
        raise TypeError("active_parameters must be an epcsaft.ActiveParameterSet")
    specification = problem.specification
    if isinstance(specification, FixedTPFeed):
        if specification.topology is None:
            if problem.hints:
                raise ValueError("phase hints require a declared topology")
            return _discovery_result(model, problem, specification, active_parameters)
        if not isinstance(specification.topology, DeclaredTopology):
            raise TypeError("fixed-T,P topology must be DeclaredTopology or None")
        topology = specification.topology
        if not 1 <= len(topology.phases) <= 3:
            raise ValueError("fixed-T,P declared topology requires one to three phases")
        if any(type(phase.support) is not AllComponents for phase in topology.phases):
            raise ValueError("fixed-T,P declared topology currently requires AllComponents")
        if sum(phase.role == "vapor" for phase in topology.phases) > 1:
            raise ValueError("fixed-T,P declared topology permits at most one vapor phase")
        temperature_k = _quantity(
            specification.temperature, "kelvin", "temperature", "phase_equilibrium"
        )
        pressure_pa = _quantity(specification.pressure, "pascal", "pressure", "phase_equilibrium")
        feed = _phase_feed(specification.feed_amounts_mol, len(model.component_ids))
        hints = _validated_phase_hints(
            problem.hints,
            len(topology.phases),
            tuple(len(model.component_ids) for _ in topology.phases),
        )
        native = cast(
            Mapping[str, object],
            _equilibrium._solve_declared_fixed_tp(
                _native_transport(model, active_parameters),
                temperature_k,
                pressure_pa,
                feed,
                tuple(phase.role for phase in topology.phases),
                tuple(_phase_support_name(phase.support) for phase in topology.phases),
                *hints,
            ),
        )
        diagnostics = _native_phase_diagnostics(cast(Mapping[str, object], native["diagnostics"]))
        if str(native["parameter_fingerprint"]) != model.parameter_fingerprint:
            raise ValueError("native declared phase result has the wrong parameter fingerprint")
        native_phases = tuple(
            _native_equilibrium_phase(cast(Mapping[str, object], value))
            for value in cast(Sequence[object], native["phases"])
        )
        if diagnostics.physical_status != "accepted":
            raise PhaseEquilibriumError(
                diagnostics.failure_reason or "declared topology was not certified",
                diagnostics,
                native_phases,
            )
        extent = FeedSplit(
            feed_amounts_mol=feed,
            phase_amounts_mol=_float_tuple(native["phase_amounts_mol"]),
            phase_fractions=_float_tuple(native["phase_fractions"]),
            total_free_energy_over_rt=float(cast(float, native["total_free_energy_over_rt"])),
        )
        return PhaseEquilibriumResult(
            problem=problem,
            phases=native_phases,
            extent=extent,
            parameter_fingerprint=str(native["parameter_fingerprint"]),
            diagnostics=diagnostics,
            sensitivity=_phase_equilibrium_sensitivity(
                native.get("sensitivity"), active_parameters
            ),
        )
    if not isinstance(specification, FixedTCoexistence):
        raise TypeError("phase-equilibrium specification type is unsupported")
    if not isinstance(specification.topology, DeclaredTopology):
        raise TypeError("fixed-T coexistence topology must be DeclaredTopology")
    topology = specification.topology
    if len(topology.phases) != 2 or specification.reference_phase_index not in {0, 1}:
        raise ValueError("fixed-T coexistence requires exactly two phases and one reference index")
    temperature_k = _quantity(
        specification.temperature, "kelvin", "temperature", "phase_equilibrium"
    )
    reference = tuple(float(value) for value in specification.reference_mole_fractions)
    if len(reference) != len(model.component_ids):
        raise ValueError("reference composition has the wrong component count")
    reference_support = topology.phases[specification.reference_phase_index].support
    if type(reference_support) is not AllComponents:
        raise ValueError("coexistence reference phase must use AllComponents")
    partner_index = 1 - specification.reference_phase_index
    partner_support = topology.phases[partner_index].support
    if type(partner_support) is NeutralComponents:
        charges = tuple(model._model.charge_number)
        neutral_ids = tuple(
            component
            for component, charge in zip(model.component_ids, charges, strict=True)
            if charge == 0
        )
        if not neutral_ids:
            raise ValueError("neutral-component vapor support is empty")
        partner_parameters = model._parameters.select(neutral_ids)
        partner_model: epcsaft.Pure | epcsaft.Mixture = (
            epcsaft.Pure(partner_parameters)
            if len(neutral_ids) == 1
            else epcsaft.Mixture(partner_parameters)
        )
    elif type(partner_support) is AllComponents:
        partner_model = model
    else:
        raise TypeError("coexistence partner support is invalid")
    hints = _validated_phase_hints(
        problem.hints,
        2,
        tuple(
            len(model.component_ids)
            if index == specification.reference_phase_index
            else len(partner_model.component_ids)
            for index in range(2)
        ),
    )
    native = cast(
        Mapping[str, object],
        _equilibrium._solve_fixed_t_coexistence(
            _native_transport(model, active_parameters),
            _native_transport(
                partner_model,
                _support_active_parameters(partner_model._parameters, active_parameters),
            ),
            temperature_k,
            tuple(phase.role for phase in topology.phases),
            tuple(_phase_support_name(phase.support) for phase in topology.phases),
            specification.reference_phase_index,
            reference,
            *hints,
        ),
    )
    diagnostics = _native_phase_diagnostics(cast(Mapping[str, object], native["diagnostics"]))
    if str(native["parameter_fingerprint"]) != model.parameter_fingerprint:
        raise ValueError("native coexistence result has the wrong parameter fingerprint")
    native_phases = tuple(
        _native_equilibrium_phase(cast(Mapping[str, object], value))
        for value in cast(Sequence[object], native["phases"])
    )
    if diagnostics.physical_status != "accepted":
        raise PhaseEquilibriumError(
            diagnostics.failure_reason or "coexistence boundary was not certified",
            diagnostics,
            native_phases,
        )
    return PhaseEquilibriumResult(
        problem=problem,
        phases=native_phases,
        extent=IncipientBoundary(
            reference_phase_index=specification.reference_phase_index,
            pressure_pa=float(cast(float, native["pressure_pa"])),
        ),
        parameter_fingerprint=str(native["parameter_fingerprint"]),
        diagnostics=diagnostics,
        sensitivity=(
            None
            if active_parameters is None
            else _phase_equilibrium_sensitivity(native.get("sensitivity"), active_parameters)
        ),
    )


def _empty_chemical_search(status: str = "not_evaluated") -> ChemicalEquilibriumSearch:
    return ChemicalEquilibriumSearch(
        status=status,
        primary_budget=0,
        primary_attempt_count=0,
        generated_start_count=0,
        budget_truncated_start_count=0,
        duplicate_start_count=0,
        infeasible_start_count=0,
        evaluated_start_count=0,
        domain_rejected_start_count=0,
        construction_rejected_start_count=0,
        attempts=(),
        basins=(),
        budget_prefixes=(),
        selected_basin_ordinal=None,
        selected_objective=None,
        selection_label="lowest_observed_certified_local_value",
    )


_CHEMICAL_SEARCH_STATUSES = {
    "certified_local_minimum",
    "saddle_observed",
    "second_order_inconclusive",
    "boundary_unadjudicated",
    "domain_rejected",
    "no_feasible_start_found",
    "search_exhausted_no_certified_candidate",
    "infeasible_certified",
    "not_evaluated",
}
_CHEMICAL_TERMINAL_STATUSES = _CHEMICAL_SEARCH_STATUSES | {"solver_failed"}


def _validate_chemical_search(search: ChemicalEquilibriumSearch) -> None:
    if search.status not in _CHEMICAL_SEARCH_STATUSES:
        raise ValueError("native chemical search has an invalid status")
    if search.selection_label != "lowest_observed_certified_local_value":
        raise ValueError("native chemical search has an invalid selection label")
    if not 0 <= search.primary_attempt_count <= search.primary_budget <= 25:
        raise ValueError("native chemical search has an invalid primary budget")
    if (
        (search.status != "not_evaluated" and search.primary_attempt_count == 0)
        or search.generated_start_count
        != search.primary_attempt_count
        + search.budget_truncated_start_count
        + search.duplicate_start_count
        + search.infeasible_start_count
        or search.budget_truncated_start_count < 0
        or not 0 <= search.evaluated_start_count <= search.primary_attempt_count
        or not 0
        <= search.domain_rejected_start_count
        <= search.primary_attempt_count - search.evaluated_start_count
        or not 0 <= search.construction_rejected_start_count
        or search.evaluated_start_count
        + search.domain_rejected_start_count
        + search.construction_rejected_start_count
        != search.primary_attempt_count
    ):
        raise ValueError("native chemical search start accounting is inconsistent")
    if tuple(attempt.ordinal for attempt in search.attempts) != tuple(range(len(search.attempts))):
        raise ValueError("native chemical search attempt ordinals are inconsistent")
    if len(search.attempts) < search.primary_attempt_count or any(
        attempt.kind != "primary"
        or attempt.primary_ordinal != attempt.ordinal
        or attempt.parent_ordinal is not None
        for attempt in search.attempts[: search.primary_attempt_count]
    ):
        raise ValueError("native chemical search primary lineage is inconsistent")
    for attempt in search.attempts:
        if attempt.kind not in {"primary", "recovery"}:
            raise ValueError("native chemical search has an invalid attempt kind")
        if not 0 <= attempt.primary_ordinal < search.primary_attempt_count:
            raise ValueError("native chemical search primary ordinal is out of range")
        if attempt.kind == "recovery" and (
            attempt.parent_ordinal != attempt.primary_ordinal
            or attempt.parent_ordinal >= search.primary_attempt_count
        ):
            raise ValueError("native chemical search recovery lineage is inconsistent")
        if attempt.start_construction_status not in {
            "accepted",
            "rejected",
            "not_evaluated",
        }:
            raise ValueError("native chemical search has an invalid start status")
        if attempt.retraction_status not in {
            "not_needed",
            "passed",
            "failed",
            "balance_preserved_by_reaction_extent",
        }:
            raise ValueError("native chemical search has an invalid retraction status")
        if attempt.provider_domain_status not in {
            "passed",
            "failed",
            "not_applicable",
            "not_adjudicated",
        }:
            raise ValueError("native chemical attempt has an invalid Provider-domain status")
        if attempt.terminal_status not in _CHEMICAL_TERMINAL_STATUSES:
            raise ValueError("native chemical search has an invalid terminal status")
        if attempt.local_minimum_status not in {
            "passed",
            "saddle_observed",
            "second_order_inconclusive",
            "not_adjudicated",
        }:
            raise ValueError("native chemical attempt has an invalid curvature status")
        if attempt.trace_status not in {
            "interior",
            "at_or_below_floor",
            "not_adjudicated",
        }:
            raise ValueError("native chemical attempt has an invalid trace status")
        optional_values = (
            attempt.volume_m3,
            attempt.objective,
            attempt.balance_inf_norm,
            attempt.charge_inf_norm,
            attempt.pressure_relative_residual,
            attempt.reaction_affinity_inf_norm,
            attempt.kkt_stationarity_inf_norm,
            attempt.complementarity_inf_norm,
            attempt.condition_number_inf,
        )
        if any(value is not None and not math.isfinite(value) for value in optional_values):
            raise ValueError("native chemical attempt contains a nonfinite value")
        if any(not math.isfinite(value) or value < 0.0 for value in attempt.amounts_mol):
            raise ValueError("native chemical attempt contains invalid amounts")
        if attempt.volume_m3 is not None and attempt.volume_m3 <= 0.0:
            raise ValueError("native chemical attempt contains an invalid volume")
        if (
            attempt.kkt_dimension < 0
            or attempt.kkt_rank < 0
            or attempt.kkt_rank > attempt.kkt_dimension
            or attempt.recovery_seed_count < 0
            or attempt.recovery_seed_count < attempt.recovery_solve_count
            or attempt.recovery_solve_count < 0
        ):
            raise ValueError("native chemical attempt counters are inconsistent")
        if attempt.basin_ordinal is not None and not (
            0 <= attempt.basin_ordinal < len(search.basins)
        ):
            raise ValueError("native chemical attempt basin reference is out of range")
    if tuple(basin.ordinal for basin in search.basins) != tuple(range(len(search.basins))):
        raise ValueError("native chemical basin ordinals are inconsistent")
    for basin in search.basins:
        if not 0 <= basin.representative_attempt_ordinal < len(search.attempts):
            raise ValueError("native chemical basin representative is out of range")
        representative = search.attempts[basin.representative_attempt_ordinal]
        if (
            representative.terminal_status != "certified_local_minimum"
            or representative.basin_ordinal != basin.ordinal
            or not basin.amounts_mol
            or any(not math.isfinite(value) or value <= 0.0 for value in basin.amounts_mol)
            or not math.isfinite(basin.volume_m3)
            or basin.volume_m3 <= 0.0
            or not math.isfinite(basin.objective)
        ):
            raise ValueError("native chemical basin evidence is inconsistent")
    if (search.selected_basin_ordinal is None) != (search.selected_objective is None):
        raise ValueError("native chemical search selection is incomplete")
    if search.selected_basin_ordinal is None:
        if search.status == "certified_local_minimum":
            raise ValueError("native chemical search omitted its certified selection")
    else:
        if not 0 <= search.selected_basin_ordinal < len(search.basins):
            raise ValueError("native chemical search selection is out of range")
        selected = search.basins[search.selected_basin_ordinal]
        if (
            search.status != "certified_local_minimum"
            or search.selected_objective != selected.objective
        ):
            raise ValueError("native chemical search selection is inconsistent")
    expected_prefix_budgets = tuple(
        dict.fromkeys(
            min(requested, search.primary_attempt_count)
            for requested in (1, 5, 13, 25)
            if search.primary_attempt_count > 0
        )
    )
    if tuple(prefix.primary_budget for prefix in search.budget_prefixes) != (
        expected_prefix_budgets
    ):
        raise ValueError("native chemical search budget prefixes are incomplete")
    previous_budget = 0
    previous_selection: int | None = None
    for index, prefix in enumerate(search.budget_prefixes):
        if (
            prefix.primary_budget <= previous_budget
            or prefix.primary_budget > search.primary_attempt_count
            or prefix.attempted_primary_ordinals != tuple(range(prefix.primary_budget))
            or len(set(prefix.basin_ordinals)) != len(prefix.basin_ordinals)
            or any(not 0 <= ordinal < len(search.basins) for ordinal in prefix.basin_ordinals)
            or (
                prefix.selected_basin_ordinal is not None
                and prefix.selected_basin_ordinal not in prefix.basin_ordinals
            )
            or prefix.selection_changed
            != (index > 0 and prefix.selected_basin_ordinal != previous_selection)
        ):
            raise ValueError("native chemical search budget prefix is inconsistent")
        previous_budget = prefix.primary_budget
        previous_selection = prefix.selected_basin_ordinal


def _chemical_search(native: Mapping[str, object]) -> ChemicalEquilibriumSearch:
    payload = cast(Mapping[str, object], native["search"])
    attempts = tuple(
        ChemicalEquilibriumAttempt(
            ordinal=int(cast(int, record["ordinal"])),
            primary_ordinal=int(cast(int, record["primary_ordinal"])),
            kind=str(record["kind"]),
            parent_ordinal=(
                int(cast(int, record["parent_ordinal"]))
                if record["parent_ordinal"] is not None
                else None
            ),
            start_identity=str(record["start_identity"]),
            start_construction_status=str(record["start_construction_status"]),
            retraction_status=str(record["retraction_status"]),
            provider_domain_status=str(record["provider_domain_status"]),
            solver_status=str(record["solver_status"]),
            callback_error=str(record["callback_error"]),
            terminal_status=str(record["terminal_status"]),
            amounts_mol=tuple(float(value) for value in cast(Sequence[float], record["amounts"])),
            volume_m3=_optional_float(record, "volume_m3"),
            objective=_optional_float(record, "objective"),
            balance_inf_norm=_optional_float(record, "balance_inf_norm"),
            charge_inf_norm=_optional_float(record, "charge_inf_norm"),
            pressure_relative_residual=_optional_float(record, "pressure_relative_residual"),
            reaction_affinity_inf_norm=_optional_float(record, "reaction_affinity_inf_norm"),
            kkt_stationarity_inf_norm=_optional_float(record, "kkt_stationarity_inf_norm"),
            complementarity_inf_norm=_optional_float(record, "complementarity_inf_norm"),
            kkt_dimension=int(cast(int, record["kkt_dimension"])),
            kkt_rank=int(cast(int, record["kkt_rank"])),
            condition_number_inf=_optional_float(record, "condition_number_inf"),
            local_minimum_status=str(record["local_minimum_status"]),
            trace_status=str(record["trace_status"]),
            basin_ordinal=(
                int(cast(int, record["basin_ordinal"]))
                if record["basin_ordinal"] is not None
                else None
            ),
            recovery_seed_count=int(cast(int, record["recovery_seed_count"])),
            recovery_solve_count=int(cast(int, record["recovery_solve_count"])),
        )
        for record in cast(Sequence[Mapping[str, object]], payload["attempts"])
    )
    basins = tuple(
        ChemicalEquilibriumBasin(
            ordinal=int(cast(int, record["ordinal"])),
            representative_attempt_ordinal=int(cast(int, record["representative_attempt_ordinal"])),
            amounts_mol=tuple(float(value) for value in cast(Sequence[float], record["amounts"])),
            volume_m3=float(cast(float, record["volume_m3"])),
            objective=float(cast(float, record["objective"])),
        )
        for record in cast(Sequence[Mapping[str, object]], payload["basins"])
    )
    prefixes = tuple(
        ChemicalEquilibriumBudgetPrefix(
            primary_budget=int(cast(int, record["primary_budget"])),
            attempted_primary_ordinals=tuple(
                int(value) for value in cast(Sequence[int], record["attempted_primary_ordinals"])
            ),
            basin_ordinals=tuple(
                int(value) for value in cast(Sequence[int], record["basin_ordinals"])
            ),
            selected_basin_ordinal=(
                int(cast(int, record["selected_basin_ordinal"]))
                if record["selected_basin_ordinal"] is not None
                else None
            ),
            selection_changed=bool(record["selection_changed"]),
        )
        for record in cast(Sequence[Mapping[str, object]], payload["budget_prefixes"])
    )
    search = ChemicalEquilibriumSearch(
        status=str(payload["status"]),
        primary_budget=int(cast(int, payload["primary_budget"])),
        primary_attempt_count=int(cast(int, payload["primary_attempt_count"])),
        generated_start_count=int(cast(int, payload["generated_start_count"])),
        budget_truncated_start_count=int(cast(int, payload["budget_truncated_start_count"])),
        duplicate_start_count=int(cast(int, payload["duplicate_start_count"])),
        infeasible_start_count=int(cast(int, payload["infeasible_start_count"])),
        evaluated_start_count=int(cast(int, payload["evaluated_start_count"])),
        domain_rejected_start_count=int(cast(int, payload["domain_rejected_start_count"])),
        construction_rejected_start_count=int(
            cast(int, payload["construction_rejected_start_count"])
        ),
        attempts=attempts,
        basins=basins,
        budget_prefixes=prefixes,
        selected_basin_ordinal=(
            int(cast(int, payload["selected_basin_ordinal"]))
            if payload["selected_basin_ordinal"] is not None
            else None
        ),
        selected_objective=_optional_float(payload, "selected_objective"),
        selection_label=str(payload["selection_label"]),
    )
    _validate_chemical_search(search)
    return search


def _failed_chemical_diagnostics(
    status: str, failure_reason: str
) -> ChemicalEquilibriumDiagnostics:
    return ChemicalEquilibriumDiagnostics(
        solver_status=status,
        callback_error="",
        failure_kind="input_or_native_error",
        failure_reason=failure_reason,
        chemical_certification_level="NOT_EVALUATED",
        boundary_status="not_adjudicated",
        structural_zero_species_indices=(),
        numerical_status="not_adjudicated",
        physical_status="not_adjudicated",
        predictive_status="not_adjudicated",
        provider_domain_status="not_adjudicated",
        local_minimum_status="not_adjudicated",
        negative_curvature_recovery_status="not_needed",
        negative_curvature_recovery_attempts=0,
        negative_curvature_recovery_selected_sign=0,
        trace_status="not_adjudicated",
        globality_status="not_adjudicated",
        search=_empty_chemical_search(),
        amounts_mol=(),
        volume_m3=None,
        balance_inf_norm=None,
        charge_inf_norm=None,
        pressure_relative_residual=None,
        reaction_affinity_inf_norm=None,
        reaction_affinity_residuals=(),
        packing_fraction=None,
        packing_fraction_min=None,
        packing_fraction_max=None,
        total_ion_fraction=None,
        total_ion_fraction_max=None,
        minimum_amount_mol=None,
        trace_floor_mol=None,
        kkt_stationarity_inf_norm=None,
        physical_stationarity_residuals=(),
        physical_equality_multipliers=(),
        physical_constraint_jacobian=(),
        physical_constraint_shape=(0, 0),
        physical_lagrangian_gradient=(),
        physical_objective_gradient=(),
        physical_to_chart_jacobian=(),
        physical_to_chart_shape=(0, 0),
        chart_physical_pullback_residual_inf_norm=None,
        complementarity_inf_norm=None,
        numerical_criteria=(),
        physical_criteria=(),
        first_failed_numerical_criterion=None,
        first_failed_physical_criterion=None,
        kkt_dimension=0,
        kkt_rank=0,
        condition_number_inf=None,
        active_lower_bounds=(),
        active_upper_bounds=(),
        active_constraint_bounds=(),
        active_trace_species=(),
        chart_topology="not_adjudicated",
        reduced_hessian_status="not_adjudicated",
        reduced_hessian=(),
        reduced_hessian_nullspace_basis=(),
        reduced_hessian_nullspace_shape=(0, 0),
        reduced_hessian_eigenvalues=(),
        reduced_hessian_spectrum_status="not_adjudicated",
        reduced_hessian_raw_inertia=(0, 0, 0),
        reduced_hessian_inertia=(0, 0, 0),
        reduced_hessian_scale=None,
        reduced_hessian_eigenvalue_tolerance=None,
        objective_gradient=(),
        constraint_values=(),
        constraint_jacobian=(),
        lagrangian_gradient=(),
        equality_multipliers=(),
        chart_stationarity_inf_norm=None,
        lagrangian_hessian=(),
        covariant_lagrangian_hessian=(),
        derivative_check_step=None,
        objective_gradient_check_relative_error=None,
        physical_objective_gradient_check_relative_error=None,
        constraint_jacobian_check_relative_error=None,
        lagrangian_hessian_check_relative_error=None,
        reduced_hessian_check_relative_error=None,
        derivative_check_worst_entry="",
        derivative_check_worst_relative_error=None,
        derivative_check_worst_analytic_value=None,
        derivative_check_worst_finite_difference_value=None,
        derivative_check_worst_step=None,
        derivative_coordinate_order=(),
        derivative_objective_basis="dimensionless_fixed_TP_A_plus_PV_plus_reference_over_RT",
        derivative_constraint_basis="ordered_per_row_in_derivative_constraint_order",
        derivative_constraint_order=(),
        objective_gradient_shape=(0,),
        constraint_jacobian_shape=(0, 0),
        lagrangian_hessian_shape=(0, 0),
        kkt_root_jacobian=(),
        kkt_root_shape=(0, 0),
        kkt_root_status="not_evaluated",
        kkt_root_jacobian_check_relative_error=None,
        reference_representation_residual_inf_norm=None,
        reference_convergence_error=None,
        reference_derivative_availability=None,
    )


def _classified_source_reference_failure(
    failure_reason: str,
) -> ChemicalEquilibriumDiagnostics | None:
    classifications = (
        (
            "unsupported-source-reference-transfer",
            "unsupported_source_reference_transfer",
        ),
        ("invalid-source-composition", "invalid_source_reference"),
        ("nonfinite-source-reference-state", "invalid_source_reference"),
        ("component-order-mismatch", "component_order_mismatch"),
        ("source-reference component order", "component_order_mismatch"),
        ("missing-parameter-or-topology", "missing_parameter_or_topology"),
        ("provider-domain-rejection", "provider_domain_rejection"),
        (
            "reference-limit-convergence-failure",
            "reference_limit_convergence_failure",
        ),
        ("derivative-unavailable", "source_reference_derivative_unavailable"),
        (
            "source-reference derivative is unavailable",
            "source_reference_derivative_unavailable",
        ),
    )
    for marker, failure_kind in classifications:
        if marker in failure_reason:
            return replace(
                _failed_chemical_diagnostics("source_reference_preflight_failed", failure_reason),
                failure_kind=failure_kind,
                provider_domain_status=(
                    "failed" if failure_kind == "provider_domain_rejection" else "not_adjudicated"
                ),
            )
    if "source-reference" in failure_reason or "source " in failure_reason:
        return replace(
            _failed_chemical_diagnostics("source_reference_preflight_failed", failure_reason),
            failure_kind="source_reference_receipt_mismatch",
        )
    return None


def _chemical_diagnostics(
    native: Mapping[str, object],
) -> ChemicalEquilibriumDiagnostics:
    def criteria(name: str) -> tuple[ChemicalEquilibriumCriterion, ...]:
        def passes(record: ChemicalEquilibriumCriterion) -> bool:
            if not math.isfinite(record.value):
                return False
            if record.comparison == "<=":
                return record.value <= record.limit
            if record.comparison == ">":
                return record.value > record.limit
            if record.comparison == ">=":
                return record.value >= record.limit
            raise ValueError(f"native {name} comparison is invalid")

        records = tuple(
            ChemicalEquilibriumCriterion(
                name=str(cast(Mapping[str, object], record)["name"]),
                status=str(cast(Mapping[str, object], record)["status"]),
                value=float(cast(float, cast(Mapping[str, object], record)["value"])),
                limit=float(cast(float, cast(Mapping[str, object], record)["limit"])),
                comparison=str(cast(Mapping[str, object], record)["comparison"]),
            )
            for record in cast(Sequence[object], native[name])
        )
        if any(record.status not in {"passed", "failed"} for record in records):
            raise ValueError(f"native {name} status is invalid")
        if any(
            record.comparison not in {"<=", ">", ">="}
            or not math.isfinite(record.limit)
            or record.limit < 0.0
            or (record.status == "passed") != passes(record)
            for record in records
        ):
            raise ValueError(f"native {name} evidence is inconsistent")
        return records

    numerical_criteria = criteria("numerical_criteria")
    physical_criteria = criteria("physical_criteria")
    first_failed_numerical = next(
        (record.name for record in numerical_criteria if record.status == "failed"), None
    )
    first_failed_physical = next(
        (record.name for record in physical_criteria if record.status == "failed"), None
    )
    reduced_inertia = tuple(
        int(value) for value in cast(Sequence[int], native["reduced_hessian_inertia"])
    )
    if len(reduced_inertia) != 3 or any(value < 0 for value in reduced_inertia):
        raise ValueError("native reduced-Hessian inertia is invalid")
    reduced_raw_inertia = tuple(
        int(value) for value in cast(Sequence[int], native["reduced_hessian_raw_inertia"])
    )
    if len(reduced_raw_inertia) != 3 or any(value < 0 for value in reduced_raw_inertia):
        raise ValueError("native raw reduced-Hessian inertia is invalid")
    failure_kind = str(native["failure_kind"])
    if failure_kind not in {
        "none",
        "genuine_saddle",
        "derivative_inconsistency",
        "rank_deficiency",
        "ill_conditioning",
        "physical_domain_failure",
        "exhausted_multistart_search",
        "unsupported_derivative_capability",
    }:
        raise ValueError("native chemical failure kind is invalid")
    if bool(native["accepted"]) != (failure_kind == "none"):
        raise ValueError("native chemical acceptance and failure kind disagree")
    accepted = bool(native["accepted"])
    certification = str(native["chemical_certification_level"])
    numerical_status = str(native["numerical_status"])
    physical_status = str(native["physical_status"])
    local_minimum_status = str(native["local_minimum_status"])
    search = _chemical_search(native)
    if accepted != (certification == "LOCAL_EQUILIBRIUM") or (
        accepted
        and (
            numerical_status != "passed"
            or physical_status != "passed"
            or local_minimum_status != "passed"
            or search.status != "certified_local_minimum"
            or first_failed_numerical is not None
            or first_failed_physical is not None
        )
    ):
        raise ValueError("native chemical acceptance and certification evidence disagree")
    kkt_dimension = int(cast(int, native["kkt_dimension"]))
    kkt_rank = int(cast(int, native["kkt_rank"]))
    if kkt_dimension < 0 or not 0 <= kkt_rank <= kkt_dimension:
        raise ValueError("native chemical KKT rank evidence is invalid")
    objective_gradient = _float_tuple(native["objective_gradient"])
    reaction_affinity_residuals = _float_tuple(native["reaction_affinity_residuals"])
    physical_stationarity_residuals = _float_tuple(native["physical_stationarity_residuals"])
    physical_equality_multipliers = _float_tuple(native["physical_equality_multipliers"])
    physical_constraint_jacobian = _float_tuple(native["physical_constraint_jacobian"])
    physical_constraint_shape = tuple(
        int(value) for value in cast(Sequence[int], native["physical_constraint_shape"])
    )
    constraint_values = _float_tuple(native["constraint_values"])
    constraint_jacobian = _float_tuple(native["constraint_jacobian"])
    lagrangian_gradient = _float_tuple(native["lagrangian_gradient"])
    equality_multipliers = _float_tuple(native["equality_multipliers"])
    lagrangian_hessian = _float_tuple(native["lagrangian_hessian"])
    covariant_lagrangian_hessian = _float_tuple(native["covariant_lagrangian_hessian"])
    reduced_hessian = _float_tuple(native["reduced_hessian"])
    reduced_hessian_eigenvalues = _float_tuple(native["reduced_hessian_eigenvalues"])
    reduced_hessian_nullspace_basis = _float_tuple(native["reduced_hessian_nullspace_basis"])
    reduced_hessian_nullspace_shape = tuple(
        int(value) for value in cast(Sequence[int], native["reduced_hessian_nullspace_shape"])
    )
    spectrum_status = str(native["reduced_hessian_spectrum_status"])
    derivative_dimension = len(objective_gradient)
    derivative_coordinate_order = tuple(
        str(value) for value in cast(Sequence[str], native["derivative_coordinate_order"])
    )
    objective_gradient_shape = tuple(
        int(value) for value in cast(Sequence[int], native["objective_gradient_shape"])
    )
    constraint_jacobian_shape = tuple(
        int(value) for value in cast(Sequence[int], native["constraint_jacobian_shape"])
    )
    lagrangian_hessian_shape = tuple(
        int(value) for value in cast(Sequence[int], native["lagrangian_hessian_shape"])
    )
    derivative_constraint_order = tuple(
        str(value) for value in cast(Sequence[str], native["derivative_constraint_order"])
    )
    physical_lagrangian_gradient = _float_tuple(native["physical_lagrangian_gradient"])
    physical_objective_gradient = _float_tuple(native["physical_objective_gradient"])
    physical_to_chart_jacobian = _float_tuple(native["physical_to_chart_jacobian"])
    physical_to_chart_shape = tuple(
        int(value) for value in cast(Sequence[int], native["physical_to_chart_shape"])
    )
    kkt_root_jacobian = _float_tuple(native["kkt_root_jacobian"])
    kkt_root_shape = tuple(int(value) for value in cast(Sequence[int], native["kkt_root_shape"]))
    kkt_root_status = str(native["kkt_root_status"])
    reduced_dimension = math.isqrt(len(reduced_hessian))
    tolerance = _optional_float(native, "reduced_hessian_eigenvalue_tolerance")
    if (
        len(lagrangian_gradient) != derivative_dimension
        or len(constraint_jacobian) != len(constraint_values) * derivative_dimension
        or len(lagrangian_hessian) != derivative_dimension * derivative_dimension
        or len(covariant_lagrangian_hessian) != derivative_dimension * derivative_dimension
        or len(derivative_coordinate_order) != derivative_dimension
        or objective_gradient_shape != (derivative_dimension,)
        or constraint_jacobian_shape != (len(constraint_values), derivative_dimension)
        or lagrangian_hessian_shape != (derivative_dimension, derivative_dimension)
        or len(derivative_constraint_order) != len(constraint_values)
        or len(equality_multipliers) != len(constraint_values)
        or len(reduced_hessian) != reduced_dimension * reduced_dimension
        or len(reduced_hessian_nullspace_shape) != 2
        or len(kkt_root_shape) != 2
        or len(physical_to_chart_shape) != 2
        or len(physical_lagrangian_gradient) != physical_to_chart_shape[0]
        or len(physical_objective_gradient) != physical_to_chart_shape[0]
        or len(physical_to_chart_jacobian)
        != physical_to_chart_shape[0] * physical_to_chart_shape[1]
        or len(kkt_root_jacobian) != kkt_root_shape[0] * kkt_root_shape[1]
        or len(physical_constraint_shape) != 2
        or len(physical_equality_multipliers) != physical_constraint_shape[0]
        or len(physical_stationarity_residuals) != physical_constraint_shape[1]
        or len(physical_constraint_jacobian)
        != physical_constraint_shape[0] * physical_constraint_shape[1]
    ):
        raise ValueError("native chemical derivative evidence has inconsistent dimensions")
    if reduced_dimension:
        if (
            reduced_hessian_nullspace_shape != (reduced_dimension, derivative_dimension)
            or len(reduced_hessian_nullspace_basis) != reduced_dimension * derivative_dimension
        ):
            raise ValueError("native reduced-Hessian basis dimensions are inconsistent")
    elif (
        reduced_hessian_nullspace_shape
        not in {
            (0, 0),
            (0, derivative_dimension),
        }
        or reduced_hessian_nullspace_basis
    ):
        raise ValueError("native empty reduced-Hessian basis is inconsistent")
    if kkt_root_status not in {
        "not_evaluated",
        "interior_no_active_bounds",
        "unavailable_active_bounds",
        "evaluation_failed",
    }:
        raise ValueError("native KKT root-system status is invalid")
    if kkt_root_status == "interior_no_active_bounds":
        if (
            kkt_root_shape != (kkt_dimension, kkt_dimension)
            or len(kkt_root_jacobian) != kkt_dimension * kkt_dimension
        ):
            raise ValueError("native KKT root-system evidence is inconsistent")
    elif kkt_root_shape != (0, 0) or kkt_root_jacobian:
        raise ValueError("native unavailable KKT root-system evidence is not empty")
    if spectrum_status not in {
        "not_adjudicated",
        "not_evaluated",
        "converged",
        "not_converged",
        "nonfinite",
    }:
        raise ValueError("native reduced-Hessian spectrum status is invalid")
    if reduced_dimension and spectrum_status in {"converged", "not_converged"}:
        if len(reduced_hessian_eigenvalues) != reduced_dimension or tolerance is None:
            raise ValueError("native reduced-Hessian spectrum dimensions are invalid")
        classified = (
            sum(value > tolerance for value in reduced_hessian_eigenvalues),
            sum(abs(value) <= tolerance for value in reduced_hessian_eigenvalues),
            sum(value < -tolerance for value in reduced_hessian_eigenvalues),
        )
        if reduced_raw_inertia != classified or sum(reduced_inertia) != reduced_dimension:
            raise ValueError("native reduced-Hessian spectrum and inertia disagree")
    elif reduced_dimension and spectrum_status == "nonfinite":
        if (
            reduced_hessian_eigenvalues
            or reduced_raw_inertia != (0, 0, 0)
            or reduced_inertia != (0, 0, 0)
        ):
            raise ValueError("native nonfinite reduced-Hessian evidence is inconsistent")
    elif reduced_dimension:
        raise ValueError("native reduced-Hessian spectrum was not evaluated")
    elif reduced_raw_inertia != (0, 0, 0) or reduced_inertia != (0, 0, 0):
        raise ValueError("native empty reduced-Hessian inertia is invalid")
    reported_affinity_norm = float(cast(float, native["reaction_affinity_inf_norm"]))
    reported_stationarity_norm = float(cast(float, native["kkt_stationarity_inf_norm"]))
    if (
        reaction_affinity_residuals
        and not math.isclose(
            max(abs(value) for value in reaction_affinity_residuals),
            reported_affinity_norm,
            rel_tol=2.0e-12,
            abs_tol=2.0e-12,
        )
    ) or (
        physical_stationarity_residuals
        and not math.isclose(
            max(abs(value) for value in physical_stationarity_residuals),
            reported_stationarity_norm,
            rel_tol=2.0e-12,
            abs_tol=2.0e-12,
        )
    ):
        raise ValueError("native chemical residual vectors and norms disagree")
    sensitivities = cast(Mapping[str, object], native["sensitivities"])
    return ChemicalEquilibriumDiagnostics(
        solver_status=str(native["solver_status"]),
        callback_error=str(native["callback_error"]),
        failure_kind=failure_kind,
        failure_reason=str(native["failure_reason"]),
        chemical_certification_level=certification,
        boundary_status=str(native["boundary_status"]),
        structural_zero_species_indices=tuple(
            int(index) for index in cast(Sequence[int], native["structural_zero_species_indices"])
        ),
        numerical_status=numerical_status,
        physical_status=physical_status,
        predictive_status="not_adjudicated",
        provider_domain_status=str(native["provider_domain_status"]),
        local_minimum_status=local_minimum_status,
        negative_curvature_recovery_status=str(
            native.get("negative_curvature_recovery_status", "not_needed")
        ),
        negative_curvature_recovery_attempts=int(
            cast(int, native.get("negative_curvature_recovery_attempts", 0))
        ),
        negative_curvature_recovery_selected_sign=int(
            cast(int, native.get("negative_curvature_recovery_selected_sign", 0))
        ),
        trace_status=str(native["trace_status"]),
        globality_status=str(native["globality_certificate"]),
        search=search,
        amounts_mol=_float_tuple(native["amounts"]),
        volume_m3=_optional_float(native, "volume_m3"),
        balance_inf_norm=float(cast(float, native["balance_inf_norm"])),
        charge_inf_norm=float(cast(float, native["charge_inf_norm"])),
        pressure_relative_residual=float(cast(float, native["pressure_relative_residual"])),
        reaction_affinity_inf_norm=reported_affinity_norm,
        reaction_affinity_residuals=reaction_affinity_residuals,
        packing_fraction=float(cast(float, native["packing_fraction"])),
        packing_fraction_min=_optional_float(native, "packing_fraction_min"),
        packing_fraction_max=_optional_float(native, "packing_fraction_max"),
        total_ion_fraction=_optional_float(native, "total_ion_fraction"),
        total_ion_fraction_max=_optional_float(native, "total_ion_fraction_max"),
        minimum_amount_mol=_optional_float(native, "minimum_amount_mol"),
        trace_floor_mol=_optional_float(native, "trace_floor_mol"),
        kkt_stationarity_inf_norm=reported_stationarity_norm,
        physical_stationarity_residuals=physical_stationarity_residuals,
        physical_equality_multipliers=physical_equality_multipliers,
        physical_constraint_jacobian=physical_constraint_jacobian,
        physical_constraint_shape=physical_constraint_shape,
        physical_lagrangian_gradient=physical_lagrangian_gradient,
        physical_objective_gradient=physical_objective_gradient,
        physical_to_chart_jacobian=physical_to_chart_jacobian,
        physical_to_chart_shape=physical_to_chart_shape,
        chart_physical_pullback_residual_inf_norm=_optional_float(
            native, "chart_physical_pullback_residual_inf_norm"
        ),
        complementarity_inf_norm=float(cast(float, native["complementarity_inf_norm"])),
        numerical_criteria=numerical_criteria,
        physical_criteria=physical_criteria,
        first_failed_numerical_criterion=first_failed_numerical,
        first_failed_physical_criterion=first_failed_physical,
        kkt_dimension=kkt_dimension,
        kkt_rank=kkt_rank,
        condition_number_inf=_optional_float(native, "condition_number_inf"),
        active_lower_bounds=_index_tuple(native["active_lower_bounds"], "active lower bounds"),
        active_upper_bounds=_index_tuple(native["active_upper_bounds"], "active upper bounds"),
        active_constraint_bounds=_index_tuple(
            native["active_constraint_bounds"], "active constraint bounds"
        ),
        active_trace_species=_index_tuple(
            sensitivities["active_trace_species"], "active trace species"
        ),
        chart_topology=str(sensitivities["chart_topology"]),
        reduced_hessian_status=str(native["reduced_hessian_status"]),
        reduced_hessian=reduced_hessian,
        reduced_hessian_nullspace_basis=reduced_hessian_nullspace_basis,
        reduced_hessian_nullspace_shape=reduced_hessian_nullspace_shape,
        reduced_hessian_eigenvalues=reduced_hessian_eigenvalues,
        reduced_hessian_spectrum_status=spectrum_status,
        reduced_hessian_raw_inertia=reduced_raw_inertia,
        reduced_hessian_inertia=reduced_inertia,
        reduced_hessian_scale=_optional_float(native, "reduced_hessian_scale"),
        reduced_hessian_eigenvalue_tolerance=tolerance,
        objective_gradient=objective_gradient,
        constraint_values=constraint_values,
        constraint_jacobian=constraint_jacobian,
        lagrangian_gradient=lagrangian_gradient,
        equality_multipliers=equality_multipliers,
        chart_stationarity_inf_norm=_optional_float(native, "chart_stationarity_inf_norm"),
        lagrangian_hessian=lagrangian_hessian,
        covariant_lagrangian_hessian=covariant_lagrangian_hessian,
        derivative_check_step=_optional_float(native, "derivative_check_step"),
        objective_gradient_check_relative_error=_optional_float(
            native, "objective_gradient_check_relative_error"
        ),
        physical_objective_gradient_check_relative_error=_optional_float(
            native, "physical_objective_gradient_check_relative_error"
        ),
        constraint_jacobian_check_relative_error=_optional_float(
            native, "constraint_jacobian_check_relative_error"
        ),
        lagrangian_hessian_check_relative_error=_optional_float(
            native, "lagrangian_hessian_check_relative_error"
        ),
        reduced_hessian_check_relative_error=_optional_float(
            native, "reduced_hessian_check_relative_error"
        ),
        derivative_check_worst_entry=str(native["derivative_check_worst_entry"]),
        derivative_check_worst_relative_error=_optional_float(
            native, "derivative_check_worst_relative_error"
        ),
        derivative_check_worst_analytic_value=_optional_float(
            native, "derivative_check_worst_analytic_value"
        ),
        derivative_check_worst_finite_difference_value=_optional_float(
            native, "derivative_check_worst_finite_difference_value"
        ),
        derivative_check_worst_step=_optional_float(native, "derivative_check_worst_step"),
        derivative_coordinate_order=derivative_coordinate_order,
        derivative_objective_basis=str(native["derivative_objective_basis"]),
        derivative_constraint_basis=str(native["derivative_constraint_basis"]),
        derivative_constraint_order=derivative_constraint_order,
        objective_gradient_shape=objective_gradient_shape,
        constraint_jacobian_shape=constraint_jacobian_shape,
        lagrangian_hessian_shape=lagrangian_hessian_shape,
        kkt_root_jacobian=kkt_root_jacobian,
        kkt_root_shape=kkt_root_shape,
        kkt_root_status=kkt_root_status,
        kkt_root_jacobian_check_relative_error=_optional_float(
            native, "kkt_root_jacobian_check_relative_error"
        ),
        reference_representation_residual_inf_norm=_optional_float(
            native, "reference_representation_residual_inf_norm"
        )
        if "reference_representation_residual_inf_norm" in native
        else None,
        reference_convergence_error=_optional_float(native, "reference_convergence_error")
        if "reference_convergence_error" in native
        else None,
        reference_derivative_availability=int(
            cast(int, native["reference_derivative_availability"])
        )
        if "reference_derivative_availability" in native
        else None,
    )


@cache
def _distribution_identity(distribution_name: str) -> tuple[str, str]:
    distribution = metadata.distribution(distribution_name)
    record = distribution.read_text("RECORD")
    if record is None or not record.strip():
        raise ValueError(f"installed {distribution_name} distribution has no RECORD identity")
    digest = hashlib.sha256(record.encode("utf-8")).hexdigest()
    return distribution.version, f"sha256:{digest}"


def _chemical_artifact_identity(
    native: Mapping[str, object],
    *,
    uses_provider: bool,
) -> ChemicalArtifactIdentity:
    engine_version, engine_record = _distribution_identity("epcsaft")
    if not uses_provider:
        return ChemicalArtifactIdentity(
            equilibrium_distribution="epcsaft",
            equilibrium_version=engine_version,
            equilibrium_record_sha256=engine_record,
            provider_distribution=None,
            provider_version=None,
            provider_record_sha256=None,
            provider_transport_capsule_name=None,
        )
    capsule_name = str(native["provider_transport_capsule_name"])
    if capsule_name != "epcsaft.native_transport":
        raise ValueError("native result has an incompatible Provider transport identity")
    return ChemicalArtifactIdentity(
        equilibrium_distribution="epcsaft",
        equilibrium_version=engine_version,
        equilibrium_record_sha256=engine_record,
        provider_distribution="epcsaft",
        provider_version=engine_version,
        provider_record_sha256=engine_record,
        provider_transport_capsule_name=capsule_name,
    )


def _active_parameter_name(parameter: epcsaft.ActiveParameter) -> str:
    return f"provider_parameter[{parameter.identity}]"


def _sensitivity_parameter(
    name: str,
    active_parameters: tuple[epcsaft.ActiveParameter, ...] = (),
) -> ChemicalEquilibriumSensitivityParameter:
    if name == "pressure_pa":
        return ChemicalEquilibriumSensitivityParameter(
            name=name,
            kind="pressure",
            source_index=None,
            input_unit="Pa",
            amount_derivative_unit="mol/Pa",
            volume_derivative_unit="m^3/Pa",
        )
    for prefix, kind, input_unit, amount_unit, volume_unit in (
        (
            "balance_total[",
            "compiled_balance_total",
            "mol",
            "mol/mol",
            "m^3/mol",
        ),
        (
            "ln_k_provider_basis[",
            "provider_basis_ln_k",
            "dimensionless",
            "mol",
            "m^3",
        ),
    ):
        if name.startswith(prefix) and name.endswith("]"):
            index_text = name[len(prefix) : -1]
            if index_text.isdecimal():
                return ChemicalEquilibriumSensitivityParameter(
                    name=name,
                    kind=kind,
                    source_index=int(index_text),
                    input_unit=input_unit,
                    amount_derivative_unit=amount_unit,
                    volume_derivative_unit=volume_unit,
                )
    for parameter in active_parameters:
        if name == _active_parameter_name(parameter):
            applicability = parameter.applicability
            component_ids = () if applicability is None else applicability.component_ids
            return ChemicalEquilibriumSensitivityParameter(
                name=name,
                kind=f"provider_{parameter.family}",
                source_index=None,
                input_unit=parameter.unit,
                amount_derivative_unit=f"mol/{parameter.unit}",
                volume_derivative_unit=f"m^3/{parameter.unit}",
                identity=parameter.identity,
                component_ids=component_ids,
                value=parameter.value,
            )
    raise ValueError(f"native sensitivity parameter identity is unsupported: {name}")


def _index_tuple(payload: object, name: str) -> tuple[int, ...]:
    values = tuple(int(value) for value in cast(Sequence[int], payload))
    if any(value < 0 for value in values):
        raise ValueError(f"native sensitivity {name} contains a negative index")
    return values


def _chemical_sensitivity(
    native: Mapping[str, object],
    *,
    species_ids: tuple[str, ...],
    model_fingerprint: str,
    active_parameters: tuple[epcsaft.ActiveParameter, ...],
) -> ChemicalEquilibriumSensitivity:
    payload = cast(Mapping[str, object], native["sensitivities"])
    status = str(payload["status"])
    if status not in {"available", "unavailable"}:
        raise ValueError("native sensitivity status is invalid")
    failure_reason = str(payload["failure_reason"])
    parameter_names = tuple(
        str(value) for value in cast(Sequence[object], payload["parameter_order"])
    )
    parameters = tuple(_sensitivity_parameter(name, active_parameters) for name in parameter_names)
    amount_derivatives = tuple(
        _vector(row, len(species_ids), "chemical sensitivity amount row")
        for row in cast(Sequence[object], payload["amount_derivatives"])
    )
    volume_derivatives = _float_tuple(payload["volume_derivatives"])
    if len(amount_derivatives) != len(parameters) or len(volume_derivatives) != len(parameters):
        raise ValueError("native sensitivity parameter and derivative dimensions disagree")
    if status == "available":
        if failure_reason or not parameters:
            raise ValueError("available native sensitivity has incomplete status evidence")
        balance_parameters = tuple(
            parameter for parameter in parameters if parameter.kind == "compiled_balance_total"
        )
        reaction_parameters = tuple(
            parameter for parameter in parameters if parameter.kind == "provider_basis_ln_k"
        )
        expected_parameters = (
            balance_parameters
            + reaction_parameters
            + (_sensitivity_parameter("pressure_pa"),)
            + tuple(
                _sensitivity_parameter(_active_parameter_name(parameter), active_parameters)
                for parameter in active_parameters
            )
        )
        if (
            parameters != expected_parameters
            or tuple(parameter.source_index for parameter in balance_parameters)
            != tuple(range(len(balance_parameters)))
            or tuple(parameter.source_index for parameter in reaction_parameters)
            != tuple(range(len(reaction_parameters)))
        ):
            raise ValueError("available native sensitivity parameter order is incomplete")
        values = (value for row in amount_derivatives for value in row)
        if not all(math.isfinite(value) for value in (*values, *volume_derivatives)):
            raise ValueError("available native sensitivity contains a non-finite derivative")
    elif parameters or amount_derivatives or volume_derivatives or not failure_reason:
        raise ValueError("unavailable native sensitivity exposes derivative columns")
    kkt_dimension = int(cast(int, payload["kkt_dimension"]))
    kkt_rank = int(cast(int, payload["kkt_rank"]))
    condition_number = float(cast(float, payload["condition_number_inf"]))
    if kkt_dimension < 0 or not 0 <= kkt_rank <= kkt_dimension:
        raise ValueError("native sensitivity KKT rank evidence is invalid")
    if math.isnan(condition_number) or condition_number < 0.0:
        raise ValueError("native sensitivity KKT conditioning evidence is invalid")
    if status == "available" and (
        kkt_dimension == 0 or kkt_rank != kkt_dimension or not math.isfinite(condition_number)
    ):
        raise ValueError("available native sensitivity lacks a conditioned square KKT system")
    if status == "available" and len(parameters) - len(active_parameters) != kkt_dimension:
        raise ValueError("available native sensitivity omits a KKT input column")
    chart_topology = str(payload["chart_topology"])
    parameter_fingerprint = str(payload["parameter_fingerprint"])
    if not chart_topology or parameter_fingerprint != model_fingerprint:
        raise ValueError("native sensitivity topology or parameter fingerprint is invalid")
    active_trace_species = _index_tuple(payload["active_trace_species"], "active trace species")
    if any(index >= len(species_ids) for index in active_trace_species):
        raise ValueError("native sensitivity trace-species index is out of range")
    provider_parameter_status = str(payload["provider_parameter_status"])
    reference_parameter_status = str(payload["reference_parameter_status"])
    if provider_parameter_status not in {"available", "not_applicable", "unavailable"} or (
        reference_parameter_status not in {"available", "not_applicable", "unavailable"}
    ):
        raise ValueError("native sensitivity parameter-support status is invalid")
    provider_parameter_failure = str(payload["provider_parameter_failure_reason"])
    reference_parameter_failure = str(payload["reference_parameter_failure_reason"])
    if (provider_parameter_status == "unavailable") != bool(provider_parameter_failure) or (
        reference_parameter_status == "unavailable"
    ) != bool(reference_parameter_failure):
        raise ValueError("native sensitivity parameter-support evidence is inconsistent")
    active_lower_bounds = _index_tuple(payload["active_lower_bounds"], "active lower bounds")
    active_upper_bounds = _index_tuple(payload["active_upper_bounds"], "active upper bounds")
    active_constraint_bounds = _index_tuple(
        payload["active_constraint_bounds"], "active constraint bounds"
    )
    if status == "available" and (
        active_lower_bounds
        or active_upper_bounds
        or active_constraint_bounds
        or active_trace_species
    ):
        raise ValueError("available native sensitivity has active-set evidence")
    return ChemicalEquilibriumSensitivity(
        status=status,
        failure_reason=failure_reason,
        parameters=parameters,
        amount_state_order=species_ids,
        amount_derivatives=amount_derivatives,
        volume_derivatives=volume_derivatives,
        kkt_dimension=kkt_dimension,
        kkt_rank=kkt_rank,
        condition_number_inf=condition_number,
        active_lower_bounds=active_lower_bounds,
        active_upper_bounds=active_upper_bounds,
        active_constraint_bounds=active_constraint_bounds,
        active_trace_species=active_trace_species,
        chart_topology=chart_topology,
        parameter_fingerprint=parameter_fingerprint,
        provider_parameter_status=provider_parameter_status,
        provider_parameter_failure_reason=provider_parameter_failure,
        reference_parameter_status=reference_parameter_status,
        reference_parameter_failure_reason=reference_parameter_failure,
    )


def _source_reference_transfer(
    phase: ProviderPhase,
    standard_state: ChemicalStandardState,
    *,
    species_ids: tuple[str, ...],
    temperature_k: float,
    pressure_pa: float,
    sensitivities_requested: bool,
    active_parameters_requested: bool,
) -> epcsaft.SourceReferenceTransfer:
    if standard_state.source_reference_component_ids != species_ids:
        raise ValueError("source-reference component order does not match the reaction system")
    if tuple(phase.model.component_ids) != species_ids:
        raise ValueError("source-reference component order does not match the installed Provider")
    if len(standard_state.log_activity_scale_factors) != len(species_ids):
        raise ValueError("source activity-scale factors do not match the species order")
    required_derivatives = (
        ()
        if not sensitivities_requested
        else (("pressure", "parameters") if active_parameters_requested else ("pressure",))
    )
    transfer = epcsaft.source_reference_transfer(
        phase.model,
        epcsaft.SourceReferenceDeclaration(
            reference_state_id=standard_state.id,
            component_ids=standard_state.source_reference_component_ids,
            solvent_mole_fractions=(standard_state.source_reference_solvent_composition),
            ion_pairs=standard_state.source_reference_ion_pairs,
            phase=standard_state.source_reference_phase,
            reference_convention=standard_state.source_reference_convention,
            activity_convention_id=(standard_state.source_reference_activity_convention_id),
            standard_molality_mol_per_kg=(
                standard_state.source_reference_standard_molality_mol_per_kg
            ),
            reference_pressure_pa=standard_state.reference_pressure_pa,
            required_derivatives=required_derivatives,
        ),
        T=temperature_k * epcsaft.unit_registry.kelvin,
        P=pressure_pa * epcsaft.unit_registry.pascal,
    )
    if (
        transfer.status != "ok"
        or transfer.domain_status != "admitted-domain"
        or transfer.convergence_status != "converged-limit"
        or transfer.component_ids != species_ids
        or transfer.reference_state_id != standard_state.id
        or transfer.activity_convention_id != standard_state.source_reference_activity_convention_id
        or transfer.parameter_fingerprint != phase.model.parameter_fingerprint
        or transfer.temperature_k != temperature_k
        or transfer.pressure_pa != pressure_pa
        or transfer.source_reference_pressure_pa != standard_state.reference_pressure_pa
    ):
        raise ValueError("installed Provider source-reference receipt is incompatible")
    return transfer


def _source_reference_payload(
    transfer: epcsaft.SourceReferenceTransfer,
) -> dict[str, object]:
    return {
        "status": transfer.status,
        "domain_status": transfer.domain_status,
        "convergence_status": transfer.convergence_status,
        "reference_state_id": transfer.reference_state_id,
        "activity_convention_id": transfer.activity_convention_id,
        "component_ids": transfer.component_ids,
        "neutral_basis_row_count": len(transfer.neutral_basis),
        "neutral_basis": tuple(value for row in transfer.neutral_basis for value in row),
        "log_fugacity_contractions": transfer.log_fugacity_contractions,
        "activity_scale_log_contractions": (transfer.activity_scale_log_contractions),
        "transfer_log_contractions": transfer.transfer_log_contractions,
        "source_composition": transfer.source_composition,
        "derivative_availability": transfer.derivative_availability,
        "temperature_k": transfer.temperature_k,
        "pressure_pa": transfer.pressure_pa,
        "source_reference_pressure_pa": transfer.source_reference_pressure_pa,
        "standard_molality_mol_per_kg": (transfer.standard_molality_mol_per_kg),
        "reference_convergence_error": transfer.reference_convergence_error,
        "source_temperature_min_k": transfer.source_temperature_interval_k[0],
        "source_temperature_max_k": transfer.source_temperature_interval_k[1],
        "source_pressure_min_pa": transfer.source_pressure_interval_pa[0],
        "source_pressure_max_pa": transfer.source_pressure_interval_pa[1],
        "parameter_fingerprint": transfer.parameter_fingerprint,
        "topology_fingerprint": transfer.topology_fingerprint,
        "component_order_fingerprint": transfer.component_order_fingerprint,
        "reference_state_fingerprint": transfer.reference_state_fingerprint,
        "domain_fingerprint": transfer.domain_fingerprint,
        "artifact_fingerprint": transfer.artifact_fingerprint,
        "helmholtz_basis_id": transfer.helmholtz_basis_id,
    }


def chemical_equilibrium(
    phase: IdealGasPhase | ProviderPhase,
    temperature: Quantity[Any],
    pressure: Quantity[Any],
    problem: ChemicalEquilibriumProblem,
    *,
    sensitivity_request: ChemicalEquilibriumSensitivityRequest | None = None,
) -> ChemicalEquilibriumResult:
    """Solve and certify one local fixed-T,P homogeneous reacting phase."""

    try:
        if not isinstance(problem, ChemicalEquilibriumProblem):
            raise TypeError("chemical_equilibrium requires a typed problem")
        if sensitivity_request is not None and not isinstance(
            sensitivity_request, ChemicalEquilibriumSensitivityRequest
        ):
            raise TypeError("sensitivity_request must be a typed sensitivity request")
        active_parameter_set = (
            None if sensitivity_request is None else sensitivity_request.active_parameters
        )
        if active_parameter_set is not None and not isinstance(
            active_parameter_set, epcsaft.ActiveParameterSet
        ):
            raise TypeError("active_parameters must be a canonical ActiveParameterSet")
        active_parameters = () if active_parameter_set is None else active_parameter_set.coordinates
        temperature_k = _quantity(temperature, "kelvin", "temperature", "chemical_equilibrium")
        pressure_pa = _quantity(pressure, "pascal", "pressure", "chemical_equilibrium")
        standard_state = problem.source_standard_state
        if standard_state is not None:
            if (
                not math.isfinite(standard_state.reference_pressure_pa)
                or standard_state.reference_pressure_pa <= 0.0
            ):
                raise ValueError(
                    "source-standard-state reference pressure must be positive and finite"
                )
            if (
                not standard_state.id
                or not standard_state.activity_scale_id
                or not standard_state.source_reference_activity_convention_id
                or not math.isfinite(standard_state.source_reference_standard_molality_mol_per_kg)
                or standard_state.source_reference_standard_molality_mol_per_kg <= 0.0
            ):
                raise ValueError("source-standard-state declaration is incomplete")
        if not math.isfinite(problem.strict_interior_amount_floor_mol) or (
            problem.strict_interior_amount_floor_mol <= 0.0
        ):
            raise ValueError("minimum admitted amount must be positive and finite")
        records = tuple(
            {
                "source_id": record.source_id,
                "reference_id": record.reference_id,
                "reaction_orientation": record.reaction_orientation,
                "conversion_id": record.conversion_id,
                "dimensionless": record.dimensionless,
                "temperature_k": temperature_k,
                "pressure_pa": (
                    standard_state.reference_pressure_pa
                    if standard_state is not None
                    else pressure_pa
                ),
            }
            for record in problem.equilibrium_constants
        )
        spec = {
            "species_ids": problem.species_ids,
            "charges": problem.charges,
            "molar_masses_kg_per_mol": problem.molar_masses_kg_per_mol,
            "balance_matrix": problem.balance_matrix,
            "conserved_totals": problem.conserved_totals,
            "reaction_matrix": problem.reaction_matrix,
            "feed_amounts": problem.feed_amounts_mol,
            "ln_k": tuple(record.ln_value for record in problem.equilibrium_constants),
            "equilibrium_constant_records": records,
            "temperature_k": temperature_k,
            "pressure_pa": pressure_pa,
        }
        packing_bounds: tuple[float, float] | None
        source_transfer: epcsaft.SourceReferenceTransfer | None = None
        if isinstance(phase, IdealGasPhase):
            if active_parameter_set is not None:
                raise ValueError("ideal-gas sensitivity cannot request Provider coordinates")
            if problem.source_standard_state is not None:
                raise ValueError("ideal-gas problems cannot request Provider reference transport")
            if not phase.model_fingerprint or not phase.reference_id:
                raise ValueError("ideal-gas model identity is incomplete")
            if any(
                record.reference_id != phase.reference_id
                for record in problem.equilibrium_constants
            ):
                raise ValueError("reaction constants do not match the ideal-gas reference")
            spec["provider_fingerprint"] = phase.model_fingerprint
            capsule = None
            packing_bounds = None
            thermodynamic_model = "ideal_gas"
            provider_ids = None
            provider_fingerprint = None
            model_fingerprint = phase.model_fingerprint
        elif isinstance(phase, ProviderPhase):
            if not isinstance(phase.model, epcsaft.Mixture):
                raise TypeError("ProviderPhase requires an epcsaft.Mixture")
            packing_bounds = (
                float(phase.admissible_packing_fraction_interval[0]),
                float(phase.admissible_packing_fraction_interval[1]),
            )
            if (
                not all(math.isfinite(value) for value in packing_bounds)
                or not 0.0 < packing_bounds[0] < packing_bounds[1]
            ):
                raise ValueError("packing-fraction bounds must be finite and increasing")
            model_fingerprint = phase.model.parameter_fingerprint
            capsule = _native_transport(phase.model, active_parameter_set)
            thermodynamic_model = "engine_model"
            provider_ids = tuple(phase.model.component_ids)
            provider_fingerprint = model_fingerprint
            if standard_state is not None:
                source_transfer = _source_reference_transfer(
                    phase,
                    standard_state,
                    species_ids=problem.species_ids,
                    temperature_k=temperature_k,
                    pressure_pa=pressure_pa,
                    sensitivities_requested=sensitivity_request is not None,
                    active_parameters_requested=bool(active_parameters),
                )
        else:
            raise TypeError("chemical_equilibrium requires a typed phase model")
        if standard_state is not None and source_transfer is None:
            raise RuntimeError("source-reference transfer was not evaluated")
        native_standard_state: dict[str, object] | None = (
            None
            if standard_state is None
            else {
                "id": standard_state.id,
                "activity_scale_id": standard_state.activity_scale_id,
                "log_activity_scale_factors": standard_state.log_activity_scale_factors,
                "reference_pressure_pa": standard_state.reference_pressure_pa,
                "source_reference_activity_convention_id": (
                    standard_state.source_reference_activity_convention_id
                ),
                "source_reference_standard_molality_mol_per_kg": (
                    standard_state.source_reference_standard_molality_mol_per_kg
                ),
                "provider_transfer": _source_reference_payload(
                    cast(epcsaft.SourceReferenceTransfer, source_transfer)
                ),
            }
        )
        native = cast(
            Mapping[str, object],
            _equilibrium._chemical_equilibrium(
                capsule,
                spec,
                native_standard_state,
                packing_bounds,
                problem.strict_interior_amount_floor_mol,
                active_parameter_set is not None,
            ),
        )
        diagnostics = _chemical_diagnostics(native)
    except (
        AttributeError,
        IndexError,
        KeyError,
        RuntimeError,
        TypeError,
        ValueError,
    ) as error:
        reason = str(error)
        classified = _classified_source_reference_failure(reason)
        diagnostics = classified or _failed_chemical_diagnostics("input_or_native_error", reason)
        raise ChemicalEquilibriumError(
            reason,
            diagnostics,
        ) from error

    try:
        if not bool(native["accepted"]):
            reason = diagnostics.failure_reason or str(native["callback_error"])
            raise ChemicalEquilibriumError(reason, diagnostics)
        if provider_fingerprint is not None and (
            str(native["parameter_fingerprint"]) != provider_fingerprint
        ):
            raise ValueError("native result has the wrong Provider fingerprint")
        if source_transfer is not None and (
            str(native.get("source_reference_state_fingerprint", ""))
            != source_transfer.reference_state_fingerprint
            or str(native.get("provider_reference_id", "")) != source_transfer.helmholtz_basis_id
        ):
            raise ValueError("native result has the wrong source-reference receipt")
        amounts = _vector(native["amounts"], len(problem.species_ids), "chemical amounts")
        total_amount = math.fsum(amounts)
        if total_amount <= 0.0 or not math.isfinite(total_amount):
            raise ValueError("native chemical amounts are invalid")
        volume_m3 = float(cast(float, native["volume_m3"]))
        provider_reference_id = (
            str(native["provider_reference_id"]) if "provider_reference_id" in native else None
        )
        standard_offsets = (
            _float_tuple(native["standard_offsets"]) if "standard_offsets" in native else None
        )
        ln_k_provider_basis = (
            _float_tuple(native["ln_k_provider_basis"]) if "ln_k_provider_basis" in native else None
        )
        artifact_identity = _chemical_artifact_identity(
            native,
            uses_provider=provider_fingerprint is not None,
        )
        sensitivity = (
            None
            if sensitivity_request is None
            else _chemical_sensitivity(
                native,
                species_ids=problem.species_ids,
                model_fingerprint=model_fingerprint,
                active_parameters=active_parameters,
            )
        )
    except ChemicalEquilibriumError:
        raise
    except (KeyError, metadata.PackageNotFoundError, RuntimeError, TypeError, ValueError) as error:
        raise ChemicalEquilibriumError(
            str(error), _failed_chemical_diagnostics("payload_error", str(error))
        ) from error
    return ChemicalEquilibriumResult(
        temperature_k=temperature_k,
        pressure_pa=pressure_pa,
        species_ids=problem.species_ids,
        charges=problem.charges,
        amounts_mol=amounts,
        mole_fractions=tuple(amount / total_amount for amount in amounts),
        volume_m3=volume_m3,
        thermodynamic_model=thermodynamic_model,
        model_fingerprint=model_fingerprint,
        provider_component_ids=provider_ids,
        provider_parameter_fingerprint=provider_fingerprint,
        equilibrium_constants=problem.equilibrium_constants,
        source_standard_state=standard_state,
        source_reference_transfer=source_transfer,
        provider_reference_id=provider_reference_id,
        standard_offsets=standard_offsets,
        ln_k_provider_basis=ln_k_provider_basis,
        local_scope="fixed_TP_single_homogeneous_phase",
        diagnostics=diagnostics,
        response_kind=(
            "value_only"
            if sensitivity is None
            else (
                "value_plus_jacobian"
                if sensitivity.status == "available"
                else "value_with_unavailable_jacobian"
            )
        ),
        artifact_identity=artifact_identity,
        sensitivity=sensitivity,
    )
