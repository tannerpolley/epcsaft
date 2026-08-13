from __future__ import annotations

import hashlib
import json
import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, replace
from typing import Any, Literal, cast

from pint import Quantity

import epcsaft

from .._trials import NonEvaluableTrial, NonEvaluableTrialCode
from ..eos import _native_topology_fingerprint
from ._api import (
    AllComponents,
    ChemicalEquilibriumDiagnostics,
    ChemicalEquilibriumError,
    ChemicalEquilibriumProblem,
    ChemicalEquilibriumResult,
    ChemicalEquilibriumSearch,
    ChemicalEquilibriumSensitivityRequest,
    DeclaredTopology,
    FixedTPFeed,
    IdealGasPhase,
    PhaseEquilibriumError,
    PhaseEquilibriumProblem,
    PhaseSpecification,
    ProviderPhase,
    _chemical_equilibrium_impl,
    chemical_equilibrium,
    phase_equilibrium,
)

type HomogeneousReactiveQuantity = Literal[
    "amount_mol",
    "mole_fraction",
    "molarity_mol_per_m3",
    "molality_mol_per_kg",
    "molar_density_mol_per_m3",
    "compressibility_factor",
    "pressure_pa",
    "fugacity_pa",
    "log_activity_contraction",
    "log_fugacity_contraction",
    "chemical_potential_over_rt_contraction",
]
type ObservationSupport = Literal[
    "real",
    "positive",
    "nonnegative",
    "left_censored",
    "right_censored",
]
type HomogeneousPhaseRole = Literal["liquid", "vapor"]
type HomogeneousBranchPolicy = Literal["local_certified_role_selected_state"]


_QUANTITY_UNITS: dict[HomogeneousReactiveQuantity, str] = {
    "amount_mol": "mole",
    "mole_fraction": "dimensionless",
    "molarity_mol_per_m3": "mole / meter**3",
    "molality_mol_per_kg": "mole / kilogram",
    "molar_density_mol_per_m3": "mole / meter**3",
    "compressibility_factor": "dimensionless",
    "pressure_pa": "pascal",
    "fugacity_pa": "pascal",
    "log_activity_contraction": "dimensionless",
    "log_fugacity_contraction": "dimensionless",
    "chemical_potential_over_rt_contraction": "dimensionless",
}
_COEFFICIENT_QUANTITIES = {
    "amount_mol",
    "mole_fraction",
    "molarity_mol_per_m3",
    "molality_mol_per_kg",
    "log_activity_contraction",
    "log_fugacity_contraction",
    "chemical_potential_over_rt_contraction",
    "fugacity_pa",
}
_POTENTIAL_QUANTITIES = {
    "log_activity_contraction",
    "log_fugacity_contraction",
    "chemical_potential_over_rt_contraction",
    "fugacity_pa",
}


def _canonical_fingerprint(schema: str, value: object) -> str:
    payload = json.dumps(
        {"schema": schema, "value": value},
        allow_nan=False,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("utf-8")
    return f"sha256:{hashlib.sha256(payload).hexdigest()}"


def _canonical_json_bytes(value: object) -> bytes:
    return json.dumps(
        value,
        ensure_ascii=False,
        allow_nan=False,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("utf-8")


@dataclass(frozen=True, slots=True)
class HomogeneousReactiveObservationRow:
    """One generic solved-phase quantity and its measurement support contract."""

    identity: str
    quantity: HomogeneousReactiveQuantity
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

    @property
    def unit(self) -> str:
        return _QUANTITY_UNITS[self.quantity]


@dataclass(frozen=True, slots=True)
class HomogeneousContinuationReference:
    """Immutable, role-certified local branch anchor for pure trial evaluation."""

    mole_fractions: tuple[float, ...]
    molar_volume_m3_per_mol: float
    maximum_log_composition_distance: float
    maximum_log_volume_distance: float
    phase_identity: str
    phase_role: HomogeneousPhaseRole
    anchor_state_identity: str
    anchor_phase_role_receipt_fingerprint: str
    anchor_parameter_fingerprint: str
    anchor_reaction_fingerprint: str
    anchor_topology_fingerprint: str
    anchor_basin_states: tuple[tuple[tuple[float, ...], float], ...]

    def __post_init__(self) -> None:
        object.__setattr__(self, "mole_fractions", tuple(self.mole_fractions))
        object.__setattr__(
            self,
            "anchor_basin_states",
            tuple(
                (tuple(composition), float(molar_volume))
                for composition, molar_volume in self.anchor_basin_states
            ),
        )
        matching_anchor_states = 0
        for composition, molar_volume in self.anchor_basin_states:
            if (
                len(composition) != len(self.mole_fractions)
                or any(not math.isfinite(value) or value <= 0.0 for value in composition)
                or not math.isclose(
                    math.fsum(composition), 1.0, rel_tol=0.0, abs_tol=1.0e-12
                )
                or not math.isfinite(molar_volume)
                or molar_volume <= 0.0
            ):
                raise ValueError("continuation basin catalog is invalid")
            composition_distance = max(
                abs(math.log(actual / anchor))
                for actual, anchor in zip(
                    composition, self.mole_fractions, strict=True
                )
            )
            volume_distance = abs(
                math.log(molar_volume / self.molar_volume_m3_per_mol)
            )
            matching_anchor_states += int(
                composition_distance <= 2.0e-10
                and volume_distance <= 2.0e-10
            )
        if (
            not self.mole_fractions
            or any(not math.isfinite(value) or value <= 0.0 for value in self.mole_fractions)
            or not math.isclose(math.fsum(self.mole_fractions), 1.0, rel_tol=0.0, abs_tol=1.0e-12)
            or not math.isfinite(self.molar_volume_m3_per_mol)
            or self.molar_volume_m3_per_mol <= 0.0
            or not math.isfinite(self.maximum_log_composition_distance)
            or self.maximum_log_composition_distance <= 0.0
            or not math.isfinite(self.maximum_log_volume_distance)
            or self.maximum_log_volume_distance <= 0.0
            or not self.phase_identity
            or self.phase_role not in {"liquid", "vapor"}
            or not self.anchor_state_identity.startswith("sha256:")
            or not self.anchor_phase_role_receipt_fingerprint.startswith("sha256:")
            or not self.anchor_parameter_fingerprint.startswith("sha256:")
            or not self.anchor_reaction_fingerprint.startswith("sha256:")
            or not self.anchor_topology_fingerprint.startswith("sha256:")
            or not self.anchor_basin_states
            or matching_anchor_states != 1
        ):
            raise ValueError(
                "continuation reference requires a positive, role-certified local state"
            )

    def to_mapping(self) -> dict[str, object]:
        return {
            "mole_fractions": self.mole_fractions,
            "molar_volume_m3_per_mol": self.molar_volume_m3_per_mol,
            "maximum_log_composition_distance": self.maximum_log_composition_distance,
            "maximum_log_volume_distance": self.maximum_log_volume_distance,
            "phase_identity": self.phase_identity,
            "phase_role": self.phase_role,
            "anchor_state_identity": self.anchor_state_identity,
            "anchor_phase_role_receipt_fingerprint": (
                self.anchor_phase_role_receipt_fingerprint
            ),
            "anchor_parameter_fingerprint": self.anchor_parameter_fingerprint,
            "anchor_reaction_fingerprint": self.anchor_reaction_fingerprint,
            "anchor_topology_fingerprint": self.anchor_topology_fingerprint,
            "anchor_basin_states": tuple(
                {
                    "mole_fractions": composition,
                    "molar_volume_m3_per_mol": molar_volume,
                }
                for composition, molar_volume in self.anchor_basin_states
            ),
        }

    @property
    def fingerprint(self) -> str:
        return _canonical_fingerprint(
            "epcsaft.homogeneous-continuation-reference",
            self.to_mapping(),
        )


@dataclass(frozen=True, slots=True)
class HomogeneousReactiveObservationProblem:
    """One explicitly identified homogeneous phase and reacting system."""

    identity: str
    phase_identity: str
    phase_role: HomogeneousPhaseRole
    phase: IdealGasPhase | ProviderPhase
    reaction_system: ChemicalEquilibriumProblem
    continuation_identity: str
    branch_policy: HomogeneousBranchPolicy
    continuation_reference: HomogeneousContinuationReference | None = None

    def __post_init__(self) -> None:
        required = (self.identity, self.phase_identity, self.continuation_identity)
        if any(not value for value in required):
            raise ValueError("homogeneous problem and phase identities must be nonempty")
        if self.phase_role not in {"liquid", "vapor"}:
            raise ValueError("homogeneous phase role must be liquid or vapor")
        if self.branch_policy != "local_certified_role_selected_state":
            raise ValueError(
                "CAP-12 homogeneous observations require a locally certified role-selected state"
            )
        if isinstance(self.phase, IdealGasPhase) and self.phase_role != "vapor":
            raise ValueError("an ideal-gas homogeneous phase must have the vapor role")
        if not isinstance(self.reaction_system, ChemicalEquilibriumProblem):
            raise TypeError("homogeneous reaction_system must be a ChemicalEquilibriumProblem")
        if self.continuation_reference is not None and len(
            self.continuation_reference.mole_fractions
        ) != len(self.reaction_system.species_ids):
            raise ValueError("continuation reference does not match the component order")


@dataclass(frozen=True, slots=True)
class HomogeneousReactiveObservationDescriptor:
    """Closed Equilibrium-owned descriptor for one fixed-T,P reacting phase."""

    identity: Literal["equilibrium.homogeneous-reactive"]
    version: int
    component_ids: tuple[str, ...]
    charges: tuple[int, ...]
    problem_identity: str
    phase_identity: str
    phase_role: HomogeneousPhaseRole
    continuation_identity: str
    branch_policy: HomogeneousBranchPolicy
    rows: tuple[HomogeneousReactiveObservationRow, ...]
    row_identities: tuple[str, ...]
    row_units: tuple[str, ...]
    active_parameter_identities: tuple[str, ...]
    active_parameter_units: tuple[str, ...]
    parameter_fingerprint: str
    active_parameter_selection_fingerprint: str
    topology_fingerprint: str
    reference_fingerprint: str
    domain_fingerprint: str
    problem_fingerprint: str
    reaction_fingerprint: str
    local_scope: Literal["fixed_TP_single_homogeneous_phase"]
    jacobian_shape: tuple[int, int]
    fingerprint: str


@dataclass(frozen=True, slots=True)
class HomogeneousReactiveObservationRowResult:
    identity: str
    status: Literal["evaluated", "non_evaluable"]
    value: float | None
    jacobian: tuple[float, ...] | None
    reason: str


@dataclass(frozen=True, slots=True)
class ActiveParameterTrialReceipt:
    parameter_fingerprint: str
    identities: tuple[str, ...]
    values: tuple[float, ...]
    units: tuple[str, ...]
    active_parameter_fingerprint: str
    fingerprint: str


@dataclass(frozen=True, slots=True)
class HomogeneousChemicalSearchCertificate:
    """Compact finite-search evidence bound to one certified local branch."""

    status: str
    selection_label: str
    primary_budget: int
    primary_attempt_count: int
    generated_start_count: int
    budget_truncated_start_count: int
    duplicate_start_count: int
    infeasible_start_count: int
    evaluated_start_count: int
    domain_rejected_start_count: int
    construction_rejected_start_count: int
    attempt_count: int
    basin_count: int
    selected_basin_ordinal: int
    selected_basin_representative_attempt_ordinal: int
    budget_prefix_selection_history: tuple[
        tuple[int, tuple[int, ...], tuple[int, ...], int | None, bool], ...
    ]
    fingerprint: str


@dataclass(frozen=True, slots=True)
class HomogeneousContinuationCertificate:
    """Pure reference-anchored selection evidence for one local branch family."""

    reference_fingerprint: str
    branch_family_identity: str
    search_basin_count: int
    in_region_basin_count: int
    selected_basin_ordinal: int
    selected_composition_distance: float
    selected_volume_distance: float
    fingerprint: str


@dataclass(frozen=True, slots=True)
class HomogeneousReactiveOwnerObservation:
    """One caller-owned row identity bound into an Equilibrium owner receipt."""

    target_identity: str
    prediction_identity: str
    aggregate_identity: str | None = None
    covariance_identity: str | None = None


@dataclass(frozen=True, slots=True)
class HomogeneousPhaseRoleReceipt:
    scope: str
    phase_identity: str
    declared_role: HomogeneousPhaseRole
    normalized_chemical_mole_fractions: tuple[float, ...]
    audited_mole_fractions: tuple[float, ...]
    audited_molar_volume_m3_per_mol: float
    audited_pressure_pa: float
    mechanical_class: str
    route: str
    native_certificate_scope: str
    phase_selection_claim: str
    globality_claim: str
    search_scope: str
    search_budget: str
    search_exclusions: tuple[str, ...]
    completeness_criterion: str
    role_evidence: tuple[str, ...]
    provider_evaluation_count: int
    continuation_reference_fingerprint: str
    active_trial_fingerprint: str
    active_parameter_fingerprint: str
    problem_fingerprint: str
    reaction_fingerprint: str
    topology_fingerprint: str
    fingerprint: str


@dataclass(frozen=True, slots=True)
class CertifiedHomogeneousState:
    phase_identity: str
    phase_role: HomogeneousPhaseRole
    temperature_k: float
    pressure_pa: float
    eos_pressure_pa: float
    amounts_mol: tuple[float, ...]
    mole_fractions: tuple[float, ...]
    volume_m3: float
    molar_density_mol_per_m3: float
    state_identity: str
    problem_fingerprint: str
    reaction_fingerprint: str
    topology_fingerprint: str
    chart_topology: str
    chemical_search_certificate: HomogeneousChemicalSearchCertificate
    continuation_certificate: HomogeneousContinuationCertificate
    phase_role_receipt: HomogeneousPhaseRoleReceipt
    fingerprint: str


@dataclass(frozen=True, slots=True)
class PhysicalBranchReceipt:
    identity: str
    branch_family_identity: str
    continuation_identity: str
    state_identity: str
    topology_fingerprint: str
    problem_fingerprint: str
    reaction_fingerprint: str
    phase_role_certificate_fingerprint: str
    chemical_search_certificate_fingerprint: str
    continuation_certificate_fingerprint: str
    phase_identity: str
    phase_role: HomogeneousPhaseRole
    amounts_mol: tuple[float, ...]
    volume_m3: float


@dataclass(frozen=True, slots=True)
class HomogeneousReactiveObservationResult:
    """Values and exact total parameter Jacobian, or one typed failed trial."""

    descriptor: HomogeneousReactiveObservationDescriptor
    status: Literal["evaluated", "non_evaluable"]
    failure: NonEvaluableTrial | None
    rows: tuple[HomogeneousReactiveObservationRowResult, ...]
    trial: ActiveParameterTrialReceipt | None
    state: CertifiedHomogeneousState | None
    physical_branch: PhysicalBranchReceipt | None
    support_domain_status: str
    equilibrium: ChemicalEquilibriumResult | None
    diagnostics: ChemicalEquilibriumDiagnostics | None

    @property
    def values(self) -> tuple[float | None, ...]:
        return tuple(row.value for row in self.rows)

    @property
    def jacobian(self) -> tuple[tuple[float, ...] | None, ...]:
        return tuple(row.jacobian for row in self.rows)


def _problem_payload(problem: ChemicalEquilibriumProblem) -> dict[str, object]:
    standard = problem.source_standard_state
    return {
        "species_ids": problem.species_ids,
        "charges": problem.charges,
        "molar_masses_kg_per_mol": problem.molar_masses_kg_per_mol,
        "balance_matrix": problem.balance_matrix,
        "conserved_totals": problem.conserved_totals,
        "reaction_matrix": problem.reaction_matrix,
        "feed_amounts_mol": problem.feed_amounts_mol,
        "equilibrium_constants": tuple(
            (
                record.ln_value,
                record.source_id,
                record.reference_id,
                record.reaction_orientation,
                record.conversion_id,
                record.dimensionless,
            )
            for record in problem.equilibrium_constants
        ),
        "strict_interior_amount_floor_mol": problem.strict_interior_amount_floor_mol,
        "source_standard_state": None
        if standard is None
        else {
            "id": standard.id,
            "activity_scale_id": standard.activity_scale_id,
            "log_activity_scale_factors": standard.log_activity_scale_factors,
            "reference_pressure_pa": standard.reference_pressure_pa,
            "component_ids": standard.source_reference_component_ids,
            "solvent_composition": standard.source_reference_solvent_composition,
            "ion_pairs": tuple(
                (pair.component_ids, pair.stoichiometry)
                for pair in standard.source_reference_ion_pairs
            ),
            "activity_convention_id": standard.source_reference_activity_convention_id,
            "standard_molality_mol_per_kg": (
                standard.source_reference_standard_molality_mol_per_kg
            ),
        },
    }


def _validate_row(
    row: HomogeneousReactiveObservationRow,
    problem: ChemicalEquilibriumProblem,
) -> None:
    count = len(problem.species_ids)
    if not row.identity or row.quantity not in _QUANTITY_UNITS:
        raise ValueError("reactive observation identity or quantity is invalid")
    if row.quantity in _COEFFICIENT_QUANTITIES:
        if (
            len(row.coefficients) != count
            or not all(math.isfinite(value) for value in row.coefficients)
            or not any(value != 0.0 for value in row.coefficients)
        ):
            raise ValueError("reactive observation coefficients do not match component order")
    elif row.coefficients:
        raise ValueError("scalar reactive observation cannot declare species coefficients")
    linear_species_aggregate = (
        row.quantity
        in {
            "amount_mol",
            "mole_fraction",
            "molarity_mol_per_m3",
            "molality_mol_per_kg",
        }
        and sum(value != 0.0 for value in row.coefficients) > 1
    )
    if bool(row.aggregate_identity) != bool(row.covariance_identity):
        raise ValueError(
            "aggregate identity and covariance grouping identity must be declared together"
        )
    if linear_species_aggregate and not row.aggregate_identity:
        raise ValueError(
            "linear species aggregates require aggregate and covariance grouping identities"
        )
    if row.quantity == "molality_mol_per_kg":
        if (
            len(row.solvent_mass_coefficients_kg_per_mol) != count
            or not all(
                math.isfinite(value) and value >= 0.0
                for value in row.solvent_mass_coefficients_kg_per_mol
            )
            or not any(value > 0.0 for value in row.solvent_mass_coefficients_kg_per_mol)
        ):
            raise ValueError("molality observation requires a positive declared solvent-mass basis")
    elif row.solvent_mass_coefficients_kg_per_mol:
        raise ValueError("only molality observations declare a solvent-mass basis")
    if row.quantity in _POTENTIAL_QUANTITIES:
        charge = math.fsum(
            coefficient * species_charge
            for coefficient, species_charge in zip(row.coefficients, problem.charges, strict=True)
        )
        if abs(charge) > 1.0e-12:
            raise ValueError("potential contractions must be charge neutral")
    if row.quantity == "fugacity_pa" and (
        sum(value != 0.0 for value in row.coefficients) != 1 or 1.0 not in row.coefficients
    ):
        raise ValueError("fugacity observations require one neutral unit-selector coefficient")
    if row.quantity == "log_activity_contraction" and problem.source_standard_state is None:
        raise ValueError("activity contractions require a declared source standard state")
    if row.support in {"left_censored", "right_censored"}:
        if row.censor_limit is None or not math.isfinite(row.censor_limit):
            raise ValueError("censored observation requires a finite declared limit")
    elif row.censor_limit is not None:
        raise ValueError("uncensored observation cannot declare a censor limit")


def homogeneous_reactive_observation_descriptor(
    problem: HomogeneousReactiveObservationProblem,
    rows: tuple[HomogeneousReactiveObservationRow, ...],
    *,
    active_parameters: epcsaft.ActiveParameterSet | None = None,
) -> HomogeneousReactiveObservationDescriptor:
    """Validate and freeze one source-neutral homogeneous observation family."""

    if not isinstance(problem, HomogeneousReactiveObservationProblem):
        raise TypeError("homogeneous observations require a typed homogeneous problem")
    phase = problem.phase
    reaction_system = problem.reaction_system
    frozen_rows = tuple(rows)
    if not frozen_rows or any(
        not isinstance(row, HomogeneousReactiveObservationRow) for row in frozen_rows
    ):
        raise ValueError("homogeneous observations require at least one typed row")
    if len({row.identity for row in frozen_rows}) != len(frozen_rows):
        raise ValueError("homogeneous observation row identities must be unique")
    for row in frozen_rows:
        _validate_row(row, reaction_system)
    if active_parameters is not None and not isinstance(
        active_parameters, epcsaft.ActiveParameterSet
    ):
        raise TypeError("active_parameters must be an ActiveParameterSet")
    if isinstance(phase, ProviderPhase):
        if tuple(phase.model.component_ids) != reaction_system.species_ids:
            raise ValueError("homogeneous observation component order does not match the Provider")
        if active_parameters is not None:
            phase.model._validate_active_parameters(active_parameters)
        parameter_fingerprint = phase.model.parameter_fingerprint
        topology_fingerprint = _native_topology_fingerprint(phase.model)
        domain_payload: object = {
            "packing_fraction": phase.admissible_packing_fraction_interval,
            "model": phase.model._model.topology,
            "temperature": (
                phase.model._model.temperature_min_k,
                phase.model._model.temperature_max_k,
            ),
        }
    elif isinstance(phase, IdealGasPhase):
        if active_parameters is not None:
            raise ValueError("ideal-gas observations cannot select Provider parameters")
        parameter_fingerprint = phase.model_fingerprint
        topology_fingerprint = _canonical_fingerprint(
            "epcsaft.ideal-gas-topology", reaction_system.species_ids
        )
        domain_payload = {"model": "ideal_gas"}
    else:
        raise TypeError("homogeneous observations require a typed phase model")
    standard = reaction_system.source_standard_state
    reference_fingerprint = _canonical_fingerprint(
        "epcsaft.homogeneous-reactive-reference",
        None if standard is None else _problem_payload(reaction_system)["source_standard_state"],
    )
    problem_fingerprint = _canonical_fingerprint(
        "epcsaft.homogeneous-reactive-problem",
        {
            "identity": problem.identity,
            "phase_identity": problem.phase_identity,
            "phase_role": problem.phase_role,
            "continuation_identity": problem.continuation_identity,
            "branch_policy": problem.branch_policy,
            "continuation_reference": (
                None
                if problem.continuation_reference is None
                else problem.continuation_reference.to_mapping()
            ),
            "reaction_system": _problem_payload(reaction_system),
        },
    )
    reaction_fingerprint = _canonical_fingerprint(
        "epcsaft.homogeneous-reactive-reaction-system",
        _problem_payload(reaction_system),
    )
    active_ids = () if active_parameters is None else active_parameters.identities
    active_units = (
        ()
        if active_parameters is None
        else tuple(parameter.unit for parameter in active_parameters.coordinates)
    )
    descriptor_payload = {
        "component_ids": reaction_system.species_ids,
        "charges": reaction_system.charges,
        "problem_identity": problem.identity,
        "phase_identity": problem.phase_identity,
        "phase_role": problem.phase_role,
        "continuation_identity": problem.continuation_identity,
        "branch_policy": problem.branch_policy,
        "rows": tuple(
            (
                row.identity,
                row.quantity,
                row.coefficients,
                row.solvent_mass_coefficients_kg_per_mol,
                row.support,
                row.censor_limit,
                row.aggregate_identity,
                row.covariance_identity,
            )
            for row in frozen_rows
        ),
        "active_parameter_identities": active_ids,
        "parameter_fingerprint": parameter_fingerprint,
        "topology_fingerprint": topology_fingerprint,
        "reference_fingerprint": reference_fingerprint,
        "problem_fingerprint": problem_fingerprint,
        "reaction_fingerprint": reaction_fingerprint,
    }
    selection_fingerprint = _canonical_fingerprint(
        "epcsaft.active-parameter-selection",
        {
            "parameter_fingerprint": parameter_fingerprint,
            "identities": active_ids,
            "units": active_units,
        },
    )
    return HomogeneousReactiveObservationDescriptor(
        identity="equilibrium.homogeneous-reactive",
        version=2,
        component_ids=reaction_system.species_ids,
        charges=reaction_system.charges,
        problem_identity=problem.identity,
        phase_identity=problem.phase_identity,
        phase_role=problem.phase_role,
        continuation_identity=problem.continuation_identity,
        branch_policy=problem.branch_policy,
        rows=frozen_rows,
        row_identities=tuple(row.identity for row in frozen_rows),
        row_units=tuple(row.unit for row in frozen_rows),
        active_parameter_identities=active_ids,
        active_parameter_units=active_units,
        parameter_fingerprint=parameter_fingerprint,
        active_parameter_selection_fingerprint=selection_fingerprint,
        topology_fingerprint=topology_fingerprint,
        reference_fingerprint=reference_fingerprint,
        domain_fingerprint=_canonical_fingerprint(
            "epcsaft.homogeneous-reactive-domain", domain_payload
        ),
        problem_fingerprint=problem_fingerprint,
        reaction_fingerprint=reaction_fingerprint,
        local_scope="fixed_TP_single_homogeneous_phase",
        jacobian_shape=(len(frozen_rows), len(active_ids)),
        fingerprint=_canonical_fingerprint(
            "epcsaft.homogeneous-reactive-descriptor", descriptor_payload
        ),
    )


def _failed_result(
    descriptor: HomogeneousReactiveObservationDescriptor,
    code: NonEvaluableTrialCode,
    diagnostic: str,
    *,
    trial: ActiveParameterTrialReceipt | None = None,
    state: CertifiedHomogeneousState | None = None,
    equilibrium: ChemicalEquilibriumResult | None = None,
    diagnostics: ChemicalEquilibriumDiagnostics | None = None,
    physical_branch: PhysicalBranchReceipt | None = None,
    support_domain_status: str = "not_evaluated",
) -> HomogeneousReactiveObservationResult:
    rows = tuple(
        HomogeneousReactiveObservationRowResult(
            row.identity, "non_evaluable", None, None, diagnostic
        )
        for row in descriptor.rows
    )
    return HomogeneousReactiveObservationResult(
        descriptor=descriptor,
        status="non_evaluable",
        failure=NonEvaluableTrial(code, diagnostic),
        rows=rows,
        trial=trial,
        state=state,
        physical_branch=physical_branch,
        support_domain_status=support_domain_status,
        equilibrium=equilibrium,
        diagnostics=diagnostics,
    )


_CHEMICAL_FAILURE_CODES: dict[str, NonEvaluableTrialCode] = {
    "genuine_saddle": "second_order_failure",
    "derivative_inconsistency": "first_order_failure",
    "rank_deficiency": "singular_sensitivity",
    "ill_conditioning": "ill_conditioned_sensitivity",
    "physical_domain_failure": "physical_rejection",
    "exhausted_multistart_search": "numerical_convergence_failure",
    "unsupported_derivative_capability": "unsupported_derivative",
    "unsupported_source_reference_transfer": "source_reference_failure",
    "invalid_source_reference": "source_reference_failure",
    "component_order_mismatch": "component_order_mismatch",
    "missing_parameter_or_topology": "parameter_topology_mismatch",
    "provider_domain_rejection": "provider_domain_rejection",
    "reference_limit_convergence_failure": "source_reference_limit_failure",
    "source_reference_derivative_unavailable": ("source_reference_derivative_unavailable"),
    "source_reference_receipt_mismatch": "source_reference_failure",
    "input_or_native_error": "numerical_convergence_failure",
}


def _chemical_failure_code(
    diagnostics: ChemicalEquilibriumDiagnostics,
) -> NonEvaluableTrialCode:
    if diagnostics.provider_domain_status == "failed":
        return "provider_domain_rejection"
    return _CHEMICAL_FAILURE_CODES.get(diagnostics.failure_kind, "numerical_convergence_failure")


def _sensitivity_failure_code(
    sensitivity: object | None,
) -> NonEvaluableTrialCode:
    if sensitivity is None:
        return "first_order_failure"
    active_bounds = tuple(
        index
        for attribute in (
            "active_lower_bounds",
            "active_upper_bounds",
            "active_constraint_bounds",
        )
        for index in getattr(sensitivity, attribute, ())
    )
    reason = str(getattr(sensitivity, "failure_reason", ""))
    if active_bounds or reason in {
        "active_set_change_not_differentiable",
        "provider_domain_boundary_active",
    }:
        return "active_bound_ambiguity"
    if "singular" in reason or "rank" in reason:
        return "singular_sensitivity"
    if "ill_condition" in reason or "condition" in reason:
        return "ill_conditioned_sensitivity"
    if "source_reference" in reason:
        return "source_reference_derivative_unavailable"
    return "first_order_failure"


def _trial_receipt(
    descriptor: HomogeneousReactiveObservationDescriptor,
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
    return ActiveParameterTrialReceipt(
        parameter_fingerprint=descriptor.parameter_fingerprint,
        identities=identities,
        values=values,
        units=units,
        active_parameter_fingerprint=active_fingerprint,
        fingerprint=_canonical_fingerprint(
            "epcsaft.active-parameter-trial",
            {
                "parameter_fingerprint": descriptor.parameter_fingerprint,
                "identities": identities,
                "values": values,
                "units": units,
                "active_parameter_fingerprint": active_fingerprint,
            },
        ),
    )


def _neutral_basis_coordinates(
    coefficients: tuple[float, ...],
    basis: tuple[tuple[float, ...], ...],
) -> tuple[float, ...]:
    size = len(basis)
    matrix = [
        [
            math.fsum(left * right for left, right in zip(basis[row], basis[column], strict=True))
            for column in range(size)
        ]
        for row in range(size)
    ]
    right_hand_side = [
        math.fsum(value * coefficient for value, coefficient in zip(row, coefficients, strict=True))
        for row in basis
    ]
    for pivot in range(size):
        selected = max(range(pivot, size), key=lambda row: abs(matrix[row][pivot]))
        if abs(matrix[selected][pivot]) <= 1.0e-14:
            raise ValueError("source-reference neutral basis is singular")
        matrix[pivot], matrix[selected] = matrix[selected], matrix[pivot]
        right_hand_side[pivot], right_hand_side[selected] = (
            right_hand_side[selected],
            right_hand_side[pivot],
        )
        diagonal = matrix[pivot][pivot]
        for column in range(pivot, size):
            matrix[pivot][column] /= diagonal
        right_hand_side[pivot] /= diagonal
        for row in range(size):
            if row == pivot:
                continue
            factor = matrix[row][pivot]
            for column in range(pivot, size):
                matrix[row][column] -= factor * matrix[pivot][column]
            right_hand_side[row] -= factor * right_hand_side[pivot]
    reconstructed = tuple(
        math.fsum(
            coordinate * basis_row[component]
            for coordinate, basis_row in zip(right_hand_side, basis, strict=True)
        )
        for component in range(len(coefficients))
    )
    if (
        max(
            abs(value - coefficient)
            for value, coefficient in zip(reconstructed, coefficients, strict=True)
        )
        > 1.0e-10
    ):
        raise ValueError("contraction is outside the charge-neutral reference span")
    return tuple(right_hand_side)


def _chemical_search_certificate_payload(
    certificate: HomogeneousChemicalSearchCertificate,
) -> dict[str, object]:
    return {
        "status": certificate.status,
        "selection_label": certificate.selection_label,
        "primary_budget": certificate.primary_budget,
        "primary_attempt_count": certificate.primary_attempt_count,
        "generated_start_count": certificate.generated_start_count,
        "budget_truncated_start_count": certificate.budget_truncated_start_count,
        "duplicate_start_count": certificate.duplicate_start_count,
        "infeasible_start_count": certificate.infeasible_start_count,
        "evaluated_start_count": certificate.evaluated_start_count,
        "domain_rejected_start_count": certificate.domain_rejected_start_count,
        "construction_rejected_start_count": certificate.construction_rejected_start_count,
        "attempt_count": certificate.attempt_count,
        "basin_count": certificate.basin_count,
        "selected_basin_ordinal": certificate.selected_basin_ordinal,
        "selected_basin_representative_attempt_ordinal": (
            certificate.selected_basin_representative_attempt_ordinal
        ),
        "budget_prefix_selection_history": certificate.budget_prefix_selection_history,
    }


def _chemical_search_certificate(
    search: ChemicalEquilibriumSearch,
) -> HomogeneousChemicalSearchCertificate:
    selected_ordinal = search.selected_basin_ordinal
    if selected_ordinal is None or not 0 <= selected_ordinal < len(search.basins):
        raise ValueError("certified chemical search has no selected local basin")
    selected = search.basins[selected_ordinal]
    certificate = HomogeneousChemicalSearchCertificate(
        status=search.status,
        selection_label=search.selection_label,
        primary_budget=search.primary_budget,
        primary_attempt_count=search.primary_attempt_count,
        generated_start_count=search.generated_start_count,
        budget_truncated_start_count=search.budget_truncated_start_count,
        duplicate_start_count=search.duplicate_start_count,
        infeasible_start_count=search.infeasible_start_count,
        evaluated_start_count=search.evaluated_start_count,
        domain_rejected_start_count=search.domain_rejected_start_count,
        construction_rejected_start_count=search.construction_rejected_start_count,
        attempt_count=len(search.attempts),
        basin_count=len(search.basins),
        selected_basin_ordinal=selected_ordinal,
        selected_basin_representative_attempt_ordinal=(selected.representative_attempt_ordinal),
        budget_prefix_selection_history=tuple(
            (
                prefix.primary_budget,
                prefix.attempted_primary_ordinals,
                prefix.basin_ordinals,
                prefix.selected_basin_ordinal,
                prefix.selection_changed,
            )
            for prefix in search.budget_prefixes
        ),
        fingerprint="",
    )
    return replace(
        certificate,
        fingerprint=_canonical_fingerprint(
            "epcsaft.homogeneous-chemical-search-certificate",
            _chemical_search_certificate_payload(certificate),
        ),
    )


def _validate_chemical_search_certificate(
    certificate: HomogeneousChemicalSearchCertificate,
) -> None:
    counts = (
        certificate.primary_budget,
        certificate.primary_attempt_count,
        certificate.generated_start_count,
        certificate.budget_truncated_start_count,
        certificate.duplicate_start_count,
        certificate.infeasible_start_count,
        certificate.evaluated_start_count,
        certificate.domain_rejected_start_count,
        certificate.construction_rejected_start_count,
        certificate.attempt_count,
        certificate.basin_count,
    )
    previous_budget = 0
    previous_selection: int | None = None
    history_valid = True
    expected_prefix_budgets = tuple(
        dict.fromkeys(
            min(requested, certificate.primary_attempt_count) for requested in (1, 5, 13, 25)
        )
    )
    for index, record in enumerate(certificate.budget_prefix_selection_history):
        budget, attempts, basins, selection, changed = record
        if (
            budget <= previous_budget
            or budget > certificate.primary_attempt_count
            or attempts != tuple(range(budget))
            or len(set(basins)) != len(basins)
            or any(not 0 <= basin < certificate.basin_count for basin in basins)
            or (selection is not None and selection not in basins)
            or changed != (index > 0 and selection != previous_selection)
        ):
            history_valid = False
            break
        previous_budget = budget
        previous_selection = selection
    expected_fingerprint = _canonical_fingerprint(
        "epcsaft.homogeneous-chemical-search-certificate",
        _chemical_search_certificate_payload(certificate),
    )
    if (
        certificate.status != "certified_local_minimum"
        or certificate.selection_label != "lowest_observed_certified_local_value"
        or any(count < 0 for count in counts)
        or not 0 < certificate.primary_attempt_count <= certificate.primary_budget <= 25
        or certificate.generated_start_count
        != certificate.primary_attempt_count
        + certificate.budget_truncated_start_count
        + certificate.duplicate_start_count
        + certificate.infeasible_start_count
        or certificate.evaluated_start_count
        + certificate.domain_rejected_start_count
        + certificate.construction_rejected_start_count
        != certificate.primary_attempt_count
        or certificate.attempt_count < certificate.primary_attempt_count
        or certificate.basin_count < 1
        or not 0 <= certificate.selected_basin_ordinal < certificate.basin_count
        or not 0
        <= certificate.selected_basin_representative_attempt_ordinal
        < certificate.attempt_count
        or not certificate.budget_prefix_selection_history
        or tuple(record[0] for record in certificate.budget_prefix_selection_history)
        != expected_prefix_budgets
        or certificate.budget_prefix_selection_history[-1][3]
        != certificate.selected_basin_ordinal
        or not history_valid
        or certificate.fingerprint != expected_fingerprint
    ):
        raise ValueError("homogeneous chemical search certificate is inconsistent")


def _continuation_certificate_payload(
    certificate: HomogeneousContinuationCertificate,
) -> dict[str, object]:
    return {
        "reference_fingerprint": certificate.reference_fingerprint,
        "branch_family_identity": certificate.branch_family_identity,
        "search_basin_count": certificate.search_basin_count,
        "in_region_basin_count": certificate.in_region_basin_count,
        "selected_basin_ordinal": certificate.selected_basin_ordinal,
        "selected_composition_distance": certificate.selected_composition_distance,
        "selected_volume_distance": certificate.selected_volume_distance,
    }


def _certify_continuation(
    problem: HomogeneousReactiveObservationProblem,
    descriptor: HomogeneousReactiveObservationDescriptor,
    equilibrium: ChemicalEquilibriumResult,
) -> HomogeneousContinuationCertificate:
    search = equilibrium.diagnostics.search
    selected = search.selected_basin_ordinal
    if selected is None or not 0 <= selected < len(search.basins):
        raise _ObservationFailure("branch_loss", "chemical search selected no certified basin")
    reference = problem.continuation_reference
    if reference is None:
        if len(search.basins) != 1:
            raise _ObservationFailure(
                "multiple_unresolved_branches",
                "unanchored evaluation found more than one certified local basin",
            )
        selected_composition_distance = 0.0
        selected_volume_distance = 0.0
        reference_fingerprint = _canonical_fingerprint(
            "epcsaft.homogeneous-unanchored-state-reference",
            {
                "amounts_mol": search.basins[selected].amounts_mol,
                "volume_m3": search.basins[selected].volume_m3,
            },
        )
        in_region = (selected,)
    else:
        anchor_catalog_metrics = tuple(
            (
                max(
                    abs(math.log(actual / anchor))
                    for actual, anchor in zip(
                        composition,
                        reference.mole_fractions,
                        strict=True,
                    )
                ),
                abs(
                    math.log(
                        molar_volume / reference.molar_volume_m3_per_mol
                    )
                ),
            )
            for composition, molar_volume in reference.anchor_basin_states
        )
        catalog_selected = tuple(
            index
            for index, (composition_distance, volume_distance) in enumerate(
                anchor_catalog_metrics
            )
            if composition_distance <= 2.0e-10
            and volume_distance <= 2.0e-10
        )
        if len(catalog_selected) != 1:
            raise _ObservationFailure(
                "branch_loss",
                "immutable continuation catalog lost its selected anchor basin",
            )
        if any(
            index != catalog_selected[0]
            and composition_distance <= reference.maximum_log_composition_distance
            and volume_distance <= reference.maximum_log_volume_distance
            for index, (composition_distance, volume_distance) in enumerate(
                anchor_catalog_metrics
            )
        ):
            raise _ObservationFailure(
                "multiple_unresolved_branches",
                "an alternative anchor-certified basin lies inside the continuation region",
            )
        metrics: list[tuple[int, float, float]] = []
        for basin in search.basins:
            total_amount = math.fsum(basin.amounts_mol)
            if total_amount <= 0.0 or basin.volume_m3 <= 0.0:
                continue
            mole_fractions = tuple(amount / total_amount for amount in basin.amounts_mol)
            if any(actual <= 0.0 for actual in mole_fractions):
                continue
            composition_distance = max(
                abs(math.log(actual / anchor))
                for actual, anchor in zip(
                    mole_fractions,
                    reference.mole_fractions,
                    strict=True,
                )
            )
            volume_distance = abs(
                math.log(
                    basin.volume_m3
                    / total_amount
                    / reference.molar_volume_m3_per_mol
                )
            )
            metrics.append((basin.ordinal, composition_distance, volume_distance))
        in_region_metrics = tuple(
            item
            for item in metrics
            if item[1] <= reference.maximum_log_composition_distance
            and item[2] <= reference.maximum_log_volume_distance
        )
        if not in_region_metrics:
            raise _ObservationFailure(
                "branch_loss",
                "no certified basin remains inside the immutable continuation region",
            )
        if len(in_region_metrics) != 1:
            raise _ObservationFailure(
                "multiple_unresolved_branches",
                "multiple certified basins lie inside the immutable continuation region",
            )
        selected_ordinal, selected_composition_distance, selected_volume_distance = (
            in_region_metrics[0]
        )
        if selected_ordinal != selected:
            raise _ObservationFailure(
                "branch_loss",
                "the reference-selected basin differs from the chemical search selection",
            )
        reference_fingerprint = reference.fingerprint
        in_region = (selected_ordinal,)
    branch_family_identity = _canonical_fingerprint(
        "epcsaft.homogeneous-reactive-branch-family",
        {
            "continuation_identity": descriptor.continuation_identity,
            "reference_fingerprint": reference_fingerprint,
            "problem_fingerprint": descriptor.problem_fingerprint,
            "reaction_fingerprint": descriptor.reaction_fingerprint,
            "topology_fingerprint": descriptor.topology_fingerprint,
            "phase_identity": descriptor.phase_identity,
            "phase_role": descriptor.phase_role,
            "branch_policy": descriptor.branch_policy,
        },
    )
    certificate = HomogeneousContinuationCertificate(
        reference_fingerprint=reference_fingerprint,
        branch_family_identity=branch_family_identity,
        search_basin_count=len(search.basins),
        in_region_basin_count=len(in_region),
        selected_basin_ordinal=selected,
        selected_composition_distance=selected_composition_distance,
        selected_volume_distance=selected_volume_distance,
        fingerprint="",
    )
    return replace(
        certificate,
        fingerprint=_canonical_fingerprint(
            "epcsaft.homogeneous-continuation-certificate",
            _continuation_certificate_payload(certificate),
        ),
    )


def _immutable_continuation_start(
    problem: HomogeneousReactiveObservationProblem,
) -> tuple[tuple[float, ...], float] | None:
    reference = problem.continuation_reference
    if reference is None:
        return None
    projected_balances = tuple(
        math.fsum(
            coefficient * mole_fraction
            for coefficient, mole_fraction in zip(
                row,
                reference.mole_fractions,
                strict=True,
            )
        )
        for row in problem.reaction_system.balance_matrix
    )
    denominator = math.fsum(value * value for value in projected_balances)
    if not math.isfinite(denominator) or denominator <= 0.0:
        raise ValueError("continuation composition does not determine the conserved scale")
    total_amount = math.fsum(
        coefficient * total
        for coefficient, total in zip(
            projected_balances,
            problem.reaction_system.conserved_totals,
            strict=True,
        )
    ) / denominator
    amounts = tuple(total_amount * value for value in reference.mole_fractions)
    balance_residual = max(
        abs(
            math.fsum(
                coefficient * amount
                for coefficient, amount in zip(row, amounts, strict=True)
            )
            - expected
        )
        for row, expected in zip(
            problem.reaction_system.balance_matrix,
            problem.reaction_system.conserved_totals,
            strict=True,
        )
    )
    charge_residual = abs(
        math.fsum(
            charge * amount
            for charge, amount in zip(
                problem.reaction_system.charges,
                amounts,
                strict=True,
            )
        )
    )
    scale = max(1.0, *(abs(value) for value in problem.reaction_system.conserved_totals))
    volume = total_amount * reference.molar_volume_m3_per_mol
    if (
        not math.isfinite(total_amount)
        or total_amount <= 0.0
        or any(
            not math.isfinite(amount)
            or amount <= problem.reaction_system.strict_interior_amount_floor_mol
            for amount in amounts
        )
        or not math.isfinite(volume)
        or volume <= 0.0
        or balance_residual > 2.0e-10 * scale
        or charge_residual > 2.0e-10 * scale
    ):
        raise ValueError(
            "immutable continuation state is inconsistent with the reaction constraints"
        )
    return amounts, volume


def _validate_continuation_certificate(
    descriptor: HomogeneousReactiveObservationDescriptor,
    certificate: HomogeneousContinuationCertificate,
) -> None:
    expected_family_identity = _canonical_fingerprint(
        "epcsaft.homogeneous-reactive-branch-family",
        {
            "continuation_identity": descriptor.continuation_identity,
            "reference_fingerprint": certificate.reference_fingerprint,
            "problem_fingerprint": descriptor.problem_fingerprint,
            "reaction_fingerprint": descriptor.reaction_fingerprint,
            "topology_fingerprint": descriptor.topology_fingerprint,
            "phase_identity": descriptor.phase_identity,
            "phase_role": descriptor.phase_role,
            "branch_policy": descriptor.branch_policy,
        },
    )
    expected_fingerprint = _canonical_fingerprint(
        "epcsaft.homogeneous-continuation-certificate",
        _continuation_certificate_payload(certificate),
    )
    if (
        certificate.branch_family_identity != expected_family_identity
        or certificate.search_basin_count < 1
        or certificate.in_region_basin_count != 1
        or not 0 <= certificate.selected_basin_ordinal < certificate.search_basin_count
        or certificate.selected_composition_distance < 0.0
        or certificate.selected_volume_distance < 0.0
        or not math.isfinite(certificate.selected_composition_distance)
        or not math.isfinite(certificate.selected_volume_distance)
        or certificate.fingerprint != expected_fingerprint
    ):
        raise ValueError("homogeneous continuation certificate is inconsistent")


def _phase_role_payload(receipt: HomogeneousPhaseRoleReceipt) -> dict[str, object]:
    return {
        "scope": receipt.scope,
        "phase_identity": receipt.phase_identity,
        "declared_role": receipt.declared_role,
        "normalized_chemical_mole_fractions": receipt.normalized_chemical_mole_fractions,
        "audited_mole_fractions": receipt.audited_mole_fractions,
        "audited_molar_volume_m3_per_mol": receipt.audited_molar_volume_m3_per_mol,
        "audited_pressure_pa": receipt.audited_pressure_pa,
        "mechanical_class": receipt.mechanical_class,
        "route": receipt.route,
        "native_certificate_scope": receipt.native_certificate_scope,
        "phase_selection_claim": receipt.phase_selection_claim,
        "globality_claim": receipt.globality_claim,
        "search_scope": receipt.search_scope,
        "search_budget": receipt.search_budget,
        "search_exclusions": receipt.search_exclusions,
        "completeness_criterion": receipt.completeness_criterion,
        "role_evidence": receipt.role_evidence,
        "provider_evaluation_count": receipt.provider_evaluation_count,
        "continuation_reference_fingerprint": receipt.continuation_reference_fingerprint,
        "active_trial_fingerprint": receipt.active_trial_fingerprint,
        "active_parameter_fingerprint": receipt.active_parameter_fingerprint,
        "problem_fingerprint": receipt.problem_fingerprint,
        "reaction_fingerprint": receipt.reaction_fingerprint,
        "topology_fingerprint": receipt.topology_fingerprint,
    }


def _certify_phase_role(
    problem: HomogeneousReactiveObservationProblem,
    descriptor: HomogeneousReactiveObservationDescriptor,
    equilibrium: ChemicalEquilibriumResult,
    active_parameters: epcsaft.ActiveParameterSet | None,
    trial: ActiveParameterTrialReceipt,
    pressure_derivative_log_volume_pa: float | None,
    continuation_reference_fingerprint: str,
) -> HomogeneousPhaseRoleReceipt:
    total_amount = math.fsum(equilibrium.amounts_mol)
    if (
        not math.isfinite(total_amount)
        or total_amount <= 0.0
        or not math.isfinite(equilibrium.volume_m3)
        or equilibrium.volume_m3 <= 0.0
        or any(not math.isfinite(amount) or amount <= 0.0 for amount in equilibrium.amounts_mol)
    ):
        raise ValueError("phase-role certification requires finite strictly positive amounts")
    normalized = tuple(amount / total_amount for amount in equilibrium.amounts_mol)
    certified_molar_volume = equilibrium.volume_m3 / total_amount

    def receipt(
        *,
        scope: str,
        audited_mole_fractions: tuple[float, ...],
        audited_molar_volume_m3_per_mol: float,
        audited_pressure_pa: float,
        mechanical_class: str,
        route: str,
        native_certificate_scope: str,
        phase_selection_claim: str,
        globality_claim: str,
        search_scope: str,
        search_budget: str,
        search_exclusions: tuple[str, ...],
        completeness_criterion: str,
        role_evidence: tuple[str, ...],
        provider_evaluation_count: int,
    ) -> HomogeneousPhaseRoleReceipt:
        value = HomogeneousPhaseRoleReceipt(
            scope=scope,
            phase_identity=descriptor.phase_identity,
            declared_role=descriptor.phase_role,
            normalized_chemical_mole_fractions=normalized,
            audited_mole_fractions=audited_mole_fractions,
            audited_molar_volume_m3_per_mol=audited_molar_volume_m3_per_mol,
            audited_pressure_pa=audited_pressure_pa,
            mechanical_class=mechanical_class,
            route=route,
            native_certificate_scope=native_certificate_scope,
            phase_selection_claim=phase_selection_claim,
            globality_claim=globality_claim,
            search_scope=search_scope,
            search_budget=search_budget,
            search_exclusions=search_exclusions,
            completeness_criterion=completeness_criterion,
            role_evidence=role_evidence,
            provider_evaluation_count=provider_evaluation_count,
            continuation_reference_fingerprint=continuation_reference_fingerprint,
            active_trial_fingerprint=trial.fingerprint,
            active_parameter_fingerprint=trial.active_parameter_fingerprint,
            problem_fingerprint=descriptor.problem_fingerprint,
            reaction_fingerprint=descriptor.reaction_fingerprint,
            topology_fingerprint=descriptor.topology_fingerprint,
            fingerprint="",
        )
        return replace(
            value,
            fingerprint=_canonical_fingerprint(
                "epcsaft.homogeneous-phase-role-certificate",
                _phase_role_payload(value),
            ),
        )

    if isinstance(problem.phase, IdealGasPhase):
        return receipt(
            scope="typed_ideal_gas_fixed_composition_role",
            audited_mole_fractions=normalized,
            audited_molar_volume_m3_per_mol=certified_molar_volume,
            audited_pressure_pa=equilibrium.pressure_pa,
            mechanical_class="strict_stable",
            route="typed_ideal_gas_model",
            native_certificate_scope="analytical_ideal_gas_state",
            phase_selection_claim="typed_ideal_gas_vapor_role",
            globality_claim="not_established",
            search_scope="exact_certified_chemical_state_only",
            search_budget="not_applicable",
            search_exclusions=("other_compositions", "other_phase_counts"),
            completeness_criterion="no_global_phase_stability_or_equilibrium_claim",
            role_evidence=("ideal_gas_pressure_volume_identity",),
            provider_evaluation_count=0,
        )
    reference = problem.continuation_reference
    if reference is not None and pressure_derivative_log_volume_pa is not None:
        mechanical_scale = max(1.0, abs(equilibrium.pressure_pa))
        if (
            reference.phase_identity != descriptor.phase_identity
            or reference.phase_role != descriptor.phase_role
            or reference.anchor_parameter_fingerprint != trial.parameter_fingerprint
            or reference.anchor_reaction_fingerprint != descriptor.reaction_fingerprint
            or reference.anchor_topology_fingerprint != descriptor.topology_fingerprint
            or not math.isfinite(pressure_derivative_log_volume_pa)
            or pressure_derivative_log_volume_pa >= -1.0e-6 * mechanical_scale
        ):
            raise ValueError(
                "immutable role anchor does not certify this exact local Provider state"
            )
        return receipt(
            scope="local_reference_anchored_fixed_composition_role",
            audited_mole_fractions=normalized,
            audited_molar_volume_m3_per_mol=certified_molar_volume,
            audited_pressure_pa=equilibrium.pressure_pa,
            mechanical_class="strict_stable",
            route="exact_fixed_TnV_provider_hessian",
            native_certificate_scope=(
                "exact_state_local_mechanical_stability_plus_role_certified_anchor"
            ),
            phase_selection_claim="reference_anchored_locally_certified_declared_role",
            globality_claim="not_established",
            search_scope="immutable_role_anchor_and_exact_certified_chemical_state",
            search_budget="one_exact_state_provider_derivative_evaluation",
            search_exclusions=(
                "global_phase_stability",
                "other_compositions",
                "other_phase_counts",
                "root_completeness_proof",
            ),
            completeness_criterion=(
                "one_certified_chemical_basin_inside_declared_reference_region;"
                "global_phase_stability_not_established"
            ),
            role_evidence=(
                "anchor_state_identity:" + reference.anchor_state_identity,
                "anchor_phase_role_receipt:"
                + reference.anchor_phase_role_receipt_fingerprint,
                "exact_state_negative_pressure_derivative_log_volume",
                "exact_state_matches_reference_selected_chemical_branch",
            ),
            provider_evaluation_count=1,
        )
    role_result = phase_equilibrium(
        problem.phase.model,
        PhaseEquilibriumProblem(
            FixedTPFeed(
                equilibrium.temperature_k * epcsaft.unit_registry.kelvin,
                equilibrium.pressure_pa * epcsaft.unit_registry.pascal,
                equilibrium.amounts_mol,
                DeclaredTopology((PhaseSpecification(problem.phase_role, AllComponents()),)),
            )
        ),
        active_parameters=active_parameters,
    )
    if len(role_result.phases) != 1:
        raise ValueError("declared phase-role route did not return exactly one phase")
    audited = role_result.phases[0]
    composition_error = max(
        abs(actual - expected)
        for actual, expected in zip(
            audited.mole_fractions,
            normalized,
            strict=True,
        )
    )
    volume_error = abs(audited.molar_volume_m3_per_mol - certified_molar_volume) / max(
        certified_molar_volume, 1.0e-300
    )
    pressure_error = abs(audited.pressure_pa - equilibrium.pressure_pa) / max(
        equilibrium.pressure_pa, 1.0
    )
    if (
        audited.role != problem.phase_role
        or audited.support != "all_components"
        or audited.supported_component_ids != descriptor.component_ids
        or audited.mechanical_class != "strict_stable"
        or composition_error > 2.0e-12
        or volume_error > 2.0e-8
        or pressure_error > 2.0e-8
        or role_result.diagnostics.physical_status != "accepted"
        or role_result.diagnostics.phase_selection_claim
        != "mechanically_certified_declared_single_phase"
        or role_result.diagnostics.globality_claim != "not_established"
        or role_result.diagnostics.completeness_criterion
        != "detected_roots_only;pressure_envelope_root_completeness_not_proven"
        or "detected_strict_stable_volume_extreme_for_declared_role"
        not in role_result.diagnostics.role_evidence
        or (
            "minimum_detected_strict_stable_volume_selected"
            if problem.phase_role == "liquid"
            else "maximum_detected_strict_stable_volume_selected"
        )
        not in role_result.diagnostics.role_evidence
    ):
        raise ValueError(
            "declared Provider phase role does not certify the exact chemical state branch"
        )
    return receipt(
        scope="local_detected_fixed_composition_role",
        audited_mole_fractions=audited.mole_fractions,
        audited_molar_volume_m3_per_mol=audited.molar_volume_m3_per_mol,
        audited_pressure_pa=audited.pressure_pa,
        mechanical_class=audited.mechanical_class,
        route=role_result.diagnostics.route,
        native_certificate_scope=role_result.diagnostics.certificate_scope,
        phase_selection_claim=role_result.diagnostics.phase_selection_claim,
        globality_claim=role_result.diagnostics.globality_claim,
        search_scope=role_result.diagnostics.search_scope,
        search_budget=role_result.diagnostics.search_budget,
        search_exclusions=role_result.diagnostics.search_exclusions,
        completeness_criterion=role_result.diagnostics.completeness_criterion,
        role_evidence=role_result.diagnostics.role_evidence,
        provider_evaluation_count=role_result.diagnostics.provider_evaluation_count,
    )


def certify_homogeneous_continuation_reference(
    problem: HomogeneousReactiveObservationProblem,
    temperature: Quantity[Any],
    pressure: Quantity[Any],
    *,
    maximum_log_composition_distance: float,
    maximum_log_volume_distance: float,
    active_parameters: epcsaft.ActiveParameterSet | None = None,
) -> HomogeneousContinuationReference:
    """Certify one immutable local role anchor using the full detected-root route."""

    if problem.continuation_reference is not None:
        raise ValueError("continuation anchoring requires an unanchored homogeneous problem")
    anchor_row = HomogeneousReactiveObservationRow(
        "continuation-anchor-pressure", "pressure_pa"
    )
    descriptor = homogeneous_reactive_observation_descriptor(
        problem,
        (anchor_row,),
        active_parameters=active_parameters,
    )
    trial = _trial_receipt(descriptor, active_parameters)
    equilibrium = chemical_equilibrium(
        problem.phase,
        temperature,
        pressure,
        problem.reaction_system,
        sensitivity_request=ChemicalEquilibriumSensitivityRequest(
            active_parameters
        ),
    )
    role_receipt = _certify_phase_role(
        problem,
        descriptor,
        equilibrium,
        active_parameters,
        trial,
        None,
        "not_applicable_anchor_construction",
    )
    anchor_state_identity = _canonical_fingerprint(
        "epcsaft.homogeneous-continuation-anchor-state",
        {
            "phase_identity": descriptor.phase_identity,
            "phase_role": descriptor.phase_role,
            "temperature_k": equilibrium.temperature_k,
            "pressure_pa": equilibrium.pressure_pa,
            "amounts_mol": equilibrium.amounts_mol,
            "mole_fractions": equilibrium.mole_fractions,
            "volume_m3": equilibrium.volume_m3,
            "parameter_fingerprint": trial.parameter_fingerprint,
            "reaction_fingerprint": descriptor.reaction_fingerprint,
            "topology_fingerprint": descriptor.topology_fingerprint,
            "phase_role_receipt_fingerprint": role_receipt.fingerprint,
        },
    )
    return HomogeneousContinuationReference(
        equilibrium.mole_fractions,
        equilibrium.volume_m3 / math.fsum(equilibrium.amounts_mol),
        maximum_log_composition_distance,
        maximum_log_volume_distance,
        descriptor.phase_identity,
        descriptor.phase_role,
        anchor_state_identity,
        role_receipt.fingerprint,
        trial.parameter_fingerprint,
        descriptor.reaction_fingerprint,
        descriptor.topology_fingerprint,
        tuple(
            (
                tuple(
                    amount / math.fsum(basin.amounts_mol)
                    for amount in basin.amounts_mol
                ),
                basin.volume_m3 / math.fsum(basin.amounts_mol),
            )
            for basin in equilibrium.diagnostics.search.basins
        ),
    )


def _state_identity_payload(
    descriptor: HomogeneousReactiveObservationDescriptor,
    *,
    temperature_k: float,
    pressure_pa: float,
    eos_pressure_pa: float,
    amounts_mol: tuple[float, ...],
    mole_fractions: tuple[float, ...],
    volume_m3: float,
    molar_density_mol_per_m3: float,
    chart_topology: str,
    phase_role_certificate_fingerprint: str,
) -> dict[str, object]:
    return {
        "phase_identity": descriptor.phase_identity,
        "phase_role": descriptor.phase_role,
        "problem_fingerprint": descriptor.problem_fingerprint,
        "reaction_fingerprint": descriptor.reaction_fingerprint,
        "topology_fingerprint": descriptor.topology_fingerprint,
        "temperature_k": temperature_k,
        "pressure_pa": pressure_pa,
        "eos_pressure_pa": eos_pressure_pa,
        "amounts_mol": amounts_mol,
        "mole_fractions": mole_fractions,
        "volume_m3": volume_m3,
        "molar_density_mol_per_m3": molar_density_mol_per_m3,
        "chart_topology": chart_topology,
        "phase_role_certificate_fingerprint": phase_role_certificate_fingerprint,
    }


def _make_certified_receipts(
    descriptor: HomogeneousReactiveObservationDescriptor,
    trial: ActiveParameterTrialReceipt,
    equilibrium: ChemicalEquilibriumResult,
    eos_pressure_pa: float,
    phase_role_receipt: HomogeneousPhaseRoleReceipt,
    continuation_certificate: HomogeneousContinuationCertificate,
) -> tuple[CertifiedHomogeneousState, PhysicalBranchReceipt]:
    if equilibrium.diagnostics.search.selected_basin_ordinal is None:
        raise ValueError("certified chemical result has no selected physical state")
    search_certificate = _chemical_search_certificate(equilibrium.diagnostics.search)
    _validate_chemical_search_certificate(search_certificate)
    _validate_continuation_certificate(descriptor, continuation_certificate)
    amounts = equilibrium.amounts_mol
    total_amount = math.fsum(amounts)
    state_payload = _state_identity_payload(
        descriptor,
        temperature_k=equilibrium.temperature_k,
        pressure_pa=equilibrium.pressure_pa,
        eos_pressure_pa=eos_pressure_pa,
        amounts_mol=amounts,
        mole_fractions=equilibrium.mole_fractions,
        volume_m3=equilibrium.volume_m3,
        molar_density_mol_per_m3=total_amount / equilibrium.volume_m3,
        chart_topology=equilibrium.diagnostics.chart_topology,
        phase_role_certificate_fingerprint=phase_role_receipt.fingerprint,
    )
    state_identity = _canonical_fingerprint(
        "epcsaft.homogeneous-reactive-state-point", state_payload
    )
    state = CertifiedHomogeneousState(
        phase_identity=descriptor.phase_identity,
        phase_role=descriptor.phase_role,
        temperature_k=equilibrium.temperature_k,
        pressure_pa=equilibrium.pressure_pa,
        eos_pressure_pa=eos_pressure_pa,
        amounts_mol=amounts,
        mole_fractions=equilibrium.mole_fractions,
        volume_m3=equilibrium.volume_m3,
        molar_density_mol_per_m3=total_amount / equilibrium.volume_m3,
        state_identity=state_identity,
        problem_fingerprint=descriptor.problem_fingerprint,
        reaction_fingerprint=descriptor.reaction_fingerprint,
        topology_fingerprint=descriptor.topology_fingerprint,
        chart_topology=equilibrium.diagnostics.chart_topology,
        chemical_search_certificate=search_certificate,
        continuation_certificate=continuation_certificate,
        phase_role_receipt=phase_role_receipt,
        fingerprint=_canonical_fingerprint(
            "epcsaft.certified-homogeneous-state",
            {
                "state_identity": state_identity,
                "trial_fingerprint": trial.fingerprint,
                "chemical_search_certificate_fingerprint": search_certificate.fingerprint,
                "continuation_certificate_fingerprint": continuation_certificate.fingerprint,
            },
        ),
    )
    branch_payload = {
        "phase_identity": descriptor.phase_identity,
        "phase_role": descriptor.phase_role,
        "problem_fingerprint": descriptor.problem_fingerprint,
        "reaction_fingerprint": descriptor.reaction_fingerprint,
        "topology_fingerprint": descriptor.topology_fingerprint,
        "continuation_identity": descriptor.continuation_identity,
        "chart_topology": equilibrium.diagnostics.chart_topology,
        "canonical_thermodynamic_state": state_identity,
        "phase_role_certificate_fingerprint": phase_role_receipt.fingerprint,
        "branch_family_identity": continuation_certificate.branch_family_identity,
    }
    branch = PhysicalBranchReceipt(
        identity=_canonical_fingerprint(
            "epcsaft.homogeneous-reactive-physical-branch", branch_payload
        ),
        branch_family_identity=continuation_certificate.branch_family_identity,
        continuation_identity=descriptor.continuation_identity,
        state_identity=state_identity,
        topology_fingerprint=descriptor.topology_fingerprint,
        problem_fingerprint=descriptor.problem_fingerprint,
        reaction_fingerprint=descriptor.reaction_fingerprint,
        phase_role_certificate_fingerprint=phase_role_receipt.fingerprint,
        chemical_search_certificate_fingerprint=search_certificate.fingerprint,
        continuation_certificate_fingerprint=continuation_certificate.fingerprint,
        phase_identity=descriptor.phase_identity,
        phase_role=descriptor.phase_role,
        amounts_mol=amounts,
        volume_m3=equilibrium.volume_m3,
    )
    _validate_homogeneous_receipts(descriptor, trial, state, branch)
    return state, branch


def _validate_phase_role_receipt(
    descriptor: HomogeneousReactiveObservationDescriptor,
    trial: ActiveParameterTrialReceipt,
    state: CertifiedHomogeneousState,
) -> None:
    receipt = state.phase_role_receipt
    total_amount = math.fsum(state.amounts_mol)
    if (
        not math.isfinite(total_amount)
        or total_amount <= 0.0
        or any(not math.isfinite(amount) or amount <= 0.0 for amount in state.amounts_mol)
    ):
        raise ValueError("homogeneous phase-role receipt has invalid chemical amounts")
    normalized = tuple(amount / total_amount for amount in state.amounts_mol)
    composition_error = max(
        abs(actual - expected)
        for actual, expected in zip(
            receipt.audited_mole_fractions,
            normalized,
            strict=True,
        )
    )
    volume_error = abs(
        receipt.audited_molar_volume_m3_per_mol - state.volume_m3 / total_amount
    ) / max(state.volume_m3 / total_amount, 1.0e-300)
    pressure_error = abs(receipt.audited_pressure_pa - state.pressure_pa) / max(
        state.pressure_pa, 1.0
    )
    expected_fingerprint = _canonical_fingerprint(
        "epcsaft.homogeneous-phase-role-certificate",
        _phase_role_payload(receipt),
    )
    local_detected = receipt.scope == "local_detected_fixed_composition_role"
    local_anchored = receipt.scope == "local_reference_anchored_fixed_composition_role"
    ideal_gas = receipt.scope == "typed_ideal_gas_fixed_composition_role"
    if (
        receipt.normalized_chemical_mole_fractions != normalized
        or composition_error > 2.0e-12
        or volume_error > 2.0e-8
        or pressure_error > 2.0e-8
        or receipt.mechanical_class != "strict_stable"
        or receipt.phase_identity != descriptor.phase_identity
        or receipt.declared_role != descriptor.phase_role
        or receipt.active_trial_fingerprint != trial.fingerprint
        or receipt.active_parameter_fingerprint != trial.active_parameter_fingerprint
        or receipt.problem_fingerprint != descriptor.problem_fingerprint
        or receipt.reaction_fingerprint != descriptor.reaction_fingerprint
        or receipt.topology_fingerprint != descriptor.topology_fingerprint
        or receipt.globality_claim != "not_established"
        or receipt.fingerprint != expected_fingerprint
        or receipt.provider_evaluation_count < 0
        or receipt.continuation_reference_fingerprint
        != state.continuation_certificate.reference_fingerprint
        or not (local_detected or local_anchored or ideal_gas)
        or (
            local_detected
            and (
                receipt.phase_selection_claim != "mechanically_certified_declared_single_phase"
                or receipt.completeness_criterion
                != "detected_roots_only;pressure_envelope_root_completeness_not_proven"
                or "detected_strict_stable_volume_extreme_for_declared_role"
                not in receipt.role_evidence
                or (
                    "minimum_detected_strict_stable_volume_selected"
                    if receipt.declared_role == "liquid"
                    else "maximum_detected_strict_stable_volume_selected"
                )
                not in receipt.role_evidence
                or "root_completeness_proof" not in receipt.search_exclusions
                or receipt.provider_evaluation_count <= 0
            )
        )
        or (
            local_anchored
            and (
                receipt.phase_selection_claim
                != "reference_anchored_locally_certified_declared_role"
                or receipt.native_certificate_scope
                != "exact_state_local_mechanical_stability_plus_role_certified_anchor"
                or "exact_state_negative_pressure_derivative_log_volume"
                not in receipt.role_evidence
                or "exact_state_matches_reference_selected_chemical_branch"
                not in receipt.role_evidence
                or not any(
                    evidence.startswith("anchor_state_identity:sha256:")
                    for evidence in receipt.role_evidence
                )
                or not any(
                    evidence.startswith("anchor_phase_role_receipt:sha256:")
                    for evidence in receipt.role_evidence
                )
                or receipt.provider_evaluation_count != 1
            )
        )
        or (ideal_gas and receipt.provider_evaluation_count != 0)
    ):
        raise ValueError("homogeneous phase-role receipt is inconsistent")


def _validate_homogeneous_receipts(
    descriptor: HomogeneousReactiveObservationDescriptor,
    trial: ActiveParameterTrialReceipt,
    state: CertifiedHomogeneousState,
    branch: PhysicalBranchReceipt,
) -> None:
    _validate_phase_role_receipt(descriptor, trial, state)
    _validate_chemical_search_certificate(state.chemical_search_certificate)
    _validate_continuation_certificate(descriptor, state.continuation_certificate)
    state_payload = _state_identity_payload(
        descriptor,
        temperature_k=state.temperature_k,
        pressure_pa=state.pressure_pa,
        eos_pressure_pa=state.eos_pressure_pa,
        amounts_mol=state.amounts_mol,
        mole_fractions=state.mole_fractions,
        volume_m3=state.volume_m3,
        molar_density_mol_per_m3=state.molar_density_mol_per_m3,
        chart_topology=state.chart_topology,
        phase_role_certificate_fingerprint=state.phase_role_receipt.fingerprint,
    )
    expected_state_identity = _canonical_fingerprint(
        "epcsaft.homogeneous-reactive-state-point", state_payload
    )
    expected_state_fingerprint = _canonical_fingerprint(
        "epcsaft.certified-homogeneous-state",
        {
            "state_identity": expected_state_identity,
            "trial_fingerprint": trial.fingerprint,
            "chemical_search_certificate_fingerprint": (
                state.chemical_search_certificate.fingerprint
            ),
            "continuation_certificate_fingerprint": (
                state.continuation_certificate.fingerprint
            ),
        },
    )
    expected_branch_identity = _canonical_fingerprint(
        "epcsaft.homogeneous-reactive-physical-branch",
        {
            "phase_identity": descriptor.phase_identity,
            "phase_role": descriptor.phase_role,
            "problem_fingerprint": descriptor.problem_fingerprint,
            "reaction_fingerprint": descriptor.reaction_fingerprint,
            "topology_fingerprint": descriptor.topology_fingerprint,
            "continuation_identity": descriptor.continuation_identity,
            "chart_topology": state.chart_topology,
            "canonical_thermodynamic_state": expected_state_identity,
            "phase_role_certificate_fingerprint": state.phase_role_receipt.fingerprint,
            "branch_family_identity": state.continuation_certificate.branch_family_identity,
        },
    )
    if (
        state.phase_identity != descriptor.phase_identity
        or state.phase_role != descriptor.phase_role
        or state.state_identity != expected_state_identity
        or state.fingerprint != expected_state_fingerprint
        or branch.identity != expected_branch_identity
        or branch.branch_family_identity
        != state.continuation_certificate.branch_family_identity
        or branch.phase_identity != descriptor.phase_identity
        or branch.phase_role != descriptor.phase_role
        or branch.continuation_identity != descriptor.continuation_identity
        or branch.state_identity != state.state_identity
        or branch.amounts_mol != state.amounts_mol
        or branch.volume_m3 != state.volume_m3
        or branch.problem_fingerprint != descriptor.problem_fingerprint
        or state.problem_fingerprint != descriptor.problem_fingerprint
        or branch.reaction_fingerprint != descriptor.reaction_fingerprint
        or state.reaction_fingerprint != descriptor.reaction_fingerprint
        or branch.topology_fingerprint != descriptor.topology_fingerprint
        or state.topology_fingerprint != descriptor.topology_fingerprint
        or branch.phase_role_certificate_fingerprint != state.phase_role_receipt.fingerprint
        or branch.chemical_search_certificate_fingerprint
        != state.chemical_search_certificate.fingerprint
        or branch.continuation_certificate_fingerprint
        != state.continuation_certificate.fingerprint
    ):
        raise ValueError("homogeneous state and physical-branch receipt are inconsistent")


def _receipt_float(value: object, field: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise TypeError(f"Equilibrium owner receipt {field} must be numeric")
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"Equilibrium owner receipt {field} must be finite")
    return result


def _receipt_int(value: object, field: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"Equilibrium owner receipt {field} must be an integer")
    return value


def _receipt_bool(value: object, field: str) -> bool:
    if not isinstance(value, bool):
        raise TypeError(f"Equilibrium owner receipt {field} must be a boolean")
    return value


def _receipt_float_tuple(value: object, field: str) -> tuple[float, ...]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        raise TypeError(f"Equilibrium owner receipt {field} must be a numeric sequence")
    return tuple(_receipt_float(item, field) for item in value)


def _receipt_int_tuple(value: object, field: str) -> tuple[int, ...]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        raise TypeError(f"Equilibrium owner receipt {field} must be an integer sequence")
    return tuple(_receipt_int(item, field) for item in value)


def _receipt_string_tuple(value: object, field: str) -> tuple[str, ...]:
    if (
        not isinstance(value, Sequence)
        or isinstance(value, (str, bytes))
        or not all(isinstance(item, str) for item in value)
    ):
        raise TypeError(f"Equilibrium owner receipt {field} must be a string sequence")
    return tuple(cast(str, item) for item in value)


def _chemical_search_certificate_mapping(
    certificate: HomogeneousChemicalSearchCertificate,
) -> dict[str, object]:
    return {
        **_chemical_search_certificate_payload(certificate),
        "budget_prefix_selection_history": tuple(
            {
                "primary_budget": budget,
                "attempted_primary_ordinals": attempts,
                "basin_ordinals": basins,
                "selected_basin_ordinal": selection,
                "selection_changed": changed,
            }
            for budget, attempts, basins, selection, changed in (
                certificate.budget_prefix_selection_history
            )
        ),
        "fingerprint": certificate.fingerprint,
    }


def _state_receipt_mapping(state: CertifiedHomogeneousState) -> dict[str, object]:
    role = state.phase_role_receipt
    return {
        "phase_identity": state.phase_identity,
        "phase_role": state.phase_role,
        "temperature_k": state.temperature_k,
        "pressure_pa": state.pressure_pa,
        "eos_pressure_pa": state.eos_pressure_pa,
        "amounts_mol": state.amounts_mol,
        "mole_fractions": state.mole_fractions,
        "volume_m3": state.volume_m3,
        "molar_density_mol_per_m3": state.molar_density_mol_per_m3,
        "state_identity": state.state_identity,
        "problem_fingerprint": state.problem_fingerprint,
        "reaction_fingerprint": state.reaction_fingerprint,
        "topology_fingerprint": state.topology_fingerprint,
        "chart_topology": state.chart_topology,
        "chemical_search_certificate": _chemical_search_certificate_mapping(
            state.chemical_search_certificate
        ),
        "continuation_certificate": {
            **_continuation_certificate_payload(state.continuation_certificate),
            "fingerprint": state.continuation_certificate.fingerprint,
        },
        "phase_role_receipt": {**_phase_role_payload(role), "fingerprint": role.fingerprint},
        "fingerprint": state.fingerprint,
    }


def _branch_receipt_mapping(branch: PhysicalBranchReceipt) -> dict[str, object]:
    return {
        "identity": branch.identity,
        "branch_family_identity": branch.branch_family_identity,
        "continuation_identity": branch.continuation_identity,
        "state_identity": branch.state_identity,
        "topology_fingerprint": branch.topology_fingerprint,
        "problem_fingerprint": branch.problem_fingerprint,
        "reaction_fingerprint": branch.reaction_fingerprint,
        "phase_role_certificate_fingerprint": branch.phase_role_certificate_fingerprint,
        "chemical_search_certificate_fingerprint": (branch.chemical_search_certificate_fingerprint),
        "continuation_certificate_fingerprint": branch.continuation_certificate_fingerprint,
        "phase_identity": branch.phase_identity,
        "phase_role": branch.phase_role,
        "amounts_mol": branch.amounts_mol,
        "volume_m3": branch.volume_m3,
    }


def _owner_observation_mapping(
    observation: HomogeneousReactiveOwnerObservation,
) -> dict[str, object]:
    return {
        "target_identity": observation.target_identity,
        "prediction_identity": observation.prediction_identity,
        "aggregate_identity": observation.aggregate_identity,
        "covariance_identity": observation.covariance_identity,
    }


def _requested_state(temperature: Quantity[Any], pressure: Quantity[Any]) -> tuple[float, float]:
    try:
        temperature_k = float(temperature.to("kelvin").magnitude)
        pressure_pa = float(pressure.to("pascal").magnitude)
    except (AttributeError, TypeError, ValueError) as error:
        raise ValueError("Equilibrium owner receipt state has invalid units") from error
    if not math.isfinite(temperature_k) or not math.isfinite(pressure_pa):
        raise ValueError("Equilibrium owner receipt state must be finite")
    return temperature_k, pressure_pa


def serialize_homogeneous_reactive_owner_receipt(
    result: HomogeneousReactiveObservationResult,
    *,
    batch_identity: str,
    requested_temperature: Quantity[Any],
    requested_pressure: Quantity[Any],
    observations: tuple[HomogeneousReactiveOwnerObservation, ...],
) -> str:
    """Serialize one canonical final-state owner receipt inside Equilibrium."""

    selected = tuple(observations)
    if (
        result.status != "evaluated"
        or result.trial is None
        or result.state is None
        or result.physical_branch is None
        or not batch_identity
        or not selected
    ):
        raise ValueError("Equilibrium cannot serialize an incomplete owner receipt")
    temperature_k, pressure_pa = _requested_state(requested_temperature, requested_pressure)
    if result.state.temperature_k != temperature_k or result.state.pressure_pa != pressure_pa:
        raise ValueError("certified state does not match the requested owner batch state")
    rows = {row.identity: row for row in result.descriptor.rows}
    if len({observation.target_identity for observation in selected}) != len(selected):
        raise ValueError("Equilibrium owner receipt target identities must be unique")
    for observation in selected:
        row = rows.get(observation.prediction_identity)
        if (
            not observation.target_identity
            or row is None
            or row.aggregate_identity != observation.aggregate_identity
            or row.covariance_identity != observation.covariance_identity
        ):
            raise ValueError("Equilibrium owner receipt observation selection is inconsistent")
    _validate_homogeneous_receipts(
        result.descriptor,
        result.trial,
        result.state,
        result.physical_branch,
    )
    return _canonical_json_bytes(
        {
            "schema": "epcsaft.equilibrium-owner-batch-receipt",
            "schema_version": 4,
            "batch_identity": batch_identity,
            "requested_temperature_k": temperature_k,
            "requested_pressure_pa": pressure_pa,
            "descriptor_fingerprint": result.descriptor.fingerprint,
            "problem_fingerprint": result.descriptor.problem_fingerprint,
            "reaction_fingerprint": result.descriptor.reaction_fingerprint,
            "trial": {
                "parameter_fingerprint": result.trial.parameter_fingerprint,
                "identities": result.trial.identities,
                "values": result.trial.values,
                "units": result.trial.units,
                "active_parameter_fingerprint": result.trial.active_parameter_fingerprint,
                "fingerprint": result.trial.fingerprint,
            },
            "certified_state": _state_receipt_mapping(result.state),
            "physical_branch": _branch_receipt_mapping(result.physical_branch),
            "observations": tuple(_owner_observation_mapping(item) for item in selected),
        }
    ).decode("utf-8")


def _parse_chemical_search_certificate(
    value: object,
) -> HomogeneousChemicalSearchCertificate:
    if not isinstance(value, Mapping):
        raise TypeError("Equilibrium owner receipt chemical search certificate is missing")
    raw_history = value.get("budget_prefix_selection_history")
    if not isinstance(raw_history, Sequence) or isinstance(raw_history, (str, bytes)):
        raise TypeError("Equilibrium owner receipt search history must be a sequence")
    history = []
    for raw_record in raw_history:
        if not isinstance(raw_record, Mapping):
            raise TypeError("Equilibrium owner receipt search history record is invalid")
        raw_selection = raw_record.get("selected_basin_ordinal")
        history.append(
            (
                _receipt_int(raw_record.get("primary_budget"), "search prefix budget"),
                _receipt_int_tuple(
                    raw_record.get("attempted_primary_ordinals"),
                    "search attempted ordinals",
                ),
                _receipt_int_tuple(raw_record.get("basin_ordinals"), "search basin ordinals"),
                (
                    None
                    if raw_selection is None
                    else _receipt_int(raw_selection, "search prefix selection")
                ),
                _receipt_bool(raw_record.get("selection_changed"), "selection changed"),
            )
        )
    certificate = HomogeneousChemicalSearchCertificate(
        status=str(value.get("status")),
        selection_label=str(value.get("selection_label")),
        primary_budget=_receipt_int(value.get("primary_budget"), "primary budget"),
        primary_attempt_count=_receipt_int(
            value.get("primary_attempt_count"), "primary attempt count"
        ),
        generated_start_count=_receipt_int(
            value.get("generated_start_count"), "generated start count"
        ),
        budget_truncated_start_count=_receipt_int(
            value.get("budget_truncated_start_count"), "truncated start count"
        ),
        duplicate_start_count=_receipt_int(
            value.get("duplicate_start_count"), "duplicate start count"
        ),
        infeasible_start_count=_receipt_int(
            value.get("infeasible_start_count"), "infeasible start count"
        ),
        evaluated_start_count=_receipt_int(
            value.get("evaluated_start_count"), "evaluated start count"
        ),
        domain_rejected_start_count=_receipt_int(
            value.get("domain_rejected_start_count"), "domain rejected start count"
        ),
        construction_rejected_start_count=_receipt_int(
            value.get("construction_rejected_start_count"),
            "construction rejected start count",
        ),
        attempt_count=_receipt_int(value.get("attempt_count"), "attempt count"),
        basin_count=_receipt_int(value.get("basin_count"), "basin count"),
        selected_basin_ordinal=_receipt_int(
            value.get("selected_basin_ordinal"), "selected basin ordinal"
        ),
        selected_basin_representative_attempt_ordinal=_receipt_int(
            value.get("selected_basin_representative_attempt_ordinal"),
            "selected representative attempt ordinal",
        ),
        budget_prefix_selection_history=tuple(history),
        fingerprint=str(value.get("fingerprint")),
    )
    _validate_chemical_search_certificate(certificate)
    return certificate


def _parse_continuation_certificate(
    value: object,
    descriptor: HomogeneousReactiveObservationDescriptor,
) -> HomogeneousContinuationCertificate:
    if not isinstance(value, Mapping):
        raise TypeError("Equilibrium owner receipt continuation certificate must be an object")
    certificate = HomogeneousContinuationCertificate(
        reference_fingerprint=str(value.get("reference_fingerprint")),
        branch_family_identity=str(value.get("branch_family_identity")),
        search_basin_count=_receipt_int(value.get("search_basin_count"), "search basin count"),
        in_region_basin_count=_receipt_int(
            value.get("in_region_basin_count"), "in-region basin count"
        ),
        selected_basin_ordinal=_receipt_int(
            value.get("selected_basin_ordinal"), "continuation selected basin ordinal"
        ),
        selected_composition_distance=_receipt_float(
            value.get("selected_composition_distance"), "selected composition distance"
        ),
        selected_volume_distance=_receipt_float(
            value.get("selected_volume_distance"), "selected volume distance"
        ),
        fingerprint=str(value.get("fingerprint")),
    )
    _validate_continuation_certificate(descriptor, certificate)
    return certificate


def validate_homogeneous_reactive_owner_receipt(
    receipt_json: str,
    *,
    descriptor: HomogeneousReactiveObservationDescriptor,
    batch_identity: str,
    requested_temperature: Quantity[Any],
    requested_pressure: Quantity[Any],
    active_parameters: epcsaft.ActiveParameterSet,
    observations: tuple[HomogeneousReactiveOwnerObservation, ...],
) -> str:
    """Validate opaque canonical owner bytes and return their SHA-256 identity."""

    if not receipt_json:
        raise ValueError("usable Equilibrium owner receipt must be nonempty")
    receipt = json.loads(receipt_json)
    if (
        not isinstance(receipt, Mapping)
        or _canonical_json_bytes(receipt).decode("utf-8") != receipt_json
    ):
        raise ValueError("Equilibrium owner receipt is not canonical JSON")
    trial = receipt.get("trial")
    state = receipt.get("certified_state")
    branch = receipt.get("physical_branch")
    role = state.get("phase_role_receipt") if isinstance(state, Mapping) else None
    search = state.get("chemical_search_certificate") if isinstance(state, Mapping) else None
    continuation = state.get("continuation_certificate") if isinstance(state, Mapping) else None
    raw_observations = receipt.get("observations")
    if (
        not isinstance(trial, Mapping)
        or not isinstance(state, Mapping)
        or not isinstance(branch, Mapping)
        or not isinstance(role, Mapping)
        or not isinstance(raw_observations, Sequence)
        or isinstance(raw_observations, (str, bytes))
    ):
        raise TypeError("Equilibrium owner receipt is incomplete")
    trial_receipt = ActiveParameterTrialReceipt(
        parameter_fingerprint=str(trial.get("parameter_fingerprint")),
        identities=_receipt_string_tuple(trial.get("identities"), "trial identities"),
        values=_receipt_float_tuple(trial.get("values"), "trial values"),
        units=_receipt_string_tuple(trial.get("units"), "trial units"),
        active_parameter_fingerprint=str(trial.get("active_parameter_fingerprint")),
        fingerprint=str(trial.get("fingerprint")),
    )
    role_receipt = HomogeneousPhaseRoleReceipt(
        scope=str(role.get("scope")),
        phase_identity=str(role.get("phase_identity")),
        declared_role=cast(Any, role.get("declared_role")),
        normalized_chemical_mole_fractions=_receipt_float_tuple(
            role.get("normalized_chemical_mole_fractions"), "normalized chemical composition"
        ),
        audited_mole_fractions=_receipt_float_tuple(
            role.get("audited_mole_fractions"), "audited composition"
        ),
        audited_molar_volume_m3_per_mol=_receipt_float(
            role.get("audited_molar_volume_m3_per_mol"), "audited molar volume"
        ),
        audited_pressure_pa=_receipt_float(role.get("audited_pressure_pa"), "audited pressure"),
        mechanical_class=str(role.get("mechanical_class")),
        route=str(role.get("route")),
        native_certificate_scope=str(role.get("native_certificate_scope")),
        phase_selection_claim=str(role.get("phase_selection_claim")),
        globality_claim=str(role.get("globality_claim")),
        search_scope=str(role.get("search_scope")),
        search_budget=str(role.get("search_budget")),
        search_exclusions=_receipt_string_tuple(
            role.get("search_exclusions"), "phase-role search exclusions"
        ),
        completeness_criterion=str(role.get("completeness_criterion")),
        role_evidence=_receipt_string_tuple(role.get("role_evidence"), "phase-role evidence"),
        provider_evaluation_count=_receipt_int(
            role.get("provider_evaluation_count"), "Provider evaluation count"
        ),
        continuation_reference_fingerprint=str(
            role.get("continuation_reference_fingerprint")
        ),
        active_trial_fingerprint=str(role.get("active_trial_fingerprint")),
        active_parameter_fingerprint=str(role.get("active_parameter_fingerprint")),
        problem_fingerprint=str(role.get("problem_fingerprint")),
        reaction_fingerprint=str(role.get("reaction_fingerprint")),
        topology_fingerprint=str(role.get("topology_fingerprint")),
        fingerprint=str(role.get("fingerprint")),
    )
    state_receipt = CertifiedHomogeneousState(
        phase_identity=str(state.get("phase_identity")),
        phase_role=cast(Any, state.get("phase_role")),
        temperature_k=_receipt_float(state.get("temperature_k"), "state temperature"),
        pressure_pa=_receipt_float(state.get("pressure_pa"), "state pressure"),
        eos_pressure_pa=_receipt_float(state.get("eos_pressure_pa"), "EOS pressure"),
        amounts_mol=_receipt_float_tuple(state.get("amounts_mol"), "state amounts"),
        mole_fractions=_receipt_float_tuple(state.get("mole_fractions"), "state mole fractions"),
        volume_m3=_receipt_float(state.get("volume_m3"), "state volume"),
        molar_density_mol_per_m3=_receipt_float(
            state.get("molar_density_mol_per_m3"), "state molar density"
        ),
        state_identity=str(state.get("state_identity")),
        problem_fingerprint=str(state.get("problem_fingerprint")),
        reaction_fingerprint=str(state.get("reaction_fingerprint")),
        topology_fingerprint=str(state.get("topology_fingerprint")),
        chart_topology=str(state.get("chart_topology")),
        chemical_search_certificate=_parse_chemical_search_certificate(search),
        continuation_certificate=_parse_continuation_certificate(continuation, descriptor),
        phase_role_receipt=role_receipt,
        fingerprint=str(state.get("fingerprint")),
    )
    branch_receipt = PhysicalBranchReceipt(
        identity=str(branch.get("identity")),
        branch_family_identity=str(branch.get("branch_family_identity")),
        continuation_identity=str(branch.get("continuation_identity")),
        state_identity=str(branch.get("state_identity")),
        topology_fingerprint=str(branch.get("topology_fingerprint")),
        problem_fingerprint=str(branch.get("problem_fingerprint")),
        reaction_fingerprint=str(branch.get("reaction_fingerprint")),
        phase_role_certificate_fingerprint=str(branch.get("phase_role_certificate_fingerprint")),
        chemical_search_certificate_fingerprint=str(
            branch.get("chemical_search_certificate_fingerprint")
        ),
        continuation_certificate_fingerprint=str(
            branch.get("continuation_certificate_fingerprint")
        ),
        phase_identity=str(branch.get("phase_identity")),
        phase_role=cast(Any, branch.get("phase_role")),
        amounts_mol=_receipt_float_tuple(branch.get("amounts_mol"), "branch amounts"),
        volume_m3=_receipt_float(branch.get("volume_m3"), "branch volume"),
    )
    temperature_k, pressure_pa = _requested_state(requested_temperature, requested_pressure)
    expected_trial = _trial_receipt(descriptor, active_parameters)
    expected_observations = tuple(_owner_observation_mapping(item) for item in observations)
    if (
        receipt.get("schema") != "epcsaft.equilibrium-owner-batch-receipt"
        or receipt.get("schema_version") != 4
        or receipt.get("batch_identity") != batch_identity
        or receipt.get("requested_temperature_k") != temperature_k
        or receipt.get("requested_pressure_pa") != pressure_pa
        or receipt.get("descriptor_fingerprint") != descriptor.fingerprint
        or receipt.get("problem_fingerprint") != descriptor.problem_fingerprint
        or receipt.get("reaction_fingerprint") != descriptor.reaction_fingerprint
        or trial_receipt != expected_trial
        or tuple(raw_observations) != expected_observations
        or state_receipt.temperature_k != temperature_k
        or state_receipt.pressure_pa != pressure_pa
    ):
        raise ValueError("Equilibrium owner receipt does not bind the accepted batch trial")
    _validate_homogeneous_receipts(
        descriptor,
        trial_receipt,
        state_receipt,
        branch_receipt,
    )
    return f"sha256:{hashlib.sha256(receipt_json.encode('utf-8')).hexdigest()}"


class _ObservationFailure(ValueError):
    def __init__(self, code: NonEvaluableTrialCode, diagnostic: str) -> None:
        super().__init__(diagnostic)
        self.code = code


def evaluate_homogeneous_reactive_observations(
    problem: HomogeneousReactiveObservationProblem,
    temperature: Quantity[Any],
    pressure: Quantity[Any],
    descriptor: HomogeneousReactiveObservationDescriptor,
    *,
    active_parameters: epcsaft.ActiveParameterSet | None = None,
) -> HomogeneousReactiveObservationResult:
    """Evaluate one descriptor without constructing residuals or failure penalties."""

    if not isinstance(descriptor, HomogeneousReactiveObservationDescriptor):
        raise TypeError("descriptor must be a homogeneous reactive descriptor")
    if not isinstance(problem, HomogeneousReactiveObservationProblem):
        raise TypeError("problem must be a homogeneous reactive observation problem")
    phase = problem.phase
    reaction_system = problem.reaction_system
    try:
        expected = homogeneous_reactive_observation_descriptor(
            problem,
            descriptor.rows,
            active_parameters=active_parameters,
        )
    except (TypeError, ValueError) as error:
        code: NonEvaluableTrialCode = (
            "component_order_mismatch"
            if "component order" in str(error)
            else "parameter_topology_mismatch"
        )
        return _failed_result(descriptor, code, str(error))
    trial = _trial_receipt(descriptor, active_parameters)
    if (
        expected.fingerprint != descriptor.fingerprint
        or expected.active_parameter_selection_fingerprint
        != descriptor.active_parameter_selection_fingerprint
    ):
        return _failed_result(
            descriptor,
            "parameter_topology_mismatch",
            "descriptor identity does not match the trial problem",
            trial=trial,
        )
    try:
        equilibrium = _chemical_equilibrium_impl(
            phase,
            temperature,
            pressure,
            reaction_system,
            sensitivity_request=ChemicalEquilibriumSensitivityRequest(active_parameters),
            continuation_start=_immutable_continuation_start(problem),
        )
    except ChemicalEquilibriumError as error:
        return _failed_result(
            descriptor,
            _chemical_failure_code(error.diagnostics),
            str(error),
            trial=trial,
            diagnostics=error.diagnostics,
            support_domain_status=error.diagnostics.provider_domain_status,
        )
    except ValueError as error:
        return _failed_result(
            descriptor,
            "branch_loss",
            str(error),
            trial=trial,
        )
    try:
        continuation_certificate = _certify_continuation(
            problem,
            descriptor,
            equilibrium,
        )
    except _ObservationFailure as error:
        return _failed_result(
            descriptor,
            error.code,
            str(error),
            trial=trial,
            equilibrium=equilibrium,
            diagnostics=equilibrium.diagnostics,
            support_domain_status=equilibrium.diagnostics.provider_domain_status,
        )
    sensitivity = equilibrium.sensitivity
    if sensitivity is None or sensitivity.status != "available":
        reason = "chemical equilibrium did not return an exact conditioned sensitivity"
        if sensitivity is not None and sensitivity.failure_reason:
            reason = sensitivity.failure_reason
        return _failed_result(
            descriptor,
            _sensitivity_failure_code(sensitivity),
            reason,
            trial=trial,
            equilibrium=equilibrium,
            diagnostics=equilibrium.diagnostics,
            support_domain_status=equilibrium.diagnostics.provider_domain_status,
        )
    active_ids = descriptor.active_parameter_identities
    parameter_rows = {
        parameter.identity: (amounts, volume)
        for parameter, amounts, volume in zip(
            sensitivity.parameters,
            sensitivity.amount_derivatives,
            sensitivity.volume_derivatives,
            strict=True,
        )
        if parameter.identity is not None
    }
    if tuple(parameter_rows) != active_ids:
        return _failed_result(
            descriptor,
            "parameter_topology_mismatch",
            "chemical sensitivity did not preserve active parameter order",
            trial=trial,
            equilibrium=equilibrium,
            diagnostics=equilibrium.diagnostics,
            support_domain_status=equilibrium.diagnostics.provider_domain_status,
        )
    amounts = equilibrium.amounts_mol
    volume = equilibrium.volume_m3
    total_amount = math.fsum(amounts)
    parameter_count = len(active_ids)
    amount_derivatives = tuple(parameter_rows[identity][0] for identity in active_ids)
    volume_derivatives = tuple(parameter_rows[identity][1] for identity in active_ids)
    total_amount_derivatives = tuple(math.fsum(row) for row in amount_derivatives)
    mole_fraction_derivatives = tuple(
        tuple(
            (
                amount_derivatives[parameter][component] * total_amount
                - amounts[component] * total_amount_derivatives[parameter]
            )
            / (total_amount * total_amount)
            for parameter in range(parameter_count)
        )
        for component in range(len(amounts))
    )
    potential_values: tuple[float, ...]
    potential_derivatives: tuple[tuple[float, ...], ...]
    provider_phase = phase if isinstance(phase, ProviderPhase) else None
    if provider_phase is not None:
        pressure_derivative_log_volume_pa: float | None = None
        try:
            if active_parameters is None:
                fixed_state = provider_phase.model._state(
                    T=equilibrium.temperature_k * epcsaft.unit_registry.kelvin,
                    rho=(
                        total_amount
                        / volume
                        * epcsaft.unit_registry.mole
                        / epcsaft.unit_registry.meter**3
                    ),
                    x=equilibrium.mole_fractions,
                )
                potential_values = tuple(
                    math.log(amount / volume) + residual
                    for amount, residual in zip(
                        amounts,
                        fixed_state.residual_chemical_potential_over_rt,
                        strict=True,
                    )
                )
                potential_derivatives = tuple(() for _ in amounts)
                eos_pressure_pa = float(fixed_state.pressure.to("pascal").magnitude)
            else:
                coordinates = provider_phase.model._active_coordinates(active_parameters)
                block = provider_phase.model._native.evaluate_phase_derivatives(
                    equilibrium.temperature_k,
                    amounts,
                    volume,
                    coordinates,
                    2,
                )
                width = len(amounts) + 1 + parameter_count
                hessian = tuple(float(value) for value in block.hessian)
                pressure_derivative_log_volume_pa = (
                    -8.31446261815324
                    * equilibrium.temperature_k
                    * volume
                    * hessian[len(amounts) * width + len(amounts)]
                )
                potential_values = tuple(float(value) for value in block.gradient[: len(amounts)])
                potential_derivatives = tuple(
                    tuple(
                        math.fsum(
                            hessian[component * width + species]
                            * amount_derivatives[parameter][species]
                            for species in range(len(amounts))
                        )
                        + hessian[component * width + len(amounts)] * volume_derivatives[parameter]
                        + hessian[component * width + len(amounts) + 1 + parameter]
                        for parameter in range(parameter_count)
                    )
                    for component in range(len(amounts))
                )
                eos_pressure_pa = float(block.pressure_pa)
        except (RuntimeError, ValueError) as error:
            return _failed_result(
                descriptor,
                "provider_domain_rejection",
                str(error),
                trial=trial,
                equilibrium=equilibrium,
                diagnostics=equilibrium.diagnostics,
                support_domain_status=equilibrium.diagnostics.provider_domain_status,
            )
    else:
        pressure_derivative_log_volume_pa = None
        eos_pressure_pa = total_amount * 8.31446261815324 * equilibrium.temperature_k / volume
        potential_values = tuple(math.log(amount / volume) for amount in amounts)
        potential_derivatives = tuple(() for _ in amounts)
    pressure_error = abs(eos_pressure_pa - equilibrium.pressure_pa) / max(
        1.0, abs(equilibrium.pressure_pa)
    )
    if not math.isfinite(eos_pressure_pa) or pressure_error > 2.0e-8:
        return _failed_result(
            descriptor,
            "physical_rejection",
            "exact certified-state EOS pressure does not close the declared pressure",
            trial=trial,
            equilibrium=equilibrium,
            diagnostics=equilibrium.diagnostics,
            support_domain_status=equilibrium.diagnostics.provider_domain_status,
        )
    try:
        role_receipt = _certify_phase_role(
            problem,
            descriptor,
            equilibrium,
            active_parameters,
            trial,
            pressure_derivative_log_volume_pa,
            continuation_certificate.reference_fingerprint,
        )
        state, physical_branch = _make_certified_receipts(
            descriptor,
            trial,
            equilibrium,
            eos_pressure_pa,
            role_receipt,
            continuation_certificate,
        )
    except PhaseEquilibriumError as error:
        return _failed_result(
            descriptor,
            "physical_rejection",
            f"declared phase-role certification failed: {error}",
            trial=trial,
            equilibrium=equilibrium,
            diagnostics=equilibrium.diagnostics,
            support_domain_status=equilibrium.diagnostics.provider_domain_status,
        )
    except (RuntimeError, TypeError, ValueError) as error:
        return _failed_result(
            descriptor,
            "physical_rejection",
            f"declared phase-role receipt is invalid: {error}",
            trial=trial,
            equilibrium=equilibrium,
            diagnostics=equilibrium.diagnostics,
            support_domain_status=equilibrium.diagnostics.provider_domain_status,
        )
    row_results: list[HomogeneousReactiveObservationRowResult] = []
    row_failure_codes: list[NonEvaluableTrialCode] = []
    for row in descriptor.rows:
        try:
            coefficients = row.coefficients
            if row.quantity == "amount_mol":
                value = math.fsum(c * n for c, n in zip(coefficients, amounts, strict=True))
                jacobian = tuple(
                    math.fsum(c * dn for c, dn in zip(coefficients, derivative, strict=True))
                    for derivative in amount_derivatives
                )
            elif row.quantity == "mole_fraction":
                value = math.fsum(
                    c * x for c, x in zip(coefficients, equilibrium.mole_fractions, strict=True)
                )
                jacobian = tuple(
                    math.fsum(
                        coefficients[component] * mole_fraction_derivatives[component][parameter]
                        for component in range(len(amounts))
                    )
                    for parameter in range(parameter_count)
                )
            elif row.quantity == "molarity_mol_per_m3":
                numerator = math.fsum(c * n for c, n in zip(coefficients, amounts, strict=True))
                value = numerator / volume
                jacobian = tuple(
                    math.fsum(c * dn for c, dn in zip(coefficients, derivative, strict=True))
                    / volume
                    - numerator * dvolume / (volume * volume)
                    for derivative, dvolume in zip(
                        amount_derivatives, volume_derivatives, strict=True
                    )
                )
            elif row.quantity == "molality_mol_per_kg":
                numerator = math.fsum(c * n for c, n in zip(coefficients, amounts, strict=True))
                denominator = math.fsum(
                    c * n
                    for c, n in zip(
                        row.solvent_mass_coefficients_kg_per_mol,
                        amounts,
                        strict=True,
                    )
                )
                if denominator <= 0.0:
                    raise _ObservationFailure(
                        "observation_out_of_support",
                        "molality solvent-mass denominator is not positive",
                    )
                value = numerator / denominator
                jacobian = tuple(
                    (
                        math.fsum(c * dn for c, dn in zip(coefficients, derivative, strict=True))
                        * denominator
                        - numerator
                        * math.fsum(
                            c * dn
                            for c, dn in zip(
                                row.solvent_mass_coefficients_kg_per_mol,
                                derivative,
                                strict=True,
                            )
                        )
                    )
                    / (denominator * denominator)
                    for derivative in amount_derivatives
                )
            elif row.quantity == "molar_density_mol_per_m3":
                value = total_amount / volume
                jacobian = tuple(
                    dtotal / volume - total_amount * dvolume / (volume * volume)
                    for dtotal, dvolume in zip(
                        total_amount_derivatives, volume_derivatives, strict=True
                    )
                )
            elif row.quantity == "compressibility_factor":
                value = (
                    equilibrium.pressure_pa
                    * volume
                    / (total_amount * 8.31446261815324 * equilibrium.temperature_k)
                )
                jacobian = tuple(
                    value * (dvolume / volume - dtotal / total_amount)
                    for dvolume, dtotal in zip(
                        volume_derivatives, total_amount_derivatives, strict=True
                    )
                )
            elif row.quantity == "pressure_pa":
                value = equilibrium.pressure_pa
                jacobian = tuple(0.0 for _ in active_ids)
            else:
                value = math.fsum(
                    c * potential
                    for c, potential in zip(coefficients, potential_values, strict=True)
                )
                jacobian = tuple(
                    math.fsum(
                        coefficients[component] * potential_derivatives[component][parameter]
                        for component in range(len(amounts))
                    )
                    for parameter in range(parameter_count)
                )
                if row.quantity in {"log_fugacity_contraction", "fugacity_pa"}:
                    value += math.fsum(coefficients) * math.log(
                        8.31446261815324 * equilibrium.temperature_k
                    )
                    if row.quantity == "fugacity_pa":
                        value = math.exp(value)
                        jacobian = tuple(value * derivative for derivative in jacobian)
                elif row.quantity == "log_activity_contraction":
                    transfer = equilibrium.source_reference_transfer
                    if transfer is None:
                        raise _ObservationFailure(
                            "source_reference_failure",
                            "source-reference transfer is unavailable",
                        )
                    if parameter_count > 0 and "parameters" not in transfer.derivative_availability:
                        raise _ObservationFailure(
                            "source_reference_derivative_unavailable",
                            "source-reference parameter derivative is unavailable",
                        )
                    try:
                        reference_coordinates = _neutral_basis_coordinates(
                            coefficients, transfer.neutral_basis
                        )
                    except ValueError as error:
                        raise _ObservationFailure("source_reference_failure", str(error)) from error
                    value -= math.fsum(
                        coordinate * reference
                        for coordinate, reference in zip(
                            reference_coordinates,
                            transfer.transfer_log_contractions,
                            strict=True,
                        )
                    )
                    jacobian = tuple(
                        derivative
                        - math.fsum(
                            reference_coordinates[basis]
                            * transfer.parameter_derivatives[basis][parameter]
                            for basis in range(len(reference_coordinates))
                        )
                        for parameter, derivative in enumerate(jacobian)
                    )
            if not math.isfinite(value) or not all(math.isfinite(item) for item in jacobian):
                raise _ObservationFailure(
                    "unsupported_observable",
                    "reactive observation value or Jacobian is nonfinite",
                )
            if row.support == "positive" and value <= 0.0:
                raise _ObservationFailure(
                    "observation_out_of_support",
                    "positive observation left its support",
                )
            if row.support == "nonnegative" and value < 0.0:
                raise _ObservationFailure(
                    "observation_out_of_support",
                    "nonnegative observation left its support",
                )
            row_results.append(
                HomogeneousReactiveObservationRowResult(
                    row.identity, "evaluated", value, jacobian, ""
                )
            )
        except _ObservationFailure as error:
            row_failure_codes.append(error.code)
            row_results.append(
                HomogeneousReactiveObservationRowResult(
                    row.identity, "non_evaluable", None, None, str(error)
                )
            )
        except (ArithmeticError, ValueError) as error:
            row_failure_codes.append("unsupported_observable")
            row_results.append(
                HomogeneousReactiveObservationRowResult(
                    row.identity, "non_evaluable", None, None, str(error)
                )
            )
    if any(row.status == "non_evaluable" for row in row_results):
        reason = next(row.reason for row in row_results if row.status == "non_evaluable")
        return HomogeneousReactiveObservationResult(
            descriptor=descriptor,
            status="non_evaluable",
            failure=NonEvaluableTrial(row_failure_codes[0], reason),
            rows=tuple(row_results),
            trial=trial,
            state=state,
            physical_branch=physical_branch,
            support_domain_status=equilibrium.diagnostics.provider_domain_status,
            equilibrium=equilibrium,
            diagnostics=equilibrium.diagnostics,
        )
    return HomogeneousReactiveObservationResult(
        descriptor=descriptor,
        status="evaluated",
        failure=None,
        rows=tuple(row_results),
        trial=trial,
        state=state,
        physical_branch=physical_branch,
        support_domain_status=equilibrium.diagnostics.provider_domain_status,
        equilibrium=equilibrium,
        diagnostics=equilibrium.diagnostics,
    )


__all__ = (
    "ActiveParameterTrialReceipt",
    "CertifiedHomogeneousState",
    "HomogeneousBranchPolicy",
    "HomogeneousChemicalSearchCertificate",
    "HomogeneousContinuationCertificate",
    "HomogeneousContinuationReference",
    "HomogeneousPhaseRole",
    "HomogeneousPhaseRoleReceipt",
    "HomogeneousReactiveObservationDescriptor",
    "HomogeneousReactiveObservationProblem",
    "HomogeneousReactiveObservationResult",
    "HomogeneousReactiveObservationRow",
    "HomogeneousReactiveObservationRowResult",
    "HomogeneousReactiveOwnerObservation",
    "PhysicalBranchReceipt",
    "certify_homogeneous_continuation_reference",
    "evaluate_homogeneous_reactive_observations",
    "homogeneous_reactive_observation_descriptor",
    "serialize_homogeneous_reactive_owner_receipt",
    "validate_homogeneous_reactive_owner_receipt",
)
