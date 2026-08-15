from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass
from typing import Any, Literal, cast

from pint import Quantity

import epcsaft

from .._trials import NonEvaluableTrial
from ._api import (
    AllComponents,
    ChemicalEquilibriumProblem,
    DeclaredNeutralComponents,
    EquilibriumPhaseState,
    IdealGasPhase,
    NeutralComponents,
    PhaseSupport,
    ProviderPhase,
)
from ._observations import (
    ActiveParameterTrialReceipt,
    HomogeneousBranchPolicy,
    HomogeneousContinuationReference,
    HomogeneousReactiveObservationDescriptor,
    HomogeneousReactiveObservationProblem,
    HomogeneousReactiveObservationRow,
    HomogeneousReactiveOwnerObservation,
    evaluate_homogeneous_reactive_observations,
    homogeneous_reactive_observation_descriptor,
    serialize_homogeneous_reactive_owner_receipt,
    validate_homogeneous_reactive_owner_receipt,
)
from ._reactive_vle import (
    IdealGasVaporComparison,
    ProviderNonidealVapor,
    ReactiveBubbleContinuationReference,
    ReactiveBubbleObservationDescriptor,
    ReactiveBubbleObservationRow,
    ReactiveBubbleVLEProblem,
    certify_reactive_bubble_continuation_reference,
    evaluate_reactive_bubble_observations,
    reactive_bubble_observation_descriptor,
    serialize_reactive_bubble_owner_receipt,
    validate_reactive_bubble_owner_receipt,
)

type AmountRole = Literal["finite", "incipient"]
type OutputSupport = Literal[
    "real",
    "positive",
    "nonnegative",
    "left_censored",
    "right_censored",
]
type EquilibriumSelector = Literal[
    "system.pressure",
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

_EQUILIBRIUM_SELECTORS = frozenset(
    {
        "system.pressure",
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
class ReactivePhase:
    """One declared fluid phase, amount role, support, model, and local anchor."""

    identity: str
    fluid_role: Literal["liquid", "vapor"]
    amount_role: AmountRole
    support: PhaseSupport
    model: ReactivePhaseModel
    continuation_identity: str | None = None
    branch_policy: HomogeneousBranchPolicy | None = None
    continuation_reference: HomogeneousContinuationReference | None = None

    def __post_init__(self) -> None:
        if (
            not self.identity
            or self.fluid_role not in {"liquid", "vapor"}
            or self.amount_role not in {"finite", "incipient"}
            or not isinstance(
                self.support, (AllComponents, NeutralComponents, DeclaredNeutralComponents)
            )
            or not isinstance(self.model, (ProviderModel, IdealGasModel))
        ):
            raise ValueError("reactive phase declaration is invalid")
        if self.amount_role == "finite" and (
            not self.continuation_identity
            or self.branch_policy != "local_certified_role_selected_state"
        ):
            raise ValueError("finite reactive phases require a local branch declaration")
        if self.amount_role == "incipient" and any(
            value is not None
            for value in (
                self.continuation_identity,
                self.branch_policy,
                self.continuation_reference,
            )
        ):
            raise ValueError("incipient phase continuation belongs to the coupled request")
        if isinstance(self.model, IdealGasModel) and self.fluid_role != "vapor":
            raise ValueError("an ideal-gas phase must have the vapor role")


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
        if self.selector == "system.pressure":
            if self.phase_identity is not None or self.coefficients:
                raise ValueError("system pressure has no phase or contraction coefficients")
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
class GeneralReactiveEquilibriumProblem:
    """One role-declared GREPE condition over the current admitted formulations."""

    identity: str
    temperature: Fixed | Solved
    pressure: Fixed | Solved
    phases: tuple[ReactivePhase, ...]
    reaction_system: ChemicalEquilibriumProblem
    reaction_phase_ids: tuple[str, ...]
    outputs: tuple[EquilibriumOutput, ...]
    continuation_identity: str | None = None
    continuation_reference: ReactiveBubbleContinuationReference | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "phases", tuple(self.phases))
        object.__setattr__(self, "reaction_phase_ids", tuple(self.reaction_phase_ids))
        object.__setattr__(self, "outputs", tuple(self.outputs))
        phase_ids = tuple(phase.identity for phase in self.phases)
        output_ids = tuple(output.identity for output in self.outputs)
        if (
            not self.identity
            or not self.phases
            or len(set(phase_ids)) != len(phase_ids)
            or not self.outputs
            or len(set(output_ids)) != len(output_ids)
            or not self.reaction_phase_ids
            or len(set(self.reaction_phase_ids)) != len(self.reaction_phase_ids)
            or any(identity not in phase_ids for identity in self.reaction_phase_ids)
            or not isinstance(self.reaction_system, ChemicalEquilibriumProblem)
        ):
            raise ValueError("general reactive equilibrium request is incomplete")
        finite_ids = {
            phase.identity for phase in self.phases if phase.amount_role == "finite"
        }
        if any(identity not in finite_ids for identity in self.reaction_phase_ids):
            raise ValueError("reactions may be declared only in finite phases")
        for output in self.outputs:
            if output.phase_identity is not None and output.phase_identity not in phase_ids:
                raise ValueError("equilibrium output names an undeclared phase")
        incipient_count = sum(phase.amount_role == "incipient" for phase in self.phases)
        if incipient_count:
            if not self.continuation_identity:
                raise ValueError("an incipient formulation requires a coupled continuation ID")
        elif self.continuation_identity is not None or self.continuation_reference is not None:
            raise ValueError("single finite-phase requests have no coupled continuation")


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
    _owner_kind: Literal["homogeneous", "reactive_bubble"]
    _owner_descriptor: (
        HomogeneousReactiveObservationDescriptor | ReactiveBubbleObservationDescriptor
    )


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
    """One generic result over a locally certified admitted GREPE formulation."""

    problem: GeneralReactiveEquilibriumProblem
    descriptor: GeneralReactiveEquilibriumDescriptor | None
    status: Literal["evaluated", "non_evaluable"]
    rows: tuple[GeneralEquilibriumRowResult, ...]
    phases: tuple[EquilibriumPhaseState, ...]
    failure: NonEvaluableTrial | None
    certificate_scope: str
    globality_claim: Literal["not_established"]
    trial: ActiveParameterTrialReceipt | None
    branch_signature: tuple[object, ...] | None
    owner_result: object | None

    @property
    def values(self) -> tuple[float | None, ...]:
        return tuple(row.value for row in self.rows)

    @property
    def jacobian(self) -> tuple[tuple[float, ...] | None, ...]:
        return tuple(row.jacobian for row in self.rows)


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
    reaction = problem.reaction_system
    component_ids = model.component_ids
    if reaction.species_ids != component_ids:
        raise ValueError("GREPE species order must equal the Provider component order")
    component_count = len(component_ids)
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
    if standard is not None:
        eligible = {
            phase.identity: phase for phase in problem.phases if phase.identity in problem.reaction_phase_ids
        }
        if standard.source_reference_phase != "liquid" or any(
            phase.fluid_role != "liquid" for phase in eligible.values()
        ):
            raise ValueError("reaction standard state is inconsistent with its eligible phase")


def _phase_model(
    model: epcsaft.Mixture,
    phase: ReactivePhase,
) -> ProviderPhase | IdealGasPhase:
    if isinstance(phase.model, ProviderModel):
        if phase.model.reference_id != "installed-provider-eos":
            raise ValueError("GREPE provider model reference is not registered")
        return ProviderPhase(model, phase.model.admissible_packing_fraction_interval)
    return IdealGasPhase(model.parameter_fingerprint, phase.model.reference_id)


def _homogeneous_quantity(output: EquilibriumOutput) -> str:
    quantities = {
        "phase.amount": "amount_mol",
        "phase.mole_fraction": "mole_fraction",
        "phase.molality": "molality_mol_per_kg",
        "phase.molar_density": "molar_density_mol_per_m3",
        "phase.compressibility_factor": "compressibility_factor",
        "phase.pressure": "pressure_pa",
        "system.pressure": "pressure_pa",
        "phase.fugacity": "fugacity_pa",
        "phase.log_activity": "log_activity_contraction",
        "phase.log_fugacity": "log_fugacity_contraction",
        "phase.chemical_potential_over_rt": "chemical_potential_over_rt_contraction",
    }
    try:
        return quantities[output.selector]
    except KeyError as error:
        raise ValueError("output selector is unavailable for one finite phase") from error


def _bubble_quantity(output: EquilibriumOutput, finite: ReactivePhase, incipient: ReactivePhase) -> str:
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
                else "" if cast(Any, row).failure is None else cast(Any, row).failure.diagnostic
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
    if isinstance(phase.support, DeclaredNeutralComponents):
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
        "continuation_identity": phase.continuation_identity,
        "branch_policy": phase.branch_policy,
    }


def _descriptor(
    model: epcsaft.Mixture,
    problem: GeneralReactiveEquilibriumProblem,
    owner_kind: Literal["homogeneous", "reactive_bubble"],
    owner: HomogeneousReactiveObservationDescriptor | ReactiveBubbleObservationDescriptor,
    certificate_scope: str,
) -> GeneralReactiveEquilibriumDescriptor:
    outputs = tuple(problem.outputs)
    units_match = all(
        epcsaft.unit_registry.Quantity(1.0, output.unit).to(cast(Any, row).unit).magnitude
        == 1.0
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
        "reaction_phase_ids": problem.reaction_phase_ids,
        "coupled_continuation_identity": problem.continuation_identity,
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
        "owner_kind": owner_kind,
        "owner_descriptor_fingerprint": owner.fingerprint,
        "certificate_scope": certificate_scope,
    }
    return GeneralReactiveEquilibriumDescriptor(
        "equilibrium.general-reactive",
        1,
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
        _fingerprint("epcsaft.grepe-descriptor.v1", payload),
        owner_kind,
        owner,
    )


def _homogeneous_branch_signature(owner: object) -> tuple[object, ...] | None:
    state = cast(Any, owner).state
    branch = cast(Any, owner).physical_branch
    chemical = cast(Any, owner).equilibrium
    if state is None or branch is None or chemical is None:
        return None
    return (
        branch.branch_family_identity,
        branch.continuation_identity,
        state.phase_identity,
        state.phase_role,
        state.chart_topology,
        state.phase_role_receipt.phase_selection_claim,
        state.chemical_search_certificate.selection_label,
    )


def _bubble_branch_signature(owner: object) -> tuple[object, ...] | None:
    state = cast(Any, owner).liquid_state
    branch = cast(Any, owner).physical_branch
    if state is None or branch is None:
        return None
    return (
        branch.branch_family_identity,
        branch.continuation_identity,
        state.phase_identity,
        state.phase_role,
        state.chart_topology,
        state.phase_role_receipt.phase_selection_claim,
        state.chemical_search_certificate.selection_label,
    )


def _unsupported(
    problem: GeneralReactiveEquilibriumProblem,
    diagnostic: str,
) -> GeneralReactiveEquilibriumResult:
    failure = NonEvaluableTrial("unsupported_topology", diagnostic)
    return GeneralReactiveEquilibriumResult(
        problem,
        None,
        "non_evaluable",
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
        (),
        failure,
        "unsupported",
        "not_established",
        None,
        None,
        None,
    )


@dataclass(frozen=True, slots=True)
class _HomogeneousCompilation:
    temperature: Quantity[Any]
    pressure: Quantity[Any]
    owner_problem: HomogeneousReactiveObservationProblem
    owner_descriptor: HomogeneousReactiveObservationDescriptor
    descriptor: GeneralReactiveEquilibriumDescriptor


@dataclass(frozen=True, slots=True)
class _BubbleCompilation:
    temperature: Quantity[Any]
    owner_problem: ReactiveBubbleVLEProblem
    owner_descriptor: ReactiveBubbleObservationDescriptor
    descriptor: GeneralReactiveEquilibriumDescriptor


def _compile_general_reactive_equilibrium(
    model: epcsaft.Mixture,
    problem: GeneralReactiveEquilibriumProblem,
    *,
    active_parameters: epcsaft.ActiveParameterSet | None = None,
) -> _HomogeneousCompilation | _BubbleCompilation | str:
    if not isinstance(model, epcsaft.Mixture):
        raise TypeError("GREPE currently requires an epcsaft.Mixture")
    if not isinstance(problem, GeneralReactiveEquilibriumProblem):
        raise TypeError("problem must be a GeneralReactiveEquilibriumProblem")
    _preflight(model, problem)
    if not isinstance(problem.temperature, Fixed):
        return "solved temperature is not admitted"
    temperature_k = _quantity(problem.temperature.value, "kelvin", "temperature")
    temperature = temperature_k * epcsaft.unit_registry.kelvin
    finite = tuple(phase for phase in problem.phases if phase.amount_role == "finite")
    incipient = tuple(phase for phase in problem.phases if phase.amount_role == "incipient")

    if len(finite) == 1 and not incipient and isinstance(problem.pressure, Fixed):
        phase = finite[0]
        if problem.reaction_phase_ids != (phase.identity,) or not isinstance(
            phase.support, AllComponents
        ):
            return "single-phase reaction/support topology is unsupported"
        pressure_pa = _quantity(problem.pressure.value, "pascal", "pressure")
        pressure = pressure_pa * epcsaft.unit_registry.pascal
        homogeneous_owner_problem = HomogeneousReactiveObservationProblem(
            problem.identity,
            phase.identity,
            phase.fluid_role,
            _phase_model(model, phase),
            problem.reaction_system,
            cast(str, phase.continuation_identity),
            cast(HomogeneousBranchPolicy, phase.branch_policy),
            phase.continuation_reference,
        )
        owner_rows = tuple(
            HomogeneousReactiveObservationRow(
                output.identity,
                cast(Any, _homogeneous_quantity(output)),
                output.coefficients,
                output.solvent_mass_coefficients_kg_per_mol,
                output.support,
                output.censor_limit,
                output.aggregate_identity,
                output.covariance_identity,
            )
            for output in problem.outputs
        )
        homogeneous_descriptor = homogeneous_reactive_observation_descriptor(
            homogeneous_owner_problem, owner_rows, active_parameters=active_parameters
        )
        descriptor = _descriptor(
            model,
            problem,
            "homogeneous",
            homogeneous_descriptor,
            "fixed_TP_single_finite_reacting_phase",
        )
        return _HomogeneousCompilation(
            temperature,
            pressure,
            homogeneous_owner_problem,
            homogeneous_descriptor,
            descriptor,
        )

    if (
        len(finite) == 1
        and len(incipient) == 1
        and isinstance(problem.pressure, Solved)
        and finite[0].fluid_role == "liquid"
        and incipient[0].fluid_role == "vapor"
        and problem.reaction_phase_ids == (finite[0].identity,)
        and isinstance(finite[0].support, AllComponents)
        and isinstance(incipient[0].support, DeclaredNeutralComponents)
    ):
        liquid, vapor = finite[0], incipient[0]
        liquid_model = _phase_model(model, liquid)
        if not isinstance(liquid_model, ProviderPhase):
            return "the admitted finite reacting liquid uses Provider EOS"
        liquid_problem = HomogeneousReactiveObservationProblem(
            f"{problem.identity}/finite/{liquid.identity}",
            liquid.identity,
            "liquid",
            liquid_model,
            problem.reaction_system,
            cast(str, liquid.continuation_identity),
            cast(HomogeneousBranchPolicy, liquid.branch_policy),
            liquid.continuation_reference,
        )
        if isinstance(vapor.model, ProviderModel):
            bubble_vapor_model: ProviderNonidealVapor | IdealGasVaporComparison = (
                ProviderNonidealVapor(vapor.model.reference_id)
            )
        else:
            bubble_vapor_model = IdealGasVaporComparison(vapor.model.reference_id)
        pressure_bounds = tuple(
            _quantity(value, "pascal", "pressure bound") for value in problem.pressure.bounds
        )
        pressure_starts = tuple(
            _quantity(value, "pascal", "pressure start") for value in problem.pressure.starts
        )
        initial_pressure = _quantity(problem.pressure.initial, "pascal", "pressure initial")
        if not pressure_bounds[0] <= initial_pressure <= pressure_bounds[1]:
            raise ValueError("GREPE pressure initial lies outside its bounds")
        bubble_owner_problem = ReactiveBubbleVLEProblem(
            problem.identity,
            liquid_problem,
            vapor.identity,
            cast(DeclaredNeutralComponents, vapor.support).component_ids,
            bubble_vapor_model,
            cast(tuple[float, float], pressure_bounds),
            pressure_starts,
            cast(str, problem.continuation_identity),
            problem.continuation_reference,
        )
        bubble_owner_rows = tuple(
            ReactiveBubbleObservationRow(
                output.identity,
                cast(Any, _bubble_quantity(output, liquid, vapor)),
                output.coefficients,
                output.solvent_mass_coefficients_kg_per_mol,
                cast(Any, output.support),
                output.censor_limit,
                output.aggregate_identity,
                output.covariance_identity,
            )
            for output in problem.outputs
        )
        bubble_descriptor = reactive_bubble_observation_descriptor(
            bubble_owner_problem,
            bubble_owner_rows,
            active_parameters=active_parameters,
        )
        descriptor = _descriptor(
            model,
            problem,
            "reactive_bubble",
            bubble_descriptor,
            "fixed_T_one_finite_liquid_one_incipient_vapor",
        )
        return _BubbleCompilation(
            temperature,
            bubble_owner_problem,
            bubble_descriptor,
            descriptor,
        )

    return "declared reactive phase/state topology is not admitted"


def general_reactive_equilibrium_descriptor(
    model: epcsaft.Mixture,
    problem: GeneralReactiveEquilibriumProblem,
    *,
    active_parameters: epcsaft.ActiveParameterSet | None = None,
) -> GeneralReactiveEquilibriumDescriptor:
    """Describe one admitted GREPE request without evaluating thermodynamics."""

    compiled = _compile_general_reactive_equilibrium(
        model, problem, active_parameters=active_parameters
    )
    if isinstance(compiled, str):
        raise ValueError(compiled)
    return compiled.descriptor


def certify_general_reactive_equilibrium_continuation_reference(
    model: epcsaft.Mixture,
    problem: GeneralReactiveEquilibriumProblem,
    *,
    active_parameters: epcsaft.ActiveParameterSet | None = None,
    maximum_log_pressure_distance: float,
    maximum_log_composition_distance: float,
    maximum_log_volume_distance: float,
) -> ReactiveBubbleContinuationReference:
    """Certify the local coupled branch of an admitted incipient GREPE request."""

    if problem.continuation_reference is not None:
        raise ValueError("GREPE continuation certification requires an unanchored request")
    compiled = _compile_general_reactive_equilibrium(
        model, problem, active_parameters=active_parameters
    )
    if isinstance(compiled, str):
        raise ValueError(compiled)
    if not isinstance(compiled, _BubbleCompilation):
        raise ValueError("this admitted GREPE formulation has no coupled continuation")
    return certify_reactive_bubble_continuation_reference(
        compiled.owner_problem,
        compiled.temperature,
        active_parameters=active_parameters,
        maximum_log_pressure_distance=maximum_log_pressure_distance,
        maximum_log_composition_distance=maximum_log_composition_distance,
        maximum_log_volume_distance=maximum_log_volume_distance,
    )


def solve_general_reactive_equilibrium(
    model: epcsaft.Mixture,
    problem: GeneralReactiveEquilibriumProblem,
    *,
    active_parameters: epcsaft.ActiveParameterSet | None = None,
) -> GeneralReactiveEquilibriumResult:
    """Compile and solve one currently admitted role-declared GREPE condition."""

    compiled = _compile_general_reactive_equilibrium(
        model, problem, active_parameters=active_parameters
    )
    if isinstance(compiled, str):
        return _unsupported(problem, compiled)
    if isinstance(compiled, _HomogeneousCompilation):
        homogeneous_owner = evaluate_homogeneous_reactive_observations(
            compiled.owner_problem,
            compiled.temperature,
            compiled.pressure,
            compiled.owner_descriptor,
            active_parameters=active_parameters,
        )
        return GeneralReactiveEquilibriumResult(
            problem,
            compiled.descriptor,
            homogeneous_owner.status,
            _generic_rows(problem.outputs, cast(tuple[object, ...], homogeneous_owner.rows)),
            (),
            homogeneous_owner.failure,
            compiled.descriptor.certificate_scope,
            "not_established",
            homogeneous_owner.trial,
            _homogeneous_branch_signature(homogeneous_owner),
            homogeneous_owner,
        )
    bubble_owner = evaluate_reactive_bubble_observations(
        compiled.owner_problem,
        compiled.temperature,
        compiled.owner_descriptor,
        active_parameters=active_parameters,
    )
    return GeneralReactiveEquilibriumResult(
        problem,
        compiled.descriptor,
        bubble_owner.status,
        _generic_rows(problem.outputs, cast(tuple[object, ...], bubble_owner.rows)),
        bubble_owner.phases,
        bubble_owner.failure,
        compiled.descriptor.certificate_scope,
        "not_established",
        bubble_owner.trial,
        _bubble_branch_signature(bubble_owner),
        bubble_owner,
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


def serialize_general_reactive_equilibrium_receipt(
    result: GeneralReactiveEquilibriumResult,
    *,
    batch_identity: str,
    observations: tuple[GeneralEquilibriumOwnerObservation, ...],
) -> str:
    """Serialize one canonical opaque GREPE owner receipt."""

    descriptor = result.descriptor
    if result.status != "evaluated" or descriptor is None or result.owner_result is None:
        raise ValueError("only evaluated GREPE results have owner receipts")
    selected = _validated_observation_selection(descriptor, observations)
    if descriptor._owner_kind == "homogeneous":
        if not isinstance(result.problem.temperature, Fixed) or not isinstance(
            result.problem.pressure, Fixed
        ):
            raise ValueError("homogeneous GREPE receipt requires fixed T and P")
        owner_json = serialize_homogeneous_reactive_owner_receipt(
            cast(Any, result.owner_result),
            batch_identity=batch_identity,
            requested_temperature=result.problem.temperature.value,
            requested_pressure=result.problem.pressure.value,
            observations=tuple(
                HomogeneousReactiveOwnerObservation(
                    item.target_identity,
                    item.prediction_identity,
                    item.aggregate_identity,
                    item.covariance_identity,
                )
                for item in selected
            ),
        )
    else:
        owner_json = serialize_reactive_bubble_owner_receipt(
            cast(Any, result.owner_result),
            batch_identity=batch_identity,
        )
    return _canonical_json(
        {
            "schema": "epcsaft.grepe-owner-receipt.v1",
            "batch_identity": batch_identity,
            "descriptor_fingerprint": descriptor.fingerprint,
            "owner_kind": descriptor._owner_kind,
            "owner_receipt": json.loads(owner_json),
            "observations": tuple(_owner_observation_payload(item) for item in selected),
        }
    )


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
    expected_observations = tuple(_owner_observation_payload(item) for item in selected)
    if (
        payload.get("schema") != "epcsaft.grepe-owner-receipt.v1"
        or payload.get("batch_identity") != batch_identity
        or payload.get("descriptor_fingerprint") != descriptor.fingerprint
        or payload.get("owner_kind") != descriptor._owner_kind
        or payload.get("observations") != list(expected_observations)
    ):
        raise ValueError("GREPE owner receipt identity or observations do not match")
    owner_payload = payload.get("owner_receipt")
    if not isinstance(owner_payload, dict):
        raise TypeError("GREPE owner receipt payload is missing")
    owner_json = _canonical_json(owner_payload)
    if descriptor._owner_kind == "homogeneous":
        if not isinstance(problem.temperature, Fixed) or not isinstance(problem.pressure, Fixed):
            raise ValueError("homogeneous GREPE receipt requires fixed T and P")
        validate_homogeneous_reactive_owner_receipt(
            owner_json,
            descriptor=cast(HomogeneousReactiveObservationDescriptor, descriptor._owner_descriptor),
            batch_identity=batch_identity,
            requested_temperature=problem.temperature.value,
            requested_pressure=problem.pressure.value,
            active_parameters=active_parameters,
            observations=tuple(
                HomogeneousReactiveOwnerObservation(
                    item.target_identity,
                    item.prediction_identity,
                    item.aggregate_identity,
                    item.covariance_identity,
                )
                for item in selected
            ),
        )
    else:
        if not isinstance(problem.temperature, Fixed):
            raise ValueError("reactive bubble GREPE receipt requires fixed temperature")
        validate_reactive_bubble_owner_receipt(
            owner_json,
            descriptor=cast(ReactiveBubbleObservationDescriptor, descriptor._owner_descriptor),
            batch_identity=batch_identity,
            requested_temperature=problem.temperature.value,
            active_parameters=active_parameters,
        )
    return f"sha256:{hashlib.sha256(receipt_json.encode('utf-8')).hexdigest()}"


__all__ = (
    "EquilibriumOutput",
    "Fixed",
    "GeneralEquilibriumRowResult",
    "GeneralEquilibriumOwnerObservation",
    "GeneralReactiveEquilibriumDescriptor",
    "GeneralReactiveEquilibriumProblem",
    "GeneralReactiveEquilibriumResult",
    "IdealGasModel",
    "ProviderModel",
    "ReactivePhase",
    "Solved",
    "certify_general_reactive_equilibrium_continuation_reference",
    "general_reactive_equilibrium_descriptor",
    "serialize_general_reactive_equilibrium_receipt",
    "solve_general_reactive_equilibrium",
    "validate_general_reactive_equilibrium_receipt",
)
