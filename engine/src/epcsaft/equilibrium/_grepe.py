from __future__ import annotations

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
    HomogeneousBranchPolicy,
    HomogeneousContinuationReference,
    HomogeneousReactiveObservationProblem,
    HomogeneousReactiveObservationRow,
    evaluate_homogeneous_reactive_observations,
    homogeneous_reactive_observation_descriptor,
)
from ._reactive_vle import (
    IdealGasVaporComparison,
    ProviderNonidealVapor,
    ReactiveBubbleContinuationReference,
    ReactiveBubbleObservationRow,
    ReactiveBubbleVLEProblem,
    evaluate_reactive_bubble_observations,
    reactive_bubble_observation_descriptor,
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
        if not self.identity or self.selector not in _EQUILIBRIUM_SELECTORS:
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
    status: Literal["evaluated", "non_evaluable"]
    value: float | None
    jacobian: tuple[float, ...] | None
    reason: str


@dataclass(frozen=True, slots=True)
class GeneralReactiveEquilibriumResult:
    """One generic result over a locally certified admitted GREPE formulation."""

    problem: GeneralReactiveEquilibriumProblem
    status: Literal["evaluated", "non_evaluable"]
    rows: tuple[GeneralEquilibriumRowResult, ...]
    phases: tuple[EquilibriumPhaseState, ...]
    failure: NonEvaluableTrial | None
    certificate_scope: str
    globality_claim: Literal["not_established"]
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


def _generic_rows(rows: tuple[object, ...]) -> tuple[GeneralEquilibriumRowResult, ...]:
    return tuple(
        GeneralEquilibriumRowResult(
            identity=cast(Any, row).identity,
            status=cast(Any, row).status,
            value=cast(Any, row).value,
            jacobian=cast(Any, row).jacobian,
            reason=(
                cast(Any, row).reason
                if hasattr(row, "reason")
                else "" if cast(Any, row).failure is None else cast(Any, row).failure.diagnostic
            ),
        )
        for row in rows
    )


def _unsupported(
    problem: GeneralReactiveEquilibriumProblem,
    diagnostic: str,
) -> GeneralReactiveEquilibriumResult:
    failure = NonEvaluableTrial("unsupported_topology", diagnostic)
    return GeneralReactiveEquilibriumResult(
        problem,
        "non_evaluable",
        tuple(
            GeneralEquilibriumRowResult(output.identity, "non_evaluable", None, None, diagnostic)
            for output in problem.outputs
        ),
        (),
        failure,
        "unsupported",
        "not_established",
        None,
    )


def solve_general_reactive_equilibrium(
    model: epcsaft.Mixture,
    problem: GeneralReactiveEquilibriumProblem,
    *,
    active_parameters: epcsaft.ActiveParameterSet | None = None,
) -> GeneralReactiveEquilibriumResult:
    """Compile and solve one currently admitted role-declared GREPE condition."""

    if not isinstance(model, epcsaft.Mixture):
        raise TypeError("GREPE currently requires an epcsaft.Mixture")
    if not isinstance(problem, GeneralReactiveEquilibriumProblem):
        raise TypeError("problem must be a GeneralReactiveEquilibriumProblem")
    _preflight(model, problem)
    if not isinstance(problem.temperature, Fixed):
        return _unsupported(problem, "solved temperature is not admitted")
    temperature_k = _quantity(problem.temperature.value, "kelvin", "temperature")
    temperature = temperature_k * epcsaft.unit_registry.kelvin
    finite = tuple(phase for phase in problem.phases if phase.amount_role == "finite")
    incipient = tuple(phase for phase in problem.phases if phase.amount_role == "incipient")

    if len(finite) == 1 and not incipient and isinstance(problem.pressure, Fixed):
        phase = finite[0]
        if problem.reaction_phase_ids != (phase.identity,) or not isinstance(
            phase.support, AllComponents
        ):
            return _unsupported(problem, "single-phase reaction/support topology is unsupported")
        pressure_pa = _quantity(problem.pressure.value, "pascal", "pressure")
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
        homogeneous_owner = evaluate_homogeneous_reactive_observations(
            homogeneous_owner_problem,
            temperature,
            pressure_pa * epcsaft.unit_registry.pascal,
            homogeneous_descriptor,
            active_parameters=active_parameters,
        )
        return GeneralReactiveEquilibriumResult(
            problem,
            homogeneous_owner.status,
            _generic_rows(cast(tuple[object, ...], homogeneous_owner.rows)),
            (),
            homogeneous_owner.failure,
            "fixed_TP_single_finite_reacting_phase",
            "not_established",
            homogeneous_owner,
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
            return _unsupported(problem, "the admitted finite reacting liquid uses Provider EOS")
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
        bubble_owner = evaluate_reactive_bubble_observations(
            bubble_owner_problem,
            temperature,
            bubble_descriptor,
            active_parameters=active_parameters,
        )
        return GeneralReactiveEquilibriumResult(
            problem,
            bubble_owner.status,
            _generic_rows(cast(tuple[object, ...], bubble_owner.rows)),
            bubble_owner.phases,
            bubble_owner.failure,
            "fixed_T_one_finite_liquid_one_incipient_vapor",
            "not_established",
            bubble_owner,
        )

    return _unsupported(problem, "declared reactive phase/state topology is not admitted")


__all__ = (
    "EquilibriumOutput",
    "Fixed",
    "GeneralEquilibriumRowResult",
    "GeneralReactiveEquilibriumProblem",
    "GeneralReactiveEquilibriumResult",
    "IdealGasModel",
    "ProviderModel",
    "ReactivePhase",
    "Solved",
    "solve_general_reactive_equilibrium",
)
