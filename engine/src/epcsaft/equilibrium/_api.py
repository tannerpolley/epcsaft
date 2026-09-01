from __future__ import annotations

import hashlib
import json
import math
from collections.abc import Callable, Mapping, Sequence
from dataclasses import asdict, dataclass, replace
from functools import lru_cache
from itertools import pairwise
from typing import TYPE_CHECKING, Any, Literal, cast, overload

from pint import Quantity
from pint.errors import DimensionalityError

import epcsaft

from .. import _core
from ..eos import _MOLAR_GAS_CONSTANT


class _ScalarBoundarySolveError(ValueError):
    """A bounded scalar boundary could not produce an accepted root."""

    def __init__(self, reason: str, value: float, state: object) -> None:
        self.reason = reason
        self.state = state
        super().__init__(f"scalar boundary solve {reason} at {value!r}")


def _solve_bracketed_scalar[ScalarState](
    evaluate: Callable[[float], tuple[ScalarState, float, float]],
    lower: float,
    upper: float,
    *,
    initial: float | None,
    bracket_points: Sequence[float] = (),
    tolerance: float,
    maximum_iterations: int,
) -> tuple[float, ScalarState, float]:
    """Solve one bounded scalar residual with safeguarded exact-derivative Newton."""

    if not lower < upper:
        raise ValueError("scalar boundary bounds must be ordered")
    if not math.isfinite(tolerance) or tolerance <= 0.0:
        raise ValueError("scalar boundary tolerance must be positive and finite")
    if maximum_iterations <= 0:
        raise ValueError("scalar boundary iteration limit must be positive")
    if initial is not None and not lower <= initial <= upper:
        raise ValueError("scalar boundary initial must lie inside its bounds")

    def checked(value: float) -> tuple[ScalarState, float, float]:
        state, residual, derivative = evaluate(value)
        if not math.isfinite(residual):
            raise _ScalarBoundarySolveError("invalid_residual", value, state)
        return state, residual, derivative

    anchor = 0.5 * (lower + upper) if initial is None else initial
    evaluations: dict[float, tuple[ScalarState, float, float]] = {}
    if initial is not None:
        evaluations[initial] = checked(initial)
        if abs(evaluations[initial][1]) <= tolerance:
            return initial, evaluations[initial][0], evaluations[initial][1]
    points = tuple(sorted({lower, *bracket_points, upper}))
    if any(not math.isfinite(point) or not lower <= point <= upper for point in points):
        raise ValueError("scalar boundary bracket points must lie inside its bounds")
    for point in points:
        if point not in evaluations:
            evaluations[point] = checked(point)
    roots = tuple(point for point in points if abs(evaluations[point][1]) <= tolerance)
    if roots:
        point = min(roots, key=lambda value: abs(value - anchor))
        state, residual, _ = evaluations[point]
        return point, state, residual
    brackets = tuple(
        (left, right)
        for left, right in pairwise(points)
        if math.copysign(1.0, evaluations[left][1]) != math.copysign(1.0, evaluations[right][1])
    )
    if not brackets:
        raise _ScalarBoundarySolveError("not_bracketed", lower, evaluations[lower][0])
    left, right = min(
        brackets,
        key=lambda pair: abs(anchor - 0.5 * (pair[0] + pair[1])),
    )
    left_residual = evaluations[left][1]

    if initial is None:
        current = 0.5 * (left + right)
        current_state, current_residual, current_derivative = checked(current)
    elif left <= initial <= right:
        current = initial
        current_state, current_residual, current_derivative = evaluations[initial]
    else:
        current = 0.5 * (left + right)
        if current in evaluations:
            current_state, current_residual, current_derivative = evaluations[current]
        else:
            current_state, current_residual, current_derivative = checked(current)

    for _ in range(maximum_iterations):
        if abs(current_residual) <= tolerance:
            return current, current_state, current_residual
        if math.copysign(1.0, current_residual) == math.copysign(1.0, left_residual):
            left, left_residual = current, current_residual
        else:
            right = current
        candidate = 0.5 * (left + right)
        if not math.isfinite(current_derivative):
            raise _ScalarBoundarySolveError("invalid_derivative", current, current_state)
        if current_derivative != 0.0:
            newton = current - current_residual / current_derivative
            if left < newton < right:
                candidate = newton
        if candidate == current:
            raise _ScalarBoundarySolveError("iteration_limit", current, current_state)
        current = candidate
        current_state, current_residual, current_derivative = checked(current)

    if abs(current_residual) <= tolerance:
        return current, current_state, current_residual
    raise _ScalarBoundarySolveError("iteration_limit", current, current_state)


if TYPE_CHECKING:
    from ._grepe import GeneralReactiveEquilibriumProblem, GeneralReactiveEquilibriumResult


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
    eos_evaluation_count: int
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
    eos_evaluation_count: int
    failure_reason: str


PhaseRole = Literal["liquid", "vapor"]


@dataclass(frozen=True, slots=True)
class AllComponents:
    """Every model component is admitted as a positive phase coordinate."""


@dataclass(frozen=True, slots=True)
class NeutralComponents:
    """Only the model's neutral components are present in this phase."""


@dataclass(frozen=True, slots=True)
class DeclaredNeutralComponents:
    """An ordered caller-declared nonempty subset of neutral model components."""

    component_ids: tuple[str, ...]

    def __post_init__(self) -> None:
        component_ids = tuple(self.component_ids)
        if not component_ids or any(not identity for identity in component_ids):
            raise ValueError("declared neutral component support must be nonempty")
        if len(set(component_ids)) != len(component_ids):
            raise ValueError("declared neutral component support contains duplicates")
        object.__setattr__(self, "component_ids", component_ids)


@dataclass(frozen=True, slots=True)
class DeclaredComponents:
    """An ordered caller-declared nonempty finite-phase species support."""

    component_ids: tuple[str, ...]

    def __post_init__(self) -> None:
        component_ids = tuple(self.component_ids)
        if not component_ids or any(not identity for identity in component_ids):
            raise ValueError("declared component support must be nonempty")
        if len(set(component_ids)) != len(component_ids):
            raise ValueError("declared component support contains duplicates")
        object.__setattr__(self, "component_ids", component_ids)


PhaseSupport = AllComponents | NeutralComponents | DeclaredNeutralComponents


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


@dataclass(frozen=True, slots=True)
class FixedPCoexistence:
    """Fixed-pressure two-phase boundary on one declared temperature bracket."""

    pressure: Quantity[Any]
    temperature_bounds: tuple[Quantity[Any], Quantity[Any]]
    topology: DeclaredTopology
    reference_phase_index: int
    reference_mole_fractions: tuple[float, ...]
    temperature_initial: Quantity[Any] | None = None
    packing_fraction_intervals: tuple[tuple[float, float], ...] | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "temperature_bounds", tuple(self.temperature_bounds))
        object.__setattr__(self, "reference_mole_fractions", tuple(self.reference_mole_fractions))
        if self.packing_fraction_intervals is not None:
            object.__setattr__(
                self,
                "packing_fraction_intervals",
                tuple(tuple(interval) for interval in self.packing_fraction_intervals),
            )


PhaseEquilibriumSpecification = FixedTPFeed | FixedTCoexistence | FixedPCoexistence


@dataclass(frozen=True, slots=True)
class PhaseHint:
    """Non-authoritative start-order hint bound to one phase index."""

    phase_index: int
    mole_fractions: tuple[float, ...]
    molar_density_mol_m3: float

    def __post_init__(self) -> None:
        object.__setattr__(self, "mole_fractions", tuple(self.mole_fractions))


@dataclass(frozen=True, slots=True)
class DeclaredPhaseInitialGuess:
    """Authoritative initial variables for one declared finite phase."""

    phase_index: int
    mole_fractions: tuple[float, ...]
    molar_density_mol_m3: float
    phase_fraction: float

    def __post_init__(self) -> None:
        object.__setattr__(self, "mole_fractions", tuple(self.mole_fractions))


@dataclass(frozen=True, slots=True)
class PhaseEquilibriumProblem:
    """One typed HELD-lineage feed or coexistence specification."""

    specification: PhaseEquilibriumSpecification
    hints: tuple[PhaseHint, ...] = ()
    initial_guesses: tuple[DeclaredPhaseInitialGuess, ...] = ()

    def __post_init__(self) -> None:
        object.__setattr__(self, "hints", tuple(self.hints))
        object.__setattr__(self, "initial_guesses", tuple(self.initial_guesses))


@dataclass(frozen=True, slots=True)
class EquilibriumPhaseState:
    """One intensive phase state on its declared component support."""

    role: PhaseRole
    support: Literal[
        "all_components",
        "neutral_components",
        "declared_neutral_components",
        "declared_components",
    ]
    supported_component_ids: tuple[str, ...]
    mole_fractions: tuple[float, ...]
    molar_volume_m3_per_mol: float
    molar_density_mol_m3: float
    packing_fraction: float
    pressure_pa: float
    chemical_potential_component_ids: tuple[str, ...]
    chemical_potential_over_rt: tuple[float, ...]
    mechanical_class: str


@dataclass(frozen=True, slots=True)
class EquilibriumContinuationState:
    """Converged generic phase state reusable by a neighboring solve."""

    phase_identities: tuple[str, ...]
    phases: tuple[EquilibriumPhaseState, ...]
    pressure_pa: float
    problem_fingerprint: str
    reaction_fingerprint: str
    topology_fingerprint: str
    parameter_fingerprint: str

    def __post_init__(self) -> None:
        object.__setattr__(self, "phase_identities", tuple(self.phase_identities))
        object.__setattr__(self, "phases", tuple(self.phases))
        if (
            not self.phases
            or len(self.phase_identities) != len(self.phases)
            or len(set(self.phase_identities)) != len(self.phase_identities)
            or any(not identity for identity in self.phase_identities)
            or not math.isfinite(self.pressure_pa)
            or self.pressure_pa <= 0.0
            or any(
                not value.startswith("sha256:")
                for value in (
                    self.problem_fingerprint,
                    self.reaction_fingerprint,
                    self.topology_fingerprint,
                    self.parameter_fingerprint,
                )
            )
        ):
            raise ValueError("equilibrium continuation state is incomplete")

    def phase(self, identity: str) -> EquilibriumPhaseState:
        try:
            return self.phases[self.phase_identities.index(identity)]
        except ValueError as error:
            raise ValueError(f"continuation state has no phase {identity!r}") from error

    def to_mapping(self) -> dict[str, object]:
        return {
            "phase_identities": self.phase_identities,
            "phases": tuple(asdict(phase) for phase in self.phases),
            "pressure_pa": self.pressure_pa,
            "problem_fingerprint": self.problem_fingerprint,
            "reaction_fingerprint": self.reaction_fingerprint,
            "topology_fingerprint": self.topology_fingerprint,
            "parameter_fingerprint": self.parameter_fingerprint,
        }

    @property
    def fingerprint(self) -> str:
        payload = json.dumps(
            self.to_mapping(),
            allow_nan=False,
            separators=(",", ":"),
            sort_keys=True,
        ).encode()
        return "sha256:" + hashlib.sha256(payload).hexdigest()


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
    eos_evaluation_count: int
    hint_count: int
    hint_reordered_start_count: int
    initial_guess_count: int
    accepted_start_count: int
    distinct_root_count: int
    final_correction_attempted: bool
    final_correction_accepted: bool
    final_correction_iteration_count: int
    pre_correction_residual_inf: float
    post_correction_residual_inf: float
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
class _ReferenceDirectionSensitivity:
    """Exact coexistence response to admissible reference-composition directions."""

    status: Literal["available", "unavailable"]
    failure_reason: str | None
    identities: tuple[str, ...]
    residual_identities: tuple[str, ...]
    state_identities: tuple[str, ...]
    output_identities: tuple[str, ...]
    residual_state_jacobian: tuple[tuple[float, ...], ...]
    residual_direction_jacobian: tuple[tuple[float, ...], ...]
    state_direction_jacobian: tuple[tuple[float, ...], ...]
    output_direction_jacobian: tuple[tuple[float, ...], ...]
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


@dataclass(frozen=True, slots=True)
class IncipientTemperatureBoundary:
    """Two-phase boundary at fixed pressure with solved temperature."""

    reference_phase_index: int
    pressure_pa: float
    temperature_k: float


EquilibriumExtent = FeedSplit | IncipientBoundary | IncipientTemperatureBoundary


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
    reaction_id: str | None = None
    correlation_kind: (
        Literal["ln-k-a-plus-b-over-t", "negative-log10-temperature-polynomial"] | None
    ) = None
    coefficient_identities: tuple[str, ...] = ()
    coefficient_values: tuple[float, ...] = ()

    def __post_init__(self) -> None:
        object.__setattr__(self, "coefficient_identities", tuple(self.coefficient_identities))
        object.__setattr__(self, "coefficient_values", tuple(self.coefficient_values))
        declarations = (self.reaction_id, self.correlation_kind)
        if any(value is None for value in declarations):
            if not all(value is None for value in declarations) or self.coefficient_identities:
                raise ValueError("reaction correlation declaration is partial")
            return
        expected = 2 if self.correlation_kind == "ln-k-a-plus-b-over-t" else 3
        if (
            self.correlation_kind
            not in {
                "ln-k-a-plus-b-over-t",
                "negative-log10-temperature-polynomial",
            }
            or len(self.coefficient_identities) != expected
            or len(self.coefficient_values) != expected
        ):
            raise ValueError("reaction correlation declaration is incompatible")


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
    source_reference_convention: Literal[
        "molality-infinite-dilution", "mole-fraction-mixed-reference"
    ]
    source_reference_activity_convention_id: str
    source_reference_standard_molality_mol_per_kg: float
    source_reference_pure_components: tuple[epcsaft.SourceReferencePureComponent, ...] = ()

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
        object.__setattr__(
            self,
            "source_reference_pure_components",
            tuple(self.source_reference_pure_components),
        )


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


def _continuation_start(
    reaction: ChemicalEquilibriumProblem,
    reference: EquilibriumPhaseState | None,
) -> tuple[tuple[float, ...], float] | None:
    if reference is None:
        return None
    projected = tuple(
        math.fsum(
            coefficient * fraction
            for coefficient, fraction in zip(row, reference.mole_fractions, strict=True)
        )
        for row in reaction.balance_matrix
    )
    denominator = math.fsum(value * value for value in projected)
    if not math.isfinite(denominator) or denominator <= 0.0:
        raise ValueError("continuation composition does not determine the conserved scale")
    total = (
        math.fsum(
            value * target
            for value, target in zip(projected, reaction.conserved_totals, strict=True)
        )
        / denominator
    )
    return (
        tuple(total * fraction for fraction in reference.mole_fractions),
        total * reference.molar_volume_m3_per_mol,
    )


def _contract_phase_output(
    selector: str,
    coefficients: tuple[float, ...],
    solvent_coefficients: tuple[float, ...],
    amounts: tuple[float, ...],
    volume: float,
    amount_derivatives: tuple[tuple[float, ...], ...],
    volume_derivatives: tuple[float, ...],
) -> tuple[float, tuple[float, ...]]:
    total = math.fsum(amounts)
    total_derivatives = tuple(math.fsum(row) for row in amount_derivatives)
    if selector == "phase.molar_density":
        return total / volume, tuple(
            dtotal / volume - total * dvolume / volume**2
            for dtotal, dvolume in zip(total_derivatives, volume_derivatives, strict=True)
        )
    numerator = math.fsum(
        coefficient * amount for coefficient, amount in zip(coefficients, amounts, strict=True)
    )
    numerator_derivatives = tuple(
        math.fsum(
            coefficient * derivative
            for coefficient, derivative in zip(coefficients, row, strict=True)
        )
        for row in amount_derivatives
    )
    if selector == "phase.amount":
        return numerator, numerator_derivatives
    if selector == "phase.mole_fraction":
        return numerator / total, tuple(
            (derivative * total - numerator * dtotal) / total**2
            for derivative, dtotal in zip(numerator_derivatives, total_derivatives, strict=True)
        )
    if selector == "phase.molality":
        denominator = math.fsum(
            coefficient * amount
            for coefficient, amount in zip(solvent_coefficients, amounts, strict=True)
        )
        if denominator <= 0.0:
            raise ValueError("phase molality denominator is not positive")
        denominator_derivatives = tuple(
            math.fsum(
                coefficient * derivative
                for coefficient, derivative in zip(solvent_coefficients, row, strict=True)
            )
            for row in amount_derivatives
        )
        return numerator / denominator, tuple(
            (derivative * denominator - numerator * ddenominator) / denominator**2
            for derivative, ddenominator in zip(
                numerator_derivatives, denominator_derivatives, strict=True
            )
        )
    raise ValueError(f"unsupported generic phase output {selector!r}")


@dataclass(frozen=True, slots=True)
class ActiveParameterTrialReceipt:
    parameter_fingerprint: str
    identities: tuple[str, ...]
    values: tuple[float, ...]
    units: tuple[str, ...]
    active_parameter_fingerprint: str
    fingerprint: str


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
        raise ValueError(f"solve requires exactly {component_count} components")
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
        eos_evaluation_count=int(cast(int, payload["eos_evaluation_count"])),
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
        eos_evaluation_count=0,
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
        eos_evaluation_count=int(cast(int, payload["eos_evaluation_count"])),
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
        raise ValueError("native HELD2 result has the wrong EOS fingerprint")
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


def _reference_direction_sensitivity(
    value: object,
    identities: tuple[str, ...],
) -> _ReferenceDirectionSensitivity:
    if not isinstance(value, Mapping):
        raise TypeError("native reference-direction sensitivity is missing")
    status = str(value["status"])
    if status not in {"available", "unavailable"}:
        raise ValueError("native reference-direction sensitivity status is invalid")
    unavailable = status == "unavailable"
    failure = value["failure_reason"]
    if (unavailable and not isinstance(failure, str)) or (not unavailable and failure is not None):
        raise ValueError("native reference-direction sensitivity failure state is invalid")
    native_identities = tuple(
        str(item) for item in cast(Sequence[object], value["parameter_identities"])
    )
    if native_identities != identities:
        raise ValueError("native reference-direction identity or order changed")
    dimension = int(cast(int, value["dimension"]))
    direction_count = int(cast(int, value["parameter_count"]))
    output_count = int(cast(int, value["output_count"]))
    rank = int(cast(int, value["rank"]))
    if direction_count != len(identities) or min(dimension, output_count, rank) < 0:
        raise ValueError("native reference-direction sensitivity dimensions are invalid")
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
        raise ValueError("native reference-direction residual/state identities are incomplete")
    if len(output_identities) != output_count:
        raise ValueError("native reference-direction output identities are incomplete")
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
        raise ValueError("native reference-direction linear solve is uncertified")
    return _ReferenceDirectionSensitivity(
        status=cast(Literal["available", "unavailable"], status),
        failure_reason=None if failure is None else str(failure),
        identities=identities,
        residual_identities=residual_identities,
        state_identities=state_identities,
        output_identities=output_identities,
        residual_state_jacobian=_row_major_matrix(
            value["residual_state_jacobian"],
            dimension,
            dimension,
            "reference-direction residual-state Jacobian",
            unavailable=unavailable,
        ),
        residual_direction_jacobian=_row_major_matrix(
            value["residual_parameter_jacobian"],
            dimension,
            direction_count,
            "reference-direction residual-direction Jacobian",
            unavailable=unavailable,
        ),
        state_direction_jacobian=_row_major_matrix(
            value["state_parameter_jacobian"],
            dimension,
            direction_count,
            "reference-direction state-direction Jacobian",
            unavailable=unavailable,
        ),
        output_direction_jacobian=_row_major_matrix(
            value["output_parameter_jacobian"],
            output_count,
            direction_count,
            "reference-direction output-direction Jacobian",
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
        native_model = model._native_model(active_parameters)
    except (TypeError, ValueError) as error:
        diagnostics = _failed_held_diagnostics("invalid_input", "input_rejected", str(error))
        raise _HeldDiscoveryError(str(error), diagnostics) from error

    try:
        native = cast(
            Mapping[str, object],
            _core._solve_held_discovery(
                native_model,
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
            raise ValueError("native HELD result has the wrong eos fingerprint")
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
        eos_evaluation_count=int(cast(int, payload["eos_evaluation_count"])),
        hint_count=int(cast(int, payload["hint_count"])),
        hint_reordered_start_count=int(cast(int, payload["hint_reordered_start_count"])),
        initial_guess_count=int(cast(int, payload["initial_guess_count"])),
        accepted_start_count=int(cast(int, payload["accepted_start_count"])),
        distinct_root_count=int(cast(int, payload["distinct_root_count"])),
        final_correction_attempted=bool(payload["final_correction_attempted"]),
        final_correction_accepted=bool(payload["final_correction_accepted"]),
        final_correction_iteration_count=int(
            cast(int, payload["final_correction_iteration_count"])
        ),
        pre_correction_residual_inf=float(cast(float, payload["pre_correction_residual_inf"])),
        post_correction_residual_inf=float(cast(float, payload["post_correction_residual_inf"])),
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
        "declared_neutral_components",
        "declared_components",
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
        support=cast(
            Literal[
                "all_components",
                "neutral_components",
                "declared_neutral_components",
                "declared_components",
            ],
            support,
        ),
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
    if type(support) is DeclaredNeutralComponents:
        return "declared_neutral_components"
    raise TypeError(
        "phase support must be AllComponents, NeutralComponents, or DeclaredNeutralComponents"
    )


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


def _validated_phase_initial_guesses(
    guesses: tuple[DeclaredPhaseInitialGuess, ...],
    phase_count: int,
    component_counts: tuple[int, ...],
) -> tuple[
    tuple[int, ...],
    tuple[tuple[float, ...], ...],
    tuple[float, ...],
    tuple[float, ...],
]:
    if any(type(guess) is not DeclaredPhaseInitialGuess for guess in guesses):
        raise TypeError("initial_guesses must contain DeclaredPhaseInitialGuess values")
    if guesses and len(guesses) != phase_count:
        raise ValueError("initial_guesses must cover every declared phase")
    indices = tuple(guess.phase_index for guess in guesses)
    if len(indices) != len(set(indices)) or any(
        index < 0 or index >= phase_count for index in indices
    ):
        raise ValueError("declared phase initial-guess indices must be unique and in range")
    for guess in guesses:
        if len(guess.mole_fractions) != component_counts[guess.phase_index]:
            raise ValueError(
                "declared phase initial-guess composition has the wrong component count"
            )
        if not all(math.isfinite(value) and value > 0.0 for value in guess.mole_fractions):
            raise ValueError("declared phase initial-guess composition must be positive and finite")
        if not math.isclose(math.fsum(guess.mole_fractions), 1.0, rel_tol=0.0, abs_tol=1.0e-12):
            raise ValueError(
                "declared phase initial-guess composition must sum to one within 1e-12"
            )
        if not math.isfinite(guess.molar_density_mol_m3) or guess.molar_density_mol_m3 <= 0.0:
            raise ValueError("declared phase initial-guess density must be positive and finite")
        if not math.isfinite(guess.phase_fraction) or guess.phase_fraction <= 0.0:
            raise ValueError("declared phase initial-guess fraction must be positive and finite")
    if guesses and not math.isclose(
        math.fsum(guess.phase_fraction for guess in guesses),
        1.0,
        rel_tol=0.0,
        abs_tol=1.0e-12,
    ):
        raise ValueError("declared phase initial-guess fractions must sum to one within 1e-12")
    return (
        indices,
        tuple(guess.mole_fractions for guess in guesses),
        tuple(guess.molar_density_mol_m3 for guess in guesses),
        tuple(guess.phase_fraction for guess in guesses),
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
        eos_evaluation_count=retained.eos_evaluation_count,
        hint_count=hint_count,
        hint_reordered_start_count=0,
        initial_guess_count=0,
        accepted_start_count=0,
        distinct_root_count=0,
        final_correction_attempted=False,
        final_correction_accepted=False,
        final_correction_iteration_count=0,
        pre_correction_residual_inf=0.0,
        post_correction_residual_inf=0.0,
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
    temperature_k = _quantity(specification.temperature, "kelvin", "temperature", "solve")
    pressure_pa = _quantity(specification.pressure, "pascal", "pressure", "solve")
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
    native_model = model._native_model(active_parameters)
    roles: list[PhaseRole] = []
    for phase in retained.phases:
        molar_volume = phase.volume_m3 / phase.amount_mol
        compressibility = pressure_pa * molar_volume / (_MOLAR_GAS_CONSTANT * temperature_k)
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
                    _core._audit_phase_state(
                        native_model,
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


@overload
def solve(
    model: epcsaft.Mixture,
    problem: PhaseEquilibriumProblem,
    *,
    active_parameters: epcsaft.ActiveParameterSet | None = None,
) -> PhaseEquilibriumResult: ...


@overload
def solve(
    model: epcsaft.Mixture,
    problem: GeneralReactiveEquilibriumProblem,
    *,
    active_parameters: epcsaft.ActiveParameterSet | None = None,
) -> GeneralReactiveEquilibriumResult: ...


def solve(
    model: epcsaft.Mixture,
    problem: PhaseEquilibriumProblem | GeneralReactiveEquilibriumProblem,
    *,
    active_parameters: epcsaft.ActiveParameterSet | None = None,
) -> PhaseEquilibriumResult | GeneralReactiveEquilibriumResult:
    """Solve one typed nonreactive or GREPE equilibrium specification."""

    if not isinstance(model, epcsaft.Mixture):
        raise TypeError("solve requires an epcsaft.Mixture")
    from ._grepe import GeneralReactiveEquilibriumProblem, _solve_general_reactive_equilibrium

    if isinstance(problem, GeneralReactiveEquilibriumProblem):
        return _solve_general_reactive_equilibrium(
            model, problem, active_parameters=active_parameters
        )
    if type(problem) is not PhaseEquilibriumProblem:
        raise TypeError(
            "problem must be a PhaseEquilibriumProblem or GeneralReactiveEquilibriumProblem"
        )
    if active_parameters is not None and type(active_parameters) is not epcsaft.ActiveParameterSet:
        raise TypeError("active_parameters must be an epcsaft.ActiveParameterSet")
    specification = problem.specification
    if isinstance(specification, FixedTPFeed):
        if specification.topology is None:
            if problem.hints or problem.initial_guesses:
                raise ValueError("phase hints and initial guesses require a declared topology")
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
        if problem.hints and problem.initial_guesses:
            raise ValueError("phase hints and initial_guesses are mutually exclusive")
        temperature_k = _quantity(specification.temperature, "kelvin", "temperature", "solve")
        pressure_pa = _quantity(specification.pressure, "pascal", "pressure", "solve")
        feed = _phase_feed(specification.feed_amounts_mol, len(model.component_ids))
        hints = _validated_phase_hints(
            problem.hints,
            len(topology.phases),
            tuple(len(model.component_ids) for _ in topology.phases),
        )
        initial_guesses = _validated_phase_initial_guesses(
            problem.initial_guesses,
            len(topology.phases),
            tuple(len(model.component_ids) for _ in topology.phases),
        )
        native = cast(
            Mapping[str, object],
            _core._solve_declared_fixed_tp(
                model._native_model(active_parameters),
                temperature_k,
                pressure_pa,
                feed,
                tuple(phase.role for phase in topology.phases),
                tuple(_phase_support_name(phase.support) for phase in topology.phases),
                hints[0],
                hints[1],
                hints[2],
                initial_guesses[0],
                initial_guesses[1],
                initial_guesses[2],
                initial_guesses[3],
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
    if isinstance(specification, FixedPCoexistence):
        return _fixed_p_coexistence(model, problem, specification, active_parameters)
    if not isinstance(specification, FixedTCoexistence):
        raise TypeError("phase-equilibrium specification type is unsupported")
    if problem.initial_guesses:
        raise ValueError("declared phase initial guesses currently require fixed T,P")
    result, directions, _, _, _ = _fixed_t_coexistence_core(
        model,
        problem,
        active_parameters=active_parameters,
        direction_identities=(),
        reference_composition_directions=(),
    )
    if directions is not None:
        raise ValueError("public coexistence unexpectedly returned direction sensitivities")
    return result


def _fixed_t_coexistence_with_reference_directions(
    model: epcsaft.Mixture,
    problem: PhaseEquilibriumProblem,
    *,
    active_parameters: epcsaft.ActiveParameterSet | None,
    direction_identities: tuple[str, ...],
    reference_composition_directions: tuple[tuple[float, ...], ...],
) -> tuple[
    PhaseEquilibriumResult,
    _ReferenceDirectionSensitivity,
    str,
    tuple[str, ...],
]:
    """Solve coexistence and differentiate exact admissible liquid-composition directions."""

    result, directions, _, partner_fingerprint, partner_ids = _fixed_t_coexistence_core(
        model,
        problem,
        active_parameters=active_parameters,
        direction_identities=direction_identities,
        reference_composition_directions=reference_composition_directions,
    )
    if directions is None:
        raise TypeError("reference-direction coexistence returned no direction sensitivity")
    return result, directions, partner_fingerprint, partner_ids


def _fixed_t_coexistence_core(
    model: epcsaft.Mixture,
    problem: PhaseEquilibriumProblem,
    *,
    active_parameters: epcsaft.ActiveParameterSet | None,
    direction_identities: tuple[str, ...],
    reference_composition_directions: tuple[tuple[float, ...], ...],
) -> tuple[
    PhaseEquilibriumResult,
    _ReferenceDirectionSensitivity | None,
    _ReferenceDirectionSensitivity,
    str,
    tuple[str, ...],
]:
    """One fixed-T setup/result boundary with optional exact composition directions."""

    if not isinstance(model, epcsaft.Mixture):
        raise TypeError("reference-direction coexistence requires an epcsaft.Mixture")
    if type(problem) is not PhaseEquilibriumProblem or not isinstance(
        problem.specification, FixedTCoexistence
    ):
        raise TypeError("reference-direction solve requires FixedTCoexistence")
    specification = problem.specification
    topology = specification.topology
    if len(topology.phases) != 2 or specification.reference_phase_index not in {0, 1}:
        raise ValueError("fixed-T coexistence requires exactly two phases and one reference index")
    if len(direction_identities) != len(reference_composition_directions):
        raise ValueError("reference-direction identities and rows must be aligned")
    temperature_k = _quantity(specification.temperature, "kelvin", "temperature", "solve")
    reference = tuple(float(value) for value in specification.reference_mole_fractions)
    if len(reference) != len(model.component_ids):
        raise ValueError("reference composition has the wrong component count")
    reference_support = topology.phases[specification.reference_phase_index].support
    if type(reference_support) is not AllComponents:
        raise ValueError("coexistence reference phase must use AllComponents")
    partner_index = 1 - specification.reference_phase_index
    partner_support = topology.phases[partner_index].support
    if type(partner_support) in {NeutralComponents, DeclaredNeutralComponents}:
        charges = tuple(model._model.charge_number)
        all_neutral_ids = tuple(
            component
            for component, charge in zip(model.component_ids, charges, strict=True)
            if charge == 0
        )
        partner_ids = (
            all_neutral_ids
            if type(partner_support) is NeutralComponents
            else cast(DeclaredNeutralComponents, partner_support).component_ids
        )
        if not partner_ids or any(component not in all_neutral_ids for component in partner_ids):
            raise ValueError("declared vapor support contains a charged or unknown component")
        partner_parameters = model._parameters.select(partner_ids)
        partner_model: epcsaft.Pure | epcsaft.Mixture = (
            epcsaft.Pure(partner_parameters)
            if len(partner_ids) == 1
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
        _core._solve_fixed_t_coexistence(
            model._native_model(active_parameters),
            partner_model._native_model(
                _support_active_parameters(partner_model._parameters, active_parameters)
            ),
            temperature_k,
            tuple(phase.role for phase in topology.phases),
            tuple(_phase_support_name(phase.support) for phase in topology.phases),
            specification.reference_phase_index,
            reference,
            direction_identities,
            reference_composition_directions,
            *hints,
        ),
    )
    diagnostics = _native_phase_diagnostics(cast(Mapping[str, object], native["diagnostics"]))
    if str(native["parameter_fingerprint"]) != model.parameter_fingerprint:
        raise ValueError("native coexistence result has the wrong reference fingerprint")
    if str(native["partner_parameter_fingerprint"]) != partner_model.parameter_fingerprint:
        raise ValueError("native coexistence result has the wrong partner fingerprint")
    if tuple(cast(Sequence[str], native["partner_component_ids"])) != partner_model.component_ids:
        raise ValueError("native coexistence result changed ordered partner support")
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
    result = PhaseEquilibriumResult(
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
    directions = (
        None
        if not direction_identities
        else _reference_direction_sensitivity(
            native.get("reference_direction_sensitivity"), direction_identities
        )
    )
    temperature_directions = _reference_direction_sensitivity(
        native.get("temperature_sensitivity"), ("state/temperature_k",)
    )
    return (
        result,
        directions,
        temperature_directions,
        partner_model.parameter_fingerprint,
        partner_model.component_ids,
    )


def _fixed_p_coexistence_sensitivity(
    fixed: PhaseEquilibriumSensitivity | None,
    temperature: _ReferenceDirectionSensitivity,
) -> PhaseEquilibriumSensitivity | None:
    if fixed is None:
        return None
    if fixed.status != "available" or temperature.status != "available":
        return replace(
            fixed,
            status="unavailable",
            failure_reason="fixed_pressure_temperature_or_parameter_sensitivity_unavailable",
        )
    if fixed.output_identities != temperature.output_identities:
        return replace(
            fixed,
            status="unavailable",
            failure_reason="fixed_pressure_output_layout_changed",
        )
    pressure_temperature = temperature.output_direction_jacobian[0][0]
    if not math.isfinite(pressure_temperature) or pressure_temperature == 0.0:
        return replace(
            fixed,
            status="unavailable",
            failure_reason="fixed_pressure_temperature_closure_singular",
        )
    pressure_parameters = fixed.output_parameter_jacobian[0]
    temperature_parameters = tuple(-value / pressure_temperature for value in pressure_parameters)
    output_rows: list[tuple[float, ...]] = [temperature_parameters]
    for fixed_row, temperature_row in zip(
        fixed.output_parameter_jacobian[1:],
        temperature.output_direction_jacobian[1:],
        strict=True,
    ):
        output_rows.append(
            tuple(
                direct + temperature_row[0] * temperature_parameters[index]
                for index, direct in enumerate(fixed_row)
            )
        )
    return PhaseEquilibriumSensitivity(
        status="available",
        failure_reason=None,
        derivative_order=1,
        higher_order_status="unavailable",
        higher_order_reason="fixed-pressure coexistence provides exact first derivatives only",
        tensor_layout="row-major",
        parameters=fixed.parameters,
        residual_identities=("boundary.pressure_pa_minus_target",),
        state_identities=("boundary.temperature_k",),
        output_identities=("boundary.temperature_k", *fixed.output_identities[1:]),
        residual_state_jacobian=((pressure_temperature,),),
        residual_parameter_jacobian=(pressure_parameters,),
        state_parameter_jacobian=(temperature_parameters,),
        output_parameter_jacobian=tuple(output_rows),
        rank=1,
        dimension=1,
        scaled_condition_number_2=1.0,
        linear_solve_backward_error=0.0,
        reference_phase_index=fixed.reference_phase_index,
    )


def _fixed_p_coexistence(
    model: epcsaft.Mixture,
    problem: PhaseEquilibriumProblem,
    specification: FixedPCoexistence,
    active_parameters: epcsaft.ActiveParameterSet | None,
) -> PhaseEquilibriumResult:
    if problem.initial_guesses:
        raise ValueError("fixed-pressure coexistence accepts phase hints, not finite-phase starts")
    pressure_pa = _quantity(specification.pressure, "pascal", "pressure", "solve")
    if pressure_pa <= 0.0:
        raise ValueError("fixed-pressure coexistence pressure must be positive")
    if len(specification.temperature_bounds) != 2:
        raise ValueError("fixed-pressure coexistence requires two temperature bounds")
    lower = _quantity(specification.temperature_bounds[0], "kelvin", "temperature", "solve")
    upper = _quantity(specification.temperature_bounds[1], "kelvin", "temperature", "solve")
    if not lower < upper:
        raise ValueError("fixed-pressure coexistence temperature bounds must be ordered")
    intervals = specification.packing_fraction_intervals
    if intervals is not None and (
        len(intervals) != 2
        or any(
            len(interval) != 2 or not 0.0 <= interval[0] < interval[1] < 1.0
            for interval in intervals
        )
    ):
        raise ValueError(
            "fixed-pressure coexistence requires one physical packing interval per phase"
        )

    evaluation_count = 0
    eos_evaluations = 0

    def evaluate(
        temperature_k: float,
    ) -> tuple[
        PhaseEquilibriumResult,
        _ReferenceDirectionSensitivity,
        float,
    ]:
        nonlocal evaluation_count, eos_evaluations
        fixed_problem = PhaseEquilibriumProblem(
            FixedTCoexistence(
                temperature_k * epcsaft.unit_registry.kelvin,
                specification.topology,
                specification.reference_phase_index,
                specification.reference_mole_fractions,
            ),
            problem.hints,
        )
        solved, _, temperature_direction, _, _ = _fixed_t_coexistence_core(
            model,
            fixed_problem,
            active_parameters=active_parameters,
            direction_identities=(),
            reference_composition_directions=(),
        )
        if not isinstance(solved.extent, IncipientBoundary):
            raise TypeError("fixed-temperature coexistence returned the wrong extent")
        evaluation_count += 1
        eos_evaluations += solved.diagnostics.eos_evaluation_count
        return solved, temperature_direction, solved.extent.pressure_pa - pressure_pa

    initial = (
        None
        if specification.temperature_initial is None
        else _quantity(
            specification.temperature_initial,
            "kelvin",
            "temperature initial",
            "solve",
        )
    )
    if initial is not None and not lower < initial < upper:
        raise ValueError(
            "fixed-pressure coexistence temperature initial must lie inside its bounds"
        )

    def evaluate_scalar(
        temperature_k: float,
    ) -> tuple[tuple[PhaseEquilibriumResult, _ReferenceDirectionSensitivity], float, float]:
        solved, direction, residual = evaluate(temperature_k)
        return (solved, direction), residual, direction.output_direction_jacobian[0][0]

    try:
        temperature_k, selected_state, residual = _solve_bracketed_scalar(
            evaluate_scalar,
            lower,
            upper,
            initial=initial,
            tolerance=1.0e-10 * max(1.0, pressure_pa),
            maximum_iterations=64,
        )
    except _ScalarBoundarySolveError as error:
        solved, _ = cast(tuple[PhaseEquilibriumResult, _ReferenceDirectionSensitivity], error.state)
        failure_reason = (
            "fixed_pressure_temperature_not_bracketed"
            if error.reason == "not_bracketed"
            else "fixed_pressure_temperature_iteration_limit"
        )
        diagnostics = replace(
            solved.diagnostics,
            route="held_fixed_p_coexistence",
            solver_status="failed",
            numerical_status="not_certified",
            physical_status="rejected",
            failure_stage="intensive_boundary",
            failure_reason=failure_reason,
            certificate_scope="fixed_p_declared_coexistence",
            search_budget=(
                "declared_temperature_bracket_endpoints"
                if error.reason == "not_bracketed"
                else "declared_temperature_bracket;64_iterations"
            ),
            optimizer_iteration_count=(0 if error.reason == "not_bracketed" else 64),
            eos_evaluation_count=eos_evaluations,
        )
        raise PhaseEquilibriumError(
            diagnostics.failure_reason or "fixed-pressure boundary did not converge",
            diagnostics,
            solved.phases,
        )

    solved, temperature_direction = selected_state
    if intervals is not None and any(
        not lower_packing < phase.packing_fraction < upper_packing
        for phase, (lower_packing, upper_packing) in zip(solved.phases, intervals, strict=True)
    ):
        diagnostics = replace(
            solved.diagnostics,
            route="held_fixed_p_coexistence",
            solver_status="converged",
            numerical_status="certified",
            physical_status="rejected",
            failure_stage="eos_domain",
            failure_reason="eos_domain_rejection",
            certificate_scope="fixed_p_declared_coexistence",
            eos_evaluation_count=eos_evaluations,
        )
        raise PhaseEquilibriumError(
            "fixed-pressure coexistence root is outside a declared EOS domain",
            diagnostics,
            solved.phases,
        )
    if temperature_direction.status != "available":
        diagnostics = replace(
            solved.diagnostics,
            route="held_fixed_p_coexistence",
            solver_status="converged",
            numerical_status="uncertified",
            physical_status="rejected",
            failure_stage="intensive_boundary_sensitivity",
            failure_reason=temperature_direction.failure_reason,
            certificate_scope="fixed_p_declared_coexistence",
            eos_evaluation_count=eos_evaluations,
        )
        raise PhaseEquilibriumError(
            diagnostics.failure_reason or "fixed-pressure derivative is unavailable",
            diagnostics,
            solved.phases,
        )
    return replace(
        solved,
        problem=problem,
        extent=IncipientTemperatureBoundary(
            specification.reference_phase_index,
            pressure_pa,
            temperature_k,
        ),
        diagnostics=replace(
            solved.diagnostics,
            route="held_fixed_p_coexistence",
            certificate_scope="fixed_p_declared_coexistence",
            search_scope="declared_temperature_bracket",
            search_budget=f"2_endpoints;64_iteration_limit;{evaluation_count}_evaluations",
            optimizer_solve_count=evaluation_count,
            optimizer_iteration_count=max(0, evaluation_count - 2),
            eos_evaluation_count=eos_evaluations,
            pressure_residual_inf=abs(residual) / max(1.0, pressure_pa),
        ),
        sensitivity=_fixed_p_coexistence_sensitivity(solved.sensitivity, temperature_direction),
    )


@lru_cache(maxsize=128)
def _source_reference_transfer(
    model: epcsaft.Mixture,
    standard_state: ChemicalStandardState,
    *,
    species_ids: tuple[str, ...],
    temperature_k: float,
    pressure_pa: float,
    sensitivities_requested: bool,
    active_parameters: epcsaft.ActiveParameterSet | None,
) -> epcsaft.SourceReferenceTransfer:
    if standard_state.source_reference_component_ids != species_ids:
        raise ValueError("source-reference component order does not match the reaction system")
    if tuple(model.component_ids) != species_ids:
        raise ValueError("source-reference component order does not match the installed EOS")
    if len(standard_state.log_activity_scale_factors) != len(species_ids):
        raise ValueError("source activity-scale factors do not match the species order")
    required_derivatives = (
        ()
        if not sensitivities_requested
        else (("pressure", "parameters") if active_parameters is not None else ("pressure",))
    )
    transfer = epcsaft.source_reference_transfer(
        model,
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
            pure_components=standard_state.source_reference_pure_components,
        ),
        T=temperature_k * epcsaft.unit_registry.kelvin,
        P=pressure_pa * epcsaft.unit_registry.pascal,
        active_parameters=active_parameters,
    )
    if (
        transfer.status != "ok"
        or transfer.domain_status != "admitted-domain"
        or transfer.convergence_status != "converged-limit"
        or transfer.component_ids != species_ids
        or transfer.reference_state_id != standard_state.id
        or transfer.activity_convention_id != standard_state.source_reference_activity_convention_id
        or transfer.parameter_fingerprint != model.parameter_fingerprint
        or transfer.temperature_k != temperature_k
        or transfer.pressure_pa != pressure_pa
        or transfer.source_reference_pressure_pa != standard_state.reference_pressure_pa
    ):
        raise ValueError("installed EOS source-reference receipt is incompatible")
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
        "pressure_derivatives_per_pa": transfer.pressure_derivatives_per_pa,
        "active_parameter_identities": transfer.active_parameter_identities,
        "parameter_derivatives": tuple(
            value for row in transfer.parameter_derivatives for value in row
        ),
        "temperature_k": transfer.temperature_k,
        "pressure_pa": transfer.pressure_pa,
        "source_reference_pressure_pa": transfer.source_reference_pressure_pa,
        "standard_molality_mol_per_kg": (transfer.standard_molality_mol_per_kg),
        "reference_convergence_error": transfer.reference_convergence_error,
        "source_temperature_min_k": transfer.source_temperature_interval_k[0],
        "source_temperature_max_k": transfer.source_temperature_interval_k[1],
        "source_pressure_min_pa": transfer.source_pressure_interval_pa[0],
        "source_pressure_max_pa": transfer.source_pressure_interval_pa[1],
        "candidate_temperature_min_k": transfer.candidate_temperature_interval_k[0],
        "candidate_temperature_max_k": transfer.candidate_temperature_interval_k[1],
        "candidate_pressure_min_pa": transfer.candidate_pressure_interval_pa[0],
        "candidate_pressure_max_pa": transfer.candidate_pressure_interval_pa[1],
        "candidate_qualification": transfer.candidate_qualification,
        "candidate_source_fingerprint": transfer.candidate_source_fingerprint,
        "parameter_fingerprint": transfer.parameter_fingerprint,
        "topology_fingerprint": transfer.topology_fingerprint,
        "component_order_fingerprint": transfer.component_order_fingerprint,
        "reference_state_fingerprint": transfer.reference_state_fingerprint,
        "domain_fingerprint": transfer.domain_fingerprint,
        "artifact_fingerprint": transfer.artifact_fingerprint,
        "helmholtz_basis_id": transfer.helmholtz_basis_id,
    }
