from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any, Literal, cast

from pint import Quantity

from .. import _core
from .._trials import NON_EVALUABLE_TRIAL_CODES, NonEvaluableTrialCode
from ..eos import Mixture, _python_native_artifact_fingerprint
from ..equilibrium import (
    ChemicalEquilibriumConstant,
    ChemicalEquilibriumProblem,
    ChemicalStandardState,
    HomogeneousReactiveObservationDescriptor,
    HomogeneousReactiveObservationProblem,
    HomogeneousReactiveObservationRow,
    HomogeneousReactiveOwnerObservation,
    IdealGasPhase,
    IdealGasVaporComparison,
    ProviderNonidealVapor,
    ProviderPhase,
    ReactiveBubbleContinuationReference,
    ReactiveBubbleObservationDescriptor,
    ReactiveBubbleObservationRow,
    ReactiveBubbleVLEProblem,
    evaluate_homogeneous_reactive_observations,
    evaluate_reactive_bubble_observations,
    homogeneous_reactive_observation_descriptor,
    reactive_bubble_observation_descriptor,
    serialize_homogeneous_reactive_owner_receipt,
    serialize_reactive_bubble_owner_receipt,
    validate_homogeneous_reactive_owner_receipt,
    validate_reactive_bubble_owner_receipt,
)
from ..parameters import ActiveParameterSet, Parameters
from ..reference import SourceReferenceIonPair
from ..units import ParameterError, unit_registry
from ._fit import (
    _FUNCTION_TOLERANCE,
    _GRADIENT_TOLERANCE,
    _MAXIMUM_ITERATIONS,
    _PARAMETER_TOLERANCE,
    FitPreflightError,
    _package_fingerprint,
)
from ._model import EquilibriumFitObservation
from ._resolve import ResolvedFit
from ._result import (
    FitCoordinateReceipt,
    FitResult,
    FitSpecificationReceipt,
    FitStart,
)

type EquilibriumResidual = Literal[
    "log_ratio",
    "scaled_difference",
    "left_censored",
    "right_censored",
]


def _is_sha256(value: str) -> bool:
    if not value.startswith("sha256:") or len(value) != 71:
        return False
    return all(character in "0123456789abcdef" for character in value[7:])


@dataclass(frozen=True, slots=True)
class _ReactiveCoordinate:
    """One affine optimizer coordinate bound to one or more Engine parameters."""

    identity: str
    parameter_identities: tuple[str, ...]
    unit: str
    scale: float
    lower: float
    upper: float
    origin: float

    def __post_init__(self) -> None:
        object.__setattr__(self, "parameter_identities", tuple(self.parameter_identities))
        if not self.identity or not self.unit or not self.parameter_identities:
            raise ValueError("regression coordinate identity, parameters, and unit are required")
        if len(set(self.parameter_identities)) != len(self.parameter_identities):
            raise ValueError("a regression coordinate cannot repeat a parameter identity")
        if not (
            math.isfinite(self.scale)
            and self.scale > 0.0
            and math.isfinite(self.lower)
            and math.isfinite(self.upper)
            and self.lower < self.upper
            and math.isfinite(self.origin)
            and self.lower <= self.origin <= self.upper
        ):
            raise ValueError("regression coordinate scale, bounds, and origin are inconsistent")

    @property
    def affine_bounds(self) -> tuple[float, float]:
        return (
            (self.lower - self.origin) / self.scale,
            (self.upper - self.origin) / self.scale,
        )


@dataclass(frozen=True, slots=True)
class _ReactiveStart:
    identity: str
    values: tuple[float, ...]

    def __post_init__(self) -> None:
        object.__setattr__(self, "values", tuple(self.values))
        if (
            not self.identity
            or not self.values
            or not all(math.isfinite(value) for value in self.values)
        ):
            raise ValueError("regression starts require an identity and finite physical values")


@dataclass(frozen=True, slots=True)
class _ReactiveTarget:
    """One measurement row and its Regression-owned residual/weight contract.

    ``covariance_identity`` is an immutable grouping/provenance identity, not a
    numeric covariance matrix. Regression retains ownership of weights.
    """

    identity: str
    prediction_identity: str
    observed: float
    unit: str
    residual: EquilibriumResidual
    scale: float
    multiplier: float
    role: str
    measurement_classification: str
    source_identity: str
    source_hash: str
    included: bool = True
    drop_reason: str = ""
    aggregate_identity: str | None = None
    covariance_identity: str | None = None

    def __post_init__(self) -> None:
        required = (
            self.identity,
            self.prediction_identity,
            self.unit,
            self.role,
            self.measurement_classification,
            self.source_identity,
            self.source_hash,
        )
        if any(not value for value in required):
            raise ValueError("regression target identity, basis, role, and provenance are required")
        if not _is_sha256(self.source_hash):
            raise ValueError("regression target source_hash must be a lowercase SHA-256 identity")
        if self.residual not in {
            "log_ratio",
            "scaled_difference",
            "left_censored",
            "right_censored",
        }:
            raise ValueError("regression target residual is unsupported")
        if not (
            math.isfinite(self.observed)
            and math.isfinite(self.scale)
            and self.scale > 0.0
            and math.isfinite(self.multiplier)
            and self.multiplier > 0.0
        ):
            raise ValueError("regression target value, scale, and multiplier must be finite")
        if self.residual == "log_ratio" and self.observed <= 0.0:
            raise ValueError("log-ratio regression targets must be strictly positive")
        if self.included == bool(self.drop_reason):
            raise ValueError("only excluded regression targets declare a drop reason")
        if bool(self.aggregate_identity) != bool(self.covariance_identity):
            raise ValueError(
                "regression aggregate and covariance grouping identities must be paired"
            )


@dataclass(frozen=True, slots=True)
class _ReactiveBatch:
    """Compatible measurements evaluated by one Equilibrium-owned solve."""

    identity: str
    observation_problem: HomogeneousReactiveObservationProblem | ReactiveBubbleVLEProblem
    temperature: Quantity[Any]
    pressure: Quantity[Any] | None
    descriptor: HomogeneousReactiveObservationDescriptor | ReactiveBubbleObservationDescriptor
    targets: tuple[_ReactiveTarget, ...]

    def __post_init__(self) -> None:
        object.__setattr__(self, "targets", tuple(self.targets))
        if not self.identity or not self.targets:
            raise ValueError("reactive regression batch identity and targets are required")
        if len({target.identity for target in self.targets}) != len(self.targets):
            raise ValueError("reactive regression target identities must be unique within a batch")
        if isinstance(self.descriptor, HomogeneousReactiveObservationDescriptor):
            if not isinstance(self.observation_problem, HomogeneousReactiveObservationProblem):
                raise TypeError("homogeneous descriptor has the wrong Equilibrium problem")
            if self.descriptor.problem_identity != self.observation_problem.identity:
                raise ValueError("reactive regression descriptor and problem identities differ")
            if self.descriptor.phase_identity != self.observation_problem.phase_identity:
                raise ValueError("reactive regression descriptor and phase identities differ")
            if self.descriptor.phase_role != self.observation_problem.phase_role:
                raise ValueError("reactive regression descriptor and phase roles differ")
            if (
                self.descriptor.continuation_identity
                != self.observation_problem.continuation_identity
                or self.descriptor.branch_policy != self.observation_problem.branch_policy
            ):
                raise ValueError("reactive regression continuation declaration is inconsistent")
            if self.pressure is None:
                raise ValueError("homogeneous reactive regression requires fixed pressure")
        elif not isinstance(self.observation_problem, ReactiveBubbleVLEProblem):
            raise TypeError("bubble descriptor has the wrong Equilibrium problem")
        rows = {row.identity: row for row in self.descriptor.rows}
        for target in self.targets:
            row = rows.get(target.prediction_identity)
            if row is None:
                raise ValueError("regression target is absent from its Equilibrium descriptor")
            if target.unit != row.unit:
                raise ValueError("regression target unit does not match its Equilibrium row")
            if (
                target.aggregate_identity != row.aggregate_identity
                or target.covariance_identity != row.covariance_identity
            ):
                raise ValueError(
                    "regression target aggregate/covariance identities do not match its row"
                )
            if row.support == "left_censored":
                if target.residual != "left_censored" or target.observed != row.censor_limit:
                    raise ValueError("left-censored target does not match its declared limit")
            elif row.support == "right_censored":
                if target.residual != "right_censored" or target.observed != row.censor_limit:
                    raise ValueError("right-censored target does not match its declared limit")
            elif target.residual in {"left_censored", "right_censored"}:
                raise ValueError("censored residual requires a censored Equilibrium row")
        try:
            float(self.temperature.to("kelvin").magnitude)
            if self.pressure is not None:
                float(self.pressure.to("pascal").magnitude)
        except (AttributeError, TypeError, ValueError) as error:
            raise ValueError("reactive regression state has invalid units") from error


@dataclass(frozen=True, slots=True)
class _ReactiveControls:
    starts: tuple[_ReactiveStart, ...]
    affine_bound_margin: float = 1.0e-7
    maximum_iterations: int = _MAXIMUM_ITERATIONS
    function_tolerance: float = _FUNCTION_TOLERANCE
    gradient_tolerance: float = _GRADIENT_TOLERANCE
    parameter_tolerance: float = _PARAMETER_TOLERANCE

    def __post_init__(self) -> None:
        object.__setattr__(self, "starts", tuple(self.starts))
        if (
            not self.starts
            or len({start.identity for start in self.starts}) != len(self.starts)
            or self.maximum_iterations <= 0
            or not math.isfinite(self.affine_bound_margin)
            or self.affine_bound_margin < 0.0
        ):
            raise ValueError("reactive regression start and iteration controls are invalid")
        tolerances = (
            self.function_tolerance,
            self.gradient_tolerance,
            self.parameter_tolerance,
        )
        if any(not math.isfinite(value) or value <= 0.0 for value in tolerances):
            raise ValueError("reactive regression tolerances must be positive and finite")


@dataclass(frozen=True, slots=True)
class _ResolvedReactiveFit:
    """Optimizer-independent specification for generic Equilibrium observations."""

    parameters: Parameters
    coordinates: tuple[_ReactiveCoordinate, ...]
    batches: tuple[_ReactiveBatch, ...]
    controls: _ReactiveControls

    def __post_init__(self) -> None:
        object.__setattr__(self, "coordinates", tuple(self.coordinates))
        object.__setattr__(self, "batches", tuple(self.batches))
        if not isinstance(self.parameters, Parameters):
            raise TypeError("equilibrium regression parameters must be Parameters")
        if not self.coordinates or not self.batches:
            raise ValueError("equilibrium regression requires coordinates and batches")
        if len({coordinate.identity for coordinate in self.coordinates}) != len(self.coordinates):
            raise ValueError("equilibrium regression coordinate identities must be unique")
        parameter_ids = tuple(
            identity
            for coordinate in self.coordinates
            for identity in coordinate.parameter_identities
        )
        if len(set(parameter_ids)) != len(parameter_ids):
            raise ValueError("each active parameter must bind to exactly one coordinate")
        active = ActiveParameterSet(self.parameters, parameter_ids)
        units = {coordinate.identity: coordinate.unit for coordinate in self.coordinates}
        coordinate_for_parameter = {
            parameter: coordinate.identity
            for coordinate in self.coordinates
            for parameter in coordinate.parameter_identities
        }
        for parameter in active.coordinates:
            if parameter.unit != units[coordinate_for_parameter[parameter.identity]]:
                raise ValueError("tied active parameters must use their canonical common unit")
        coordinate_by_parameter = {
            parameter: coordinate
            for coordinate in self.coordinates
            for parameter in coordinate.parameter_identities
        }
        for label, values in (
            (
                "lower",
                tuple(coordinate_by_parameter[identity].lower for identity in parameter_ids),
            ),
            (
                "upper",
                tuple(coordinate_by_parameter[identity].upper for identity in parameter_ids),
            ),
            (
                "origin",
                tuple(coordinate_by_parameter[identity].origin for identity in parameter_ids),
            ),
        ):
            try:
                ActiveParameterSet(self.parameters, parameter_ids, values=values)
            except ParameterError as error:
                raise ValueError(
                    f"equilibrium regression {label} coordinates leave a parameter domain"
                ) from error
        if len({batch.identity for batch in self.batches}) != len(self.batches):
            raise ValueError("equilibrium regression batch identities must be unique")
        targets = tuple(target.identity for batch in self.batches for target in batch.targets)
        if len(set(targets)) != len(targets):
            raise ValueError("equilibrium regression target identities must be globally unique")
        for batch in self.batches:
            if batch.descriptor.active_parameter_identities != parameter_ids:
                raise ValueError("Equilibrium descriptor does not preserve active parameter order")
        for start in self.controls.starts:
            if len(start.values) != len(self.coordinates):
                raise ValueError("every regression start must cover every physical coordinate")
            if any(
                value < coordinate.lower or value > coordinate.upper
                for value, coordinate in zip(start.values, self.coordinates, strict=True)
            ):
                raise ValueError("an explicit physical regression start lies outside its bounds")
            start_values_by_identity = {
                identity: value
                for coordinate, value in zip(self.coordinates, start.values, strict=True)
                for identity in coordinate.parameter_identities
            }
            try:
                ActiveParameterSet(
                    self.parameters,
                    parameter_ids,
                    values=tuple(start_values_by_identity[item] for item in parameter_ids),
                )
            except ParameterError as error:
                raise ValueError(
                    "an explicit regression start leaves a parameter domain"
                ) from error


def _exact_mapping(
    value: object,
    fields: frozenset[str],
    label: str,
) -> Mapping[str, object]:
    if not isinstance(value, Mapping) or set(value) != fields:
        raise ValueError(f"{label} must contain exactly {tuple(sorted(fields))}")
    return cast(Mapping[str, object], value)


def _sequence(value: object, label: str) -> Sequence[object]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        raise TypeError(f"{label} must be an array")
    return value


def _number(value: object, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise TypeError(f"{label} must be numeric")
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{label} must be finite")
    return result


def _integer(value: object, label: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"{label} must be an integer")
    return value


def _float_tuple(value: object, label: str) -> tuple[float, ...]:
    return tuple(_number(item, label) for item in _sequence(value, label))


def _float_rows(value: object, label: str) -> tuple[tuple[float, ...], ...]:
    return tuple(
        _float_tuple(row, f"{label}[{index}]") for index, row in enumerate(_sequence(value, label))
    )


def _reaction_system(value: object) -> ChemicalEquilibriumProblem:
    item = _exact_mapping(
        value,
        frozenset(
            {
                "species_ids",
                "charges",
                "molar_masses_kg_per_mol",
                "balance_matrix",
                "conserved_totals",
                "reaction_matrix",
                "feed_amounts_mol",
                "equilibrium_constants",
                "strict_interior_amount_floor_mol",
                "source_standard_state",
            }
        ),
        "reactive reaction system",
    )
    constants_list: list[ChemicalEquilibriumConstant] = []
    for raw_record in _sequence(item["equilibrium_constants"], "equilibrium constants"):
        record = _sequence(raw_record, "equilibrium constant")
        if len(record) != 6 or not isinstance(record[5], bool):
            raise ValueError("equilibrium constant requires six typed fields")
        constants_list.append(
            ChemicalEquilibriumConstant(
                _number(record[0], "equilibrium constant"),
                str(record[1]),
                str(record[2]),
                str(record[3]),
                str(record[4]),
                record[5],
            )
        )
    constants = tuple(constants_list)
    standard_value = item["source_standard_state"]
    standard = None
    if standard_value is not None:
        source = _exact_mapping(
            standard_value,
            frozenset(
                {
                    "id",
                    "activity_scale_id",
                    "log_activity_scale_factors",
                    "reference_pressure_pa",
                    "component_ids",
                    "solvent_composition",
                    "ion_pairs",
                    "phase",
                    "reference_convention",
                    "activity_convention_id",
                    "standard_molality_mol_per_kg",
                }
            ),
            "source standard state",
        )
        ion_pair_values: list[SourceReferenceIonPair] = []
        for raw_record in _sequence(source["ion_pairs"], "source ion pairs"):
            record = _sequence(raw_record, "source ion pair")
            if len(record) != 2:
                raise ValueError("source ion pair requires components and stoichiometry")
            components = tuple(
                str(component) for component in _sequence(record[0], "ion components")
            )
            stoichiometry = tuple(
                _integer(coefficient, "ion stoichiometry")
                for coefficient in _sequence(record[1], "ion stoichiometry")
            )
            if len(components) != 2 or len(stoichiometry) != 2:
                raise ValueError("source ion-pair rows require two components")
            ion_pair_values.append(
                SourceReferenceIonPair(
                    components,
                    stoichiometry,
                )
            )
        ion_pairs = tuple(ion_pair_values)
        standard = ChemicalStandardState(
            str(source["id"]),
            str(source["activity_scale_id"]),
            _float_tuple(source["log_activity_scale_factors"], "activity scale factors"),
            _number(source["reference_pressure_pa"], "reference pressure"),
            tuple(str(item) for item in _sequence(source["component_ids"], "reference components")),
            _float_tuple(source["solvent_composition"], "reference solvent composition"),
            ion_pairs,
            cast(Any, source["phase"]),
            cast(Any, source["reference_convention"]),
            str(source["activity_convention_id"]),
            _number(source["standard_molality_mol_per_kg"], "standard molality"),
        )
    return ChemicalEquilibriumProblem(
        tuple(str(item) for item in _sequence(item["species_ids"], "species identities")),
        tuple(_integer(value, "charge") for value in _sequence(item["charges"], "charges")),
        _float_tuple(item["molar_masses_kg_per_mol"], "molar masses"),
        _float_rows(item["balance_matrix"], "balance matrix"),
        _float_tuple(item["conserved_totals"], "conserved totals"),
        _float_rows(item["reaction_matrix"], "reaction matrix"),
        _float_tuple(item["feed_amounts_mol"], "feed amounts"),
        constants,
        _number(item["strict_interior_amount_floor_mol"], "interior amount floor"),
        standard,
    )


def _homogeneous_problem_from_mapping(
    identity: str,
    phase_value: object,
    continuation_value: object,
    reaction_value: object,
    resolved: ResolvedFit,
) -> HomogeneousReactiveObservationProblem:
    phase = _exact_mapping(
        phase_value,
        frozenset({"identity", "role", "model", "admissible_packing_fraction_interval"}),
        "reactive phase",
    )
    model_kind = str(phase["model"])
    phase_model: IdealGasPhase | ProviderPhase
    if model_kind == "provider":
        bounds = _float_tuple(phase["admissible_packing_fraction_interval"], "packing interval")
        if len(bounds) != 2:
            raise ValueError("packing interval requires two bounds")
        phase_model = ProviderPhase(Mixture(resolved.parameters), cast(Any, bounds))
    elif model_kind == "ideal_gas":
        phase_model = IdealGasPhase(resolved.parameters.fingerprint, "parameter-fit-ideal-gas")
    else:
        raise ValueError("reactive phase model must be provider or ideal_gas")
    continuation = _exact_mapping(
        continuation_value,
        frozenset({"identity", "branch_policy", "reference"}),
        "reactive continuation",
    )
    raw_reference = continuation["reference"]
    reference = None
    if raw_reference is not None:
        from ..equilibrium import HomogeneousContinuationReference

        reference_item = _exact_mapping(
            raw_reference,
            frozenset(
                {
                    "mole_fractions",
                    "molar_volume_m3_per_mol",
                    "maximum_log_composition_distance",
                    "maximum_log_volume_distance",
                    "phase_identity",
                    "phase_role",
                    "anchor_state_identity",
                    "anchor_phase_role_receipt_fingerprint",
                    "anchor_parameter_fingerprint",
                    "anchor_reaction_fingerprint",
                    "anchor_topology_fingerprint",
                    "anchor_basin_states",
                }
            ),
            "reactive continuation reference",
        )
        reference = HomogeneousContinuationReference(
            _float_tuple(reference_item["mole_fractions"], "reference composition"),
            _number(reference_item["molar_volume_m3_per_mol"], "reference molar volume"),
            _number(
                reference_item["maximum_log_composition_distance"],
                "maximum composition distance",
            ),
            _number(
                reference_item["maximum_log_volume_distance"],
                "maximum volume distance",
            ),
            str(reference_item["phase_identity"]),
            cast(Any, reference_item["phase_role"]),
            str(reference_item["anchor_state_identity"]),
            str(reference_item["anchor_phase_role_receipt_fingerprint"]),
            str(reference_item["anchor_parameter_fingerprint"]),
            str(reference_item["anchor_reaction_fingerprint"]),
            str(reference_item["anchor_topology_fingerprint"]),
            tuple(
                (
                    _float_tuple(
                        cast(Mapping[str, object], record)["mole_fractions"],
                        "anchor basin composition",
                    ),
                    _number(
                        cast(Mapping[str, object], record)[
                            "molar_volume_m3_per_mol"
                        ],
                        "anchor basin molar volume",
                    ),
                )
                for record in cast(
                    Sequence[object], reference_item["anchor_basin_states"]
                )
            ),
        )
    return HomogeneousReactiveObservationProblem(
        identity,
        str(phase["identity"]),
        cast(Any, phase["role"]),
        phase_model,
        _reaction_system(reaction_value),
        str(continuation["identity"]),
        cast(Any, continuation["branch_policy"]),
        reference,
    )


_TARGET_FIELDS = frozenset(
    {
        "identity",
        "prediction_identity",
        "observed",
        "unit",
        "residual",
        "scale",
        "multiplier",
        "role",
        "measurement_classification",
        "source_identity",
        "source_hash",
        "included",
        "drop_reason",
        "aggregate_identity",
        "covariance_identity",
    }
)


def _targets_from_mapping(value: object) -> tuple[_ReactiveTarget, ...]:
    return tuple(
        _ReactiveTarget(
            str(record["identity"]),
            str(record["prediction_identity"]),
            _number(record["observed"], "observed value"),
            str(record["unit"]),
            cast(Any, record["residual"]),
            _number(record["scale"], "residual scale"),
            _number(record["multiplier"], "target multiplier"),
            str(record["role"]),
            str(record["measurement_classification"]),
            str(record["source_identity"]),
            str(record["source_hash"]),
            bool(record["included"]),
            "" if record["drop_reason"] is None else str(record["drop_reason"]),
            None if record["aggregate_identity"] is None else str(record["aggregate_identity"]),
            None if record["covariance_identity"] is None else str(record["covariance_identity"]),
        )
        for record in (
            _exact_mapping(raw, _TARGET_FIELDS, "reactive regression target")
            for raw in _sequence(value, "reactive targets")
        )
    )


def _reactive_bubble_batch(
    item: Mapping[str, object],
    resolved: ResolvedFit,
    active: ActiveParameterSet,
) -> _ReactiveBatch:
    identity = str(item["identity"])
    liquid = _homogeneous_problem_from_mapping(
        str(item["liquid_problem_identity"]),
        item["liquid_phase"],
        item["liquid_continuation"],
        item["reaction_system"],
        resolved,
    )
    vapor = _exact_mapping(
        item["vapor"],
        frozenset({"identity", "component_ids", "model", "reference_id"}),
        "reactive bubble vapor",
    )
    model_kind = str(vapor["model"])
    vapor_model: ProviderNonidealVapor | IdealGasVaporComparison
    if model_kind == "provider_nonideal":
        vapor_model = ProviderNonidealVapor(str(vapor["reference_id"]))
    elif model_kind == "ideal_gas_comparison":
        vapor_model = IdealGasVaporComparison(str(vapor["reference_id"]))
    else:
        raise ValueError("reactive bubble vapor model is unsupported")
    continuation = _exact_mapping(
        item["bubble_continuation"],
        frozenset({"identity", "reference"}),
        "reactive bubble continuation",
    )
    raw_reference = continuation["reference"]
    reference = None
    if raw_reference is not None:
        record = _exact_mapping(
            raw_reference,
            frozenset(
                {
                    "pressure_pa",
                    "vapor_mole_fractions",
                    "liquid_molar_volume_m3_per_mol",
                    "vapor_molar_volume_m3_per_mol",
                    "maximum_log_pressure_distance",
                    "maximum_log_composition_distance",
                    "maximum_log_volume_distance",
                    "anchor_state_identity",
                    "anchor_problem_fingerprint",
                    "anchor_topology_fingerprint",
                    "anchor_parameter_fingerprint",
                }
            ),
            "reactive bubble continuation reference",
        )
        reference = ReactiveBubbleContinuationReference(
            _number(record["pressure_pa"], "bubble reference pressure"),
            _float_tuple(record["vapor_mole_fractions"], "bubble reference vapor composition"),
            _number(record["liquid_molar_volume_m3_per_mol"], "bubble liquid molar volume"),
            _number(record["vapor_molar_volume_m3_per_mol"], "bubble vapor molar volume"),
            _number(record["maximum_log_pressure_distance"], "bubble pressure distance"),
            _number(record["maximum_log_composition_distance"], "bubble composition distance"),
            _number(record["maximum_log_volume_distance"], "bubble volume distance"),
            str(record["anchor_state_identity"]),
            str(record["anchor_problem_fingerprint"]),
            str(record["anchor_topology_fingerprint"]),
            str(record["anchor_parameter_fingerprint"]),
        )
    interval = _float_tuple(item["pressure_interval_pa"], "bubble pressure interval")
    if len(interval) != 2:
        raise ValueError("bubble pressure interval requires two bounds")
    problem = ReactiveBubbleVLEProblem(
        identity,
        liquid,
        str(vapor["identity"]),
        tuple(str(value) for value in _sequence(vapor["component_ids"], "vapor components")),
        vapor_model,
        interval,
        _float_tuple(item["pressure_starts_pa"], "bubble pressure starts"),
        str(continuation["identity"]),
        reference,
    )
    rows = tuple(
        ReactiveBubbleObservationRow(
            identity=str(record["identity"]),
            quantity=cast(Any, record["quantity"]),
            coefficients=_float_tuple(record["coefficients"], "bubble row coefficients"),
            solvent_mass_coefficients_kg_per_mol=_float_tuple(
                record["solvent_mass_coefficients_kg_per_mol"], "bubble solvent mass basis"
            ),
            support=cast(Any, record["support"]),
            censor_limit=(
                None
                if record["censor_limit"] is None
                else _number(record["censor_limit"], "bubble censor limit")
            ),
            aggregate_identity=(
                None if record["aggregate_identity"] is None else str(record["aggregate_identity"])
            ),
            covariance_identity=(
                None
                if record["covariance_identity"] is None
                else str(record["covariance_identity"])
            ),
        )
        for record in (
            _exact_mapping(
                raw,
                frozenset(
                    {
                        "identity",
                        "quantity",
                        "coefficients",
                        "solvent_mass_coefficients_kg_per_mol",
                        "support",
                        "censor_limit",
                        "aggregate_identity",
                        "covariance_identity",
                    }
                ),
                "reactive bubble observation row",
            )
            for raw in _sequence(item["rows"], "reactive bubble rows")
        )
    )
    descriptor = reactive_bubble_observation_descriptor(problem, rows, active_parameters=active)
    return _ReactiveBatch(
        identity,
        problem,
        _number(item["temperature_k"], "temperature") * unit_registry.kelvin,
        None,
        descriptor,
        _targets_from_mapping(item["targets"]),
    )


def _reactive_fit_problem(resolved: ResolvedFit) -> _ResolvedReactiveFit:
    observations = tuple(
        item for item in resolved.fit.observations if isinstance(item, EquilibriumFitObservation)
    )
    if len(observations) != len(resolved.fit.observations):
        raise TypeError("reactive Parameter Fit cannot mix owner families")
    active = ActiveParameterSet(resolved.parameters, resolved.selected_identities)
    batches: list[_ReactiveBatch] = []
    for specification in observations:
        item = specification.to_mapping()
        if specification.family == "reactive_bubble_vle":
            batches.append(_reactive_bubble_batch(item, resolved, active))
            continue
        problem = _homogeneous_problem_from_mapping(
            str(item["identity"]),
            item["phase"],
            item["continuation"],
            item["reaction_system"],
            resolved,
        )
        rows = tuple(
            HomogeneousReactiveObservationRow(
                identity=str(record["identity"]),
                quantity=cast(Any, record["quantity"]),
                coefficients=_float_tuple(record["coefficients"], "row coefficients"),
                solvent_mass_coefficients_kg_per_mol=_float_tuple(
                    record["solvent_mass_coefficients_kg_per_mol"], "solvent mass basis"
                ),
                support=cast(Any, record["support"]),
                censor_limit=(
                    None
                    if record["censor_limit"] is None
                    else _number(record["censor_limit"], "censor limit")
                ),
                aggregate_identity=(
                    None
                    if record["aggregate_identity"] is None
                    else str(record["aggregate_identity"])
                ),
                covariance_identity=(
                    None
                    if record["covariance_identity"] is None
                    else str(record["covariance_identity"])
                ),
            )
            for record in (
                _exact_mapping(
                    raw,
                    frozenset(
                        {
                            "identity",
                            "quantity",
                            "coefficients",
                            "solvent_mass_coefficients_kg_per_mol",
                            "support",
                            "censor_limit",
                            "aggregate_identity",
                            "covariance_identity",
                        }
                    ),
                    "reactive observation row",
                )
                for raw in _sequence(item["rows"], "reactive rows")
            )
        )
        descriptor = homogeneous_reactive_observation_descriptor(
            problem,
            rows,
            active_parameters=active,
        )
        batches.append(
            _ReactiveBatch(
                str(item["identity"]),
                problem,
                _number(item["temperature_k"], "temperature") * unit_registry.kelvin,
                _number(item["pressure_pa"], "pressure") * unit_registry.pascal,
                descriptor,
                _targets_from_mapping(item["targets"]),
            )
        )
    return _ResolvedReactiveFit(
        resolved.parameters,
        tuple(
            _ReactiveCoordinate(
                coordinate.identity,
                coordinate.tied_identities,
                coordinate.unit,
                coordinate.scale,
                coordinate.bounds[0],
                coordinate.bounds[1],
                coordinate.origin,
            )
            for coordinate in resolved.coordinates
        ),
        tuple(batches),
        _ReactiveControls(
            starts=tuple(
                _ReactiveStart(start.identity, start.values)
                for start in resolved.fit.solver.multistart.starts
            ),
            affine_bound_margin=resolved.fit.solver.multistart.affine_bound_margin,
        ),
    )


def _validated_owner_receipts(
    problem: _ResolvedReactiveFit,
    start: FitStart,
) -> tuple[Mapping[str, Any], ...]:
    if not start.usable:
        if start.owner_receipts:
            raise ValueError("unusable equilibrium start must not retain owner receipts")
        return ()
    batches = tuple(
        batch for batch in problem.batches if any(target.included for target in batch.targets)
    )
    if len(start.owner_receipts) != len(batches):
        raise ValueError("usable equilibrium start has incomplete owner batch receipts")
    parameter_ids = tuple(
        identity
        for coordinate in problem.coordinates
        for identity in coordinate.parameter_identities
    )
    coordinate_values = {
        coordinate.identity: value
        for coordinate, value in zip(problem.coordinates, start.final_physical, strict=True)
    }
    coordinate_by_parameter = {
        parameter: coordinate.identity
        for coordinate in problem.coordinates
        for parameter in coordinate.parameter_identities
    }
    active_values = tuple(
        coordinate_values[coordinate_by_parameter[identity]] for identity in parameter_ids
    )
    active = ActiveParameterSet(problem.parameters, parameter_ids, values=active_values)
    validated: list[Mapping[str, Any]] = []
    for batch, opaque_receipt in zip(batches, start.owner_receipts, strict=True):
        receipt_json = opaque_receipt.canonical_bytes_utf8
        if isinstance(batch.descriptor, ReactiveBubbleObservationDescriptor):
            receipt_digest = validate_reactive_bubble_owner_receipt(
                receipt_json,
                descriptor=batch.descriptor,
                batch_identity=batch.identity,
                requested_temperature=batch.temperature,
                active_parameters=active,
            )
        else:
            observations = tuple(
                HomogeneousReactiveOwnerObservation(
                    target.identity,
                    target.prediction_identity,
                    target.aggregate_identity,
                    target.covariance_identity,
                )
                for target in batch.targets
                if target.included
            )
            assert batch.pressure is not None
            receipt_digest = validate_homogeneous_reactive_owner_receipt(
                receipt_json,
                descriptor=batch.descriptor,
                batch_identity=batch.identity,
                requested_temperature=batch.temperature,
                requested_pressure=batch.pressure,
                active_parameters=active,
                observations=observations,
            )
        if receipt_digest != opaque_receipt.digest:
            raise ValueError("Equilibrium owner receipt digest does not match its opaque bytes")
        validated.append(
            {
                "canonical_bytes_utf8": receipt_json,
                "digest": receipt_digest,
            }
        )
    return tuple(validated)


def _target_payload(target: _ReactiveTarget) -> dict[str, object]:
    return {
        "identity": target.identity,
        "prediction_identity": target.prediction_identity,
        "observed": target.observed,
        "unit": target.unit,
        "residual": target.residual,
        "scale": target.scale,
        "multiplier": target.multiplier,
        "role": target.role,
        "measurement_classification": target.measurement_classification,
        "source_identity": target.source_identity,
        "source_hash": target.source_hash,
        "included": target.included,
        "drop_reason": target.drop_reason,
        "aggregate_identity": target.aggregate_identity,
        "covariance_identity": target.covariance_identity,
    }


def _batch_payload(batch: _ReactiveBatch) -> dict[str, object]:
    common = {
        "identity": batch.identity,
        "family": batch.descriptor.identity,
        "descriptor_fingerprint": batch.descriptor.fingerprint,
        "problem_fingerprint": batch.descriptor.problem_fingerprint,
        "reaction_fingerprint": batch.descriptor.reaction_fingerprint,
        "topology_fingerprint": batch.descriptor.topology_fingerprint,
        "reference_fingerprint": batch.descriptor.reference_fingerprint,
        "domain_fingerprint": batch.descriptor.domain_fingerprint,
        "temperature_k": float(batch.temperature.to("kelvin").magnitude),
        "targets": tuple(_target_payload(target) for target in batch.targets),
    }
    if isinstance(batch.observation_problem, ReactiveBubbleVLEProblem):
        return {
            **common,
            "observation_problem_identity": batch.observation_problem.identity,
            "liquid_problem_identity": batch.observation_problem.liquid_problem.identity,
            "vapor_phase_identity": batch.observation_problem.vapor_phase_identity,
            "vapor_component_ids": batch.observation_problem.vapor_component_ids,
            "continuation_identity": batch.observation_problem.continuation_identity,
            "continuation_reference": (
                None
                if batch.observation_problem.continuation_reference is None
                else batch.observation_problem.continuation_reference.to_mapping()
            ),
            "pressure_interval_pa": batch.observation_problem.pressure_interval_pa,
            "pressure_starts_pa": batch.observation_problem.pressure_starts_pa,
        }
    assert batch.pressure is not None
    return {
        **common,
        "observation_problem_identity": batch.observation_problem.identity,
        "phase_identity": batch.observation_problem.phase_identity,
        "phase_role": batch.observation_problem.phase_role,
        "continuation_identity": batch.observation_problem.continuation_identity,
        "branch_policy": batch.observation_problem.branch_policy,
        "continuation_reference": (
            None
            if batch.observation_problem.continuation_reference is None
            else batch.observation_problem.continuation_reference.to_mapping()
        ),
        "pressure_pa": float(batch.pressure.to("pascal").magnitude),
    }


def _problem_payload(problem: _ResolvedReactiveFit) -> dict[str, object]:
    return {
        "parameter_fingerprint": problem.parameters.fingerprint,
        "coordinates": tuple(
            {
                "identity": coordinate.identity,
                "parameter_identities": coordinate.parameter_identities,
                "unit": coordinate.unit,
                "scale": coordinate.scale,
                "lower": coordinate.lower,
                "upper": coordinate.upper,
                "origin": coordinate.origin,
                "affine_bounds": coordinate.affine_bounds,
            }
            for coordinate in problem.coordinates
        ),
        "batches": tuple(_batch_payload(batch) for batch in problem.batches),
        "controls": {
            "starts": tuple(
                {"identity": start.identity, "values": start.values}
                for start in problem.controls.starts
            ),
            "affine_bound_margin": problem.controls.affine_bound_margin,
            "maximum_iterations": problem.controls.maximum_iterations,
            "function_tolerance": problem.controls.function_tolerance,
            "gradient_tolerance": problem.controls.gradient_tolerance,
            "parameter_tolerance": problem.controls.parameter_tolerance,
        },
    }


def _specification_receipt(
    problem: _ResolvedReactiveFit,
    canonical_fit_mapping: Mapping[str, object] | None = None,
) -> FitSpecificationReceipt:
    active_batches = tuple(
        batch for batch in problem.batches if any(target.included for target in batch.targets)
    )
    return FitSpecificationReceipt.create(
        _problem_payload(problem) if canonical_fit_mapping is None else canonical_fit_mapping,
        parameters=problem.parameters,
        coordinates=tuple(
            FitCoordinateReceipt(
                coordinate.identity,
                coordinate.parameter_identities,
                coordinate.unit,
                coordinate.origin,
                coordinate.scale,
            )
            for coordinate in problem.coordinates
        ),
        row_identities=tuple(
            target.identity
            for batch in active_batches
            for target in batch.targets
            if target.included
        ),
        block_sizes=tuple(
            sum(target.included for target in batch.targets) for batch in active_batches
        ),
        input_row_count=sum(len(batch.targets) for batch in problem.batches),
        dropped_row_count=sum(
            not target.included for batch in problem.batches for target in batch.targets
        ),
        allows_parameter_promotion=False,
    )


def _active_trial(
    problem: _ResolvedReactiveFit,
    values: Sequence[float],
) -> ActiveParameterSet:
    active_ids = tuple(
        identity
        for coordinate in problem.coordinates
        for identity in coordinate.parameter_identities
    )
    physical_by_parameter = {
        identity: value
        for coordinate, value in zip(problem.coordinates, values, strict=True)
        for identity in coordinate.parameter_identities
    }
    return ActiveParameterSet(
        problem.parameters,
        active_ids,
        values=tuple(physical_by_parameter[identity] for identity in active_ids),
    )


def _evaluate_batch(batch: _ReactiveBatch, active: ActiveParameterSet) -> Any:
    if isinstance(batch.observation_problem, ReactiveBubbleVLEProblem):
        assert isinstance(batch.descriptor, ReactiveBubbleObservationDescriptor)
        return evaluate_reactive_bubble_observations(
            batch.observation_problem,
            batch.temperature,
            batch.descriptor,
            active_parameters=active,
        )
    assert isinstance(batch.descriptor, HomogeneousReactiveObservationDescriptor)
    assert batch.pressure is not None
    return evaluate_homogeneous_reactive_observations(
        batch.observation_problem,
        batch.temperature,
        batch.pressure,
        batch.descriptor,
        active_parameters=active,
    )


def _preflight(problem: _ResolvedReactiveFit) -> None:
    """Fail before Ceres unless every exact start has complete owner derivatives."""

    branch_signatures: dict[str, tuple[object, ...]] = {}
    for start in problem.controls.starts:
        active = _active_trial(problem, start.values)
        for batch in problem.batches:
            selected_targets = tuple(target for target in batch.targets if target.included)
            if not selected_targets:
                continue
            result = _evaluate_batch(batch, active)
            if result.status != "evaluated":
                failure = result.failure
                if failure is None:
                    raise FitPreflightError(
                        "owner_contract_failure",
                        start.identity,
                        batch.identity,
                        "missing typed failure",
                    )
                # A scientific rejection at an exact caller start is part of
                # multistart accounting. The pure owner will reproduce it for
                # native solve accounting; it must not suppress usable starts.
                continue
            certified_state = getattr(result, "state", None) or getattr(
                result, "liquid_state", None
            )
            if certified_state is None or result.physical_branch is None:
                raise FitPreflightError(
                    "owner_contract_failure",
                    start.identity,
                    batch.identity,
                    "Equilibrium returned no certified branch",
                )
            rows = {row.identity: row for row in result.rows}
            for target in selected_targets:
                row = rows[target.prediction_identity]
                if (
                    row.value is None
                    or row.jacobian is None
                    or len(row.jacobian) != len(active.identities)
                    or not math.isfinite(row.value)
                    or any(not math.isfinite(value) for value in row.jacobian)
                ):
                    raise FitPreflightError(
                        "owner_contract_failure",
                        start.identity,
                        batch.identity,
                        "Equilibrium returned an incomplete selected row",
                    )
            signature = (
                result.physical_branch.branch_family_identity,
                result.physical_branch.continuation_identity,
                certified_state.phase_identity,
                certified_state.phase_role,
                certified_state.chart_topology,
                certified_state.phase_role_receipt.phase_selection_claim,
                certified_state.chemical_search_certificate.selection_label,
            )
            prior = branch_signatures.setdefault(batch.identity, signature)
            if signature != prior:
                raise FitPreflightError(
                    "branch_loss",
                    start.identity,
                    batch.identity,
                    "cross-start physical branch signature changed",
                )


class _EquilibriumOwnerEvaluator:
    __slots__ = ("_active_ids", "_batches", "_parameters")

    def __init__(
        self,
        problem: _ResolvedReactiveFit,
        batches: Sequence[_ReactiveBatch],
    ) -> None:
        self._parameters = problem.parameters
        self._active_ids = tuple(
            identity
            for coordinate in problem.coordinates
            for identity in coordinate.parameter_identities
        )
        self._batches = tuple(batches)

    @staticmethod
    def _rejected(code: NonEvaluableTrialCode, reason: str) -> dict[str, object]:
        return {"status": "rejected", "failure": {"code": code, "reason": reason}}

    def __call__(
        self,
        _start_index: int,
        index: int,
        values: Sequence[float],
        receipt_requested: bool,
    ) -> dict[str, object]:
        batch = self._batches[index]
        try:
            active = ActiveParameterSet(self._parameters, self._active_ids, values=values)
        except (ParameterError, TypeError, ValueError) as error:
            return self._rejected("parameter_topology_mismatch", str(error))
        result = _evaluate_batch(batch, active)
        if result.status != "evaluated":
            failure = result.failure
            if failure is None:
                return self._rejected(
                    "owner_contract_failure",
                    "Equilibrium returned a non-evaluable trial without a typed failure",
                )
            return self._rejected(failure.code, failure.diagnostic)
        if (
            result.descriptor.active_parameter_identities != self._active_ids
            or result.trial is None
            or result.trial.identities != self._active_ids
            or result.trial.values != active.values
            or result.trial.units != tuple(parameter.unit for parameter in active.coordinates)
            or result.trial.active_parameter_fingerprint != active.fingerprint
        ):
            return self._rejected(
                "owner_contract_failure",
                "Equilibrium did not preserve the actual active-parameter trial receipt",
            )
        certified_state = getattr(result, "state", None) or getattr(result, "liquid_state", None)
        equilibrium_state = getattr(result, "equilibrium", None)
        if (
            result.physical_branch is None
            or certified_state is None
            or (
                isinstance(batch.observation_problem, HomogeneousReactiveObservationProblem)
                and equilibrium_state is None
            )
        ):
            return self._rejected(
                "owner_contract_failure",
                "Equilibrium accepted a trial without a physical branch and certified state",
            )
        rows = {row.identity: row for row in result.rows}
        selected = tuple(
            rows[target.prediction_identity] for target in batch.targets if target.included
        )
        if any(row.value is None or row.jacobian is None for row in selected):
            return self._rejected(
                "owner_contract_failure",
                "Equilibrium did not evaluate every selected measurement row",
            )
        response: dict[str, object] = {
            "status": "accepted",
            "predictions": [row.value for row in selected],
            "jacobian": [value for row in selected for value in row.jacobian],
        }
        if receipt_requested:
            if isinstance(batch.observation_problem, ReactiveBubbleVLEProblem):
                response["receipt_json"] = serialize_reactive_bubble_owner_receipt(
                    result, batch_identity=batch.identity
                )
            else:
                observations = tuple(
                    HomogeneousReactiveOwnerObservation(
                        target.identity,
                        target.prediction_identity,
                        target.aggregate_identity,
                        target.covariance_identity,
                    )
                    for target in batch.targets
                    if target.included
                )
                assert batch.pressure is not None
                response["receipt_json"] = serialize_homogeneous_reactive_owner_receipt(
                    result,
                    batch_identity=batch.identity,
                    requested_temperature=batch.temperature,
                    requested_pressure=batch.pressure,
                    observations=observations,
                )
        return response


def _fit_equilibrium_observations(
    problem: _ResolvedReactiveFit,
    *,
    canonical_fit_mapping: Mapping[str, object] | None = None,
) -> FitResult:
    """Fit generic Equilibrium observations with exact owner-supplied Jacobians."""

    if not isinstance(problem, _ResolvedReactiveFit):
        raise TypeError("problem must be a resolved reactive Parameter Fit")
    active_ids = tuple(
        identity
        for coordinate in problem.coordinates
        for identity in coordinate.parameter_identities
    )
    coordinate_indices = {
        identity: index
        for index, coordinate in enumerate(problem.coordinates)
        for identity in coordinate.parameter_identities
    }
    active_batches = tuple(
        batch for batch in problem.batches if any(target.included for target in batch.targets)
    )
    if not active_batches:
        raise ValueError("equilibrium regression has no included target rows")
    if any(
        batch.observation_problem.continuation_reference is None
        or (
            isinstance(batch.observation_problem, ReactiveBubbleVLEProblem)
            and batch.observation_problem.liquid_problem.continuation_reference is None
        )
        for batch in active_batches
    ):
        raise ValueError(
            "equilibrium regression requires an immutable continuation reference per batch"
        )
    _preflight(problem)
    residual_blocks = [
        {
            "multiplier": 1.0,
            "targets": [
                {
                    "observed": target.observed,
                    "residual": target.residual,
                    "scale": target.scale,
                    "multiplier": target.multiplier,
                }
                for target in batch.targets
                if target.included
            ],
        }
        for batch in active_batches
    ]
    native = _core._fit_resolved(
        _EquilibriumOwnerEvaluator(problem, active_batches),
        len(active_ids),
        [coordinate_indices[identity] for identity in active_ids],
        [
            {
                "origin": coordinate.origin,
                "scale": coordinate.scale,
                "lower_affine": coordinate.affine_bounds[0],
                "upper_affine": coordinate.affine_bounds[1],
            }
            for coordinate in problem.coordinates
        ],
        residual_blocks,
        {
            "starts": [
                {
                    "identity": start.identity,
                    "affine": [
                        (value - coordinate.origin) / coordinate.scale
                        for value, coordinate in zip(
                            start.values,
                            problem.coordinates,
                            strict=True,
                        )
                    ],
                }
                for start in problem.controls.starts
            ],
            "affine_bound_margin": problem.controls.affine_bound_margin,
            "workers": 1,
            "maximum_iterations": problem.controls.maximum_iterations,
            "function_tolerance": problem.controls.function_tolerance,
            "gradient_tolerance": problem.controls.gradient_tolerance,
            "parameter_tolerance": problem.controls.parameter_tolerance,
        },
        tuple(sorted(NON_EVALUABLE_TRIAL_CODES)),
    )
    starts = cast(Mapping[str, Any], native).get("starts")
    if not isinstance(starts, Sequence) or isinstance(starts, (str, bytes)):
        raise TypeError("native equilibrium fit starts must be a sequence")
    result = FitResult(
        _specification_receipt(problem, canonical_fit_mapping),
        tuple(FitStart.from_native(item) for item in starts),
        artifact_fingerprint=_python_native_artifact_fingerprint(),
        package_fingerprint=_package_fingerprint(),
    )
    for start in result.starts:
        _validated_owner_receipts(problem, start)
    return result


def _fit_equilibrium_parameter_fit(resolved: ResolvedFit) -> FitResult:
    return _fit_equilibrium_observations(
        _reactive_fit_problem(resolved),
        canonical_fit_mapping=resolved.to_mapping(),
    )


__all__: tuple[str, ...] = ()
