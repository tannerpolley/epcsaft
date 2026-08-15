from __future__ import annotations

import json
import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any, Literal, cast

from .. import _core
from .._trials import NON_EVALUABLE_TRIAL_CODES, NonEvaluableTrialCode
from ..eos import Mixture, _python_native_artifact_fingerprint
from ..equilibrium import (
    GeneralEquilibriumOwnerObservation,
    GeneralReactiveEquilibriumDescriptor,
    GeneralReactiveEquilibriumProblem,
    general_reactive_equilibrium_descriptor,
    general_reactive_equilibrium_problem_from_mapping,
    phase_equilibrium,
    serialize_general_reactive_equilibrium_receipt,
    validate_general_reactive_equilibrium_receipt,
)
from ..parameters import ActiveParameterSet, Parameters
from ..units import ParameterError
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
    FitRowReceipt,
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
    basis: str
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
            self.basis,
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
    model: Mixture
    problem: GeneralReactiveEquilibriumProblem
    descriptor: GeneralReactiveEquilibriumDescriptor
    targets: tuple[_ReactiveTarget, ...]

    def __post_init__(self) -> None:
        object.__setattr__(self, "targets", tuple(self.targets))
        if not self.identity or not self.targets:
            raise ValueError("reactive regression batch identity and targets are required")
        if len({target.identity for target in self.targets}) != len(self.targets):
            raise ValueError("reactive regression target identities must be unique within a batch")
        if self.descriptor.problem_identity != self.problem.identity:
            raise ValueError("equilibrium regression descriptor and problem identities differ")
        rows = {row.identity: row for row in self.descriptor.outputs}
        for target in self.targets:
            row = rows.get(target.prediction_identity)
            if row is None:
                raise ValueError("regression target is absent from its Equilibrium descriptor")
            if target.unit != row.unit:
                raise ValueError("regression target unit does not match its Equilibrium row")
            if target.basis != row.basis:
                raise ValueError("regression target basis does not match its Equilibrium row")
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
            if target.residual == "log_ratio" and row.support != "positive":
                raise ValueError("log-ratio residual requires a positive Equilibrium row")


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
    if (
        not isinstance(value, Mapping)
        or any(type(key) is not str for key in value)
        or set(value) != fields
    ):
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


def _text(value: object, label: str) -> str:
    if not isinstance(value, str) or not value or value != value.strip():
        raise TypeError(f"{label} must be a nonempty string")
    return value


def _boolean(value: object, label: str) -> bool:
    if not isinstance(value, bool):
        raise TypeError(f"{label} must be a boolean")
    return value


def _optional_text(value: object, label: str) -> str | None:
    return None if value is None else _text(value, label)


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
        "basis",
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
            _text(record["identity"], "target identity"),
            _text(record["prediction_identity"], "prediction identity"),
            _number(record["observed"], "observed value"),
            _text(record["unit"], "target unit"),
            cast(Any, _text(record["residual"], "residual transform")),
            _number(record["scale"], "residual scale"),
            _number(record["multiplier"], "target multiplier"),
            _text(record["role"], "target role"),
            _text(record["measurement_classification"], "measurement classification"),
            _text(record["basis"], "target basis"),
            _text(record["source_identity"], "source identity"),
            _text(record["source_hash"], "source hash"),
            _boolean(record["included"], "included flag"),
            "" if record["drop_reason"] is None else _text(record["drop_reason"], "drop reason"),
            _optional_text(record["aggregate_identity"], "aggregate identity"),
            _optional_text(record["covariance_identity"], "covariance identity"),
        )
        for record in (
            _exact_mapping(raw, _TARGET_FIELDS, "reactive regression target")
            for raw in _sequence(value, "reactive targets")
        )
    )


def _reactive_fit_problem(resolved: ResolvedFit) -> _ResolvedReactiveFit:
    observations = tuple(
        item for item in resolved.fit.observations if isinstance(item, EquilibriumFitObservation)
    )
    if len(observations) != len(resolved.fit.observations):
        raise TypeError("reactive Parameter Fit cannot mix owner families")
    active = ActiveParameterSet(resolved.parameters, resolved.selected_identities)
    model = Mixture(resolved.parameters)
    batches = []
    for specification in observations:
        item = specification.to_mapping()
        request = general_reactive_equilibrium_problem_from_mapping(item["request"])
        batches.append(
            _ReactiveBatch(
                _text(item["identity"], "equilibrium batch identity"),
                model,
                request,
                general_reactive_equilibrium_descriptor(model, request, active_parameters=active),
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
        observations = tuple(
            GeneralEquilibriumOwnerObservation(
                target.identity,
                target.prediction_identity,
                target.aggregate_identity,
                target.covariance_identity,
            )
            for target in batch.targets
            if target.included
        )
        receipt_digest = validate_general_reactive_equilibrium_receipt(
            receipt_json,
            descriptor=batch.descriptor,
            problem=batch.problem,
            batch_identity=batch.identity,
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
        "basis": target.basis,
        "source_identity": target.source_identity,
        "source_hash": target.source_hash,
        "included": target.included,
        "drop_reason": target.drop_reason,
        "aggregate_identity": target.aggregate_identity,
        "covariance_identity": target.covariance_identity,
    }


def _batch_payload(batch: _ReactiveBatch) -> dict[str, object]:
    return {
        "identity": batch.identity,
        "family": batch.descriptor.identity,
        "descriptor_fingerprint": batch.descriptor.fingerprint,
        "problem_fingerprint": batch.descriptor.problem_fingerprint,
        "reaction_fingerprint": batch.descriptor.reaction_fingerprint,
        "topology_fingerprint": batch.descriptor.topology_fingerprint,
        "reference_fingerprint": batch.descriptor.reference_fingerprint,
        "domain_fingerprint": batch.descriptor.domain_fingerprint,
        "problem_identity": batch.problem.identity,
        "phase_identities": batch.descriptor.phase_identities,
        "certificate_scope": batch.descriptor.certificate_scope,
        "outputs": tuple(
            {
                "identity": output.identity,
                "selector": output.selector,
                "unit": output.unit,
                "basis": output.basis,
            }
            for output in batch.problem.outputs
        ),
        "targets": tuple(_target_payload(target) for target in batch.targets),
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
        rows=tuple(
            FitRowReceipt(
                target.identity,
                target.source_identity,
                target.source_hash,
                target.prediction_identity,
                json.dumps(
                    {
                        "aggregate_identity": target.aggregate_identity,
                        "batch_identity": batch.identity,
                        "covariance_identity": target.covariance_identity,
                        "measurement_classification": target.measurement_classification,
                        "role": target.role,
                    },
                    separators=(",", ":"),
                    sort_keys=True,
                ),
                target.observed,
                target.unit,
                target.basis,
                target.residual,
                target.scale,
                target.multiplier,
            )
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
    return phase_equilibrium(batch.model, batch.problem, active_parameters=active)


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
            if result.branch_signature is None:
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
            signature = result.branch_signature
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
            result.descriptor is None
            or result.descriptor.active_parameter_identities != self._active_ids
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
        if result.branch_signature is None:
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
            response["receipt_json"] = serialize_general_reactive_equilibrium_receipt(
                result,
                batch_identity=batch.identity,
                observations=tuple(
                    GeneralEquilibriumOwnerObservation(
                        target.identity,
                        target.prediction_identity,
                        target.aggregate_identity,
                        target.covariance_identity,
                    )
                    for target in batch.targets
                    if target.included
                ),
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
        any(
            phase.amount_role == "finite"
            and phase.continuation_reference is None
            and phase.start is None
            for phase in batch.problem.phases
        )
        or (
            any(phase.amount_role == "incipient" for phase in batch.problem.phases)
            and batch.problem.continuation_reference is None
        )
        for batch in active_batches
    ):
        raise ValueError(
            "equilibrium regression requires an immutable continuation reference or start per batch"
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
