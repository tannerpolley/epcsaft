from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any, Literal, cast

from .._canonical import _canonical_json
from ..eos import Mixture
from ..equilibrium import (
    GeneralEquilibriumOwnerObservation,
    GeneralReactiveEquilibriumDescriptor,
    GeneralReactiveEquilibriumProblem,
    general_reactive_equilibrium_descriptor,
    general_reactive_equilibrium_problem_from_mapping,
    serialize_general_reactive_equilibrium_receipt,
    validate_general_reactive_equilibrium_receipt,
)
from ..equilibrium._grepe import (
    Fixed,
    Solved,
    _fingerprint,
    _map_native_general_reactive_equilibrium,
    _native_general_reactive_equilibrium_inputs,
)
from ..parameters import ActiveParameterSet
from ._model import EquilibriumFitObservation
from ._resolve import ResolvedFit
from ._result import (
    FitCoordinateReceipt,
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


def _native_intensive(value: Fixed | Solved, unit: str) -> Mapping[str, object]:
    if isinstance(value, Fixed):
        return {
            "solved": False,
            "value": float(value.value.to(unit).magnitude),
            "bounds": (),
            "starts": (),
        }
    return {
        "solved": True,
        "value": float(value.initial.to(unit).magnitude),
        "bounds": tuple(float(bound.to(unit).magnitude) for bound in value.bounds),
        "starts": tuple(float(start.to(unit).magnitude) for start in value.starts),
    }


def _is_sha256(value: str) -> bool:
    if not value.startswith("sha256:") or len(value) != 71:
        return False
    return all(character in "0123456789abcdef" for character in value[7:])


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


def _reactive_batches(resolved: ResolvedFit) -> tuple[_ReactiveBatch, ...]:
    observations = tuple(
        item for item in resolved.fit.observations if isinstance(item, EquilibriumFitObservation)
    )
    active = ActiveParameterSet(resolved.parameters, resolved.selected_identities)
    model = Mixture(resolved.parameters)
    batches: list[_ReactiveBatch] = []
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
    result = tuple(batches)
    if len({batch.identity for batch in result}) != len(result):
        raise ValueError("equilibrium regression batch identities must be unique")
    targets = tuple(target.identity for batch in result for target in batch.targets)
    if len(set(targets)) != len(targets):
        raise ValueError("equilibrium regression target identities must be globally unique")
    if any(
        batch.descriptor.active_parameter_identities != resolved.selected_identities
        for batch in result
    ):
        raise ValueError("Equilibrium descriptor does not preserve active parameter order")
    return result


def _validated_owner_receipts(
    resolved: ResolvedFit,
    batches: Sequence[_ReactiveBatch],
    start: FitStart,
) -> tuple[Mapping[str, Any], ...]:
    if not start.usable:
        if start.owner_receipts:
            raise ValueError("unusable equilibrium start must not retain owner receipts")
        return ()
    active_batches = tuple(
        batch for batch in batches if any(target.included for target in batch.targets)
    )
    if len(start.owner_receipts) != len(active_batches):
        raise ValueError("usable equilibrium start has incomplete owner batch receipts")
    parameter_ids = resolved.selected_identities
    coordinate_values = {
        coordinate.identity: value
        for coordinate, value in zip(resolved.coordinates, start.final_physical, strict=True)
    }
    coordinate_by_parameter = {
        parameter: coordinate.identity
        for coordinate in resolved.coordinates
        for parameter in coordinate.tied_identities
    }
    active_values = tuple(
        coordinate_values[coordinate_by_parameter[identity]] for identity in parameter_ids
    )
    active = ActiveParameterSet(resolved.parameters, parameter_ids, values=active_values)
    validated: list[Mapping[str, Any]] = []
    for batch, opaque_receipt in zip(active_batches, start.owner_receipts, strict=True):
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


def _specification_receipt(
    resolved: ResolvedFit,
    batches: Sequence[_ReactiveBatch],
) -> FitSpecificationReceipt:
    active_batches = tuple(
        batch for batch in batches if any(target.included for target in batch.targets)
    )
    return FitSpecificationReceipt.create(
        resolved.to_mapping(),
        parameters=resolved.parameters,
        coordinates=tuple(
            FitCoordinateReceipt(
                coordinate.identity,
                coordinate.tied_identities,
                coordinate.unit,
                coordinate.origin,
                coordinate.scale,
            )
            for coordinate in resolved.coordinates
        ),
        rows=tuple(
            FitRowReceipt(
                target.identity,
                target.source_identity,
                target.source_hash,
                target.prediction_identity,
                _canonical_json(
                    {
                        "aggregate_identity": target.aggregate_identity,
                        "batch_identity": batch.identity,
                        "covariance_identity": target.covariance_identity,
                        "measurement_classification": target.measurement_classification,
                        "role": target.role,
                    },
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
        input_row_count=sum(len(batch.targets) for batch in batches),
        dropped_row_count=sum(not target.included for batch in batches for target in batch.targets),
    )


def _equilibrium_owner_batches(
    resolved: ResolvedFit,
) -> tuple[
    tuple[_ReactiveBatch, ...],
    tuple[_ReactiveBatch, ...],
    tuple[Mapping[str, object], ...],
    tuple[Mapping[str, object], ...],
    FitSpecificationReceipt,
]:
    batches = _reactive_batches(resolved)
    active_batches = tuple(
        batch for batch in batches if any(target.included for target in batch.targets)
    )
    if not active_batches:
        raise ValueError("equilibrium regression has no included target rows")
    if any(
        not batch.problem.intensive_boundaries
        and batch.descriptor.formulation != "temperature_boundary"
        and (
            (
                not any(phase.amount_role == "incipient" for phase in batch.problem.phases)
                and any(
                    phase.amount_role == "finite"
                    and batch.problem.continuation_state is None
                    and phase.start is None
                    for phase in batch.problem.phases
                )
            )
            or (
                any(phase.amount_role == "incipient" for phase in batch.problem.phases)
                and batch.problem.continuation_state is None
            )
        )
        for batch in active_batches
    ):
        raise ValueError("equilibrium regression requires a continuation state or start per batch")
    residual_blocks = tuple(
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
    )
    active = ActiveParameterSet(resolved.parameters, resolved.selected_identities)
    experiments = []
    for batch in active_batches:
        incipient_model, native_problem, reference_fingerprint = (
            _native_general_reactive_equilibrium_inputs(
                batch.model, batch.problem, active_parameters=active
            )
        )
        output_index = {
            output.identity: index for index, output in enumerate(batch.problem.outputs)
        }
        experiments.append(
            {
                "kind": "grepe",
                "structure_fingerprint": _fingerprint(
                    "epcsaft.regression-grepe-structure.v1",
                    {
                        "parameter_fingerprint": batch.descriptor.parameter_fingerprint,
                        "active_parameter_selection_fingerprint": (
                            batch.descriptor.active_parameter_selection_fingerprint
                        ),
                        "active_parameter_identities": (
                            batch.descriptor.active_parameter_identities
                        ),
                        "reaction_fingerprint": batch.descriptor.reaction_fingerprint,
                        "topology_fingerprint": batch.descriptor.topology_fingerprint,
                        "reference_fingerprint": batch.descriptor.reference_fingerprint,
                        "domain_fingerprint": batch.descriptor.domain_fingerprint,
                        "source_reference_fingerprint": reference_fingerprint,
                        "temperature": _native_intensive(batch.problem.temperature, "kelvin"),
                        "pressure": _native_intensive(batch.problem.pressure, "pascal"),
                        "outputs": tuple(
                            (
                                output.identity,
                                output.selector,
                                output.unit,
                                output.basis,
                                output.phase_identity,
                                output.coefficients,
                                output.solvent_mass_coefficients_kg_per_mol,
                                output.support,
                                output.censor_limit,
                                output.aggregate_identity,
                                output.covariance_identity,
                            )
                            for output in batch.problem.outputs
                        ),
                    },
                ),
                "model": batch.model._native_model(active),
                "incipient_model": incipient_model,
                "problem": native_problem,
                "temperature": _native_intensive(batch.problem.temperature, "kelvin"),
                "pressure": _native_intensive(batch.problem.pressure, "pascal"),
                "source_reference_fingerprint": reference_fingerprint,
                "selected_rows": [
                    output_index[target.prediction_identity]
                    for target in batch.targets
                    if target.included
                ],
            }
        )
    return (
        batches,
        active_batches,
        tuple(experiments),
        residual_blocks,
        _specification_receipt(resolved, batches),
    )


def _serialize_native_owner_result(
    batch: _ReactiveBatch,
    active: ActiveParameterSet,
    native: object,
) -> str:
    result = _map_native_general_reactive_equilibrium(batch.model, batch.problem, active, native)
    return serialize_general_reactive_equilibrium_receipt(
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


__all__: tuple[str, ...] = ()
