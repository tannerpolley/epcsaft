from __future__ import annotations

import json
import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any, Literal, cast

from .._trials import NonEvaluableTrialCode
from .. import _core
from ..eos import EosError, Mixture
from ..equilibrium import (
    AllComponents,
    GeneralEquilibriumOwnerObservation,
    GeneralReactiveEquilibriumDescriptor,
    GeneralReactiveEquilibriumProblem,
    ProviderModel,
    general_reactive_equilibrium_descriptor,
    general_reactive_equilibrium_problem_from_mapping,
    serialize_general_reactive_equilibrium_receipt,
    solve,
    validate_general_reactive_equilibrium_receipt,
)
from ..equilibrium._api import _continuation_start
from ..equilibrium._grepe import (
    _compile_general_reactive_equilibrium,
    _effective_reaction_system,
    _finite_phase_spec,
    _FinitePhaseCompilation,
    _reactive_finite_source_standard_state,
    _solve_general_reactive_equilibrium,
)
from ..parameters import ActiveParameterSet
from ..units import ParameterError
from ._fit import FitPreflightError
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
        input_row_count=sum(len(batch.targets) for batch in batches),
        dropped_row_count=sum(not target.included for batch in batches for target in batch.targets),
    )


def _active_trial(
    resolved: ResolvedFit,
    values: Sequence[float],
) -> ActiveParameterSet:
    active_ids = resolved.selected_identities
    physical_by_parameter = {
        identity: value
        for coordinate, value in zip(resolved.coordinates, values, strict=True)
        for identity in coordinate.tied_identities
    }
    return ActiveParameterSet(
        resolved.parameters,
        active_ids,
        values=tuple(physical_by_parameter[identity] for identity in active_ids),
    )


def _evaluate_batch(
    batch: _ReactiveBatch,
    active: ActiveParameterSet,
    continuation_start: tuple[tuple[float, ...], float] | None = None,
) -> Any:
    if continuation_start is None:
        return solve(batch.model, batch.problem, active_parameters=active)
    return _solve_general_reactive_equilibrium(
        batch.model,
        batch.problem,
        active_parameters=active,
        continuation_start=continuation_start,
    )


def _preflight(
    resolved: ResolvedFit,
    batches: Sequence[_ReactiveBatch],
) -> dict[tuple[str, str, tuple[str, ...], tuple[float, ...]], Any]:
    """Fail before optimization unless every exact start has complete owner derivatives."""

    branch_signatures: dict[str, tuple[object, ...]] = {}
    exact_start_results: dict[tuple[str, str, tuple[str, ...], tuple[float, ...]], Any] = {}
    for start in resolved.fit.solver.multistart.starts:
        active = _active_trial(resolved, start.values)
        for batch in batches:
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
            exact_start_results[
                (start.identity, batch.identity, active.identities, active.values)
            ] = result
    return exact_start_results


class _EquilibriumOwnerEvaluator:
    __slots__ = (
        "_active_ids",
        "_batches",
        "_exact_start_results",
        "_full_space_starts",
        "_parameters",
        "_resolved",
        "_starts",
    )

    def __init__(
        self,
        resolved: ResolvedFit,
        batches: Sequence[_ReactiveBatch],
    ) -> None:
        self._parameters = resolved.parameters
        self._resolved = resolved
        self._starts = tuple(start.identity for start in resolved.fit.solver.multistart.starts)
        self._active_ids = resolved.selected_identities
        self._batches = tuple(batches)
        self._exact_start_results: (
            dict[tuple[str, str, tuple[str, ...], tuple[float, ...]], Any] | None
        ) = None
        self._full_space_starts: dict[int, tuple[tuple[tuple[float, ...], float], ...]] = {}

    def set_full_space_starts(
        self,
        start_index: int,
        starts: tuple[tuple[tuple[float, ...], float], ...],
    ) -> None:
        self._full_space_starts[start_index] = starts

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
        if not receipt_requested and self._exact_start_results is None:
            self._exact_start_results = _preflight(self._resolved, self._batches)
        batch = self._batches[index]
        try:
            active = ActiveParameterSet(self._parameters, self._active_ids, values=values)
        except (ParameterError, TypeError, ValueError) as error:
            return self._rejected("parameter_topology_mismatch", str(error))
        result = None
        if not receipt_requested:
            assert self._exact_start_results is not None
            result = self._exact_start_results.get(
                (self._starts[_start_index], batch.identity, self._active_ids, active.values)
            )
        if result is None:
            starts = self._full_space_starts.get(_start_index)
            result = _evaluate_batch(
                batch,
                active,
                None if starts is None else starts[index],
            )
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


def _equilibrium_owner_batches(
    resolved: ResolvedFit,
) -> tuple[
    tuple[_ReactiveBatch, ...],
    tuple[_ReactiveBatch, ...],
    _EquilibriumOwnerEvaluator,
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
    return (
        batches,
        active_batches,
        _EquilibriumOwnerEvaluator(resolved, active_batches),
        residual_blocks,
        _specification_receipt(resolved, batches),
    )


def _full_space_batches(
    resolved: ResolvedFit,
    batches: Sequence[_ReactiveBatch],
    owner: _EquilibriumOwnerEvaluator | None = None,
) -> tuple[Any, tuple[Mapping[str, object], ...]] | None:
    active = ActiveParameterSet(resolved.parameters, resolved.selected_identities)
    compiled_batches: list[tuple[Any, ...]] = []
    for batch in batches:
        compiled = _compile_general_reactive_equilibrium(
            batch.model, batch.problem, active_parameters=active
        )
        phase = batch.problem.phases[0]
        outputs = {output.identity: output for output in batch.problem.outputs}
        selected = tuple(
            outputs[target.prediction_identity] for target in batch.targets if target.included
        )
        reaction = _effective_reaction_system(batch.problem)
        if (
            not isinstance(compiled, _FinitePhaseCompilation)
            or len(batch.problem.phases) != 1
            or not isinstance(phase.model, ProviderModel)
            or not isinstance(phase.support, AllComponents)
            or batch.problem.continuation_state is None
            or phase.start is not None
            or any(
                output.selector not in {"phase.amount", "phase.mole_fraction"}
                for output in selected
            )
        ):
            return None
        temperature_k = compiled.temperature_k
        pressure_pa = compiled.pressure_pa
        spec = _finite_phase_spec(batch.problem, temperature_k, pressure_pa)
        reference = batch.problem.continuation_state.phase(phase.identity)
        if (continuation := _continuation_start(reaction, reference)) is not None:
            spec["continuation_start"] = {
                "amounts": continuation[0],
                "volume_m3": continuation[1],
            }
        compiled_batches.append(
            (
                batch,
                spec,
                phase.model.admissible_packing_fraction_interval,
                reaction.strict_interior_amount_floor_mol,
                selected,
            )
        )

    coordinate_count = len(resolved.coordinates)

    def source_state(
        batch: _ReactiveBatch,
        spec: Mapping[str, object],
        trial: ActiveParameterSet,
    ) -> dict[str, object] | None:
        return _reactive_finite_source_standard_state(
            batch.model,
            batch.problem,
            batch.problem.phases[0],
            float(cast(Any, spec["temperature_k"])),
            float(cast(Any, spec["pressure_pa"])),
            trial,
        )

    setups: list[Mapping[str, object]] = []
    start_blocks: list[tuple[Mapping[str, object], ...]] = []
    expected_widths: tuple[tuple[int, int], ...] | None = None
    try:
        for start in resolved.fit.solver.multistart.starts:
            trial = _active_trial(resolved, start.values)
            initial: list[float] = []
            lower: list[float] = []
            upper: list[float] = []
            constraint_lower: list[float] = []
            constraint_upper: list[float] = []
            blocks: list[Mapping[str, object]] = []
            for batch, spec, packing, trace_floor, _ in compiled_batches:
                native = cast(
                    Mapping[str, object],
                    _core._chemical_equilibrium(
                        batch.model._native_model(trial),
                        spec,
                        source_state(batch, spec, trial),
                        packing,
                        trace_floor,
                        True,
                        True,
                        0,
                    ),
                )
                local = tuple(cast(Sequence[float], native["local_variables"]))
                local_lower = tuple(cast(Sequence[float], native["local_lower_bounds"]))
                local_upper = tuple(cast(Sequence[float], native["local_upper_bounds"]))
                if not local or len(local_lower) != len(local) or len(local_upper) != len(local):
                    return None
                blocks.append(
                    {
                        "local_offset": len(initial),
                        "local_width": len(local),
                        "constraint_offset": len(constraint_lower),
                        "constraint_width": len(local),
                        "lower": local_lower,
                        "upper": local_upper,
                    }
                )
                initial.extend(local)
                lower.extend(local_lower)
                upper.extend(local_upper)
                constraint_lower.extend([0.0] * len(local))
                constraint_upper.extend([0.0] * len(local))
            widths = tuple(
                (cast(int, block["local_width"]), cast(int, block["constraint_width"]))
                for block in blocks
            )
            if expected_widths is not None and widths != expected_widths:
                return None
            expected_widths = widths
            rows: list[int] = []
            columns: list[int] = []
            for block in blocks:
                for row in range(cast(int, block["constraint_width"])):
                    rows.extend([cast(int, block["constraint_offset"]) + row] * coordinate_count)
                    columns.extend(range(coordinate_count))
                    rows.extend(
                        [cast(int, block["constraint_offset"]) + row]
                        * cast(int, block["local_width"])
                    )
                    columns.extend(
                        coordinate_count + cast(int, block["local_offset"]) + column
                        for column in range(cast(int, block["local_width"]))
                    )
            setups.append(
                {
                    "initial": initial,
                    "lower": lower,
                    "upper": upper,
                    "constraint_lower": constraint_lower,
                    "constraint_upper": constraint_upper,
                    "jacobian_rows": rows,
                    "jacobian_columns": columns,
                }
            )
            start_blocks.append(tuple(blocks))
    except (EosError, ParameterError, RuntimeError, TypeError, ValueError):
        return None

    coordinate_by_parameter = {
        identity: index
        for index, coordinate in enumerate(resolved.coordinates)
        for identity in coordinate.tied_identities
    }

    def evaluate(start_index: int, values: Sequence[float]) -> Mapping[str, object]:
        physical = tuple(
            coordinate.origin + float(values[index]) * coordinate.scale
            for index, coordinate in enumerate(resolved.coordinates)
        )
        trial = ActiveParameterSet(
            resolved.parameters,
            resolved.selected_identities,
            values=tuple(
                physical[coordinate_by_parameter[identity]]
                for identity in resolved.selected_identities
            ),
        )
        variable_count = coordinate_count + len(
            cast(Sequence[object], setups[start_index]["initial"])
        )
        residuals: list[float] = []
        residual_jacobian: list[float] = []
        constraints: list[float] = []
        constraint_jacobian: list[float] = []
        continuation_starts: list[tuple[tuple[float, ...], float]] = []
        for compiled, bounds in zip(compiled_batches, start_blocks[start_index], strict=True):
            batch, spec, packing, _, outputs = compiled
            local_width = cast(int, bounds["local_width"])
            local_offset = coordinate_count + cast(int, bounds["local_offset"])
            native = cast(
                Mapping[str, object],
                _core._evaluate_provider_reaction_local(
                    batch.model._native_model(trial),
                    spec,
                    source_state(batch, spec, trial),
                    packing,
                    values[local_offset : local_offset + local_width],
                    cast(Sequence[float], bounds["lower"]),
                    cast(Sequence[float], bounds["upper"]),
                    True,
                ),
            )
            amounts = tuple(cast(Sequence[float], native["amounts"]))
            total = math.fsum(amounts)
            volume = float(cast(float, native["volume_m3"]))
            continuation_starts.append((amounts, volume))
            chart = cast(Sequence[float], native["physical_to_chart_jacobian"])
            for target, output in zip(
                (target for target in batch.targets if target.included), outputs, strict=True
            ):
                numerator = math.fsum(
                    coefficient * amount
                    for coefficient, amount in zip(output.coefficients, amounts, strict=True)
                )
                if output.selector == "phase.amount":
                    prediction = numerator
                    output_gradient = (*output.coefficients, 0.0)
                else:
                    prediction = numerator / total
                    output_gradient = (
                        *(
                            (coefficient * total - numerator) / total**2
                            for coefficient in output.coefficients
                        ),
                        0.0,
                    )
                local_jacobian = tuple(
                    math.fsum(
                        output_gradient[row] * chart[row * local_width + column]
                        for row in range(len(output_gradient))
                    )
                    for column in range(local_width)
                )
                multiplier = target.multiplier
                factor = multiplier / target.scale
                if target.residual == "log_ratio":
                    factor = multiplier / prediction
                    residual = multiplier * math.log(prediction / target.observed)
                elif (
                    target.residual == "left_censored"
                    and prediction <= target.observed
                    or target.residual == "right_censored"
                    and prediction >= target.observed
                ):
                    factor = residual = 0.0
                else:
                    residual = multiplier * (prediction - target.observed) / target.scale
                residuals.append(residual)
                row_jacobian = [0.0] * variable_count
                row_jacobian[local_offset : local_offset + local_width] = (
                    factor * derivative for derivative in local_jacobian
                )
                residual_jacobian.extend(row_jacobian)
            block_constraints = tuple(cast(Sequence[float], native["kkt_root_values"]))
            constraints.extend(block_constraints)
            state_jacobian = cast(Sequence[float], native["kkt_root_jacobian"])
            order = tuple(cast(Sequence[str], native["local_parameter_order"]))
            parameter_jacobian = cast(Sequence[float], native["local_parameter_jacobian"])
            parameter_columns = {
                identity: order.index(f"provider_parameter[{identity}]")
                for identity in resolved.selected_identities
            }
            for row in range(len(block_constraints)):
                for coordinate in resolved.coordinates:
                    constraint_jacobian.append(
                        coordinate.scale
                        * math.fsum(
                            parameter_jacobian[row * len(order) + parameter_columns[identity]]
                            for identity in coordinate.tied_identities
                            if identity in parameter_columns
                        )
                    )
                constraint_jacobian.extend(
                    state_jacobian[row * local_width : (row + 1) * local_width]
                )
        if owner is not None:
            owner.set_full_space_starts(start_index, tuple(continuation_starts))
        return {
            "residuals": residuals,
            "residual_jacobian": residual_jacobian,
            "constraints": constraints,
            "constraint_jacobian": constraint_jacobian,
        }

    return evaluate, tuple(setups)


__all__: tuple[str, ...] = ()
