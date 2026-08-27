from __future__ import annotations

import hashlib
import json
import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from os import PathLike
from pathlib import Path
from typing import Any, cast

from .._canonical import _canonical_json, _canonical_json_bytes
from .._trials import NON_EVALUABLE_TRIAL_CODES, NonEvaluableTrialCode
from ..parameters import Parameters
from ..units import unit_registry
from ._resolve import ResolvedFit


def _finite_or_none(value: float) -> float | None:
    return value if math.isfinite(value) else None


def _finite_values(values: Sequence[float]) -> list[float | None]:
    return [_finite_or_none(value) for value in values]


@dataclass(frozen=True, slots=True)
class _ReceiptRow:
    values: tuple[tuple[str, float], ...]

    def value(self, name: str) -> float:
        for identity, value in self.values:
            if identity == name:
                return value
        raise KeyError(name)


@dataclass(frozen=True, slots=True)
class OwnerReceipt:
    canonical_bytes_utf8: str
    digest: str

    @classmethod
    def from_canonical_bytes(cls, value: str) -> OwnerReceipt:
        if not value:
            raise ValueError("native owner receipts must be nonempty when present")
        return cls(value, f"sha256:{hashlib.sha256(value.encode('utf-8')).hexdigest()}")

    def to_mapping(self) -> dict[str, str]:
        return {
            "canonical_bytes_utf8": self.canonical_bytes_utf8,
            "digest": self.digest,
        }


@dataclass(frozen=True, slots=True)
class FitCoordinateReceipt:
    identity: str
    parameter_identities: tuple[str, ...]
    unit: str
    origin: float
    scale: float


@dataclass(frozen=True, slots=True)
class FitRowReceipt:
    """Immutable measurement metadata aligned with one optimizer row."""

    identity: str
    source_identity: str
    source_hash: str
    selector: str
    series_metadata_json: str
    observed: float
    unit: str
    basis: str
    residual: str
    scale: float
    multiplier: float

    def __post_init__(self) -> None:
        if (
            not self.identity
            or not self.source_identity
            or not self.source_hash.startswith("sha256:")
            or not self.selector
            or not self.series_metadata_json
            or not self.unit
            or not self.basis
            or not self.residual
            or not math.isfinite(self.observed)
            or not math.isfinite(self.scale)
            or self.scale <= 0.0
            or not math.isfinite(self.multiplier)
            or self.multiplier <= 0.0
        ):
            raise ValueError("fit row receipt is incomplete")
        metadata = json.loads(self.series_metadata_json)
        if not isinstance(metadata, Mapping):
            raise TypeError("fit row series metadata must be a canonical object")
        if _canonical_json(metadata) != self.series_metadata_json:
            raise ValueError("fit row series metadata is not canonical JSON")

    def to_mapping(self) -> dict[str, object]:
        return {
            "identity": self.identity,
            "source_identity": self.source_identity,
            "source_hash": self.source_hash,
            "selector": self.selector,
            "series_metadata": json.loads(self.series_metadata_json),
            "observed": self.observed,
            "unit": self.unit,
            "basis": self.basis,
            "residual": self.residual,
            "scale": self.scale,
            "multiplier": self.multiplier,
        }


@dataclass(frozen=True, slots=True)
class FitSpecificationReceipt:
    canonical_bytes_utf8: str
    digest: str
    parameters: Parameters
    coordinates: tuple[FitCoordinateReceipt, ...]
    rows: tuple[FitRowReceipt, ...]
    row_identities: tuple[str, ...]
    block_sizes: tuple[int, ...]
    input_row_count: int
    dropped_row_count: int

    @classmethod
    def create(
        cls,
        mapping: Mapping[str, object],
        *,
        parameters: Parameters,
        coordinates: Sequence[FitCoordinateReceipt],
        rows: Sequence[FitRowReceipt],
        block_sizes: Sequence[int],
        input_row_count: int | None = None,
        dropped_row_count: int = 0,
    ) -> FitSpecificationReceipt:
        frozen_rows = tuple(rows)
        row_identities = tuple(row.identity for row in frozen_rows)
        frozen_blocks = tuple(block_sizes)
        frozen_coordinates = tuple(coordinates)
        declared_input_rows = len(row_identities) if input_row_count is None else input_row_count
        if (
            sum(frozen_blocks) != len(row_identities)
            or declared_input_rows != len(row_identities) + dropped_row_count
            or dropped_row_count < 0
            or len(set(row_identities)) != len(row_identities)
        ):
            raise ValueError("fit specification row accounting is inconsistent")
        receipt_mapping = {
            "fit": mapping,
            "input_parameters": parameters.to_mapping(),
            "coordinates": [
                {
                    "identity": coordinate.identity,
                    "parameter_identities": list(coordinate.parameter_identities),
                    "unit": coordinate.unit,
                    "origin": coordinate.origin,
                    "scale": coordinate.scale,
                }
                for coordinate in frozen_coordinates
            ],
            "rows": [row.to_mapping() for row in frozen_rows],
            "row_identities": list(row_identities),
            "block_sizes": list(frozen_blocks),
            "input_row_count": declared_input_rows,
            "dropped_row_count": dropped_row_count,
        }
        canonical = _canonical_json(receipt_mapping)
        return cls(
            canonical,
            f"sha256:{hashlib.sha256(canonical.encode('utf-8')).hexdigest()}",
            parameters,
            frozen_coordinates,
            frozen_rows,
            row_identities,
            frozen_blocks,
            declared_input_rows,
            dropped_row_count,
        )

    @classmethod
    def from_resolved_fit(cls, resolved: ResolvedFit) -> FitSpecificationReceipt:
        rows: list[FitRowReceipt] = []
        block_sizes: list[int] = []
        for observation_index, observation in enumerate(resolved.observations):
            specification = observation.specification
            canonical_units = dict(observation.canonical_units)
            bases = dict(specification.basis)
            residuals = dict(specification.residuals)
            source = specification.source
            source_identity = source.citation or source.path
            if observation.rows is not None:
                resolved_rows = tuple(
                    (_ReceiptRow(row.values), specification.multiplier) for row in observation.rows
                )
            elif observation.quadrature is not None:
                rule = observation.quadrature.rules[-1]
                target = specification.targets[0]
                resolved_rows = tuple(
                    (
                        _ReceiptRow((("temperature_k", temperature), (target, observed))),
                        specification.multiplier * math.sqrt(weight),
                    )
                    for temperature, weight, observed in zip(
                        rule.nodes, rule.weights, rule.observed_values, strict=True
                    )
                )
            else:
                raise ValueError("resolved observation has no measurement rows")
            state_fields = tuple(name for name, _ in source.domain.bounds)
            for state_index, (row, multiplier) in enumerate(resolved_rows):
                block_sizes.append(len(observation.specification.targets))
                metadata = {
                    "observation_index": observation_index,
                    "state_index": state_index,
                    "state": {name: row.value(name) for name in state_fields},
                    "composition": specification.composition,
                }
                for target in specification.targets:
                    rows.append(
                        FitRowReceipt(
                            f"observation/{observation_index}/state/{state_index}/{target}",
                            source_identity,
                            observation.source_check.actual_sha256,
                            target,
                            _canonical_json(metadata),
                            row.value(target),
                            canonical_units[target],
                            bases[target],
                            residuals[target],
                            row.value(target),
                            multiplier,
                        )
                    )
        return cls.create(
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
            rows=rows,
            block_sizes=block_sizes,
            input_row_count=len(rows),
        )

    def to_mapping(self) -> Mapping[str, object]:
        value = json.loads(self.canonical_bytes_utf8)
        if not isinstance(value, Mapping):
            raise TypeError("canonical fit specification receipt is not an object")
        return cast(Mapping[str, object], value)


@dataclass(frozen=True, slots=True)
class BatchOutcome:
    index: int
    row_count: int
    status: str
    failure_code: NonEvaluableTrialCode | None

    @classmethod
    def from_native(cls, value: Mapping[str, Any]) -> BatchOutcome:
        status = str(value["status"])
        if status not in {"evaluated", "rejected", "not_evaluated"}:
            raise ValueError("native batch outcome status is outside the public vocabulary")
        raw_code = str(value["failure_code"])
        if raw_code and raw_code not in NON_EVALUABLE_TRIAL_CODES:
            raise ValueError("native batch failure code is outside the public vocabulary")
        return cls(
            int(value["index"]),
            int(value["row_count"]),
            status,
            None if not raw_code else cast(NonEvaluableTrialCode, raw_code),
        )

    def to_mapping(self) -> dict[str, object]:
        return {
            "index": self.index,
            "row_count": self.row_count,
            "status": self.status,
            "failure_code": self.failure_code,
        }


@dataclass(frozen=True, slots=True)
class FitAccounting:
    input_rows: int
    included_rows: int
    dropped_rows: int
    attempted_starts: int
    accepted_starts: int
    failed_starts: int
    evaluated_row_trials: int
    rejected_row_trials: int
    failed_row_trials: int
    non_evaluable_trials: int

    def to_mapping(self) -> dict[str, int]:
        return {field: int(getattr(self, field)) for field in self.__dataclass_fields__}


@dataclass(frozen=True, slots=True)
class OwnerFailure:
    code: NonEvaluableTrialCode
    reason: str
    count: int

    @classmethod
    def from_native(cls, value: Mapping[str, Any]) -> OwnerFailure:
        code = str(value["code"])
        if code not in NON_EVALUABLE_TRIAL_CODES:
            raise ValueError("native owner failure code is outside the public vocabulary")
        return cls(
            cast(NonEvaluableTrialCode, code),
            str(value["reason"]),
            int(value["count"]),
        )

    def to_mapping(self) -> dict[str, object]:
        return {"code": self.code, "reason": self.reason, "count": self.count}


@dataclass(frozen=True, slots=True)
class FitStart:
    index: int
    start_identity: str
    initial_affine: tuple[float, ...]
    initial_physical: tuple[float, ...]
    final_physical: tuple[float, ...]
    final_affine: tuple[float, ...]
    predictions: tuple[float, ...]
    residuals: tuple[float, ...]
    physical_jacobian: tuple[float, ...]
    affine_jacobian: tuple[float, ...]
    singular_values: tuple[float, ...]
    rank: int
    free_subspace_rank: int
    condition_number: float
    affine_covariance: tuple[float, ...]
    degrees_of_freedom: int
    information_row_count: int
    free_coordinate_indices: tuple[int, ...]
    residual_variance: float
    covariance_status: str
    active_bounds: tuple[int, ...]
    solver_status: str
    failure_reason: str
    trial_failure_code: NonEvaluableTrialCode | None
    iterations: int
    non_evaluable_count: int
    owner_failures: tuple[OwnerFailure, ...]
    owner_receipts: tuple[OwnerReceipt, ...]
    batch_outcomes: tuple[BatchOutcome, ...]
    cost: float
    usable: bool
    optimizer_basin: int | None

    @property
    def outcomes(self) -> dict[str, dict[str, str]]:
        solver_converged = self.solver_status in {
            "solve_succeeded",
            "solved_to_acceptable_level",
        }
        if self.usable and all(outcome.status == "evaluated" for outcome in self.batch_outcomes):
            physical_status = "certified"
        elif any(outcome.status == "rejected" for outcome in self.batch_outcomes):
            physical_status = "rejected"
        else:
            physical_status = "not_evaluated"
        return {
            "solver_convergence": {
                "status": "converged" if solver_converged else "not_converged",
                "reason": self.solver_status,
            },
            "numerical_convergence": {
                "status": "not_assessed" if solver_converged else "not_converged",
                "reason": (
                    "no independent tolerance or initialization sensitivity study"
                    if solver_converged
                    else self.solver_status
                ),
            },
            "physical_certification": {
                "status": physical_status,
                "reason": (
                    "all owner batches evaluated"
                    if physical_status == "certified"
                    else self.failure_reason or physical_status
                ),
            },
            "scientific_fit_acceptance": {
                "status": "not_declared",
                "reason": "the fit specification declares no scientific acceptance criterion",
            },
        }

    @classmethod
    def from_native(cls, value: Mapping[str, Any]) -> FitStart:
        diagnostics = value["diagnostics"]
        if not isinstance(diagnostics, Mapping):
            raise TypeError("native fit diagnostics must be a mapping")
        raw_trial_failure_code = str(value["trial_failure_code"])
        if raw_trial_failure_code and raw_trial_failure_code not in NON_EVALUABLE_TRIAL_CODES:
            raise ValueError("native trial failure code is outside the public vocabulary")
        owner_receipts = tuple(
            OwnerReceipt.from_canonical_bytes(str(item)) for item in value["owner_receipts_json"]
        )
        return cls(
            index=int(value["index"]),
            start_identity=str(value["start_identity"]),
            initial_affine=tuple(float(item) for item in value["initial_affine"]),
            initial_physical=tuple(float(item) for item in value["initial_physical"]),
            final_physical=tuple(float(item) for item in value["final_physical"]),
            final_affine=tuple(float(item) for item in value["final_affine"]),
            predictions=tuple(float(item) for item in value["predictions"]),
            residuals=tuple(float(item) for item in value["residuals"]),
            physical_jacobian=tuple(float(item) for item in value["physical_jacobian"]),
            affine_jacobian=tuple(float(item) for item in value["affine_jacobian"]),
            singular_values=tuple(float(item) for item in diagnostics["singular_values"]),
            rank=int(diagnostics["rank"]),
            free_subspace_rank=int(diagnostics["free_subspace_rank"]),
            condition_number=float(diagnostics["condition_number"]),
            affine_covariance=tuple(float(item) for item in diagnostics["affine_covariance"]),
            degrees_of_freedom=int(diagnostics["degrees_of_freedom"]),
            information_row_count=int(diagnostics["information_row_count"]),
            free_coordinate_indices=tuple(
                int(item) for item in diagnostics["free_coordinate_indices"]
            ),
            residual_variance=float(diagnostics["residual_variance"]),
            covariance_status=str(diagnostics["covariance_status"]),
            active_bounds=tuple(int(item) for item in value["active_bounds"]),
            solver_status=str(value["solver_status"]),
            failure_reason=str(value["failure_reason"]),
            trial_failure_code=(
                None
                if not raw_trial_failure_code
                else cast(NonEvaluableTrialCode, raw_trial_failure_code)
            ),
            iterations=int(value["iterations"]),
            non_evaluable_count=int(value["non_evaluable_count"]),
            owner_failures=tuple(
                OwnerFailure.from_native(item) for item in value["owner_failures"]
            ),
            owner_receipts=owner_receipts,
            batch_outcomes=tuple(
                BatchOutcome.from_native(item) for item in value["batch_outcomes"]
            ),
            cost=float(value["cost"]),
            usable=bool(value["usable"]),
            optimizer_basin=(
                None if value["optimizer_basin"] is None else int(value["optimizer_basin"])
            ),
        )

    def to_mapping(self, block_sizes: Sequence[int]) -> dict[str, object]:
        contributions: list[float | None] = []
        offset = 0
        for size in block_sizes:
            block = self.residuals[offset : offset + size]
            contributions.append(
                0.5 * sum(value * value for value in block)
                if all(math.isfinite(value) for value in block)
                else None
            )
            offset += size
        return {
            "index": self.index,
            "start_identity": self.start_identity,
            "initial_affine": _finite_values(self.initial_affine),
            "initial_physical": _finite_values(self.initial_physical),
            "final_physical": _finite_values(self.final_physical),
            "final_affine": _finite_values(self.final_affine),
            "predictions": _finite_values(self.predictions),
            "residuals": _finite_values(self.residuals),
            "block_contributions": contributions,
            "physical_jacobian": _finite_values(self.physical_jacobian),
            "affine_jacobian": _finite_values(self.affine_jacobian),
            "singular_values": _finite_values(self.singular_values),
            "rank": self.rank,
            "free_subspace_rank": self.free_subspace_rank,
            "condition_number": _finite_or_none(self.condition_number),
            "affine_covariance": _finite_values(self.affine_covariance),
            "degrees_of_freedom": self.degrees_of_freedom,
            "information_row_count": self.information_row_count,
            "free_coordinate_indices": list(self.free_coordinate_indices),
            "residual_variance": _finite_or_none(self.residual_variance),
            "covariance_status": self.covariance_status,
            "active_bounds": list(self.active_bounds),
            "solver_status": self.solver_status,
            "failure_reason": self.failure_reason,
            "trial_failure_code": self.trial_failure_code,
            "iterations": self.iterations,
            "non_evaluable_count": self.non_evaluable_count,
            "owner_failures": [failure.to_mapping() for failure in self.owner_failures],
            "owner_receipts": [receipt.to_mapping() for receipt in self.owner_receipts],
            "batch_outcomes": [outcome.to_mapping() for outcome in self.batch_outcomes],
            "cost": _finite_or_none(self.cost),
            "usable": self.usable,
            "optimizer_basin": self.optimizer_basin,
            "outcomes": self.outcomes,
        }


@dataclass(frozen=True, slots=True)
class EvaluatedFitRow:
    """One exact optimizer-row result; plotting consumes this table directly."""

    identity: str
    source_identity: str
    source_hash: str
    selector: str
    series_metadata_json: str
    observed: float
    predicted: float
    raw_difference: float
    scaled_residual: float
    unit: str
    basis: str
    evaluation_status: str
    candidate_parameter_identity: str

    def to_mapping(self) -> dict[str, object]:
        return {
            "identity": self.identity,
            "source_identity": self.source_identity,
            "source_hash": self.source_hash,
            "selector": self.selector,
            "series_metadata": json.loads(self.series_metadata_json),
            "observed": self.observed,
            "predicted": self.predicted,
            "raw_difference": self.raw_difference,
            "scaled_residual": self.scaled_residual,
            "unit": self.unit,
            "basis": self.basis,
            "evaluation_status": self.evaluation_status,
            "candidate_parameter_identity": self.candidate_parameter_identity,
        }


@dataclass(frozen=True, slots=True, init=False)
class FitResult:
    specification: FitSpecificationReceipt
    starts: tuple[FitStart, ...]
    artifact_fingerprint: str
    package_fingerprint: str
    accounting: FitAccounting
    digest: str

    def __init__(
        self,
        specification: FitSpecificationReceipt,
        starts: Sequence[FitStart],
        *,
        artifact_fingerprint: str,
        package_fingerprint: str,
    ) -> None:
        if not isinstance(specification, FitSpecificationReceipt):
            raise TypeError("specification must be a FitSpecificationReceipt")
        frozen_starts = tuple(starts)
        if not frozen_starts or any(not isinstance(item, FitStart) for item in frozen_starts):
            raise ValueError("starts must contain every native start")
        if tuple(item.index for item in frozen_starts) != tuple(range(len(frozen_starts))):
            raise ValueError("starts must retain deterministic native order")
        if any(
            tuple(outcome.index for outcome in start.batch_outcomes)
            != tuple(range(len(specification.block_sizes)))
            for start in frozen_starts
        ):
            raise ValueError("every start must retain ordered final batch outcomes")
        object.__setattr__(self, "specification", specification)
        object.__setattr__(self, "starts", frozen_starts)
        object.__setattr__(self, "artifact_fingerprint", artifact_fingerprint)
        object.__setattr__(self, "package_fingerprint", package_fingerprint)
        evaluated = sum(
            outcome.row_count
            for start in frozen_starts
            for outcome in start.batch_outcomes
            if outcome.status == "evaluated"
        )
        rejected = sum(
            outcome.row_count
            for start in frozen_starts
            for outcome in start.batch_outcomes
            if outcome.status == "rejected"
        )
        failed = sum(
            outcome.row_count
            for start in frozen_starts
            for outcome in start.batch_outcomes
            if outcome.status == "not_evaluated"
        )
        object.__setattr__(
            self,
            "accounting",
            FitAccounting(
                input_rows=specification.input_row_count,
                included_rows=len(specification.row_identities),
                dropped_rows=specification.dropped_row_count,
                attempted_starts=len(frozen_starts),
                accepted_starts=sum(start.usable for start in frozen_starts),
                failed_starts=sum(not start.usable for start in frozen_starts),
                evaluated_row_trials=evaluated,
                rejected_row_trials=rejected,
                failed_row_trials=failed,
                non_evaluable_trials=sum(start.non_evaluable_count for start in frozen_starts),
            ),
        )
        object.__setattr__(self, "digest", "")
        object.__setattr__(
            self,
            "digest",
            f"sha256:{hashlib.sha256(_canonical_json_bytes(self._payload())).hexdigest()}",
        )

    @classmethod
    def from_native(
        cls,
        specification: FitSpecificationReceipt,
        native: Mapping[str, Any],
        *,
        artifact_fingerprint: str,
        package_fingerprint: str,
    ) -> FitResult:
        starts = native.get("starts")
        if not isinstance(starts, Sequence) or isinstance(starts, (str, bytes)):
            raise TypeError("native fit starts must be a sequence")
        return cls(
            specification,
            tuple(FitStart.from_native(item) for item in starts),
            artifact_fingerprint=artifact_fingerprint,
            package_fingerprint=package_fingerprint,
        )

    @property
    def best_usable_start(self) -> FitStart | None:
        usable = tuple(start for start in self.starts if start.usable)
        return min(usable, key=lambda start: (start.cost, start.start_identity)) if usable else None

    @property
    def candidate_parameter_identity(self) -> str:
        accepted = self.best_usable_start
        if accepted is None:
            return ""
        payload = {
            "input_parameter_fingerprint": self.specification.parameters.fingerprint,
            "coordinates": [
                {
                    "parameter_identities": coordinate.parameter_identities,
                    "unit": coordinate.unit,
                    "value": value,
                }
                for coordinate, value in zip(
                    self.specification.coordinates, accepted.final_physical, strict=True
                )
            ],
        }
        return f"sha256:{hashlib.sha256(_canonical_json_bytes(payload)).hexdigest()}"

    @property
    def evaluated_rows(self) -> tuple[EvaluatedFitRow, ...]:
        accepted = self.best_usable_start
        if accepted is None:
            return ()
        row_count = len(self.specification.rows)
        if (
            len(accepted.predictions) != row_count
            or len(accepted.residuals) != row_count
            or len(accepted.physical_jacobian) != row_count * len(self.specification.coordinates)
            or any(outcome.status != "evaluated" for outcome in accepted.batch_outcomes)
            or any(
                not math.isfinite(value)
                for value in (
                    *accepted.predictions,
                    *accepted.residuals,
                    *accepted.physical_jacobian,
                )
            )
        ):
            return ()
        candidate = self.candidate_parameter_identity
        return tuple(
            EvaluatedFitRow(
                identity=row.identity,
                source_identity=row.source_identity,
                source_hash=row.source_hash,
                selector=row.selector,
                series_metadata_json=row.series_metadata_json,
                observed=row.observed,
                predicted=prediction,
                raw_difference=prediction - row.observed,
                scaled_residual=residual,
                unit=row.unit,
                basis=row.basis,
                evaluation_status="evaluated",
                candidate_parameter_identity=candidate,
            )
            for row, prediction, residual in zip(
                self.specification.rows,
                accepted.predictions,
                accepted.residuals,
                strict=True,
            )
        )

    @property
    def evaluated_row_table_hash(self) -> str:
        rows = self.evaluated_rows
        if not rows:
            return ""
        return (
            "sha256:"
            + hashlib.sha256(_canonical_json_bytes([row.to_mapping() for row in rows])).hexdigest()
        )

    @property
    def covariance_coordinate_identities(self) -> tuple[str, ...]:
        accepted = self.best_usable_start
        if accepted is None or accepted.covariance_status != "available":
            return ()
        return tuple(
            self.specification.coordinates[index].identity
            for index in accepted.free_coordinate_indices
        )

    @property
    def physical_covariance(self) -> tuple[float, ...]:
        accepted = self.best_usable_start
        if accepted is None or accepted.covariance_status != "available":
            return ()
        scales = tuple(
            self.specification.coordinates[index].scale
            for index in accepted.free_coordinate_indices
        )
        size = len(scales)
        return tuple(
            accepted.affine_covariance[row * size + column] * scales[row] * scales[column]
            for row in range(size)
            for column in range(size)
        )

    def fitted_values(self) -> tuple[tuple[str, float, str], ...]:
        accepted = self.best_usable_start
        if accepted is None:
            return ()
        values: list[tuple[str, float, str]] = []
        for coordinate, value in zip(
            self.specification.coordinates, accepted.final_physical, strict=True
        ):
            values.extend(
                (identity, value, coordinate.unit) for identity in coordinate.parameter_identities
            )
        return tuple(values)

    @property
    def outcomes(self) -> dict[str, dict[str, object]] | None:
        accepted = self.best_usable_start
        if accepted is None:
            return None
        outcomes: dict[str, dict[str, object]] = {
            name: dict(outcome) for name, outcome in accepted.outcomes.items()
        }
        basins = {
            start.optimizer_basin for start in self.starts if start.optimizer_basin is not None
        }
        if len(self.starts) > 1:
            usable_starts = sum(start.usable for start in self.starts)
            converged = usable_starts == len(self.starts) and len(basins) == 1
            outcomes["numerical_convergence"] = {
                "status": "converged" if converged else "not_converged",
                "reason": (
                    "all declared starts reached one optimizer basin"
                    if converged
                    else "declared starts did not all reach one optimizer basin"
                ),
                "attempted_starts": len(self.starts),
                "usable_starts": usable_starts,
                "distinct_basins": len(basins),
            }

        resolved = cast(Mapping[str, object], self.specification.to_mapping()["fit"])
        fit = cast(Mapping[str, object], resolved["fit"])
        acceptance = fit.get("acceptance")
        if isinstance(acceptance, Mapping):
            threshold = float(acceptance["maximum_cost"])
            passed = accepted.cost <= threshold
            outcomes["scientific_fit_acceptance"] = {
                "status": "accepted" if passed else "rejected",
                "reason": "best usable cost compared with declared maximum_cost",
                "metric": "half_sum_of_squared_scaled_residuals",
                "achieved": accepted.cost,
                "threshold": threshold,
            }
        return outcomes

    def _payload(self) -> dict[str, object]:
        accepted = self.best_usable_start
        return {
            "schema": "epcsaft.fit-result",
            "schema_version": 2,
            "specification_digest": self.specification.digest,
            "specification": self.specification.to_mapping(),
            "artifact_fingerprint": self.artifact_fingerprint,
            "package_fingerprint": self.package_fingerprint,
            "accounting": self.accounting.to_mapping(),
            "covariance_coordinate_identities": list(self.covariance_coordinate_identities),
            "physical_covariance": _finite_values(self.physical_covariance),
            "fitted_values": [
                {"identity": identity, "value": _finite_or_none(value), "unit": unit}
                for identity, value, unit in self.fitted_values()
            ],
            "candidate_parameter_identity": self.candidate_parameter_identity or None,
            "evaluated_row_table_hash": self.evaluated_row_table_hash or None,
            "evaluated_rows": [row.to_mapping() for row in self.evaluated_rows],
            "best_usable_start": None if accepted is None else accepted.index,
            "outcomes": self.outcomes,
            "starts": [start.to_mapping(self.specification.block_sizes) for start in self.starts],
        }

    def to_mapping(self) -> dict[str, object]:
        return {**self._payload(), "digest": self.digest}

    def to_json_bytes(self) -> bytes:
        return _canonical_json_bytes(self.to_mapping())

    def to_json(self, path: str | PathLike[str]) -> Path:
        target = Path(path)
        target.write_bytes(self.to_json_bytes() + b"\n")
        return target

    def export_candidate_parameters(self) -> Parameters:
        """Create a new provenance-linked candidate without granting adoption authority."""

        if not self.evaluated_rows or not self.evaluated_row_table_hash:
            raise RuntimeError("fit has no complete certified candidate point")
        values = {
            identity: unit_registry.Quantity(value, unit)
            for identity, value, unit in self.fitted_values()
        }
        if not values:
            raise RuntimeError("fit has no usable start")
        suffix = self.candidate_parameter_identity.removeprefix("sha256:")[:16]
        return self.specification.parameters.with_values(
            values,
            source_id=f"fit-candidate-{suffix}",
            locator=(
                f"FitResult {self.digest}; evaluated-row-table {self.evaluated_row_table_hash}"
            ),
        )


__all__ = (
    "BatchOutcome",
    "EvaluatedFitRow",
    "FitAccounting",
    "FitCoordinateReceipt",
    "FitResult",
    "FitRowReceipt",
    "FitSpecificationReceipt",
    "FitStart",
    "OwnerReceipt",
)
