from __future__ import annotations

import hashlib
import json
import math
from collections.abc import Iterable
from dataclasses import dataclass
from enum import StrEnum

from epcsaft import regression as _package

import epcsaft.regression._native as _native
from .parameter_regression import (
    ComponentParameterIdentity,
    FittedParameterDiagnostic,
    GeneralJacobianDiagnostics,
    ModelParameterIdentity,
    ObservationPartition,
    PairParameterIdentity,
    ParameterCoordinate,
    SourceDescriptor,
    _complete_failure_reasons,
)


def _nonempty(value: str, field: str) -> None:
    if type(value) is not str or not value.strip():
        raise ValueError(f"{field} must be a nonempty string")


def _finite(value: float, field: str, *, positive: bool = False) -> None:
    if type(value) not in (int, float) or not math.isfinite(value):
        raise ValueError(f"{field} must be finite")
    if positive and value <= 0.0:
        raise ValueError(f"{field} must be positive")


def _sha256(value: str, field: str) -> None:
    _nonempty(value, field)
    if not value.startswith("sha256:"):
        raise ValueError(f"{field} must begin with 'sha256:'")
    body = value.removeprefix("sha256:")
    if len(body) != 64 or any(character not in "0123456789abcdef" for character in body):
        raise ValueError(f"{field} must be a lowercase SHA-256 identity")


def _installed_artifact_identity(value: str, field: str) -> None:
    _nonempty(value, field)
    fields = value.split(";")
    for label in ("RECORD", "HEADER"):
        matches = tuple(
            item.removeprefix(f"{label}=") for item in fields if item.startswith(f"{label}=")
        )
        if len(matches) != 1:
            raise ValueError(f"{field} must contain exactly one {label} identity")
        _sha256(matches[0], f"{field} {label}")


def _parameter_id(parameter: ParameterCoordinate) -> str:
    if isinstance(parameter.identity, ComponentParameterIdentity):
        identity = "component"
    elif isinstance(parameter.identity, PairParameterIdentity):
        identity = "unordered_component_pair"
    elif isinstance(parameter.identity, ModelParameterIdentity):
        identity = "model"
    else:
        raise TypeError("parameter identity is unsupported")
    components = ",".join(parameter.identity.canonical_component_ids)
    return f"{parameter.family.value};{identity};{components}"


class PositiveObservationTransform(StrEnum):
    IDENTITY = "identity"
    NATURAL_LOG = "natural_log"


@dataclass(frozen=True, slots=True)
class PositiveScalarObservation:
    row_id: str
    state_id: str
    state_schema_id: str
    source_id: str
    source_locator: str
    primitive_id: str
    primitive_unit: str
    transform: PositiveObservationTransform
    reference_id: str
    reference_fingerprint: str
    observed_value: float
    residual_scale: float
    residual_scale_unit: str
    partition: ObservationPartition

    def __post_init__(self) -> None:
        for field in (
            "row_id",
            "state_id",
            "state_schema_id",
            "source_id",
            "source_locator",
            "primitive_id",
            "primitive_unit",
            "reference_id",
            "residual_scale_unit",
        ):
            _nonempty(getattr(self, field), field)
        if not isinstance(self.transform, PositiveObservationTransform):
            raise TypeError("transform must be a PositiveObservationTransform")
        if not isinstance(self.partition, ObservationPartition):
            raise TypeError("partition must be an ObservationPartition")
        _sha256(self.reference_fingerprint, "reference fingerprint")
        _finite(self.observed_value, "observed value", positive=True)
        _finite(self.residual_scale, "residual scale", positive=True)
        expected_scale_unit = (
            "1"
            if self.transform is PositiveObservationTransform.NATURAL_LOG
            else self.primitive_unit
        )
        if self.residual_scale_unit != expected_scale_unit:
            raise ValueError("residual scale unit is incompatible with the observation transform")


def canonical_positive_dataset_sha256(
    observations: Iterable[PositiveScalarObservation],
) -> str:
    rows = sorted(
        (
            {
                "row_id": row.row_id,
                "state_id": row.state_id,
                "state_schema_id": row.state_schema_id,
                "source_id": row.source_id,
                "source_locator": row.source_locator,
                "primitive_id": row.primitive_id,
                "primitive_unit": row.primitive_unit,
                "transform": row.transform.value,
                "reference_id": row.reference_id,
                "reference_fingerprint": row.reference_fingerprint,
                "observed_value": row.observed_value,
                "residual_scale": row.residual_scale,
                "residual_scale_unit": row.residual_scale_unit,
                "partition": row.partition.value,
            }
            for row in observations
        ),
        key=lambda row: row["row_id"],
    )
    encoded = json.dumps(rows, allow_nan=False, separators=(",", ":"), sort_keys=True).encode(
        "utf-8"
    )
    return hashlib.sha256(encoded).hexdigest()


@dataclass(frozen=True, slots=True)
class EvaluatorContract:
    evaluator_identity: str
    capability_id: str
    capability_fingerprint: str
    provider_artifact_identity: str
    owner_artifact_identity: str
    contract_fingerprint: str
    model_fingerprint: str
    artifact_identity: str

    def __post_init__(self) -> None:
        for field in (
            "evaluator_identity",
            "capability_id",
            "provider_artifact_identity",
            "owner_artifact_identity",
        ):
            _nonempty(getattr(self, field), field)
        for field in (
            "capability_fingerprint",
            "contract_fingerprint",
            "model_fingerprint",
            "artifact_identity",
        ):
            _sha256(getattr(self, field), field.replace("_", " "))
        for field in ("provider_artifact_identity", "owner_artifact_identity"):
            _installed_artifact_identity(getattr(self, field), field.replace("_", " "))


@dataclass(frozen=True, slots=True)
class PositiveEvaluatorProblem:
    sources: tuple[SourceDescriptor, ...]
    parameters: tuple[ParameterCoordinate, ...]
    parameter_ids: tuple[str, ...]
    start_vectors: tuple[tuple[float, ...], ...]
    observations: tuple[PositiveScalarObservation, ...]
    evaluator: EvaluatorContract
    maximum_condition_number: float
    maximum_iterations: int
    maximum_solver_time_seconds: float
    function_tolerance: float
    gradient_tolerance: float
    parameter_tolerance: float
    confirmation_parameter_scaled_max_delta: float
    confirmation_cost_relative_delta: float

    def __post_init__(self) -> None:
        if not self.sources or not self.parameters or not self.observations:
            raise ValueError("sources, parameters, and observations must be nonempty")
        if not isinstance(self.evaluator, EvaluatorContract):
            raise TypeError("evaluator must be an EvaluatorContract")
        if (
            type(self.parameter_ids) is not tuple
            or len(self.parameter_ids) != len(self.parameters)
            or len(set(self.parameter_ids)) != len(self.parameter_ids)
        ):
            raise ValueError("parameter_ids must uniquely match the ordered parameter count")
        for parameter_id in self.parameter_ids:
            _nonempty(parameter_id, "parameter_id")
        expected_parameter_ids = tuple(_parameter_id(parameter) for parameter in self.parameters)
        if self.parameter_ids != expected_parameter_ids:
            raise ValueError("parameter_ids must match the ordered typed parameter identities")
        if len(self.start_vectors) < 2:
            raise ValueError("start_vectors must contain a primary and at least one confirmation")
        for vector in self.start_vectors:
            if type(vector) is not tuple or len(vector) != len(self.parameters):
                raise ValueError("every start vector must match the parameter count")
            for coordinate, value in zip(self.parameters, vector, strict=True):
                _finite(value, "parameter start")
                if not coordinate.lower_bound <= value <= coordinate.upper_bound:
                    raise ValueError("every parameter start must lie within its bounds")
        row_ids = tuple(row.row_id for row in self.observations)
        if len(set(row_ids)) != len(row_ids):
            raise ValueError("positive-observation row identities must be unique")
        source_ids = tuple(source.source_id for source in self.sources)
        if len(set(source_ids)) != len(source_ids):
            raise ValueError("source identities must be unique")
        source_map = {source.source_id: source for source in self.sources}
        for row in self.observations:
            if row.partition is not ObservationPartition.TRAINING:
                raise ValueError(
                    "positive evaluator observations must all be training rows; "
                    "reserved partitions require a separately admitted workflow"
                )
            if row.source_id not in source_map:
                raise ValueError(f"row {row.row_id!r} references an unknown source")
        for source in self.sources:
            rows = tuple(row for row in self.observations if row.source_id == source.source_id)
            if not rows:
                raise ValueError(f"source {source.source_id!r} has no observations")
            if canonical_positive_dataset_sha256(rows) != source.canonical_dataset_sha256:
                raise ValueError(
                    f"source {source.source_id!r} canonical dataset SHA-256 does not match its rows"
                )
        if not self.training_observations:
            raise ValueError("at least one training observation is required")
        if len(self.training_observations) < len(self.parameters):
            raise ValueError("training residual count must be at least the parameter count")
        if len({parameter.provider_parameter_fingerprint for parameter in self.parameters}) != 1:
            raise ValueError("parameters must bind one Provider parameter fingerprint")
        if len({parameter.provider_topology_fingerprint for parameter in self.parameters}) != 1:
            raise ValueError("parameters must bind one Provider topology fingerprint")
        if any(
            parameter.capability_id != self.evaluator.capability_id for parameter in self.parameters
        ):
            raise ValueError("parameter capability_id must match the evaluator capability_id")
        _finite(
            self.maximum_condition_number,
            "maximum condition number",
            positive=True,
        )
        if type(self.maximum_iterations) is not int or self.maximum_iterations <= 0:
            raise ValueError("maximum_iterations must be a positive integer")
        for field in (
            "maximum_solver_time_seconds",
            "function_tolerance",
            "gradient_tolerance",
            "parameter_tolerance",
            "confirmation_parameter_scaled_max_delta",
            "confirmation_cost_relative_delta",
        ):
            _finite(getattr(self, field), field, positive=True)

    @property
    def training_observations(self) -> tuple[PositiveScalarObservation, ...]:
        return tuple(
            row for row in self.observations if row.partition is ObservationPartition.TRAINING
        )

    @property
    def solver_start_vectors(self) -> tuple[tuple[float, ...], ...]:
        return tuple(
            tuple(
                coordinate.transform.to_solver(value)
                for coordinate, value in zip(self.parameters, vector, strict=True)
            )
            for vector in self.start_vectors
        )


@dataclass(frozen=True, slots=True)
class PositiveEvaluatorCapability:
    evaluator_identity: str
    capability_id: str
    capability_fingerprint: str
    provider_artifact_identity: str
    owner_artifact_identity: str
    contract_fingerprint: str
    model_fingerprint: str
    artifact_identity: str
    provider_sdk_capsule_name: str
    provider_sdk_abi_version: int
    provider_sdk_table_size: int
    provider_sdk_result_size: int
    provider_sdk_mixture_result_size: int
    provider_sdk_neutral_reference_result_size: int
    provider_sdk_neutral_reference_derivative_result_size: int
    provider_sdk_reacting_phase_parameter_result_size: int
    single_thread_non_reentrant: bool
    value_only_avoids_derivative_work: bool


@dataclass(frozen=True, slots=True)
class ComposedPositiveRowDiagnostic:
    row_id: str
    partition: str
    state_id: str
    state_schema_id: str
    source_id: str
    source_locator: str
    primitive_id: str
    primitive_unit: str
    transform: str
    reference_id: str
    reference_fingerprint: str
    observed_value: float
    modeled_value: float
    scaled_residual: float
    scaled_solver_jacobian: tuple[float, ...]
    solver_status: str
    numerical_status: str
    physical_status: str
    derivative_status: str
    chart_topology: str
    provider_topology_fingerprint: str
    kkt_dimension: int
    kkt_rank: int
    kkt_condition_number_inf: float


def _native_payload(problem: PositiveEvaluatorProblem) -> tuple[object, ...]:
    contract = problem.evaluator
    return (
        (
            contract.evaluator_identity,
            contract.capability_id,
            contract.capability_fingerprint,
            contract.provider_artifact_identity,
            contract.owner_artifact_identity,
            contract.contract_fingerprint,
            contract.model_fingerprint,
            problem.parameters[0].provider_parameter_fingerprint,
            problem.parameters[0].provider_topology_fingerprint,
            contract.artifact_identity,
        ),
        tuple(
            (
                problem.parameter_ids[index],
                parameter.unit,
                parameter.transform.origin,
                parameter.transform.scale,
                parameter.lower_bound,
                parameter.upper_bound,
            )
            for index, parameter in enumerate(problem.parameters)
        ),
        problem.start_vectors,
        tuple(
            (
                row.row_id,
                row.partition.value,
                row.state_id,
                row.state_schema_id,
                row.source_id,
                row.source_locator,
                row.primitive_id,
                row.primitive_unit,
                row.transform.value,
                row.reference_id,
                row.reference_fingerprint,
                row.observed_value,
                row.residual_scale,
            )
            for row in problem.observations
        ),
        (
            problem.maximum_condition_number,
            problem.maximum_iterations,
            problem.maximum_solver_time_seconds,
            problem.function_tolerance,
            problem.gradient_tolerance,
            problem.parameter_tolerance,
            problem.confirmation_parameter_scaled_max_delta,
            problem.confirmation_cost_relative_delta,
        ),
    )


def fit_positive_observations(
    problem: PositiveEvaluatorProblem,
    evaluator_handle: object,
) -> _package.RegressionResult:
    from .result import RegressionResult

    if not isinstance(problem, PositiveEvaluatorProblem):
        raise TypeError("problem must be a PositiveEvaluatorProblem")
    native = _native.solve_evaluator(evaluator_handle, _native_payload(problem))
    if native[21]:
        raise RuntimeError(str(native[21]))
    physical_parameters = tuple(float(value) for value in native[5])
    bound_distances = tuple(float(value) for value in native[6])
    active_bounds = tuple(str(value) for value in native[7])
    row_records = native[20]
    rows = tuple(
        ComposedPositiveRowDiagnostic(
            row_id=observation.row_id,
            partition=observation.partition.value,
            state_id=observation.state_id,
            state_schema_id=observation.state_schema_id,
            source_id=observation.source_id,
            source_locator=observation.source_locator,
            primitive_id=observation.primitive_id,
            primitive_unit=observation.primitive_unit,
            transform=observation.transform.value,
            reference_id=observation.reference_id,
            reference_fingerprint=observation.reference_fingerprint,
            observed_value=observation.observed_value,
            modeled_value=float(record[0]),
            scaled_residual=float(record[1]),
            scaled_solver_jacobian=tuple(record[2]),
            solver_status=str(record[3]),
            numerical_status=str(record[4]),
            physical_status=str(record[5]),
            derivative_status=str(record[6]),
            chart_topology=str(record[7]),
            provider_topology_fingerprint=str(record[8]),
            kkt_dimension=int(record[9]),
            kkt_rank=int(record[10]),
            kkt_condition_number_inf=float(record[11]),
        )
        for observation, record in zip(problem.observations, row_records, strict=True)
    )
    parameter_diagnostics = tuple(
        FittedParameterDiagnostic(
            family=coordinate.family,
            component_ids=coordinate.identity.canonical_component_ids,
            unit=coordinate.unit,
            transform_origin=coordinate.transform.origin,
            transform_scale=coordinate.transform.scale,
            start=problem.start_vectors[0][index],
            final=physical_parameters[index],
            movement=physical_parameters[index] - problem.start_vectors[0][index],
            lower_bound=coordinate.lower_bound,
            upper_bound=coordinate.upper_bound,
            active_bound_distance=bound_distances[index],
            active_bound=active_bounds[index] or None,
        )
        for index, coordinate in enumerate(problem.parameters)
    )
    jacobian = GeneralJacobianDiagnostics(
        residual_count=int(native[23]),
        variable_count=int(native[22]),
        full_singular_values=tuple(native[10]),
        full_rank=int(native[11]),
        full_condition_number=float(native[12]),
        projected_parameter_singular_values=tuple(native[13]),
        projected_parameter_rank=int(native[14]),
        projected_parameter_condition_number=float(native[15]),
    )
    solver_converged = native[0] == "CONVERGENCE" and bool(native[1]) and not native[21]
    confirmations_usable = bool(native[19])
    numerically_converged = (
        solver_converged
        and confirmations_usable
        and native[17] <= problem.confirmation_parameter_scaled_max_delta
        and native[18] <= problem.confirmation_cost_relative_delta
        and jacobian.full_rank == jacobian.variable_count
        and jacobian.projected_parameter_rank == len(problem.parameters)
        and math.isfinite(jacobian.full_condition_number)
        and jacobian.full_condition_number <= problem.maximum_condition_number
        and math.isfinite(jacobian.projected_parameter_condition_number)
        and jacobian.projected_parameter_condition_number <= problem.maximum_condition_number
    )
    workflow_valid = len(rows) == len(problem.observations)
    failure_reasons = _complete_failure_reasons(
        (),
        termination=str(native[0]),
        solution_usable=bool(native[1]),
        solver_converged=solver_converged,
        confirmations_usable=confirmations_usable,
        confirmation_parameter_delta=float(native[17]),
        confirmation_parameter_limit=problem.confirmation_parameter_scaled_max_delta,
        confirmation_cost_delta=float(native[18]),
        confirmation_cost_limit=problem.confirmation_cost_relative_delta,
        jacobian=jacobian,
        parameter_count=len(problem.parameters),
        maximum_condition_number=problem.maximum_condition_number,
        workflow_valid=workflow_valid,
    )
    capability_record = native[26]
    capability = PositiveEvaluatorCapability(
        evaluator_identity=problem.evaluator.evaluator_identity,
        capability_id=problem.evaluator.capability_id,
        capability_fingerprint=problem.evaluator.capability_fingerprint,
        provider_artifact_identity=problem.evaluator.provider_artifact_identity,
        owner_artifact_identity=problem.evaluator.owner_artifact_identity,
        contract_fingerprint=problem.evaluator.contract_fingerprint,
        model_fingerprint=problem.evaluator.model_fingerprint,
        artifact_identity=problem.evaluator.artifact_identity,
        provider_sdk_capsule_name=str(capability_record[0]),
        provider_sdk_abi_version=int(capability_record[1]),
        provider_sdk_table_size=int(capability_record[2]),
        provider_sdk_result_size=int(capability_record[3]),
        provider_sdk_mixture_result_size=int(capability_record[4]),
        provider_sdk_neutral_reference_result_size=int(capability_record[5]),
        provider_sdk_neutral_reference_derivative_result_size=int(capability_record[6]),
        provider_sdk_reacting_phase_parameter_result_size=int(capability_record[7]),
        single_thread_non_reentrant=bool(capability_record[8]),
        value_only_avoids_derivative_work=bool(capability_record[9]),
    )
    return RegressionResult(
        problem=problem,
        capabilities=(capability,),
        preparation_fingerprint=None,
        provider_parameter_fingerprint=str(native[27]),
        provider_topology_fingerprint=(problem.parameters[0].provider_topology_fingerprint),
        solver_converged=solver_converged,
        numerically_converged=numerically_converged,
        workflow_valid=workflow_valid,
        physical_status=(
            "UPSTREAM_CERTIFIED_ALL_ROWS" if workflow_valid else "UPSTREAM_ROW_CERTIFICATE_REJECTED"
        ),
        scientific_status="NOT_ADJUDICATED_SOURCE_BOUND_FIT_ONLY",
        predictive_status="NOT_ADJUDICATED_NO_APPROVED_HELD_OUT_CUTOFF",
        authority_status="NO_AUTHORITY_CHANGE_REGRESSION_RESULT_ONLY",
        termination=str(native[0]),
        solution_usable=bool(native[1]),
        initial_cost=float(native[2]),
        final_cost=float(native[3]),
        iterations=int(native[4]),
        residual_evaluation_count=int(native[24]),
        jacobian_evaluation_count=int(native[25]),
        parameters=parameter_diagnostics,
        jacobian=jacobian,
        rows=rows,
        confirmation_count=int(native[16]),
        confirmation_parameter_scaled_max_delta=float(native[17]),
        confirmation_cost_relative_max_delta=float(native[18]),
        confirmations_usable=confirmations_usable,
        training_row_count=len(problem.training_observations),
        held_out_row_count=0,
        stress_row_count=0,
        evaluated_row_count=len(rows),
        skipped_row_count=0,
        failed_row_count=0,
        failure_reasons=failure_reasons,
    )


__all__ = (
    "ComposedPositiveRowDiagnostic",
    "EvaluatorContract",
    "PositiveEvaluatorCapability",
    "PositiveEvaluatorProblem",
    "PositiveObservationTransform",
    "PositiveScalarObservation",
    "canonical_positive_dataset_sha256",
    "fit_positive_observations",
)
