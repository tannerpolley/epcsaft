from __future__ import annotations

import hashlib
import json
import math
from collections.abc import Iterable, Mapping
from dataclasses import asdict, dataclass, fields, is_dataclass, replace
from enum import StrEnum

from . import _native
from .parameter_regression import (
    AffineParameterTransform,
    AqueousKijMeanIonicActivityObservation,
    AssociationParameterIdentity,
    ComponentParameterIdentity,
    FixedCompositionMixtureDensityObservation,
    FixedCompositionVleObservation,
    FixedTopologyAssociationCapability,
    GeneralJacobianDiagnostics,
    IonSolvationKijObservation,
    MeanIonicActivityObservation,
    ModelParameterIdentity,
    ObservationPartition,
    PairParameterIdentity,
    ParameterCapability,
    ParameterCoordinate,
    ParameterFamily,
    PureDensityObservation,
    PureSaturationObservation,
    PureVaporPressureObservation,
    RegressionObservation,
    RegressionProblem,
    RelativePermittivityRatioObservation,
    SolvationGibbsObservation,
    SourceDescriptor,
    _evaluate_parameters,
    _evaluate_parameters_at_start,
    _require_finite,
    _require_nonempty_string,
    _require_sha256,
    _uses_fixed_topology_descriptor,
    canonical_dataset_sha256,
    fit_parameters,
    parameter_capabilities,
)
from .result import RegressionResult, _installed_eos_artifact_identity


class AcquisitionClass(StrEnum):
    DIRECT_MEASUREMENT = "direct_measurement"
    DIGITIZED_VALUE = "digitized_value"
    DATABASE_RECORD = "database_record"
    AUTHOR_CORRELATION = "author_correlation"
    RECONSTRUCTED_CORRELATION = "reconstructed_correlation"


@dataclass(frozen=True, slots=True)
class CorrelationProvenance:
    equation: str
    coefficients: tuple[tuple[str, float], ...]
    units: str
    validity_interval: str
    sampling_fields: tuple[str, ...]
    sampling_grid: tuple[tuple[float, ...], ...]
    transformation_record: str

    def __post_init__(self) -> None:
        for name in ("equation", "units", "validity_interval", "transformation_record"):
            _require_nonempty_string(getattr(self, name), name)
        if not self.coefficients or not self.sampling_fields or not self.sampling_grid:
            raise ValueError(
                "correlation coefficients, sampling_fields, and sampling_grid must be nonempty"
            )
        if any(type(name) is not str or not name.strip() for name, _ in self.coefficients) or len(
            {name for name, _ in self.coefficients}
        ) != len(self.coefficients):
            raise ValueError("correlation coefficient names must be unique and nonempty")
        if any(
            type(field) is not str or not field.strip() for field in self.sampling_fields
        ) or len(set(self.sampling_fields)) != len(self.sampling_fields):
            raise ValueError("correlation sampling_fields must be unique and nonempty")
        if any(
            type(value) not in (int, float) or not math.isfinite(value)
            for _, value in self.coefficients
        ) or any(
            type(point) is not tuple
            or len(point) != len(self.sampling_fields)
            or any(type(value) not in (int, float) or not math.isfinite(value) for value in point)
            for point in self.sampling_grid
        ):
            raise ValueError("correlation coefficients and sampling_grid must be finite")


@dataclass(frozen=True, slots=True)
class RowProvenance:
    acquisition: AcquisitionClass
    duplicate_decision: str
    exclusion_decision: str
    critical_region_decision: str
    censoring_decision: str
    outlier_decision: str
    correlation: CorrelationProvenance | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.acquisition, AcquisitionClass):
            raise TypeError("acquisition must be an AcquisitionClass")
        for name in (
            "duplicate_decision",
            "exclusion_decision",
            "critical_region_decision",
            "censoring_decision",
            "outlier_decision",
        ):
            _require_nonempty_string(getattr(self, name), name)
        correlation_backed = self.acquisition in (
            AcquisitionClass.AUTHOR_CORRELATION,
            AcquisitionClass.RECONSTRUCTED_CORRELATION,
        )
        if self.correlation is not None and not isinstance(self.correlation, CorrelationProvenance):
            raise TypeError("correlation must be a CorrelationProvenance")
        if correlation_backed != (self.correlation is not None):
            raise ValueError(
                "correlation acquisition classes require correlation provenance, "
                "and other acquisition classes forbid it"
            )


@dataclass(frozen=True, slots=True)
class ObjectiveContract:
    residual_family: str
    interpretation: str
    row_weighting: str
    covariance_interpretation: str
    loss: str
    loss_parameters: tuple[tuple[str, float], ...]
    failed_row_treatment: str

    def __post_init__(self) -> None:
        for name in (
            "residual_family",
            "interpretation",
            "row_weighting",
            "covariance_interpretation",
            "failed_row_treatment",
        ):
            _require_nonempty_string(getattr(self, name), name)
        if self.loss != "squared":
            raise ValueError("the current Ceres problem supports only squared loss")
        if self.loss_parameters:
            raise ValueError("squared loss has no loss parameters")
        if self.interpretation != "native_scaled_least_squares":
            raise ValueError("interpretation must be 'native_scaled_least_squares'")
        if self.row_weighting != "observation_residual_scales":
            raise ValueError("row_weighting must be 'observation_residual_scales'")
        if self.covariance_interpretation != "independent_no_covariance":
            raise ValueError("covariance_interpretation must be 'independent_no_covariance'")
        if self.failed_row_treatment != "fail_fit":
            raise ValueError("failed_row_treatment must be 'fail_fit'")


@dataclass(frozen=True, slots=True)
class SourceInput:
    source_id: str
    citation: str
    durable_locator: str
    source_artifact_sha256: str
    transformation_record: str
    units_and_bases: str
    use_basis: str
    residual_scale_rationale: str

    def __post_init__(self) -> None:
        for name in (
            "source_id",
            "citation",
            "durable_locator",
            "transformation_record",
            "units_and_bases",
            "use_basis",
            "residual_scale_rationale",
        ):
            _require_nonempty_string(getattr(self, name), name)
        _require_sha256(self.source_artifact_sha256, "source_artifact_sha256")


def _json_value(value: object) -> object:
    if isinstance(value, StrEnum):
        return value.value
    if is_dataclass(value):
        return {key: _json_value(item) for key, item in asdict(value).items()}
    if isinstance(value, Mapping):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json_value(item) for item in value]
    return value


def _canonical_sha256(value: object) -> str:
    encoded = json.dumps(
        _json_value(value), sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode()
    return hashlib.sha256(encoded).hexdigest()


_TUPLE_FIELDS = {
    "component_ids",
    "mole_fractions",
    "active_pair_component_ids",
    "fixed_k_ij",
    "chemical_potential_scales",
    "liquid_volume_bounds_m3_per_mol",
    "vapor_volume_bounds_m3_per_mol",
    "volume_bounds_m3_per_mol",
}
_OBSERVATION_TYPES = (
    FixedCompositionVleObservation,
    PureSaturationObservation,
    PureVaporPressureObservation,
    PureDensityObservation,
    FixedCompositionMixtureDensityObservation,
    MeanIonicActivityObservation,
    AqueousKijMeanIonicActivityObservation,
    IonSolvationKijObservation,
    SolvationGibbsObservation,
    RelativePermittivityRatioObservation,
)


@dataclass(frozen=True, slots=True)
class ObservationDataset:
    source: SourceInput
    observations: tuple[RegressionObservation, ...]
    objective: ObjectiveContract
    row_provenance: tuple[tuple[str, RowProvenance], ...]
    provenance_sha256: str

    def __post_init__(self) -> None:
        if not isinstance(self.source, SourceInput):
            raise TypeError("source must be a SourceInput")
        if not isinstance(self.objective, ObjectiveContract):
            raise TypeError("objective must be an ObjectiveContract")
        if type(self.observations) is not tuple or any(
            type(row) not in _OBSERVATION_TYPES for row in self.observations
        ):
            raise TypeError("observations must be a tuple of fit-ready direct-EOS rows")
        if type(self.row_provenance) is not tuple:
            raise TypeError("row_provenance must be a tuple")
        if not self.observations:
            raise ValueError("dataset observations must be nonempty")
        _validate_objective(self)
        row_ids = tuple(row.row_id for row in self.observations)
        if len(set(row_ids)) != len(row_ids):
            raise ValueError("dataset row IDs must be unique")
        provenance_ids = tuple(row_id for row_id, _ in self.row_provenance)
        if len(set(provenance_ids)) != len(provenance_ids) or set(provenance_ids) != set(row_ids):
            raise ValueError("row provenance must match every observation row exactly")
        if any(not isinstance(provenance, RowProvenance) for _, provenance in self.row_provenance):
            raise TypeError("every row provenance value must be a RowProvenance")
        provenance_by_id = dict(self.row_provenance)
        correlations = {
            provenance.correlation
            for provenance in provenance_by_id.values()
            if provenance.correlation is not None
        }
        for correlation in correlations:
            correlated_rows = tuple(
                row
                for row in self.observations
                if provenance_by_id[row.row_id].correlation == correlation
            )
            if any(
                not hasattr(row, field)
                for row in correlated_rows
                for field in correlation.sampling_fields
            ):
                raise ValueError(
                    "correlation sampling_fields must name fields on every transformed row"
                )
            observed_grid = tuple(
                tuple(getattr(row, field) for field in correlation.sampling_fields)
                for row in correlated_rows
            )
            if any(
                type(value) not in (int, float) or not math.isfinite(value)
                for point in observed_grid
                for value in point
            ):
                raise ValueError("correlation sampling_fields must name finite numeric row fields")
            if observed_grid != correlation.sampling_grid:
                raise ValueError(
                    "correlation sampling_grid must exactly match the ordered "
                    "temperatures of its transformed rows"
                )
        expected = _canonical_sha256(
            {
                "source": self.source,
                "objective": self.objective,
                "rows": tuple((row_id, provenance) for row_id, provenance in self.row_provenance),
                "canonical_dataset_sha256": canonical_dataset_sha256(self.observations),
            }
        )
        if self.provenance_sha256 != expected:
            raise ValueError("provenance_sha256 does not match the validated dataset")

    @classmethod
    def from_records(
        cls,
        observation_type: type,
        records: Iterable[Mapping[str, object]],
        *,
        source: SourceInput,
        objective: ObjectiveContract,
        row_provenance: Mapping[str, RowProvenance],
    ) -> ObservationDataset:
        if observation_type not in _OBSERVATION_TYPES:
            raise ValueError("observation_type is not a fit-ready direct-EOS contract")
        allowed = {field.name for field in fields(observation_type)}
        observations: list[RegressionObservation] = []
        for index, record in enumerate(records):
            if not isinstance(record, Mapping):
                raise TypeError(f"row {index}: record must be a mapping")
            row_label = str(record.get("row_id", f"row {index}"))
            extra = sorted(set(record) - allowed)
            missing = sorted(allowed - set(record))
            if extra or missing:
                detail = []
                if missing:
                    detail.append(f"missing {', '.join(missing)}")
                if extra:
                    detail.append(f"unexpected {', '.join(extra)}")
                raise ValueError(f"{row_label}: {'; '.join(detail)}")
            values = dict(record)
            for name in _TUPLE_FIELDS.intersection(values):
                if isinstance(values[name], list):
                    values[name] = tuple(values[name])
            if isinstance(values.get("partition"), str):
                try:
                    values["partition"] = ObservationPartition(values["partition"])
                except ValueError as exc:
                    raise ValueError(f"{row_label}: invalid partition") from exc
            try:
                observations.append(observation_type(**values))
            except (TypeError, ValueError) as exc:
                raise ValueError(f"{row_label}: {exc}") from exc
        rows = tuple(observations)
        mismatched_sources = tuple(row.row_id for row in rows if row.source_id != source.source_id)
        if mismatched_sources:
            raise ValueError(
                "rows reference a source_id other than the dataset source: "
                + ", ".join(mismatched_sources)
            )
        provenance = tuple(sorted(row_provenance.items()))
        payload = {
            "source": source,
            "objective": objective,
            "rows": provenance,
            "canonical_dataset_sha256": canonical_dataset_sha256(rows),
        }
        return cls(source, rows, objective, provenance, _canonical_sha256(payload))


@dataclass(frozen=True, slots=True)
class SolverControls:
    maximum_iterations: int
    maximum_solver_time_seconds: float
    function_tolerance: float
    gradient_tolerance: float
    parameter_tolerance: float

    def __post_init__(self) -> None:
        if type(self.maximum_iterations) is not int or self.maximum_iterations <= 0:
            raise ValueError("maximum_iterations must be positive")
        for name in (
            "maximum_solver_time_seconds",
            "function_tolerance",
            "gradient_tolerance",
            "parameter_tolerance",
        ):
            _require_finite(getattr(self, name), name, positive=True)


@dataclass(frozen=True, slots=True)
class RankControls:
    maximum_condition_number: float

    def __post_init__(self) -> None:
        _require_finite(
            self.maximum_condition_number,
            "maximum_condition_number",
            positive=True,
        )


@dataclass(frozen=True, slots=True)
class ConfirmationControls:
    parameter_scaled_max_delta: float
    cost_relative_delta: float

    def __post_init__(self) -> None:
        _require_finite(
            self.parameter_scaled_max_delta,
            "parameter_scaled_max_delta",
            positive=True,
        )
        _require_finite(
            self.cost_relative_delta,
            "cost_relative_delta",
            positive=True,
        )


@dataclass(frozen=True, slots=True)
class ParameterRequest:
    family: ParameterFamily
    identity: (
        PairParameterIdentity
        | ComponentParameterIdentity
        | ModelParameterIdentity
        | AssociationParameterIdentity
    )
    transform: AffineParameterTransform
    lower_bound: float
    upper_bound: float

    def __post_init__(self) -> None:
        if not isinstance(self.family, ParameterFamily):
            raise TypeError("family must be a ParameterFamily")
        if type(self.lower_bound) not in (int, float) or type(self.upper_bound) not in (
            int,
            float,
        ):
            raise ValueError("parameter bounds must be finite")
        if not (
            math.isfinite(self.lower_bound)
            and math.isfinite(self.upper_bound)
            and self.lower_bound < self.upper_bound
        ):
            raise ValueError("parameter bounds must be finite and increasing")


def _unit(value: str) -> str:
    units = {"dimensionless": "1", "angstrom": "angstrom", "kelvin": "K"}
    if value not in units:
        raise ValueError(f"installed capability uses unsupported unit {value!r}")
    return units[value]


def _identity_matches(request: ParameterRequest, capability: ParameterCapability) -> bool:
    identity = request.identity.canonical_component_ids
    if isinstance(request.identity, ModelParameterIdentity):
        return capability.identity_shape == "model"
    if isinstance(request.identity, PairParameterIdentity):
        return identity == tuple(sorted(capability.active_component_ids))
    active = capability.active_component_ids or capability.component_ids
    return identity == active


def _observation_contract(observations: tuple[RegressionObservation, ...]) -> str:
    observation = observations[0]
    if isinstance(
        observation,
        (
            FixedCompositionVleObservation,
            PureSaturationObservation,
            PureVaporPressureObservation,
            PureDensityObservation,
            FixedCompositionMixtureDensityObservation,
        ),
    ):
        return "fixed_composition_helmholtz_phase"
    if isinstance(
        observation,
        (MeanIonicActivityObservation, AqueousKijMeanIonicActivityObservation),
    ):
        return "aqueous_mean_ionic_activity"
    if isinstance(observation, (SolvationGibbsObservation, IonSolvationKijObservation)):
        return "ion_solvation_gibbs"
    if isinstance(observation, RelativePermittivityRatioObservation):
        return "relative_permittivity_ratio"
    raise TypeError("unsupported regression observation contract")


def _resolve_coordinates(
    model: object,
    requests: tuple[ParameterRequest, ...],
    observations: tuple[RegressionObservation, ...],
) -> tuple[ParameterCoordinate, ...]:
    advertised = tuple(
        item for item in parameter_capabilities(model) if isinstance(item, ParameterCapability)
    )
    derivative_ready = tuple(item for item in advertised if item.installed_ready)
    matches = tuple(
        item
        for item in parameter_capabilities(model)
        if isinstance(item, FixedTopologyAssociationCapability) and item.installed_ready
    )
    use_fixed_topology = any(
        isinstance(request.identity, AssociationParameterIdentity) for request in requests
    ) or (
        len(matches) == 1
        and all(
            any(
                slot.family is request.family and slot.identity == request.identity
                for slot in matches[0].slots
            )
            for request in requests
        )
    )
    if use_fixed_topology:
        if len(matches) != 1:
            raise ValueError(
                "installed EOS does not advertise exactly one derivative-ready "
                "fixed-topology association descriptor"
            )
        capability = matches[0]
        coordinates = []
        for request in requests:
            identities = (
                tuple(
                    AssociationParameterIdentity(request.family, (pair,))
                    for pair in request.identity.site_pairs
                )
                if isinstance(request.identity, AssociationParameterIdentity)
                else (request.identity,)
            )
            slots = tuple(
                slot
                for identity in identities
                for slot in capability.slots
                if slot.family is request.family and slot.identity == identity
            )
            if len(slots) != len(identities) or any(slot.unit != slots[0].unit for slot in slots):
                raise ValueError(
                    "association request identity does not match the installed "
                    "fixed-topology descriptor"
                )
            coordinates.append(
                ParameterCoordinate(
                    family=request.family,
                    identity=request.identity,
                    capability_id=None,
                    provider_parameter_fingerprint=capability.parameter_fingerprint,
                    provider_topology_fingerprint=capability.topology_fingerprint,
                    unit=slots[0].unit,
                    transform=request.transform,
                    lower_bound=request.lower_bound,
                    upper_bound=request.upper_bound,
                    provider_artifact_fingerprint=capability.artifact_fingerprint,
                )
            )
        return tuple(coordinates)
    coordinates = []
    observation_contract = _observation_contract(observations)
    for request in requests:
        matches = tuple(
            item
            for item in derivative_ready
            if item.family is request.family
            and item.observation_contract == observation_contract
            and _identity_matches(request, item)
        )
        if len(matches) != 1:
            raise ValueError(
                f"installed EOS capability for {request.family.value} is "
                f"{'ambiguous' if len(matches) > 1 else 'unavailable'}"
            )
        capability = matches[0]
        coordinates.append(
            ParameterCoordinate(
                request.family,
                request.identity,
                capability.capability_id,
                capability.parameter_fingerprint,
                capability.topology_fingerprint,
                _unit(capability.coordinate_units[-1]),
                request.transform,
                request.lower_bound,
                request.upper_bound,
            )
        )
    return tuple(coordinates)


_OBJECTIVE_FAMILIES = {
    FixedCompositionVleObservation: "fixed_composition_vle",
    PureSaturationObservation: "pure_saturation",
    PureVaporPressureObservation: "pure_vapor_pressure",
    PureDensityObservation: "pure_density",
    FixedCompositionMixtureDensityObservation: "fixed_composition_mixture_density",
    MeanIonicActivityObservation: "mean_ionic_activity",
    AqueousKijMeanIonicActivityObservation: "aqueous_kij_mean_ionic_activity",
    IonSolvationKijObservation: "ion_solvation_kij",
    SolvationGibbsObservation: "solvation_gibbs",
    RelativePermittivityRatioObservation: "relative_permittivity_ratio",
}


def _validate_objective(dataset: ObservationDataset) -> None:
    families = {_OBJECTIVE_FAMILIES[type(row)] for row in dataset.observations}
    allowed = (
        families
        if len(families) == 1
        else {"fixed_topology_association_mixed"}
        if families.issubset(
            {
                "pure_saturation",
                "pure_vapor_pressure",
                "pure_density",
            }
        )
        else set()
    )
    if dataset.objective.residual_family not in allowed:
        raise ValueError(
            "objective residual_family does not match the validated observation "
            f"contract; expected one of {sorted(allowed)}"
        )


@dataclass(frozen=True, slots=True)
class PreflightStart:
    start: tuple[float, ...]
    residual_count: int
    variable_count: int
    full_rank: int
    full_condition_number: float
    projected_parameter_rank: int
    projected_parameter_condition_number: float
    fitted_bounds_status: tuple[str, ...]
    lifted_bounds_status: tuple[str, ...]
    derivative_complete: bool
    failure_reason: str


@dataclass(frozen=True, slots=True)
class PreflightReport:
    fitted_parameter_count: int
    lifted_variable_count: int
    residual_count: int
    training_row_count: int
    held_out_row_count: int
    stress_row_count: int
    partition_counts: tuple[tuple[str, int], ...]
    starts: tuple[PreflightStart, ...]
    local_evidence_only: bool
    practical_identifiability_gate: str | None
    ready: bool
    reasons: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class ResidualIdentity:
    row_id: str
    component: str


@dataclass(frozen=True, slots=True)
class PreparedFitEvaluation:
    physical_parameter_point: tuple[float, ...]
    solver_parameter_point: tuple[float, ...]
    lifted_physical_point: tuple[float, ...]
    lifted_solver_point: tuple[float, ...]
    fitted_variable_identities: tuple[ParameterCoordinate, ...]
    lifted_variable_ids: tuple[str, ...]
    residual_identities: tuple[ResidualIdentity, ...]
    residual_vector: tuple[float, ...]
    jacobian: tuple[float, ...]
    jacobian_layout: str
    jacobian_diagnostics: GeneralJacobianDiagnostics
    parameter_fingerprints: tuple[str, ...]
    topology_fingerprints: tuple[str, ...]
    capability_artifact_fingerprints: tuple[str, ...]
    installed_eos_artifact_fingerprint: str
    preparation_fingerprint: str
    dataset_fingerprints: tuple[str, ...]
    source_fingerprints: tuple[str, ...]
    observation_order_fingerprint: str


def _shape(rows: tuple[RegressionObservation, ...]) -> tuple[int, int]:
    residuals = lifted = 0
    for row in rows:
        if isinstance(row, (FixedCompositionVleObservation, PureSaturationObservation)):
            residuals += 4
            lifted += 2
        elif isinstance(row, PureVaporPressureObservation):
            residuals += 3
            lifted += 2
        elif isinstance(
            row,
            (PureDensityObservation, FixedCompositionMixtureDensityObservation),
        ):
            residuals += 2
            lifted += 1
        else:
            residuals += 1
    return residuals, lifted


def _residual_identities(
    rows: tuple[RegressionObservation, ...],
) -> tuple[ResidualIdentity, ...]:
    identities = []
    for row in rows:
        components = (
            (
                "liquid_pressure",
                "vapor_pressure",
                "chemical_potential_0",
                "chemical_potential_1",
            )
            if isinstance(row, FixedCompositionVleObservation)
            else (
                "liquid_pressure",
                "vapor_pressure",
                "chemical_potential",
                "liquid_density",
            )
            if isinstance(row, PureSaturationObservation)
            else ("liquid_pressure", "vapor_pressure", "chemical_potential")
            if isinstance(row, PureVaporPressureObservation)
            else ("pressure", "density")
            if isinstance(
                row,
                (PureDensityObservation, FixedCompositionMixtureDensityObservation),
            )
            else (_OBJECTIVE_FAMILIES[type(row)],)
        )
        identities.extend(ResidualIdentity(row.row_id, component) for component in components)
    return tuple(identities)


def _lifted_variable_ids(
    rows: tuple[RegressionObservation, ...],
) -> tuple[str, ...]:
    identities = []
    for row in rows:
        if isinstance(
            row,
            (PureDensityObservation, FixedCompositionMixtureDensityObservation),
        ):
            identities.append(f"{row.row_id}:volume_m3_per_mol")
        elif isinstance(
            row,
            (
                FixedCompositionVleObservation,
                PureSaturationObservation,
                PureVaporPressureObservation,
            ),
        ):
            identities.extend(
                (
                    f"{row.row_id}:liquid_volume_m3_per_mol",
                    f"{row.row_id}:vapor_volume_m3_per_mol",
                )
            )
    return tuple(identities)


def _lifted_physical_point(
    rows: tuple[RegressionObservation, ...],
    solver_point: tuple[float, ...],
) -> tuple[float, ...]:
    values = []
    index = 0
    for row in rows:
        if isinstance(
            row,
            (PureDensityObservation, FixedCompositionMixtureDensityObservation),
        ):
            values.append(row.volume_origin_m3_per_mol * math.exp(solver_point[index]))
            index += 1
        elif isinstance(
            row,
            (
                FixedCompositionVleObservation,
                PureSaturationObservation,
                PureVaporPressureObservation,
            ),
        ):
            values.extend(
                (
                    row.liquid_volume_origin_m3_per_mol * math.exp(solver_point[index]),
                    row.vapor_volume_origin_m3_per_mol * math.exp(solver_point[index + 1]),
                )
            )
            index += 2
    if index != len(solver_point):
        raise ValueError("lifted solver point has the wrong dimension")
    return tuple(values)


def _lifted_start_variables(
    rows: tuple[RegressionObservation, ...],
) -> tuple[float, ...]:
    values: list[float] = []
    for row in rows:
        if isinstance(
            row,
            (PureDensityObservation, FixedCompositionMixtureDensityObservation),
        ):
            values.append(math.log(row.volume_start_m3_per_mol / row.volume_origin_m3_per_mol))
        elif isinstance(
            row,
            (
                FixedCompositionVleObservation,
                PureSaturationObservation,
                PureVaporPressureObservation,
            ),
        ):
            values.extend(
                (
                    math.log(
                        row.liquid_volume_start_m3_per_mol / row.liquid_volume_origin_m3_per_mol
                    ),
                    math.log(
                        row.vapor_volume_start_m3_per_mol / row.vapor_volume_origin_m3_per_mol
                    ),
                )
            )
    return tuple(values)


def _requires_resolved_lifted_start(problem: RegressionProblem) -> bool:
    return len(problem.parameters) > 1 and all(
        isinstance(
            row,
            (
                PureDensityObservation,
                FixedCompositionMixtureDensityObservation,
                PureSaturationObservation,
                PureVaporPressureObservation,
            ),
        )
        for row in problem.training_observations
    )


def _bound_status(value: float, bounds: tuple[float, float]) -> str:
    return "lower" if value == bounds[0] else "upper" if value == bounds[1] else "interior"


def _lifted_bounds_status(
    rows: tuple[RegressionObservation, ...],
) -> tuple[str, ...]:
    status: list[str] = []
    for row in rows:
        if isinstance(
            row,
            (PureDensityObservation, FixedCompositionMixtureDensityObservation),
        ):
            status.append(
                _bound_status(
                    row.volume_start_m3_per_mol,
                    row.volume_bounds_m3_per_mol,
                )
            )
        elif isinstance(
            row,
            (
                FixedCompositionVleObservation,
                PureSaturationObservation,
                PureVaporPressureObservation,
            ),
        ):
            status.extend(
                (
                    _bound_status(
                        row.liquid_volume_start_m3_per_mol,
                        row.liquid_volume_bounds_m3_per_mol,
                    ),
                    _bound_status(
                        row.vapor_volume_start_m3_per_mol,
                        row.vapor_volume_bounds_m3_per_mol,
                    ),
                )
            )
    return tuple(status)


def _lifted_bounds_status_at_variables(
    rows: tuple[RegressionObservation, ...],
    variables: tuple[float, ...],
) -> tuple[str, ...]:
    status: list[str] = []
    variable = 0
    for row in rows:
        if isinstance(
            row,
            (PureDensityObservation, FixedCompositionMixtureDensityObservation),
        ):
            volume = row.volume_origin_m3_per_mol * math.exp(variables[variable])
            status.append(_bound_status(volume, row.volume_bounds_m3_per_mol))
            variable += 1
        elif isinstance(
            row,
            (
                FixedCompositionVleObservation,
                PureSaturationObservation,
                PureVaporPressureObservation,
            ),
        ):
            liquid = row.liquid_volume_origin_m3_per_mol * math.exp(variables[variable])
            vapor = row.vapor_volume_origin_m3_per_mol * math.exp(variables[variable + 1])
            status.extend(
                (
                    _bound_status(liquid, row.liquid_volume_bounds_m3_per_mol),
                    _bound_status(vapor, row.vapor_volume_bounds_m3_per_mol),
                )
            )
            variable += 2
    if variable != len(variables):
        raise ValueError("resolved lifted start has the wrong dimension")
    return tuple(status)


def _evaluation_failure_reason(error: Exception) -> str:
    message = str(error)
    normalized = message.casefold()
    if "nonfinite" in normalized:
        return f"nonfinite_evaluation: {message}"
    if "derivative" in normalized or any(
        phrase in normalized
        for phrase in (
            "does not advertise",
            "does not expose",
            "capability unavailable",
            "capabilities do not expose",
        )
    ):
        return f"unavailable_derivatives: {message}"
    if any(
        phrase in normalized
        for phrase in (
            "violates bounds",
            "outside fitted bounds",
            "domain error",
            "not mechanically stable",
        )
    ):
        return f"eos_domain_failure: {message}"
    if any(
        token in normalized
        for token in (
            "component order",
            "capabilit",
            "descriptor",
            "fingerprint",
            "identity",
            "topology",
            "unit",
            "coordinate contract",
        )
    ):
        return f"incompatible_contract: {message}"
    if isinstance(error, TypeError) or any(
        token in normalized
        for token in (
            "field count",
            "payload",
            "must be a sequence",
            "wrong finite dimension",
            "inconsistent parameter",
        )
    ):
        return f"malformed_rows: {message}"
    return f"eos_domain_failure: {message}"


@dataclass(frozen=True, slots=True)
class PreparedFit:
    model: object
    problem: RegressionProblem
    datasets: tuple[ObservationDataset, ...]

    @property
    def preparation_fingerprint(self) -> str:
        return _canonical_sha256(
            {
                "problem": self.problem,
                "datasets": tuple(dataset.provenance_sha256 for dataset in self.datasets),
            }
        )

    def fit(self) -> RegressionResult:
        return replace(
            fit_parameters(self.problem, self.model),
            preparation_fingerprint=self.preparation_fingerprint,
        )

    def evaluate(
        self,
        physical_parameter_point: tuple[float, ...],
        *,
        lifted_solver_point: tuple[float, ...] | None = None,
    ) -> PreparedFitEvaluation:
        parameter_count = len(self.problem.parameters)
        residual_count, lifted_count = _shape(self.problem.training_observations)
        if (
            type(physical_parameter_point) is not tuple
            or len(physical_parameter_point) != parameter_count
        ):
            raise ValueError("physical parameter point has the wrong dimension")
        for value, coordinate in zip(
            physical_parameter_point, self.problem.parameters, strict=True
        ):
            _require_finite(value, "physical parameter point")
            if not coordinate.lower_bound <= value <= coordinate.upper_bound:
                raise ValueError("physical parameter point is outside fitted bounds")

        if lifted_solver_point is not None:
            if type(lifted_solver_point) is not tuple or len(lifted_solver_point) != lifted_count:
                raise ValueError("lifted solver point has the wrong dimension")
            for value in lifted_solver_point:
                _require_finite(value, "lifted solver point")

        solver_parameter_point = tuple(
            coordinate.transform.to_solver(value)
            for coordinate, value in zip(
                self.problem.parameters,
                physical_parameter_point,
                strict=True,
            )
        )
        try:
            if lifted_solver_point is None and _requires_resolved_lifted_start(self.problem):
                variables, residuals, jacobian = _evaluate_parameters_at_start(
                    self.problem,
                    self.model,
                    physical_parameter_point,
                )
                solver_parameter_point = variables[:parameter_count]
                lifted_solver_point = variables[parameter_count:]
            else:
                lifted_solver_point = (
                    _lifted_start_variables(self.problem.training_observations)
                    if lifted_solver_point is None
                    else lifted_solver_point
                )
                residuals, jacobian = _evaluate_parameters(
                    self.problem,
                    self.model,
                    (*solver_parameter_point, *lifted_solver_point),
                )
        except (ArithmeticError, RuntimeError, TypeError, ValueError) as exc:
            raise RuntimeError(_evaluation_failure_reason(exc)) from exc

        variable_count = parameter_count + lifted_count
        try:
            residual_vector = tuple(float(value) for value in residuals)
            jacobian_vector = tuple(float(value) for value in jacobian)
        except (TypeError, ValueError) as exc:
            raise RuntimeError(f"malformed_evaluation: {exc}") from exc
        if (
            len(residual_vector) != residual_count
            or len(jacobian_vector) != residual_count * variable_count
        ):
            raise RuntimeError(
                "unavailable_derivatives: incomplete exact residual or Jacobian payload"
            )
        if not all(
            math.isfinite(value)
            for value in (
                *solver_parameter_point,
                *lifted_solver_point,
                *residual_vector,
                *jacobian_vector,
            )
        ):
            raise RuntimeError("nonfinite_evaluation: exact evaluation payload")

        try:
            full, projected = _native.diagnose_jacobian(
                parameter_count,
                lifted_count,
                residual_count,
                jacobian_vector,
            )
            diagnostics = GeneralJacobianDiagnostics(
                residual_count=residual_count,
                variable_count=variable_count,
                full_singular_values=tuple(float(value) for value in full[0]),
                full_rank=int(full[1]),
                full_condition_number=float(full[2]),
                projected_parameter_singular_values=tuple(float(value) for value in projected[0]),
                projected_parameter_rank=int(projected[1]),
                projected_parameter_condition_number=float(projected[2]),
            )
        except (ArithmeticError, RuntimeError, TypeError, ValueError) as exc:
            raise RuntimeError(f"diagnostic_failure: {exc}") from exc
        if not all(
            math.isfinite(value)
            for value in (
                *diagnostics.full_singular_values,
                diagnostics.full_condition_number,
                *diagnostics.projected_parameter_singular_values,
                diagnostics.projected_parameter_condition_number,
            )
        ):
            raise RuntimeError("nonfinite_evaluation: Jacobian diagnostics")

        residual_identities = _residual_identities(self.problem.training_observations)
        try:
            lifted_physical_point = _lifted_physical_point(
                self.problem.training_observations,
                lifted_solver_point,
            )
        except (ArithmeticError, ValueError) as exc:
            raise RuntimeError(f"eos_domain_failure: {exc}") from exc
        try:
            installed_identity = dict(_installed_eos_artifact_identity())
        except (RuntimeError, TypeError, ValueError) as exc:
            raise RuntimeError(f"artifact_identity_failure: {exc}") from exc
        return PreparedFitEvaluation(
            physical_parameter_point=physical_parameter_point,
            solver_parameter_point=solver_parameter_point,
            lifted_physical_point=lifted_physical_point,
            lifted_solver_point=lifted_solver_point,
            fitted_variable_identities=self.problem.parameters,
            lifted_variable_ids=_lifted_variable_ids(self.problem.training_observations),
            residual_identities=residual_identities,
            residual_vector=residual_vector,
            jacobian=jacobian_vector,
            jacobian_layout="row_major",
            jacobian_diagnostics=diagnostics,
            parameter_fingerprints=tuple(
                dict.fromkeys(
                    item.provider_parameter_fingerprint for item in self.problem.parameters
                )
            ),
            topology_fingerprints=tuple(
                dict.fromkeys(
                    item.provider_topology_fingerprint for item in self.problem.parameters
                )
            ),
            capability_artifact_fingerprints=tuple(
                dict.fromkeys(
                    item.provider_artifact_fingerprint
                    for item in self.problem.parameters
                    if item.provider_artifact_fingerprint is not None
                )
            ),
            installed_eos_artifact_fingerprint=(f"sha256:{installed_identity['wheel_sha256']}"),
            preparation_fingerprint=self.preparation_fingerprint,
            dataset_fingerprints=tuple(dataset.provenance_sha256 for dataset in self.datasets),
            source_fingerprints=tuple(
                dataset.source.source_artifact_sha256 for dataset in self.datasets
            ),
            observation_order_fingerprint=(f"sha256:{_canonical_sha256(residual_identities)}"),
        )

    def preflight(self) -> PreflightReport:
        parameter_count = len(self.problem.parameters)
        residual_count, lifted_count = _shape(self.problem.training_observations)
        lifted_start = _lifted_start_variables(self.problem.training_observations)
        if len(lifted_start) != lifted_count:
            raise ValueError("preflight lifted start accounting is incomplete")
        lifted_bounds = _lifted_bounds_status(self.problem.training_observations)
        if len(lifted_bounds) != lifted_count:
            raise ValueError("preflight lifted bound accounting is incomplete")
        starts = []
        reasons: list[str] = []
        variable_count = parameter_count + lifted_count
        for physical_start, solver_start in zip(
            self.problem.start_vectors, self.problem.solver_start_vectors, strict=True
        ):
            bounds = tuple(
                _bound_status(
                    value,
                    (coordinate.lower_bound, coordinate.upper_bound),
                )
                for value, coordinate in zip(physical_start, self.problem.parameters, strict=True)
            )
            if residual_count < variable_count:
                reason = (
                    "structural_insufficiency: residual count is smaller than "
                    "the fitted-plus-lifted variable count"
                )
                reasons.append(reason)
                starts.append(
                    PreflightStart(
                        physical_start,
                        residual_count,
                        variable_count,
                        0,
                        math.inf,
                        0,
                        math.inf,
                        bounds,
                        lifted_bounds,
                        False,
                        reason,
                    )
                )
                continue
            try:
                if _requires_resolved_lifted_start(self.problem):
                    variables, residuals, jacobian = _evaluate_parameters_at_start(
                        self.problem,
                        self.model,
                        physical_start,
                    )
                    start_lifted_bounds = _lifted_bounds_status_at_variables(
                        self.problem.training_observations,
                        variables[parameter_count:],
                    )
                else:
                    residuals, jacobian = _evaluate_parameters(
                        self.problem,
                        self.model,
                        (*solver_start, *lifted_start),
                    )
                    start_lifted_bounds = lifted_bounds
                if (
                    len(residuals) != residual_count
                    or len(jacobian) != residual_count * variable_count
                ):
                    raise RuntimeError(
                        "unavailable_derivatives: incomplete exact derivative columns"
                    )
                if not all(math.isfinite(value) for value in (*residuals, *jacobian)):
                    raise ArithmeticError("nonfinite exact residual or derivative payload")
                full, projected = _native.diagnose_jacobian(
                    parameter_count,
                    lifted_count,
                    residual_count,
                    jacobian,
                )
                full_rank, full_condition = int(full[1]), float(full[2])
                projected_rank, projected_condition = (
                    int(projected[1]),
                    float(projected[2]),
                )
                failure = ""
                if full_rank < variable_count:
                    failure = "rank_deficient_full_jacobian"
                elif projected_rank < parameter_count:
                    failure = "rank_deficient_projected_parameter_jacobian"
                elif (
                    full_condition > self.problem.maximum_condition_number
                    or projected_condition > self.problem.maximum_condition_number
                ):
                    failure = "poor_conditioning"
                if failure:
                    reasons.append(failure)
                starts.append(
                    PreflightStart(
                        physical_start,
                        residual_count,
                        variable_count,
                        full_rank,
                        full_condition,
                        projected_rank,
                        projected_condition,
                        bounds,
                        start_lifted_bounds,
                        True,
                        failure,
                    )
                )
            except (ArithmeticError, RuntimeError, TypeError, ValueError) as exc:
                reason = (
                    str(exc)
                    if str(exc).startswith("unavailable_derivatives:")
                    else _evaluation_failure_reason(exc)
                )
                reasons.append(reason)
                starts.append(
                    PreflightStart(
                        physical_start,
                        residual_count,
                        parameter_count + lifted_count,
                        0,
                        math.inf,
                        0,
                        math.inf,
                        bounds,
                        lifted_bounds,
                        False,
                        reason,
                    )
                )
        gate = (
            "topology-specific identifiability and accepted-region evidence required"
            if _uses_fixed_topology_descriptor(self.problem.parameters)
            else None
        )
        unique_reasons = tuple(dict.fromkeys(reasons))
        return PreflightReport(
            parameter_count,
            lifted_count,
            residual_count,
            len(self.problem.training_observations),
            len(self.problem.held_out_observations),
            len(self.problem.stress_observations),
            (
                ("training", len(self.problem.training_observations)),
                ("held_out", len(self.problem.held_out_observations)),
                ("stress", len(self.problem.stress_observations)),
            ),
            tuple(starts),
            True,
            gate,
            not unique_reasons,
            unique_reasons,
        )


def prepare_fit(
    model: object,
    *,
    datasets: tuple[ObservationDataset, ...],
    parameters: tuple[ParameterRequest, ...],
    parameter_slot_indices: tuple[int, ...],
    start_vectors: tuple[tuple[float, ...], ...],
    solver: SolverControls,
    rank: RankControls,
    confirmation: ConfirmationControls,
) -> PreparedFit:
    if not datasets:
        raise ValueError("datasets must be nonempty")
    for dataset in datasets:
        _validate_objective(dataset)
    observations = tuple(row for dataset in datasets for row in dataset.observations)
    coordinates = _resolve_coordinates(model, parameters, observations)
    observation_types = {type(row) for row in observations}
    if len(observation_types) > 1 and not _uses_fixed_topology_descriptor(coordinates):
        raise ValueError(
            "mixed observation contracts are supported only by the installed "
            "fixed-topology association descriptor"
        )
    sources = tuple(
        SourceDescriptor(
            source_id=dataset.source.source_id,
            citation=dataset.source.citation,
            durable_locator=dataset.source.durable_locator,
            source_artifact_sha256=dataset.source.source_artifact_sha256,
            canonical_dataset_sha256=canonical_dataset_sha256(dataset.observations),
            transformation_record=dataset.source.transformation_record,
            units_and_bases=dataset.source.units_and_bases,
            use_basis=dataset.source.use_basis,
            residual_scale_rationale=dataset.source.residual_scale_rationale,
        )
        for dataset in datasets
    )
    problem = RegressionProblem(
        sources=sources,
        parameters=coordinates,
        parameter_slot_indices=parameter_slot_indices,
        start_vectors=start_vectors,
        observations=observations,
        maximum_condition_number=rank.maximum_condition_number,
        maximum_iterations=solver.maximum_iterations,
        maximum_solver_time_seconds=solver.maximum_solver_time_seconds,
        function_tolerance=solver.function_tolerance,
        gradient_tolerance=solver.gradient_tolerance,
        parameter_tolerance=solver.parameter_tolerance,
        confirmation_parameter_scaled_max_delta=confirmation.parameter_scaled_max_delta,
        confirmation_cost_relative_delta=confirmation.cost_relative_delta,
    )
    return PreparedFit(model, problem, datasets)


__all__ = (
    "AcquisitionClass",
    "ConfirmationControls",
    "CorrelationProvenance",
    "ObjectiveContract",
    "ObservationDataset",
    "ParameterRequest",
    "PreflightReport",
    "PreflightStart",
    "PreparedFit",
    "PreparedFitEvaluation",
    "RankControls",
    "ResidualIdentity",
    "RowProvenance",
    "SolverControls",
    "SourceInput",
    "prepare_fit",
)
