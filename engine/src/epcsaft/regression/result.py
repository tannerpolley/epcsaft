from __future__ import annotations

import hashlib
import json
import math
from collections.abc import Mapping
from dataclasses import dataclass, replace
from enum import StrEnum
from functools import cache
from importlib.metadata import distribution
from pathlib import Path
from urllib.parse import unquote, urlparse

from epcsaft import Parameters, SourceRecord, ValidityDomain
from epcsaft import regression as _package
from epcsaft.parameters import _ParameterBundle
from epcsaft.records import (
    AssociationParameterRecord,
    ConstantCorrelation,
    ConstantPlusSumOfExponentialsCorrelation,
    ModelParameterRecord,
    PairParameterRecord,
    SingleParameterRecord,
    SiteRecord,
)

from .evaluator_regression import (
    ComposedPositiveRowDiagnostic,
    PositiveEvaluatorCapability,
    PositiveEvaluatorProblem,
)
from .parameter_regression import (
    AssociationParameterIdentity,
    ComponentParameterIdentity,
    DensityRowDiagnostic,
    DirectObservationRowDiagnostic,
    FittedParameterDiagnostic,
    FixedTopologyAssociationCapability,
    GeneralJacobianDiagnostics,
    GeneralRowDiagnostic,
    ParameterCapability,
    PairParameterIdentity,
    PureSaturationRowDiagnostic,
    PureVaporPressureRowDiagnostic,
    RegressionProblem,
    ModelParameterIdentity,
    _require_nonempty_string,
    _require_sha256,
    canonical_dataset_sha256,
)


@dataclass(frozen=True, slots=True)
class RegressionResult:
    problem: RegressionProblem | PositiveEvaluatorProblem
    capabilities: tuple[
        ParameterCapability | FixedTopologyAssociationCapability | PositiveEvaluatorCapability,
        ...,
    ]
    preparation_fingerprint: str | None
    provider_parameter_fingerprint: str
    provider_topology_fingerprint: str
    solver_converged: bool
    numerically_converged: bool
    workflow_valid: bool
    physical_status: str
    scientific_status: str
    predictive_status: str
    authority_status: str
    termination: str
    solution_usable: bool
    initial_cost: float
    final_cost: float
    iterations: int
    residual_evaluation_count: int
    jacobian_evaluation_count: int
    parameters: tuple[FittedParameterDiagnostic, ...]
    jacobian: GeneralJacobianDiagnostics
    rows: tuple[
        GeneralRowDiagnostic
        | PureSaturationRowDiagnostic
        | PureVaporPressureRowDiagnostic
        | DensityRowDiagnostic
        | DirectObservationRowDiagnostic
        | ComposedPositiveRowDiagnostic,
        ...,
    ]
    confirmation_count: int
    confirmation_parameter_scaled_max_delta: float
    confirmation_cost_relative_max_delta: float
    confirmations_usable: bool
    training_row_count: int
    held_out_row_count: int
    stress_row_count: int
    evaluated_row_count: int
    skipped_row_count: int
    failed_row_count: int
    failure_reasons: tuple[str, ...]

    def to_record(
        self,
        *,
        prepared: _package.PreparedFit | None = None,
        context: ResultContext | None = None,
    ) -> dict[str, object]:
        return _result_record(self, prepared=prepared, context=context)

    def to_json_bytes(
        self,
        *,
        prepared: _package.PreparedFit | None = None,
        context: ResultContext | None = None,
    ) -> bytes:
        return _result_json_bytes(self, prepared=prepared, context=context)

    def to_parameters(
        self,
        parameters: Parameters,
        *,
        bundle_id: str,
        bundle_version: int,
        source: SourceRecord,
        domain: ValidityDomain,
        locator: str,
    ) -> Parameters:
        """Apply a valid fit to its exact parameter subject with explicit provenance."""
        return _parameters_from_result(
            self,
            parameters,
            bundle_id=bundle_id,
            bundle_version=bundle_version,
            source=source,
            domain=domain,
            locator=locator,
        )


def _parameters_from_result(
    result: RegressionResult,
    parameters: Parameters,
    *,
    bundle_id: str,
    bundle_version: int,
    source: SourceRecord,
    domain: ValidityDomain,
    locator: str,
) -> Parameters:
    if not result.workflow_valid or not result.solution_usable:
        raise ValueError("only a valid usable regression result can update parameters")
    if result.provider_parameter_fingerprint != parameters.fingerprint:
        raise ValueError("regression result does not match the supplied parameter fingerprint")
    if not isinstance(source, SourceRecord) or not isinstance(domain, ValidityDomain):
        raise TypeError("source and domain must be typed parameter provenance records")
    if type(locator) is not str or not locator.strip():
        raise ValueError("locator must be a nonempty string")
    if not isinstance(result.problem, RegressionProblem):
        raise TypeError("only typed parameter-regression results can update parameters")
    if len(result.problem.parameters) != len(result.parameters):
        raise ValueError("regression result has incomplete fitted-parameter diagnostics")

    records = list(parameters._bundle.records)
    for coordinate, diagnostic in zip(
        result.problem.parameters,
        result.parameters,
        strict=True,
    ):
        if (
            diagnostic.family is not coordinate.family
            or diagnostic.component_ids != coordinate.identity.canonical_component_ids
        ):
            raise ValueError("fitted-parameter diagnostics do not match their coordinates")

        matches: list[int] = []
        identity = coordinate.identity
        for index, record in enumerate(records):
            if record.family != coordinate.family.value:
                continue
            if isinstance(identity, ComponentParameterIdentity) and isinstance(
                record, SingleParameterRecord
            ):
                matched = record.component_id == identity.component_id
            elif isinstance(identity, PairParameterIdentity) and isinstance(
                record, PairParameterRecord
            ):
                matched = (
                    record.component_id_a,
                    record.component_id_b,
                ) == identity.component_ids
            elif isinstance(identity, ModelParameterIdentity) and isinstance(
                record, ModelParameterRecord
            ):
                matched = True
            elif isinstance(identity, AssociationParameterIdentity) and isinstance(
                record, AssociationParameterRecord
            ):
                matched = any(
                    (
                        record.component_id_a,
                        record.site_id_a,
                        record.component_id_b,
                        record.site_id_b,
                    )
                    == (
                        pair.component_id_a,
                        pair.site_id_a,
                        pair.component_id_b,
                        pair.site_id_b,
                    )
                    for pair in identity.site_pairs
                )
            else:
                matched = False
            if matched:
                matches.append(index)

        expected_matches = (
            len(identity.site_pairs) if isinstance(identity, AssociationParameterIdentity) else 1
        )
        if len(matches) != expected_matches:
            raise ValueError("fitted parameter does not map to exactly one typed record identity")
        for index in matches:
            records[index] = replace(
                records[index],
                value=diagnostic.final,
                source_id=source.source_id,
                locator=locator,
                domain_id=domain.domain_id,
            )

    sources = {item.source_id: item for item in parameters._bundle.sources}
    if source.source_id in sources and sources[source.source_id] != source:
        raise ValueError("source_id conflicts with existing parameter provenance")
    sources[source.source_id] = source
    domains = {item.domain_id: item for item in parameters._bundle.domains}
    if domain.domain_id in domains and domains[domain.domain_id] != domain:
        raise ValueError("domain_id conflicts with existing parameter provenance")
    domains[domain.domain_id] = domain

    bundle = _ParameterBundle.from_records(
        bundle_id=bundle_id,
        bundle_version=bundle_version,
        purpose="regression-output",
        sources=sources.values(),
        domains=domains.values(),
        components=parameters._bundle.components,
        singles=(item for item in records if isinstance(item, SingleParameterRecord)),
        pairs=(item for item in records if isinstance(item, PairParameterRecord)),
        sites=(item for item in records if isinstance(item, SiteRecord)),
        associations=(item for item in records if isinstance(item, AssociationParameterRecord)),
        correlations=(
            item
            for item in records
            if isinstance(item, (ConstantCorrelation, ConstantPlusSumOfExponentialsCorrelation))
        ),
        models=(item for item in records if isinstance(item, ModelParameterRecord)),
    )
    return bundle.select(parameters.component_ids)


class ReproductionClass(StrEnum):
    EXACT_AUTHOR_METHOD_REPLAY = "EXACT_AUTHOR_METHOD_REPLAY"
    SOURCE_FAITHFUL_RECONSTRUCTION = "SOURCE_FAITHFUL_RECONSTRUCTION"
    PUBLISHED_TUPLE_PROPERTY_REPLAY = "PUBLISHED_TUPLE_PROPERTY_REPLAY"
    MODERN_REFIT = "MODERN_REFIT"


@dataclass(frozen=True, slots=True)
class LiteratureModelIdentity:
    formulation: str
    phase_reference_convention: str
    association_scheme: str
    site_multiplicities: str
    mixing_rule: str
    combining_rule: str
    mixture_parameter_treatment: str

    def __post_init__(self) -> None:
        for name in self.__dataclass_fields__:
            _require_nonempty_string(getattr(self, name), name)


@dataclass(frozen=True, slots=True)
class ArtifactReference:
    artifact_id: str
    sha256: str

    def __post_init__(self) -> None:
        _require_nonempty_string(self.artifact_id, "artifact_id")
        _require_sha256(self.sha256, "artifact SHA-256")


@dataclass(frozen=True, slots=True)
class ResultContext:
    reproduction_class: ReproductionClass | None = None
    model_identity: LiteratureModelIdentity | None = None
    reconstruction_decisions: tuple[str, ...] = ()
    source_printed_parameters: tuple[tuple[str, str], ...] = ()
    profile_artifact: ArtifactReference | None = None
    bootstrap_artifact: ArtifactReference | None = None
    uncertainty_artifact: ArtifactReference | None = None
    validation_campaign_artifact: ArtifactReference | None = None

    def __post_init__(self) -> None:
        if self.reproduction_class is not None and not isinstance(
            self.reproduction_class, ReproductionClass
        ):
            raise TypeError("reproduction_class must be a ReproductionClass")
        if self.reproduction_class is not None and self.model_identity is None:
            raise ValueError("a literature reproduction class requires complete model_identity")
        if self.model_identity is not None and self.reproduction_class is None:
            raise ValueError("literature model_identity requires exactly one reproduction class")
        if type(self.reconstruction_decisions) is not tuple or any(
            type(decision) is not str or not decision.strip()
            for decision in self.reconstruction_decisions
        ):
            raise ValueError("reconstruction_decisions must be a tuple of nonempty strings")
        reconstructed = self.reproduction_class in (
            ReproductionClass.SOURCE_FAITHFUL_RECONSTRUCTION,
            ReproductionClass.MODERN_REFIT,
        )
        if reconstructed != bool(self.reconstruction_decisions):
            raise ValueError(
                "source-faithful reconstructions and modern refits require "
                "explicit reconstruction decisions; exact and tuple replays "
                "forbid them"
            )
        if self.model_identity is not None and not isinstance(
            self.model_identity, LiteratureModelIdentity
        ):
            raise TypeError("model_identity must be a LiteratureModelIdentity")
        if any(
            not isinstance(value, ArtifactReference)
            for value in (
                self.profile_artifact,
                self.bootstrap_artifact,
                self.uncertainty_artifact,
                self.validation_campaign_artifact,
            )
            if value is not None
        ):
            raise TypeError("artifact references must be ArtifactReference values")
        if type(self.source_printed_parameters) is not tuple or any(
            type(item) is not tuple or len(item) != 2 for item in self.source_printed_parameters
        ):
            raise TypeError("source_printed_parameters must be a tuple of name/value pairs")
        if any(
            type(name) is not str or not name.strip() or type(value) is not str or not value.strip()
            for name, value in self.source_printed_parameters
        ) or len({name for name, _ in self.source_printed_parameters}) != len(
            self.source_printed_parameters
        ):
            raise ValueError("source_printed_parameters must contain unique nonempty pairs")
        if self.source_printed_parameters and self.reproduction_class is None:
            raise ValueError("source-printed parameters require a literature reproduction class")


@cache
def _installed_eos_artifact_identity() -> tuple[tuple[str, str], ...]:
    installed = distribution("epcsaft")
    files = tuple(installed.files or ())

    def digest(suffix: str) -> str:
        matches = tuple(file for file in files if str(file).endswith(suffix))
        if len(matches) != 1:
            raise RuntimeError(f"installed epcsaft artifact must contain exactly one {suffix}")
        return hashlib.sha256(installed.locate_file(matches[0]).read_bytes()).hexdigest()

    direct_url_matches = tuple(
        file for file in files if str(file).endswith(".dist-info/direct_url.json")
    )
    if len(direct_url_matches) != 1:
        raise RuntimeError("installed epcsaft artifact must contain exactly one direct_url.json")
    direct_url = json.loads(
        installed.locate_file(direct_url_matches[0]).read_text(encoding="utf-8")
    )
    parsed = urlparse(direct_url["url"])
    if parsed.scheme != "file":
        raise RuntimeError("installed epcsaft wheel must have a local immutable artifact URL")
    wheel = Path(unquote(parsed.path))
    if not wheel.is_file():
        raise RuntimeError("installed epcsaft wheel artifact is unavailable")
    wheel_sha256 = hashlib.sha256(wheel.read_bytes()).hexdigest()
    header_sha256 = digest("epcsaft/include/epcsaft/native_sdk_v1.h")

    return (
        ("distribution", f"{installed.metadata['Name']}=={installed.version}"),
        ("wheel_sha256", wheel_sha256),
        ("record_sha256", digest(".dist-info/RECORD")),
        ("public_header_sha256", header_sha256),
    )


def _ensure_finite(value: object, path: str = "record") -> None:
    if isinstance(value, float) and not math.isfinite(value):
        raise ValueError(f"{path} must contain only finite numbers")
    if isinstance(value, Mapping):
        for key, item in value.items():
            _ensure_finite(item, f"{path}.{key}")
    elif isinstance(value, (tuple, list)):
        for index, item in enumerate(value):
            _ensure_finite(item, f"{path}[{index}]")


def _without_unavailable(value: object, path: str = "") -> tuple[object | None, tuple[str, ...]]:
    if isinstance(value, float) and not math.isfinite(value):
        return None, (path,)
    if isinstance(value, Mapping):
        available: dict[object, object] = {}
        unavailable: list[str] = []
        for key, item in value.items():
            item_path = f"{path}.{key}" if path else str(key)
            retained, missing = _without_unavailable(item, item_path)
            unavailable.extend(missing)
            if retained is not None:
                available[key] = retained
        return available, tuple(unavailable)
    if isinstance(value, list):
        retained_items = []
        unavailable = []
        for index, item in enumerate(value):
            item_path = f"{path}[{index}]"
            retained, missing = _without_unavailable(item, item_path)
            unavailable.extend(missing)
            if retained is not None:
                retained_items.append(retained)
        if unavailable:
            return None, tuple(unavailable)
        return retained_items, tuple(unavailable)
    return value, ()


def _jacobian_record(result: RegressionResult, json_value) -> dict[str, object]:
    record = json_value(result.jacobian)
    assert isinstance(record, dict)
    retained, unavailable = _without_unavailable(record)
    assert isinstance(retained, dict)
    retained["unavailable_fields"] = list(unavailable)
    return retained


def _row_records(result: RegressionResult, json_value) -> list[object]:
    records = []
    for row in result.rows:
        record = json_value(row)
        if getattr(row, "evaluated", True) is not False:
            records.append(record)
            continue
        retained, unavailable = _without_unavailable(record)
        assert isinstance(retained, dict)
        retained["unavailable_fields"] = list(unavailable)
        records.append(retained)
    return records


def _capability_records(result: RegressionResult, json_value) -> list[object]:
    records = json_value(result.capabilities)
    assert isinstance(records, list)
    for capability, record in zip(result.capabilities, records, strict=True):
        if not isinstance(capability, FixedTopologyAssociationCapability):
            continue
        assert isinstance(record, dict)
        slots = record["slots"]
        assert isinstance(slots, list)
        for slot, slot_record in zip(capability.slots, slots, strict=True):
            assert isinstance(slot_record, dict)
            if not math.isfinite(slot.lower_bound_exclusive):
                raise ValueError("fixed-topology lower_bound_exclusive must be finite")
            if math.isfinite(slot.upper_bound_exclusive):
                continue
            if slot.upper_bound_exclusive != math.inf:
                raise ValueError("fixed-topology upper_bound_exclusive must be finite or open")
            del slot_record["upper_bound_exclusive"]
            slot_record["unavailable_fields"] = ["upper_bound_exclusive"]
            slot_record["open_bounds"] = [
                {
                    "direction": "upper",
                    "field": "upper_bound_exclusive",
                }
            ]
    return records


def _resolved_data_identity(result: RegressionResult, json_value) -> object:
    return {
        source.source_id: {
            "source": json_value(source),
            "rows": [
                json_value(row)
                for row in result.problem.observations
                if row.source_id == source.source_id
            ],
            "row_provenance": {"status": "not_available_in_low_level_problem"},
        }
        for source in result.problem.sources
    }


def _resolved_objective_identity(result: RegressionResult, json_value) -> object:
    from .usability import _OBJECTIVE_FAMILIES, ObjectiveContract

    if isinstance(result.problem, PositiveEvaluatorProblem):
        residual_family = "positive_scalar"
    else:
        families = {_OBJECTIVE_FAMILIES[type(row)] for row in result.problem.observations}
        residual_family = (
            next(iter(families)) if len(families) == 1 else "fixed_topology_association_mixed"
        )
    return {
        "resolution": "canonical_problem_contract",
        "contract": json_value(
            ObjectiveContract(
                residual_family=residual_family,
                interpretation="native_scaled_least_squares",
                row_weighting="observation_residual_scales",
                covariance_interpretation="independent_no_covariance",
                loss="squared",
                loss_parameters=(),
                failed_row_treatment="fail_fit",
            )
        ),
    }


def _result_record(
    result: RegressionResult,
    *,
    prepared: _package.PreparedFit | None = None,
    context: ResultContext | None = None,
) -> dict[str, object]:
    from .usability import _json_value

    context = context or ResultContext()
    if context.reproduction_class is not None and prepared is None:
        raise ValueError(
            "literature result export requires PreparedFit data and objective provenance"
        )
    if context.reproduction_class is ReproductionClass.PUBLISHED_TUPLE_PROPERTY_REPLAY:
        raise ValueError(
            "PUBLISHED_TUPLE_PROPERTY_REPLAY requires a no-fit property result, "
            "not RegressionResult"
        )
    if prepared is None and result.preparation_fingerprint is not None:
        raise ValueError("a prepared RegressionResult requires its exact PreparedFit provenance")
    if prepared is not None and (
        prepared.problem != result.problem
        or result.preparation_fingerprint != prepared.preparation_fingerprint
    ):
        raise ValueError("prepared fit provenance does not own this RegressionResult")
    positive = isinstance(result.problem, PositiveEvaluatorProblem)
    problem = {
        "kind": type(result.problem).__name__,
        "resolved": _json_value(result.problem),
        "datasets": (
            [
                {
                    "source": _json_value(dataset.source),
                    "canonical_dataset_sha256": canonical_dataset_sha256(dataset.observations),
                    "provenance_sha256": dataset.provenance_sha256,
                    "objective": _json_value(dataset.objective),
                    "rows": [
                        {
                            "row_id": row_id,
                            **_json_value(provenance),
                        }
                        for row_id, provenance in dataset.row_provenance
                    ],
                }
                for dataset in prepared.datasets
            ]
            if prepared is not None
            else None
        ),
    }
    record = {
        "schema_id": "epcsaft.regression.result.v1",
        "schema_version": 1,
        "problem": problem,
        "preparation_fingerprint": result.preparation_fingerprint,
        "data_identity": (
            _resolved_data_identity(result, _json_value)
            if prepared is None
            else _json_value(
                {
                    dataset.source.source_id: {
                        row_id: provenance for row_id, provenance in dataset.row_provenance
                    }
                    for dataset in prepared.datasets
                }
            )
        ),
        "objective_identity": (
            _resolved_objective_identity(result, _json_value)
            if prepared is None
            else _json_value(
                {dataset.source.source_id: dataset.objective for dataset in prepared.datasets}
            )
        ),
        "capabilities": _capability_records(result, _json_value),
        "installed_artifacts": (
            {
                "provider": result.capabilities[0].provider_artifact_identity,
                "owner": result.capabilities[0].owner_artifact_identity,
                "evaluator": result.capabilities[0].artifact_identity,
            }
            if positive
            else dict(_installed_eos_artifact_identity())
        ),
        "provider_parameter_fingerprint": result.provider_parameter_fingerprint,
        "provider_topology_fingerprint": result.provider_topology_fingerprint,
        "status": {
            "solver_converged": result.solver_converged,
            "numerically_converged": result.numerically_converged,
            "workflow_valid": result.workflow_valid,
            "physical_status": result.physical_status,
            "scientific_status": result.scientific_status,
            "predictive_status": result.predictive_status,
            "authority_status": result.authority_status,
            "termination": result.termination,
            "solution_usable": result.solution_usable,
            "failure_reasons": list(result.failure_reasons),
        },
        "solver": {
            "initial_cost": result.initial_cost,
            "final_cost": result.final_cost,
            "iterations": result.iterations,
            "residual_evaluation_count": result.residual_evaluation_count,
            "jacobian_evaluation_count": result.jacobian_evaluation_count,
        },
        "parameters": _json_value(result.parameters),
        "jacobian": _jacobian_record(result, _json_value),
        "rows": _row_records(result, _json_value),
        "confirmation": {
            "count": result.confirmation_count,
            "parameter_scaled_max_delta": (result.confirmation_parameter_scaled_max_delta),
            "cost_relative_max_delta": result.confirmation_cost_relative_max_delta,
            "usable": result.confirmations_usable,
        },
        "row_accounting": {
            "training": result.training_row_count,
            "held_out": result.held_out_row_count,
            "stress": result.stress_row_count,
            "evaluated": result.evaluated_row_count,
            "skipped": result.skipped_row_count,
            "failed": result.failed_row_count,
        },
        "literature": {
            "reproduction_class": (
                context.reproduction_class.value if context.reproduction_class is not None else None
            ),
            "model_identity": _json_value(context.model_identity),
            "reconstruction_decisions": list(context.reconstruction_decisions),
            "source_printed_parameters": dict(context.source_printed_parameters),
            "profile_artifact": _json_value(context.profile_artifact),
            "bootstrap_artifact": _json_value(context.bootstrap_artifact),
            "uncertainty_artifact": _json_value(context.uncertainty_artifact),
            "validation_campaign_artifact": _json_value(context.validation_campaign_artifact),
        },
    }
    _ensure_finite(record)
    return record


def _result_json_bytes(
    result: RegressionResult,
    *,
    prepared: _package.PreparedFit | None = None,
    context: ResultContext | None = None,
) -> bytes:
    return json.dumps(
        _result_record(result, prepared=prepared, context=context),
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode("utf-8")


__all__ = (
    "ArtifactReference",
    "LiteratureModelIdentity",
    "RegressionResult",
    "ReproductionClass",
    "ResultContext",
)
