from __future__ import annotations

import csv
import hashlib
import json
import math
from collections.abc import Sequence
from dataclasses import dataclass
from os import PathLike
from pathlib import Path
from typing import SupportsFloat, cast

import pint

from ..eos import EosError, Mixture, ParameterCapability, Pure, SaturationError
from ..parameters import ActiveParameterSet, Parameters
from ..units import unit_registry
from ._model import (
    SCHEMA_REGISTRY,
    ActiveSelection,
    ObservationDefinition,
    ObservationSpec,
    ParameterFit,
    PublishedObservationCorrelation,
    QuadratureSpec,
    SourceDomain,
    TableRow,
    TableSource,
)


class FitResolutionError(ValueError):
    def __init__(self, code: str, path: str, detail: object) -> None:
        self.code = code
        self.path = path
        self.detail = detail
        super().__init__(f"{code} at {path}: {detail}")


@dataclass(frozen=True, slots=True)
class SourceCheck:
    path: str
    expected_sha256: str
    actual_sha256: str


@dataclass(frozen=True, slots=True)
class ResolvedCoordinate:
    identity: str
    tied_identities: tuple[str, ...]
    bounds: tuple[float, float]
    scale: float
    unit: str


@dataclass(frozen=True, slots=True)
class QuadratureRule:
    order: int
    nodes: tuple[float, ...]
    weights: tuple[float, ...]
    observed_values: tuple[float, ...]
    observed_derivatives: tuple[float, ...]


@dataclass(frozen=True, slots=True)
class QuadraturePlan:
    correlation: str
    domain: SourceDomain
    measure: str
    controls: QuadratureSpec
    rules: tuple[QuadratureRule, ...]
    accepted_order: int
    normalized_objective: float
    objective_jacobian: tuple[float, ...]
    objective_change: float
    jacobian_change: float


@dataclass(frozen=True, slots=True)
class ResolvedObservation:
    specification: ObservationSpec
    capability: ParameterCapability
    canonical_units: tuple[tuple[str, str], ...]
    rows: tuple[TableRow, ...] | None
    source_check: SourceCheck
    quadrature: QuadraturePlan | None


@dataclass(frozen=True, slots=True)
class ResolvedFit:
    fit: ParameterFit
    parameters: Parameters
    active_parameters: ActiveParameterSet
    selected_identities: tuple[str, ...]
    independent_identities: tuple[str, ...]
    coordinates: tuple[ResolvedCoordinate, ...]
    observations: tuple[ResolvedObservation, ...]
    parameters_source_check: SourceCheck
    digest: str

    @property
    def source_checks(self) -> tuple[SourceCheck, ...]:
        return (
            self.parameters_source_check,
            *(observation.source_check for observation in self.observations),
        )

    def to_mapping(self) -> dict[str, object]:
        fit_observations = self.fit.to_mapping()["observations"]
        assert isinstance(fit_observations, list)
        return {
            "fit": self.fit.to_mapping(),
            "parameters_fingerprint": self.parameters.fingerprint,
            "parameters_source_check": _source_check_mapping(self.parameters_source_check),
            "active_parameter_fingerprint": self.active_parameters.fingerprint,
            "selected_identities": list(self.selected_identities),
            "independent_identities": list(self.independent_identities),
            "coordinates": [
                {
                    "identity": coordinate.identity,
                    "tied_identities": list(coordinate.tied_identities),
                    "bounds": list(coordinate.bounds),
                    "scale": coordinate.scale,
                    "unit": coordinate.unit,
                }
                for coordinate in self.coordinates
            ],
            "observations": [
                {
                    "specification": fit_observations[index],
                    "capability": _capability_mapping(observation.capability),
                    "canonical_units": dict(observation.canonical_units),
                    "rows": None
                    if observation.rows is None
                    else [{name: value for name, value in row.values} for row in observation.rows],
                    "source_check": _source_check_mapping(observation.source_check),
                    "quadrature": _quadrature_mapping(observation.quadrature),
                }
                for index, observation in enumerate(self.observations)
            ],
        }


def _canonical(value: object) -> bytes:
    return json.dumps(value, ensure_ascii=False, separators=(",", ":"), sort_keys=True).encode(
        "utf-8"
    )


def _hash(value: object) -> str:
    return f"sha256:{hashlib.sha256(_canonical(value)).hexdigest()}"


def _path_root(base_path: str | PathLike[str] | None) -> Path:
    if base_path is None:
        return Path.cwd()
    candidate = Path(base_path)
    if candidate.is_dir():
        return candidate
    return candidate.parent


def _local_path(path: str, root: Path) -> Path:
    target = Path(path)
    return target if target.is_absolute() else root / target


def _file_check(
    path: str,
    expected: str,
    *,
    root: Path,
    inline: object | None = None,
) -> SourceCheck:
    if inline is not None:
        actual = _hash(inline)
    else:
        target = _local_path(path, root)
        if not target.is_file():
            raise FitResolutionError("source-missing", path, str(target))
        try:
            actual = f"sha256:{hashlib.sha256(target.read_bytes()).hexdigest()}"
        except OSError as exc:
            raise FitResolutionError("source-read", path, str(exc)) from exc
    if actual != expected:
        raise FitResolutionError(
            "source-digest-mismatch", path, {"expected": expected, "actual": actual}
        )
    return SourceCheck(path, expected, actual)


def _verify_parameters_reference(fit: ParameterFit, parameters: Parameters, root: Path) -> None:
    target = _local_path(fit.parameters.path, root)
    try:
        referenced = Parameters.from_json(target, components=parameters.component_ids)
    except (OSError, ValueError) as exc:
        raise FitResolutionError(
            "invalid-parameters-source", fit.parameters.path, str(exc)
        ) from exc
    if referenced.fingerprint != parameters.fingerprint:
        raise FitResolutionError(
            "parameters-fingerprint-mismatch",
            fit.parameters.path,
            {
                "declared": referenced.fingerprint,
                "provided": parameters.fingerprint,
            },
        )


def _source_check_mapping(check: SourceCheck) -> dict[str, object]:
    return {
        "path": check.path,
        "expected_sha256": check.expected_sha256,
        "actual_sha256": check.actual_sha256,
    }


def _capability_mapping(capability: ParameterCapability) -> dict[str, object]:
    return {
        "observation_family": capability.observation_family,
        "parameter_identities": list(capability.parameter_identities),
        "parameter_units": list(capability.parameter_units),
        "model_domain": capability.model_domain,
        "temperature_domain_k": list(capability.temperature_domain_k),
        "value_labels": list(capability.value_labels),
        "value_units": list(capability.value_units),
        "value_shape": list(capability.value_shape),
        "derivative_shape": list(capability.derivative_shape),
        "parameter_fingerprint": capability.parameter_fingerprint,
        "active_parameter_fingerprint": capability.active_parameter_fingerprint,
        "topology_fingerprint": capability.topology_fingerprint,
    }


def _quadrature_mapping(quadrature: QuadraturePlan | None) -> dict[str, object] | None:
    if quadrature is None:
        return None
    return {
        "correlation": quadrature.correlation,
        "domain": {name: list(bounds) for name, bounds in quadrature.domain.bounds},
        "measure": quadrature.measure,
        "controls": {
            "initial_order": quadrature.controls.initial_order,
            "max_order": quadrature.controls.max_order,
            "tolerance": quadrature.controls.tolerance,
        },
        "rules": [
            {
                "order": rule.order,
                "nodes": list(rule.nodes),
                "weights": list(rule.weights),
                "observed_values": list(rule.observed_values),
                "observed_derivatives": list(rule.observed_derivatives),
            }
            for rule in quadrature.rules
        ],
        "accepted_order": quadrature.accepted_order,
        "normalized_objective": quadrature.normalized_objective,
        "objective_jacobian": list(quadrature.objective_jacobian),
        "objective_change": quadrature.objective_change,
        "jacobian_change": quadrature.jacobian_change,
    }


def _selected_identities(fit: ParameterFit, parameters: Parameters) -> tuple[str, ...]:
    specs = tuple(parameters.parameter_specs)
    available = {spec.identity: spec for spec in specs}
    if len(available) != len(specs):
        raise FitResolutionError("duplicate-parameter", "parameters.parameter_specs", None)
    selected = fit.active.identities
    missing = tuple(identity for identity in selected if identity not in available)
    if missing:
        raise FitResolutionError("missing-parameter", "$.active.identities", missing)
    for identity in selected:
        spec = available[identity]
        if spec.continuous is not True:
            raise FitResolutionError("discrete-parameter", identity, spec.family)
    try:
        ActiveParameterSet(parameters, selected)
    except ValueError as exc:
        raise FitResolutionError(
            "invalid-active-parameter-set", "$.active.identities", str(exc)
        ) from exc
    return selected


def _tie_groups(
    selection: ActiveSelection, selected: tuple[str, ...]
) -> tuple[tuple[str, ...], ...]:
    selected_set = set(selected)
    memberships: dict[str, int] = {}
    for index, tie in enumerate(selection.ties):
        for identity in tie.identities:
            if identity not in selected_set:
                raise FitResolutionError("tie-outside-selection", "$.active.ties", identity)
            if identity in memberships:
                raise FitResolutionError("overlapping-tie", "$.active.ties", identity)
            memberships[identity] = index
    groups: list[tuple[str, ...]] = []
    emitted: set[int] = set()
    for identity in selected:
        tie_index = memberships.get(identity)
        if tie_index is None:
            groups.append((identity,))
            continue
        if tie_index in emitted:
            continue
        members = tuple(member for member in selected if memberships.get(member) == tie_index)
        groups.append(members)
        emitted.add(tie_index)
    return tuple(groups)


def _resolve_coordinates(
    fit: ParameterFit,
    parameters: Parameters,
    selected: tuple[str, ...],
) -> tuple[tuple[str, ...], tuple[ResolvedCoordinate, ...]]:
    specs = {spec.identity: spec for spec in parameters.parameter_specs}
    declared = dict(fit.active.coordinates)
    groups = _tie_groups(fit.active, selected)
    representatives = {group[0] for group in groups}
    unexpected = tuple(sorted(set(declared) - representatives))
    if unexpected:
        raise FitResolutionError("coordinate-not-independent", "$.active.coordinates", unexpected)
    independent: list[str] = []
    coordinates: list[ResolvedCoordinate] = []
    for group in groups:
        representative = group[0]
        coordinate = declared.get(representative)
        if coordinate is None:
            raise FitResolutionError("missing-coordinate", "$.active.coordinates", representative)
        if coordinate.bounds[0] >= coordinate.bounds[1]:
            raise FitResolutionError(
                "invalid-bounds", f"$.active.coordinates.{representative}", coordinate.bounds
            )
        parameter_units = tuple(str(specs[identity].unit) for identity in group)
        if len(set(parameter_units)) != 1:
            raise FitResolutionError("incompatible-tie-units", representative, parameter_units)
        unit = coordinate.unit or parameter_units[0]
        try:
            unit_registry.Quantity(1, unit).to(parameter_units[0])
        except (AttributeError, TypeError, ValueError, pint.errors.PintError) as exc:
            raise FitResolutionError(
                "unit-mismatch",
                f"$.active.coordinates.{representative}.unit",
                {"declared": unit, "parameter": parameter_units[0]},
            ) from exc
        if fit.solver.multistart.include_input:
            value = specs[representative].value
            try:
                magnitude = (
                    float(cast(SupportsFloat, value.to(unit).magnitude))
                    if isinstance(value, pint.Quantity)
                    else float(cast(SupportsFloat, value))
                )
            except (TypeError, ValueError, pint.errors.PintError) as exc:
                raise FitResolutionError("invalid-parameter-value", representative, value) from exc
            if (
                not math.isfinite(magnitude)
                or not coordinate.bounds[0] <= magnitude <= coordinate.bounds[1]
            ):
                raise FitResolutionError("input-outside-bounds", representative, magnitude)
        independent.append(representative)
        coordinates.append(
            ResolvedCoordinate(representative, group, coordinate.bounds, coordinate.scale, unit)
        )
    return tuple(independent), tuple(coordinates)


def _table_rows(
    source: TableSource,
    *,
    root: Path,
    definition: ObservationDefinition,
    targets: tuple[str, ...],
) -> tuple[TableRow, ...]:
    if source.rows is not None:
        return source.rows
    target = Path(source.path)
    if not target.is_absolute():
        target = root / target
    if not target.is_file():
        raise FitResolutionError("source-missing", source.path, str(target))
    try:
        with target.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            fields = tuple(reader.fieldnames or ())
            expected = (*definition.state_fields, *targets)
            if len(fields) != len(set(fields)) or set(fields) != set(expected):
                raise FitResolutionError(
                    "table-columns", source.path, {"expected": expected, "actual": fields}
                )
            rows: list[TableRow] = []
            for index, raw in enumerate(reader):
                if None in raw:
                    raise FitResolutionError(
                        "table-columns",
                        f"{source.path}[{index}]",
                        {"unexpected": raw[None]},
                    )
                values: list[tuple[str, float]] = []
                for name in sorted(expected):
                    value = raw.get(name)
                    if value is None or value == "":
                        raise FitResolutionError(
                            "table-value", f"{source.path}[{index}].{name}", "missing"
                        )
                    try:
                        number = float(value)
                    except ValueError as exc:
                        raise FitResolutionError(
                            "table-value", f"{source.path}[{index}].{name}", value
                        ) from exc
                    if not math.isfinite(number):
                        raise FitResolutionError(
                            "table-value", f"{source.path}[{index}].{name}", value
                        )
                    values.append((name, number))
                rows.append(TableRow(tuple(values)))
    except OSError as exc:
        raise FitResolutionError("source-read", source.path, str(exc)) from exc
    if not rows:
        raise FitResolutionError("table-empty", source.path, None)
    return tuple(rows)


def _validate_rows(rows: Sequence[TableRow], source: TableSource, targets: tuple[str, ...]) -> None:
    states: set[tuple[float, ...]] = set()
    state_fields = tuple(name for name, _ in source.domain.bounds)
    for index, row in enumerate(rows):
        state = tuple(row.value(name) for name in state_fields)
        if state in states:
            raise FitResolutionError("duplicate-state", f"{source.path}[{index}]", state)
        states.add(state)
        for name, (lower, upper) in source.domain.bounds:
            value = row.value(name)
            if not lower <= value <= upper:
                raise FitResolutionError(
                    "row-outside-domain", f"{source.path}[{index}].{name}", value
                )
        for target in targets:
            value = row.value(target)
            if SCHEMA_REGISTRY.targets[target].positive and value <= 0:
                raise FitResolutionError(
                    "nonpositive-target", f"{source.path}[{index}].{target}", value
                )


def _validate_observation(
    observation: ObservationSpec, capability: ParameterCapability
) -> tuple[tuple[tuple[str, str], ...], tuple[tuple[str, float], ...]]:
    capability_labels = {
        label: unit
        for label, unit in zip(capability.value_labels, capability.value_units, strict=True)
    }
    unsupported_targets = tuple(
        target
        for target in observation.targets
        if SCHEMA_REGISTRY.targets[target].capability_label not in capability_labels
    )
    if unsupported_targets:
        raise FitResolutionError(
            "unsupported-observation-target", observation.family, unsupported_targets
        )
    for name, (lower, upper) in observation.source.domain.bounds:
        if lower >= upper or (name == "temperature_k" and lower <= 0):
            raise FitResolutionError(
                "invalid-domain", f"{observation.source.path}.{name}", (lower, upper)
            )
        if name == "temperature_k":
            model_lower, model_upper = capability.temperature_domain_k
            if (model_lower is not None and lower < model_lower) or (
                model_upper is not None and upper > model_upper
            ):
                raise FitResolutionError(
                    "source-outside-model-domain",
                    f"{observation.source.path}.{name}",
                    {
                        "source": (lower, upper),
                        "model": capability.temperature_domain_k,
                    },
                )
    canonical_units: list[tuple[str, str]] = []
    factors: list[tuple[str, float]] = []
    for target, unit in observation.units:
        capability_unit = capability_labels[SCHEMA_REGISTRY.targets[target].capability_label]
        try:
            zero = float(unit_registry.Quantity(0, unit).to(capability_unit).magnitude)
            factor = float(unit_registry.Quantity(1, unit).to(capability_unit).magnitude)
        except (AttributeError, TypeError, ValueError, pint.errors.PintError) as exc:
            raise FitResolutionError("invalid-unit", target, unit) from exc
        if zero != 0.0 or not math.isfinite(factor) or factor <= 0.0:
            raise FitResolutionError("invalid-unit", target, unit)
        canonical_units.append((target, capability_unit))
        factors.append((target, factor))
    return tuple(canonical_units), tuple(factors)


def _canonical_rows(rows: Sequence[TableRow], factors: dict[str, float]) -> tuple[TableRow, ...]:
    converted: list[TableRow] = []
    for index, row in enumerate(rows):
        values = tuple((name, value * factors.get(name, 1.0)) for name, value in row.values)
        if not all(math.isfinite(value) for _, value in values):
            raise FitResolutionError("unit-conversion-nonfinite", f"rows[{index}]", None)
        converted.append(TableRow(values))
    return tuple(converted)


def _gauss_legendre(order: int) -> tuple[tuple[float, ...], tuple[float, ...]]:
    nodes = [0.0] * order
    weights = [0.0] * order
    for index in range((order + 1) // 2):
        root = math.cos(math.pi * (index + 0.75) / (order + 0.5))
        for _ in range(64):
            p_previous, p_current = 1.0, root
            for degree in range(2, order + 1):
                p_previous, p_current = (
                    p_current,
                    ((2 * degree - 1) * root * p_current - (degree - 1) * p_previous) / degree,
                )
            derivative = order * (root * p_current - p_previous) / (root * root - 1.0)
            step = p_current / derivative
            root -= step
            if abs(step) <= 2.0e-15:
                break
        p_previous, p_current = 1.0, root
        for degree in range(2, order + 1):
            p_previous, p_current = (
                p_current,
                ((2 * degree - 1) * root * p_current - (degree - 1) * p_previous) / degree,
            )
        derivative = order * (root * p_current - p_previous) / (root * root - 1.0)
        nodes[index] = -root
        nodes[order - index - 1] = root
        weight = 2.0 / ((1.0 - root * root) * derivative * derivative)
        weights[index] = weight
        weights[order - index - 1] = weight
    return tuple(nodes), tuple(weights)


def _transformed_input(kind: str, temperature_k: float) -> tuple[float, float]:
    if kind == "identity":
        return temperature_k, 1.0
    if kind == "log":
        return math.log(temperature_k), 1.0 / temperature_k
    if kind == "log10":
        return math.log10(temperature_k), 1.0 / (temperature_k * math.log(10.0))
    raise RuntimeError(f"unhandled input transformation: {kind}")


def _transformed_output(kind: str, value: float) -> tuple[float, float]:
    if kind == "identity":
        return value, 1.0
    if kind == "log":
        result = math.exp(value)
        return result, result
    if kind == "log10":
        result = 10.0**value
        return result, math.log(10.0) * result
    raise RuntimeError(f"unhandled output transformation: {kind}")


def _correlation_value_derivative(
    source: PublishedObservationCorrelation,
    observation: ObservationSpec,
    temperature_k: float,
) -> tuple[float, float]:
    definition = SCHEMA_REGISTRY.correlations[source.correlation]
    transformation = dict(observation.transformation)[observation.targets[0]]
    try:
        independent, independent_derivative = _transformed_input(
            transformation.input, temperature_k
        )
        if definition.form != "power_series":
            raise RuntimeError(f"unhandled correlation form: {definition.form}")
        value = source.coefficients[-1]
        derivative = 0.0
        for coefficient in reversed(source.coefficients[:-1]):
            derivative = derivative * independent + value
            value = value * independent + coefficient
        observed, output_derivative = _transformed_output(transformation.output, value)
        derivative *= independent_derivative * output_derivative
    except (OverflowError, ValueError, ZeroDivisionError) as exc:
        raise FitResolutionError(
            "correlation-numerical", source.correlation, temperature_k
        ) from exc
    if not math.isfinite(observed) or not math.isfinite(derivative):
        raise FitResolutionError("correlation-nonfinite", source.correlation, temperature_k)
    target = SCHEMA_REGISTRY.targets[observation.targets[0]]
    if target.positive and observed <= 0:
        raise FitResolutionError("correlation-nonpositive", source.correlation, temperature_k)
    return observed, derivative


def _quadrature_objective(
    rule: QuadratureRule,
    observation: ObservationSpec,
    model: Pure,
    active_parameters: ActiveParameterSet,
) -> tuple[float, tuple[float, ...]]:
    target = observation.targets[0]
    label = SCHEMA_REGISTRY.targets[target].capability_label
    residual_kind = dict(observation.residuals)[target]
    if residual_kind != "log_ratio":
        raise RuntimeError(f"unhandled residual: {residual_kind}")
    residuals: list[float] = []
    jacobian: list[tuple[float, ...]] = []
    for temperature_k, observed in zip(rule.nodes, rule.observed_values, strict=True):
        try:
            result = model.saturation_parameter_derivatives(
                active_parameters,
                T=temperature_k * unit_registry.kelvin,
            )
            index = result.capability.value_labels.index(label)
        except (EosError, SaturationError, ValueError) as exc:
            raise FitResolutionError(
                "quadrature-evaluation", observation.source.path, temperature_k
            ) from exc
        predicted = result.values[index]
        if not math.isfinite(predicted) or predicted <= 0:
            raise FitResolutionError(
                "quadrature-evaluation", observation.source.path, temperature_k
            )
        residuals.append(observation.multiplier * math.log(predicted / observed))
        jacobian.append(
            tuple(
                observation.multiplier * derivative / predicted
                for derivative in result.derivatives[index]
            )
        )
    objective = sum(
        weight * residual * residual
        for weight, residual in zip(rule.weights, residuals, strict=True)
    )
    objective_jacobian = tuple(
        2.0
        * sum(
            weight * residual * row[column]
            for weight, residual, row in zip(rule.weights, residuals, jacobian, strict=True)
        )
        for column in range(len(active_parameters.identities))
    )
    if not math.isfinite(objective) or not all(
        math.isfinite(value) for value in objective_jacobian
    ):
        raise FitResolutionError("quadrature-numerical", observation.source.path, rule.order)
    return objective, objective_jacobian


def _relative_change(current: float, previous: float) -> float:
    return abs(current - previous) / max(1.0, abs(current), abs(previous))


def _quadrature(
    observation: ObservationSpec,
    model: Pure,
    active_parameters: ActiveParameterSet,
    unit_factor: float,
) -> QuadraturePlan | None:
    source = observation.source
    controls = observation.quadrature
    if not isinstance(source, PublishedObservationCorrelation) or controls is None:
        return None
    if controls.initial_order > controls.max_order:
        raise FitResolutionError(
            "invalid-quadrature", source.correlation, "initial order exceeds maximum"
        )
    lower, upper = source.domain.value("temperature_k")
    midpoint = (lower + upper) / 2.0
    half_width = (upper - lower) / 2.0
    rules: list[QuadratureRule] = []
    previous_objective: float | None = None
    previous_jacobian: tuple[float, ...] | None = None
    last_changes: tuple[float, float] | None = None
    order = controls.initial_order
    while order <= controls.max_order:
        reference_nodes, reference_weights = _gauss_legendre(order)
        nodes = tuple(midpoint + half_width * node for node in reference_nodes)
        values_and_derivatives = tuple(
            _correlation_value_derivative(source, observation, node) for node in nodes
        )
        rule = QuadratureRule(
            order,
            nodes,
            tuple(weight / 2.0 for weight in reference_weights),
            tuple(value * unit_factor for value, _ in values_and_derivatives),
            tuple(derivative * unit_factor for _, derivative in values_and_derivatives),
        )
        rules.append(rule)
        objective, objective_jacobian = _quadrature_objective(
            rule, observation, model, active_parameters
        )
        if previous_objective is not None and previous_jacobian is not None:
            objective_change = _relative_change(objective, previous_objective)
            jacobian_change = max(
                (
                    _relative_change(current, previous)
                    for current, previous in zip(objective_jacobian, previous_jacobian, strict=True)
                ),
                default=0.0,
            )
            last_changes = (objective_change, jacobian_change)
            if max(last_changes) <= controls.tolerance:
                return QuadraturePlan(
                    source.correlation,
                    source.domain,
                    "uniform_temperature",
                    controls,
                    tuple(rules),
                    order,
                    objective,
                    objective_jacobian,
                    objective_change,
                    jacobian_change,
                )
        previous_objective = objective
        previous_jacobian = objective_jacobian
        order *= 2
    raise FitResolutionError(
        "quadrature-not-converged",
        source.correlation,
        {
            "max_order": controls.max_order,
            "tolerance": controls.tolerance,
            "changes": last_changes,
        },
    )


def resolve_fit(
    fit: ParameterFit,
    parameters: Parameters,
    *,
    base_path: str | PathLike[str] | None = None,
) -> ResolvedFit:
    if not isinstance(fit, ParameterFit):
        raise FitResolutionError("invalid-fit", "fit", type(fit).__name__)
    if not isinstance(parameters, Parameters):
        raise FitResolutionError("invalid-parameters", "parameters", type(parameters).__name__)
    component_count = len(parameters.component_ids)
    if fit.model.kind == "pure" and component_count != 1:
        raise FitResolutionError("model-cardinality", "$.model.kind", component_count)
    if fit.model.kind == "mixture" and component_count < 2:
        raise FitResolutionError("model-cardinality", "$.model.kind", component_count)

    root = _path_root(base_path)
    parameters_source_check = _file_check(
        fit.parameters.path,
        fit.parameters.sha256,
        root=root,
    )
    _verify_parameters_reference(fit, parameters, root)
    selected = _selected_identities(fit, parameters)
    independent, coordinates = _resolve_coordinates(fit, parameters, selected)
    active_parameters = ActiveParameterSet(parameters, selected)
    model = Pure(parameters) if fit.model.kind == "pure" else Mixture(parameters)
    resolved_observations: list[ResolvedObservation] = []
    for index, observation in enumerate(fit.observations):
        if fit.model.kind == "mixture":
            if observation.composition is None or len(observation.composition) != component_count:
                raise FitResolutionError(
                    "invalid-composition", f"$.observations[{index}].composition", component_count
                )
        elif observation.composition is not None:
            raise FitResolutionError(
                "invalid-composition", f"$.observations[{index}].composition", "pure model"
            )
        definition = SCHEMA_REGISTRY.observations[observation.family]
        try:
            capability = model.parameter_capability(
                active_parameters, observation_family=observation.family
            )
        except EosError as exc:
            raise FitResolutionError(
                "unsupported-parameter-derivative", observation.family, str(exc)
            ) from exc
        canonical_units, unit_factors = _validate_observation(observation, capability)
        factors = dict(unit_factors)
        source = observation.source
        if isinstance(source, PublishedObservationCorrelation) and not isinstance(model, Pure):
            raise FitResolutionError(
                "unsupported-source", f"$.observations[{index}].source", "mixture correlation"
            )
        if isinstance(source, TableSource):
            inline_rows = (
                None
                if source.rows is None
                else [{name: value for name, value in row.values} for row in source.rows]
            )
            check = _file_check(
                source.path,
                source.sha256,
                root=root,
                inline=inline_rows,
            )
            rows = _table_rows(
                source,
                root=root,
                definition=definition,
                targets=observation.targets,
            )
            _validate_rows(rows, source, observation.targets)
            rows = _canonical_rows(rows, factors)
        else:
            correlation_payload = {
                "correlation": source.correlation,
                "coefficients": list(source.coefficients),
                "domain": {name: list(bounds) for name, bounds in source.domain.bounds},
            }
            check = _file_check(
                source.path,
                source.sha256,
                root=root,
                inline=correlation_payload,
            )
            rows = None
        resolved_observations.append(
            ResolvedObservation(
                observation,
                capability,
                canonical_units,
                rows,
                check,
                _quadrature(
                    observation,
                    cast(Pure, model),
                    active_parameters,
                    factors[observation.targets[0]],
                ),
            )
        )

    provisional = ResolvedFit(
        fit,
        parameters,
        active_parameters,
        selected,
        independent,
        coordinates,
        tuple(resolved_observations),
        parameters_source_check,
        "",
    )
    return ResolvedFit(
        provisional.fit,
        provisional.parameters,
        provisional.active_parameters,
        provisional.selected_identities,
        provisional.independent_identities,
        provisional.coordinates,
        provisional.observations,
        provisional.parameters_source_check,
        _hash(provisional.to_mapping()),
    )
