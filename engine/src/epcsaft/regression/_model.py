from __future__ import annotations

import hashlib
import json
import math
import re
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from itertools import combinations
from os import PathLike
from pathlib import Path
from types import MappingProxyType
from typing import TYPE_CHECKING, Literal, cast

from .._canonical import _canonical_json, _canonical_json_bytes
from ..eos import OBSERVABLE_FAMILIES, OBSERVABLE_TARGETS, ObservableFamily
from ..parameters import Parameters

if TYPE_CHECKING:
    from ._resolve import ResolvedFit


SCHEMA_VERSION = 2
_INTEGER_MAX = 2_147_483_647
_ERROR_DETAIL_LIMIT = 512
_DIGEST = re.compile(r"sha256:[0-9a-f]{64}\Z")
_IDENTIFIER_PATTERN = r"[a-z0-9]+(?:[-_][a-z0-9]+)*"
_IDENTITY_PATTERN = (
    rf"(?:{_IDENTIFIER_PATTERN}(?:/{_IDENTIFIER_PATTERN})+"
    r"|reaction:R[0-9]+:correlation:[a-z][a-z0-9_]*)"
)
_IDENTITY = re.compile(rf"{_IDENTITY_PATTERN}\Z")


class FitSchemaError(ValueError):
    def __init__(self, code: str, path: str, detail: object) -> None:
        try:
            rendered = str(detail)
        except (OverflowError, ValueError):
            detail = f"<{type(detail).__name__} omitted>"
            rendered = str(detail)
        if len(rendered) > _ERROR_DETAIL_LIMIT:
            detail = f"<{type(detail).__name__} omitted>"
            rendered = str(detail)
        self.code = code
        self.path = path
        self.detail = detail
        super().__init__(f"{code} at {path}: {rendered}")


@dataclass(frozen=True, slots=True)
class _Fields:
    required: tuple[str, ...]
    optional: tuple[str, ...] = ()


@dataclass(frozen=True, slots=True)
class ResidualDefinition:
    positive: bool


@dataclass(frozen=True, slots=True)
class CorrelationDefinition:
    form: Literal["power_series"]
    coefficient_count: int | None
    independent_fields: tuple[str, ...]


_CORRELATIONS = MappingProxyType(
    {
        "power_series_temperature": CorrelationDefinition("power_series", None, ("temperature_k",)),
    }
)

_FIELDS = MappingProxyType(
    {
        "fit": _Fields(
            ("schema_version", "parameters", "model", "active", "observations", "solver"),
            ("acceptance",),
        ),
        "parameters": _Fields(("path", "sha256")),
        "model": _Fields(("kind",)),
        "active": _Fields(("identities", "coordinates", "ties")),
        "coordinate": _Fields(("bounds", "origin", "scale"), ("unit",)),
        "tie": _Fields(("identities",)),
        "observation": _Fields(
            (
                "family",
                "source",
                "targets",
                "residuals",
                "units",
                "basis",
                "reference_convention",
                "transformation",
                "multiplier",
            ),
            ("quadrature", "composition"),
        ),
        "reactive_observation": _Fields(
            (
                "family",
                "identity",
                "request",
                "targets",
                "multiplier",
            )
        ),
        "table_source": _Fields(("kind", "path", "sha256", "citation", "domain"), ("rows",)),
        "correlation_source": _Fields(
            (
                "kind",
                "correlation",
                "coefficients",
                "path",
                "sha256",
                "citation",
                "domain",
            )
        ),
        "transformation": _Fields(("input", "output")),
        "quadrature": _Fields(("initial_order", "max_order", "tolerance")),
        "solver": _Fields(("multistart",), ("controls",)),
        "acceptance": _Fields(("maximum_cost",)),
        "solver_controls": _Fields(
            (),
            (
                "maximum_iterations",
                "optimizer_tolerance",
            ),
        ),
        "multistart": _Fields(("starts", "affine_bound_margin")),
        "start": _Fields(("identity", "values")),
    }
)
_RESIDUALS = MappingProxyType(
    {
        "log_ratio": ResidualDefinition(True),
        "scaled_difference": ResidualDefinition(False),
    }
)
_BASES = frozenset(base for target in OBSERVABLE_TARGETS.values() for base in target.bases)
_REFERENCES = frozenset(
    reference for target in OBSERVABLE_TARGETS.values() for reference in target.references
)
_TRANSFORMATIONS = frozenset({"identity", "log", "log10"})
_MODEL_KINDS = ("pure", "mixture")
_SOURCE_KINDS = ("table", "published_observation_correlation")
_QUADRATURE_ORDERS = (2, 4, 8, 16, 32, 64)


def _digest(value: object) -> str:
    return f"sha256:{hashlib.sha256(_canonical_json_bytes(value)).hexdigest()}"


def _mapping(value: object, path: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping) or any(type(key) is not str for key in value):
        raise FitSchemaError("invalid-type", path, "object with string keys required")
    return cast(Mapping[str, object], value)


def _keys(value: Mapping[str, object], name: str, path: str) -> None:
    fields = _FIELDS[name]
    allowed = frozenset((*fields.required, *fields.optional))
    unknown = tuple(sorted(set(value) - allowed))
    missing = tuple(field for field in fields.required if field not in value)
    if unknown:
        raise FitSchemaError("unknown-key", path, unknown)
    if missing:
        raise FitSchemaError("missing-key", path, missing)


def _text(value: object, path: str) -> str:
    if not isinstance(value, str) or not value or value != value.strip():
        raise FitSchemaError("invalid-type", path, "non-empty string required")
    return value


def _identity(value: object, path: str) -> str:
    result = _text(value, path)
    if _IDENTITY.fullmatch(result) is None:
        raise FitSchemaError("invalid-identity", path, result)
    return result


def _number(value: object, path: str) -> float:
    if type(value) not in (int, float):
        raise FitSchemaError("invalid-type", path, "finite number required")
    try:
        result = float(cast(int | float, value))
    except OverflowError as exc:
        raise FitSchemaError("invalid-number", path, "outside finite float range") from exc
    if not math.isfinite(result):
        raise FitSchemaError("invalid-number", path, value)
    return result


def _integer(value: object, path: str) -> int:
    if type(value) is not int:
        raise FitSchemaError("invalid-type", path, "integer required")
    result = value
    if abs(result) > _INTEGER_MAX:
        raise FitSchemaError("invalid-number", path, "integer outside supported range")
    return result


def _boolean(value: object, path: str) -> bool:
    if type(value) is not bool:
        raise FitSchemaError("invalid-type", path, "boolean required")
    return value


def _digest_text(value: object, path: str) -> str:
    result = _text(value, path)
    if _DIGEST.fullmatch(result) is None:
        raise FitSchemaError("invalid-digest", path, result)
    return result


def _list(value: object, path: str) -> list[object]:
    if not isinstance(value, list):
        raise FitSchemaError("invalid-type", path, "array required")
    return value


def _enum(value: object, allowed: Sequence[str], path: str) -> str:
    result = _text(value, path)
    if result not in allowed:
        raise FitSchemaError("unsupported-value", path, result)
    return result


@dataclass(frozen=True, slots=True)
class ParameterReference:
    path: str
    sha256: str


@dataclass(frozen=True, slots=True)
class ModelSpec:
    kind: Literal["pure", "mixture"]


@dataclass(frozen=True, slots=True)
class FitCoordinate:
    bounds: tuple[float, float]
    origin: float
    scale: float
    unit: str | None


@dataclass(frozen=True, slots=True)
class ParameterTie:
    identities: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class SourceDomain:
    bounds: tuple[tuple[str, tuple[float, float]], ...]

    def value(self, name: str) -> tuple[float, float]:
        for key, bounds in self.bounds:
            if key == name:
                return bounds
        raise KeyError(name)


@dataclass(frozen=True, slots=True)
class Transformation:
    input: str
    output: str


@dataclass(frozen=True, slots=True)
class QuadratureSpec:
    initial_order: int
    max_order: int
    tolerance: float


@dataclass(frozen=True, slots=True)
class TableRow:
    values: tuple[tuple[str, float], ...]

    def value(self, name: str) -> float:
        for key, value in self.values:
            if key == name:
                return value
        raise KeyError(name)


@dataclass(frozen=True, slots=True)
class TableSource:
    path: str
    sha256: str
    citation: str
    domain: SourceDomain
    rows: tuple[TableRow, ...] | None


@dataclass(frozen=True, slots=True)
class PublishedObservationCorrelation:
    correlation: str
    coefficients: tuple[float, ...]
    path: str
    sha256: str
    citation: str
    domain: SourceDomain


Source = TableSource | PublishedObservationCorrelation


@dataclass(frozen=True, slots=True)
class ObservationSpec:
    family: str
    composition: tuple[float, ...] | None
    source: Source
    targets: tuple[str, ...]
    residuals: tuple[tuple[str, str], ...]
    units: tuple[tuple[str, str], ...]
    basis: tuple[tuple[str, str], ...]
    reference_convention: tuple[tuple[str, str], ...]
    transformation: tuple[tuple[str, Transformation], ...]
    multiplier: float
    quadrature: QuadratureSpec | None


@dataclass(frozen=True, slots=True)
class EquilibriumFitObservation:
    """Canonical JSON-backed declaration validated by its Equilibrium owner."""

    family: Literal["equilibrium"]
    canonical_bytes_utf8: str

    def to_mapping(self) -> dict[str, object]:
        value = json.loads(self.canonical_bytes_utf8)
        if not isinstance(value, dict):
            raise TypeError("reactive fit observation is not a canonical object")
        return cast(dict[str, object], value)


@dataclass(frozen=True, slots=True)
class MultistartSpec:
    starts: tuple[FitStartSpec, ...]
    affine_bound_margin: float


@dataclass(frozen=True, slots=True)
class FitStartSpec:
    identity: str
    values: tuple[float, ...]


@dataclass(frozen=True, slots=True)
class SolverControlsSpec:
    maximum_iterations: int = 25
    optimizer_tolerance: float = 1.0e-6


@dataclass(frozen=True, slots=True)
class SolverSpec:
    multistart: MultistartSpec
    controls: SolverControlsSpec


@dataclass(frozen=True, slots=True)
class FitAcceptanceSpec:
    maximum_cost: float


@dataclass(frozen=True, slots=True)
class ActiveSelection:
    identities: tuple[str, ...]
    coordinates: tuple[tuple[str, FitCoordinate], ...]
    ties: tuple[ParameterTie, ...]


@dataclass(frozen=True, slots=True, init=False)
class ParameterFit:
    schema_version: int
    parameters: ParameterReference
    model: ModelSpec
    active: ActiveSelection
    observations: tuple[ObservationSpec | EquilibriumFitObservation, ...]
    solver: SolverSpec
    acceptance: FitAcceptanceSpec | None

    def __init__(self) -> None:
        raise TypeError("use ParameterFit.from_mapping or ParameterFit.from_json")

    @classmethod
    def from_mapping(cls, value: Mapping[str, object]) -> ParameterFit:
        root = _mapping(value, "$")
        _keys(root, "fit", "$")
        version = _integer(root["schema_version"], "$.schema_version")
        if version != SCHEMA_VERSION:
            raise FitSchemaError(
                "unsupported-schema", "$.schema_version", f"expected {SCHEMA_VERSION}"
            )
        model = _parse_model(root["model"], "$.model")
        observations = tuple(
            _parse_fit_observation(item, f"$.observations[{index}]")
            for index, item in enumerate(_list(root["observations"], "$.observations"))
        )
        if not observations:
            raise FitSchemaError("missing-value", "$.observations", "at least one required")
        families = tuple(observation.family for observation in observations)
        equilibrium_families = {"equilibrium"}
        if any(item in equilibrium_families for item in families) and model.kind != "mixture":
            raise FitSchemaError(
                "unsupported-model-domain",
                "$.observations",
                "Equilibrium observations require mixture-only owner families",
            )
        for index, observation in enumerate(observations):
            if isinstance(observation, EquilibriumFitObservation):
                continue
            if model.kind not in OBSERVABLE_FAMILIES[observation.family].model_kinds:
                raise FitSchemaError(
                    "unsupported-model-domain",
                    f"$.observations[{index}].family",
                    {"family": observation.family, "model": model.kind},
                )
        instance = object.__new__(cls)
        for name, item in (
            ("schema_version", version),
            ("parameters", _parse_parameters(root["parameters"], "$.parameters")),
            ("model", model),
            ("active", _parse_active(root["active"], "$.active")),
            ("observations", observations),
            ("solver", _parse_solver(root["solver"], "$.solver")),
            (
                "acceptance",
                _parse_acceptance(root["acceptance"], "$.acceptance")
                if "acceptance" in root
                else None,
            ),
        ):
            object.__setattr__(instance, name, item)
        return instance

    @classmethod
    def from_json(cls, path: str | PathLike[str]) -> ParameterFit:
        source = Path(path)

        def reject_constant(value: str) -> object:
            raise FitSchemaError("invalid-number", str(source), value)

        try:
            with source.open(encoding="utf-8") as handle:
                value = json.load(
                    handle,
                    object_pairs_hook=_unique_object,
                    parse_constant=reject_constant,
                )
        except json.JSONDecodeError as exc:
            raise FitSchemaError("invalid-json", str(source), exc.msg) from exc
        except ValueError as exc:
            raise FitSchemaError(
                "invalid-json", str(source), "numeric token exceeds supported size"
            ) from exc
        except (OSError, UnicodeError) as exc:
            raise FitSchemaError("invalid-json", str(source), str(exc)) from exc
        return cls.from_mapping(_mapping(value, str(source)))

    @classmethod
    def json_schema(cls) -> dict[str, object]:
        return parameter_fit_json_schema()

    @property
    def digest(self) -> str:
        return _digest(self.to_mapping())

    def to_mapping(self) -> dict[str, object]:
        result: dict[str, object] = {
            "schema_version": self.schema_version,
            "parameters": {"path": self.parameters.path, "sha256": self.parameters.sha256},
            "model": {"kind": self.model.kind},
            "active": {
                "identities": list(self.active.identities),
                "coordinates": {
                    identity: _coordinate_mapping(coordinate)
                    for identity, coordinate in self.active.coordinates
                },
                "ties": [{"identities": list(tie.identities)} for tie in self.active.ties],
            },
            "observations": [
                item.to_mapping()
                if isinstance(item, EquilibriumFitObservation)
                else _observation_mapping(item)
                for item in self.observations
            ],
            "solver": _solver_mapping(self.solver),
        }
        if self.acceptance is not None:
            result["acceptance"] = {"maximum_cost": self.acceptance.maximum_cost}
        return result

    def resolve(
        self,
        parameters: Parameters,
        *,
        base_path: str | PathLike[str] | None = None,
    ) -> ResolvedFit:
        from ._resolve import resolve_fit

        return resolve_fit(
            self,
            parameters,
            base_path=base_path,
        )


def _unique_object(pairs: list[tuple[str, object]]) -> dict[str, object]:
    result: dict[str, object] = {}
    for key, value in pairs:
        if key in result:
            raise FitSchemaError("duplicate-key", "$", key)
        result[key] = value
    return result


def _parse_parameters(value: object, path: str) -> ParameterReference:
    mapping = _mapping(value, path)
    _keys(mapping, "parameters", path)
    return ParameterReference(
        _text(mapping["path"], f"{path}.path"),
        _digest_text(mapping["sha256"], f"{path}.sha256"),
    )


def _parse_model(value: object, path: str) -> ModelSpec:
    mapping = _mapping(value, path)
    _keys(mapping, "model", path)
    kind = _enum(mapping["kind"], _MODEL_KINDS, f"{path}.kind")
    return ModelSpec(cast(Literal["pure", "mixture"], kind))


def _parse_coordinate(value: object, path: str) -> FitCoordinate:
    mapping = _mapping(value, path)
    _keys(mapping, "coordinate", path)
    bounds = _list(mapping["bounds"], f"{path}.bounds")
    if len(bounds) != 2:
        raise FitSchemaError("invalid-value", f"{path}.bounds", "two values required")
    lower, upper = (_number(item, f"{path}.bounds[{index}]") for index, item in enumerate(bounds))
    scale = _number(mapping["scale"], f"{path}.scale")
    if scale <= 0:
        raise FitSchemaError("invalid-value", f"{path}.scale", "positive required")
    unit = None if "unit" not in mapping else _text(mapping["unit"], f"{path}.unit")
    origin = _number(mapping["origin"], f"{path}.origin")
    return FitCoordinate((lower, upper), origin, scale, unit)


def _parse_active(value: object, path: str) -> ActiveSelection:
    mapping = _mapping(value, path)
    _keys(mapping, "active", path)
    identities = tuple(
        _identity(item, f"{path}.identities[{index}]")
        for index, item in enumerate(_list(mapping["identities"], f"{path}.identities"))
    )
    if not identities or len(set(identities)) != len(identities):
        raise FitSchemaError(
            "invalid-value", f"{path}.identities", "unique non-empty list required"
        )
    coordinates_mapping = _mapping(mapping["coordinates"], f"{path}.coordinates")
    coordinates = tuple(
        (
            _identity(identity, f"{path}.coordinates.{identity}"),
            _parse_coordinate(item, f"{path}.coordinates.{identity}"),
        )
        for identity, item in sorted(coordinates_mapping.items())
    )
    ties = tuple(
        _parse_tie(item, f"{path}.ties[{index}]")
        for index, item in enumerate(_list(mapping["ties"], f"{path}.ties"))
    )
    return ActiveSelection(identities, coordinates, ties)


def _parse_tie(value: object, path: str) -> ParameterTie:
    mapping = _mapping(value, path)
    _keys(mapping, "tie", path)
    identities = tuple(
        _identity(item, f"{path}.identities[{index}]")
        for index, item in enumerate(_list(mapping["identities"], f"{path}.identities"))
    )
    if len(identities) < 2 or len(set(identities)) != len(identities):
        raise FitSchemaError(
            "invalid-value", f"{path}.identities", "two unique identities required"
        )
    return ParameterTie(tuple(sorted(identities)))


def _parse_domain(value: object, path: str, allowed: Sequence[str]) -> SourceDomain:
    mapping = _mapping(value, path)
    expected = set(allowed)
    unknown = tuple(sorted(set(mapping) - expected))
    missing = tuple(sorted(expected - set(mapping)))
    if unknown or missing:
        raise FitSchemaError("invalid-keys", path, {"unknown": unknown, "missing": missing})
    bounds: list[tuple[str, tuple[float, float]]] = []
    for name in sorted(mapping):
        values = _list(mapping[name], f"{path}.{name}")
        if len(values) != 2:
            raise FitSchemaError("invalid-value", f"{path}.{name}", "two values required")
        lower, upper = (
            _number(item, f"{path}.{name}[{index}]") for index, item in enumerate(values)
        )
        bounds.append((name, (lower, upper)))
    return SourceDomain(tuple(bounds))


def _row(
    value: object,
    path: str,
    definition: ObservableFamily,
    targets: tuple[str, ...],
    domain: SourceDomain,
) -> TableRow:
    mapping = _mapping(value, path)
    expected = frozenset((*definition.state_fields, *targets))
    unknown = tuple(sorted(set(mapping) - expected))
    missing = tuple(sorted(expected - set(mapping)))
    if unknown or missing:
        raise FitSchemaError("invalid-row", path, {"unknown": unknown, "missing": missing})
    values = tuple((name, _number(mapping[name], f"{path}.{name}")) for name in sorted(mapping))
    row = TableRow(values)
    for field, (lower, upper) in domain.bounds:
        current = row.value(field)
        if not lower <= current <= upper:
            raise FitSchemaError("row-outside-domain", f"{path}.{field}", current)
    for target in targets:
        if OBSERVABLE_TARGETS[target].positive and row.value(target) <= 0:
            raise FitSchemaError("nonpositive-target", f"{path}.{target}", row.value(target))
    return row


def _parse_table_source(
    value: object,
    path: str,
    definition: ObservableFamily,
    targets: tuple[str, ...],
) -> TableSource:
    mapping = _mapping(value, path)
    _keys(mapping, "table_source", path)
    domain = _parse_domain(mapping["domain"], f"{path}.domain", definition.state_fields)
    rows = None
    if "rows" in mapping:
        parsed = tuple(
            _row(item, f"{path}.rows[{index}]", definition, targets, domain)
            for index, item in enumerate(_list(mapping["rows"], f"{path}.rows"))
        )
        if not parsed:
            raise FitSchemaError("missing-value", f"{path}.rows", "non-empty")
        rows = parsed
    return TableSource(
        _text(mapping["path"], f"{path}.path"),
        _digest_text(mapping["sha256"], f"{path}.sha256"),
        _text(mapping["citation"], f"{path}.citation"),
        domain,
        rows,
    )


def _parse_correlation_source(
    value: object,
    path: str,
    definition: ObservableFamily,
    targets: tuple[str, ...],
) -> PublishedObservationCorrelation:
    mapping = _mapping(value, path)
    _keys(mapping, "correlation_source", path)
    name = _enum(mapping["correlation"], tuple(sorted(_CORRELATIONS)), f"{path}.correlation")
    correlation = _CORRELATIONS[name]
    if len(targets) != 1:
        raise FitSchemaError(
            "invalid-value", f"{path}.correlation", "correlation source requires one target"
        )
    if correlation.independent_fields != definition.state_fields:
        raise FitSchemaError("unsupported-correlation", f"{path}.correlation", name)
    coefficients = tuple(
        _number(item, f"{path}.coefficients[{index}]")
        for index, item in enumerate(_list(mapping["coefficients"], f"{path}.coefficients"))
    )
    if not coefficients:
        raise FitSchemaError("missing-value", f"{path}.coefficients", "non-empty")
    if (
        correlation.coefficient_count is not None
        and len(coefficients) != correlation.coefficient_count
    ):
        raise FitSchemaError("invalid-value", f"{path}.coefficients", correlation.coefficient_count)
    return PublishedObservationCorrelation(
        name,
        coefficients,
        _text(mapping["path"], f"{path}.path"),
        _digest_text(mapping["sha256"], f"{path}.sha256"),
        _text(mapping["citation"], f"{path}.citation"),
        _parse_domain(mapping["domain"], f"{path}.domain", definition.state_fields),
    )


def _parse_source(
    value: object,
    path: str,
    definition: ObservableFamily,
    targets: tuple[str, ...],
) -> Source:
    mapping = _mapping(value, path)
    kind = _enum(mapping.get("kind"), _SOURCE_KINDS, f"{path}.kind")
    if kind == _SOURCE_KINDS[0]:
        return _parse_table_source(value, path, definition, targets)
    return _parse_correlation_source(value, path, definition, targets)


def _parse_target_map(
    value: object,
    path: str,
    targets: tuple[str, ...],
    allowed: Sequence[str],
) -> tuple[tuple[str, str], ...]:
    mapping = _mapping(value, path)
    if set(mapping) != set(targets):
        raise FitSchemaError(
            "invalid-keys", path, {"expected": targets, "actual": tuple(sorted(mapping))}
        )
    return tuple(
        (target, _enum(mapping[target], allowed, f"{path}.{target}")) for target in targets
    )


def _parse_units(value: object, path: str, targets: tuple[str, ...]) -> tuple[tuple[str, str], ...]:
    mapping = _mapping(value, path)
    if set(mapping) != set(targets):
        raise FitSchemaError(
            "invalid-keys", path, {"expected": targets, "actual": tuple(sorted(mapping))}
        )
    units: list[tuple[str, str]] = []
    for target in targets:
        unit = _text(mapping[target], f"{path}.{target}")
        units.append((target, unit))
    return tuple(units)


def _parse_transformations(
    value: object, path: str, targets: tuple[str, ...]
) -> tuple[tuple[str, Transformation], ...]:
    mapping = _mapping(value, path)
    if set(mapping) != set(targets):
        raise FitSchemaError(
            "invalid-keys", path, {"expected": targets, "actual": tuple(sorted(mapping))}
        )
    transformations: list[tuple[str, Transformation]] = []
    allowed = tuple(sorted(_TRANSFORMATIONS))
    for target in targets:
        item = _mapping(mapping[target], f"{path}.{target}")
        _keys(item, "transformation", f"{path}.{target}")
        transformations.append(
            (
                target,
                Transformation(
                    _enum(item["input"], allowed, f"{path}.{target}.input"),
                    _enum(item["output"], allowed, f"{path}.{target}.output"),
                ),
            )
        )
    return tuple(transformations)


def _parse_fit_observation(
    value: object,
    path: str,
) -> ObservationSpec | EquilibriumFitObservation:
    mapping = _mapping(value, path)
    if mapping.get("family") == "equilibrium":
        _keys(mapping, "reactive_observation", path)
        if _text(mapping["identity"], f"{path}.identity") == "equilibrium":
            raise FitSchemaError("invalid-value", f"{path}.identity", "batch identity required")
        multiplier = _number(mapping["multiplier"], f"{path}.multiplier")
        if multiplier <= 0.0:
            raise FitSchemaError("invalid-value", f"{path}.multiplier", "positive required")
        _mapping(mapping["request"], f"{path}.request")
        if not _list(mapping["targets"], f"{path}.targets"):
            raise FitSchemaError("missing-value", f"{path}.targets", "non-empty required")
        return EquilibriumFitObservation("equilibrium", _canonical_json(mapping))
    return _parse_observation(value, path)


def _parse_observation(value: object, path: str) -> ObservationSpec:
    mapping = _mapping(value, path)
    _keys(mapping, "observation", path)
    family = _enum(mapping["family"], tuple(sorted(OBSERVABLE_FAMILIES)), f"{path}.family")
    definition = OBSERVABLE_FAMILIES[family]
    targets = tuple(
        _enum(item, definition.targets, f"{path}.targets[{index}]")
        for index, item in enumerate(_list(mapping["targets"], f"{path}.targets"))
    )
    if not targets or len(set(targets)) != len(targets):
        raise FitSchemaError("invalid-value", f"{path}.targets", "unique non-empty list required")
    canonical_targets = tuple(target for target in definition.targets if target in targets)
    if targets != canonical_targets:
        raise FitSchemaError(
            "invalid-value", f"{path}.targets", "canonical registry order required"
        )
    source = _parse_source(mapping["source"], f"{path}.source", definition, targets)
    composition = None
    if "composition" in mapping:
        composition = tuple(
            _number(item, f"{path}.composition[{index}]")
            for index, item in enumerate(_list(mapping["composition"], f"{path}.composition"))
        )
        if (
            not composition
            or any(value <= 0.0 for value in composition)
            or not math.isclose(sum(composition), 1.0, rel_tol=0.0, abs_tol=1.0e-12)
        ):
            raise FitSchemaError(
                "invalid-value", f"{path}.composition", "positive normalized values required"
            )
    residuals = _parse_target_map(
        mapping["residuals"],
        f"{path}.residuals",
        targets,
        tuple(sorted(_RESIDUALS)),
    )
    units = _parse_units(mapping["units"], f"{path}.units", targets)
    basis = _parse_target_map(mapping["basis"], f"{path}.basis", targets, tuple(sorted(_BASES)))
    references = _parse_target_map(
        mapping["reference_convention"],
        f"{path}.reference_convention",
        targets,
        tuple(sorted(_REFERENCES)),
    )
    for target, declared in basis:
        if declared not in OBSERVABLE_TARGETS[target].bases:
            raise FitSchemaError("invalid-basis", f"{path}.basis.{target}", declared)
    for target, declared in references:
        if declared not in OBSERVABLE_TARGETS[target].references:
            raise FitSchemaError(
                "invalid-reference", f"{path}.reference_convention.{target}", declared
            )
    transformations = _parse_transformations(
        mapping["transformation"], f"{path}.transformation", targets
    )
    if isinstance(source, TableSource) and any(
        transform.input != "identity" or transform.output != "identity"
        for _, transform in transformations
    ):
        raise FitSchemaError(
            "invalid-transformation",
            f"{path}.transformation",
            "table sources require identity transformations",
        )
    multiplier = _number(mapping["multiplier"], f"{path}.multiplier")
    if multiplier <= 0:
        raise FitSchemaError("invalid-value", f"{path}.multiplier", "positive required")
    quadrature = None
    if "quadrature" in mapping:
        item = _mapping(mapping["quadrature"], f"{path}.quadrature")
        _keys(item, "quadrature", f"{path}.quadrature")
        initial = _integer(item["initial_order"], f"{path}.quadrature.initial_order")
        maximum = _integer(item["max_order"], f"{path}.quadrature.max_order")
        tolerance = _number(item["tolerance"], f"{path}.quadrature.tolerance")
        if initial not in _QUADRATURE_ORDERS or maximum not in _QUADRATURE_ORDERS or tolerance <= 0:
            raise FitSchemaError("invalid-value", f"{path}.quadrature", "orders or tolerance")
        quadrature = QuadratureSpec(initial, maximum, tolerance)
    if isinstance(source, PublishedObservationCorrelation) and quadrature is None:
        raise FitSchemaError("missing-key", path, "quadrature required for correlation source")
    if isinstance(source, TableSource) and quadrature is not None:
        raise FitSchemaError(
            "invalid-value", f"{path}.quadrature", "table source cannot be integrated"
        )
    return ObservationSpec(
        family,
        composition,
        source,
        targets,
        residuals,
        units,
        basis,
        references,
        transformations,
        multiplier,
        quadrature,
    )


def _parse_solver(value: object, path: str) -> SolverSpec:
    mapping = _mapping(value, path)
    _keys(mapping, "solver", path)
    item = _mapping(mapping["multistart"], f"{path}.multistart")
    _keys(item, "multistart", f"{path}.multistart")
    starts = []
    for index, raw_start in enumerate(_list(item["starts"], f"{path}.multistart.starts")):
        start_path = f"{path}.multistart.starts[{index}]"
        start = _mapping(raw_start, start_path)
        _keys(start, "start", start_path)
        values = tuple(
            _number(raw, f"{start_path}.values[{coordinate}]")
            for coordinate, raw in enumerate(_list(start["values"], f"{start_path}.values"))
        )
        if not values:
            raise FitSchemaError("invalid-value", f"{start_path}.values", "non-empty required")
        starts.append(FitStartSpec(_text(start["identity"], f"{start_path}.identity"), values))
    identities = tuple(start.identity for start in starts)
    if not starts or len(set(identities)) != len(identities):
        raise FitSchemaError(
            "invalid-value", f"{path}.multistart.starts", "unique non-empty identities required"
        )
    margin = _number(item["affine_bound_margin"], f"{path}.multistart.affine_bound_margin")
    if margin < 0.0:
        raise FitSchemaError(
            "invalid-value", f"{path}.multistart.affine_bound_margin", "nonnegative required"
        )
    controls_item = _mapping(mapping.get("controls", {}), f"{path}.controls")
    _keys(controls_item, "solver_controls", f"{path}.controls")
    defaults = SolverControlsSpec()
    controls = SolverControlsSpec(
        maximum_iterations=_integer(
            controls_item.get("maximum_iterations", defaults.maximum_iterations),
            f"{path}.controls.maximum_iterations",
        ),
        optimizer_tolerance=_number(
            controls_item.get("optimizer_tolerance", defaults.optimizer_tolerance),
            f"{path}.controls.optimizer_tolerance",
        ),
    )
    if controls.maximum_iterations <= 0:
        raise FitSchemaError(
            "invalid-value",
            f"{path}.controls.maximum_iterations",
            "positive required",
        )
    if controls.optimizer_tolerance <= 0.0:
        raise FitSchemaError(
            "invalid-value", f"{path}.controls.optimizer_tolerance", "positive required"
        )
    return SolverSpec(MultistartSpec(tuple(starts), margin), controls)


def _parse_acceptance(value: object, path: str) -> FitAcceptanceSpec:
    mapping = _mapping(value, path)
    _keys(mapping, "acceptance", path)
    maximum_cost = _number(mapping["maximum_cost"], f"{path}.maximum_cost")
    if maximum_cost < 0.0:
        raise FitSchemaError("invalid-value", f"{path}.maximum_cost", "nonnegative required")
    return FitAcceptanceSpec(maximum_cost)


def _coordinate_mapping(coordinate: FitCoordinate) -> dict[str, object]:
    result: dict[str, object] = {
        "bounds": list(coordinate.bounds),
        "origin": coordinate.origin,
        "scale": coordinate.scale,
    }
    if coordinate.unit is not None:
        result["unit"] = coordinate.unit
    return result


def _domain_mapping(domain: SourceDomain) -> dict[str, object]:
    return {name: list(bounds) for name, bounds in domain.bounds}


def _source_mapping(source: Source) -> dict[str, object]:
    if isinstance(source, TableSource):
        result: dict[str, object] = {
            "kind": _SOURCE_KINDS[0],
            "path": source.path,
            "sha256": source.sha256,
            "citation": source.citation,
            "domain": _domain_mapping(source.domain),
        }
        if source.rows is not None:
            result["rows"] = [{name: value for name, value in row.values} for row in source.rows]
        return result
    return {
        "kind": _SOURCE_KINDS[1],
        "correlation": source.correlation,
        "coefficients": list(source.coefficients),
        "path": source.path,
        "sha256": source.sha256,
        "citation": source.citation,
        "domain": _domain_mapping(source.domain),
    }


def _observation_mapping(observation: ObservationSpec) -> dict[str, object]:
    result: dict[str, object] = {
        "family": observation.family,
        "source": _source_mapping(observation.source),
        "targets": list(observation.targets),
        "residuals": dict(observation.residuals),
        "units": dict(observation.units),
        "basis": dict(observation.basis),
        "reference_convention": dict(observation.reference_convention),
        "transformation": {
            target: {"input": transform.input, "output": transform.output}
            for target, transform in observation.transformation
        },
        "multiplier": observation.multiplier,
    }
    if observation.quadrature is not None:
        result["quadrature"] = {
            "initial_order": observation.quadrature.initial_order,
            "max_order": observation.quadrature.max_order,
            "tolerance": observation.quadrature.tolerance,
        }
    if observation.composition is not None:
        result["composition"] = list(observation.composition)
    return result


def _solver_mapping(solver: SolverSpec) -> dict[str, object]:
    return {
        "multistart": {
            "starts": [
                {"identity": start.identity, "values": list(start.values)}
                for start in solver.multistart.starts
            ],
            "affine_bound_margin": solver.multistart.affine_bound_margin,
        },
        "controls": {
            "maximum_iterations": solver.controls.maximum_iterations,
            "optimizer_tolerance": solver.controls.optimizer_tolerance,
        },
    }


def _object_schema(
    properties: dict[str, object], required: Sequence[str], *, additional: bool = False
) -> dict[str, object]:
    return {
        "type": "object",
        "properties": properties,
        "required": list(required),
        "additionalProperties": additional,
    }


def _exact_map_schema(
    targets: Sequence[str], properties: Mapping[str, object]
) -> dict[str, object]:
    return _object_schema(
        {target: properties[target] for target in targets},
        targets,
    )


def _domain_schema(state_fields: Sequence[str]) -> dict[str, object]:
    bounds = {
        "type": "array",
        "items": {"type": "number"},
        "minItems": 2,
        "maxItems": 2,
    }
    return _object_schema({field: bounds for field in state_fields}, state_fields)


def _source_base_schema(state_fields: Sequence[str]) -> dict[str, object]:
    return {
        "path": {"type": "string", "minLength": 1},
        "sha256": {"type": "string", "pattern": r"^sha256:[0-9a-f]{64}$"},
        "citation": {"type": "string", "minLength": 1},
        "domain": _domain_schema(state_fields),
    }


def _correlation_schema(
    name: str, definition: CorrelationDefinition, base: dict[str, object]
) -> dict[str, object]:
    coefficients: dict[str, object] = {
        "type": "array",
        "items": {"type": "number"},
        "minItems": 1,
    }
    if definition.coefficient_count is not None:
        coefficients["maxItems"] = definition.coefficient_count
        coefficients["minItems"] = definition.coefficient_count
    return _object_schema(
        {
            "kind": {"const": _SOURCE_KINDS[1]},
            "correlation": {"const": name},
            "coefficients": coefficients,
            **base,
        },
        ("kind", "correlation", "coefficients", "path", "sha256", "citation", "domain"),
    )


def _target_array_schema(targets: Sequence[str]) -> dict[str, object]:
    return {
        "type": "array",
        "prefixItems": [{"const": target} for target in targets],
        "items": False,
        "minItems": len(targets),
        "maxItems": len(targets),
    }


def _table_source_schema(
    definition: ObservableFamily,
    targets: Sequence[str],
    number: dict[str, object],
) -> dict[str, object]:
    row_fields = (*definition.state_fields, *targets)
    return _object_schema(
        {
            "kind": {"const": _SOURCE_KINDS[0]},
            **_source_base_schema(definition.state_fields),
            "rows": {
                "type": "array",
                "items": _object_schema(
                    {field: number for field in row_fields},
                    row_fields,
                ),
                "minItems": 1,
            },
        },
        _FIELDS["table_source"].required,
    )


def _quadrature_schema() -> dict[str, object]:
    return _object_schema(
        {
            "initial_order": {"enum": list(_QUADRATURE_ORDERS)},
            "max_order": {"enum": list(_QUADRATURE_ORDERS)},
            "tolerance": {"type": "number", "exclusiveMinimum": 0},
        },
        _FIELDS["quadrature"].required,
    )


def _observation_schema(
    family: str,
    definition: ObservableFamily,
    targets: tuple[str, ...],
    source: dict[str, object],
    *,
    correlated: bool,
    strings: dict[str, object],
    transformation: dict[str, object],
) -> dict[str, object]:
    residual = {target: {"enum": sorted(_RESIDUALS)} for target in targets}
    unit = {target: strings for target in targets}
    basis = {target: {"enum": list(OBSERVABLE_TARGETS[target].bases)} for target in targets}
    reference = {
        target: {"enum": list(OBSERVABLE_TARGETS[target].references)} for target in targets
    }
    transformations = {target: transformation for target in targets}
    properties: dict[str, object] = {
        "family": {"const": family},
        "source": source,
        "targets": _target_array_schema(targets),
        "residuals": _exact_map_schema(targets, residual),
        "units": _exact_map_schema(targets, unit),
        "basis": _exact_map_schema(targets, basis),
        "reference_convention": _exact_map_schema(targets, reference),
        "transformation": _exact_map_schema(targets, transformations),
        "multiplier": {"type": "number", "exclusiveMinimum": 0},
        "composition": {
            "type": "array",
            "items": {"type": "number", "exclusiveMinimum": 0},
            "minItems": 1,
        },
    }
    required = list(_FIELDS["observation"].required)
    if correlated:
        properties["quadrature"] = _quadrature_schema()
        required.append("quadrature")
    return _object_schema(properties, required)


def _observation_schemas(
    strings: dict[str, object],
    number: dict[str, object],
    transformation: dict[str, object],
) -> tuple[dict[str, object], ...]:
    schemas: list[dict[str, object]] = []
    table_transformation = _object_schema(
        {"input": {"const": "identity"}, "output": {"const": "identity"}},
        _FIELDS["transformation"].required,
    )
    for family, definition in OBSERVABLE_FAMILIES.items():
        for count in range(1, len(definition.targets) + 1):
            for targets in combinations(definition.targets, count):
                schemas.append(
                    _observation_schema(
                        family,
                        definition,
                        targets,
                        _table_source_schema(definition, targets, number),
                        correlated=False,
                        strings=strings,
                        transformation=table_transformation,
                    )
                )
        for target in definition.targets:
            for name, correlation in _CORRELATIONS.items():
                if correlation.independent_fields != definition.state_fields:
                    continue
                schemas.append(
                    _observation_schema(
                        family,
                        definition,
                        (target,),
                        _correlation_schema(
                            name,
                            correlation,
                            _source_base_schema(definition.state_fields),
                        ),
                        correlated=True,
                        strings=strings,
                        transformation=transformation,
                    )
                )
    return tuple(schemas)


def parameter_fit_json_schema() -> dict[str, object]:
    strings: dict[str, object] = {"type": "string", "minLength": 1}
    identity: dict[str, object] = {
        "type": "string",
        "pattern": rf"^{_IDENTITY_PATTERN}$",
    }
    number: dict[str, object] = {"type": "number"}
    digest: dict[str, object] = {
        "type": "string",
        "pattern": r"^sha256:[0-9a-f]{64}$",
    }
    coordinate = _object_schema(
        {
            "bounds": {"type": "array", "items": number, "minItems": 2, "maxItems": 2},
            "origin": number,
            "scale": {"type": "number", "exclusiveMinimum": 0},
            "unit": strings,
        },
        _FIELDS["coordinate"].required,
    )
    transformation = _object_schema(
        {
            "input": {"enum": sorted(_TRANSFORMATIONS)},
            "output": {"enum": sorted(_TRANSFORMATIONS)},
        },
        _FIELDS["transformation"].required,
    )
    observation_schemas = (
        *_observation_schemas(strings, number, transformation),
        _object_schema(
            {
                "family": {"const": "equilibrium"},
                "identity": strings,
                "request": {"type": "object"},
                "targets": {
                    "type": "array",
                    "items": {"type": "object"},
                    "minItems": 1,
                },
                "multiplier": {"type": "number", "exclusiveMinimum": 0},
            },
            _FIELDS["reactive_observation"].required,
        ),
    )
    valid_model_kinds = sorted(
        {
            "mixture",
            *(
                kind
                for definition in OBSERVABLE_FAMILIES.values()
                for kind in definition.model_kinds
            ),
        }
    )
    return {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "$id": "epcsaft.parameter-fit.v2",
        **_object_schema(
            {
                "schema_version": {"const": SCHEMA_VERSION, "type": "integer"},
                "parameters": _object_schema(
                    {"path": strings, "sha256": digest},
                    _FIELDS["parameters"].required,
                ),
                "model": _object_schema(
                    {"kind": {"enum": valid_model_kinds}},
                    _FIELDS["model"].required,
                ),
                "active": _object_schema(
                    {
                        "identities": {
                            "type": "array",
                            "items": identity,
                            "minItems": 1,
                            "uniqueItems": True,
                        },
                        "coordinates": {
                            "type": "object",
                            "propertyNames": identity,
                            "additionalProperties": coordinate,
                            "minProperties": 1,
                        },
                        "ties": {
                            "type": "array",
                            "items": _object_schema(
                                {
                                    "identities": {
                                        "type": "array",
                                        "items": identity,
                                        "minItems": 2,
                                        "uniqueItems": True,
                                    }
                                },
                                _FIELDS["tie"].required,
                            ),
                            "uniqueItems": True,
                        },
                    },
                    _FIELDS["active"].required,
                ),
                "observations": {
                    "type": "array",
                    "items": {"oneOf": list(observation_schemas)},
                    "minItems": 1,
                },
                "solver": _object_schema(
                    {
                        "multistart": _object_schema(
                            {
                                "starts": {
                                    "type": "array",
                                    "items": _object_schema(
                                        {
                                            "identity": strings,
                                            "values": {
                                                "type": "array",
                                                "items": number,
                                                "minItems": 1,
                                            },
                                        },
                                        _FIELDS["start"].required,
                                    ),
                                    "minItems": 1,
                                },
                                "affine_bound_margin": {"type": "number", "minimum": 0},
                            },
                            _FIELDS["multistart"].required,
                        ),
                        "controls": _object_schema(
                            {
                                "maximum_iterations": {
                                    "type": "integer",
                                    "minimum": 1,
                                },
                                "optimizer_tolerance": {
                                    "type": "number",
                                    "exclusiveMinimum": 0,
                                },
                            },
                            _FIELDS["solver_controls"].required,
                        ),
                    },
                    _FIELDS["solver"].required,
                ),
                "acceptance": _object_schema(
                    {"maximum_cost": {"type": "number", "minimum": 0}},
                    _FIELDS["acceptance"].required,
                ),
            },
            _FIELDS["fit"].required,
        ),
    }
