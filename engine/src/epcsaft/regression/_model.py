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

from ..parameters import Parameters

if TYPE_CHECKING:
    from ._resolve import ResolvedFit


SCHEMA_VERSION = 1
_INTEGER_MAX = 2_147_483_647
_ERROR_DETAIL_LIMIT = 512
_DIGEST = re.compile(r"sha256:[0-9a-f]{64}\Z")
_IDENTIFIER_PATTERN = r"[a-z0-9]+(?:[-_][a-z0-9]+)*"
_IDENTITY_PATTERN = rf"{_IDENTIFIER_PATTERN}(?:/{_IDENTIFIER_PATTERN})+"
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
class TargetDefinition:
    positive: bool
    observable_row_identity: str
    unit: str
    bases: tuple[str, ...]
    references: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class ObservationDefinition:
    targets: tuple[str, ...]
    state_fields: tuple[str, ...]
    model_kinds: tuple[Literal["pure", "mixture"], ...]


@dataclass(frozen=True, slots=True)
class ResidualDefinition:
    positive: bool


@dataclass(frozen=True, slots=True)
class CorrelationDefinition:
    form: Literal["power_series"]
    coefficient_count: int | None
    independent_fields: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class _Registry:
    fields: Mapping[str, _Fields]
    targets: Mapping[str, TargetDefinition]
    observations: Mapping[str, ObservationDefinition]
    residuals: Mapping[str, ResidualDefinition]
    correlations: Mapping[str, CorrelationDefinition]
    bases: frozenset[str]
    references: frozenset[str]
    transformations: frozenset[str]
    model_kinds: tuple[str, ...]
    source_kinds: tuple[str, ...]
    quadrature_orders: tuple[int, ...]


_TARGETS = MappingProxyType(
    {
        "pressure": TargetDefinition(True, "pressure_pa", "pascal", ("absolute",), ("absolute",)),
        "saturation_pressure": TargetDefinition(
            True, "saturation_pressure_pa", "pascal", ("absolute",), ("absolute",)
        ),
        "saturated_liquid_mass_density": TargetDefinition(
            True,
            "liquid_mass_density_kg_per_m3",
            "kilogram / meter**3",
            ("mass",),
            ("saturated_liquid",),
        ),
        "saturated_liquid_molar_density": TargetDefinition(
            True,
            "liquid_molar_density_mol_per_m3",
            "mole / meter**3",
            ("molar",),
            ("saturated_liquid",),
        ),
        "saturated_vapor_mass_density": TargetDefinition(
            True,
            "vapor_mass_density_kg_per_m3",
            "kilogram / meter**3",
            ("mass",),
            ("saturated_vapor",),
        ),
        "saturated_vapor_molar_density": TargetDefinition(
            True,
            "vapor_molar_density_mol_per_m3",
            "mole / meter**3",
            ("molar",),
            ("saturated_vapor",),
        ),
    }
)

_OBSERVATIONS = MappingProxyType(
    {
        "fixed_state_pressure": ObservationDefinition(
            ("pressure",),
            ("temperature_k", "molar_density_mol_per_m3"),
            ("pure", "mixture"),
        ),
        "pure_saturation": ObservationDefinition(
            tuple(target for target in _TARGETS if target != "pressure"),
            ("temperature_k",),
            ("pure",),
        ),
    }
)

_CORRELATIONS = MappingProxyType(
    {
        "power_series_temperature": CorrelationDefinition("power_series", None, ("temperature_k",)),
    }
)

SCHEMA_REGISTRY = _Registry(
    fields=MappingProxyType(
        {
            "fit": _Fields(
                ("schema_version", "parameters", "model", "active", "observations", "solver")
            ),
            "parameters": _Fields(("path", "sha256")),
            "model": _Fields(("kind",)),
            "active": _Fields(("identities", "coordinates", "ties")),
            "coordinate": _Fields(("bounds", "scale"), ("unit",)),
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
            "solver": _Fields(("multistart",)),
            "multistart": _Fields(("count", "include_input")),
        }
    ),
    targets=_TARGETS,
    observations=_OBSERVATIONS,
    residuals=MappingProxyType({"log_ratio": ResidualDefinition(True)}),
    correlations=_CORRELATIONS,
    bases=frozenset({base for target in _TARGETS.values() for base in target.bases}),
    references=frozenset(
        {reference for target in _TARGETS.values() for reference in target.references}
    ),
    transformations=frozenset({"identity", "log", "log10"}),
    model_kinds=("pure", "mixture"),
    source_kinds=("table", "published_observation_correlation"),
    quadrature_orders=(2, 4, 8, 16, 32, 64),
)


def _canonical(value: object) -> bytes:
    return json.dumps(value, ensure_ascii=False, separators=(",", ":"), sort_keys=True).encode(
        "utf-8"
    )


def _digest(value: object) -> str:
    return f"sha256:{hashlib.sha256(_canonical(value)).hexdigest()}"


def _mapping(value: object, path: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping) or any(type(key) is not str for key in value):
        raise FitSchemaError("invalid-type", path, "object with string keys required")
    return cast(Mapping[str, object], value)


def _keys(value: Mapping[str, object], name: str, path: str) -> None:
    fields = SCHEMA_REGISTRY.fields[name]
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
class MultistartSpec:
    count: int
    include_input: bool


@dataclass(frozen=True, slots=True)
class SolverSpec:
    multistart: MultistartSpec


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
    observations: tuple[ObservationSpec, ...]
    solver: SolverSpec

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
            _parse_observation(item, f"$.observations[{index}]")
            for index, item in enumerate(_list(root["observations"], "$.observations"))
        )
        if not observations:
            raise FitSchemaError("missing-value", "$.observations", "at least one required")
        for index, observation in enumerate(observations):
            if model.kind not in SCHEMA_REGISTRY.observations[observation.family].model_kinds:
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
        return {
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
            "observations": [_observation_mapping(item) for item in self.observations],
            "solver": _solver_mapping(self.solver),
        }

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
    kind = _enum(mapping["kind"], SCHEMA_REGISTRY.model_kinds, f"{path}.kind")
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
    return FitCoordinate((lower, upper), scale, unit)


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
    definition: ObservationDefinition,
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
        if SCHEMA_REGISTRY.targets[target].positive and row.value(target) <= 0:
            raise FitSchemaError("nonpositive-target", f"{path}.{target}", row.value(target))
    return row


def _parse_table_source(
    value: object,
    path: str,
    definition: ObservationDefinition,
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
    definition: ObservationDefinition,
    targets: tuple[str, ...],
) -> PublishedObservationCorrelation:
    mapping = _mapping(value, path)
    _keys(mapping, "correlation_source", path)
    name = _enum(
        mapping["correlation"], tuple(sorted(SCHEMA_REGISTRY.correlations)), f"{path}.correlation"
    )
    correlation = SCHEMA_REGISTRY.correlations[name]
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
    definition: ObservationDefinition,
    targets: tuple[str, ...],
) -> Source:
    mapping = _mapping(value, path)
    kind = _enum(mapping.get("kind"), SCHEMA_REGISTRY.source_kinds, f"{path}.kind")
    if kind == SCHEMA_REGISTRY.source_kinds[0]:
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
    allowed = tuple(sorted(SCHEMA_REGISTRY.transformations))
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


def _parse_observation(value: object, path: str) -> ObservationSpec:
    mapping = _mapping(value, path)
    _keys(mapping, "observation", path)
    family = _enum(mapping["family"], tuple(sorted(SCHEMA_REGISTRY.observations)), f"{path}.family")
    definition = SCHEMA_REGISTRY.observations[family]
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
        tuple(sorted(SCHEMA_REGISTRY.residuals)),
    )
    units = _parse_units(mapping["units"], f"{path}.units", targets)
    basis = _parse_target_map(
        mapping["basis"], f"{path}.basis", targets, tuple(sorted(SCHEMA_REGISTRY.bases))
    )
    references = _parse_target_map(
        mapping["reference_convention"],
        f"{path}.reference_convention",
        targets,
        tuple(sorted(SCHEMA_REGISTRY.references)),
    )
    for target, declared in basis:
        if declared not in SCHEMA_REGISTRY.targets[target].bases:
            raise FitSchemaError("invalid-basis", f"{path}.basis.{target}", declared)
    for target, declared in references:
        if declared not in SCHEMA_REGISTRY.targets[target].references:
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
        if (
            initial not in SCHEMA_REGISTRY.quadrature_orders
            or maximum not in SCHEMA_REGISTRY.quadrature_orders
            or tolerance <= 0
        ):
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
    count = _integer(item["count"], f"{path}.multistart.count")
    if count < 1:
        raise FitSchemaError("invalid-value", f"{path}.multistart.count", "positive required")
    return SolverSpec(
        MultistartSpec(
            count,
            _boolean(item["include_input"], f"{path}.multistart.include_input"),
        ),
    )


def _coordinate_mapping(coordinate: FitCoordinate) -> dict[str, object]:
    result: dict[str, object] = {"bounds": list(coordinate.bounds), "scale": coordinate.scale}
    if coordinate.unit is not None:
        result["unit"] = coordinate.unit
    return result


def _domain_mapping(domain: SourceDomain) -> dict[str, object]:
    return {name: list(bounds) for name, bounds in domain.bounds}


def _source_mapping(source: Source) -> dict[str, object]:
    if isinstance(source, TableSource):
        result: dict[str, object] = {
            "kind": SCHEMA_REGISTRY.source_kinds[0],
            "path": source.path,
            "sha256": source.sha256,
            "citation": source.citation,
            "domain": _domain_mapping(source.domain),
        }
        if source.rows is not None:
            result["rows"] = [{name: value for name, value in row.values} for row in source.rows]
        return result
    return {
        "kind": SCHEMA_REGISTRY.source_kinds[1],
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
            "count": solver.multistart.count,
            "include_input": solver.multistart.include_input,
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
            "kind": {"const": SCHEMA_REGISTRY.source_kinds[1]},
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
    definition: ObservationDefinition,
    targets: Sequence[str],
    number: dict[str, object],
) -> dict[str, object]:
    row_fields = (*definition.state_fields, *targets)
    return _object_schema(
        {
            "kind": {"const": SCHEMA_REGISTRY.source_kinds[0]},
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
        SCHEMA_REGISTRY.fields["table_source"].required,
    )


def _quadrature_schema() -> dict[str, object]:
    return _object_schema(
        {
            "initial_order": {"enum": list(SCHEMA_REGISTRY.quadrature_orders)},
            "max_order": {"enum": list(SCHEMA_REGISTRY.quadrature_orders)},
            "tolerance": {"type": "number", "exclusiveMinimum": 0},
        },
        SCHEMA_REGISTRY.fields["quadrature"].required,
    )


def _observation_schema(
    family: str,
    definition: ObservationDefinition,
    targets: tuple[str, ...],
    source: dict[str, object],
    *,
    correlated: bool,
    strings: dict[str, object],
    transformation: dict[str, object],
) -> dict[str, object]:
    residual = {target: {"enum": sorted(SCHEMA_REGISTRY.residuals)} for target in targets}
    unit = {target: strings for target in targets}
    basis = {target: {"enum": list(SCHEMA_REGISTRY.targets[target].bases)} for target in targets}
    reference = {
        target: {"enum": list(SCHEMA_REGISTRY.targets[target].references)} for target in targets
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
    required = list(SCHEMA_REGISTRY.fields["observation"].required)
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
        SCHEMA_REGISTRY.fields["transformation"].required,
    )
    for family, definition in SCHEMA_REGISTRY.observations.items():
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
            for name, correlation in SCHEMA_REGISTRY.correlations.items():
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
            "scale": {"type": "number", "exclusiveMinimum": 0},
            "unit": strings,
        },
        SCHEMA_REGISTRY.fields["coordinate"].required,
    )
    transformation = _object_schema(
        {
            "input": {"enum": sorted(SCHEMA_REGISTRY.transformations)},
            "output": {"enum": sorted(SCHEMA_REGISTRY.transformations)},
        },
        SCHEMA_REGISTRY.fields["transformation"].required,
    )
    observation_schemas = _observation_schemas(strings, number, transformation)
    valid_model_kinds = sorted(
        {
            kind
            for definition in SCHEMA_REGISTRY.observations.values()
            for kind in definition.model_kinds
        }
    )
    return {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "$id": "epcsaft.parameter-fit.v1",
        **_object_schema(
            {
                "schema_version": {"const": SCHEMA_VERSION, "type": "integer"},
                "parameters": _object_schema(
                    {"path": strings, "sha256": digest},
                    SCHEMA_REGISTRY.fields["parameters"].required,
                ),
                "model": _object_schema(
                    {"kind": {"enum": valid_model_kinds}},
                    SCHEMA_REGISTRY.fields["model"].required,
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
                                SCHEMA_REGISTRY.fields["tie"].required,
                            ),
                            "uniqueItems": True,
                        },
                    },
                    SCHEMA_REGISTRY.fields["active"].required,
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
                                "count": {
                                    "type": "integer",
                                    "minimum": 1,
                                    "maximum": _INTEGER_MAX,
                                },
                                "include_input": {"type": "boolean"},
                            },
                            SCHEMA_REGISTRY.fields["multistart"].required,
                        ),
                    },
                    SCHEMA_REGISTRY.fields["solver"].required,
                ),
            },
            SCHEMA_REGISTRY.fields["fit"].required,
        ),
    }
