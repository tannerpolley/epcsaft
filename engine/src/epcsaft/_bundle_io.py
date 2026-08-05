from __future__ import annotations

import csv
import json
import tomllib
from decimal import Decimal, InvalidOperation
from os import PathLike
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any, Literal, Protocol, TypedDict, cast

import pint

from .records import (
    _DIMENSIONAL_UNITS,
    AssociationParameterRecord,
    ComponentRecord,
    ConstantCorrelation,
    ConstantPlusSumOfExponentialsCorrelation,
    CorrelationFamily,
    CorrelationRecord,
    ExponentialTerm,
    ModelParameterRecord,
    PairParameterRecord,
    ParameterRecord,
    SingleParameterRecord,
    SiteRecord,
    SourceRecord,
    ValidityDomain,
)
from .units import ParameterError, _error, unit_registry


class _LoadedBundle(TypedDict):
    bundle_id: str
    bundle_version: int
    purpose: Literal[
        "external-bundle", "package-test-fixture", "regression-output", "user-provided"
    ]
    sources: tuple[SourceRecord, ...]
    domains: tuple[ValidityDomain, ...]
    components: tuple[ComponentRecord, ...]
    singles: tuple[SingleParameterRecord, ...]
    pairs: tuple[PairParameterRecord, ...]
    sites: tuple[SiteRecord, ...]
    associations: tuple[AssociationParameterRecord, ...]
    correlations: tuple[ConstantCorrelation | ConstantPlusSumOfExponentialsCorrelation, ...]
    models: tuple[ModelParameterRecord, ...]


class _WritableBundle(Protocol):
    @property
    def bundle_id(self) -> str: ...

    @property
    def bundle_version(self) -> int: ...

    @property
    def purpose(self) -> str: ...

    @property
    def sources(self) -> tuple[SourceRecord, ...]: ...

    @property
    def domains(self) -> tuple[ValidityDomain, ...]: ...

    @property
    def components(self) -> tuple[ComponentRecord, ...]: ...

    @property
    def records(self) -> tuple[ParameterRecord, ...]: ...

_OPTIONAL_FILES = frozenset(
    {
        "single.csv",
        "pair.csv",
        "sites.csv",
        "association.csv",
        "model.csv",
        "correlations.toml",
    }
)
_COLUMNS = {
    "components.csv": ("component_id", "name", "cas", "formula", "aliases"),
    "single.csv": (
        "record_id",
        "component_id",
        "family",
        "value",
        "unit",
        "source_id",
        "locator",
        "domain_id",
    ),
    "pair.csv": (
        "record_id",
        "component_id_a",
        "component_id_b",
        "family",
        "value",
        "unit",
        "source_id",
        "locator",
        "domain_id",
    ),
    "sites.csv": (
        "record_id",
        "component_id",
        "site_id",
        "site_class",
        "multiplicity",
        "source_id",
        "locator",
        "domain_id",
    ),
    "association.csv": (
        "record_id",
        "component_id_a",
        "site_id_a",
        "component_id_b",
        "site_id_b",
        "family",
        "value",
        "unit",
        "source_id",
        "locator",
        "domain_id",
    ),
    "model.csv": (
        "record_id",
        "family",
        "value",
        "unit",
        "source_id",
        "locator",
        "domain_id",
    ),
}

_TEMPLATE_FILES = {
    "neutral": ("components.csv", "model.csv", "single.csv"),
    "associating": (
        "association.csv",
        "components.csv",
        "model.csv",
        "single.csv",
        "sites.csv",
    ),
    "electrolyte": (
        "components.csv",
        "model.csv",
        "pair.csv",
        "single.csv",
    ),
    "associating-electrolyte": (
        "association.csv",
        "components.csv",
        "model.csv",
        "pair.csv",
        "single.csv",
        "sites.csv",
    ),
}


def _exact_keys(
    value: object, required: set[str], optional: set[str] = set()
) -> dict[str, object]:
    if (
        not isinstance(value, dict)
        or not required <= value.keys()
        or not value.keys() <= required | optional
    ):
        raise _error("invalid-schema", field="keys", value="unexpected-or-missing")
    return value


def _toml(path: Path) -> dict[str, object]:
    try:
        return tomllib.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, tomllib.TOMLDecodeError) as exc:
        raise _error("invalid-schema", field="file", value=path.name) from exc


def _decimal(value: object, *, field: str) -> Decimal:
    try:
        result = Decimal(str(value))
    except (InvalidOperation, ValueError) as exc:
        raise _error("invalid-value", field=field, value=value) from exc
    if not result.is_finite():
        raise _error("invalid-value", field=field, value=value)
    return result


def _quantity(value: object, unit: object, *, field: str) -> pint.Quantity[Any]:
    if not isinstance(unit, str) or not unit:
        raise _error("invalid-unit", field=field, expected="declared-unit")
    try:
        return unit_registry.Quantity(_decimal(value, field=field), unit)
    except Exception as exc:
        if isinstance(exc, ParameterError):
            raise
        raise _error("invalid-unit", field=field, expected=unit) from exc


def _domain(value: object) -> ValidityDomain:
    required = {"domain_id", "kind"}
    optional = {
        "temperature_min",
        "temperature_max",
        "temperature_unit",
        "pressure_min",
        "pressure_max",
        "pressure_unit",
        "ion_mole_fraction_max",
    }
    data = _exact_keys(value, required, optional)
    kwargs: dict[str, object] = {"domain_id": data["domain_id"], "kind": data["kind"]}
    for prefix in ("temperature", "pressure"):
        bounds = [name for name in (f"{prefix}_min", f"{prefix}_max") if name in data]
        unit_name = f"{prefix}_unit"
        if bounds and unit_name not in data:
            raise _error("invalid-unit", field=unit_name, expected="declared-unit")
        if not bounds and unit_name in data:
            raise _error("invalid-schema", field=unit_name, value="unused")
        for name in bounds:
            kwargs[name] = _quantity(data[name], data[unit_name], field=name)
    if "ion_mole_fraction_max" in data:
        kwargs["ion_mole_fraction_max"] = _decimal(
            data["ion_mole_fraction_max"], field="ion_mole_fraction_max"
        )
    return ValidityDomain(**kwargs)  # type: ignore[arg-type]


def _rows(path: Path, name: str) -> list[dict[str, str]]:
    try:
        stream = path.open("r", encoding="utf-8", newline="")
    except OSError as exc:
        raise _error("missing-record", field="file", value=name) from exc
    with stream:
        reader = csv.DictReader(stream)
        if tuple(reader.fieldnames or ()) != _COLUMNS[name]:
            raise _error("invalid-schema", field="columns", value=name)
        rows: list[dict[str, str]] = []
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise _error("invalid-schema", field="row", value=name)
            if not any(row.values()):
                continue
            rows.append(row)
    if not rows:
        raise _error("invalid-schema", field="file", value=f"empty:{name}")
    return rows


def _aliases(value: str) -> tuple[str, ...]:
    if value == "":
        return ()
    try:
        aliases = json.loads(value)
    except json.JSONDecodeError as exc:
        raise _error("invalid-value", field="aliases", value=value) from exc
    if not isinstance(aliases, list) or not all(isinstance(alias, str) for alias in aliases):
        raise _error("invalid-value", field="aliases", value=value)
    return tuple(aliases)


def _parameter_value(row: dict[str, str]) -> object:
    if row["value"] == "":
        raise _error("missing-record", field="value", value=row["record_id"])
    family = row["family"]
    if family in {
        "electrolyte_formulation",
        "relative_permittivity_formulation",
        "polar_formulation",
    }:
        if row["unit"]:
            raise _error("invalid-unit", field="unit", expected="blank")
        return row["value"]
    if family == "charge_number":
        try:
            value = int(row["value"])
        except ValueError as exc:
            raise _error("invalid-value", field="value", value=row["value"]) from exc
        if str(value) != row["value"] and f"+{value}" != row["value"]:
            raise _error("invalid-value", field="value", value=row["value"])
        if row["unit"]:
            raise _error("invalid-unit", field="unit", expected="blank")
        return value
    if family in _DIMENSIONAL_UNITS:
        return _quantity(row["value"], row["unit"], field="value")
    if row["unit"]:
        raise _error("invalid-unit", field="unit", expected="blank")
    return _decimal(row["value"], field="value")


def _components(path: Path) -> tuple[ComponentRecord, ...]:
    return tuple(
        ComponentRecord(
            component_id=row["component_id"],
            name=row["name"] or None,
            cas=row["cas"] or None,
            formula=row["formula"] or None,
            aliases=_aliases(row["aliases"]),
        )
        for row in _rows(path, "components.csv")
    )


def _parameters(
    path: Path, name: str, record_type: type[ParameterRecord]
) -> tuple[ParameterRecord, ...]:
    result = []
    for row in _rows(path, name):
        data: dict[str, object] = {
            key: value for key, value in row.items() if key not in {"unit", "value"}
        }
        data["value"] = _parameter_value(row)
        result.append(cast(ParameterRecord, cast(Any, record_type)(**data)))
    return tuple(result)


def _sites(path: Path) -> tuple[SiteRecord, ...]:
    result = []
    for row in _rows(path, "sites.csv"):
        try:
            multiplicity = int(row["multiplicity"])
        except ValueError as exc:
            raise _error("invalid-value", field="multiplicity", value=row["multiplicity"]) from exc
        result.append(
            SiteRecord(
                record_id=row["record_id"],
                component_id=row["component_id"],
                site_id=row["site_id"],
                site_class=row["site_class"],
                multiplicity=multiplicity,
                source_id=row["source_id"],
                locator=row["locator"],
                domain_id=row["domain_id"],
            )
        )
    return tuple(result)


def _correlations(
    path: Path,
) -> tuple[ConstantCorrelation | ConstantPlusSumOfExponentialsCorrelation, ...]:
    document = _toml(path)
    _exact_keys(document, {"correlations"})
    if not isinstance(document["correlations"], list) or not document["correlations"]:
        raise _error("invalid-schema", field="correlations", value="empty")
    result: list[ConstantCorrelation | ConstantPlusSumOfExponentialsCorrelation] = []
    common = {
        "record_id",
        "component_id",
        "family",
        "form",
        "constant",
        "output_unit",
        "source_id",
        "locator",
        "domain_id",
    }
    for raw in document["correlations"]:
        if not isinstance(raw, dict):
            raise _error("invalid-schema", field="correlation", value="not-table")
        form = raw.get("form")
        if form == "constant":
            data = _exact_keys(raw, common)
            constant = _correlation_value(data["family"], data["constant"], data["output_unit"])
            result.append(
                ConstantCorrelation(
                    record_id=cast(str, data["record_id"]),
                    component_id=cast(str, data["component_id"]),
                    family=cast(CorrelationFamily, data["family"]),
                    constant=constant,
                    source_id=cast(str, data["source_id"]),
                    locator=cast(str, data["locator"]),
                    domain_id=cast(str, data["domain_id"]),
                )
            )
        elif form == "constant-plus-sum-of-exponentials":
            data = _exact_keys(raw, common | {"coefficients"}, {"input_unit"})
            input_unit = data.get("input_unit", "kelvin")
            coefficients = data["coefficients"]
            if not isinstance(coefficients, list) or not coefficients:
                raise _error("invalid-correlation", field="coefficients", value="empty")
            terms = []
            for raw_term in coefficients:
                term = _exact_keys(raw_term, {"amplitude", "exponent_coefficient"})
                terms.append(
                    ExponentialTerm(
                        _quantity(term["amplitude"], data["output_unit"], field="amplitude"),
                        _quantity(
                            term["exponent_coefficient"],
                            f"1 / ({input_unit})",
                            field="exponent_coefficient",
                        ),
                    )
                )
            result.append(
                ConstantPlusSumOfExponentialsCorrelation(
                    record_id=cast(str, data["record_id"]),
                    component_id=cast(str, data["component_id"]),
                    family="segment_diameter",
                    constant=_correlation_value(
                        data["family"], data["constant"], data["output_unit"]
                    ),
                    terms=tuple(terms),
                    source_id=cast(str, data["source_id"]),
                    locator=cast(str, data["locator"]),
                    domain_id=cast(str, data["domain_id"]),
                )
            )
        else:
            raise _error("invalid-correlation", field="form", value=form)
    return tuple(result)


def _correlation_value(family: object, value: object, output_unit: object) -> object:
    if family in _DIMENSIONAL_UNITS:
        return _quantity(value, output_unit, field="constant")
    if output_unit not in {"", "dimensionless"}:
        raise _error("invalid-unit", field="output_unit", expected="dimensionless")
    return _decimal(value, field="constant")


def _load_path(path: str | PathLike[str]) -> _LoadedBundle:
    try:
        root = Path(path)
    except TypeError as exc:
        raise _error("invalid-value", field="path", value=path) from exc
    if not root.is_dir():
        raise _error("missing-record", field="path", value=root)
    manifest_path = root / "bundle.toml"
    if not manifest_path.is_file():
        raise _error("missing-record", field="file", value="bundle.toml")
    document = _toml(manifest_path)
    data = _exact_keys(
        document,
        {
            "schema_version",
            "bundle_id",
            "bundle_version",
            "purpose",
            "files",
            "sources",
            "domains",
        },
    )
    if type(data["schema_version"]) is not int or data["schema_version"] != 1:
        raise _error("invalid-schema", field="schema_version", value=data["schema_version"])
    files = data["files"]
    allowed = _OPTIONAL_FILES | {"components.csv"}
    if (
        not isinstance(files, list)
        or not all(isinstance(name, str) for name in files)
        or len(files) != len(set(files))
        or "components.csv" not in files
        or not set(files) <= allowed
    ):
        raise _error("invalid-schema", field="files", value="invalid-inventory")
    actual = {entry.name for entry in root.iterdir()}
    if actual != {"bundle.toml", *files}:
        raise _error("invalid-schema", field="files", value="inventory-mismatch")

    if not isinstance(data["sources"], list) or not isinstance(data["domains"], list):
        raise _error("invalid-schema", field="metadata", value="not-array")
    sources_list: list[SourceRecord] = []
    for raw in data["sources"]:
        source = _exact_keys(raw, {"source_id", "citation", "use_basis"}, {"doi"})
        sources_list.append(
            SourceRecord(
                source_id=cast(str, source["source_id"]),
                citation=cast(str, source["citation"]),
                use_basis=cast(str, source["use_basis"]),
                doi=cast(str | None, source.get("doi")),
            )
        )
    sources = tuple(sources_list)
    domains = tuple(_domain(raw) for raw in data["domains"])
    result: _LoadedBundle = {
        "bundle_id": cast(str, data["bundle_id"]),
        "bundle_version": cast(int, data["bundle_version"]),
        "purpose": cast(
            Literal[
                "external-bundle",
                "package-test-fixture",
                "regression-output",
                "user-provided",
            ],
            data["purpose"],
        ),
        "sources": sources,
        "domains": domains,
        "components": _components(root / "components.csv"),
        "singles": (),
        "pairs": (),
        "sites": (),
        "associations": (),
        "correlations": (),
        "models": (),
    }
    if "single.csv" in files:
        result["singles"] = cast(
            tuple[SingleParameterRecord, ...],
            _parameters(root / "single.csv", "single.csv", SingleParameterRecord),
        )
    if "pair.csv" in files:
        result["pairs"] = cast(
            tuple[PairParameterRecord, ...],
            _parameters(root / "pair.csv", "pair.csv", PairParameterRecord),
        )
    if "sites.csv" in files:
        result["sites"] = _sites(root / "sites.csv")
    if "association.csv" in files:
        result["associations"] = cast(
            tuple[AssociationParameterRecord, ...],
            _parameters(
                root / "association.csv", "association.csv", AssociationParameterRecord
            ),
        )
    if "model.csv" in files:
        result["models"] = cast(
            tuple[ModelParameterRecord, ...],
            _parameters(root / "model.csv", "model.csv", ModelParameterRecord),
        )
    if "correlations.toml" in files:
        result["correlations"] = _correlations(root / "correlations.toml")
    return result


def _create_template(
    path: str | PathLike[str],
    *,
    components: tuple[str, ...],
    profile: str,
    correlations: bool,
) -> None:
    """Create one non-overwriting authoring skeleton for the strict bundle schema."""
    if profile not in _TEMPLATE_FILES:
        raise _error("invalid-value", field="profile", value=profile)
    try:
        root = Path(path)
    except TypeError as exc:
        raise _error("invalid-value", field="path", value=path) from exc
    if root.exists():
        raise _error("invalid-value", field="path", value="already-exists")
    if not components or len(components) != len(set(components)):
        raise _error("invalid-selection", field="components", value=components)

    files = list(_TEMPLATE_FILES[profile])
    if correlations:
        files.append("correlations.toml")
    files.sort()
    root.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(prefix=f".{root.name}.", dir=root.parent) as temporary:
        staging = Path(temporary) / root.name
        staging.mkdir()
        manifest = (
            "schema_version = 1\n"
            'bundle_id = "replace-with-bundle-id"\n'
            "bundle_version = 1\n"
            'purpose = "user-provided"\n'
            f"files = [{', '.join(_toml_string(name) for name in files)}]\n\n"
            "[[sources]]\n"
            'source_id = "replace-with-source-id"\n'
            'citation = "replace with complete citation"\n'
            'use_basis = "replace with parameter use basis"\n\n'
            "[[domains]]\n"
            'domain_id = "replace-with-domain-id"\n'
            'kind = "replace-with-domain-kind"\n'
        )
        (staging / "bundle.toml").write_text(manifest, encoding="utf-8")
        _write_csv(
            staging,
            "components.csv",
            [(component, "", "", "", "") for component in components],
        )
        for name in files:
            if name == "components.csv":
                continue
            if name == "correlations.toml":
                (staging / name).write_text(
                    "# Add one or more [[correlations]] records.\n",
                    encoding="utf-8",
                )
            else:
                _write_csv(staging, name, [])
        staging.replace(root)


def _number_text(value: object) -> str:
    number = Decimal(str(value))
    if number.is_zero():
        return "0"
    return format(number.normalize(), "f")


def _toml_string(value: str) -> str:
    return json.dumps(value, ensure_ascii=False)


def _parameter_text(family: str, value: object) -> tuple[str, str]:
    if family in {
        "electrolyte_formulation",
        "relative_permittivity_formulation",
        "polar_formulation",
    }:
        assert isinstance(value, str)
        return value, ""
    if family == "charge_number":
        return str(value), ""
    unit = _DIMENSIONAL_UNITS.get(family, "")
    if unit:
        if not isinstance(value, pint.Quantity):
            raise TypeError(f"expected quantity for {family}")
        return _number_text(value.to(unit).magnitude), unit
    return _number_text(value), ""


def _write_manifest(root: Path, bundle: _WritableBundle, files: tuple[str, ...]) -> None:
    lines = [
        "schema_version = 1",
        f"bundle_id = {_toml_string(bundle.bundle_id)}",
        f"bundle_version = {bundle.bundle_version}",
        f"purpose = {_toml_string(bundle.purpose)}",
        f"files = [{', '.join(_toml_string(name) for name in files)}]",
    ]
    sources = bundle.sources
    if not sources:
        lines.append("sources = []")
    domains = bundle.domains
    if not domains:
        lines.append("domains = []")
    for source in sources:
        lines.extend(
            (
                "",
                "[[sources]]",
                f"source_id = {_toml_string(source.source_id)}",
                f"citation = {_toml_string(source.citation)}",
                f"use_basis = {_toml_string(source.use_basis)}",
            )
        )
        if source.doi is not None:
            lines.append(f"doi = {_toml_string(source.doi)}")

    for domain in domains:
        lines.extend(
            (
                "",
                "[[domains]]",
                f"domain_id = {_toml_string(domain.domain_id)}",
                f"kind = {_toml_string(domain.kind)}",
            )
        )
        temperature_bounds = tuple(
            (name, getattr(domain, name))
            for name in ("temperature_min", "temperature_max")
            if getattr(domain, name) is not None
        )
        for name, value in temperature_bounds:
            if not isinstance(value, pint.Quantity):
                raise TypeError("expected temperature quantity")
            magnitude = _number_text(value.to("kelvin").magnitude)
            lines.append(f"{name} = {_toml_string(magnitude)}")
        if temperature_bounds:
            lines.append('temperature_unit = "kelvin"')
        pressure_bounds = tuple(
            (name, getattr(domain, name))
            for name in ("pressure_min", "pressure_max")
            if getattr(domain, name) is not None
        )
        for name, value in pressure_bounds:
            if not isinstance(value, pint.Quantity):
                raise TypeError("expected pressure quantity")
            magnitude = _number_text(value.to("pascal").magnitude)
            lines.append(f"{name} = {_toml_string(magnitude)}")
        if pressure_bounds:
            lines.append('pressure_unit = "pascal"')
        if domain.ion_mole_fraction_max is not None:
            ion_max = _number_text(domain.ion_mole_fraction_max)
            lines.append(f"ion_mole_fraction_max = {_toml_string(ion_max)}")
    (root / "bundle.toml").write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_csv(root: Path, name: str, rows: list[tuple[object, ...]]) -> None:
    with (root / name).open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(_COLUMNS[name])
        writer.writerows(rows)


def _write_components(root: Path, components: tuple[ComponentRecord, ...]) -> None:
    _write_csv(
        root,
        "components.csv",
        [
            (
                record.component_id,
                record.name or "",
                record.cas or "",
                record.formula or "",
                json.dumps(record.aliases, ensure_ascii=False, separators=(",", ":"))
                if record.aliases
                else "",
            )
            for record in components
        ],
    )


def _parameter_row(
    record: SingleParameterRecord
    | PairParameterRecord
    | AssociationParameterRecord
    | ModelParameterRecord,
    fields: tuple[str, ...],
) -> tuple[object, ...]:
    value, unit = _parameter_text(record.family, record.value)
    return (
        record.record_id,
        *(getattr(record, field) for field in fields),
        record.family,
        value,
        unit,
        record.source_id,
        record.locator,
        record.domain_id,
    )


def _converted_magnitude(value: object, unit: str) -> object:
    if not isinstance(value, pint.Quantity):
        raise TypeError(f"expected quantity convertible to {unit}")
    return value.to(unit).magnitude


def _write_correlations(root: Path, records: tuple[CorrelationRecord, ...]) -> None:
    lines: list[str] = []
    for record in records:
        constant, output_unit = _parameter_text(record.family, record.constant)
        if lines:
            lines.append("")
        lines.extend(
            (
                "[[correlations]]",
                f"record_id = {_toml_string(record.record_id)}",
                f"component_id = {_toml_string(record.component_id)}",
                f"family = {_toml_string(record.family)}",
                f"form = {_toml_string(record.form)}",
                f"constant = {_toml_string(constant)}",
                f"output_unit = {_toml_string(output_unit or 'dimensionless')}",
            )
        )
        if isinstance(record, ConstantPlusSumOfExponentialsCorrelation):
            coefficients = ", ".join(
                "{ amplitude = "
                f"{_toml_string(_number_text(_converted_magnitude(term.amplitude, 'angstrom')))}, "
                "exponent_coefficient = "
                f"{_toml_string(_number_text(_converted_magnitude(term.exponent_coefficient, '1 / kelvin')))} }}"
                for term in record.terms
            )
            lines.extend(
                (
                    'input_unit = "kelvin"',
                    f"coefficients = [{coefficients}]",
                )
            )
        lines.extend(
            (
                f"source_id = {_toml_string(record.source_id)}",
                f"locator = {_toml_string(record.locator)}",
                f"domain_id = {_toml_string(record.domain_id)}",
            )
        )
    (root / "correlations.toml").write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_bundle(root: Path, bundle: _WritableBundle) -> None:
    records = bundle.records
    singles = tuple(record for record in records if isinstance(record, SingleParameterRecord))
    pairs = tuple(record for record in records if isinstance(record, PairParameterRecord))
    sites = tuple(record for record in records if isinstance(record, SiteRecord))
    associations = tuple(
        record for record in records if isinstance(record, AssociationParameterRecord)
    )
    models = tuple(record for record in records if isinstance(record, ModelParameterRecord))
    correlations = tuple(
        record
        for record in records
        if isinstance(record, (ConstantCorrelation, ConstantPlusSumOfExponentialsCorrelation))
    )
    groups = (
        ("single.csv", singles),
        ("pair.csv", pairs),
        ("sites.csv", sites),
        ("association.csv", associations),
        ("model.csv", models),
        ("correlations.toml", correlations),
    )
    files = ("components.csv", *(name for name, group in groups if group))
    _write_manifest(root, bundle, files)
    _write_components(root, bundle.components)
    if singles:
        _write_csv(
            root, "single.csv", [_parameter_row(record, ("component_id",)) for record in singles]
        )
    if pairs:
        _write_csv(
            root,
            "pair.csv",
            [_parameter_row(record, ("component_id_a", "component_id_b")) for record in pairs],
        )
    if sites:
        _write_csv(
            root,
            "sites.csv",
            [
                (
                    record.record_id,
                    record.component_id,
                    record.site_id,
                    record.site_class,
                    record.multiplicity,
                    record.source_id,
                    record.locator,
                    record.domain_id,
                )
                for record in sites
            ],
        )
    if associations:
        _write_csv(
            root,
            "association.csv",
            [
                _parameter_row(
                    record,
                    ("component_id_a", "site_id_a", "component_id_b", "site_id_b"),
                )
                for record in associations
            ],
        )
    if models:
        _write_csv(root, "model.csv", [_parameter_row(record, ()) for record in models])
    if correlations:
        _write_correlations(root, correlations)


def _dump_path(bundle: _WritableBundle, path: str | PathLike[str]) -> None:
    try:
        root = Path(path)
    except TypeError as exc:
        raise _error("invalid-value", field="path", value=path) from exc
    if root.exists() or root.is_symlink():
        raise _error("invalid-value", field="path", value=root, reason="already-exists")
    if not root.parent.is_dir():
        raise _error("invalid-value", field="path", value=root, reason="missing-parent")
    if not bundle.components:
        raise _error("invalid-value", field="components", value="empty")

    try:
        with TemporaryDirectory(prefix=f".{root.name}.", dir=root.parent) as temporary:
            staged = Path(temporary)
            _write_bundle(staged, bundle)
            if root.exists() or root.is_symlink():
                raise _error("invalid-value", field="path", value=root, reason="already-exists")
            staged.rename(root)
    except ParameterError:
        raise
    except (OSError, UnicodeError) as exc:
        raise _error("invalid-value", field="path", value=root, reason="write-failed") from exc
