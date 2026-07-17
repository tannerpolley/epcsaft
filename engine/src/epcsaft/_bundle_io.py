from __future__ import annotations

import csv
from decimal import Decimal, InvalidOperation
import json
from os import PathLike
from pathlib import Path
import tomllib

from .records import (
    AssociationParameterRecord,
    ComponentRecord,
    ConstantCorrelation,
    ConstantPlusSumOfExponentialsCorrelation,
    ExponentialTerm,
    ModelParameterRecord,
    PairParameterRecord,
    SingleParameterRecord,
    SiteRecord,
    SourceRecord,
    ValidityDomain,
    _DIMENSIONAL_UNITS,
)
from .units import ParameterError, _error, unit_registry


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


def _exact_keys(value: object, required: set[str], optional: set[str] = set()) -> dict:
    if (
        not isinstance(value, dict)
        or not required <= value.keys()
        or not value.keys() <= required | optional
    ):
        raise _error("invalid-schema", field="keys", value="unexpected-or-missing")
    return value


def _toml(path: Path) -> dict:
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


def _quantity(value: object, unit: object, *, field: str):
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


def _parameters(path: Path, name: str, record_type: type) -> tuple:
    result = []
    for row in _rows(path, name):
        data: dict[str, object] = {
            key: value for key, value in row.items() if key not in {"unit", "value"}
        }
        data["value"] = _parameter_value(row)
        result.append(record_type(**data))
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


def _correlations(path: Path) -> tuple:
    document = _toml(path)
    _exact_keys(document, {"correlations"})
    if not isinstance(document["correlations"], list) or not document["correlations"]:
        raise _error("invalid-schema", field="correlations", value="empty")
    result = []
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
                    record_id=data["record_id"],
                    component_id=data["component_id"],
                    family=data["family"],
                    constant=constant,
                    source_id=data["source_id"],
                    locator=data["locator"],
                    domain_id=data["domain_id"],
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
                    record_id=data["record_id"],
                    component_id=data["component_id"],
                    family=data["family"],
                    constant=_correlation_value(
                        data["family"], data["constant"], data["output_unit"]
                    ),
                    terms=tuple(terms),
                    source_id=data["source_id"],
                    locator=data["locator"],
                    domain_id=data["domain_id"],
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


def _load_path(path: str | PathLike[str]) -> dict[str, object]:
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
    sources = tuple(
        SourceRecord(**_exact_keys(raw, {"source_id", "citation", "use_basis"}, {"doi"}))
        for raw in data["sources"]
    )
    domains = tuple(_domain(raw) for raw in data["domains"])
    result: dict[str, object] = {
        "bundle_id": data["bundle_id"],
        "bundle_version": data["bundle_version"],
        "purpose": data["purpose"],
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
        result["singles"] = _parameters(root / "single.csv", "single.csv", SingleParameterRecord)
    if "pair.csv" in files:
        result["pairs"] = _parameters(root / "pair.csv", "pair.csv", PairParameterRecord)
    if "sites.csv" in files:
        result["sites"] = _sites(root / "sites.csv")
    if "association.csv" in files:
        result["associations"] = _parameters(
            root / "association.csv", "association.csv", AssociationParameterRecord
        )
    if "model.csv" in files:
        result["models"] = _parameters(root / "model.csv", "model.csv", ModelParameterRecord)
    if "correlations.toml" in files:
        result["correlations"] = _correlations(root / "correlations.toml")
    return result
