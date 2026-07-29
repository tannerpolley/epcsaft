from __future__ import annotations

import csv
from decimal import Decimal, InvalidOperation
import json
from os import PathLike
from pathlib import Path
from tempfile import TemporaryDirectory
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
    if family in {"electrolyte_formulation", "relative_permittivity_formulation"}:
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
    if family in {"electrolyte_formulation", "relative_permittivity_formulation"}:
        assert isinstance(value, str)
        return value, ""
    if family == "charge_number":
        return str(value), ""
    unit = _DIMENSIONAL_UNITS.get(family, "")
    if unit:
        return _number_text(value.to(unit).magnitude), unit  # type: ignore[attr-defined]
    return _number_text(value), ""


def _write_manifest(root: Path, bundle: object, files: tuple[str, ...]) -> None:
    lines = [
        "schema_version = 1",
        f"bundle_id = {_toml_string(bundle.bundle_id)}",  # type: ignore[attr-defined]
        f"bundle_version = {bundle.bundle_version}",  # type: ignore[attr-defined]
        f"purpose = {_toml_string(bundle.purpose)}",  # type: ignore[attr-defined]
        f"files = [{', '.join(_toml_string(name) for name in files)}]",
    ]
    sources = bundle.sources  # type: ignore[attr-defined]
    if not sources:
        lines.append("sources = []")
    domains = bundle.domains  # type: ignore[attr-defined]
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
            magnitude = _number_text(value.to("kelvin").magnitude)  # type: ignore[attr-defined]
            lines.append(f"{name} = {_toml_string(magnitude)}")
        if temperature_bounds:
            lines.append('temperature_unit = "kelvin"')
        pressure_bounds = tuple(
            (name, getattr(domain, name))
            for name in ("pressure_min", "pressure_max")
            if getattr(domain, name) is not None
        )
        for name, value in pressure_bounds:
            magnitude = _number_text(value.to("pascal").magnitude)  # type: ignore[attr-defined]
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


def _parameter_row(record: object, fields: tuple[str, ...]) -> tuple[object, ...]:
    value, unit = _parameter_text(record.family, record.value)  # type: ignore[attr-defined]
    return (
        record.record_id,  # type: ignore[attr-defined]
        *(getattr(record, field) for field in fields),
        record.family,  # type: ignore[attr-defined]
        value,
        unit,
        record.source_id,  # type: ignore[attr-defined]
        record.locator,  # type: ignore[attr-defined]
        record.domain_id,  # type: ignore[attr-defined]
    )


def _write_correlations(root: Path, records: tuple[object, ...]) -> None:
    lines: list[str] = []
    for record in records:
        constant, output_unit = _parameter_text(record.family, record.constant)  # type: ignore[attr-defined]
        if lines:
            lines.append("")
        lines.extend(
            (
                "[[correlations]]",
                f"record_id = {_toml_string(record.record_id)}",  # type: ignore[attr-defined]
                f"component_id = {_toml_string(record.component_id)}",  # type: ignore[attr-defined]
                f"family = {_toml_string(record.family)}",  # type: ignore[attr-defined]
                f"form = {_toml_string(record.form)}",  # type: ignore[attr-defined]
                f"constant = {_toml_string(constant)}",
                f"output_unit = {_toml_string(output_unit or 'dimensionless')}",
            )
        )
        if isinstance(record, ConstantPlusSumOfExponentialsCorrelation):
            coefficients = ", ".join(
                "{ amplitude = "
                f"{_toml_string(_number_text(term.amplitude.to('angstrom').magnitude))}, "
                "exponent_coefficient = "
                f"{_toml_string(_number_text(term.exponent_coefficient.to('1 / kelvin').magnitude))} }}"
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
                f"source_id = {_toml_string(record.source_id)}",  # type: ignore[attr-defined]
                f"locator = {_toml_string(record.locator)}",  # type: ignore[attr-defined]
                f"domain_id = {_toml_string(record.domain_id)}",  # type: ignore[attr-defined]
            )
        )
    (root / "correlations.toml").write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_bundle(root: Path, bundle: object) -> None:
    records = bundle.records  # type: ignore[attr-defined]
    groups = {
        "single.csv": tuple(
            record for record in records if isinstance(record, SingleParameterRecord)
        ),
        "pair.csv": tuple(record for record in records if isinstance(record, PairParameterRecord)),
        "sites.csv": tuple(record for record in records if isinstance(record, SiteRecord)),
        "association.csv": tuple(
            record for record in records if isinstance(record, AssociationParameterRecord)
        ),
        "model.csv": tuple(
            record for record in records if isinstance(record, ModelParameterRecord)
        ),
        "correlations.toml": tuple(
            record
            for record in records
            if isinstance(record, (ConstantCorrelation, ConstantPlusSumOfExponentialsCorrelation))
        ),
    }
    files = ("components.csv", *(name for name, group in groups.items() if group))
    _write_manifest(root, bundle, files)
    _write_components(root, bundle.components)  # type: ignore[attr-defined]
    if singles := groups["single.csv"]:
        _write_csv(
            root, "single.csv", [_parameter_row(record, ("component_id",)) for record in singles]
        )
    if pairs := groups["pair.csv"]:
        _write_csv(
            root,
            "pair.csv",
            [_parameter_row(record, ("component_id_a", "component_id_b")) for record in pairs],
        )
    if sites := groups["sites.csv"]:
        _write_csv(
            root,
            "sites.csv",
            [
                (
                    record.record_id,  # type: ignore[attr-defined]
                    record.component_id,  # type: ignore[attr-defined]
                    record.site_id,  # type: ignore[attr-defined]
                    record.site_class,  # type: ignore[attr-defined]
                    record.multiplicity,  # type: ignore[attr-defined]
                    record.source_id,  # type: ignore[attr-defined]
                    record.locator,  # type: ignore[attr-defined]
                    record.domain_id,  # type: ignore[attr-defined]
                )
                for record in sites
            ],
        )
    if associations := groups["association.csv"]:
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
    if models := groups["model.csv"]:
        _write_csv(root, "model.csv", [_parameter_row(record, ()) for record in models])
    if correlations := groups["correlations.toml"]:
        _write_correlations(root, correlations)


def _dump_path(bundle: object, path: str | PathLike[str]) -> None:
    try:
        root = Path(path)
    except TypeError as exc:
        raise _error("invalid-value", field="path", value=path) from exc
    if root.exists() or root.is_symlink():
        raise _error("invalid-value", field="path", value=root, reason="already-exists")
    if not root.parent.is_dir():
        raise _error("invalid-value", field="path", value=root, reason="missing-parent")
    if not bundle.components:  # type: ignore[attr-defined]
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
