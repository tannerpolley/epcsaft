from __future__ import annotations

from decimal import Decimal
import hashlib
import json
import unicodedata

from .records import (
    AssociationParameterRecord,
    ComponentRecord,
    ConstantCorrelation,
    ConstantPlusSumOfExponentialsCorrelation,
    ModelParameterRecord,
    PairParameterRecord,
    SingleParameterRecord,
    SiteRecord,
    SourceRecord,
    ValidityDomain,
    _DIMENSIONAL_UNITS,
)


_TAG = "epcsaft-parameters-v1"


def _string(value: str) -> str:
    return unicodedata.normalize("NFC", value.replace("\r\n", "\n").replace("\r", "\n"))


def _decimal(value: object) -> str:
    number = Decimal(str(value))
    if number.is_zero():
        return "0"
    return format(number.normalize(), "f")


def _quantity(value: object, unit: str) -> dict[str, str]:
    return {"magnitude": _decimal(value.magnitude), "unit": unit}  # type: ignore[attr-defined]


def _optional(target: dict[str, object], name: str, value: object) -> None:
    if value is not None:
        target[name] = value


def _source(record: SourceRecord) -> dict[str, object]:
    result: dict[str, object] = {
        "source_id": record.source_id,
        "citation": record.citation,
        "use_basis": record.use_basis,
    }
    _optional(result, "doi", record.doi)
    return result


def _domain(record: ValidityDomain) -> dict[str, object]:
    result: dict[str, object] = {"domain_id": record.domain_id, "kind": record.kind}
    for name in ("temperature_min", "temperature_max"):
        value = getattr(record, name)
        if value is not None:
            result[name] = _quantity(value, "kelvin")
    for name in ("pressure_min", "pressure_max"):
        value = getattr(record, name)
        if value is not None:
            result[name] = _quantity(value, "pascal")
    if record.ion_mole_fraction_max is not None:
        result["ion_mole_fraction_max"] = _decimal(
            record.ion_mole_fraction_max
        )
    return result


def _component(record: ComponentRecord) -> dict[str, object]:
    result: dict[str, object] = {
        "component_id": record.component_id,
        "aliases": list(record.aliases),
    }
    for name in ("name", "cas", "formula"):
        _optional(result, name, getattr(record, name))
    return result


def _value(family: str, value: object) -> object:
    if family == "charge_number":
        return value
    if family in _DIMENSIONAL_UNITS:
        return _quantity(value, _DIMENSIONAL_UNITS[family])
    return _decimal(value)


def _provenance(record: object) -> dict[str, object]:
    return {
        "record_id": record.record_id,  # type: ignore[attr-defined]
        "source_id": record.source_id,  # type: ignore[attr-defined]
        "locator": record.locator,  # type: ignore[attr-defined]
        "domain_id": record.domain_id,  # type: ignore[attr-defined]
    }


def _record(record: object) -> dict[str, object]:
    result = _provenance(record)
    if isinstance(record, ModelParameterRecord):
        result.update(
            type="model",
            family=record.family,
            value=_value(record.family, record.value),
        )
    elif isinstance(record, SingleParameterRecord):
        result.update(
            type="single",
            component_id=record.component_id,
            family=record.family,
            value=_value(record.family, record.value),
        )
    elif isinstance(record, PairParameterRecord):
        result.update(
            type="pair",
            component_id_a=record.component_id_a,
            component_id_b=record.component_id_b,
            family=record.family,
            value=_value(record.family, record.value),
        )
    elif isinstance(record, SiteRecord):
        result.update(
            type="site",
            component_id=record.component_id,
            site_id=record.site_id,
            site_class=record.site_class,
            multiplicity=record.multiplicity,
        )
    elif isinstance(record, AssociationParameterRecord):
        result.update(
            type="association",
            component_id_a=record.component_id_a,
            site_id_a=record.site_id_a,
            component_id_b=record.component_id_b,
            site_id_b=record.site_id_b,
            family=record.family,
            value=_value(record.family, record.value),
        )
    elif isinstance(record, ConstantCorrelation):
        result.update(
            type="correlation",
            form=record.form,
            component_id=record.component_id,
            family=record.family,
            constant=_value(record.family, record.constant),
        )
    elif isinstance(record, ConstantPlusSumOfExponentialsCorrelation):
        result.update(
            type="correlation",
            form=record.form,
            component_id=record.component_id,
            family=record.family,
            constant=_value(record.family, record.constant),
            coefficients=[
                {
                    "amplitude": _quantity(term.amplitude, "angstrom"),
                    "exponent_coefficient": _quantity(
                        term.exponent_coefficient, "1 / kelvin"
                    ),
                }
                for term in record.terms
            ],
        )
    else:
        raise TypeError(f"unsupported parameter record: {type(record).__name__}")
    return result


def _normalize(value: object) -> object:
    if isinstance(value, str):
        return _string(value)
    if isinstance(value, list):
        return [_normalize(item) for item in value]
    if isinstance(value, dict):
        return {key: _normalize(item) for key, item in value.items()}
    return value


def _json(value: object) -> str:
    return json.dumps(
        _normalize(value),
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    )


def _bundle_json(bundle: object) -> str:
    projection = {
        "tag": _TAG,
        "bundle": {
            "bundle_id": bundle.bundle_id,  # type: ignore[attr-defined]
            "bundle_version": bundle.bundle_version,  # type: ignore[attr-defined]
            "purpose": bundle.purpose,  # type: ignore[attr-defined]
            "sources": [
                _source(record)
                for record in sorted(
                    bundle.sources, key=lambda item: item.source_id  # type: ignore[attr-defined]
                )
            ],
            "domains": [
                _domain(record)
                for record in sorted(
                    bundle.domains, key=lambda item: item.domain_id  # type: ignore[attr-defined]
                )
            ],
            "components": [
                _component(record)
                for record in sorted(
                    bundle.components, key=lambda item: item.component_id  # type: ignore[attr-defined]
                )
            ],
            "records": [
                _record(record)
                for record in sorted(
                    bundle.records, key=lambda item: item.record_id  # type: ignore[attr-defined]
                )
            ],
        },
    }
    return _json(projection)


def _fingerprint_json(value: str) -> str:
    return f"sha256:{hashlib.sha256(value.encode('utf-8')).hexdigest()}"


def _bundle_fingerprint(bundle: object) -> str:
    return _fingerprint_json(_bundle_json(bundle))


def _set_fingerprint(
    bundle_fingerprint: str,
    component_ids: tuple[str, ...],
    record_ids: tuple[str, ...],
) -> str:
    return _fingerprint_json(
        _json(
            {
                "tag": _TAG,
                "parameter_set": {
                    "bundle_fingerprint": bundle_fingerprint,
                    "component_ids": list(component_ids),
                    "record_ids": list(record_ids),
                },
            }
        )
    )
