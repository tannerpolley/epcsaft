from __future__ import annotations

import hashlib
import json
import unicodedata
from decimal import Decimal
from typing import cast

import pint

from .records import (
    _DIMENSIONAL_UNITS,
    AssociationParameterRecord,
    ConstantCorrelation,
    ConstantPlusSumOfExponentialsCorrelation,
    CrossAssociationCorrectionRecord,
    ModelParameterRecord,
    PairParameterRecord,
    ParameterRecord,
    SingleParameterRecord,
    SiteRecord,
    ValidityDomain,
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
    if not isinstance(value, pint.Quantity):
        raise TypeError(f"expected quantity for {unit}")
    return {"magnitude": _decimal(value.magnitude), "unit": unit}


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
        result["ion_mole_fraction_max"] = _decimal(record.ion_mole_fraction_max)
    return result


def _resolved_validity(domains: tuple[ValidityDomain, ...]) -> dict[str, object]:
    result: dict[str, object] = {"kinds": sorted({domain.kind for domain in domains})}
    for minimum, maximum, unit in (
        ("temperature_min", "temperature_max", "kelvin"),
        ("pressure_min", "pressure_max", "pascal"),
    ):
        minima = tuple(
            value for domain in domains if (value := getattr(domain, minimum)) is not None
        )
        maxima = tuple(
            value for domain in domains if (value := getattr(domain, maximum)) is not None
        )
        if minima:
            result[minimum] = _quantity(
                max(minima, key=lambda value: value.to(unit).magnitude),
                unit,
            )
        if maxima:
            result[maximum] = _quantity(
                min(maxima, key=lambda value: value.to(unit).magnitude),
                unit,
            )
    ion_maxima = tuple(
        domain.ion_mole_fraction_max
        for domain in domains
        if domain.ion_mole_fraction_max is not None
    )
    if ion_maxima:
        result["ion_mole_fraction_max"] = _decimal(
            min(cast(tuple[Decimal, ...], ion_maxima))
        )
    return result


def _value(family: str, value: object) -> object:
    if family in {
        "electrolyte_formulation",
        "relative_permittivity_formulation",
        "polar_formulation",
    }:
        assert isinstance(value, str)
        return _string(value)
    if family == "charge_number":
        return value
    if family in _DIMENSIONAL_UNITS:
        return _quantity(value, _DIMENSIONAL_UNITS[family])
    return _decimal(value)


def _provenance(record: ParameterRecord) -> dict[str, object]:
    return {
        "record_id": record.record_id,
        "source_id": record.source_id,
        "locator": record.locator,
        "domain_id": record.domain_id,
    }


def _record(record: ParameterRecord) -> dict[str, object]:
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
    elif isinstance(record, CrossAssociationCorrectionRecord):
        result.update(
            type="cross-association-correction",
            component_id_a=record.component_id_a,
            component_id_b=record.component_id_b,
            role_a=record.role_a,
            role_b=record.role_b,
            family=record.family,
            value=_value(record.family, record.value),
        )
    elif isinstance(record, SiteRecord):
        result.update(
            type="site",
            component_id=record.component_id,
            site_id=record.site_id,
            site_role=record.site_role,
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
                    "exponent_coefficient": _quantity(term.exponent_coefficient, "1 / kelvin"),
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


def _fingerprint_json(value: str) -> str:
    return f"sha256:{hashlib.sha256(value.encode('utf-8')).hexdigest()}"


def _document_fingerprint(value: object) -> str:
    return _fingerprint_json(_json({"tag": _TAG, "parameters": value}))


def _selected_document_fingerprint(
    document_fingerprint: str,
    component_ids: tuple[str, ...],
) -> str:
    return _fingerprint_json(
        _json(
            {
                "tag": _TAG,
                "document_fingerprint": document_fingerprint,
                "component_ids": list(component_ids),
            }
        )
    )


def _parameter_fingerprint(
    component_ids: tuple[str, ...],
    records: tuple[ParameterRecord, ...],
    domains: tuple[ValidityDomain, ...],
) -> str:
    domain_by_id = {domain.domain_id: domain for domain in domains}
    scientific_records = []
    for record in records:
        projection = _record(record)
        for field in ("record_id", "source_id", "locator", "domain_id"):
            del projection[field]
        if projection["type"] == "correlation" and projection["form"] == "constant":
            projection = {
                "type": "single",
                "component_id": projection["component_id"],
                "family": projection["family"],
                "value": projection["constant"],
            }
        scientific_records.append(projection)
    selected_domains = tuple(
        domain_by_id[domain_id]
        for domain_id in sorted({record.domain_id for record in records})
    )
    return _fingerprint_json(
        _json(
            {
                "tag": _TAG,
                "parameter_set": {
                    "component_ids": list(component_ids),
                    "records": sorted(scientific_records, key=_json),
                    "validity": _resolved_validity(selected_domains),
                },
            }
        )
    )
