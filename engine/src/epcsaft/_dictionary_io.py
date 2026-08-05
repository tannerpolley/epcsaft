from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from decimal import Decimal
from typing import Any, Literal

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
    _identifier,
)
from .units import _decimal_number, _error, _exact_integer, unit_registry

_MISSING = object()
_TOP_LEVEL = frozenset(
    {"schema", "schema_version", "components", "parameters", "options", "validity"}
)
_ALIASES = {
    "molar_mass": "mw",
    "segment_count": "m",
    "segment_diameter": "s",
    "dispersion_energy_over_k": "e",
    "charge_number": "z",
    "relative_permittivity": "epsilon_r",
    "born_diameter": "d_born",
    "solvation_factor": "f_solv",
    "dielectric_volume": "v_dielectric",
    "dielectric_temperature": "t_dielectric",
    "ion_specific_suppression_coefficient": "a_ion",
    "dipole_scaling": "mu2",
    "polarizability_scaling": "alpha",
    "correlation_integral_parameter": "a2",
    "dipole_moment": "mu",
    "quadrupole_moment": "Q",
    "dispersion_interaction": "k_ij",
    "diameter_interaction": "l_ij",
}
_SINGLE_UNITS = {
    "molar_mass": "kilogram / mole",
    "segment_count": None,
    "segment_diameter": "angstrom",
    "dispersion_energy_over_k": "kelvin",
    "charge_number": None,
    "relative_permittivity": None,
    "born_diameter": "angstrom",
    "solvation_factor": None,
    "dielectric_volume": "meter ** 3 / mole",
    "dielectric_temperature": "kelvin",
    "ion_specific_suppression_coefficient": None,
    "dipole_scaling": "debye ** 2",
    "polarizability_scaling": "angstrom ** 3",
    "correlation_integral_parameter": None,
    "dipole_moment": "debye",
    "quadrupole_moment": "debye * angstrom",
}
_PARAMETER_KEYS = frozenset({*_ALIASES, *_ALIASES.values(), "sites", "association", "correlations"})
_OPTION_ALIASES = {
    "electrolyte_formulation": "electrolyte",
    "relative_permittivity_formulation": "permittivity_model",
    "ionic_region_relative_permittivity": "epsilon_r_ion",
    "ion_fraction_suppression_coefficient": "a_dh",
    "polar_formulation": "polar_model",
}
_OPTION_KEYS = frozenset({*_OPTION_ALIASES, *_OPTION_ALIASES.values()})


@dataclass(frozen=True, slots=True)
class _DictionaryRecords:
    component_ids: tuple[str, ...]
    bundle_id: str
    bundle_version: int
    purpose: Literal["user-provided"]
    sources: tuple[SourceRecord, ...]
    domains: tuple[ValidityDomain, ...]
    components: tuple[ComponentRecord, ...]
    singles: tuple[SingleParameterRecord, ...]
    pairs: tuple[PairParameterRecord, ...]
    sites: tuple[SiteRecord, ...]
    associations: tuple[AssociationParameterRecord, ...]
    correlations: tuple[ConstantCorrelation | ConstantPlusSumOfExponentialsCorrelation, ...]
    models: tuple[ModelParameterRecord, ...]


def _mapping(value: object, *, field: str) -> dict[str, object]:
    if not isinstance(value, Mapping) or any(type(key) is not str for key in value):
        raise _error("invalid-schema", field=field, expected="mapping")
    return dict(value)


def _strict(
    value: dict[str, object],
    *,
    field: str,
    allowed: frozenset[str],
    required: frozenset[str] = frozenset(),
) -> None:
    unknown = sorted(set(value) - allowed)
    if unknown:
        raise _error("invalid-schema", field=field, unknown=unknown[0])
    missing = sorted(required - set(value))
    if missing:
        raise _error("invalid-schema", field=field, missing=missing[0])


def _alias(
    value: dict[str, object],
    descriptive: str,
    abbreviated: str,
    *,
    required: bool = False,
) -> object:
    if descriptive in value and abbreviated in value:
        raise _error(
            "duplicate-alias",
            field=descriptive,
            alias=abbreviated,
        )
    result = value.get(descriptive, value.get(abbreviated, _MISSING))
    if result is _MISSING and required:
        raise _error("invalid-schema", field="parameters", missing=descriptive)
    return result


def _sequence(value: object, *, field: str) -> tuple[object, ...]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise _error("invalid-shape", field=field, expected="sequence")
    return tuple(value)


def _vector(
    value: object,
    *,
    field: str,
    size: int,
    integer: bool = False,
    allow_none: bool = False,
) -> tuple[object | None, ...]:
    result = _sequence(value, field=field)
    if len(result) != size:
        raise _error("invalid-shape", field=field, expected=size, observed=len(result))
    normalized: list[object | None] = []
    for item in result:
        if item is None and allow_none:
            normalized.append(None)
        elif integer:
            normalized.append(_exact_integer(item, field=field))
        else:
            normalized.append(_decimal_number(item, field=field))
    return tuple(normalized)


def _matrix(value: object, *, field: str, size: int) -> tuple[tuple[Decimal, ...], ...]:
    rows = _sequence(value, field=field)
    if len(rows) != size:
        raise _error("invalid-shape", field=field, expected=f"{size}x{size}")
    matrix = tuple(
        tuple(_decimal_number(item, field=field) for item in _sequence(row, field=field))
        for row in rows
    )
    if any(len(row) != size for row in matrix):
        raise _error("invalid-shape", field=field, expected=f"{size}x{size}")
    for i in range(size):
        if matrix[i][i] != 0:
            raise _error("invalid-value", field=field, reason="nonzero-diagonal")
        for j in range(i):
            if matrix[i][j] != matrix[j][i]:
                raise _error("invalid-value", field=field, reason="asymmetric")
    return matrix


def _quantity(value: object, unit: str | None, *, field: str) -> object:
    number = _decimal_number(value, field=field)
    return number if unit is None else unit_registry.Quantity(number, unit)


def _validity(value: object) -> ValidityDomain:
    raw = _mapping(value, field="validity")
    allowed = frozenset(
        {
            "kind",
            "temperature_min_k",
            "temperature_max_k",
            "pressure_min_pa",
            "pressure_max_pa",
            "ion_mole_fraction_max",
        }
    )
    _strict(raw, field="validity", allowed=allowed)
    if "kind" not in raw:
        raise _error("invalid-domain", field="kind", value="missing")
    kind = raw["kind"]
    if kind not in {"unknown", "fit-range", "reported-conditions"}:
        raise _error("invalid-domain", field="kind", value=kind)
    if kind == "unknown" and len(raw) != 1:
        raise _error("invalid-domain", field="kind", value=kind)
    if kind != "unknown" and len(raw) == 1:
        raise _error("invalid-domain", field="bounds", value="missing")
    return ValidityDomain(
        "dictionary-validity",
        kind,
        temperature_min=None
        if "temperature_min_k" not in raw
        else unit_registry.Quantity(
            _decimal_number(raw["temperature_min_k"], field="temperature_min_k"), "kelvin"
        ),
        temperature_max=None
        if "temperature_max_k" not in raw
        else unit_registry.Quantity(
            _decimal_number(raw["temperature_max_k"], field="temperature_max_k"), "kelvin"
        ),
        pressure_min=None
        if "pressure_min_pa" not in raw
        else unit_registry.Quantity(
            _decimal_number(raw["pressure_min_pa"], field="pressure_min_pa"), "pascal"
        ),
        pressure_max=None
        if "pressure_max_pa" not in raw
        else unit_registry.Quantity(
            _decimal_number(raw["pressure_max_pa"], field="pressure_max_pa"), "pascal"
        ),
        ion_mole_fraction_max=raw.get("ion_mole_fraction_max"),
    )


def _single_records(
    raw: dict[str, object], component_ids: tuple[str, ...]
) -> tuple[SingleParameterRecord, ...]:
    records: list[SingleParameterRecord] = []
    for family, unit in _SINGLE_UNITS.items():
        value = _alias(
            raw,
            family,
            _ALIASES[family],
            required=family in {"segment_count", "segment_diameter", "dispersion_energy_over_k"},
        )
        if value is _MISSING:
            continue
        vector = _vector(
            value,
            field=family,
            size=len(component_ids),
            integer=family == "charge_number",
            allow_none=family not in {"segment_count", "dispersion_energy_over_k", "charge_number"},
        )
        for index, (component_id, item) in enumerate(zip(component_ids, vector, strict=True)):
            if item is None:
                continue
            records.append(
                SingleParameterRecord(
                    f"dictionary-{component_id}-{family.replace('_', '-')}",
                    component_id,
                    family,  # type: ignore[arg-type]
                    item if family == "charge_number" else _quantity(item, unit, field=family),
                    "dictionary-input",
                    f"parameters.{family}[{index}]",
                    "dictionary-validity",
                )
            )
    return tuple(records)


def _pair_records(
    raw: dict[str, object], component_ids: tuple[str, ...]
) -> tuple[PairParameterRecord, ...]:
    records: list[PairParameterRecord] = []
    for descriptive, family in (
        ("dispersion_interaction", "k_ij"),
        ("diameter_interaction", "l_ij"),
    ):
        value = _alias(raw, descriptive, family, required=family == "k_ij")
        if value is _MISSING:
            continue
        matrix = _matrix(value, field=descriptive, size=len(component_ids))
        for left in range(len(component_ids)):
            for right in range(left + 1, len(component_ids)):
                records.append(
                    PairParameterRecord(
                        f"dictionary-{family.replace('_', '-')}-{left}-{right}",
                        component_ids[left],
                        component_ids[right],
                        family,  # type: ignore[arg-type]
                        matrix[left][right],
                        "dictionary-input",
                        f"parameters.{descriptive}[{left}][{right}]",
                        "dictionary-validity",
                    )
                )
    return tuple(records)


def _model_records(raw: dict[str, object]) -> tuple[ModelParameterRecord, ...]:
    records = []
    for family, abbreviated in _OPTION_ALIASES.items():
        value = _alias(
            raw,
            family,
            abbreviated,
            required=family == "relative_permittivity_formulation",
        )
        if value is _MISSING:
            continue
        if family in {
            "ionic_region_relative_permittivity",
            "ion_fraction_suppression_coefficient",
        }:
            value = _decimal_number(value, field=family)
        records.append(
            ModelParameterRecord(
                f"dictionary-{family.replace('_', '-')}",
                family,  # type: ignore[arg-type]
                value,
                "dictionary-input",
                f"options.{family}",
                "dictionary-validity",
            )
        )
    return tuple(records)


def _sites(raw: object, component_ids: tuple[str, ...]) -> tuple[SiteRecord, ...]:
    if raw is _MISSING:
        return ()
    records = []
    for index, value in enumerate(_sequence(raw, field="sites")):
        row = _mapping(value, field=f"sites[{index}]")
        _strict(
            row,
            field=f"sites[{index}]",
            allowed=frozenset({"component_id", "site_id", "site_class", "multiplicity"}),
            required=frozenset({"component_id", "site_id", "site_class", "multiplicity"}),
        )
        component_id = _identifier(row["component_id"], field="component_id")
        if component_id not in component_ids:
            raise _error("invalid-value", field="component_id", value=component_id)
        site_id = _identifier(row["site_id"], field="site_id")
        records.append(
            SiteRecord(
                f"dictionary-site-{index}",
                component_id,
                site_id,
                row["site_class"],  # type: ignore[arg-type]
                row["multiplicity"],  # type: ignore[arg-type]
                "dictionary-input",
                f"parameters.sites[{index}]",
                "dictionary-validity",
            )
        )
    return tuple(records)


def _association(raw: object) -> tuple[AssociationParameterRecord, ...]:
    if raw is _MISSING:
        return ()
    records = []
    required = frozenset(
        {
            "component_id_a",
            "site_id_a",
            "component_id_b",
            "site_id_b",
            "association_energy_over_k",
            "association_volume",
        }
    )
    for index, value in enumerate(_sequence(raw, field="association")):
        row = _mapping(value, field=f"association[{index}]")
        _strict(row, field=f"association[{index}]", allowed=required, required=required)
        for family, unit in (
            ("association_energy_over_k", "kelvin"),
            ("association_volume", None),
        ):
            records.append(
                AssociationParameterRecord(
                    f"dictionary-association-{index}-{family.replace('_', '-')}",
                    row["component_id_a"],  # type: ignore[arg-type]
                    row["site_id_a"],  # type: ignore[arg-type]
                    row["component_id_b"],  # type: ignore[arg-type]
                    row["site_id_b"],  # type: ignore[arg-type]
                    family,  # type: ignore[arg-type]
                    _quantity(row[family], unit, field=family),
                    "dictionary-input",
                    f"parameters.association[{index}].{family}",
                    "dictionary-validity",
                )
            )
    return tuple(records)


def _correlations(
    raw: object,
) -> tuple[ConstantCorrelation | ConstantPlusSumOfExponentialsCorrelation, ...]:
    if raw is _MISSING:
        return ()
    records: list[ConstantCorrelation | ConstantPlusSumOfExponentialsCorrelation] = []
    allowed = frozenset({"component_id", "family", "form", "constant", "terms"})
    for index, value in enumerate(_sequence(raw, field="correlations")):
        row = _mapping(value, field=f"correlations[{index}]")
        _strict(
            row,
            field=f"correlations[{index}]",
            allowed=allowed,
            required=frozenset({"component_id", "family", "form", "constant"}),
        )
        family = row["family"]
        if family not in {"segment_diameter", "relative_permittivity", "solvation_factor"}:
            raise _error("invalid-value", field="correlation.family", value=family)
        unit = "angstrom" if family == "segment_diameter" else None
        common: dict[str, Any] = {
            "record_id": f"dictionary-correlation-{index}",
            "component_id": row["component_id"],
            "family": family,
            "constant": _quantity(row["constant"], unit, field="constant"),
            "source_id": "dictionary-input",
            "locator": f"parameters.correlations[{index}]",
            "domain_id": "dictionary-validity",
        }
        if row["form"] == "constant":
            if "terms" in row:
                raise _error("invalid-schema", field="correlation.terms", value="extraneous")
            records.append(ConstantCorrelation(**common))
        elif row["form"] == "constant-plus-sum-of-exponentials":
            if family != "segment_diameter" or "terms" not in row:
                raise _error("invalid-value", field="correlation.form", value=row["form"])
            terms = []
            for term_index, term_value in enumerate(
                _sequence(row["terms"], field="correlation.terms")
            ):
                term = _mapping(term_value, field=f"correlation.terms[{term_index}]")
                _strict(
                    term,
                    field=f"correlation.terms[{term_index}]",
                    allowed=frozenset({"amplitude", "exponent_coefficient_per_k"}),
                    required=frozenset({"amplitude", "exponent_coefficient_per_k"}),
                )
                terms.append(
                    ExponentialTerm(
                        unit_registry.Quantity(
                            _decimal_number(term["amplitude"], field="amplitude"), "angstrom"
                        ),
                        unit_registry.Quantity(
                            _decimal_number(
                                term["exponent_coefficient_per_k"],
                                field="exponent_coefficient_per_k",
                            ),
                            "1 / kelvin",
                        ),
                    )
                )
            records.append(ConstantPlusSumOfExponentialsCorrelation(**common, terms=tuple(terms)))
        else:
            raise _error("invalid-value", field="correlation.form", value=row["form"])
    return tuple(records)


def _records_from_dictionary(value: object) -> _DictionaryRecords:
    raw = _mapping(value, field="dictionary")
    _strict(
        raw,
        field="dictionary",
        allowed=_TOP_LEVEL,
        required=_TOP_LEVEL,
    )
    if raw["schema"] != "epcsaft.parameters" or raw["schema_version"] != 1:
        raise _error("invalid-schema", field="schema", value=raw.get("schema_version"))

    component_ids = tuple(
        _identifier(component, field="component_id")
        for component in _sequence(raw["components"], field="components")
    )
    if not component_ids or len(component_ids) != len(set(component_ids)):
        raise _error("invalid-selection", field="components", value="empty-or-duplicate")

    parameters = _mapping(raw["parameters"], field="parameters")
    options = _mapping(raw["options"], field="options")
    _strict(parameters, field="parameters", allowed=_PARAMETER_KEYS)
    _strict(options, field="options", allowed=_OPTION_KEYS)
    domain = _validity(raw["validity"])
    return _DictionaryRecords(
        component_ids=component_ids,
        bundle_id="dictionary-parameters",
        bundle_version=1,
        purpose="user-provided",
        sources=(
            SourceRecord(
                "dictionary-input",
                "Caller-supplied source-neutral parameter dictionary",
                "Explicit fixed-unit values",
            ),
        ),
        domains=(domain,),
        components=tuple(ComponentRecord(component) for component in component_ids),
        singles=_single_records(parameters, component_ids),
        pairs=_pair_records(parameters, component_ids),
        sites=_sites(parameters.get("sites", _MISSING), component_ids),
        associations=_association(parameters.get("association", _MISSING)),
        correlations=_correlations(parameters.get("correlations", _MISSING)),
        models=_model_records(options),
    )


__all__ = ("_records_from_dictionary",)
