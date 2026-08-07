from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from decimal import Decimal
from typing import Any, Literal, cast

import pint

from .records import (
    _DIMENSIONAL_UNITS,
    AssociationDerivativeSource,
    AssociationEdge,
    AssociationEdgeSource,
    AssociationParameterRecord,
    AssociationPreset,
    AssociationTopology,
    ComponentRecord,
    ConstantCorrelation,
    ConstantPlusSumOfExponentialsCorrelation,
    ExponentialTerm,
    ModelFamilyKind,
    ModelFamilySpec,
    ModelParameterRecord,
    PairParameterRecord,
    ParameterSpec,
    SingleParameterFamily,
    SingleParameterRecord,
    SiteRecord,
    SourceRecord,
    ValidityDomain,
    _identifier,
    _parameter_identity,
)
from .units import ParameterError, _decimal_number, _error, _exact_integer, unit_registry

type _Purpose = Literal[
    "external-document", "package-test-fixture", "regression-output", "user-provided"
]

_TOP_LEVEL = frozenset(
    {
        "schema",
        "schema_version",
        "document_id",
        "document_version",
        "purpose",
        "sources",
        "domains",
        "components",
        "pairs",
        "model_families",
        "model_coefficients",
        "correlations",
        "topology",
    }
)
_FIXED_FAMILIES = frozenset({"molar_mass", "charge_number", "dipole_moment", "quadrupole_moment"})
_COMPONENT_COEFFICIENTS = frozenset(
    {
        "segment_count",
        "segment_diameter",
        "dispersion_energy_over_k",
        "relative_permittivity",
        "born_diameter",
        "solvation_factor",
        "dielectric_volume",
        "dielectric_temperature",
        "ion_specific_suppression_coefficient",
        "dipole_scaling",
        "polarizability_scaling",
        "correlation_integral_parameter",
    }
)
_PAIR_COEFFICIENTS = frozenset({"k_ij", "l_ij"})
_MODEL_COEFFICIENTS = frozenset(
    {"ionic_region_relative_permittivity", "ion_fraction_suppression_coefficient"}
)
_MODEL_FAMILY_IDS: dict[str, str] = {
    "base": "model/base",
    "association": "model/association",
    "electrolyte": "model/electrolyte",
    "permittivity": "model/relative_permittivity",
    "polar": "model/polar",
}
_MODEL_CHOICES = {
    "base": frozenset({"pc-saft"}),
    "association": frozenset({"none", "general-site"}),
    "electrolyte": frozenset({"none", "ssm-ds", "original-born", "fully-dissociated-debye-huckel"}),
    "permittivity": frozenset(
        {
            "none",
            "ion-fraction-suppression",
            "component-permittivity-mixing",
            "density-volume-temperature",
            "ion-specific-suppression",
            "dipolar-density-perturbation",
            "solvent-only",
        }
    ),
    "polar": frozenset({"none", "point-multipole"}),
}
_MODEL_IDENTITIES = {
    "ionic_region_relative_permittivity": (
        "model/relative_permittivity/ionic_region_relative_permittivity"
    ),
    "ion_fraction_suppression_coefficient": (
        "model/relative_permittivity/ion_fraction_suppression_coefficient"
    ),
}
_PRESET_SITES: dict[str, tuple[tuple[str, str], ...]] = {
    "1a": (("a", "effective"),),
    "2b": (("a", "donor"), ("b", "acceptor")),
    "3b": (("a", "donor"), ("b", "acceptor"), ("c", "acceptor")),
    "4c": (
        ("a", "donor"),
        ("b", "donor"),
        ("c", "acceptor"),
        ("d", "acceptor"),
    ),
}


@dataclass(frozen=True, slots=True)
class _DocumentRecords:
    component_ids: tuple[str, ...]
    document_id: str
    document_version: int
    purpose: _Purpose
    sources: tuple[SourceRecord, ...]
    domains: tuple[ValidityDomain, ...]
    components: tuple[ComponentRecord, ...]
    singles: tuple[SingleParameterRecord, ...]
    pairs: tuple[PairParameterRecord, ...]
    sites: tuple[SiteRecord, ...]
    associations: tuple[AssociationParameterRecord, ...]
    correlations: tuple[ConstantCorrelation | ConstantPlusSumOfExponentialsCorrelation, ...]
    models: tuple[ModelParameterRecord, ...]
    model_families: tuple[ModelFamilySpec, ...]
    topology: AssociationTopology
    parameter_specs: tuple[ParameterSpec, ...]


def _mapping(value: object, *, field: str) -> dict[str, object]:
    if not isinstance(value, Mapping) or any(type(key) is not str for key in value):
        raise _error("invalid-schema", field=field, expected="mapping")
    return dict(value)


def _sequence(value: object, *, field: str) -> tuple[object, ...]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise _error("invalid-schema", field=field, expected="array")
    return tuple(value)


def _strict(
    value: dict[str, object],
    *,
    field: str,
    required: frozenset[str],
    optional: frozenset[str] = frozenset(),
) -> None:
    if unknown := sorted(set(value) - required - optional):
        raise _error("invalid-schema", field=field, unknown=unknown[0])
    if missing := sorted(required - set(value)):
        raise _error("invalid-schema", field=field, missing=missing[0])


def _positive_version(value: object, *, field: str) -> int:
    return _exact_integer(value, field=field, positive=True)


def _provenance(value: object, *, field: str) -> tuple[str, str, str]:
    raw = _mapping(value, field=field)
    _strict(
        raw,
        field=field,
        required=frozenset({"source_id", "locator", "domain_id"}),
    )
    source_id = _identifier(raw["source_id"], field=f"{field}.source_id")
    locator = raw["locator"]
    if not isinstance(locator, str) or not locator.strip():
        raise _error("invalid-value", field=f"{field}.locator", value=locator)
    domain_id = _identifier(raw["domain_id"], field=f"{field}.domain_id")
    return source_id, locator, domain_id


def _measurement(
    value: object,
    *,
    family: str,
    field: str,
) -> tuple[object, str]:
    raw = _mapping(value, field=field)
    _strict(raw, field=field, required=frozenset({"magnitude", "unit"}))
    unit = raw["unit"]
    if not isinstance(unit, str) or not unit:
        raise _error("invalid-unit", field=f"{field}.unit", expected="declared-unit")
    magnitude = raw["magnitude"]
    if family == "charge_number":
        if unit != "elementary-charge":
            raise _error("invalid-unit", field=f"{field}.unit", expected="elementary-charge")
        return _exact_integer(magnitude, field=f"{field}.magnitude"), unit
    target = (
        "1 / kelvin"
        if family == "exponential_temperature_coefficient"
        else _DIMENSIONAL_UNITS.get(family, "dimensionless")
    )
    try:
        quantity = unit_registry.Quantity(
            _decimal_number(magnitude, field=f"{field}.magnitude"), unit
        ).to(target)
    except ParameterError:
        raise
    except (AttributeError, TypeError, ValueError, pint.errors.PintError) as exc:
        raise _error("invalid-unit", field=f"{field}.unit", expected=target) from exc
    if target == "dimensionless":
        return (
            _decimal_number(quantity.magnitude, field=f"{field}.magnitude"),
            "dimensionless",
        )
    return quantity, target


def _coefficient(
    value: object,
    *,
    field: str,
    allowed_families: frozenset[str],
    expected_identity: str,
) -> tuple[str, object, ParameterSpec, tuple[str, str, str]]:
    raw = _mapping(value, field=field)
    _strict(
        raw,
        field=field,
        required=frozenset({"identity", "family", "value", "provenance"}),
    )
    family = raw["family"]
    if not isinstance(family, str) or family not in allowed_families:
        raise _error("invalid-value", field=f"{field}.family", value=family)
    identity = _parameter_identity(raw["identity"], field=f"{field}.identity")
    if identity != expected_identity:
        raise _error(
            "invalid-value",
            field=f"{field}.identity",
            expected=expected_identity,
            value=identity,
        )
    parameter_value, unit = _measurement(raw["value"], family=family, field=f"{field}.value")
    provenance = _provenance(raw["provenance"], field=f"{field}.provenance")
    return (
        family,
        parameter_value,
        ParameterSpec(identity, family, parameter_value, unit, True, *provenance),
        provenance,
    )


def _record_id(identity: str) -> str:
    return "parameter-" + identity.replace("/", "-").replace("_", "-")


def _sources(value: object) -> tuple[SourceRecord, ...]:
    result = []
    for index, item in enumerate(_sequence(value, field="sources")):
        field = f"sources[{index}]"
        raw = _mapping(item, field=field)
        _strict(
            raw,
            field=field,
            required=frozenset({"source_id", "citation", "use_basis"}),
            optional=frozenset({"doi"}),
        )
        result.append(
            SourceRecord(
                cast(str, raw["source_id"]),
                cast(str, raw["citation"]),
                cast(str, raw["use_basis"]),
                cast(str | None, raw.get("doi")),
            )
        )
    return tuple(result)


def _domains(value: object) -> tuple[ValidityDomain, ...]:
    result = []
    required = frozenset({"domain_id", "kind"})
    optional = frozenset(
        {
            "temperature_min",
            "temperature_max",
            "pressure_min",
            "pressure_max",
            "ion_mole_fraction_max",
        }
    )
    for index, item in enumerate(_sequence(value, field="domains")):
        field = f"domains[{index}]"
        raw = _mapping(item, field=field)
        _strict(raw, field=field, required=required, optional=optional)
        kwargs: dict[str, object] = {
            "domain_id": raw["domain_id"],
            "kind": raw["kind"],
        }
        for name, family in (
            ("temperature_min", "dielectric_temperature"),
            ("temperature_max", "dielectric_temperature"),
            ("pressure_min", "pressure"),
            ("pressure_max", "pressure"),
        ):
            if name in raw:
                if family == "pressure":
                    measurement = _mapping(raw[name], field=f"{field}.{name}")
                    _strict(
                        measurement,
                        field=f"{field}.{name}",
                        required=frozenset({"magnitude", "unit"}),
                    )
                    unit = measurement["unit"]
                    if not isinstance(unit, str):
                        raise _error("invalid-unit", field=f"{field}.{name}.unit")
                    try:
                        kwargs[name] = unit_registry.Quantity(
                            _decimal_number(
                                measurement["magnitude"], field=f"{field}.{name}.magnitude"
                            ),
                            unit,
                        ).to("pascal")
                    except (TypeError, ValueError, pint.errors.PintError) as exc:
                        raise _error(
                            "invalid-unit", field=f"{field}.{name}.unit", expected="pascal"
                        ) from exc
                else:
                    kwargs[name] = _measurement(raw[name], family=family, field=f"{field}.{name}")[
                        0
                    ]
        if "ion_mole_fraction_max" in raw:
            kwargs["ion_mole_fraction_max"] = _measurement(
                raw["ion_mole_fraction_max"],
                family="solvation_factor",
                field=f"{field}.ion_mole_fraction_max",
            )[0]
        result.append(ValidityDomain(**kwargs))  # type: ignore[arg-type]
    return tuple(result)


def _components(
    value: object,
) -> tuple[
    tuple[ComponentRecord, ...], tuple[SingleParameterRecord, ...], tuple[ParameterSpec, ...]
]:
    components = []
    singles = []
    specs = []
    required = frozenset({"component_id", "fixed", "coefficients"})
    optional = frozenset({"name", "cas", "formula", "aliases"})
    for index, item in enumerate(_sequence(value, field="components")):
        field = f"components[{index}]"
        raw = _mapping(item, field=field)
        _strict(raw, field=field, required=required, optional=optional)
        component_id = _identifier(raw["component_id"], field=f"{field}.component_id")
        aliases = _sequence(raw.get("aliases", ()), field=f"{field}.aliases")
        components.append(
            ComponentRecord(
                component_id,
                cast(str | None, raw.get("name")),
                cast(str | None, raw.get("cas")),
                cast(str | None, raw.get("formula")),
                cast(tuple[str, ...], aliases),
            )
        )

        fixed = _mapping(raw["fixed"], field=f"{field}.fixed")
        if unknown := sorted(set(fixed) - _FIXED_FAMILIES):
            raise _error("invalid-schema", field=f"{field}.fixed", unknown=unknown[0])
        for family, item_value in fixed.items():
            fixed_field = f"{field}.fixed.{family}"
            fixed_value = _mapping(item_value, field=fixed_field)
            _strict(
                fixed_value,
                field=fixed_field,
                required=frozenset({"value", "provenance"}),
            )
            parameter_value, _ = _measurement(
                fixed_value["value"], family=family, field=f"{fixed_field}.value"
            )
            provenance = _provenance(fixed_value["provenance"], field=f"{fixed_field}.provenance")
            singles.append(
                SingleParameterRecord(
                    _record_id(f"component/{component_id}/{family}"),
                    component_id,
                    cast(SingleParameterFamily, family),
                    parameter_value,
                    *provenance,
                )
            )

        for coefficient_index, coefficient_value in enumerate(
            _sequence(raw["coefficients"], field=f"{field}.coefficients")
        ):
            coefficient_field = f"{field}.coefficients[{coefficient_index}]"
            coefficient_raw = _mapping(coefficient_value, field=coefficient_field)
            coefficient_family = coefficient_raw.get("family")
            if not isinstance(coefficient_family, str):
                raise _error(
                    "invalid-value",
                    field=f"{coefficient_field}.family",
                    value=coefficient_family,
                )
            identity = f"component/{component_id}/{coefficient_family}"
            resolved_family, parameter_value, spec, provenance = _coefficient(
                coefficient_raw,
                field=coefficient_field,
                allowed_families=_COMPONENT_COEFFICIENTS,
                expected_identity=identity,
            )
            singles.append(
                SingleParameterRecord(
                    _record_id(identity),
                    component_id,
                    cast(SingleParameterFamily, resolved_family),
                    parameter_value,
                    *provenance,
                )
            )
            specs.append(spec)
    return tuple(components), tuple(singles), tuple(specs)


def _pairs(
    value: object,
    component_ids: frozenset[str],
) -> tuple[tuple[PairParameterRecord, ...], tuple[ParameterSpec, ...]]:
    records = []
    specs = []
    for index, item in enumerate(_sequence(value, field="pairs")):
        field = f"pairs[{index}]"
        raw = _mapping(item, field=field)
        _strict(
            raw,
            field=field,
            required=frozenset({"component_id_a", "component_id_b", "coefficients"}),
        )
        component_a = _identifier(raw["component_id_a"], field=f"{field}.component_id_a")
        component_b = _identifier(raw["component_id_b"], field=f"{field}.component_id_b")
        if component_a not in component_ids or component_b not in component_ids:
            raise _error("missing-record", field=field, value="component")
        if component_a >= component_b:
            raise _error("invalid-value", field=field, reason="noncanonical-pair-order")
        for coefficient_index, coefficient_value in enumerate(
            _sequence(raw["coefficients"], field=f"{field}.coefficients")
        ):
            coefficient_field = f"{field}.coefficients[{coefficient_index}]"
            coefficient_raw = _mapping(coefficient_value, field=coefficient_field)
            family = coefficient_raw.get("family")
            if not isinstance(family, str):
                raise _error("invalid-value", field=f"{coefficient_field}.family", value=family)
            identity = f"pair/{component_a}/{component_b}/{family}"
            family, parameter_value, spec, provenance = _coefficient(
                coefficient_raw,
                field=coefficient_field,
                allowed_families=_PAIR_COEFFICIENTS,
                expected_identity=identity,
            )
            records.append(
                PairParameterRecord(
                    _record_id(identity),
                    component_a,
                    component_b,
                    cast(Any, family),
                    parameter_value,
                    *provenance,
                )
            )
            specs.append(spec)
    return tuple(records), tuple(specs)


def _model_families(value: object) -> tuple[ModelFamilySpec, ...]:
    result = []
    seen: set[str] = set()
    for index, item in enumerate(_sequence(value, field="model_families")):
        field = f"model_families[{index}]"
        raw = _mapping(item, field=field)
        _strict(
            raw,
            field=field,
            required=frozenset({"family_id", "kind", "choice", "provenance"}),
        )
        kind = raw["kind"]
        if not isinstance(kind, str) or kind not in _MODEL_CHOICES:
            raise _error("invalid-value", field=f"{field}.kind", value=kind)
        if kind in seen:
            raise _error("duplicate-record", field="model_family.kind", value=kind)
        seen.add(kind)
        family_id = _parameter_identity(raw["family_id"], field=f"{field}.family_id")
        if family_id != _MODEL_FAMILY_IDS[kind]:
            raise _error("invalid-value", field=f"{field}.family_id", value=family_id)
        choice = raw["choice"]
        if not isinstance(choice, str) or choice not in _MODEL_CHOICES[kind]:
            raise _error("invalid-value", field=f"{field}.choice", value=choice)
        result.append(
            ModelFamilySpec(
                family_id,
                cast(ModelFamilyKind, kind),
                choice,
                *_provenance(raw["provenance"], field=f"{field}.provenance"),
            )
        )
    if seen != set(_MODEL_CHOICES):
        missing = min(set(_MODEL_CHOICES) - seen)
        raise _error("missing-record", field="model_family.kind", value=missing)
    return tuple(sorted(result, key=lambda item: item.kind))


def _model_coefficients(
    value: object,
) -> tuple[tuple[ModelParameterRecord, ...], tuple[ParameterSpec, ...]]:
    records = []
    specs = []
    for index, item in enumerate(_sequence(value, field="model_coefficients")):
        field = f"model_coefficients[{index}]"
        raw = _mapping(item, field=field)
        family = raw.get("family")
        if not isinstance(family, str) or family not in _MODEL_COEFFICIENTS:
            raise _error("invalid-value", field=f"{field}.family", value=family)
        identity = _MODEL_IDENTITIES[family]
        family, parameter_value, spec, provenance = _coefficient(
            raw,
            field=field,
            allowed_families=_MODEL_COEFFICIENTS,
            expected_identity=identity,
        )
        records.append(
            ModelParameterRecord(
                _record_id(identity), cast(Any, family), parameter_value, *provenance
            )
        )
        specs.append(spec)
    return tuple(records), tuple(specs)


def _correlation_coefficient(
    value: object,
    *,
    field: str,
    identity: str,
    family: str,
    provenance: tuple[str, str, str],
) -> tuple[object, ParameterSpec]:
    raw = _mapping(value, field=field)
    _strict(raw, field=field, required=frozenset({"identity", "value"}))
    observed = _parameter_identity(raw["identity"], field=f"{field}.identity")
    if observed != identity:
        raise _error("invalid-value", field=f"{field}.identity", expected=identity, value=observed)
    parameter_value, unit = _measurement(raw["value"], family=family, field=f"{field}.value")
    return parameter_value, ParameterSpec(
        identity, family, parameter_value, unit, True, *provenance
    )


def _correlations(
    value: object,
    component_ids: frozenset[str],
) -> tuple[
    tuple[ConstantCorrelation | ConstantPlusSumOfExponentialsCorrelation, ...],
    tuple[ParameterSpec, ...],
]:
    records: list[ConstantCorrelation | ConstantPlusSumOfExponentialsCorrelation] = []
    specs = []
    for index, item in enumerate(_sequence(value, field="correlations")):
        field = f"correlations[{index}]"
        raw = _mapping(item, field=field)
        _strict(
            raw,
            field=field,
            required=frozenset(
                {
                    "correlation_id",
                    "component_id",
                    "family",
                    "form",
                    "independent_variables",
                    "constant",
                    "provenance",
                }
            ),
            optional=frozenset({"terms"}),
        )
        component_id = _identifier(raw["component_id"], field=f"{field}.component_id")
        if component_id not in component_ids:
            raise _error("missing-record", field=f"{field}.component_id", value=component_id)
        family = raw["family"]
        if family not in {"segment_diameter", "relative_permittivity", "solvation_factor"}:
            raise _error("invalid-value", field=f"{field}.family", value=family)
        form = raw["form"]
        if form not in {"constant", "constant-plus-sum-of-exponentials"}:
            raise _error("invalid-correlation", field=f"{field}.form", value=form)
        if form == "constant-plus-sum-of-exponentials" and family != "segment_diameter":
            raise _error("invalid-correlation", field=f"{field}.form", family=family)
        correlation_id = f"component/{component_id}/{family}/{form}"
        observed_id = _parameter_identity(raw["correlation_id"], field=f"{field}.correlation_id")
        if observed_id != correlation_id:
            raise _error("invalid-value", field=f"{field}.correlation_id", value=observed_id)
        independent = _sequence(
            raw["independent_variables"], field=f"{field}.independent_variables"
        )
        expected_independent = () if form == "constant" else ("temperature",)
        if independent != expected_independent:
            raise _error(
                "invalid-correlation",
                field=f"{field}.independent_variables",
                value=independent,
            )
        provenance = _provenance(raw["provenance"], field=f"{field}.provenance")
        constant, constant_spec = _correlation_coefficient(
            raw["constant"],
            field=f"{field}.constant",
            identity=f"{correlation_id}/constant",
            family=family,
            provenance=provenance,
        )
        specs.append(constant_spec)
        if form == "constant":
            if "terms" in raw:
                raise _error("invalid-schema", field=f"{field}.terms", value="extraneous")
            records.append(
                ConstantCorrelation(
                    _record_id(correlation_id),
                    component_id,
                    family,
                    constant,
                    *provenance,
                    correlation_id=correlation_id,
                )
            )
            continue
        if "terms" not in raw:
            raise _error("invalid-schema", field=field, missing="terms")
        terms = []
        for term_index, item_term in enumerate(_sequence(raw["terms"], field=f"{field}.terms")):
            term_field = f"{field}.terms[{term_index}]"
            term = _mapping(item_term, field=term_field)
            _strict(
                term,
                field=term_field,
                required=frozenset({"amplitude", "exponent_coefficient"}),
            )
            prefix = f"{correlation_id}/term-{term_index}"
            amplitude, amplitude_spec = _correlation_coefficient(
                term["amplitude"],
                field=f"{term_field}.amplitude",
                identity=f"{prefix}/amplitude",
                family="segment_diameter",
                provenance=provenance,
            )
            exponent, exponent_spec = _correlation_coefficient(
                term["exponent_coefficient"],
                field=f"{term_field}.exponent_coefficient",
                identity=f"{prefix}/exponent_coefficient",
                family="exponential_temperature_coefficient",
                provenance=provenance,
            )
            terms.append(
                ExponentialTerm(
                    amplitude,
                    exponent,
                    amplitude_identity=amplitude_spec.identity,
                    exponent_identity=exponent_spec.identity,
                )
            )
            specs.extend((amplitude_spec, exponent_spec))
        records.append(
            ConstantPlusSumOfExponentialsCorrelation(
                _record_id(correlation_id),
                component_id,
                "segment_diameter",
                constant,
                tuple(terms),
                *provenance,
                correlation_id=correlation_id,
            )
        )
    return tuple(records), tuple(specs)


def _endpoint(
    value: object,
    *,
    field: str,
    sites: frozenset[tuple[str, str]],
) -> tuple[str, str]:
    raw = _mapping(value, field=field)
    _strict(raw, field=field, required=frozenset({"component_id", "site_id"}))
    endpoint = (
        _identifier(raw["component_id"], field=f"{field}.component_id"),
        _identifier(raw["site_id"], field=f"{field}.site_id"),
    )
    if endpoint not in sites:
        raise _error("missing-record", field=field, value=endpoint)
    return endpoint


def _provenance_list(value: object, *, field: str) -> tuple[tuple[str, str, str], ...]:
    result = tuple(
        _provenance(item, field=f"{field}[{index}]")
        for index, item in enumerate(_sequence(value, field=field))
    )
    if not result:
        raise _error("invalid-value", field=field, value="empty")
    if len(result) != len(set(result)):
        raise _error("duplicate-record", field=field)
    return tuple(sorted(result))


def _edge_value(
    value: object,
    *,
    field: str,
    identity: str,
    family: str,
    provenance: tuple[str, str, str],
) -> tuple[object, ParameterSpec]:
    raw = _mapping(value, field=field)
    _strict(raw, field=field, required=frozenset({"identity", "value"}))
    observed = _parameter_identity(raw["identity"], field=f"{field}.identity")
    if observed != identity:
        raise _error("invalid-value", field=f"{field}.identity", expected=identity, value=observed)
    parameter_value, unit = _measurement(raw["value"], family=family, field=f"{field}.value")
    return parameter_value, ParameterSpec(
        identity, family, parameter_value, unit, True, *provenance
    )


def _edge_records(
    edge: AssociationEdge,
) -> tuple[AssociationParameterRecord, AssociationParameterRecord]:
    provenance = edge.source.provenance[0]
    common = (
        edge.endpoint_a[0],
        edge.endpoint_a[1],
        edge.endpoint_b[0],
        edge.endpoint_b[1],
    )
    return (
        AssociationParameterRecord(
            _record_id(edge.energy_identity),
            *common,
            "association_energy_over_k",
            edge.energy_over_k,
            *provenance,
        ),
        AssociationParameterRecord(
            _record_id(edge.volume_identity),
            *common,
            "association_volume",
            edge.volume,
            *provenance,
        ),
    )


def _explicit_edge(
    raw: dict[str, object],
    *,
    field: str,
    endpoint_a: tuple[str, str],
    endpoint_b: tuple[str, str],
) -> tuple[AssociationEdge, tuple[ParameterSpec, ParameterSpec]]:
    _strict(
        raw,
        field=field,
        required=frozenset({"endpoint_a", "endpoint_b", "energy_over_k", "volume", "source"}),
    )
    source_raw = _mapping(raw["source"], field=f"{field}.source")
    _strict(
        source_raw,
        field=f"{field}.source",
        required=frozenset({"kind", "provenance"}),
    )
    if source_raw["kind"] != "explicit":
        raise _error("invalid-value", field=f"{field}.source.kind", value=source_raw["kind"])
    provenance = _provenance_list(source_raw["provenance"], field=f"{field}.source.provenance")
    prefix = f"association/{endpoint_a[0]}/{endpoint_a[1]}/{endpoint_b[0]}/{endpoint_b[1]}"
    energy_identity = f"{prefix}/energy_over_k"
    volume_identity = f"{prefix}/volume"
    energy, energy_spec = _edge_value(
        raw["energy_over_k"],
        field=f"{field}.energy_over_k",
        identity=energy_identity,
        family="association_energy_over_k",
        provenance=provenance[0],
    )
    volume, volume_spec = _edge_value(
        raw["volume"],
        field=f"{field}.volume",
        identity=volume_identity,
        family="association_volume",
        provenance=provenance[0],
    )
    return (
        AssociationEdge(
            endpoint_a,
            endpoint_b,
            energy_identity,
            volume_identity,
            energy,
            volume,
            AssociationEdgeSource("explicit", None, provenance, (), ()),
        ),
        (energy_spec, volume_spec),
    )


def _rule_inputs(value: object, *, field: str) -> tuple[tuple[str, str], tuple[str, str]]:
    raw = _mapping(value, field=field)
    _strict(raw, field=field, required=frozenset({"energy_over_k", "volume"}))
    result = []
    for family in ("energy_over_k", "volume"):
        identities = tuple(
            _parameter_identity(item, field=f"{field}.{family}")
            for item in _sequence(raw[family], field=f"{field}.{family}")
        )
        if len(identities) != 2 or len(set(identities)) != 2:
            raise _error("invalid-value", field=f"{field}.{family}", value=identities)
        result.append(tuple(sorted(identities)))
    return cast(tuple[tuple[str, str], tuple[str, str]], tuple(result))


def _association_input_edge(
    identity: str,
    *,
    family: str,
    field: str,
) -> tuple[tuple[str, str], tuple[str, str]]:
    parts = identity.split("/")
    if len(parts) != 6 or parts[0] != "association" or parts[5] != family:
        raise _error("invalid-value", field=field, value=identity)
    endpoints = tuple(sorted(((parts[1], parts[2]), (parts[3], parts[4]))))
    if endpoints[0][0] != endpoints[1][0]:
        raise _error(
            "invalid-value",
            field=field,
            value=identity,
            reason="combining-input-must-be-self-edge",
        )
    return cast(tuple[tuple[str, str], tuple[str, str]], endpoints)


def _validate_combining_inputs(
    *,
    field: str,
    endpoint_a: tuple[str, str],
    endpoint_b: tuple[str, str],
    energy_inputs: tuple[str, str],
    volume_inputs: tuple[str, str],
    site_roles: Mapping[tuple[str, str], str],
) -> None:
    output_endpoints = (endpoint_a, endpoint_b)
    if endpoint_a[0] == endpoint_b[0]:
        raise _error(
            "invalid-value",
            field=field,
            reason="combining-rule-requires-cross-component-edge",
        )
    if {site_roles[endpoint] for endpoint in output_endpoints} != {"donor", "acceptor"}:
        raise _error(
            "invalid-value",
            field=field,
            reason="combining-rule-requires-donor-acceptor-edge",
        )
    expected_components = {endpoint_a[0], endpoint_b[0]}
    input_edge_sets: dict[
        str,
        frozenset[tuple[tuple[str, str], tuple[str, str]]],
    ] = {}
    for family, identities in (
        ("energy_over_k", energy_inputs),
        ("volume", volume_inputs),
    ):
        input_edges = tuple(
            _association_input_edge(
                identity,
                family=family,
                field=f"{field}.source.inputs.{family}",
            )
            for identity in identities
        )
        input_edge_sets[family] = frozenset(input_edges)
        if {edge[0][0] for edge in input_edges} != expected_components:
            raise _error(
                "invalid-value",
                field=f"{field}.source.inputs.{family}",
                reason="combining-input-components-do-not-match-edge",
            )
        for output_endpoint in output_endpoints:
            matching = tuple(edge for edge in input_edges if edge[0][0] == output_endpoint[0])
            if len(matching) != 1 or output_endpoint not in matching[0]:
                raise _error(
                    "invalid-value",
                    field=f"{field}.source.inputs.{family}",
                    reason="combining-input-site-does-not-match-edge",
                )
            if {site_roles[endpoint] for endpoint in matching[0]} != {"donor", "acceptor"}:
                raise _error(
                    "invalid-value",
                    field=f"{field}.source.inputs.{family}",
                    reason="combining-input-must-be-donor-acceptor-self-edge",
                )
    if input_edge_sets["energy_over_k"] != input_edge_sets["volume"]:
        raise _error(
            "invalid-value",
            field=f"{field}.source.inputs",
            reason="combining-input-families-must-use-the-same-self-edges",
        )


def _combining_edge(
    raw: dict[str, object],
    *,
    field: str,
    endpoint_a: tuple[str, str],
    endpoint_b: tuple[str, str],
    explicit_values: Mapping[str, tuple[str, object]],
    site_roles: Mapping[tuple[str, str], str],
) -> tuple[AssociationEdge, tuple[ParameterSpec, ParameterSpec]]:
    _strict(
        raw,
        field=field,
        required=frozenset({"endpoint_a", "endpoint_b", "source"}),
    )
    source_raw = _mapping(raw["source"], field=f"{field}.source")
    _strict(
        source_raw,
        field=f"{field}.source",
        required=frozenset({"kind", "rule_id", "inputs", "provenance"}),
    )
    if source_raw["kind"] != "combining-rule":
        raise _error("invalid-value", field=f"{field}.source.kind", value=source_raw["kind"])
    rule_id = source_raw["rule_id"]
    if rule_id != "arithmetic-energy-geometric-volume":
        raise _error("invalid-value", field=f"{field}.source.rule_id", value=rule_id)
    energy_inputs, volume_inputs = _rule_inputs(
        source_raw["inputs"], field=f"{field}.source.inputs"
    )
    _validate_combining_inputs(
        field=field,
        endpoint_a=endpoint_a,
        endpoint_b=endpoint_b,
        energy_inputs=energy_inputs,
        volume_inputs=volume_inputs,
        site_roles=site_roles,
    )
    for identity in energy_inputs:
        if explicit_values.get(identity, (None,))[0] != "association_energy_over_k":
            raise _error("missing-record", field=f"{field}.source.inputs", value=identity)
    for identity in volume_inputs:
        if explicit_values.get(identity, (None,))[0] != "association_volume":
            raise _error("missing-record", field=f"{field}.source.inputs", value=identity)
    energy_values = tuple(
        cast(pint.Quantity[Any], explicit_values[item][1]) for item in energy_inputs
    )
    volume_values = tuple(cast(Decimal, explicit_values[item][1]) for item in volume_inputs)
    energy = (energy_values[0] + energy_values[1]) * Decimal("0.5")
    volume = (volume_values[0] * volume_values[1]).sqrt()
    prefix = f"association/{endpoint_a[0]}/{endpoint_a[1]}/{endpoint_b[0]}/{endpoint_b[1]}"
    energy_identity = f"{prefix}/energy_over_k"
    volume_identity = f"{prefix}/volume"
    derivatives = (
        *(
            AssociationDerivativeSource(energy_identity, identity, Decimal("0.5"))
            for identity in energy_inputs
        ),
        *(
            AssociationDerivativeSource(
                volume_identity,
                identity,
                volume / (Decimal(2) * value),
            )
            for identity, value in zip(volume_inputs, volume_values, strict=True)
        ),
    )
    provenance = _provenance_list(source_raw["provenance"], field=f"{field}.source.provenance")
    source = AssociationEdgeSource(
        "combining-rule",
        rule_id,
        provenance,
        (*energy_inputs, *volume_inputs),
        derivatives,
    )
    edge = AssociationEdge(
        endpoint_a,
        endpoint_b,
        energy_identity,
        volume_identity,
        energy,
        volume,
        source,
    )
    energy_spec = ParameterSpec(
        energy_identity,
        "association_energy_over_k",
        energy,
        "kelvin",
        False,
        *provenance[0],
    )
    volume_spec = ParameterSpec(
        volume_identity,
        "association_volume",
        volume,
        "dimensionless",
        False,
        *provenance[0],
    )
    return edge, (energy_spec, volume_spec)


def _topology(
    value: object,
    component_ids: frozenset[str],
) -> tuple[
    AssociationTopology,
    tuple[SiteRecord, ...],
    tuple[AssociationParameterRecord, ...],
    tuple[ParameterSpec, ...],
]:
    raw = _mapping(value, field="topology")
    _strict(
        raw,
        field="topology",
        required=frozenset({"presets", "sites", "edges"}),
    )

    presets: list[tuple[str, AssociationPreset]] = []
    sites: list[SiteRecord] = []
    preset_components: set[str] = set()
    for index, item in enumerate(_sequence(raw["presets"], field="topology.presets")):
        field = f"topology.presets[{index}]"
        preset_raw = _mapping(item, field=field)
        _strict(
            preset_raw,
            field=field,
            required=frozenset({"component_id", "preset", "provenance"}),
        )
        component_id = _identifier(preset_raw["component_id"], field=f"{field}.component_id")
        if component_id not in component_ids:
            raise _error("missing-record", field=f"{field}.component_id", value=component_id)
        preset = preset_raw["preset"]
        if not isinstance(preset, str) or preset not in _PRESET_SITES:
            raise _error("invalid-value", field=f"{field}.preset", value=preset)
        if component_id in preset_components:
            raise _error("duplicate-record", field="topology.preset.component_id")
        preset_components.add(component_id)
        provenance = _provenance(preset_raw["provenance"], field=f"{field}.provenance")
        presets.append((component_id, cast(AssociationPreset, preset)))
        for site_id, site_role in _PRESET_SITES[preset]:
            sites.append(
                SiteRecord(
                    f"site-{component_id}-{site_id}",
                    component_id,
                    site_id,
                    cast(Any, site_role),
                    1,
                    *provenance,
                )
            )

    for index, item in enumerate(_sequence(raw["sites"], field="topology.sites")):
        field = f"topology.sites[{index}]"
        site_raw = _mapping(item, field=field)
        _strict(
            site_raw,
            field=field,
            required=frozenset(
                {"component_id", "site_id", "site_role", "multiplicity", "provenance"}
            ),
        )
        component_id = _identifier(site_raw["component_id"], field=f"{field}.component_id")
        if component_id not in component_ids:
            raise _error("missing-record", field=f"{field}.component_id", value=component_id)
        if component_id in preset_components:
            raise _error(
                "invalid-schema",
                field=field,
                reason="preset-and-explicit-sites",
            )
        provenance = _provenance(site_raw["provenance"], field=f"{field}.provenance")
        sites.append(
            SiteRecord(
                f"site-{component_id}-{site_raw['site_id']}",
                component_id,
                cast(str, site_raw["site_id"]),
                cast(Any, site_raw["site_role"]),
                cast(int, site_raw["multiplicity"]),
                *provenance,
            )
        )

    ordered_sites = tuple(sorted(sites, key=lambda item: (item.component_id, item.site_id)))
    site_keys = frozenset((site.component_id, site.site_id) for site in ordered_sites)
    site_roles = {(site.component_id, site.site_id): site.site_role for site in ordered_sites}
    if len(site_keys) != len(ordered_sites):
        raise _error("duplicate-record", field="topology.site")

    pending: list[
        tuple[
            str,
            dict[str, object],
            tuple[str, str],
            tuple[str, str],
        ]
    ] = []
    edge_keys: set[tuple[tuple[str, str], tuple[str, str]]] = set()
    for index, item in enumerate(_sequence(raw["edges"], field="topology.edges")):
        field = f"topology.edges[{index}]"
        edge_raw = _mapping(item, field=field)
        endpoint_a = _endpoint(
            edge_raw.get("endpoint_a"), field=f"{field}.endpoint_a", sites=site_keys
        )
        endpoint_b = _endpoint(
            edge_raw.get("endpoint_b"), field=f"{field}.endpoint_b", sites=site_keys
        )
        endpoint_a, endpoint_b = sorted((endpoint_a, endpoint_b))
        edge_key = (endpoint_a, endpoint_b)
        if edge_key in edge_keys:
            raise _error("duplicate-record", field="topology.edge", value=edge_key)
        edge_keys.add(edge_key)
        source_raw = _mapping(edge_raw.get("source"), field=f"{field}.source")
        kind = source_raw.get("kind")
        if kind not in {"explicit", "combining-rule"}:
            raise _error("invalid-value", field=f"{field}.source.kind", value=kind)
        pending.append((field, edge_raw, endpoint_a, endpoint_b))

    edges: list[AssociationEdge] = []
    specs: list[ParameterSpec] = []
    explicit_values: dict[str, tuple[str, object]] = {}
    for field, edge_raw, endpoint_a, endpoint_b in pending:
        source_raw = _mapping(edge_raw["source"], field=f"{field}.source")
        if source_raw["kind"] != "explicit":
            continue
        edge, edge_specs = _explicit_edge(
            edge_raw,
            field=field,
            endpoint_a=endpoint_a,
            endpoint_b=endpoint_b,
        )
        edges.append(edge)
        specs.extend(edge_specs)
        explicit_values[edge.energy_identity] = (
            "association_energy_over_k",
            edge.energy_over_k,
        )
        explicit_values[edge.volume_identity] = ("association_volume", edge.volume)

    for field, edge_raw, endpoint_a, endpoint_b in pending:
        source_raw = _mapping(edge_raw["source"], field=f"{field}.source")
        if source_raw["kind"] != "combining-rule":
            continue
        edge, edge_specs = _combining_edge(
            edge_raw,
            field=field,
            endpoint_a=endpoint_a,
            endpoint_b=endpoint_b,
            explicit_values=explicit_values,
            site_roles=site_roles,
        )
        edges.append(edge)
        specs.extend(edge_specs)

    ordered_edges = tuple(sorted(edges, key=lambda edge: (edge.endpoint_a, edge.endpoint_b)))
    associations = tuple(record for edge in ordered_edges for record in _edge_records(edge))
    return (
        AssociationTopology(tuple(sorted(presets)), ordered_sites, ordered_edges),
        ordered_sites,
        associations,
        tuple(sorted(specs, key=lambda item: item.identity)),
    )


def _validate_document(
    *,
    sources: tuple[SourceRecord, ...],
    domains: tuple[ValidityDomain, ...],
    components: tuple[ComponentRecord, ...],
    records: tuple[
        SingleParameterRecord
        | PairParameterRecord
        | AssociationParameterRecord
        | ConstantCorrelation
        | ConstantPlusSumOfExponentialsCorrelation
        | ModelParameterRecord
        | SiteRecord,
        ...,
    ],
    model_families: tuple[ModelFamilySpec, ...],
    topology: AssociationTopology,
    parameter_specs: tuple[ParameterSpec, ...],
) -> None:
    def unique(values: Sequence[str], field: str) -> None:
        if len(values) != len(set(values)):
            raise _error("duplicate-record", field=field)

    if not sources or not domains or not components:
        raise _error("missing-record", field="sources-domains-components")
    unique(tuple(item.source_id for item in sources), "source_id")
    unique(tuple(item.domain_id for item in domains), "domain_id")
    unique(tuple(item.component_id for item in components), "component_id")
    unique(tuple(item.record_id for item in records), "record_id")
    unique(tuple(item.identity for item in parameter_specs), "parameter_identity")

    source_ids = {item.source_id for item in sources}
    domain_ids = {item.domain_id for item in domains}
    component_ids = {item.component_id for item in components}

    provenance_owners: list[object] = [*records, *model_families, *parameter_specs]
    for owner in provenance_owners:
        source_id = getattr(owner, "source_id", None)
        domain_id = getattr(owner, "domain_id", None)
        if source_id is not None and source_id not in source_ids:
            raise _error("invalid-source", source_id=source_id)
        if domain_id is not None and domain_id not in domain_ids:
            raise _error("invalid-domain", domain_id=domain_id)
    for edge in topology.edges:
        for source_id, _, domain_id in edge.source.provenance:
            if source_id not in source_ids:
                raise _error("invalid-source", source_id=source_id)
            if domain_id not in domain_ids:
                raise _error("invalid-domain", domain_id=domain_id)

    for record in records:
        for attribute in ("component_id", "component_id_a", "component_id_b"):
            component_id = getattr(record, attribute, None)
            if component_id is not None and component_id not in component_ids:
                raise _error("missing-record", field=attribute, value=component_id)

    fixed_identities = {
        record.identity
        for record in records
        if isinstance(record, SingleParameterRecord) and record.family in _FIXED_FAMILIES
    }
    if fixed_identities & {item.identity for item in parameter_specs}:
        raise _error("invalid-schema", field="parameter_specs", reason="fixed-coordinate")
    derived_identities = {
        identity
        for edge in topology.edges
        if edge.source.kind == "combining-rule"
        for identity in (edge.energy_identity, edge.volume_identity)
    }
    noncontinuous_identities = {item.identity for item in parameter_specs if not item.continuous}
    if noncontinuous_identities != derived_identities:
        raise _error(
            "invalid-schema",
            field="parameter_specs",
            reason="derived-coordinate-classification",
        )

    direct_component_families = {
        (record.component_id, record.family)
        for record in records
        if isinstance(record, SingleParameterRecord)
    }
    correlation_families = {
        (record.component_id, record.family)
        for record in records
        if isinstance(record, (ConstantCorrelation, ConstantPlusSumOfExponentialsCorrelation))
    }
    if direct_component_families & correlation_families:
        raise _error("duplicate-record", field="component-correlation-family")

    family_choices = {item.kind: item.choice for item in model_families}
    has_association = bool(topology.sites or topology.edges)
    if family_choices["association"] == "none" and has_association:
        raise _error("invalid-schema", field="topology", reason="association-disabled")
    if family_choices["association"] == "general-site" and (
        not topology.sites or not topology.edges
    ):
        raise _error("missing-record", field="topology", reason="association-enabled")


def _records_from_mapping(value: object) -> _DocumentRecords:
    raw = _mapping(value, field="parameters")
    _strict(raw, field="parameters", required=_TOP_LEVEL)
    if raw["schema"] != "epcsaft.parameters" or raw["schema_version"] != 1:
        raise _error("invalid-schema", field="schema", value=raw.get("schema_version"))
    document_id = _identifier(raw["document_id"], field="document_id")
    document_version = _positive_version(raw["document_version"], field="document_version")
    purpose = raw["purpose"]
    if purpose not in {
        "external-document",
        "package-test-fixture",
        "regression-output",
        "user-provided",
    }:
        raise _error("invalid-value", field="purpose", value=purpose)

    sources = _sources(raw["sources"])
    domains = _domains(raw["domains"])
    components, singles, component_specs = _components(raw["components"])
    component_ids = frozenset(item.component_id for item in components)
    pairs, pair_specs = _pairs(raw["pairs"], component_ids)
    model_families = _model_families(raw["model_families"])
    models, model_specs = _model_coefficients(raw["model_coefficients"])
    correlations, correlation_specs = _correlations(raw["correlations"], component_ids)
    topology, sites, associations, association_specs = _topology(raw["topology"], component_ids)
    parameter_specs = tuple(
        sorted(
            (
                *component_specs,
                *pair_specs,
                *model_specs,
                *correlation_specs,
                *association_specs,
            ),
            key=lambda item: item.identity,
        )
    )
    records = (
        *singles,
        *pairs,
        *sites,
        *associations,
        *correlations,
        *models,
    )
    _validate_document(
        sources=sources,
        domains=domains,
        components=components,
        records=records,
        model_families=model_families,
        topology=topology,
        parameter_specs=parameter_specs,
    )
    return _DocumentRecords(
        component_ids=tuple(item.component_id for item in components),
        document_id=document_id,
        document_version=document_version,
        purpose=purpose,
        sources=sources,
        domains=domains,
        components=components,
        singles=singles,
        pairs=pairs,
        sites=sites,
        associations=associations,
        correlations=correlations,
        models=models,
        model_families=model_families,
        topology=topology,
        parameter_specs=parameter_specs,
    )


__all__ = ("_records_from_mapping",)
