from __future__ import annotations

from dataclasses import dataclass
from typing import ClassVar, Literal, TypeAlias, get_args
import re

from .units import (
    _decimal_number,
    _declared_quantity,
    _error,
    _exact_integer,
    _quantity,
)


_ID_PATTERN = re.compile(r"[a-z0-9]+(?:-[a-z0-9]+)*\Z")
DomainKind: TypeAlias = Literal["fit-range", "reported-conditions", "unknown"]
SingleParameterFamily: TypeAlias = Literal[
    "molar_mass",
    "segment_count",
    "segment_diameter",
    "dispersion_energy_over_k",
    "charge_number",
    "relative_permittivity",
    "born_diameter",
    "solvation_factor",
    "schreckenberg_dielectric_volume",
    "schreckenberg_dielectric_temperature",
    "zuber_ion_suppression_coefficient",
    "rueben_dipole_scaling",
    "rueben_polarizability_scaling",
    "rueben_correlation_integral_parameter",
]
PairParameterFamily: TypeAlias = Literal["k_ij", "l_ij"]
AssociationParameterFamily: TypeAlias = Literal[
    "association_energy_over_k", "association_volume"
]
ModelParameterFamily: TypeAlias = Literal[
    "dielectric_ion_suppression_coefficient",
    "ionic_region_relative_permittivity",
]
CorrelationFamily: TypeAlias = Literal[
    "segment_diameter", "relative_permittivity", "solvation_factor"
]

_DOMAIN_KINDS = frozenset(get_args(DomainKind))
_SINGLE_FAMILIES = frozenset(get_args(SingleParameterFamily))
_PAIR_FAMILIES = frozenset(get_args(PairParameterFamily))
_ASSOCIATION_FAMILIES = frozenset(get_args(AssociationParameterFamily))
_MODEL_FAMILIES = frozenset(get_args(ModelParameterFamily))
_CORRELATION_FAMILIES = frozenset(get_args(CorrelationFamily))
_DIMENSIONAL_UNITS = {
    "molar_mass": "kilogram / mole",
    "segment_diameter": "angstrom",
    "dispersion_energy_over_k": "kelvin",
    "born_diameter": "angstrom",
    "schreckenberg_dielectric_volume": "meter ** 3 / mole",
    "schreckenberg_dielectric_temperature": "kelvin",
    "rueben_dipole_scaling": "debye ** 2",
    "rueben_polarizability_scaling": "angstrom ** 3",
    "association_energy_over_k": "kelvin",
}
_DIMENSIONLESS_FAMILIES = frozenset(
    {
        "segment_count",
        "relative_permittivity",
        "solvation_factor",
        "association_volume",
        "k_ij",
        "l_ij",
        "dielectric_ion_suppression_coefficient",
        "ionic_region_relative_permittivity",
        "zuber_ion_suppression_coefficient",
        "rueben_correlation_integral_parameter",
    }
)


def _identifier(value: object, *, field: str) -> str:
    if not isinstance(value, str) or _ID_PATTERN.fullmatch(value) is None:
        raise _error("invalid-id", field=field, value=value)
    return value


def _text(value: object, *, field: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise _error("invalid-value", field=field, value=value)
    return value


def _optional_text(value: object, *, field: str) -> str | None:
    if value is None:
        return None
    return _text(value, field=field)


def _family(value: object, *, allowed: frozenset[str], field: str = "family") -> str:
    if not isinstance(value, str) or value not in allowed:
        raise _error("invalid-value", field=field, value=value)
    return value


def _parameter_value(family: str, value: object) -> object:
    if family == "charge_number":
        return _exact_integer(value, field="value")
    if family in _DIMENSIONAL_UNITS:
        return _quantity(value, _DIMENSIONAL_UNITS[family], field="value")
    if family in _DIMENSIONLESS_FAMILIES:
        return _decimal_number(value, field="value")
    raise _error("invalid-value", field="family", value=family)


def _provenance(
    source_id: object,
    locator: object,
    domain_id: object,
) -> tuple[str, str, str]:
    return (
        _identifier(source_id, field="source_id"),
        _text(locator, field="locator"),
        _identifier(domain_id, field="domain_id"),
    )


@dataclass(frozen=True, slots=True)
class SourceRecord:
    source_id: str
    citation: str
    use_basis: str
    doi: str | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "source_id", _identifier(self.source_id, field="source_id"))
        object.__setattr__(self, "citation", _text(self.citation, field="citation"))
        object.__setattr__(self, "use_basis", _text(self.use_basis, field="use_basis"))
        object.__setattr__(self, "doi", _optional_text(self.doi, field="doi"))


@dataclass(frozen=True, slots=True)
class ValidityDomain:
    domain_id: str
    kind: DomainKind
    temperature_min: object | None = None
    temperature_max: object | None = None
    pressure_min: object | None = None
    pressure_max: object | None = None
    ion_mole_fraction_max: object | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "domain_id", _identifier(self.domain_id, field="domain_id"))
        if self.kind not in _DOMAIN_KINDS:
            raise _error("invalid-domain", field="kind", value=self.kind)

        normalized = {
            "temperature_min": None
            if self.temperature_min is None
            else _quantity(self.temperature_min, "kelvin", field="temperature_min"),
            "temperature_max": None
            if self.temperature_max is None
            else _quantity(self.temperature_max, "kelvin", field="temperature_max"),
            "pressure_min": None
            if self.pressure_min is None
            else _quantity(self.pressure_min, "pascal", field="pressure_min"),
            "pressure_max": None
            if self.pressure_max is None
            else _quantity(self.pressure_max, "pascal", field="pressure_max"),
            "ion_mole_fraction_max": None
            if self.ion_mole_fraction_max is None
            else _decimal_number(self.ion_mole_fraction_max, field="ion_mole_fraction_max"),
        }
        for field, value in normalized.items():
            object.__setattr__(self, field, value)

        bounds = (
            (
                "temperature",
                normalized["temperature_min"],
                normalized["temperature_max"],
            ),
            ("pressure", normalized["pressure_min"], normalized["pressure_max"]),
        )
        ion_max = normalized["ion_mole_fraction_max"]
        if ion_max is not None and not (0 < ion_max <= 1):
            raise _error("invalid-domain", field="ion_mole_fraction_max", value=ion_max)
        if self.kind == "unknown" and any(value is not None for value in normalized.values()):
            raise _error("invalid-domain", field="kind", value=self.kind)
        if self.kind != "unknown" and all(value is None for value in normalized.values()):
            raise _error("invalid-domain", field="bounds", value="missing")
        for field, lower, upper in bounds:
            if lower is not None and upper is not None and lower.magnitude > upper.magnitude:
                raise _error("invalid-domain", field=field, value="reversed-bounds")


@dataclass(frozen=True, slots=True)
class ComponentRecord:
    component_id: str
    name: str | None = None
    cas: str | None = None
    formula: str | None = None
    aliases: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        object.__setattr__(
            self, "component_id", _identifier(self.component_id, field="component_id")
        )
        object.__setattr__(self, "name", _optional_text(self.name, field="name"))
        object.__setattr__(self, "cas", _optional_text(self.cas, field="cas"))
        object.__setattr__(self, "formula", _optional_text(self.formula, field="formula"))
        if isinstance(self.aliases, str):
            raise _error("invalid-value", field="aliases", value=self.aliases)
        aliases = tuple(_text(alias, field="alias") for alias in self.aliases)
        if len(aliases) != len(set(aliases)):
            raise _error("invalid-value", field="aliases", value="duplicate")
        object.__setattr__(self, "aliases", aliases)


@dataclass(frozen=True, slots=True)
class SingleParameterRecord:
    record_id: str
    component_id: str
    family: SingleParameterFamily
    value: object
    source_id: str
    locator: str
    domain_id: str

    def __post_init__(self) -> None:
        object.__setattr__(self, "record_id", _identifier(self.record_id, field="record_id"))
        object.__setattr__(
            self, "component_id", _identifier(self.component_id, field="component_id")
        )
        family = _family(self.family, allowed=_SINGLE_FAMILIES)
        object.__setattr__(self, "family", family)
        object.__setattr__(self, "value", _parameter_value(family, self.value))
        source_id, locator, domain_id = _provenance(self.source_id, self.locator, self.domain_id)
        object.__setattr__(self, "source_id", source_id)
        object.__setattr__(self, "locator", locator)
        object.__setattr__(self, "domain_id", domain_id)


@dataclass(frozen=True, slots=True)
class ModelParameterRecord:
    record_id: str
    family: ModelParameterFamily
    value: object
    source_id: str
    locator: str
    domain_id: str

    def __post_init__(self) -> None:
        object.__setattr__(self, "record_id", _identifier(self.record_id, field="record_id"))
        family = _family(self.family, allowed=_MODEL_FAMILIES)
        object.__setattr__(self, "family", family)
        object.__setattr__(self, "value", _parameter_value(family, self.value))
        source_id, locator, domain_id = _provenance(self.source_id, self.locator, self.domain_id)
        object.__setattr__(self, "source_id", source_id)
        object.__setattr__(self, "locator", locator)
        object.__setattr__(self, "domain_id", domain_id)


@dataclass(frozen=True, slots=True)
class PairParameterRecord:
    record_id: str
    component_id_a: str
    component_id_b: str
    family: PairParameterFamily
    value: object
    source_id: str
    locator: str
    domain_id: str

    def __post_init__(self) -> None:
        object.__setattr__(self, "record_id", _identifier(self.record_id, field="record_id"))
        component_a = _identifier(self.component_id_a, field="component_id_a")
        component_b = _identifier(self.component_id_b, field="component_id_b")
        if component_a == component_b:
            raise _error("invalid-value", field="component_pair", value=component_a)
        component_a, component_b = sorted((component_a, component_b))
        object.__setattr__(self, "component_id_a", component_a)
        object.__setattr__(self, "component_id_b", component_b)
        family = _family(self.family, allowed=_PAIR_FAMILIES)
        object.__setattr__(self, "family", family)
        object.__setattr__(self, "value", _parameter_value(family, self.value))
        source_id, locator, domain_id = _provenance(self.source_id, self.locator, self.domain_id)
        object.__setattr__(self, "source_id", source_id)
        object.__setattr__(self, "locator", locator)
        object.__setattr__(self, "domain_id", domain_id)


@dataclass(frozen=True, slots=True)
class SiteRecord:
    record_id: str
    component_id: str
    site_id: str
    site_class: str
    multiplicity: int
    source_id: str
    locator: str
    domain_id: str

    def __post_init__(self) -> None:
        object.__setattr__(self, "record_id", _identifier(self.record_id, field="record_id"))
        object.__setattr__(
            self, "component_id", _identifier(self.component_id, field="component_id")
        )
        object.__setattr__(self, "site_id", _identifier(self.site_id, field="site_id"))
        object.__setattr__(self, "site_class", _identifier(self.site_class, field="site_class"))
        object.__setattr__(
            self,
            "multiplicity",
            _exact_integer(self.multiplicity, field="multiplicity", positive=True),
        )
        source_id, locator, domain_id = _provenance(self.source_id, self.locator, self.domain_id)
        object.__setattr__(self, "source_id", source_id)
        object.__setattr__(self, "locator", locator)
        object.__setattr__(self, "domain_id", domain_id)


@dataclass(frozen=True, slots=True)
class AssociationParameterRecord:
    record_id: str
    component_id_a: str
    site_id_a: str
    component_id_b: str
    site_id_b: str
    family: AssociationParameterFamily
    value: object
    source_id: str
    locator: str
    domain_id: str

    def __post_init__(self) -> None:
        object.__setattr__(self, "record_id", _identifier(self.record_id, field="record_id"))
        endpoint_a = (
            _identifier(self.component_id_a, field="component_id_a"),
            _identifier(self.site_id_a, field="site_id_a"),
        )
        endpoint_b = (
            _identifier(self.component_id_b, field="component_id_b"),
            _identifier(self.site_id_b, field="site_id_b"),
        )
        endpoint_a, endpoint_b = sorted((endpoint_a, endpoint_b))
        object.__setattr__(self, "component_id_a", endpoint_a[0])
        object.__setattr__(self, "site_id_a", endpoint_a[1])
        object.__setattr__(self, "component_id_b", endpoint_b[0])
        object.__setattr__(self, "site_id_b", endpoint_b[1])
        family = _family(self.family, allowed=_ASSOCIATION_FAMILIES)
        object.__setattr__(self, "family", family)
        object.__setattr__(self, "value", _parameter_value(family, self.value))
        source_id, locator, domain_id = _provenance(self.source_id, self.locator, self.domain_id)
        object.__setattr__(self, "source_id", source_id)
        object.__setattr__(self, "locator", locator)
        object.__setattr__(self, "domain_id", domain_id)


@dataclass(frozen=True, slots=True)
class ExponentialTerm:
    amplitude: object
    exponent_coefficient: object

    def __post_init__(self) -> None:
        amplitude = _declared_quantity(self.amplitude, field="amplitude")
        if amplitude.dimensionless:
            raise _error("invalid-unit", field="amplitude", expected="dimensional")
        object.__setattr__(self, "amplitude", amplitude)
        object.__setattr__(
            self,
            "exponent_coefficient",
            _quantity(
                self.exponent_coefficient,
                "1 / kelvin",
                field="exponent_coefficient",
            ),
        )


@dataclass(frozen=True, slots=True)
class ConstantCorrelation:
    form: ClassVar[str] = "constant"

    record_id: str
    component_id: str
    family: CorrelationFamily
    constant: object
    source_id: str
    locator: str
    domain_id: str

    def __post_init__(self) -> None:
        object.__setattr__(self, "record_id", _identifier(self.record_id, field="record_id"))
        object.__setattr__(
            self, "component_id", _identifier(self.component_id, field="component_id")
        )
        family = _family(self.family, allowed=_CORRELATION_FAMILIES)
        object.__setattr__(self, "family", family)
        object.__setattr__(self, "constant", _parameter_value(family, self.constant))
        source_id, locator, domain_id = _provenance(self.source_id, self.locator, self.domain_id)
        object.__setattr__(self, "source_id", source_id)
        object.__setattr__(self, "locator", locator)
        object.__setattr__(self, "domain_id", domain_id)


@dataclass(frozen=True, slots=True)
class ConstantPlusSumOfExponentialsCorrelation:
    form: ClassVar[str] = "constant-plus-sum-of-exponentials"

    record_id: str
    component_id: str
    family: Literal["segment_diameter"]
    constant: object
    terms: tuple[ExponentialTerm, ...]
    source_id: str
    locator: str
    domain_id: str

    def __post_init__(self) -> None:
        object.__setattr__(self, "record_id", _identifier(self.record_id, field="record_id"))
        object.__setattr__(
            self, "component_id", _identifier(self.component_id, field="component_id")
        )
        family = _family(self.family, allowed=frozenset({"segment_diameter"}))
        object.__setattr__(self, "family", family)
        object.__setattr__(self, "constant", _parameter_value(family, self.constant))
        terms = tuple(self.terms)
        if not terms or not all(isinstance(term, ExponentialTerm) for term in terms):
            raise _error("invalid-correlation", field="terms", value="invalid")
        normalized_terms = tuple(
            ExponentialTerm(
                amplitude=_quantity(term.amplitude, "angstrom", field="amplitude"),
                exponent_coefficient=term.exponent_coefficient,
            )
            for term in terms
        )
        object.__setattr__(self, "terms", normalized_terms)
        source_id, locator, domain_id = _provenance(self.source_id, self.locator, self.domain_id)
        object.__setattr__(self, "source_id", source_id)
        object.__setattr__(self, "locator", locator)
        object.__setattr__(self, "domain_id", domain_id)


CorrelationRecord: TypeAlias = ConstantCorrelation | ConstantPlusSumOfExponentialsCorrelation


__all__ = (
    "AssociationParameterFamily",
    "AssociationParameterRecord",
    "ComponentRecord",
    "ConstantCorrelation",
    "ConstantPlusSumOfExponentialsCorrelation",
    "CorrelationFamily",
    "CorrelationRecord",
    "DomainKind",
    "ExponentialTerm",
    "ModelParameterFamily",
    "ModelParameterRecord",
    "PairParameterFamily",
    "PairParameterRecord",
    "SingleParameterFamily",
    "SingleParameterRecord",
    "SiteRecord",
    "SourceRecord",
    "ValidityDomain",
)
