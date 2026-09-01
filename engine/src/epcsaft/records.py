from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Any, ClassVar, Literal, Protocol, cast, get_args

import pint

from .units import (
    _decimal_number,
    _declared_quantity,
    _error,
    _exact_integer,
    _quantity,
)

_ID_PATTERN = re.compile(r"[a-z0-9]+(?:-[a-z0-9]+)*\Z")
_FAMILY_PATTERN = re.compile(r"[a-z0-9]+(?:[-_][a-z0-9]+)*\Z")
_REACTION_PARAMETER_PATTERN = re.compile(r"reaction:R[0-9]+:correlation:[a-z][a-z0-9_]*\Z")
type DomainKind = Literal["fit-range", "reported-conditions", "unknown"]
type ModelQualification = Literal["candidate_extrapolation"]
type ModelFamilyKind = Literal["base", "association", "electrolyte", "permittivity"]
type AssociationSiteRole = Literal["donor", "acceptor", "effective"]
type AssociationEdgeSourceKind = Literal["explicit", "combining-rule"]
type AssociationPreset = Literal["1a", "2b", "3b", "4c"]
type AssociationCombiningRuleId = Literal["arithmetic-energy-geometric-volume",]
type RelativePermittivityFormulation = Literal[
    "none",
    "ion-fraction-suppression",
    "component-permittivity-mixing",
    "density-volume-temperature",
    "ion-specific-suppression",
    "dipolar-density-perturbation",
    "solvent-only",
]
type IonicDispersionFormulation = Literal[
    "all-ion-pairs",
    "exclude-same-sign-ion-pairs",
]
type SingleParameterFamily = Literal[
    "molar_mass",
    "segment_count",
    "segment_diameter",
    "dispersion_energy_over_k",
    "charge_number",
    "relative_permittivity",
    "packing_diameter",
    "debye_huckel_diameter",
    "born_diameter",
    "solvation_factor",
    "dielectric_volume",
    "dielectric_temperature",
    "ion_specific_suppression_coefficient",
    "dipole_scaling",
    "polarizability_scaling",
    "correlation_integral_parameter",
]
type PairParameterFamily = Literal[
    "k_ij", "l_ij", "k_ij_reciprocal_temperature_slope"
]
type CrossAssociationCorrectionFamily = Literal["k_hb_ij"]
type AssociationParameterFamily = Literal["association_energy_over_k", "association_volume"]
type ModelParameterFamily = Literal[
    "ion_fraction_suppression_coefficient",
    "ionic_region_relative_permittivity",
]
type CorrelationFamily = Literal["segment_diameter", "relative_permittivity", "solvation_factor"]


class _CandidateOwner(Protocol):
    @property
    def candidate_domain_id(self) -> str | None: ...

    @property
    def qualification(self) -> ModelQualification | None: ...

    @property
    def source_sha256(self) -> str | None: ...


_DOMAIN_KINDS = frozenset(get_args(DomainKind.__value__))
_SINGLE_FAMILIES = frozenset(get_args(SingleParameterFamily.__value__))
_PAIR_FAMILIES = frozenset(get_args(PairParameterFamily.__value__))
_ASSOCIATION_FAMILIES = frozenset(get_args(AssociationParameterFamily.__value__))
_MODEL_FAMILIES = frozenset(get_args(ModelParameterFamily.__value__))
_CORRELATION_FAMILIES = frozenset(get_args(CorrelationFamily.__value__))
_ASSOCIATION_PRESETS = frozenset(get_args(AssociationPreset.__value__))
_ASSOCIATION_COMBINING_RULES = frozenset(
    get_args(AssociationCombiningRuleId.__value__)
)
_DIMENSIONAL_UNITS = {
    "molar_mass": "kilogram / mole",
    "segment_diameter": "angstrom",
    "packing_diameter": "angstrom",
    "debye_huckel_diameter": "angstrom",
    "dispersion_energy_over_k": "kelvin",
    "born_diameter": "angstrom",
    "dielectric_volume": "meter ** 3 / mole",
    "dielectric_temperature": "kelvin",
    "dipole_scaling": "debye ** 2",
    "polarizability_scaling": "angstrom ** 3",
    "association_energy_over_k": "kelvin",
    "k_ij_reciprocal_temperature_slope": "kelvin",
}
_DIMENSIONLESS_FAMILIES = frozenset(
    {
        "segment_count",
        "relative_permittivity",
        "solvation_factor",
        "association_volume",
        "k_ij",
        "l_ij",
        "k_hb_ij",
        "ion_fraction_suppression_coefficient",
        "ionic_region_relative_permittivity",
        "ion_specific_suppression_coefficient",
        "correlation_integral_parameter",
    }
)


def _identifier(value: object, *, field: str) -> str:
    if not isinstance(value, str) or _ID_PATTERN.fullmatch(value) is None:
        raise _error("invalid-id", field=field, value=value)
    return value


def _parameter_identity(value: object, *, field: str = "identity") -> str:
    if (
        not isinstance(value, str)
        or not value
        or (
            _REACTION_PARAMETER_PATTERN.fullmatch(value) is None
            and any(_FAMILY_PATTERN.fullmatch(part) is None for part in value.split("/"))
        )
    ):
        raise _error("invalid-id", field=field, value=value)
    return value


def _family_identity(value: object, *, field: str = "family") -> str:
    if not isinstance(value, str) or _FAMILY_PATTERN.fullmatch(value) is None:
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


def _candidate_provenance(
    candidate_domain_id: object,
    qualification: object,
    source_sha256: object,
    *,
    field: str,
) -> tuple[str | None, ModelQualification | None, str | None]:
    values = (candidate_domain_id, qualification, source_sha256)
    if all(value is None for value in values):
        return None, None, None
    if any(value is None for value in values):
        raise _error("invalid-provenance", field=field, value="partial")
    domain_id = _identifier(candidate_domain_id, field=f"{field}.candidate_domain_id")
    if qualification != "candidate_extrapolation":
        raise _error("invalid-value", field=f"{field}.qualification", value=qualification)
    if (
        not isinstance(source_sha256, str)
        or re.fullmatch(r"sha256:[0-9a-f]{64}", source_sha256) is None
    ):
        raise _error("invalid-value", field=f"{field}.source_sha256", value=source_sha256)
    return domain_id, "candidate_extrapolation", source_sha256


def _normalize_candidate(owner: _CandidateOwner, *, field: str) -> None:
    candidate = _candidate_provenance(
        owner.candidate_domain_id,
        owner.qualification,
        owner.source_sha256,
        field=field,
    )
    for name, value in zip(
        ("candidate_domain_id", "qualification", "source_sha256"),
        candidate,
        strict=True,
    ):
        object.__setattr__(owner, name, value)


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
class ModelFamilySpec:
    """One closed, fixed model-family choice in a Parameters document."""

    family_id: str
    kind: ModelFamilyKind
    choice: str
    source_id: str
    locator: str
    domain_id: str
    c_shell: float | None = None
    c_dielectric: float | None = None
    ionic_dispersion: IonicDispersionFormulation | None = None
    candidate_domain_id: str | None = None
    qualification: ModelQualification | None = None
    source_sha256: str | None = None

    def __post_init__(self) -> None:
        object.__setattr__(
            self, "family_id", _parameter_identity(self.family_id, field="family_id")
        )
        if self.kind not in {"base", "association", "electrolyte", "permittivity"}:
            raise _error("invalid-value", field="model_family.kind", value=self.kind)
        object.__setattr__(self, "choice", _identifier(self.choice, field="model_family.choice"))
        corrections = (self.c_shell, self.c_dielectric)
        if self.kind == "electrolyte" and self.choice == "born":
            if corrections not in {(0.0, 0.0), (1.0, 0.0), (1.0, 1.0)}:
                raise _error(
                    "invalid-value",
                    field="model_family.born_corrections",
                    value=corrections,
                )
        elif corrections != (None, None):
            raise _error(
                "invalid-value",
                field="model_family.born_corrections",
                reason="non-born-model",
            )
        if self.kind == "electrolyte" and self.choice != "none":
            if self.ionic_dispersion is None:
                object.__setattr__(
                    self,
                    "ionic_dispersion",
                    "exclude-same-sign-ion-pairs",
                )
            if self.ionic_dispersion not in {
                "all-ion-pairs",
                "exclude-same-sign-ion-pairs",
            }:
                raise _error(
                    "invalid-value",
                    field="model_family.ionic_dispersion",
                    value=self.ionic_dispersion,
                )
        elif self.ionic_dispersion is not None:
            raise _error(
                "invalid-value",
                field="model_family.ionic_dispersion",
                reason="non-electrolyte-model",
            )
        source_id, locator, domain_id = _provenance(self.source_id, self.locator, self.domain_id)
        object.__setattr__(self, "source_id", source_id)
        object.__setattr__(self, "locator", locator)
        object.__setattr__(self, "domain_id", domain_id)
        _normalize_candidate(self, field="model_family.candidate")


@dataclass(frozen=True, slots=True)
class ParameterSpec:
    """Stable identity and selection metadata for one stored parameter value."""

    identity: str
    family: str
    value: object
    unit: str
    continuous: bool
    source_id: str | None = None
    locator: str | None = None
    domain_id: str | None = None
    candidate_domain_id: str | None = None
    qualification: ModelQualification | None = None
    source_sha256: str | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "identity", _parameter_identity(self.identity))
        object.__setattr__(self, "family", _family_identity(self.family))
        object.__setattr__(self, "unit", _text(self.unit, field="unit"))
        if self.continuous is not True and self.continuous is not False:
            raise _error("invalid-value", field="continuous", value=self.continuous)
        provenance = (self.source_id, self.locator, self.domain_id)
        if any(value is None for value in provenance):
            if not all(value is None for value in provenance):
                raise _error("invalid-provenance", field="parameter_spec", value="partial")
        else:
            source_id, locator, domain_id = _provenance(
                cast(str, self.source_id),
                cast(str, self.locator),
                cast(str, self.domain_id),
            )
            object.__setattr__(self, "source_id", source_id)
            object.__setattr__(self, "locator", locator)
            object.__setattr__(self, "domain_id", domain_id)
        _normalize_candidate(self, field="parameter_spec.candidate")


@dataclass(frozen=True, slots=True)
class ReactionCorrelationRecord:
    """One optimizer-ready coefficient in a typed reaction correlation."""

    record_id: str
    reaction_id: str
    correlation_kind: Literal["ln-k-a-plus-b-over-t", "negative-log10-temperature-polynomial"]
    coefficient_name: str
    identity: str
    value: object
    unit: str
    optimizer_lower: object
    optimizer_upper: object
    source_id: str
    locator: str
    domain_id: str
    candidate_domain_id: str
    qualification: ModelQualification
    source_sha256: str
    family: str = "reaction_correlation"

    def __post_init__(self) -> None:
        object.__setattr__(
            self, "record_id", _parameter_identity(self.record_id, field="record_id")
        )
        if (
            not isinstance(self.reaction_id, str)
            or re.fullmatch(r"R[0-9]+", self.reaction_id) is None
        ):
            raise _error("invalid-id", field="reaction_id", value=self.reaction_id)
        if self.correlation_kind not in {
            "ln-k-a-plus-b-over-t",
            "negative-log10-temperature-polynomial",
        }:
            raise _error(
                "invalid-value", field="reaction.correlation_kind", value=self.correlation_kind
            )
        object.__setattr__(
            self,
            "coefficient_name",
            _family_identity(self.coefficient_name, field="coefficient_name"),
        )
        identity = _parameter_identity(self.identity)
        if identity != f"reaction:{self.reaction_id}:correlation:{self.coefficient_name}":
            raise _error("invalid-id", field="identity", value=identity)
        object.__setattr__(self, "identity", identity)
        unit = _text(self.unit, field="unit")
        value = _quantity(self.value, unit, field="value")
        lower = _quantity(self.optimizer_lower, unit, field="optimizer_lower")
        upper = _quantity(self.optimizer_upper, unit, field="optimizer_upper")
        if lower.magnitude >= upper.magnitude or not (
            lower.magnitude <= value.magnitude <= upper.magnitude
        ):
            raise _error("invalid-domain", field="reaction.optimizer_domain", value=identity)
        object.__setattr__(self, "unit", unit)
        object.__setattr__(self, "value", value)
        object.__setattr__(self, "optimizer_lower", lower)
        object.__setattr__(self, "optimizer_upper", upper)
        source_id, locator, domain_id = _provenance(self.source_id, self.locator, self.domain_id)
        object.__setattr__(self, "source_id", source_id)
        object.__setattr__(self, "locator", locator)
        object.__setattr__(self, "domain_id", domain_id)
        _normalize_candidate(self, field="reaction.candidate")


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
            if (
                lower is not None
                and upper is not None
                and cast(pint.Quantity[Any], lower).magnitude
                > cast(pint.Quantity[Any], upper).magnitude
            ):
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
    candidate_domain_id: str | None = None
    qualification: ModelQualification | None = None
    source_sha256: str | None = None

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
        _normalize_candidate(self, field="component_parameter.candidate")

    @property
    def identity(self) -> str:
        return f"component/{self.component_id}/{self.family}"


@dataclass(frozen=True, slots=True)
class ModelParameterRecord:
    record_id: str
    family: ModelParameterFamily
    value: object
    source_id: str
    locator: str
    domain_id: str
    candidate_domain_id: str | None = None
    qualification: ModelQualification | None = None
    source_sha256: str | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "record_id", _identifier(self.record_id, field="record_id"))
        family = _family(self.family, allowed=_MODEL_FAMILIES)
        object.__setattr__(self, "family", family)
        object.__setattr__(self, "value", _parameter_value(family, self.value))
        source_id, locator, domain_id = _provenance(self.source_id, self.locator, self.domain_id)
        object.__setattr__(self, "source_id", source_id)
        object.__setattr__(self, "locator", locator)
        object.__setattr__(self, "domain_id", domain_id)
        _normalize_candidate(self, field="model_parameter.candidate")

    @property
    def identity(self) -> str:
        return f"model/{self.family}"


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
    candidate_domain_id: str | None = None
    qualification: ModelQualification | None = None
    source_sha256: str | None = None
    reference_temperature: object | None = None

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
        _normalize_candidate(self, field="pair_parameter.candidate")
        if family == "k_ij_reciprocal_temperature_slope":
            if self.reference_temperature is None:
                raise _error(
                    "invalid-correlation", field="reference_temperature", value=None
                )
            object.__setattr__(
                self,
                "reference_temperature",
                _quantity(
                    self.reference_temperature,
                    "kelvin",
                    field="reference_temperature",
                ),
            )
            if cast(pint.Quantity[Any], self.reference_temperature).magnitude <= 0.0:
                raise _error(
                    "invalid-correlation",
                    field="reference_temperature",
                    value=self.reference_temperature,
                )
        elif self.reference_temperature is not None:
            raise _error(
                "invalid-correlation",
                field="reference_temperature",
                value=self.reference_temperature,
            )

    @property
    def identity(self) -> str:
        if self.family == "k_ij_reciprocal_temperature_slope":
            return (
                f"pair/{self.component_id_a}/{self.component_id_b}/"
                "k_ij/reciprocal_temperature_slope"
            )
        return f"pair/{self.component_id_a}/{self.component_id_b}/{self.family}"


@dataclass(frozen=True, slots=True)
class CrossAssociationCorrectionRecord:
    """One unordered component- and role-pair correction to combined energy."""

    record_id: str
    component_id_a: str
    component_id_b: str
    role_a: AssociationSiteRole
    role_b: AssociationSiteRole
    family: CrossAssociationCorrectionFamily
    value: object
    source_id: str
    locator: str
    domain_id: str
    candidate_domain_id: str | None = None
    qualification: ModelQualification | None = None
    source_sha256: str | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "record_id", _identifier(self.record_id, field="record_id"))
        component_a = _identifier(self.component_id_a, field="component_id_a")
        component_b = _identifier(self.component_id_b, field="component_id_b")
        if component_a == component_b:
            raise _error("invalid-value", field="component_pair", value=component_a)
        component_a, component_b = sorted((component_a, component_b))
        object.__setattr__(self, "component_id_a", component_a)
        object.__setattr__(self, "component_id_b", component_b)
        roles = (self.role_a, self.role_b)
        if any(role not in {"donor", "acceptor", "effective"} for role in roles):
            raise _error("invalid-value", field="association_role_pair", value=roles)
        role_a, role_b = sorted(roles)
        object.__setattr__(self, "role_a", role_a)
        object.__setattr__(self, "role_b", role_b)
        if self.family != "k_hb_ij":
            raise _error("invalid-value", field="family", value=self.family)
        object.__setattr__(self, "value", _parameter_value(self.family, self.value))
        source_id, locator, domain_id = _provenance(self.source_id, self.locator, self.domain_id)
        object.__setattr__(self, "source_id", source_id)
        object.__setattr__(self, "locator", locator)
        object.__setattr__(self, "domain_id", domain_id)
        _normalize_candidate(self, field="cross_association.candidate")

    @property
    def identity(self) -> str:
        return (
            f"pair/{self.component_id_a}/{self.component_id_b}/association/"
            f"{self.role_a}/{self.role_b}/{self.family}"
        )


@dataclass(frozen=True, slots=True)
class SiteRecord:
    record_id: str
    component_id: str
    site_id: str
    site_role: AssociationSiteRole
    multiplicity: int
    source_id: str
    locator: str
    domain_id: str
    candidate_domain_id: str | None = None
    qualification: ModelQualification | None = None
    source_sha256: str | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "record_id", _identifier(self.record_id, field="record_id"))
        object.__setattr__(
            self, "component_id", _identifier(self.component_id, field="component_id")
        )
        object.__setattr__(self, "site_id", _identifier(self.site_id, field="site_id"))
        if self.site_role not in {"donor", "acceptor", "effective"}:
            raise _error("invalid-value", field="site_role", value=self.site_role)
        object.__setattr__(
            self,
            "multiplicity",
            _exact_integer(self.multiplicity, field="multiplicity", positive=True),
        )
        source_id, locator, domain_id = _provenance(self.source_id, self.locator, self.domain_id)
        object.__setattr__(self, "source_id", source_id)
        object.__setattr__(self, "locator", locator)
        object.__setattr__(self, "domain_id", domain_id)
        _normalize_candidate(self, field="association_site.candidate")

    @property
    def identity(self) -> str:
        return f"component/{self.component_id}/site/{self.site_id}"


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
    candidate_domain_id: str | None = None
    qualification: ModelQualification | None = None
    source_sha256: str | None = None

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
        _normalize_candidate(self, field="association_parameter.candidate")

    @property
    def identity(self) -> str:
        return (
            f"association/{self.component_id_a}/{self.site_id_a}/"
            f"{self.component_id_b}/{self.site_id_b}/{self.family}"
        )


@dataclass(frozen=True, slots=True)
class ExponentialTerm:
    amplitude: object
    exponent_coefficient: object
    amplitude_identity: str | None = None
    exponent_identity: str | None = None

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
        if self.amplitude_identity is not None:
            object.__setattr__(
                self, "amplitude_identity", _parameter_identity(self.amplitude_identity)
            )
        if self.exponent_identity is not None:
            object.__setattr__(
                self, "exponent_identity", _parameter_identity(self.exponent_identity)
            )


@dataclass(frozen=True, slots=True)
class TemperaturePolynomialTerm:
    power: int
    coefficient: object
    coefficient_identity: str | None = None

    def __post_init__(self) -> None:
        power = _exact_integer(self.power, field="power")
        if power < 1:
            raise _error("invalid-correlation", field="power", value=power)
        object.__setattr__(self, "power", power)
        object.__setattr__(
            self,
            "coefficient",
            _quantity(
                self.coefficient,
                f"1 / kelvin ** {power}",
                field="coefficient",
            ),
        )
        if self.coefficient_identity is not None:
            object.__setattr__(
                self,
                "coefficient_identity",
                _parameter_identity(self.coefficient_identity),
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
    candidate_domain_id: str | None = None
    qualification: ModelQualification | None = None
    source_sha256: str | None = None
    correlation_id: str | None = None

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
        _normalize_candidate(self, field="correlation.candidate")
        correlation_id = self.correlation_id or (
            f"component/{self.component_id}/{self.family}/constant"
        )
        object.__setattr__(self, "correlation_id", _parameter_identity(correlation_id))

    @property
    def identity(self) -> str:
        return self.correlation_id or ""


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
    candidate_domain_id: str | None = None
    qualification: ModelQualification | None = None
    source_sha256: str | None = None
    correlation_id: str | None = None

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
                amplitude_identity=term.amplitude_identity,
                exponent_identity=term.exponent_identity,
            )
            for term in terms
        )
        object.__setattr__(self, "terms", normalized_terms)
        source_id, locator, domain_id = _provenance(self.source_id, self.locator, self.domain_id)
        object.__setattr__(self, "source_id", source_id)
        object.__setattr__(self, "locator", locator)
        object.__setattr__(self, "domain_id", domain_id)
        _normalize_candidate(self, field="correlation.candidate")
        correlation_id = self.correlation_id or (
            f"component/{self.component_id}/{self.family}/temperature-exponential"
        )
        object.__setattr__(self, "correlation_id", _parameter_identity(correlation_id))

    @property
    def identity(self) -> str:
        return self.correlation_id or ""


@dataclass(frozen=True, slots=True)
class ConstantPlusPolynomialInTemperatureCorrelation:
    form: ClassVar[str] = "constant-plus-polynomial-in-temperature"

    record_id: str
    component_id: str
    family: Literal["relative_permittivity"]
    constant: object
    reference_temperature: object
    terms: tuple[TemperaturePolynomialTerm, ...]
    source_id: str
    locator: str
    domain_id: str
    candidate_domain_id: str | None = None
    qualification: ModelQualification | None = None
    source_sha256: str | None = None
    correlation_id: str | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "record_id", _identifier(self.record_id, field="record_id"))
        object.__setattr__(
            self, "component_id", _identifier(self.component_id, field="component_id")
        )
        family = _family(self.family, allowed=frozenset({"relative_permittivity"}))
        object.__setattr__(self, "family", family)
        object.__setattr__(self, "constant", _parameter_value(family, self.constant))
        object.__setattr__(
            self,
            "reference_temperature",
            _quantity(self.reference_temperature, "kelvin", field="reference_temperature"),
        )
        terms = tuple(self.terms)
        if not terms or not all(isinstance(term, TemperaturePolynomialTerm) for term in terms):
            raise _error("invalid-correlation", field="terms", value="invalid")
        if tuple(term.power for term in terms) != tuple(range(1, len(terms) + 1)):
            raise _error("invalid-correlation", field="terms", reason="nonconsecutive-powers")
        object.__setattr__(self, "terms", terms)
        source_id, locator, domain_id = _provenance(self.source_id, self.locator, self.domain_id)
        object.__setattr__(self, "source_id", source_id)
        object.__setattr__(self, "locator", locator)
        object.__setattr__(self, "domain_id", domain_id)
        _normalize_candidate(self, field="correlation.candidate")
        correlation_id = self.correlation_id or (
            f"component/{self.component_id}/{self.family}/{self.form}"
        )
        object.__setattr__(self, "correlation_id", _parameter_identity(correlation_id))

    @property
    def identity(self) -> str:
        return self.correlation_id or ""


@dataclass(frozen=True, slots=True)
class AssociationEdgeSource:
    """The single declared origin of one materialized association edge."""

    kind: AssociationEdgeSourceKind
    rule_id: str | None
    provenance: tuple[tuple[str, str, str], ...]
    input_identities: tuple[str, ...]
    derivative_sources: tuple[AssociationDerivativeSource, ...]

    def __post_init__(self) -> None:
        if self.kind not in {"explicit", "combining-rule"}:
            raise _error("invalid-value", field="association_edge.source.kind", value=self.kind)
        if self.kind == "combining-rule":
            if self.rule_id is None:
                raise _error(
                    "invalid-value", field="association_edge.source.rule_id", value="missing"
                )
            if self.rule_id not in _ASSOCIATION_COMBINING_RULES:
                raise _error(
                    "invalid-value", field="association_edge.source.rule_id", value=self.rule_id
                )
        elif self.rule_id is not None:
            raise _error(
                "invalid-schema", field="association_edge.source", reason="competing-sources"
            )
        if not self.provenance:
            raise _error("invalid-value", field="association_edge.source.provenance", value="empty")
        normalized = tuple(
            (
                _identifier(source_id, field="source_id"),
                _text(locator, field="locator"),
                _identifier(domain_id, field="domain_id"),
            )
            for source_id, locator, domain_id in self.provenance
        )
        object.__setattr__(self, "provenance", normalized)
        inputs = tuple(
            _parameter_identity(identity, field="association_edge.source.input_identity")
            for identity in self.input_identities
        )
        derivatives = tuple(self.derivative_sources)
        if not all(isinstance(item, AssociationDerivativeSource) for item in derivatives):
            raise _error("invalid-value", field="association_edge.source.derivative_sources")
        if self.kind == "explicit" and (inputs or derivatives):
            raise _error("invalid-schema", field="association_edge.source", reason="explicit-chain")
        if self.kind == "combining-rule" and (not inputs or not derivatives):
            raise _error("invalid-schema", field="association_edge.source", reason="missing-chain")
        if len(inputs) != len(set(inputs)):
            raise _error("duplicate-record", field="association_edge.source.input_identity")
        object.__setattr__(self, "input_identities", inputs)
        object.__setattr__(self, "derivative_sources", derivatives)


@dataclass(frozen=True, slots=True)
class AssociationDerivativeSource:
    """One exact local chain factor from a rule input to a materialized edge value."""

    output_identity: str
    input_identity: str
    factor: object

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "output_identity",
            _parameter_identity(self.output_identity, field="output_identity"),
        )
        object.__setattr__(
            self,
            "input_identity",
            _parameter_identity(self.input_identity, field="input_identity"),
        )
        object.__setattr__(self, "factor", _decimal_number(self.factor, field="factor"))


@dataclass(frozen=True, slots=True)
class AssociationEdge:
    endpoint_a: tuple[str, str]
    endpoint_b: tuple[str, str]
    energy_identity: str
    volume_identity: str
    energy_over_k: object
    volume: object
    source: AssociationEdgeSource
    candidate_domain_id: str | None = None
    qualification: ModelQualification | None = None
    source_sha256: str | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.source, AssociationEdgeSource):
            raise _error("invalid-value", field="association_edge.source", value="invalid")
        endpoint_a = (
            _identifier(self.endpoint_a[0], field="component_id_a"),
            _identifier(self.endpoint_a[1], field="site_id_a"),
        )
        endpoint_b = (
            _identifier(self.endpoint_b[0], field="component_id_b"),
            _identifier(self.endpoint_b[1], field="site_id_b"),
        )
        endpoint_a, endpoint_b = sorted((endpoint_a, endpoint_b))
        object.__setattr__(self, "endpoint_a", endpoint_a)
        object.__setattr__(self, "endpoint_b", endpoint_b)
        energy_identity = _parameter_identity(self.energy_identity, field="energy_identity")
        volume_identity = _parameter_identity(self.volume_identity, field="volume_identity")
        prefix = f"association/{endpoint_a[0]}/{endpoint_a[1]}/{endpoint_b[0]}/{endpoint_b[1]}"
        if energy_identity != f"{prefix}/energy_over_k":
            raise _error("invalid-value", field="energy_identity", value=energy_identity)
        if volume_identity != f"{prefix}/volume":
            raise _error("invalid-value", field="volume_identity", value=volume_identity)
        energy = _quantity(self.energy_over_k, "kelvin", field="energy_over_k")
        volume = _decimal_number(self.volume, field="volume")
        if energy.magnitude <= 0 or volume <= 0:
            raise _error("invalid-value", field="association_edge", reason="nonpositive")
        object.__setattr__(self, "energy_identity", energy_identity)
        object.__setattr__(self, "volume_identity", volume_identity)
        object.__setattr__(self, "energy_over_k", energy)
        object.__setattr__(self, "volume", volume)
        _normalize_candidate(self, field="association_edge.candidate")


@dataclass(frozen=True, slots=True)
class AssociationTopology:
    presets: tuple[tuple[str, AssociationPreset], ...]
    sites: tuple[SiteRecord, ...]
    edges: tuple[AssociationEdge, ...]

    def __post_init__(self) -> None:
        presets = tuple(
            (
                _identifier(component_id, field="preset.component_id"),
                preset,
            )
            for component_id, preset in self.presets
        )
        if any(preset not in _ASSOCIATION_PRESETS for _, preset in presets):
            raise _error("invalid-value", field="preset")
        if len({component_id for component_id, _ in presets}) != len(presets):
            raise _error("duplicate-record", field="preset.component_id")
        object.__setattr__(self, "presets", tuple(sorted(presets)))
        site_keys = tuple((site.component_id, site.site_id) for site in self.sites)
        if len(site_keys) != len(set(site_keys)):
            raise _error("duplicate-record", field="site")
        edge_keys = tuple((edge.endpoint_a, edge.endpoint_b) for edge in self.edges)
        if len(edge_keys) != len(set(edge_keys)):
            raise _error("duplicate-record", field="association_edge")
        site_set = set(site_keys)
        if any(
            endpoint not in site_set
            for edge in self.edges
            for endpoint in (edge.endpoint_a, edge.endpoint_b)
        ):
            raise _error("missing-record", field="association_edge.site")


type CorrelationRecord = (
    ConstantCorrelation
    | ConstantPlusSumOfExponentialsCorrelation
    | ConstantPlusPolynomialInTemperatureCorrelation
)
type ParameterRecord = (
    SingleParameterRecord
    | ModelParameterRecord
    | PairParameterRecord
    | CrossAssociationCorrectionRecord
    | SiteRecord
    | AssociationParameterRecord
    | CorrelationRecord
    | ReactionCorrelationRecord
)


__all__ = (
    "AssociationCombiningRuleId",
    "AssociationDerivativeSource",
    "AssociationEdge",
    "AssociationEdgeSource",
    "AssociationEdgeSourceKind",
    "AssociationParameterFamily",
    "AssociationParameterRecord",
    "AssociationPreset",
    "AssociationSiteRole",
    "AssociationTopology",
    "ComponentRecord",
    "ConstantCorrelation",
    "ConstantPlusPolynomialInTemperatureCorrelation",
    "ConstantPlusSumOfExponentialsCorrelation",
    "CorrelationFamily",
    "CorrelationRecord",
    "CrossAssociationCorrectionFamily",
    "CrossAssociationCorrectionRecord",
    "DomainKind",
    "ExponentialTerm",
    "IonicDispersionFormulation",
    "ModelFamilyKind",
    "ModelFamilySpec",
    "ModelParameterFamily",
    "ModelParameterRecord",
    "ModelQualification",
    "PairParameterFamily",
    "PairParameterRecord",
    "ParameterRecord",
    "ParameterSpec",
    "ReactionCorrelationRecord",
    "RelativePermittivityFormulation",
    "SingleParameterFamily",
    "SingleParameterRecord",
    "SiteRecord",
    "SourceRecord",
    "TemperaturePolynomialTerm",
    "ValidityDomain",
)
