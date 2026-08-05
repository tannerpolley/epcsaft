from __future__ import annotations

import hashlib
import json
import math
from collections.abc import Iterable
from dataclasses import dataclass, replace
from enum import StrEnum
from typing import Any

from epcsaft import regression as _package
from epcsaft.eos import Mixture, _native_transport

import epcsaft.regression._native as _native


class ParameterFamily(StrEnum):
    SEGMENT_COUNT = "segment_count"
    SEGMENT_DIAMETER = "segment_diameter"
    DISPERSION_ENERGY_OVER_K = "dispersion_energy_over_k"
    RELATIVE_PERMITTIVITY = "relative_permittivity"
    BORN_DIAMETER = "born_diameter"
    SOLVATION_FACTOR = "solvation_factor"
    SCHRECKENBERG_DIELECTRIC_VOLUME = "schreckenberg_dielectric_volume"
    SCHRECKENBERG_DIELECTRIC_TEMPERATURE = "schreckenberg_dielectric_temperature"
    ZUBER_ION_SUPPRESSION_COEFFICIENT = "zuber_ion_suppression_coefficient"
    RUEBEN_DIPOLE_SCALING = "rueben_dipole_scaling"
    RUEBEN_POLARIZABILITY_SCALING = "rueben_polarizability_scaling"
    RUEBEN_CORRELATION_INTEGRAL_PARAMETER = "rueben_correlation_integral_parameter"
    K_IJ = "k_ij"
    L_IJ = "l_ij"
    ASSOCIATION_ENERGY_OVER_K = "association_energy_over_k"
    ASSOCIATION_VOLUME = "association_volume"
    ION_FRACTION_SUPPRESSION_COEFFICIENT = "ion_fraction_suppression_coefficient"
    IONIC_REGION_RELATIVE_PERMITTIVITY = "ionic_region_relative_permittivity"


@dataclass(frozen=True, slots=True)
class ParameterCapability:
    capability_id: str
    family: ParameterFamily
    component_ids: tuple[str, ...]
    coordinate_kinds: tuple[str, ...]
    coordinate_units: tuple[str, ...]
    parameter_fingerprint: str
    topology_fingerprint: str
    derivative_order: int
    maturity: str
    authority_effect: str
    temperature_min_k: float
    temperature_max_k: float
    identity_shape: str
    observation_contract: str
    model_domain: str
    tensor_layout: str
    state_coordinate_count: int
    active_parameter_count: int
    helmholtz_basis_id: str
    unsupported_status: str
    domain_status: str
    active_component_ids: tuple[str, ...]

    @property
    def installed_ready(self) -> bool:
        required_order = (
            2 if self.observation_contract == "fixed_composition_helmholtz_phase" else 1
        )
        return self.maturity == "DERIVATIVE_READY" and self.derivative_order >= required_order

    @property
    def unsupported_reason(self) -> str:
        if self.installed_ready:
            return ""
        return f"installed maturity={self.maturity}, derivative order={self.derivative_order}"


@dataclass(frozen=True, slots=True)
class FixedTopologyParameterSlot:
    slot_index: int
    parameter_id: str
    family: ParameterFamily
    identity: ComponentParameterIdentity | AssociationParameterIdentity
    unit: str
    lower_bound_exclusive: float
    upper_bound_exclusive: float


@dataclass(frozen=True, slots=True)
class FixedTopologyAssociationCapability:
    component_ids: tuple[str, ...]
    site_ids: tuple[str, ...]
    site_component_ids: tuple[str, ...]
    site_multiplicities: tuple[int, ...]
    slots: tuple[FixedTopologyParameterSlot, ...]
    parameter_fingerprint: str
    topology_fingerprint: str
    artifact_fingerprint: str
    derivative_order: int
    maturity: str
    authority_effect: str
    temperature_min_k: float
    temperature_max_k: float
    observation_contract: str
    model_domain: str
    tensor_layout: str
    state_coordinate_count: int
    helmholtz_basis_id: str
    association_basis_id: str
    unsupported_status: str
    domain_status: str

    @property
    def installed_ready(self) -> bool:
        return self.maturity == "DERIVATIVE_READY" and self.derivative_order >= 2

    @property
    def unsupported_reason(self) -> str:
        if self.installed_ready:
            return ""
        return f"installed maturity={self.maturity}, derivative order={self.derivative_order}"


@dataclass(frozen=True, slots=True)
class UnsupportedParameterCapability:
    capability_code: int
    schema_version: int
    parameter_family_code: int

    @property
    def installed_ready(self) -> bool:
        return False

    @property
    def unsupported_reason(self) -> str:
        return (
            f"unsupported capability schema {self.schema_version}, "
            f"family code {self.parameter_family_code}"
        )


@dataclass(frozen=True, slots=True)
class FittedParameterDiagnostic:
    family: ParameterFamily
    component_ids: tuple[str, ...]
    unit: str
    transform_origin: float
    transform_scale: float
    start: float
    final: float
    movement: float
    lower_bound: float
    upper_bound: float
    active_bound_distance: float
    active_bound: str | None
    active_coordinate_ids: tuple[str, ...] = ()


@dataclass(frozen=True, slots=True)
class GeneralJacobianDiagnostics:
    residual_count: int
    variable_count: int
    full_singular_values: tuple[float, ...]
    full_rank: int
    full_condition_number: float
    projected_parameter_singular_values: tuple[float, ...]
    projected_parameter_rank: int
    projected_parameter_condition_number: float


@dataclass(frozen=True, slots=True)
class GeneralRowDiagnostic:
    row_id: str
    partition: str
    liquid_volume_m3_per_mol: float
    vapor_volume_m3_per_mol: float
    scaled_residuals: tuple[float, float, float, float]
    observed_pressure_pa: float
    liquid_model_pressure_pa: float
    vapor_model_pressure_pa: float
    chemical_potential_differences_over_rt: tuple[float, float]
    derivative_status: str
    status: str
    evaluated: bool
    failure_reason: str


@dataclass(frozen=True, slots=True)
class PureSaturationRowDiagnostic:
    row_id: str
    partition: str
    liquid_volume_m3_per_mol: float
    vapor_volume_m3_per_mol: float
    scaled_residuals: tuple[float, float, float, float]
    observed_pressure_pa: float
    liquid_model_pressure_pa: float
    vapor_model_pressure_pa: float
    chemical_potential_difference_over_rt: float
    observed_liquid_density_kg_per_m3: float
    model_liquid_density_kg_per_m3: float
    derivative_status: str
    status: str
    evaluated: bool
    failure_reason: str


@dataclass(frozen=True, slots=True)
class PureVaporPressureRowDiagnostic:
    row_id: str
    partition: str
    liquid_volume_m3_per_mol: float
    vapor_volume_m3_per_mol: float
    scaled_residuals: tuple[float, float, float]
    observed_pressure_pa: float
    liquid_model_pressure_pa: float
    vapor_model_pressure_pa: float
    chemical_potential_difference_over_rt: float
    derivative_status: str
    status: str
    evaluated: bool
    failure_reason: str


@dataclass(frozen=True, slots=True)
class DensityRowDiagnostic:
    row_id: str
    partition: str
    volume_m3_per_mol: float
    scaled_residuals: tuple[float, float]
    observed_pressure_pa: float
    model_pressure_pa: float
    observed_density_kg_per_m3: float
    model_density_kg_per_m3: float
    derivative_status: str
    status: str
    evaluated: bool
    failure_reason: str


@dataclass(frozen=True, slots=True)
class DirectObservationRowDiagnostic:
    row_id: str
    partition: str
    observable: str
    observable_unit: str
    observed_value: float
    modeled_value: float
    scaled_residual: float
    provider_derivative: float
    derivative_status: str
    status: str
    evaluated: bool
    failure_reason: str


def parameter_capabilities(
    model: Mixture,
) -> tuple[
    ParameterCapability | FixedTopologyAssociationCapability | UnsupportedParameterCapability,
    ...,
]:
    raw_capabilities = _native.parameter_capabilities(_native_transport(model))
    capabilities: list[
        ParameterCapability | FixedTopologyAssociationCapability | UnsupportedParameterCapability
    ] = []
    for raw in raw_capabilities:
        if raw[0] == "fixed_topology_association":
            component_ids = tuple(raw[1])
            site_ids = tuple(raw[2])
            site_component_indices = tuple(raw[3])
            site_component_ids = tuple(component_ids[index] for index in site_component_indices)
            slots = []
            for slot_index, item in enumerate(raw[5]):
                family = ParameterFamily(item[0])
                if item[2] < 0:
                    identity: ComponentParameterIdentity | AssociationParameterIdentity = (
                        ComponentParameterIdentity(component_ids[item[1]])
                    )
                else:
                    identity = AssociationParameterIdentity(
                        family,
                        (
                            AssociationSitePairIdentity(
                                component_ids[item[4]],
                                site_ids[item[2]],
                                component_ids[item[5]],
                                site_ids[item[3]],
                            ),
                        ),
                    )
                slots.append(
                    FixedTopologyParameterSlot(
                        slot_index=slot_index,
                        parameter_id=item[6],
                        family=family,
                        identity=identity,
                        unit=_provider_unit(item[7]),
                        lower_bound_exclusive=item[8],
                        upper_bound_exclusive=item[9],
                    )
                )
            capabilities.append(
                FixedTopologyAssociationCapability(
                    component_ids=component_ids,
                    site_ids=site_ids,
                    site_component_ids=site_component_ids,
                    site_multiplicities=tuple(raw[4]),
                    slots=tuple(slots),
                    parameter_fingerprint=raw[6],
                    topology_fingerprint=raw[7],
                    artifact_fingerprint=raw[8],
                    derivative_order=raw[9],
                    maturity=raw[10],
                    authority_effect=raw[11],
                    temperature_min_k=raw[12],
                    temperature_max_k=raw[13],
                    observation_contract=raw[14],
                    model_domain=raw[15],
                    tensor_layout=raw[16],
                    state_coordinate_count=raw[17],
                    helmholtz_basis_id=raw[18],
                    association_basis_id=raw[19],
                    unsupported_status=raw[20],
                    domain_status=raw[21],
                )
            )
            continue
        if len(raw) == 3 and type(raw[0]) is int:
            capabilities.append(
                UnsupportedParameterCapability(
                    capability_code=raw[0],
                    schema_version=raw[1],
                    parameter_family_code=raw[2],
                )
            )
            continue
        capabilities.append(
            ParameterCapability(
                capability_id=raw[0],
                family=ParameterFamily(raw[1]),
                component_ids=tuple(raw[2]),
                coordinate_kinds=tuple(raw[3]),
                coordinate_units=tuple(raw[4]),
                parameter_fingerprint=raw[5],
                topology_fingerprint=raw[6],
                derivative_order=raw[7],
                maturity=raw[8],
                authority_effect=raw[9],
                temperature_min_k=raw[10],
                temperature_max_k=raw[11],
                identity_shape=raw[12],
                observation_contract=raw[13],
                model_domain=raw[14],
                tensor_layout=raw[15],
                state_coordinate_count=raw[16],
                active_parameter_count=raw[17],
                helmholtz_basis_id=raw[18],
                unsupported_status=raw[19],
                domain_status=raw[20],
                active_component_ids=(
                    tuple(sorted(raw[21]))
                    if raw[12] == "unordered_component_pair"
                    else tuple(raw[21])
                ),
            )
        )
    return tuple(capabilities)


def _provider_unit(value: str) -> str:
    units = {"dimensionless": "1", "angstrom": "angstrom", "kelvin": "K"}
    try:
        return units[value]
    except KeyError as error:
        raise ValueError(f"installed capability uses unsupported unit {value!r}") from error


class ObservationPartition(StrEnum):
    TRAINING = "training"
    HELD_OUT = "held_out"
    STRESS = "stress"


def _require_nonempty_string(value: str, field: str) -> None:
    if type(value) is not str or not value.strip():
        raise ValueError(f"{field} must be a nonempty string")


def _require_finite(value: float, field: str, *, positive: bool = False) -> None:
    if type(value) not in (int, float) or not math.isfinite(value):
        raise ValueError(f"{field} must be finite")
    if positive and value <= 0.0:
        raise ValueError(f"{field} must be positive")


def _require_sha256(value: str, field: str, *, prefixed: bool = False) -> None:
    _require_nonempty_string(value, field)
    prefix = "sha256:" if prefixed else ""
    body = value.removeprefix(prefix)
    if prefixed and not value.startswith(prefix):
        raise ValueError(f"{field} must begin with 'sha256:'")
    if len(body) != 64 or any(character not in "0123456789abcdef" for character in body):
        raise ValueError(f"{field} must be a lowercase SHA-256 identity")


@dataclass(frozen=True, slots=True)
class PairParameterIdentity:
    component_id_a: str
    component_id_b: str

    def __post_init__(self) -> None:
        _require_nonempty_string(self.component_id_a, "component_id_a")
        _require_nonempty_string(self.component_id_b, "component_id_b")
        if self.component_id_a == self.component_id_b:
            raise ValueError("pair parameter components must be distinct")
        canonical = tuple(sorted((self.component_id_a, self.component_id_b)))
        object.__setattr__(self, "component_id_a", canonical[0])
        object.__setattr__(self, "component_id_b", canonical[1])

    @property
    def component_ids(self) -> tuple[str, str]:
        return (self.component_id_a, self.component_id_b)

    @property
    def canonical_component_ids(self) -> tuple[str, str]:
        return self.component_ids


@dataclass(frozen=True, slots=True)
class ComponentParameterIdentity:
    component_id: str

    def __post_init__(self) -> None:
        _require_nonempty_string(self.component_id, "component_id")

    @property
    def component_ids(self) -> tuple[str]:
        return (self.component_id,)

    @property
    def canonical_component_ids(self) -> tuple[str]:
        return self.component_ids


@dataclass(frozen=True, slots=True)
class ModelParameterIdentity:
    @property
    def component_ids(self) -> tuple[()]:
        return ()

    @property
    def canonical_component_ids(self) -> tuple[()]:
        return ()


@dataclass(frozen=True, slots=True, order=True)
class AssociationSitePairIdentity:
    component_id_a: str
    site_id_a: str
    component_id_b: str
    site_id_b: str

    def __post_init__(self) -> None:
        for field in (
            "component_id_a",
            "site_id_a",
            "component_id_b",
            "site_id_b",
        ):
            _require_nonempty_string(getattr(self, field), field)
        endpoints = sorted(
            (
                (self.component_id_a, self.site_id_a),
                (self.component_id_b, self.site_id_b),
            )
        )
        object.__setattr__(self, "component_id_a", endpoints[0][0])
        object.__setattr__(self, "site_id_a", endpoints[0][1])
        object.__setattr__(self, "component_id_b", endpoints[1][0])
        object.__setattr__(self, "site_id_b", endpoints[1][1])

    @property
    def component_ids(self) -> tuple[str, ...]:
        return tuple(sorted({self.component_id_a, self.component_id_b}))

    @property
    def canonical_id(self) -> str:
        return f"{self.component_id_a}/{self.site_id_a}--{self.component_id_b}/{self.site_id_b}"


@dataclass(frozen=True, slots=True)
class AssociationParameterIdentity:
    quantity: ParameterFamily
    site_pairs: tuple[AssociationSitePairIdentity, ...]

    def __post_init__(self) -> None:
        if self.quantity not in (
            ParameterFamily.ASSOCIATION_ENERGY_OVER_K,
            ParameterFamily.ASSOCIATION_VOLUME,
        ):
            raise ValueError("association parameter quantity must be association energy or volume")
        if (
            type(self.site_pairs) is not tuple
            or not self.site_pairs
            or any(not isinstance(pair, AssociationSitePairIdentity) for pair in self.site_pairs)
        ):
            raise TypeError(
                "association parameter site_pairs must be a nonempty tuple of "
                "AssociationSitePairIdentity values"
            )
        canonical = tuple(sorted(self.site_pairs))
        if len(set(canonical)) != len(canonical):
            raise ValueError("association parameter site pairs must be unique")
        object.__setattr__(self, "site_pairs", canonical)

    @property
    def component_ids(self) -> tuple[str, ...]:
        return tuple(
            sorted({component for pair in self.site_pairs for component in pair.component_ids})
        )

    @property
    def canonical_component_ids(self) -> tuple[str, ...]:
        return self.component_ids

    @property
    def canonical_ids(self) -> tuple[str, ...]:
        return tuple(pair.canonical_id for pair in self.site_pairs)


@dataclass(frozen=True, slots=True)
class AffineParameterTransform:
    origin: float
    scale: float

    def __post_init__(self) -> None:
        _require_finite(self.origin, "transform origin")
        _require_finite(self.scale, "transform scale")
        if self.scale == 0.0:
            raise ValueError("transform scale must be nonzero")

    def to_solver(self, physical_value: float) -> float:
        _require_finite(physical_value, "physical parameter value")
        return (physical_value - self.origin) / self.scale

    def to_physical(self, solver_value: float) -> float:
        _require_finite(solver_value, "solver parameter value")
        return self.origin + self.scale * solver_value


@dataclass(frozen=True, slots=True)
class ParameterCoordinate:
    family: ParameterFamily
    identity: (
        PairParameterIdentity
        | ComponentParameterIdentity
        | ModelParameterIdentity
        | AssociationParameterIdentity
    )
    capability_id: str | None
    provider_parameter_fingerprint: str
    provider_topology_fingerprint: str
    unit: str
    transform: AffineParameterTransform
    lower_bound: float
    upper_bound: float
    provider_artifact_fingerprint: str | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.family, ParameterFamily):
            raise TypeError("family must be a ParameterFamily")
        pair_families = (ParameterFamily.K_IJ, ParameterFamily.L_IJ)
        component_units = {
            ParameterFamily.SEGMENT_COUNT: "1",
            ParameterFamily.SEGMENT_DIAMETER: "angstrom",
            ParameterFamily.DISPERSION_ENERGY_OVER_K: "K",
            ParameterFamily.RELATIVE_PERMITTIVITY: "1",
            ParameterFamily.BORN_DIAMETER: "angstrom",
            ParameterFamily.SOLVATION_FACTOR: "1",
        }
        model_units = {
            ParameterFamily.ION_FRACTION_SUPPRESSION_COEFFICIENT: "1",
            ParameterFamily.IONIC_REGION_RELATIVE_PERMITTIVITY: "1",
        }
        if self.family in pair_families:
            if not isinstance(self.identity, PairParameterIdentity):
                raise TypeError("pair-parameter identity must be a PairParameterIdentity")
            expected_unit = "1"
        elif self.family in component_units:
            if not isinstance(self.identity, ComponentParameterIdentity):
                raise TypeError("component-parameter identity must be a ComponentParameterIdentity")
            expected_unit = component_units[self.family]
        elif self.family in (
            ParameterFamily.ASSOCIATION_ENERGY_OVER_K,
            ParameterFamily.ASSOCIATION_VOLUME,
        ):
            if (
                not isinstance(self.identity, AssociationParameterIdentity)
                or self.identity.quantity is not self.family
            ):
                raise TypeError(
                    "association coordinates require a matching AssociationParameterIdentity"
                )
            expected_unit = "K" if self.family is ParameterFamily.ASSOCIATION_ENERGY_OVER_K else "1"
        elif self.family in model_units:
            if not isinstance(self.identity, ModelParameterIdentity):
                raise TypeError("model-parameter identity must be a ModelParameterIdentity")
            expected_unit = model_units[self.family]
        else:
            raise ValueError(
                "the v1 coordinate contract supports only k_ij, l_ij, "
                "segment_count, segment_diameter, and "
                "dispersion_energy_over_k, relative_permittivity, "
                "born_diameter, or "
                "solvation_factor, or ion_fraction_suppression_coefficient"
                ", ionic_region_relative_permittivity, "
                "association_energy_over_k, or association_volume"
            )
        if isinstance(self.identity, AssociationParameterIdentity):
            if self.capability_id is not None:
                raise ValueError("generic association coordinates do not use capability IDs")
            if self.provider_artifact_fingerprint is None:
                raise ValueError(
                    "generic association coordinates require the installed artifact fingerprint"
                )
            _require_sha256(
                self.provider_artifact_fingerprint,
                "provider artifact fingerprint",
                prefixed=True,
            )
        else:
            if self.capability_id is None:
                if self.provider_artifact_fingerprint is None:
                    raise ValueError(
                        "descriptor-driven coordinates require the installed artifact fingerprint"
                    )
                _require_sha256(
                    self.provider_artifact_fingerprint,
                    "provider artifact fingerprint",
                    prefixed=True,
                )
            else:
                _require_nonempty_string(self.capability_id, "capability_id")
                if self.provider_artifact_fingerprint is not None:
                    _require_sha256(
                        self.provider_artifact_fingerprint,
                        "provider artifact fingerprint",
                        prefixed=True,
                    )
        _require_sha256(
            self.provider_parameter_fingerprint,
            "provider parameter fingerprint",
            prefixed=True,
        )
        _require_sha256(
            self.provider_topology_fingerprint,
            "provider topology fingerprint",
            prefixed=True,
        )
        _require_nonempty_string(self.unit, "unit")
        if self.unit != expected_unit:
            raise ValueError(f"{self.family.value} parameter unit must be {expected_unit!r}")
        if not isinstance(self.transform, AffineParameterTransform):
            raise TypeError("transform must be an AffineParameterTransform")
        _require_finite(self.lower_bound, "lower bound")
        _require_finite(self.upper_bound, "upper bound")
        if self.lower_bound >= self.upper_bound:
            raise ValueError("parameter bounds must be strictly increasing")


@dataclass(frozen=True, slots=True)
class SourceDescriptor:
    source_id: str
    citation: str
    durable_locator: str
    source_artifact_sha256: str
    canonical_dataset_sha256: str
    transformation_record: str
    units_and_bases: str
    use_basis: str
    residual_scale_rationale: str

    def __post_init__(self) -> None:
        for field in (
            "source_id",
            "citation",
            "durable_locator",
            "transformation_record",
            "units_and_bases",
            "use_basis",
            "residual_scale_rationale",
        ):
            _require_nonempty_string(getattr(self, field), field)
        _require_sha256(self.source_artifact_sha256, "source artifact SHA-256")
        _require_sha256(self.canonical_dataset_sha256, "canonical dataset SHA-256")


@dataclass(frozen=True, slots=True)
class FixedCompositionVleObservation:
    row_id: str
    source_id: str
    source_locator: str
    component_ids: tuple[str, str]
    temperature_k: float
    pressure_pa: float
    liquid_mole_fraction_first: float
    vapor_mole_fraction_first: float
    pressure_scale_pa: float
    chemical_potential_scales: tuple[float, float]
    liquid_volume_origin_m3_per_mol: float
    liquid_volume_start_m3_per_mol: float
    liquid_volume_bounds_m3_per_mol: tuple[float, float]
    vapor_volume_origin_m3_per_mol: float
    vapor_volume_start_m3_per_mol: float
    vapor_volume_bounds_m3_per_mol: tuple[float, float]
    partition: ObservationPartition

    def __post_init__(self) -> None:
        for field in ("row_id", "source_id", "source_locator"):
            _require_nonempty_string(getattr(self, field), field)
        if type(self.component_ids) is not tuple or len(self.component_ids) != 2:
            raise ValueError("component_ids must contain exactly two component identifiers")
        for component_id in self.component_ids:
            _require_nonempty_string(component_id, "component_id")
        if self.component_ids[0] == self.component_ids[1]:
            raise ValueError("component_ids must be distinct")
        _require_finite(self.temperature_k, "temperature_k", positive=True)
        _require_finite(self.pressure_pa, "pressure_pa", positive=True)
        for value, field in (
            (self.liquid_mole_fraction_first, "liquid mole fraction"),
            (self.vapor_mole_fraction_first, "vapor mole fraction"),
        ):
            _require_finite(value, field)
            if not 0.0 < value < 1.0:
                raise ValueError(f"{field} must be strictly between zero and one")
        _require_finite(self.pressure_scale_pa, "pressure_scale_pa", positive=True)
        if (
            type(self.chemical_potential_scales) is not tuple
            or len(self.chemical_potential_scales) != 2
        ):
            raise ValueError("chemical_potential_scales must contain two values")
        for scale in self.chemical_potential_scales:
            _require_finite(scale, "chemical potential scale", positive=True)
        self._validate_volume("liquid")
        self._validate_volume("vapor")
        if not isinstance(self.partition, ObservationPartition):
            raise TypeError("partition must be an ObservationPartition")

    def _validate_volume(self, phase: str) -> None:
        origin = getattr(self, f"{phase}_volume_origin_m3_per_mol")
        start = getattr(self, f"{phase}_volume_start_m3_per_mol")
        bounds = getattr(self, f"{phase}_volume_bounds_m3_per_mol")
        _require_finite(origin, f"{phase} volume origin", positive=True)
        _require_finite(start, f"{phase} volume start", positive=True)
        if type(bounds) is not tuple or len(bounds) != 2:
            raise ValueError(f"{phase} volume bounds must contain two values")
        lower, upper = bounds
        _require_finite(lower, f"{phase} volume lower bound", positive=True)
        _require_finite(upper, f"{phase} volume upper bound", positive=True)
        if lower >= upper:
            raise ValueError(f"{phase} volume bounds must be strictly increasing")
        if not lower <= start <= upper:
            raise ValueError(f"{phase} volume start must lie within its bounds")


@dataclass(frozen=True, slots=True)
class PureSaturationObservation:
    row_id: str
    source_id: str
    source_locator: str
    component_id: str
    temperature_k: float
    pressure_pa: float
    liquid_density_kg_per_m3: float
    molar_mass_kg_per_mol: float
    pressure_scale_pa: float
    chemical_potential_scale: float
    liquid_density_scale_kg_per_m3: float
    liquid_volume_origin_m3_per_mol: float
    liquid_volume_start_m3_per_mol: float
    liquid_volume_bounds_m3_per_mol: tuple[float, float]
    vapor_volume_origin_m3_per_mol: float
    vapor_volume_start_m3_per_mol: float
    vapor_volume_bounds_m3_per_mol: tuple[float, float]
    partition: ObservationPartition

    def __post_init__(self) -> None:
        for field in ("row_id", "source_id", "source_locator", "component_id"):
            _require_nonempty_string(getattr(self, field), field)
        for field in (
            "temperature_k",
            "pressure_pa",
            "liquid_density_kg_per_m3",
            "molar_mass_kg_per_mol",
            "pressure_scale_pa",
            "chemical_potential_scale",
            "liquid_density_scale_kg_per_m3",
        ):
            _require_finite(getattr(self, field), field, positive=True)
        self._validate_volume("liquid")
        self._validate_volume("vapor")
        if not isinstance(self.partition, ObservationPartition):
            raise TypeError("partition must be an ObservationPartition")

    @property
    def component_ids(self) -> tuple[str]:
        return (self.component_id,)

    def _validate_volume(self, phase: str) -> None:
        origin = getattr(self, f"{phase}_volume_origin_m3_per_mol")
        start = getattr(self, f"{phase}_volume_start_m3_per_mol")
        bounds = getattr(self, f"{phase}_volume_bounds_m3_per_mol")
        _require_finite(origin, f"{phase} volume origin", positive=True)
        _require_finite(start, f"{phase} volume start", positive=True)
        if type(bounds) is not tuple or len(bounds) != 2:
            raise ValueError(f"{phase} volume bounds must contain two values")
        lower, upper = bounds
        _require_finite(lower, f"{phase} volume lower bound", positive=True)
        _require_finite(upper, f"{phase} volume upper bound", positive=True)
        if lower >= upper:
            raise ValueError(f"{phase} volume bounds must be strictly increasing")
        if not lower <= start <= upper:
            raise ValueError(f"{phase} volume start must lie within its bounds")


@dataclass(frozen=True, slots=True)
class PureVaporPressureObservation:
    row_id: str
    source_id: str
    source_locator: str
    component_id: str
    temperature_k: float
    pressure_pa: float
    pressure_scale_pa: float
    chemical_potential_scale: float
    liquid_volume_origin_m3_per_mol: float
    liquid_volume_start_m3_per_mol: float
    liquid_volume_bounds_m3_per_mol: tuple[float, float]
    vapor_volume_origin_m3_per_mol: float
    vapor_volume_start_m3_per_mol: float
    vapor_volume_bounds_m3_per_mol: tuple[float, float]
    partition: ObservationPartition

    def __post_init__(self) -> None:
        for field in ("row_id", "source_id", "source_locator", "component_id"):
            _require_nonempty_string(getattr(self, field), field)
        for field in (
            "temperature_k",
            "pressure_pa",
            "pressure_scale_pa",
            "chemical_potential_scale",
        ):
            _require_finite(getattr(self, field), field, positive=True)
        self._validate_volume("liquid")
        self._validate_volume("vapor")
        if self.liquid_volume_bounds_m3_per_mol[1] >= self.vapor_volume_bounds_m3_per_mol[0]:
            raise ValueError(
                "pure-vapor-pressure liquid and vapor volume bounds must be strictly separated"
            )
        if not isinstance(self.partition, ObservationPartition):
            raise TypeError("partition must be an ObservationPartition")

    @property
    def component_ids(self) -> tuple[str]:
        return (self.component_id,)

    def _validate_volume(self, phase: str) -> None:
        origin = getattr(self, f"{phase}_volume_origin_m3_per_mol")
        start = getattr(self, f"{phase}_volume_start_m3_per_mol")
        bounds = getattr(self, f"{phase}_volume_bounds_m3_per_mol")
        _require_finite(origin, f"{phase} volume origin", positive=True)
        _require_finite(start, f"{phase} volume start", positive=True)
        if type(bounds) is not tuple or len(bounds) != 2:
            raise ValueError(f"{phase} volume bounds must contain two values")
        lower, upper = bounds
        _require_finite(lower, f"{phase} volume lower bound", positive=True)
        _require_finite(upper, f"{phase} volume upper bound", positive=True)
        if lower >= upper:
            raise ValueError(f"{phase} volume bounds must be strictly increasing")
        if not lower <= start <= upper:
            raise ValueError(f"{phase} volume start must lie within its bounds")


@dataclass(frozen=True, slots=True)
class PureDensityObservation:
    row_id: str
    source_id: str
    source_locator: str
    component_id: str
    temperature_k: float
    pressure_pa: float
    density_kg_per_m3: float
    molar_mass_kg_per_mol: float
    pressure_scale_pa: float
    density_scale_kg_per_m3: float
    volume_origin_m3_per_mol: float
    volume_start_m3_per_mol: float
    volume_bounds_m3_per_mol: tuple[float, float]
    partition: ObservationPartition

    def __post_init__(self) -> None:
        for field in ("row_id", "source_id", "source_locator", "component_id"):
            _require_nonempty_string(getattr(self, field), field)
        for field in (
            "temperature_k",
            "pressure_pa",
            "density_kg_per_m3",
            "molar_mass_kg_per_mol",
            "pressure_scale_pa",
            "density_scale_kg_per_m3",
            "volume_origin_m3_per_mol",
            "volume_start_m3_per_mol",
        ):
            _require_finite(getattr(self, field), field, positive=True)
        if (
            type(self.volume_bounds_m3_per_mol) is not tuple
            or len(self.volume_bounds_m3_per_mol) != 2
        ):
            raise ValueError("volume bounds must contain two values")
        lower, upper = self.volume_bounds_m3_per_mol
        _require_finite(lower, "volume lower bound", positive=True)
        _require_finite(upper, "volume upper bound", positive=True)
        if lower >= upper:
            raise ValueError("volume bounds must be strictly increasing")
        if not lower <= self.volume_start_m3_per_mol <= upper:
            raise ValueError("volume start must lie within its bounds")
        if not isinstance(self.partition, ObservationPartition):
            raise TypeError("partition must be an ObservationPartition")

    @property
    def component_ids(self) -> tuple[str]:
        return (self.component_id,)


@dataclass(frozen=True, slots=True)
class FixedCompositionMixtureDensityObservation:
    row_id: str
    source_id: str
    source_locator: str
    component_ids: tuple[str, ...]
    mole_fractions: tuple[float, ...]
    temperature_k: float
    pressure_pa: float
    density_kg_per_m3: float
    molar_mass_kg_per_mol: float
    pressure_scale_pa: float
    density_scale_kg_per_m3: float
    volume_origin_m3_per_mol: float
    volume_start_m3_per_mol: float
    volume_bounds_m3_per_mol: tuple[float, float]
    partition: ObservationPartition

    def __post_init__(self) -> None:
        for field in ("row_id", "source_id", "source_locator"):
            _require_nonempty_string(getattr(self, field), field)
        if type(self.component_ids) is not tuple or len(self.component_ids) < 2:
            raise ValueError(
                "component_ids must contain at least two ordered component identifiers"
            )
        for component_id in self.component_ids:
            _require_nonempty_string(component_id, "component_id")
        if len(set(self.component_ids)) != len(self.component_ids):
            raise ValueError("component_ids must be distinct")
        if type(self.mole_fractions) is not tuple or len(self.mole_fractions) != len(
            self.component_ids
        ):
            raise ValueError("mole_fractions must match component_ids")
        for fraction in self.mole_fractions:
            _require_finite(fraction, "mole fraction")
            if fraction <= 0.0:
                raise ValueError("mole fractions must be positive")
        if not math.isclose(sum(self.mole_fractions), 1.0, rel_tol=0.0, abs_tol=1.0e-12):
            raise ValueError("mole_fractions must sum to one")
        for field in (
            "temperature_k",
            "pressure_pa",
            "density_kg_per_m3",
            "molar_mass_kg_per_mol",
            "pressure_scale_pa",
            "density_scale_kg_per_m3",
            "volume_origin_m3_per_mol",
            "volume_start_m3_per_mol",
        ):
            _require_finite(getattr(self, field), field, positive=True)
        if (
            type(self.volume_bounds_m3_per_mol) is not tuple
            or len(self.volume_bounds_m3_per_mol) != 2
        ):
            raise ValueError("volume bounds must contain two values")
        lower, upper = self.volume_bounds_m3_per_mol
        _require_finite(lower, "volume lower bound", positive=True)
        _require_finite(upper, "volume upper bound", positive=True)
        if lower >= upper:
            raise ValueError("volume bounds must be strictly increasing")
        if not lower <= self.volume_start_m3_per_mol <= upper:
            raise ValueError("volume start must lie within its bounds")
        if not isinstance(self.partition, ObservationPartition):
            raise TypeError("partition must be an ObservationPartition")


def _validate_direct_component_identity(
    component_ids: tuple[str, ...], active_component_id: str
) -> None:
    if type(component_ids) is not tuple or len(component_ids) < 2:
        raise ValueError("component_ids must contain the ordered Provider model")
    for component_id in component_ids:
        _require_nonempty_string(component_id, "component_id")
    if len(set(component_ids)) != len(component_ids):
        raise ValueError("component_ids must be distinct")
    _require_nonempty_string(active_component_id, "active_component_id")
    if active_component_id not in component_ids:
        raise ValueError("active_component_id must occur in component_ids")


@dataclass(frozen=True, slots=True)
class MeanIonicActivityObservation:
    row_id: str
    source_id: str
    source_locator: str
    component_ids: tuple[str, ...]
    active_component_id: str
    temperature_k: float
    pressure_pa: float
    formula_unit_molality_mol_per_kg: float
    observed_mean_ionic_activity_coefficient: float
    relative_residual_scale: float
    partition: ObservationPartition

    def __post_init__(self) -> None:
        for field in ("row_id", "source_id", "source_locator"):
            _require_nonempty_string(getattr(self, field), field)
        _validate_direct_component_identity(self.component_ids, self.active_component_id)
        for field in (
            "temperature_k",
            "pressure_pa",
            "formula_unit_molality_mol_per_kg",
            "observed_mean_ionic_activity_coefficient",
            "relative_residual_scale",
        ):
            _require_finite(getattr(self, field), field, positive=True)
        if not isinstance(self.partition, ObservationPartition):
            raise TypeError("partition must be an ObservationPartition")


@dataclass(frozen=True, slots=True)
class AqueousKijMeanIonicActivityObservation:
    row_id: str
    source_id: str
    source_locator: str
    component_ids: tuple[str, str, str]
    active_pair_component_ids: tuple[str, str]
    fixed_k_ij: tuple[float, float, float]
    temperature_k: float
    pressure_pa: float
    formula_unit_molality_mol_per_kg: float
    observed_mean_ionic_activity_coefficient: float
    relative_residual_scale: float
    partition: ObservationPartition

    def __post_init__(self) -> None:
        for field in ("row_id", "source_id", "source_locator"):
            _require_nonempty_string(getattr(self, field), field)
        if type(self.component_ids) is not tuple or len(self.component_ids) != 3:
            raise ValueError("component_ids must contain solvent, cation, and anion")
        for component_id in self.component_ids:
            _require_nonempty_string(component_id, "component_id")
        if len(set(self.component_ids)) != 3:
            raise ValueError("component_ids must be distinct")
        if (
            type(self.active_pair_component_ids) is not tuple
            or len(self.active_pair_component_ids) != 2
            or len(set(self.active_pair_component_ids)) != 2
            or not set(self.active_pair_component_ids).issubset(self.component_ids)
        ):
            raise ValueError("active_pair_component_ids must identify two model components")
        if type(self.fixed_k_ij) is not tuple or len(self.fixed_k_ij) != 3:
            raise ValueError(
                "fixed_k_ij must contain water-cation, water-anion, and cation-anion values"
            )
        for value in self.fixed_k_ij:
            _require_finite(value, "fixed k_ij")
        for field in (
            "temperature_k",
            "pressure_pa",
            "formula_unit_molality_mol_per_kg",
            "observed_mean_ionic_activity_coefficient",
            "relative_residual_scale",
        ):
            _require_finite(getattr(self, field), field, positive=True)
        if not isinstance(self.partition, ObservationPartition):
            raise TypeError("partition must be an ObservationPartition")

    @property
    def canonical_active_pair_component_ids(self) -> tuple[str, str]:
        ordered = sorted(self.active_pair_component_ids)
        return ordered[0], ordered[1]


@dataclass(frozen=True, slots=True)
class IonSolvationKijObservation:
    row_id: str
    source_id: str
    source_locator: str
    component_ids: tuple[str, str, str]
    active_component_id: str
    active_pair_component_ids: tuple[str, str]
    fixed_k_ij: tuple[float, float, float]
    temperature_k: float
    pressure_pa: float
    observed_solvation_gibbs_j_per_mol: float
    residual_scale_j_per_mol: float
    partition: ObservationPartition

    def __post_init__(self) -> None:
        for field in ("row_id", "source_id", "source_locator"):
            _require_nonempty_string(getattr(self, field), field)
        _validate_direct_component_identity(self.component_ids, self.active_component_id)
        if len(self.component_ids) != 3:
            raise ValueError("component_ids must contain solvent, cation, and anion")
        if (
            type(self.active_pair_component_ids) is not tuple
            or len(self.active_pair_component_ids) != 2
            or len(set(self.active_pair_component_ids)) != 2
            or not set(self.active_pair_component_ids).issubset(self.component_ids)
        ):
            raise ValueError("active_pair_component_ids must identify two model components")
        if self.active_component_id not in self.active_pair_component_ids:
            raise ValueError("active k_ij pair must contain the active ion")
        if self.active_component_id == self.component_ids[0]:
            raise ValueError("active_component_id must identify the cation or anion")
        if type(self.fixed_k_ij) is not tuple or len(self.fixed_k_ij) != 3:
            raise ValueError(
                "fixed_k_ij must contain solvent-cation, solvent-anion, and cation-anion values"
            )
        for value in self.fixed_k_ij:
            _require_finite(value, "fixed k_ij")
        _require_finite(self.temperature_k, "temperature_k", positive=True)
        _require_finite(self.pressure_pa, "pressure_pa", positive=True)
        _require_finite(
            self.observed_solvation_gibbs_j_per_mol,
            "observed_solvation_gibbs_j_per_mol",
        )
        _require_finite(
            self.residual_scale_j_per_mol,
            "residual_scale_j_per_mol",
            positive=True,
        )
        if not isinstance(self.partition, ObservationPartition):
            raise TypeError("partition must be an ObservationPartition")

    @property
    def canonical_active_pair_component_ids(self) -> tuple[str, str]:
        ordered = sorted(self.active_pair_component_ids)
        return ordered[0], ordered[1]


@dataclass(frozen=True, slots=True)
class SolvationGibbsObservation:
    row_id: str
    source_id: str
    source_locator: str
    component_ids: tuple[str, ...]
    active_component_id: str
    temperature_k: float
    pressure_pa: float
    observed_solvation_gibbs_j_per_mol: float
    residual_scale_j_per_mol: float
    partition: ObservationPartition

    def __post_init__(self) -> None:
        for field in ("row_id", "source_id", "source_locator"):
            _require_nonempty_string(getattr(self, field), field)
        _validate_direct_component_identity(self.component_ids, self.active_component_id)
        _require_finite(self.temperature_k, "temperature_k", positive=True)
        _require_finite(self.pressure_pa, "pressure_pa", positive=True)
        _require_finite(
            self.observed_solvation_gibbs_j_per_mol,
            "observed_solvation_gibbs_j_per_mol",
        )
        _require_finite(
            self.residual_scale_j_per_mol,
            "residual_scale_j_per_mol",
            positive=True,
        )
        if not isinstance(self.partition, ObservationPartition):
            raise TypeError("partition must be an ObservationPartition")


@dataclass(frozen=True, slots=True)
class RelativePermittivityRatioObservation:
    row_id: str
    source_id: str
    source_locator: str
    solvent_id: str
    component_ids: tuple[str, str, str]
    temperature_k: float
    pressure_pa: float
    total_ion_mole_fraction: float
    observed_relative_permittivity_ratio: float
    residual_scale: float
    partition: ObservationPartition

    def __post_init__(self) -> None:
        for field in ("row_id", "source_id", "source_locator", "solvent_id"):
            _require_nonempty_string(getattr(self, field), field)
        if type(self.component_ids) is not tuple or len(self.component_ids) != 3:
            raise ValueError("component_ids must contain solvent, cation, and anion")
        for component_id in self.component_ids:
            _require_nonempty_string(component_id, "component_id")
        if len(set(self.component_ids)) != 3:
            raise ValueError("component_ids must be distinct")
        if self.component_ids[0] != self.solvent_id:
            raise ValueError("solvent_id must match the first ordered component_id")
        for field in (
            "temperature_k",
            "pressure_pa",
            "observed_relative_permittivity_ratio",
            "residual_scale",
        ):
            _require_finite(getattr(self, field), field, positive=True)
        _require_finite(
            self.total_ion_mole_fraction,
            "total_ion_mole_fraction",
            positive=True,
        )
        if self.total_ion_mole_fraction >= 1.0:
            raise ValueError("total_ion_mole_fraction must be less than one")
        if not isinstance(self.partition, ObservationPartition):
            raise TypeError("partition must be an ObservationPartition")


DirectObservation = (
    MeanIonicActivityObservation
    | AqueousKijMeanIonicActivityObservation
    | IonSolvationKijObservation
    | SolvationGibbsObservation
    | RelativePermittivityRatioObservation
)
RegressionObservation = (
    FixedCompositionVleObservation
    | PureSaturationObservation
    | PureVaporPressureObservation
    | PureDensityObservation
    | FixedCompositionMixtureDensityObservation
    | DirectObservation
)


def _canonical_row(row: RegressionObservation) -> dict[str, object]:
    if isinstance(row, RelativePermittivityRatioObservation):
        return {
            "row_id": row.row_id,
            "source_id": row.source_id,
            "source_locator": row.source_locator,
            "solvent_id": row.solvent_id,
            "component_ids": list(row.component_ids),
            "temperature_k": row.temperature_k,
            "pressure_pa": row.pressure_pa,
            "total_ion_mole_fraction": row.total_ion_mole_fraction,
            "observed_relative_permittivity_ratio": (row.observed_relative_permittivity_ratio),
            "residual_scale": row.residual_scale,
            "partition": row.partition.value,
        }
    if isinstance(row, AqueousKijMeanIonicActivityObservation):
        return {
            "row_id": row.row_id,
            "source_id": row.source_id,
            "source_locator": row.source_locator,
            "component_ids": list(row.component_ids),
            "active_pair_component_ids": list(row.active_pair_component_ids),
            "fixed_k_ij": list(row.fixed_k_ij),
            "temperature_k": row.temperature_k,
            "pressure_pa": row.pressure_pa,
            "formula_unit_molality_mol_per_kg": (row.formula_unit_molality_mol_per_kg),
            "observed_mean_ionic_activity_coefficient": (
                row.observed_mean_ionic_activity_coefficient
            ),
            "relative_residual_scale": row.relative_residual_scale,
            "partition": row.partition.value,
        }
    if isinstance(row, IonSolvationKijObservation):
        return {
            "row_id": row.row_id,
            "source_id": row.source_id,
            "source_locator": row.source_locator,
            "component_ids": list(row.component_ids),
            "active_component_id": row.active_component_id,
            "active_pair_component_ids": list(row.active_pair_component_ids),
            "fixed_k_ij": list(row.fixed_k_ij),
            "temperature_k": row.temperature_k,
            "pressure_pa": row.pressure_pa,
            "observed_solvation_gibbs_j_per_mol": (row.observed_solvation_gibbs_j_per_mol),
            "residual_scale_j_per_mol": row.residual_scale_j_per_mol,
            "partition": row.partition.value,
        }
    if isinstance(row, MeanIonicActivityObservation):
        return {
            "row_id": row.row_id,
            "source_id": row.source_id,
            "source_locator": row.source_locator,
            "component_ids": list(row.component_ids),
            "active_component_id": row.active_component_id,
            "temperature_k": row.temperature_k,
            "pressure_pa": row.pressure_pa,
            "formula_unit_molality_mol_per_kg": (row.formula_unit_molality_mol_per_kg),
            "observed_mean_ionic_activity_coefficient": (
                row.observed_mean_ionic_activity_coefficient
            ),
            "relative_residual_scale": row.relative_residual_scale,
            "partition": row.partition.value,
        }
    if isinstance(row, SolvationGibbsObservation):
        return {
            "row_id": row.row_id,
            "source_id": row.source_id,
            "source_locator": row.source_locator,
            "component_ids": list(row.component_ids),
            "active_component_id": row.active_component_id,
            "temperature_k": row.temperature_k,
            "pressure_pa": row.pressure_pa,
            "observed_solvation_gibbs_j_per_mol": (row.observed_solvation_gibbs_j_per_mol),
            "residual_scale_j_per_mol": row.residual_scale_j_per_mol,
            "partition": row.partition.value,
        }
    if isinstance(row, PureDensityObservation):
        return {
            "row_id": row.row_id,
            "source_id": row.source_id,
            "source_locator": row.source_locator,
            "component_ids": list(row.component_ids),
            "temperature_k": row.temperature_k,
            "pressure_pa": row.pressure_pa,
            "density_kg_per_m3": row.density_kg_per_m3,
            "molar_mass_kg_per_mol": row.molar_mass_kg_per_mol,
            "pressure_scale_pa": row.pressure_scale_pa,
            "density_scale_kg_per_m3": row.density_scale_kg_per_m3,
            "volume_origin_m3_per_mol": row.volume_origin_m3_per_mol,
            "volume_start_m3_per_mol": row.volume_start_m3_per_mol,
            "volume_bounds_m3_per_mol": list(row.volume_bounds_m3_per_mol),
            "partition": row.partition.value,
        }
    if isinstance(row, FixedCompositionMixtureDensityObservation):
        return {
            "row_id": row.row_id,
            "source_id": row.source_id,
            "source_locator": row.source_locator,
            "component_ids": list(row.component_ids),
            "mole_fractions": list(row.mole_fractions),
            "temperature_k": row.temperature_k,
            "pressure_pa": row.pressure_pa,
            "density_kg_per_m3": row.density_kg_per_m3,
            "molar_mass_kg_per_mol": row.molar_mass_kg_per_mol,
            "pressure_scale_pa": row.pressure_scale_pa,
            "density_scale_kg_per_m3": row.density_scale_kg_per_m3,
            "volume_origin_m3_per_mol": row.volume_origin_m3_per_mol,
            "volume_start_m3_per_mol": row.volume_start_m3_per_mol,
            "volume_bounds_m3_per_mol": list(row.volume_bounds_m3_per_mol),
            "partition": row.partition.value,
        }
    if isinstance(row, PureSaturationObservation):
        return {
            "row_id": row.row_id,
            "source_id": row.source_id,
            "source_locator": row.source_locator,
            "component_ids": list(row.component_ids),
            "temperature_k": row.temperature_k,
            "pressure_pa": row.pressure_pa,
            "liquid_density_kg_per_m3": row.liquid_density_kg_per_m3,
            "molar_mass_kg_per_mol": row.molar_mass_kg_per_mol,
            "pressure_scale_pa": row.pressure_scale_pa,
            "chemical_potential_scale": row.chemical_potential_scale,
            "liquid_density_scale_kg_per_m3": (row.liquid_density_scale_kg_per_m3),
            "liquid_volume_origin_m3_per_mol": (row.liquid_volume_origin_m3_per_mol),
            "liquid_volume_start_m3_per_mol": (row.liquid_volume_start_m3_per_mol),
            "liquid_volume_bounds_m3_per_mol": list(row.liquid_volume_bounds_m3_per_mol),
            "vapor_volume_origin_m3_per_mol": (row.vapor_volume_origin_m3_per_mol),
            "vapor_volume_start_m3_per_mol": (row.vapor_volume_start_m3_per_mol),
            "vapor_volume_bounds_m3_per_mol": list(row.vapor_volume_bounds_m3_per_mol),
            "partition": row.partition.value,
        }
    if isinstance(row, PureVaporPressureObservation):
        return {
            "row_id": row.row_id,
            "source_id": row.source_id,
            "source_locator": row.source_locator,
            "component_ids": list(row.component_ids),
            "temperature_k": row.temperature_k,
            "pressure_pa": row.pressure_pa,
            "pressure_scale_pa": row.pressure_scale_pa,
            "chemical_potential_scale": row.chemical_potential_scale,
            "liquid_volume_origin_m3_per_mol": (row.liquid_volume_origin_m3_per_mol),
            "liquid_volume_start_m3_per_mol": (row.liquid_volume_start_m3_per_mol),
            "liquid_volume_bounds_m3_per_mol": list(row.liquid_volume_bounds_m3_per_mol),
            "vapor_volume_origin_m3_per_mol": (row.vapor_volume_origin_m3_per_mol),
            "vapor_volume_start_m3_per_mol": (row.vapor_volume_start_m3_per_mol),
            "vapor_volume_bounds_m3_per_mol": list(row.vapor_volume_bounds_m3_per_mol),
            "partition": row.partition.value,
        }
    return {
        "row_id": row.row_id,
        "source_id": row.source_id,
        "source_locator": row.source_locator,
        "component_ids": list(row.component_ids),
        "temperature_k": row.temperature_k,
        "pressure_pa": row.pressure_pa,
        "liquid_mole_fraction_first": row.liquid_mole_fraction_first,
        "vapor_mole_fraction_first": row.vapor_mole_fraction_first,
        "pressure_scale_pa": row.pressure_scale_pa,
        "chemical_potential_scales": list(row.chemical_potential_scales),
        "liquid_volume_origin_m3_per_mol": row.liquid_volume_origin_m3_per_mol,
        "liquid_volume_start_m3_per_mol": row.liquid_volume_start_m3_per_mol,
        "liquid_volume_bounds_m3_per_mol": list(row.liquid_volume_bounds_m3_per_mol),
        "vapor_volume_origin_m3_per_mol": row.vapor_volume_origin_m3_per_mol,
        "vapor_volume_start_m3_per_mol": row.vapor_volume_start_m3_per_mol,
        "vapor_volume_bounds_m3_per_mol": list(row.vapor_volume_bounds_m3_per_mol),
        "partition": row.partition.value,
    }


def canonical_dataset_sha256(
    observations: Iterable[RegressionObservation],
) -> str:
    rows = sorted(
        (_canonical_row(row) for row in observations), key=lambda row: str(row["row_id"])
    )
    payload = json.dumps(rows, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


@dataclass(frozen=True, slots=True)
class RegressionProblem:
    sources: tuple[SourceDescriptor, ...]
    parameters: tuple[ParameterCoordinate, ...]
    parameter_slot_indices: tuple[int, ...]
    start_vectors: tuple[tuple[float, ...], ...]
    observations: tuple[RegressionObservation, ...]
    maximum_condition_number: float
    maximum_iterations: int
    maximum_solver_time_seconds: float
    function_tolerance: float
    gradient_tolerance: float
    parameter_tolerance: float
    confirmation_parameter_scaled_max_delta: float
    confirmation_cost_relative_delta: float

    def __post_init__(self) -> None:
        if not self.sources:
            raise ValueError("at least one source is required")
        if not self.parameters:
            raise ValueError("at least one parameter is required")
        if (
            type(self.parameter_slot_indices) is not tuple
            or not self.parameter_slot_indices
            or any(type(index) is not int for index in self.parameter_slot_indices)
            or any(
                index < 0 or index >= len(self.parameters) for index in self.parameter_slot_indices
            )
        ):
            raise ValueError(
                "parameter_slot_indices must be a nonempty ordered map into parameters"
            )
        if set(self.parameter_slot_indices) != set(range(len(self.parameters))):
            raise ValueError("every fitted parameter must own at least one evaluator slot")
        if type(self.start_vectors) is not tuple or len(self.start_vectors) < 2:
            raise ValueError(
                "start_vectors must contain a primary and at least one confirmation vector"
            )
        for vector in self.start_vectors:
            if type(vector) is not tuple or len(vector) != len(self.parameters):
                raise ValueError("every start vector must match the ordered parameter count")
            for coordinate, start in zip(self.parameters, vector, strict=True):
                _require_finite(start, "parameter start")
                if not coordinate.lower_bound <= start <= coordinate.upper_bound:
                    raise ValueError("every parameter start must lie within its declared bounds")
        if not self.observations:
            raise ValueError("at least one observation is required")
        source_ids = tuple(source.source_id for source in self.sources)
        if len(set(source_ids)) != len(source_ids):
            raise ValueError("source_id values must be unique")
        row_ids = tuple(row.row_id for row in self.observations)
        if len(set(row_ids)) != len(row_ids):
            raise ValueError("duplicate row_id values are forbidden")
        parameter_keys = tuple(
            (parameter.family, parameter.identity) for parameter in self.parameters
        )
        if len(set(parameter_keys)) != len(parameter_keys):
            raise ValueError("duplicate parameter identity values are forbidden")
        source_map = {source.source_id: source for source in self.sources}
        for row in self.observations:
            if row.source_id not in source_map:
                raise ValueError(f"observation {row.row_id!r} references an unknown source_id")
        for source in self.sources:
            rows = tuple(row for row in self.observations if row.source_id == source.source_id)
            if not rows:
                raise ValueError(f"source {source.source_id!r} has no observations")
            if canonical_dataset_sha256(rows) != source.canonical_dataset_sha256:
                raise ValueError(
                    f"source {source.source_id!r} canonical dataset SHA-256 does not match its rows"
                )
        if not self.training_observations:
            raise ValueError("at least one training observation is required")
        mixture_density_observations = tuple(
            isinstance(row, FixedCompositionMixtureDensityObservation) for row in self.observations
        )
        if any(mixture_density_observations) and not all(mixture_density_observations):
            raise ValueError(
                "fixed-composition mixture-density observations must form a "
                "homogeneous regression problem"
            )
        density_observations = tuple(
            isinstance(
                row,
                (PureDensityObservation, FixedCompositionMixtureDensityObservation),
            )
            for row in self.observations
        )
        if (
            any(density_observations)
            and not all(density_observations)
            and not all(
                isinstance(
                    row,
                    (
                        PureDensityObservation,
                        FixedCompositionMixtureDensityObservation,
                        PureSaturationObservation,
                        PureVaporPressureObservation,
                    ),
                )
                for row in self.observations
            )
        ):
            raise ValueError(
                "density and phase-equilibrium observations cannot share one regression problem"
            )
        for parameter in self.parameters:
            identity = parameter.identity.canonical_component_ids
            for row in self.observations:
                row_identity = (
                    ()
                    if isinstance(parameter.identity, ModelParameterIdentity)
                    else ()
                    if isinstance(row, RelativePermittivityRatioObservation)
                    else row.canonical_active_pair_component_ids
                    if isinstance(
                        row,
                        (
                            AqueousKijMeanIonicActivityObservation,
                            IonSolvationKijObservation,
                        ),
                    )
                    else (row.component_ids[0],)
                    if (
                        parameter.family is ParameterFamily.RELATIVE_PERMITTIVITY
                        and isinstance(row, SolvationGibbsObservation)
                    )
                    else (row.active_component_id,)
                    if isinstance(
                        row,
                        (MeanIonicActivityObservation, SolvationGibbsObservation),
                    )
                    else (
                        parameter.identity.canonical_component_ids
                        if (
                            isinstance(
                                parameter.identity,
                                (
                                    PairParameterIdentity,
                                    ComponentParameterIdentity,
                                    AssociationParameterIdentity,
                                ),
                            )
                            and isinstance(
                                row,
                                FixedCompositionMixtureDensityObservation,
                            )
                            and set(parameter.identity.canonical_component_ids).issubset(
                                row.component_ids
                            )
                        )
                        else (
                            tuple(sorted(row.component_ids))
                            if len(row.component_ids) == 2
                            else row.component_ids
                        )
                    )
                )
                if row_identity != identity:
                    raise ValueError(
                        f"observation {row.row_id!r} identity does not match the parameter identity"
                    )
        _require_finite(
            self.maximum_condition_number,
            "maximum_condition_number",
            positive=True,
        )
        if type(self.maximum_iterations) is not int or self.maximum_iterations <= 0:
            raise ValueError("maximum_iterations must be a positive integer")
        for field in (
            "maximum_solver_time_seconds",
            "function_tolerance",
            "gradient_tolerance",
            "parameter_tolerance",
            "confirmation_parameter_scaled_max_delta",
            "confirmation_cost_relative_delta",
        ):
            _require_finite(getattr(self, field), field, positive=True)

    @property
    def training_observations(self) -> tuple[RegressionObservation, ...]:
        return tuple(
            row for row in self.observations if row.partition is ObservationPartition.TRAINING
        )

    @property
    def held_out_observations(self) -> tuple[RegressionObservation, ...]:
        return tuple(
            row for row in self.observations if row.partition is ObservationPartition.HELD_OUT
        )

    @property
    def stress_observations(self) -> tuple[RegressionObservation, ...]:
        return tuple(
            row for row in self.observations if row.partition is ObservationPartition.STRESS
        )

    @property
    def solver_start_vectors(self) -> tuple[tuple[float, ...], ...]:
        return tuple(
            tuple(
                coordinate.transform.to_solver(start)
                for coordinate, start in zip(self.parameters, vector, strict=True)
            )
            for vector in self.start_vectors
        )


def _row_payload(row: RegressionObservation) -> tuple[object, ...]:
    if isinstance(row, RelativePermittivityRatioObservation):
        return (
            row.row_id,
            row.partition.value,
            row.temperature_k,
            row.pressure_pa,
            row.total_ion_mole_fraction,
            row.observed_relative_permittivity_ratio,
            row.residual_scale,
        )
    if isinstance(row, AqueousKijMeanIonicActivityObservation):
        return (
            row.row_id,
            row.partition.value,
            row.temperature_k,
            row.pressure_pa,
            row.formula_unit_molality_mol_per_kg,
            row.observed_mean_ionic_activity_coefficient,
            row.relative_residual_scale,
            *row.fixed_k_ij,
        )
    if isinstance(row, IonSolvationKijObservation):
        return (
            row.row_id,
            row.partition.value,
            row.temperature_k,
            row.pressure_pa,
            float(row.component_ids.index(row.active_component_id)),
            row.observed_solvation_gibbs_j_per_mol,
            row.residual_scale_j_per_mol,
            *row.fixed_k_ij,
        )
    if isinstance(row, MeanIonicActivityObservation):
        return (
            row.row_id,
            row.partition.value,
            row.temperature_k,
            row.pressure_pa,
            row.formula_unit_molality_mol_per_kg,
            row.observed_mean_ionic_activity_coefficient,
            row.relative_residual_scale,
        )
    if isinstance(row, SolvationGibbsObservation):
        return (
            row.row_id,
            row.partition.value,
            row.temperature_k,
            row.pressure_pa,
            0.0,
            row.observed_solvation_gibbs_j_per_mol,
            row.residual_scale_j_per_mol,
        )
    if isinstance(row, PureDensityObservation):
        return (
            row.row_id,
            row.partition.value,
            row.temperature_k,
            row.pressure_pa,
            row.pressure_scale_pa,
            row.molar_mass_kg_per_mol,
            row.density_kg_per_m3,
            row.density_scale_kg_per_m3,
            row.volume_origin_m3_per_mol,
            row.volume_start_m3_per_mol,
            row.volume_bounds_m3_per_mol[0],
            row.volume_bounds_m3_per_mol[1],
        )
    if isinstance(row, FixedCompositionMixtureDensityObservation):
        return (
            row.row_id,
            row.partition.value,
            row.temperature_k,
            row.pressure_pa,
            row.mole_fractions,
            row.pressure_scale_pa,
            row.molar_mass_kg_per_mol,
            row.density_kg_per_m3,
            row.density_scale_kg_per_m3,
            row.volume_origin_m3_per_mol,
            row.volume_start_m3_per_mol,
            row.volume_bounds_m3_per_mol[0],
            row.volume_bounds_m3_per_mol[1],
        )
    if isinstance(row, PureSaturationObservation):
        return (
            row.row_id,
            row.partition.value,
            row.temperature_k,
            row.pressure_pa,
            row.pressure_scale_pa,
            row.chemical_potential_scale,
            row.molar_mass_kg_per_mol,
            row.liquid_density_kg_per_m3,
            row.liquid_density_scale_kg_per_m3,
            row.liquid_volume_origin_m3_per_mol,
            row.liquid_volume_start_m3_per_mol,
            row.liquid_volume_bounds_m3_per_mol[0],
            row.liquid_volume_bounds_m3_per_mol[1],
            row.vapor_volume_origin_m3_per_mol,
            row.vapor_volume_start_m3_per_mol,
            row.vapor_volume_bounds_m3_per_mol[0],
            row.vapor_volume_bounds_m3_per_mol[1],
        )
    if isinstance(row, PureVaporPressureObservation):
        return (
            row.row_id,
            row.partition.value,
            row.temperature_k,
            row.pressure_pa,
            row.pressure_scale_pa,
            row.chemical_potential_scale,
            row.liquid_volume_origin_m3_per_mol,
            row.liquid_volume_start_m3_per_mol,
            row.liquid_volume_bounds_m3_per_mol[0],
            row.liquid_volume_bounds_m3_per_mol[1],
            row.vapor_volume_origin_m3_per_mol,
            row.vapor_volume_start_m3_per_mol,
            row.vapor_volume_bounds_m3_per_mol[0],
            row.vapor_volume_bounds_m3_per_mol[1],
        )
    return (
        row.row_id,
        row.partition.value,
        row.temperature_k,
        row.pressure_pa,
        row.liquid_mole_fraction_first,
        row.vapor_mole_fraction_first,
        row.pressure_scale_pa,
        row.chemical_potential_scales[0],
        row.chemical_potential_scales[1],
        row.liquid_volume_origin_m3_per_mol,
        row.liquid_volume_start_m3_per_mol,
        row.liquid_volume_bounds_m3_per_mol[0],
        row.liquid_volume_bounds_m3_per_mol[1],
        row.vapor_volume_origin_m3_per_mol,
        row.vapor_volume_start_m3_per_mol,
        row.vapor_volume_bounds_m3_per_mol[0],
        row.vapor_volume_bounds_m3_per_mol[1],
    )


def _uses_fixed_topology_descriptor(
    parameters: Iterable[ParameterCoordinate],
) -> bool:
    return any(
        parameter.capability_id is None and parameter.provider_artifact_fingerprint is not None
        for parameter in parameters
    )


def _fixed_topology_slot_match(
    problem: RegressionProblem,
    capability: FixedTopologyAssociationCapability,
) -> tuple[int, ...]:
    if (
        not capability.installed_ready
        or capability.observation_contract != "fixed_composition_helmholtz_phase"
        or capability.model_domain
        not in ("neutral_associating_pure", "neutral_associating_mixture")
        or capability.tensor_layout != "row_major"
        or capability.state_coordinate_count != len(capability.component_ids) + 1
        or not capability.helmholtz_basis_id
        or not capability.association_basis_id
    ):
        raise ValueError("installed EOS fixed-topology association descriptor is incomplete")
    matched: list[tuple[int, int]] = []
    for parameter_index, parameter in enumerate(problem.parameters):
        if (
            parameter.capability_id is not None
            or parameter.provider_parameter_fingerprint != capability.parameter_fingerprint
            or parameter.provider_topology_fingerprint != capability.topology_fingerprint
            or parameter.provider_artifact_fingerprint != capability.artifact_fingerprint
        ):
            raise ValueError("association coordinate does not match the installed EOS artifact")
        if isinstance(parameter.identity, AssociationParameterIdentity):
            expected_identities: tuple[
                ComponentParameterIdentity
                | PairParameterIdentity
                | ModelParameterIdentity
                | AssociationParameterIdentity,
                ...,
            ] = tuple(
                AssociationParameterIdentity(
                    parameter.family,
                    (pair,),
                )
                for pair in parameter.identity.site_pairs
            )
        else:
            expected_identities = (parameter.identity,)
        parameter_matches = []
        for expected_identity in expected_identities:
            candidates = tuple(
                slot
                for slot in capability.slots
                if slot.family is parameter.family
                and slot.identity == expected_identity
                and slot.unit == parameter.unit
                and parameter.lower_bound > slot.lower_bound_exclusive
                and parameter.upper_bound < slot.upper_bound_exclusive
            )
            if len(candidates) != 1:
                raise ValueError(
                    "association coordinate identity, unit, or bounds do not "
                    "match exactly one advertised EOS slot"
                )
            parameter_matches.append(candidates[0].slot_index)
        matched.extend((slot_index, parameter_index) for slot_index in parameter_matches)
    matched.sort()
    slot_indices = tuple(slot for slot, _ in matched)
    if len(set(slot_indices)) != len(slot_indices):
        raise ValueError("fixed-topology association slots cannot be selected twice")
    expected_sharing = tuple(parameter for _, parameter in matched)
    if problem.parameter_slot_indices != expected_sharing:
        raise ValueError(
            "parameter_slot_indices must map the selected advertised EOS slots "
            "to fitted parameters in descriptor order"
        )
    return slot_indices


def _native_payload(
    problem: RegressionProblem,
    capability: ParameterCapability | FixedTopologyAssociationCapability,
) -> tuple[object, ...]:
    parameter = problem.parameters[0]
    observation_shape = (
        "fixed_topology_association"
        if _uses_fixed_topology_descriptor(problem.parameters)
        and all(
            isinstance(
                row,
                (
                    PureDensityObservation,
                    PureSaturationObservation,
                    PureVaporPressureObservation,
                ),
            )
            for row in problem.observations
        )
        else "pure_density"
        if isinstance(problem.observations[0], PureDensityObservation)
        and all(isinstance(row, PureDensityObservation) for row in problem.observations)
        else "mixture_density"
        if isinstance(problem.observations[0], FixedCompositionMixtureDensityObservation)
        and all(
            isinstance(row, FixedCompositionMixtureDensityObservation)
            for row in problem.observations
        )
        else "pure_vapor_pressure"
        if isinstance(problem.observations[0], PureVaporPressureObservation)
        and all(isinstance(row, PureVaporPressureObservation) for row in problem.observations)
        else "phase_or_direct"
    )
    payload = (
        parameter.capability_id or "fixed_topology_association",
        parameter.provider_parameter_fingerprint,
        parameter.provider_topology_fingerprint,
        capability.component_ids,
        parameter.transform.origin,
        parameter.transform.scale,
        parameter.lower_bound,
        parameter.upper_bound,
        tuple(vector[0] for vector in problem.start_vectors),
        problem.maximum_condition_number,
        problem.maximum_iterations,
        problem.maximum_solver_time_seconds,
        problem.function_tolerance,
        problem.gradient_tolerance,
        problem.parameter_tolerance,
        problem.confirmation_parameter_scaled_max_delta,
        problem.confirmation_cost_relative_delta,
        tuple(_row_payload(row) for row in problem.training_observations),
        tuple(
            _row_payload(row)
            for row in (*problem.held_out_observations, *problem.stress_observations)
        ),
        observation_shape,
    )
    if len(problem.parameters) == 1 and not _uses_fixed_topology_descriptor(problem.parameters):
        return payload
    joint_payload = (
        *payload,
        tuple(item.transform.origin for item in problem.parameters),
        tuple(item.transform.scale for item in problem.parameters),
        tuple(item.lower_bound for item in problem.parameters),
        tuple(item.upper_bound for item in problem.parameters),
        problem.start_vectors,
        problem.parameter_slot_indices,
    )
    if not _uses_fixed_topology_descriptor(problem.parameters):
        return joint_payload
    if not isinstance(capability, FixedTopologyAssociationCapability):
        raise TypeError("association regression requires its structural capability")
    return (
        *joint_payload,
        capability.artifact_fingerprint,
        _fixed_topology_slot_match(problem, capability),
    )


def _matched_capability(problem: RegressionProblem, model: Mixture) -> ParameterCapability:
    capabilities = parameter_capabilities(model)
    parameter = problem.parameters[0]
    matches = tuple(
        capability
        for capability in capabilities
        if isinstance(capability, ParameterCapability)
        and capability.capability_id == parameter.capability_id
    )
    if len(matches) != 1:
        raise ValueError("installed Provider does not advertise the requested parameter capability")
    capability = matches[0]
    if (
        capability.family is not parameter.family
        or capability.parameter_fingerprint != parameter.provider_parameter_fingerprint
        or capability.topology_fingerprint != parameter.provider_topology_fingerprint
    ):
        raise ValueError("regression parameter does not match the installed Provider capability")
    provider_unit = {
        "dimensionless": "1",
        "angstrom": "angstrom",
        "kelvin": "K",
    }.get(capability.coordinate_units[-1])
    if provider_unit != parameter.unit:
        raise ValueError(
            "regression parameter unit does not match the installed Provider capability"
        )
    for row in problem.observations:
        if row.component_ids != capability.component_ids:
            raise ValueError(
                f"observation {row.row_id!r} component order does not match the Provider model"
            )
        if not (capability.temperature_min_k <= row.temperature_k <= capability.temperature_max_k):
            raise ValueError(
                f"observation {row.row_id!r} temperature is outside the Provider capability domain"
            )
        if isinstance(
            row,
            (MeanIonicActivityObservation, AqueousKijMeanIonicActivityObservation),
        ):
            if (
                row.pressure_pa != 100_000.0
                or not 0.001 <= row.formula_unit_molality_mol_per_kg <= 6.0
            ):
                raise ValueError(
                    f"observation {row.row_id!r} pressure or molality is "
                    "outside the Provider direct-observable domain"
                )
            if capability.observation_contract != "aqueous_mean_ionic_activity" or (
                capability.active_component_ids != parameter.identity.canonical_component_ids
            ):
                raise ValueError(
                    "mean-ionic-activity observation does not match the "
                    "Provider direct-observable capability"
                )
        elif isinstance(row, (SolvationGibbsObservation, IonSolvationKijObservation)):
            if row.pressure_pa != 100_000.0:
                raise ValueError(
                    f"observation {row.row_id!r} pressure is outside the "
                    "Provider direct-observable domain"
                )
            if (
                capability.observation_contract != "ion_solvation_gibbs"
                or capability.active_component_ids != parameter.identity.canonical_component_ids
            ):
                raise ValueError(
                    "solvation-Gibbs observation does not match the Provider "
                    "direct-observable capability"
                )
            if capability.capability_id in (
                "ion_solvation_ionic_region_permittivity_v1",
                "ion_solvation_solvent_permittivity_v1",
            ) and (
                len(capability.component_ids) < 2
                or row.active_component_id != capability.component_ids[1]
            ):
                raise ValueError(
                    "solvation-Gibbs observation active ion does not match "
                    "the Provider callback's fixed observable ion"
                )
        elif isinstance(row, RelativePermittivityRatioObservation):
            if row.pressure_pa != 100_000.0:
                raise ValueError(
                    f"observation {row.row_id!r} pressure is outside the "
                    "Provider direct-observable domain"
                )
            if (
                capability.observation_contract != "relative_permittivity_ratio"
                or capability.identity_shape != "model"
                or capability.active_component_ids
            ):
                raise ValueError(
                    "relative-permittivity observation does not match the "
                    "Provider direct-observable capability"
                )
        elif (
            capability.observation_contract != "fixed_composition_helmholtz_phase"
            or (
                capability.identity_shape == "unordered_component_pair"
                and capability.active_component_ids != parameter.identity.canonical_component_ids
            )
            or (
                capability.identity_shape == "component"
                and capability.active_component_ids != parameter.identity.canonical_component_ids
            )
        ):
            raise ValueError("phase observation does not match the Provider capability")
    return capability


def _matched_capabilities(
    problem: RegressionProblem, model: Mixture
) -> tuple[ParameterCapability | FixedTopologyAssociationCapability, ...]:
    if _uses_fixed_topology_descriptor(problem.parameters):
        advertised = tuple(
            capability
            for capability in parameter_capabilities(model)
            if isinstance(capability, FixedTopologyAssociationCapability)
        )
        if len(advertised) != 1:
            raise ValueError(
                "installed EOS does not advertise exactly one fixed-topology association descriptor"
            )
        capability = advertised[0]
        _fixed_topology_slot_match(problem, capability)
        expected_rows = (
            (FixedCompositionMixtureDensityObservation,)
            if capability.model_domain == "neutral_associating_mixture"
            else (
                PureDensityObservation,
                PureSaturationObservation,
                PureVaporPressureObservation,
            )
        )
        for row in problem.observations:
            if (
                not isinstance(row, expected_rows)
                or row.component_ids != capability.component_ids
                or not capability.temperature_min_k
                <= row.temperature_k
                <= capability.temperature_max_k
            ):
                raise ValueError(
                    f"observation {row.row_id!r} does not match the installed "
                    "fixed-topology association descriptor"
                )
        return (capability,)
    return tuple(
        _matched_capability(
            replace(
                problem,
                parameters=(parameter,),
                parameter_slot_indices=(0,),
                start_vectors=tuple((vector[index],) for vector in problem.start_vectors),
            ),
            model,
        )
        for index, parameter in enumerate(problem.parameters)
    )


def _supported_capability_block(
    problem: RegressionProblem,
    capabilities: tuple[ParameterCapability | FixedTopologyAssociationCapability, ...],
) -> bool:
    if _uses_fixed_topology_descriptor(problem.parameters):
        return (
            len(capabilities) == 1
            and isinstance(capabilities[0], FixedTopologyAssociationCapability)
            and capabilities[0].model_domain
            in ("neutral_associating_pure", "neutral_associating_mixture")
            and all(
                isinstance(
                    row,
                    (
                        PureDensityObservation,
                        PureSaturationObservation,
                        PureVaporPressureObservation,
                        FixedCompositionMixtureDensityObservation,
                    ),
                )
                for row in problem.observations
            )
        )
    if len(capabilities) == 1:
        return True
    expected_families = (
        ParameterFamily.SEGMENT_COUNT,
        ParameterFamily.SEGMENT_DIAMETER,
        ParameterFamily.DISPERSION_ENERGY_OVER_K,
    )
    return (
        tuple(parameter.family for parameter in problem.parameters) == expected_families
        and problem.parameter_slot_indices == (0, 1, 2)
        and all(isinstance(row, PureSaturationObservation) for row in problem.observations)
        and len(capabilities) == 3
        and len({capability.component_ids for capability in capabilities}) == 1
        and len({capability.parameter_fingerprint for capability in capabilities}) == 1
        and len({capability.topology_fingerprint for capability in capabilities}) == 1
    )


def _evaluate_parameters(
    problem: RegressionProblem, model: Mixture, variables: tuple[float, ...]
) -> tuple[tuple[float, ...], tuple[float, ...]]:
    capabilities = _matched_capabilities(problem, model)
    if not _supported_capability_block(problem, capabilities):
        raise ValueError(
            "the installed capabilities do not expose one compatible multi-active evaluator"
        )
    capability = capabilities[0]
    residuals, jacobian = _native.evaluate_general(
        _native_transport(model), _native_payload(problem, capability), variables
    )
    return tuple(residuals), tuple(jacobian)


def _evaluate_parameters_at_start(
    problem: RegressionProblem,
    model: Mixture,
    physical_start: tuple[float, ...],
) -> tuple[tuple[float, ...], tuple[float, ...], tuple[float, ...]]:
    capabilities = _matched_capabilities(problem, model)
    if not _supported_capability_block(problem, capabilities):
        raise ValueError(
            "the installed capabilities do not expose one compatible multi-active evaluator"
        )
    variables, residuals, jacobian = _native.evaluate_general_start(
        _native_transport(model),
        _native_payload(problem, capabilities[0]),
        physical_start,
    )
    return tuple(variables), tuple(residuals), tuple(jacobian)


def _complete_failure_reasons(
    existing: tuple[str, ...],
    *,
    termination: str,
    solution_usable: bool,
    solver_converged: bool,
    confirmations_usable: bool,
    confirmation_parameter_delta: float,
    confirmation_parameter_limit: float,
    confirmation_cost_delta: float,
    confirmation_cost_limit: float,
    jacobian: GeneralJacobianDiagnostics,
    parameter_count: int,
    maximum_condition_number: float,
    workflow_valid: bool,
) -> tuple[str, ...]:
    reasons = list(existing)
    if not solution_usable:
        reasons.append("solution_unusable")
    if not solver_converged and termination != "CONVERGENCE":
        reasons.append(f"solver_not_converged: termination={termination}")
    elif not solver_converged and solution_usable:
        reasons.append("solver_rejected_despite_usable_convergence")
    if not confirmations_usable:
        reasons.append("confirmation_runs_unusable")
    if confirmation_parameter_delta > confirmation_parameter_limit:
        reasons.append("confirmation_parameter_disagreement")
    if confirmation_cost_delta > confirmation_cost_limit:
        reasons.append("confirmation_cost_disagreement")
    if jacobian.full_rank < jacobian.variable_count:
        reasons.append("rank_deficient_full_jacobian")
    if jacobian.projected_parameter_rank < parameter_count:
        reasons.append("rank_deficient_projected_parameter_jacobian")
    if (
        not math.isfinite(jacobian.full_condition_number)
        or jacobian.full_condition_number > maximum_condition_number
    ):
        reasons.append("poor_full_jacobian_conditioning")
    if (
        not math.isfinite(jacobian.projected_parameter_condition_number)
        or jacobian.projected_parameter_condition_number > maximum_condition_number
    ):
        reasons.append("poor_projected_parameter_conditioning")
    if not workflow_valid:
        reasons.append("incomplete_row_evaluation")
    return tuple(dict.fromkeys(reasons))


def _direct_observed_value(observation: RegressionObservation) -> float:
    if isinstance(
        observation,
        (MeanIonicActivityObservation, AqueousKijMeanIonicActivityObservation),
    ):
        return observation.observed_mean_ionic_activity_coefficient
    if isinstance(observation, RelativePermittivityRatioObservation):
        return observation.observed_relative_permittivity_ratio
    if isinstance(observation, (IonSolvationKijObservation, SolvationGibbsObservation)):
        return observation.observed_solvation_gibbs_j_per_mol
    raise TypeError("observation is not a direct scalar observation")


def _general_row_diagnostic(row: Any, observation: RegressionObservation) -> GeneralRowDiagnostic:
    if not isinstance(observation, FixedCompositionVleObservation):
        raise TypeError("general row requires a fixed-composition VLE observation")
    return GeneralRowDiagnostic(
        row_id=row[0],
        partition=row[1],
        liquid_volume_m3_per_mol=row[2],
        vapor_volume_m3_per_mol=row[3],
        scaled_residuals=tuple(row[4]),
        observed_pressure_pa=observation.pressure_pa,
        liquid_model_pressure_pa=(
            observation.pressure_pa + row[4][0] * observation.pressure_scale_pa
        ),
        vapor_model_pressure_pa=(
            observation.pressure_pa + row[4][1] * observation.pressure_scale_pa
        ),
        chemical_potential_differences_over_rt=(
            row[4][2] * observation.chemical_potential_scales[0],
            row[4][3] * observation.chemical_potential_scales[1],
        ),
        derivative_status=("EXACT_PROVIDER_HESSIAN" if row[5] else "UNAVAILABLE"),
        status="evaluated" if row[5] else "failed",
        evaluated=row[5],
        failure_reason=row[6],
    )


def fit_parameters(problem: RegressionProblem, model: Mixture) -> _package.RegressionResult:
    from .result import RegressionResult

    capabilities = _matched_capabilities(problem, model)
    if not _supported_capability_block(problem, capabilities):
        raise ValueError(
            "the installed capabilities do not expose one compatible multi-active evaluator"
        )
    capability = capabilities[0]
    native = _native.solve_general(_native_transport(model), _native_payload(problem, capability))
    vector_parameter_result = len(problem.parameters) > 1 or _uses_fixed_topology_descriptor(
        problem.parameters
    )
    physical_parameters = (
        tuple(float(value) for value in native[5])
        if vector_parameter_result
        else (float(native[5]),)
    )
    bound_distances = (
        tuple(float(value) for value in native[6])
        if vector_parameter_result
        else (float(native[6]),)
    )
    active_bounds = (
        tuple(str(value) for value in native[7]) if vector_parameter_result else (str(native[7]),)
    )
    observations = {row.row_id: row for row in problem.observations}
    rows: tuple[
        GeneralRowDiagnostic
        | PureSaturationRowDiagnostic
        | PureVaporPressureRowDiagnostic
        | DensityRowDiagnostic
        | DirectObservationRowDiagnostic,
        ...,
    ] = tuple(
        (
            DirectObservationRowDiagnostic(
                row_id=row[0],
                partition=row[1],
                observable=(
                    "mean_ionic_activity_coefficient"
                    if isinstance(
                        observation := observations[row[0]],
                        (
                            MeanIonicActivityObservation,
                            AqueousKijMeanIonicActivityObservation,
                        ),
                    )
                    else "relative_permittivity_ratio"
                    if isinstance(observation, RelativePermittivityRatioObservation)
                    else "solvation_gibbs_energy"
                ),
                observable_unit=(
                    "1"
                    if isinstance(
                        observation,
                        (
                            MeanIonicActivityObservation,
                            AqueousKijMeanIonicActivityObservation,
                            RelativePermittivityRatioObservation,
                        ),
                    )
                    else "J/mol"
                ),
                observed_value=_direct_observed_value(observation),
                modeled_value=row[7],
                scaled_residual=row[4][0],
                provider_derivative=row[8],
                derivative_status=("EXACT_PROVIDER_FIRST_DERIVATIVE" if row[5] else "UNAVAILABLE"),
                status="evaluated" if row[5] else "failed",
                evaluated=row[5],
                failure_reason=row[6],
            )
            if isinstance(
                observations[row[0]],
                (
                    MeanIonicActivityObservation,
                    AqueousKijMeanIonicActivityObservation,
                    IonSolvationKijObservation,
                    SolvationGibbsObservation,
                    RelativePermittivityRatioObservation,
                ),
            )
            else DensityRowDiagnostic(
                row_id=row[0],
                partition=row[1],
                volume_m3_per_mol=row[2],
                scaled_residuals=tuple(row[4]),
                observed_pressure_pa=observation.pressure_pa,
                model_pressure_pa=(
                    observation.pressure_pa + row[4][0] * observation.pressure_scale_pa
                ),
                observed_density_kg_per_m3=observation.density_kg_per_m3,
                model_density_kg_per_m3=(observation.molar_mass_kg_per_mol / row[2]),
                derivative_status=("EXACT_PROVIDER_HESSIAN" if row[5] else "UNAVAILABLE"),
                status="evaluated" if row[5] else "failed",
                evaluated=row[5],
                failure_reason=row[6],
            )
            if isinstance(
                observation := observations[row[0]],
                (PureDensityObservation, FixedCompositionMixtureDensityObservation),
            )
            else PureVaporPressureRowDiagnostic(
                row_id=row[0],
                partition=row[1],
                liquid_volume_m3_per_mol=row[2],
                vapor_volume_m3_per_mol=row[3],
                scaled_residuals=tuple(row[4]),
                observed_pressure_pa=observation.pressure_pa,
                liquid_model_pressure_pa=(
                    observation.pressure_pa + row[4][0] * observation.pressure_scale_pa
                ),
                vapor_model_pressure_pa=(
                    observation.pressure_pa + row[4][1] * observation.pressure_scale_pa
                ),
                chemical_potential_difference_over_rt=(
                    row[4][2] * observation.chemical_potential_scale
                ),
                derivative_status=("EXACT_PROVIDER_HESSIAN" if row[5] else "UNAVAILABLE"),
                status="evaluated" if row[5] else "failed",
                evaluated=row[5],
                failure_reason=row[6],
            )
            if isinstance(
                observation := observations[row[0]],
                PureVaporPressureObservation,
            )
            else PureSaturationRowDiagnostic(
                row_id=row[0],
                partition=row[1],
                liquid_volume_m3_per_mol=row[2],
                vapor_volume_m3_per_mol=row[3],
                scaled_residuals=tuple(row[4]),
                observed_pressure_pa=observation.pressure_pa,
                liquid_model_pressure_pa=(
                    observation.pressure_pa + row[4][0] * observation.pressure_scale_pa
                ),
                vapor_model_pressure_pa=(
                    observation.pressure_pa + row[4][1] * observation.pressure_scale_pa
                ),
                chemical_potential_difference_over_rt=(
                    row[4][2] * observation.chemical_potential_scale
                ),
                observed_liquid_density_kg_per_m3=(observation.liquid_density_kg_per_m3),
                model_liquid_density_kg_per_m3=(observation.molar_mass_kg_per_mol / row[2]),
                derivative_status=("EXACT_PROVIDER_HESSIAN" if row[5] else "UNAVAILABLE"),
                status="evaluated" if row[5] else "failed",
                evaluated=row[5],
                failure_reason=row[6],
            )
            if isinstance(
                observation := observations[row[0]],
                PureSaturationObservation,
            )
            else _general_row_diagnostic(row, observation)
        )
        for row in native[20]
    )
    parameter_diagnostics = tuple(
        FittedParameterDiagnostic(
            family=coordinate.family,
            component_ids=coordinate.identity.canonical_component_ids,
            unit=coordinate.unit,
            transform_origin=coordinate.transform.origin,
            transform_scale=coordinate.transform.scale,
            start=problem.start_vectors[0][index],
            final=physical_parameters[index],
            movement=(physical_parameters[index] - problem.start_vectors[0][index]),
            lower_bound=coordinate.lower_bound,
            upper_bound=coordinate.upper_bound,
            active_bound_distance=bound_distances[index],
            active_bound=active_bounds[index] or None,
            active_coordinate_ids=(
                coordinate.identity.canonical_ids
                if isinstance(coordinate.identity, AssociationParameterIdentity)
                else ()
            ),
        )
        for index, coordinate in enumerate(problem.parameters)
    )
    jacobian = GeneralJacobianDiagnostics(
        residual_count=native[23],
        variable_count=native[22],
        full_singular_values=tuple(native[10]),
        full_rank=native[11],
        full_condition_number=native[12],
        projected_parameter_singular_values=tuple(native[13]),
        projected_parameter_rank=native[14],
        projected_parameter_condition_number=native[15],
    )
    solver_converged = native[0] == "CONVERGENCE" and bool(native[1]) and not native[21]
    confirmations_usable = bool(native[19])
    numerically_converged = (
        solver_converged
        and confirmations_usable
        and native[17] <= problem.confirmation_parameter_scaled_max_delta
        and native[18] <= problem.confirmation_cost_relative_delta
        and jacobian.full_rank == jacobian.variable_count
        and jacobian.projected_parameter_rank == len(problem.parameters)
        and math.isfinite(jacobian.full_condition_number)
        and jacobian.full_condition_number <= problem.maximum_condition_number
        and math.isfinite(jacobian.projected_parameter_condition_number)
        and jacobian.projected_parameter_condition_number <= problem.maximum_condition_number
    )
    workflow_valid = all(row.evaluated and not row.failure_reason for row in rows) and len(
        rows
    ) == len(problem.observations)
    failures = tuple(
        reason
        for reason in (
            native[21],
            *(f"{row.row_id}: {row.failure_reason}" for row in rows if row.failure_reason),
        )
        if reason
    )
    failures = _complete_failure_reasons(
        failures,
        termination=str(native[0]),
        solution_usable=bool(native[1]),
        solver_converged=solver_converged,
        confirmations_usable=confirmations_usable,
        confirmation_parameter_delta=float(native[17]),
        confirmation_parameter_limit=problem.confirmation_parameter_scaled_max_delta,
        confirmation_cost_delta=float(native[18]),
        confirmation_cost_limit=problem.confirmation_cost_relative_delta,
        jacobian=jacobian,
        parameter_count=len(problem.parameters),
        maximum_condition_number=problem.maximum_condition_number,
        workflow_valid=workflow_valid,
    )
    return RegressionResult(
        problem=problem,
        capabilities=capabilities,
        preparation_fingerprint=None,
        provider_parameter_fingerprint=capability.parameter_fingerprint,
        provider_topology_fingerprint=capability.topology_fingerprint,
        solver_converged=solver_converged,
        numerically_converged=numerically_converged,
        workflow_valid=workflow_valid,
        physical_status="NOT_ADJUDICATED_NO_ROW_ACCEPTANCE_CRITERIA",
        scientific_status="NOT_ADJUDICATED_NO_APPROVED_SCIENTIFIC_CUTOFF",
        predictive_status="NOT_ADJUDICATED_NO_APPROVED_HELD_OUT_CUTOFF",
        authority_status="NO_AUTHORITY_CHANGE_REGRESSION_RESULT_ONLY",
        termination=native[0],
        solution_usable=bool(native[1]),
        initial_cost=native[2],
        final_cost=native[3],
        iterations=native[4],
        residual_evaluation_count=native[24],
        jacobian_evaluation_count=native[25],
        parameters=parameter_diagnostics,
        jacobian=jacobian,
        rows=rows,
        confirmation_count=native[16],
        confirmation_parameter_scaled_max_delta=native[17],
        confirmation_cost_relative_max_delta=native[18],
        confirmations_usable=confirmations_usable,
        training_row_count=len(problem.training_observations),
        held_out_row_count=len(problem.held_out_observations),
        stress_row_count=len(problem.stress_observations),
        evaluated_row_count=sum(row.evaluated for row in rows),
        skipped_row_count=0,
        failed_row_count=sum(not row.evaluated for row in rows),
        failure_reasons=failures,
    )
