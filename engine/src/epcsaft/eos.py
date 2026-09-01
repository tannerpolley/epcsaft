from __future__ import annotations

import hashlib
import json
import math
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from types import MappingProxyType
from typing import TYPE_CHECKING, Any, Literal, SupportsFloat, cast, overload

import pint

import epcsaft._core as _native_core

from ._resolve import _resolve_model, _ResolvedModel
from ._trials import NonEvaluableTrial
from .parameters import ActiveParameterSet, Parameters
from .records import (
    ConstantCorrelation,
    ConstantPlusPolynomialInTemperatureCorrelation,
    ConstantPlusSumOfExponentialsCorrelation,
    SingleParameterRecord,
)
from .units import ParameterError, _error, _quantity, unit_registry

if TYPE_CHECKING:
    from .reference import Activity, ElectrolyteReference


type Phase = Literal["vapor", "liquid"]
type PropertyBasis = Literal["molar", "reduced"]
type ObservableBasis = Literal["absolute", "mass", "molar", "reduced"]
type ObservableIdentity = Literal[
    "eos.fixed-pressure-fugacity",
    "eos.fixed-state",
    "eos.pure-saturation",
]


@dataclass(frozen=True, slots=True)
class ObservableTarget:
    """One source-facing target resolved to an EOS observable row."""

    positive: bool
    row_identity: str
    unit: str
    bases: tuple[str, ...]
    references: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class ObservableFamily:
    """One EOS-owned observation family admitted to Parameter Fit."""

    descriptor_identity: ObservableIdentity
    targets: tuple[str, ...]
    state_fields: tuple[str, ...]
    model_kinds: tuple[Literal["pure", "mixture"], ...]


OBSERVABLE_TARGETS: Mapping[str, ObservableTarget] = MappingProxyType(
    {
        "pressure": ObservableTarget(True, "pressure_pa", "pascal", ("absolute",), ("absolute",)),
        "saturation_pressure": ObservableTarget(
            True, "saturation_pressure_pa", "pascal", ("absolute",), ("absolute",)
        ),
        "saturated_liquid_mass_density": ObservableTarget(
            True,
            "liquid_mass_density_kg_per_m3",
            "kilogram / meter**3",
            ("mass",),
            ("saturated_liquid",),
        ),
        "saturated_liquid_molar_density": ObservableTarget(
            True,
            "liquid_molar_density_mol_per_m3",
            "mole / meter**3",
            ("molar",),
            ("saturated_liquid",),
        ),
        "saturated_vapor_mass_density": ObservableTarget(
            True,
            "vapor_mass_density_kg_per_m3",
            "kilogram / meter**3",
            ("mass",),
            ("saturated_vapor",),
        ),
        "saturated_vapor_molar_density": ObservableTarget(
            True,
            "vapor_molar_density_mol_per_m3",
            "mole / meter**3",
            ("molar",),
            ("saturated_vapor",),
        ),
    }
)

OBSERVABLE_FAMILIES: Mapping[str, ObservableFamily] = MappingProxyType(
    {
        "fixed_state_pressure": ObservableFamily(
            "eos.fixed-state",
            ("pressure",),
            ("temperature_k", "molar_density_mol_per_m3"),
            ("pure", "mixture"),
        ),
        "pure_saturation": ObservableFamily(
            "eos.pure-saturation",
            tuple(target for target in OBSERVABLE_TARGETS if target != "pressure"),
            ("temperature_k",),
            ("pure",),
        ),
    }
)
_MOLAR_GAS_CONSTANT = 8.31446261815324


@lru_cache(maxsize=1)
def _python_native_artifact_fingerprint() -> str:
    return _file_fingerprint(Path(_native_core.__file__))


def _file_fingerprint(path: Path) -> str:
    return f"sha256:{hashlib.sha256(path.read_bytes()).hexdigest()}"


def _magnitude(
    record: (
        SingleParameterRecord
        | ConstantCorrelation
        | ConstantPlusSumOfExponentialsCorrelation
        | ConstantPlusPolynomialInTemperatureCorrelation
    ),
    unit: str,
) -> float:
    value = record.value if isinstance(record, SingleParameterRecord) else record.constant
    if unit != "dimensionless":
        if not isinstance(value, pint.Quantity):
            raise TypeError(f"expected quantity convertible to {unit}")
        value = value.to(unit).magnitude
    return float(cast(SupportsFloat, value))


def _quantity_magnitude(value: object, unit: str) -> float:
    if not isinstance(value, pint.Quantity):
        raise TypeError(f"expected quantity convertible to {unit}")
    return float(value.to(unit).magnitude)


class EosError(ValueError):
    """One stable public error for invalid states and EOS solver failures."""

    __slots__ = ("_code", "_context")

    def __init__(self, code: str, context: Iterable[tuple[str, object]] = ()) -> None:
        self._code = code
        self._context = tuple((str(key), str(value)) for key, value in context)
        details = ", ".join(f"{key}={value}" for key, value in self._context)
        super().__init__(f"{code} ({details})" if details else code)

    @property
    def code(self) -> str:
        return self._code

    @property
    def context(self) -> tuple[tuple[str, str], ...]:
        return self._context


def _eos_error(code: str, **context: object) -> EosError:
    return EosError(code, context.items())


@dataclass(frozen=True, slots=True)
class SaturationPhase:
    """One pure phase on a locally solved saturation boundary."""

    amount_mol: float
    volume_m3: float
    molar_density_mol_m3: float
    pressure_pa: float
    chemical_potential_over_rt: float


@dataclass(frozen=True, slots=True)
class SaturationStart:
    """Optional caller-provided vapor/liquid densities and pressure."""

    vapor_density: pint.Quantity[Any]
    liquid_density: pint.Quantity[Any]
    pressure: pint.Quantity[Any]


@dataclass(frozen=True, slots=True)
class SaturationDiagnostics:
    candidate_available: bool
    optimizer_converged: bool
    optimizer_status: str
    vapor_density_bounds_mol_m3: tuple[float, float] | None
    liquid_density_bounds_mol_m3: tuple[float, float] | None
    feasibility_measure: float | None
    feasibility_tolerance: float | None
    physical_accepted: bool
    physical_status: str
    physical_status_count: int | None
    pressure_relative_residual: float | None
    chemical_potential_absolute_residual: float | None
    phase_density_distance: float | None
    vapor_pressure_relative_residual: float | None
    liquid_pressure_relative_residual: float | None
    vapor_residual_density_derivative_m3_per_mol: float | None
    liquid_residual_density_derivative_m3_per_mol: float | None
    rejection_status: Literal["none", "solve-failure", "root-swap", "coalescence", "branch-loss"]
    exact_derivatives: bool
    sensitivity_status: str
    sensitivity_failure_reason: str
    sensitivity_rank: int | None
    sensitivity_dimension: int | None
    sensitivity_condition_number: float | None
    sensitivity_backward_error: float | None
    globality_certificate: bool
    failure_reason: str


@dataclass(frozen=True, slots=True)
class SaturationResult:
    temperature_k: float
    saturation_pressure_pa: float
    parameter_fingerprint: str
    vapor: SaturationPhase
    liquid: SaturationPhase
    diagnostics: SaturationDiagnostics
    certified: bool


@dataclass(frozen=True, slots=True)
class ObservableDescriptor:
    identity: ObservableIdentity
    version: int
    required_inputs: tuple[str, ...]
    applicability: str
    model_identity: str
    supported_topologies: tuple[str, ...]
    model_domain: str
    temperature_domain_k: tuple[float | None, float | None]
    row_identities: tuple[str, ...]
    row_units: tuple[str, ...]
    row_bases: tuple[ObservableBasis, ...]
    state_fields: tuple[str, ...]
    reference: str
    active_parameter_identities: tuple[str, ...]
    active_parameter_values: tuple[float, ...]
    active_parameter_units: tuple[str, ...]
    derivative_orders: tuple[int, ...]
    source_derivative_order: int
    jacobian_shape: tuple[int, int]
    jacobian_axes: tuple[str, str]
    jacobian_layout: Literal["row-major"]
    parameter_fingerprint: str
    active_parameter_fingerprint: str
    topology_fingerprint: str


@dataclass(frozen=True, slots=True)
class ObservableState:
    temperature_k: float
    molar_density_mol_m3: float | None
    pressure_pa: float | None
    composition: tuple[float, ...]
    phase: Phase | None
    reference: str


@dataclass(frozen=True, slots=True)
class ObservableRowStatus:
    identity: str
    status: Literal["accepted", "rejected"]
    reason: str


@dataclass(frozen=True, slots=True)
class ObservableBlockStatus:
    scope: Literal["not-required", "row", "shared-solve", "branch"]
    status: Literal["accepted", "rejected"]
    reason: str


@dataclass(frozen=True, slots=True)
class SaturationBranchEvidence:
    identity: Phase
    density_bracket_mol_per_m3: tuple[float, float]
    pressure_relative_residual: float | None
    residual_density_derivative_m3_per_mol: float | None
    rejection_status: Literal["none", "solve-failure", "root-swap", "coalescence", "branch-loss"]


@dataclass(frozen=True, slots=True)
class ObservableBlock:
    descriptor: ObservableDescriptor
    state: ObservableState
    values: tuple[float | None, ...]
    jacobian: tuple[tuple[float, ...] | None, ...]
    row_statuses: tuple[ObservableRowStatus, ...]
    status: ObservableBlockStatus
    saturation_diagnostics: SaturationDiagnostics | None = None
    saturation_branches: tuple[SaturationBranchEvidence, ...] = ()


def compile_regression_evaluator(
    specification: Mapping[str, Any],
    active_parameters: ActiveParameterSet,
) -> object:
    """Compile one EOS observation specification into a native owner.

    ``specification`` is deliberately a small mapping so source packets can
    carry it without importing the runtime. It requires canonical ``family``,
    state, and ``row_identities`` fields.
    """

    if not isinstance(specification, Mapping):
        raise TypeError("EOS regression specification must be a mapping")
    model = specification.get("model")
    if not isinstance(model, _EosModel):
        raise TypeError("EOS regression specification requires an epcsaft EOS model")
    model._validate_active_parameters(active_parameters)
    family = specification.get("family")
    if family not in {"fixed-state", "pure-saturation"}:
        raise ValueError("unsupported EOS regression family")
    rows = tuple(cast(Iterable[str], specification.get("row_identities", ())))
    if not rows:
        raise ValueError("EOS regression specification requires row_identities")
    temperature = specification.get("temperature_k")
    if temperature is None:
        raise ValueError("EOS regression specification requires temperature_k")
    temperature_k = float(temperature)
    if not math.isfinite(temperature_k):
        raise ValueError("temperature must be finite")
    from .regression._api import CompiledEvaluator

    native_spec: dict[str, object] = {
        "family": family,
        "temperature_k": temperature_k,
        "row_identities": rows,
        "active_parameter_identities": active_parameters.identities,
        "active_parameter_units": tuple(
            coordinate.unit for coordinate in active_parameters.coordinates
        ),
        "active_parameter_bases": tuple("physical" for _ in active_parameters.coordinates),
    }
    if family == "fixed-state":
        density = specification.get("molar_density_mol_m3")
        if density is None:
            raise ValueError("fixed-state EOS regression requires molar density")
        native_spec["molar_density_mol_m3"] = float(density)
        composition = specification.get("mole_fractions")
        if composition is None:
            raise ValueError("fixed-state EOS regression requires mole_fractions")
        native_spec["mole_fractions"] = tuple(
            float(value) for value in cast(Iterable[Any], composition)
        )
    source_payload = {
        "family": family,
        "temperature_k": temperature_k,
        "molar_density_mol_m3": native_spec.get("molar_density_mol_m3"),
        "mole_fractions": native_spec.get("mole_fractions"),
        "rows": rows,
    }
    native_spec["source_spec_fingerprint"] = specification.get(
        "source_spec_fingerprint",
        _native_identity_fingerprint("epcsaft.eos-regression-source-v1", source_payload),
    )
    native_spec["evaluator_fingerprint"] = specification.get(
        "evaluator_fingerprint",
        _native_identity_fingerprint(
            "epcsaft.eos-regression-evaluator-v1",
            {"source": native_spec["source_spec_fingerprint"], "rows": rows},
        ),
    )
    native = _native_core._compile_eos_regression_evaluator(
        model._native_model(active_parameters),
        native_spec,
        model._active_coordinates(active_parameters),
    )
    return CompiledEvaluator.from_native(native)


def _accepted_observable_block(
    descriptor: ObservableDescriptor,
    state: ObservableState,
    values: tuple[float, ...],
    jacobian: tuple[tuple[float, ...], ...],
    *,
    scope: Literal["not-required", "shared-solve"] = "not-required",
    saturation_diagnostics: SaturationDiagnostics | None = None,
    saturation_branches: tuple[SaturationBranchEvidence, ...] = (),
) -> ObservableBlock:
    if len(values) != len(descriptor.row_identities) or len(jacobian) != len(values):
        raise RuntimeError("observable block rows do not match its descriptor")
    row_statuses: list[ObservableRowStatus] = []
    retained_values: list[float | None] = []
    retained_jacobian: list[tuple[float, ...] | None] = []
    for identity, value, row in zip(descriptor.row_identities, values, jacobian, strict=True):
        finite = (
            math.isfinite(value)
            and len(row) == descriptor.jacobian_shape[1]
            and all(math.isfinite(derivative) for derivative in row)
        )
        row_statuses.append(
            ObservableRowStatus(
                identity, "accepted" if finite else "rejected", "" if finite else "nonfinite-row"
            )
        )
        retained_values.append(value if finite else None)
        retained_jacobian.append(row if finite else None)
    rejected = any(status.status == "rejected" for status in row_statuses)
    return ObservableBlock(
        descriptor=descriptor,
        state=state,
        values=tuple(retained_values),
        jacobian=tuple(retained_jacobian),
        row_statuses=tuple(row_statuses),
        status=ObservableBlockStatus(
            "row" if rejected else scope,
            "rejected" if rejected else "accepted",
            "one or more observable rows are nonfinite" if rejected else "",
        ),
        saturation_diagnostics=saturation_diagnostics,
        saturation_branches=saturation_branches,
    )


def _saturation_branch_evidence(
    diagnostics: SaturationDiagnostics,
) -> tuple[SaturationBranchEvidence, ...]:
    if (
        diagnostics.vapor_density_bounds_mol_m3 is None
        or diagnostics.liquid_density_bounds_mol_m3 is None
    ):
        return ()
    vapor_residual = diagnostics.vapor_pressure_relative_residual
    liquid_residual = diagnostics.liquid_pressure_relative_residual
    vapor_derivative = diagnostics.vapor_residual_density_derivative_m3_per_mol
    liquid_derivative = diagnostics.liquid_residual_density_derivative_m3_per_mol
    return (
        SaturationBranchEvidence(
            "vapor",
            diagnostics.vapor_density_bounds_mol_m3,
            vapor_residual
            if vapor_residual is not None and math.isfinite(vapor_residual)
            else None,
            vapor_derivative
            if vapor_derivative is not None and math.isfinite(vapor_derivative)
            else None,
            diagnostics.rejection_status,
        ),
        SaturationBranchEvidence(
            "liquid",
            diagnostics.liquid_density_bounds_mol_m3,
            liquid_residual
            if liquid_residual is not None and math.isfinite(liquid_residual)
            else None,
            liquid_derivative
            if liquid_derivative is not None and math.isfinite(liquid_derivative)
            else None,
            diagnostics.rejection_status,
        ),
    )


def _rejected_saturation_block(
    descriptor: ObservableDescriptor,
    temperature_k: float,
    diagnostics: SaturationDiagnostics,
) -> ObservableBlock:
    reason = diagnostics.failure_reason or diagnostics.rejection_status
    return ObservableBlock(
        descriptor=descriptor,
        state=ObservableState(
            temperature_k=temperature_k,
            molar_density_mol_m3=None,
            pressure_pa=None,
            composition=(1.0,),
            phase=None,
            reference=descriptor.reference,
        ),
        values=tuple(None for _ in descriptor.row_identities),
        jacobian=tuple(None for _ in descriptor.row_identities),
        row_statuses=tuple(
            ObservableRowStatus(identity, "rejected", reason)
            for identity in descriptor.row_identities
        ),
        status=ObservableBlockStatus(
            "shared-solve" if diagnostics.rejection_status == "solve-failure" else "branch",
            "rejected",
            reason,
        ),
        saturation_diagnostics=diagnostics,
        saturation_branches=_saturation_branch_evidence(diagnostics),
    )


class SaturationError(RuntimeError):
    """Raised when no finite in-domain saturation candidate is available."""

    def __init__(self, message: str, diagnostics: SaturationDiagnostics) -> None:
        super().__init__(message)
        self.diagnostics = diagnostics


def _native_arguments(
    model: _ResolvedModel,
) -> tuple[object, ...]:
    sigma_constant: list[float] = []
    term_component: list[int] = []
    term_amplitude: list[float] = []
    term_exponent: list[float] = []
    for component, record in enumerate(model.segment_diameter):
        if isinstance(record, SingleParameterRecord):
            sigma_constant.append(_magnitude(record, "angstrom"))
        else:
            sigma_constant.append(_quantity_magnitude(record.constant, "angstrom"))
            if isinstance(record, ConstantPlusSumOfExponentialsCorrelation):
                for term in record.terms:
                    term_component.append(component)
                    term_amplitude.append(_quantity_magnitude(term.amplitude, "angstrom"))
                    term_exponent.append(
                        _quantity_magnitude(term.exponent_coefficient, "1 / kelvin")
                    )
            elif not isinstance(record, ConstantCorrelation):
                raise _error(
                    "invalid-value",
                    family="segment_diameter",
                    reason="unsupported-correlation",
                )

    component_position = {
        component_id: index for index, component_id in enumerate(model.component_ids)
    }
    site_position = {site_id: index for index, site_id in enumerate(model.site_ids)}
    sites = len(model.site_ids)
    association_active = [0] * (sites * sites)
    association_energy = [0.0] * (sites * sites)
    association_volume = [0.0] * (sites * sites)
    for pair in model.association_pairs:
        a = site_position[pair.endpoint_a]
        b = site_position[pair.endpoint_b]
        for index in (a * sites + b, b * sites + a):
            association_active[index] = 1
            association_energy[index] = pair.energy_over_k
            association_volume[index] = pair.volume

    components = len(model.component_ids)
    permittivity_term_component: list[int] = []
    permittivity_term_coefficient: list[float] = []
    permittivity_term_power: list[int] = []
    permittivity_reference_temperature = [0.0] * components
    relative_permittivity_formulation_codes = {
        "none": 0,
        "ion-fraction-suppression": 1,
        "component-permittivity-mixing": 2,
        "density-volume-temperature": 3,
        "ion-specific-suppression": 4,
        "dipolar-density-perturbation": 5,
        "solvent-only": 6,
    }
    if model.topology in {
        "ionic-fully-dissociated-debye-huckel",
        "ionic-born",
    }:
        molar_mass = [0.0 if value is None else value for value in model.molar_mass_kg_per_mol]
        charge_number = list(model.charge_number)
        relative_permittivity = []
        for component, permittivity_record in enumerate(model.relative_permittivity):
            relative_permittivity.append(
                0.0
                if permittivity_record is None
                else _magnitude(permittivity_record, "dimensionless")
            )
            if isinstance(
                permittivity_record,
                ConstantPlusPolynomialInTemperatureCorrelation,
            ):
                permittivity_reference_temperature[component] = _quantity_magnitude(
                    permittivity_record.reference_temperature, "kelvin"
                )
                for polynomial_term in permittivity_record.terms:
                    permittivity_term_component.append(component)
                    permittivity_term_coefficient.append(
                        _quantity_magnitude(
                            polynomial_term.coefficient,
                            f"1 / kelvin ** {polynomial_term.power}",
                        )
                    )
                    permittivity_term_power.append(polynomial_term.power)
        packing_diameter = [
            0.0 if value is None else value for value in model.packing_diameter_angstrom
        ]
        debye_huckel_diameter = [
            0.0 if value is None else value for value in model.debye_huckel_diameter_angstrom
        ]
        if model.topology == "ionic-born":
            dielectric_suppression = (
                0.0
                if model.ion_fraction_suppression_coefficient is None
                else model.ion_fraction_suppression_coefficient
            )
            born_diameter = [
                0.0 if value is None else value for value in model.born_diameter_angstrom
            ]
            solvation_factor = [
                _magnitude(record, "dimensionless")
                for record in model.solvation_factor
                if record is not None
            ]
            ionic_permittivity = model.ionic_region_relative_permittivity
            electrolyte_model = 1
            born_c_shell = model.born_c_shell
            born_c_dielectric = model.born_c_dielectric
            assert ionic_permittivity is not None
            assert len(solvation_factor) == components
        elif model.topology == "ionic-fully-dissociated-debye-huckel":
            dielectric_suppression = 0.0
            born_diameter = [0.0] * components
            solvation_factor = [0.0] * components
            ionic_permittivity = 0.0
            electrolyte_model = 2
            born_c_shell = 0.0
            born_c_dielectric = 0.0
    else:
        molar_mass = [0.0] * components
        charge_number = [0] * components
        relative_permittivity = [0.0] * components
        dielectric_suppression = 0.0
        packing_diameter = [0.0] * components
        debye_huckel_diameter = [0.0] * components
        born_diameter = [0.0] * components
        solvation_factor = [0.0] * components
        ionic_permittivity = 0.0
        born_c_shell = 0.0
        born_c_dielectric = 0.0
        electrolyte_model = 0

    return (
        list(model.segment_count),
        sigma_constant,
        term_component,
        term_amplitude,
        term_exponent,
        list(model.dispersion_energy_over_k),
        list(model.k_ij),
        list(model.k_ij_reciprocal_temperature_slope),
        list(model.k_ij_reference_temperature_k),
        list(model.l_ij),
        [component_position[component_id] for component_id, _ in model.site_ids],
        list(model.site_multiplicity),
        association_active,
        association_energy,
        association_volume,
        molar_mass,
        charge_number,
        relative_permittivity,
        permittivity_term_component,
        permittivity_term_coefficient,
        permittivity_term_power,
        permittivity_reference_temperature,
        dielectric_suppression,
        packing_diameter,
        debye_huckel_diameter,
        born_diameter,
        solvation_factor,
        ionic_permittivity,
        born_c_shell,
        born_c_dielectric,
        electrolyte_model,
        list(model.dielectric_volume_m3_per_mol) or [0.0] * components,
        list(model.dielectric_temperature_k) or [0.0] * components,
        list(model.ion_specific_suppression_coefficient) or [0.0] * components,
        list(model.dipole_scaling_debye2) or [0.0] * components,
        list(model.polarizability_scaling_angstrom3) or [0.0] * components,
        list(model.correlation_integral_parameter) or [0.0] * components,
        relative_permittivity_formulation_codes[model.relative_permittivity_formulation],
    )


@dataclass(frozen=True, slots=True, init=False)
class _EosModel:
    _component_ids: tuple[str, ...]
    _parameter_fingerprint: str
    _topology_fingerprint: str
    _component_order_fingerprint: str
    _artifact_fingerprint: str
    _parameters: Parameters
    _model: _ResolvedModel
    _native: _native_core.NativeModel

    def __init__(self, parameters: Parameters) -> None:
        model = _resolve_model(parameters)
        native_arguments = _native_arguments(model)
        object.__setattr__(self, "_component_ids", parameters.component_ids)
        object.__setattr__(self, "_parameter_fingerprint", parameters.fingerprint)
        object.__setattr__(self, "_parameters", parameters)
        object.__setattr__(self, "_model", model)
        topology_fingerprint = _native_topology_fingerprint(self)
        component_order_fingerprint = _native_identity_fingerprint(
            "epcsaft-native-component-order-v1",
            {
                "component_ids": parameters.component_ids,
                "charges": tuple(model.charge_number),
            },
        )
        artifact_fingerprint = _python_native_artifact_fingerprint()
        object.__setattr__(self, "_topology_fingerprint", topology_fingerprint)
        object.__setattr__(self, "_component_order_fingerprint", component_order_fingerprint)
        object.__setattr__(self, "_artifact_fingerprint", artifact_fingerprint)
        native = _native_core.NativeModel(*native_arguments).with_parameter_identities(
            parameters.component_ids,
            tuple(f"{component}/{site}" for component, site in model.site_ids),
        )
        native = native.with_eos_metadata(
            parameters.fingerprint,
            math.nan if model.temperature_min_k is None else model.temperature_min_k,
            math.nan if model.temperature_max_k is None else model.temperature_max_k,
            math.nan if model.ion_mole_fraction_max is None else model.ion_mole_fraction_max,
            math.nan if model.pressure_min_pa is None else model.pressure_min_pa,
            math.nan if model.pressure_max_pa is None else model.pressure_max_pa,
            topology_fingerprint,
            component_order_fingerprint,
            _native_identity_fingerprint(
                "epcsaft-native-source-domain-v1",
                {
                    "qualification": model.candidate_qualification,
                    "candidate_source_fingerprint": model.candidate_source_fingerprint,
                    "temperature_min_k": model.temperature_min_k,
                    "temperature_max_k": model.temperature_max_k,
                    "pressure_min_pa": model.pressure_min_pa,
                    "pressure_max_pa": model.pressure_max_pa,
                    "ion_mole_fraction_max": model.ion_mole_fraction_max,
                    "source_temperature_min_k": model.source_temperature_min_k,
                    "source_temperature_max_k": model.source_temperature_max_k,
                    "source_pressure_min_pa": model.source_pressure_min_pa,
                    "source_pressure_max_pa": model.source_pressure_max_pa,
                },
            ),
            artifact_fingerprint,
            math.nan if model.source_temperature_min_k is None else model.source_temperature_min_k,
            math.nan if model.source_temperature_max_k is None else model.source_temperature_max_k,
            math.nan if model.source_pressure_min_pa is None else model.source_pressure_min_pa,
            math.nan if model.source_pressure_max_pa is None else model.source_pressure_max_pa,
            model.candidate_qualification or "",
            model.candidate_source_fingerprint or "",
        )
        object.__setattr__(self, "_native", native)

    @classmethod
    def from_mapping(cls, mapping: Mapping[str, object]) -> _EosModel:
        """Construct a model from one strict canonical Parameters document."""
        return cls(Parameters.from_mapping(mapping))

    @property
    def component_ids(self) -> tuple[str, ...]:
        return self._component_ids

    @property
    def parameter_fingerprint(self) -> str:
        return self._parameter_fingerprint

    def _validate_active_parameters(self, active_parameters: ActiveParameterSet) -> None:
        if not isinstance(active_parameters, ActiveParameterSet):
            raise _eos_error("invalid-active-parameter-set", expected="ActiveParameterSet")
        if active_parameters.parameters.fingerprint != self._parameter_fingerprint:
            raise _eos_error("invalid-active-parameter-set", reason="different-parameters")

    def _active_coordinates(
        self, active_parameters: ActiveParameterSet
    ) -> tuple[_native_core.NativeActiveCoordinate, ...]:
        self._validate_active_parameters(active_parameters)
        result = []
        for coordinate in active_parameters.coordinates:
            applicability = coordinate.applicability
            if applicability is None:
                raise _eos_error(
                    "unsupported-parameter-derivative",
                    observation="fixed-state-helmholtz",
                    identity=coordinate.identity,
                    reason="physical-application-unresolved",
                )
            if (
                self._model.ionic_dispersion == "exclude-same-sign-ion-pairs"
                and coordinate.family == "k_ij"
                and len(applicability.component_indices) == 2
                and self._model.charge_number
                and self._model.charge_number[applicability.component_indices[0]]
                * self._model.charge_number[applicability.component_indices[1]]
                > 0
            ):
                raise _eos_error(
                    "unsupported-parameter-derivative",
                    observation="fixed-state-helmholtz",
                    identity=coordinate.identity,
                    reason="inactive-by-selected-model",
                )
            result.append(
                _native_core.NativeActiveCoordinate(
                    coordinate.identity,
                    coordinate.family,
                    applicability.application,
                    applicability.component_indices,
                    applicability.association_edges,
                    applicability.association_energy_scales,
                    coordinate.value,
                    applicability.reference_temperature_k,
                )
            )
        return tuple(result)

    def _native_model(
        self,
        active_parameters: ActiveParameterSet | None = None,
    ) -> _native_core.NativeModel:
        if active_parameters is None:
            return self._native
        return self._native.with_active_coordinates(self._active_coordinates(active_parameters))

    def _descriptor(
        self,
        active_parameters: ActiveParameterSet,
        identity: ObservableIdentity,
        rows: tuple[str, ...],
        units: tuple[str, ...],
        bases: tuple[ObservableBasis, ...],
    ) -> ObservableDescriptor:
        self._active_coordinates(active_parameters)
        required_inputs: tuple[str, ...]
        if identity == "eos.fixed-state":
            required_inputs = ("temperature_k", "molar_density_mol_m3")
            if len(self._component_ids) > 1:
                required_inputs = (*required_inputs, "composition")
            applicability = "single-phase fixed-density state"
            reference = "residual Helmholtz energy at the declared reference amount"
        elif identity == "eos.fixed-pressure-fugacity":
            required_inputs = ("temperature_k", "pressure_pa", "phase")
            if len(self._component_ids) > 1:
                required_inputs = (*required_inputs, "composition")
            applicability = "stable single-phase fixed-pressure state"
            reference = "ideal-gas mixture at the same temperature, pressure, and composition"
        else:
            required_inputs = ("temperature_k",)
            applicability = "pure-component vapor-liquid coexistence boundary"
            reference = "shared pure-component vapor-liquid saturation state"
        return ObservableDescriptor(
            identity=identity,
            version=1,
            required_inputs=required_inputs,
            applicability=applicability,
            model_identity=self._model.topology,
            supported_topologies=(self._model.topology,),
            model_domain=self._model.topology,
            temperature_domain_k=(
                self._model.temperature_min_k,
                self._model.temperature_max_k,
            ),
            row_identities=rows,
            row_units=units,
            row_bases=bases,
            state_fields=required_inputs,
            reference=reference,
            active_parameter_identities=active_parameters.identities,
            active_parameter_values=active_parameters.values,
            active_parameter_units=tuple(
                coordinate.unit for coordinate in active_parameters.coordinates
            ),
            derivative_orders=(0, 1),
            source_derivative_order=2,
            jacobian_shape=(len(rows), len(active_parameters.coordinates)),
            jacobian_axes=("observable_row", "active_parameter"),
            jacobian_layout="row-major",
            parameter_fingerprint=self._parameter_fingerprint,
            active_parameter_fingerprint=active_parameters.fingerprint,
            topology_fingerprint=_native_topology_fingerprint(self),
        )

    def _state(
        self,
        *,
        T: object,
        rho: object | None = None,
        P: object | None = None,
        x: object,
        phase: str | None = None,
        density_anchor: object | None = None,
        native_model: _native_core.NativeModel | None = None,
    ) -> State:
        """Create one immutable state at explicit density or specified pressure.

        Dimensioned inputs must be Pint quantities. Supply exactly one of
        ``rho`` and ``P``; ``phase`` is only valid for pressure
        closure. Invalid states and failed density/EOS evaluations raise ``EosError``.
        """
        if (rho is None) == (P is None):
            raise _eos_error("invalid-state", reason="specify-exactly-one-of-density-or-pressure")
        if phase not in (None, "vapor", "liquid"):
            raise _eos_error("invalid-state", field="phase")
        if rho is not None and phase is not None:
            raise _eos_error("invalid-state", field="phase")
        if density_anchor is not None and (rho is not None or P is None or phase is None):
            raise _eos_error(
                "invalid-density-anchor",
                reason="requires-pressure-closure-and-explicit-phase",
            )
        temperature_k = _state_magnitude(T, "kelvin", "T")
        try:
            fractions = tuple(
                float(cast(SupportsFloat, value)) for value in cast(Iterable[object], x)
            )
        except (TypeError, ValueError) as exc:
            raise _eos_error("invalid-state", field="x", reason="not-numeric-sequence") from exc
        if len(fractions) != len(self._component_ids):
            raise _eos_error(
                "invalid-state",
                field="x",
                reason="component-count",
                value=len(fractions),
                expected=len(self._component_ids),
            )
        if any(not math.isfinite(value) for value in fractions):
            raise _eos_error("invalid-state", field="x", reason="nonfinite")
        if any(value < 0.0 for value in fractions):
            raise _eos_error("invalid-state", field="x", reason="negative")
        fraction_sum = math.fsum(fractions)
        if abs(fraction_sum - 1.0) > 1.0e-12:
            raise _eos_error(
                "invalid-state",
                field="x",
                reason="not-normalized",
                value=fraction_sum,
                expected=1.0,
                tolerance=1.0e-12,
            )
        if (
            self._model.temperature_min_k is not None
            and temperature_k < self._model.temperature_min_k
        ) or (
            self._model.temperature_max_k is not None
            and temperature_k > self._model.temperature_max_k
        ):
            context: dict[str, object] = {
                "field": "T",
                "reason": "outside-domain",
                "value": temperature_k,
            }
            if self._model.temperature_min_k is not None:
                context["minimum"] = self._model.temperature_min_k
            if self._model.temperature_max_k is not None:
                context["maximum"] = self._model.temperature_max_k
            context["unit"] = "kelvin"
            raise _eos_error("invalid-state", **context)

        if self._model.topology == "ionic-born":
            ion_fraction = math.fsum(
                fraction
                for fraction, charge in zip(fractions, self._model.charge_number, strict=True)
                if charge != 0
            )
            if (
                self._model.ion_mole_fraction_max is not None
                and ion_fraction > self._model.ion_mole_fraction_max + 1.0e-12
            ):
                raise _eos_error(
                    "invalid-state",
                    field="x",
                    reason="ion-fraction-limit",
                    value=ion_fraction,
                    maximum=self._model.ion_mole_fraction_max,
                )
            charge_balance = math.fsum(
                fraction * charge
                for fraction, charge in zip(fractions, self._model.charge_number, strict=True)
            )
            if abs(charge_balance) > 1.0e-12:
                raise _eos_error(
                    "invalid-state",
                    field="electroneutrality",
                    reason="charge-balance",
                    value=charge_balance,
                    tolerance=1.0e-12,
                )

        evaluator = self._native if native_model is None else native_model
        anchor: State | None = None
        if rho is not None:
            density_molar = _state_magnitude(rho, "mole / meter ** 3", "rho")
            try:
                native = evaluator.evaluate(temperature_k, density_molar, fractions)
            except (RuntimeError, ValueError) as exc:
                raise _eos_error("evaluation-failed", reason=str(exc)) from exc
        else:
            pressure_pa = _state_magnitude(P, "pascal", "P")
            if density_anchor is None:
                try:
                    native = evaluator.evaluate_pressure(
                        temperature_k, pressure_pa, fractions, phase or ""
                    )
                except (RuntimeError, ValueError) as exc:
                    raise _eos_error("evaluation-failed", reason=str(exc)) from exc
            else:
                anchor = _validate_density_anchor(
                    self,
                    density_anchor,
                    temperature_k,
                    pressure_pa,
                    cast(str, phase),
                )
                certified = cast(Mapping[str, object], anchor.certified_branch)
                try:
                    native = evaluator.evaluate_pressure_continuation(
                        temperature_k,
                        pressure_pa,
                        fractions,
                        float(anchor.molar_density.to("mole / meter ** 3").magnitude),
                        cast(str, phase),
                        int(str(certified["stable_root_index"])),
                        float.fromhex(str(certified["origin_molar_density_hex"])),
                    )
                except (RuntimeError, ValueError) as exc:
                    raise _eos_error(
                        "density-anchor-branch-unreachable",
                        reason=str(exc),
                        max_abs_log_density=math.log1p(1.0e-2),
                        centered_scan_intervals=5,
                        domain_scan_intervals=360,
                    ) from exc
        diagnostics = native.density_diagnostics
        certified_branch = (
            anchor.certified_branch if anchor is not None else _certified_branch(diagnostics, phase)
        )
        state_identity = _state_identity(
            artifact_fingerprint=self._artifact_fingerprint,
            parameter_fingerprint=self._parameter_fingerprint,
            topology_fingerprint=self._topology_fingerprint,
            component_order_fingerprint=self._component_order_fingerprint,
            component_ids=self._component_ids,
            temperature_k=temperature_k,
            pressure_pa=pressure_pa if P is not None else native.pressure,
            mole_fractions=fractions,
            molar_density=native.molar_density,
            requested_phase=phase,
            certified_branch=certified_branch,
        )
        return State(
            component_ids=self._component_ids,
            temperature=temperature_k * unit_registry.kelvin,
            mole_fractions=fractions,
            parameter_fingerprint=self._parameter_fingerprint,
            topology_fingerprint=self._topology_fingerprint,
            component_order_fingerprint=self._component_order_fingerprint,
            artifact_fingerprint=self._artifact_fingerprint,
            requested_phase=cast(Phase | None, phase),
            certified_branch=certified_branch,
            state_identity=state_identity,
            hard_chain=native.hard_chain,
            dispersion=native.dispersion,
            association=native.association,
            debye_huckel=native.debye_huckel,
            born=native.born,
            bulk_relative_permittivity=None
            if self._model.relative_permittivity_formulation == "none"
            else native.bulk_relative_permittivity,
            _residual_helmholtz=native.residual_helmholtz,
            compressibility_factor=native.compressibility_factor,
            pressure=native.pressure * unit_registry.pascal,
            residual_chemical_potential_over_rt=tuple(native.residual_chemical_potential_over_rt),
            _log_fugacity_coefficient=None
            if native.log_fugacity_coefficient is None
            else tuple(native.log_fugacity_coefficient),
            _residual_entropy=native.residual_entropy
            * unit_registry.joule
            / unit_registry.mole
            / unit_registry.kelvin,
            _residual_internal_energy=native.residual_internal_energy
            * unit_registry.joule
            / unit_registry.mole,
            _residual_enthalpy=native.residual_enthalpy * unit_registry.joule / unit_registry.mole,
            _residual_isochoric_heat_capacity=native.residual_isochoric_heat_capacity
            * unit_registry.joule
            / unit_registry.mole
            / unit_registry.kelvin,
            _residual_isobaric_heat_capacity=None
            if native.residual_isobaric_heat_capacity is None
            else native.residual_isobaric_heat_capacity
            * unit_registry.joule
            / unit_registry.mole
            / unit_registry.kelvin,
            molar_density=native.molar_density * unit_registry.mole / unit_registry.meter**3,
            density_diagnostics=None
            if diagnostics is None
            else DensityDiagnostics(
                bracket_lower=diagnostics.bracket_lower
                * unit_registry.mole
                / unit_registry.meter**3,
                bracket_upper=diagnostics.bracket_upper
                * unit_registry.mole
                / unit_registry.meter**3,
                pressure_residual=diagnostics.pressure_residual * unit_registry.pascal,
                requested_pressure=pressure_pa * unit_registry.pascal,
                pressure_density_derivative=diagnostics.pressure_density_derivative
                * unit_registry.pascal
                * unit_registry.meter**3
                / unit_registry.mole,
                iterations=diagnostics.iterations,
                stable=diagnostics.stable,
                branch=diagnostics.branch,
                stable_root_count=(None if anchor is not None else diagnostics.stable_root_count),
                selected_stable_root_index=(
                    None if anchor is not None else diagnostics.selected_stable_root_index
                ),
                anchor_status="accepted" if anchor is not None else "not_supplied",
                anchor_density=None if anchor is None else anchor.molar_density,
                anchor_source_composition=None if anchor is None else anchor.mole_fractions,
                anchor_source_state_identity=None if anchor is None else anchor.state_identity,
                requested_phase=cast(Phase | None, phase),
                certified_branch=certified_branch,
                certificate_kind=diagnostics.certificate_kind,
                local_search_outcome=diagnostics.local_search_outcome,
                local_root_count=(diagnostics.local_root_count if anchor is not None else None),
                reachability_envelope=(
                    MappingProxyType(
                        {
                            "schema": "epcsaft-density-anchor-envelope-v1",
                            "max_abs_log_density": diagnostics.reachability_log_density_limit,
                            "centered_scan_intervals": diagnostics.reachability_scan_interval_count,
                            "domain_scan_intervals": 360,
                        }
                    )
                    if anchor is not None
                    else None
                ),
            ),
            fixed_pressure_composition_derivatives=(
                _fixed_pressure_composition_derivatives(
                    self, evaluator, temperature_k, fractions, native.molar_density
                )
                if P is not None and phase is not None and len(fractions) > 1
                else None
            ),
        )

    def _fixed_state_observable_block(
        self,
        active_parameters: ActiveParameterSet,
        *,
        T: object,
        rho: object,
        x: object,
    ) -> ObservableBlock:
        """Evaluate one fixed-density EOS-owned observable block."""

        coordinates = self._active_coordinates(active_parameters)
        state = self._state(T=T, rho=rho, x=x)
        temperature_k = float(state.temperature.to("kelvin").magnitude)
        density = float(state.molar_density.to("mole / meter ** 3").magnitude)
        amounts = state.mole_fractions
        volume = 1.0 / density
        parameter_count = len(coordinates)
        component_count = len(amounts)
        width = component_count + 1 + parameter_count
        parameter_coordinates = tuple(range(component_count + 1, width))
        try:
            block = self._native.evaluate_phase_derivatives(
                temperature_k,
                amounts,
                volume,
                coordinates,
                2,
            )
        except (KeyError, RuntimeError, ValueError) as exc:
            if isinstance(exc, EosError):
                raise
            raise _eos_error(
                "unsupported-parameter-derivative",
                observation="fixed-state-helmholtz",
                reason=str(exc),
            ) from exc

        gradient = tuple(float(value) for value in block.gradient)
        hessian = tuple(float(value) for value in block.hessian)
        helmholtz_derivatives = tuple(gradient[coordinate] for coordinate in parameter_coordinates)
        pressure_derivatives = tuple(
            -_MOLAR_GAS_CONSTANT * temperature_k * hessian[component_count * width + coordinate]
            for coordinate in parameter_coordinates
        )
        chemical_derivatives = tuple(
            tuple(hessian[component * width + coordinate] for coordinate in parameter_coordinates)
            for component in range(component_count)
        )
        labels = (
            "helmholtz_over_rt_reference_amount",
            "pressure_pa",
            *(f"chemical_potential_over_rt/{component}" for component in self._component_ids),
        )
        descriptor = self._descriptor(
            active_parameters,
            "eos.fixed-state",
            labels,
            ("1", "pascal", *("1" for _ in self._component_ids)),
            cast(
                tuple[ObservableBasis, ...],
                ("reduced", "absolute", *("reduced" for _ in self._component_ids)),
            ),
        )
        return _accepted_observable_block(
            descriptor,
            ObservableState(
                temperature_k=temperature_k,
                molar_density_mol_m3=density,
                pressure_pa=float(block.pressure_pa),
                composition=state.mole_fractions,
                phase=None,
                reference=descriptor.reference,
            ),
            (
                float(block.helmholtz_over_rt_reference_amount),
                float(block.pressure_pa),
                *(gradient[component] for component in range(component_count)),
            ),
            (
                helmholtz_derivatives,
                pressure_derivatives,
                *chemical_derivatives,
            ),
        )

    def _fixed_pressure_fugacity_observable_block(
        self,
        active_parameters: ActiveParameterSet,
        *,
        T: object,
        P: object,
        x: object,
        phase: Phase,
    ) -> ObservableBlock:
        """Evaluate one fixed-pressure fugacity observable block."""

        coordinates = self._active_coordinates(active_parameters)
        active_native = self._native.with_active_coordinates(coordinates)
        state = self._state(
            T=T,
            P=P,
            x=x,
            phase=phase,
            native_model=active_native,
        )
        temperature_k = float(state.temperature.to("kelvin").magnitude)
        density = float(state.molar_density.to("mole / meter ** 3").magnitude)
        volume = 1.0 / density
        try:
            block = self._native.evaluate_phase_derivatives(
                temperature_k,
                state.mole_fractions,
                volume,
                coordinates,
                2,
            )
        except (RuntimeError, ValueError) as exc:
            raise _eos_error(
                "unsupported-parameter-derivative",
                observation="fixed-pressure-fugacity",
                reason=str(exc),
            ) from exc
        parameter_count = len(coordinates)
        component_count = len(self._component_ids)
        width = component_count + 1 + parameter_count
        volume_coordinate = component_count
        hessian = tuple(float(value) for value in block.hessian)
        pressure_volume = hessian[volume_coordinate * width + volume_coordinate]
        if not math.isfinite(pressure_volume) or pressure_volume == 0.0:
            raise _eos_error(
                "non-evaluable-parameter-derivative",
                observation="fixed-pressure-fugacity",
                reason="singular-pressure-closure",
            )
        if state._log_fugacity_coefficient is None:
            raise _eos_error(
                "non-evaluable-parameter-derivative",
                observation="fixed-pressure-fugacity",
                reason="fugacity-unavailable",
            )
        values = state._log_fugacity_coefficient
        # EqID: fixed_pressure_total_derivative
        derivatives = tuple(
            tuple(
                hessian[component * width + component_count + 1 + parameter]
                - hessian[component * width + volume_coordinate]
                * hessian[volume_coordinate * width + component_count + 1 + parameter]
                / pressure_volume
                for parameter in range(parameter_count)
            )
            for component in range(component_count)
        )
        descriptor = self._descriptor(
            active_parameters,
            "eos.fixed-pressure-fugacity",
            tuple(f"log_fugacity_coefficient/{component}" for component in self._component_ids),
            tuple("1" for _ in self._component_ids),
            tuple("reduced" for _ in self._component_ids),
        )
        return _accepted_observable_block(
            descriptor,
            ObservableState(
                temperature_k=temperature_k,
                molar_density_mol_m3=density,
                pressure_pa=float(state.pressure.to("pascal").magnitude),
                composition=state.mole_fractions,
                phase=phase,
                reference=descriptor.reference,
            ),
            values,
            derivatives,
            scope="shared-solve",
        )

    def observable_descriptor(
        self,
        active_parameters: ActiveParameterSet,
        *,
        identity: ObservableIdentity,
    ) -> ObservableDescriptor:
        if identity == "eos.fixed-state":
            labels = (
                "helmholtz_over_rt_reference_amount",
                "pressure_pa",
                *(f"chemical_potential_over_rt/{component}" for component in self._component_ids),
            )
            return self._descriptor(
                active_parameters,
                identity,
                labels,
                ("1", "pascal", *("1" for _ in self._component_ids)),
                cast(
                    tuple[ObservableBasis, ...],
                    ("reduced", "absolute", *("reduced" for _ in self._component_ids)),
                ),
            )
        if identity == "eos.fixed-pressure-fugacity":
            return self._descriptor(
                active_parameters,
                identity,
                tuple(f"log_fugacity_coefficient/{component}" for component in self._component_ids),
                tuple("1" for _ in self._component_ids),
                tuple("reduced" for _ in self._component_ids),
            )
        if identity != "eos.pure-saturation":
            raise _eos_error("unsupported-observable-descriptor", identity=identity)
        raise _eos_error("unsupported-observable-descriptor", identity=identity)

    def _reference(
        self,
        *,
        T: object,
        P: object,
        solvent_x: object,
        phase: Literal["liquid"],
        ion_component_ids: tuple[str, str],
        formula_stoichiometry: tuple[int, int],
        reference_convention: Literal["molality-infinite-dilution"],
    ) -> ElectrolyteReference:
        """Build one liquid electrolyte reference on the selected component order.

        The model must contain one source-complete cation/anion formula unit and
        solvent-only normalized reference fractions. Failures raise ``EosError``.
        """
        from .reference import build_reference_state

        return build_reference_state(
            cast(Mixture, self),
            temperature=T,
            pressure=P,
            solvent_mole_fractions=solvent_x,
            phase=phase,
            ion_component_ids=ion_component_ids,
            formula_stoichiometry=formula_stoichiometry,
            reference_convention=reference_convention,
        )


@dataclass(frozen=True, slots=True, init=False)
class Pure(_EosModel):
    """A strict one-component EOS model."""

    def __init__(self, parameters: Parameters) -> None:
        if len(parameters.component_ids) != 1:
            raise _eos_error(
                "invalid-model",
                reason="pure-requires-one-component",
                components=len(parameters.component_ids),
            )
        _EosModel.__init__(self, parameters)

    @overload
    def state(
        self,
        *,
        T: object,
        rho: pint.Quantity[Any],
        P: None = None,
        phase: None = None,
    ) -> State: ...

    @overload
    def state(
        self,
        *,
        T: pint.Quantity[Any],
        rho: None = None,
        P: pint.Quantity[Any],
        phase: Phase | None = None,
    ) -> State: ...

    def state(
        self,
        *,
        T: object,
        rho: object | None = None,
        P: object | None = None,
        phase: str | None = None,
    ) -> State:
        return _EosModel._state(self, T=T, rho=rho, P=P, x=(1.0,), phase=phase)

    def observable_block(
        self,
        active_parameters: ActiveParameterSet,
        *,
        identity: ObservableIdentity,
        T: object,
        rho: object | None = None,
        P: object | None = None,
        phase: Phase | None = None,
    ) -> ObservableBlock:
        if identity == "eos.fixed-state" and rho is not None and P is None and phase is None:
            return _EosModel._fixed_state_observable_block(
                self, active_parameters, T=T, rho=rho, x=(1.0,)
            )
        if identity == "eos.fixed-pressure-fugacity" and rho is None and P is not None and phase:
            return _EosModel._fixed_pressure_fugacity_observable_block(
                self, active_parameters, T=T, P=P, x=(1.0,), phase=phase
            )
        if identity == "eos.pure-saturation" and rho is None and P is None and phase is None:
            return self._pure_saturation_observable_block(active_parameters, T=T)
        raise _eos_error("invalid-observable-inputs", identity=identity)

    def _saturation_temperature(self, T: object) -> float:
        temperature_k = _state_magnitude(T, "kelvin", "T")
        minimum = self._model.temperature_min_k
        maximum = self._model.temperature_max_k
        if (minimum is not None and temperature_k < minimum) or (
            maximum is not None and temperature_k > maximum
        ):
            raise _eos_error(
                "invalid-state",
                field="T",
                reason="outside-domain",
                value=temperature_k,
                minimum=minimum,
                maximum=maximum,
                unit="kelvin",
            )
        return temperature_k

    def _saturation_result(
        self, native: Mapping[str, object], temperature_k: float
    ) -> SaturationResult:
        diagnostics = _saturation_diagnostics(cast(Mapping[str, object], native["diagnostics"]))
        if not diagnostics.candidate_available:
            reason = diagnostics.failure_reason or "local saturation boundary was rejected"
            raise SaturationError(reason, diagnostics)
        return SaturationResult(
            temperature_k=temperature_k,
            saturation_pressure_pa=float(cast(float, native["saturation_pressure_pa"])),
            parameter_fingerprint=self._parameter_fingerprint,
            vapor=_saturation_phase(cast(Mapping[str, object], native["vapor"])),
            liquid=_saturation_phase(cast(Mapping[str, object], native["liquid"])),
            diagnostics=diagnostics,
            certified=bool(native["accepted"]),
        )

    def saturation(
        self, *, T: pint.Quantity[Any], start: SaturationStart | None = None
    ) -> SaturationResult:
        """Solve a local pure-component vapor-liquid saturation boundary."""

        temperature_k = self._saturation_temperature(T)
        try:
            liquid_density_upper = self._native.molar_density_upper_bound(temperature_k, (1.0,))
            explicit_start = (
                None
                if start is None
                else (
                    _state_magnitude(start.vapor_density, "mole / meter ** 3", "vapor_density"),
                    _state_magnitude(start.liquid_density, "mole / meter ** 3", "liquid_density"),
                    _state_magnitude(start.pressure, "pascal", "pressure"),
                )
            )
            native = cast(
                Mapping[str, object],
                _native_core._solve_saturation(
                    self._native,
                    temperature_k,
                    liquid_density_upper,
                    explicit_start,
                ),
            )
        except (RuntimeError, ValueError) as exc:
            failure = _saturation_diagnostics(_failed_saturation_diagnostics(str(exc)))
            raise SaturationError(str(exc), failure) from exc
        return self._saturation_result(native, temperature_k)

    def _pure_saturation_observable_block(
        self,
        active_parameters: ActiveParameterSet,
        *,
        T: object,
    ) -> ObservableBlock:
        """Solve one shared saturation state and return its EOS-owned block."""

        descriptor = self.observable_descriptor(active_parameters, identity="eos.pure-saturation")
        temperature_k = self._saturation_temperature(T)
        coordinates = self._active_coordinates(active_parameters)
        try:
            native = cast(
                Mapping[str, object],
                _native_core._solve_saturation_active(
                    self._native,
                    temperature_k,
                    coordinates,
                ),
            )
        except (RuntimeError, ValueError) as exc:
            failure = _saturation_diagnostics(_failed_saturation_diagnostics(str(exc)))
            raise SaturationError(str(exc), failure) from exc
        diagnostics = _saturation_diagnostics(cast(Mapping[str, object], native["diagnostics"]))
        if not bool(native["accepted"]) or not diagnostics.exact_derivatives:
            return _rejected_saturation_block(descriptor, temperature_k, diagnostics)
        result = self._saturation_result(native, temperature_k)
        parameter_count = len(coordinates)
        flat = tuple(float(value) for value in cast(Iterable[SupportsFloat], native["derivatives"]))
        derivatives = tuple(
            flat[row * parameter_count : (row + 1) * parameter_count] for row in range(3)
        )
        values: tuple[float, ...] = (
            result.saturation_pressure_pa,
            result.vapor.molar_density_mol_m3,
            result.liquid.molar_density_mol_m3,
        )
        molar_mass = self._model.molar_mass_kg_per_mol[0]
        if molar_mass is not None:
            mass_values = (
                result.vapor.molar_density_mol_m3 * molar_mass,
                result.liquid.molar_density_mol_m3 * molar_mass,
            )
            mass_derivatives = (
                tuple(value * molar_mass for value in derivatives[1]),
                tuple(value * molar_mass for value in derivatives[2]),
            )
            values = (
                *values,
                *mass_values,
            )
            derivatives = (
                *derivatives,
                *mass_derivatives,
            )
        return _accepted_observable_block(
            descriptor,
            ObservableState(
                temperature_k=temperature_k,
                molar_density_mol_m3=None,
                pressure_pa=result.saturation_pressure_pa,
                composition=(1.0,),
                phase=None,
                reference=descriptor.reference,
            ),
            values,
            derivatives,
            scope="shared-solve",
            saturation_diagnostics=diagnostics,
            saturation_branches=_saturation_branch_evidence(diagnostics),
        )

    def observable_descriptor(
        self,
        active_parameters: ActiveParameterSet,
        *,
        identity: ObservableIdentity,
    ) -> ObservableDescriptor:
        """Return one versioned EOS-owned observable descriptor."""

        if identity != "eos.pure-saturation":
            return _EosModel.observable_descriptor(self, active_parameters, identity=identity)
        labels = [
            "saturation_pressure_pa",
            "vapor_molar_density_mol_per_m3",
            "liquid_molar_density_mol_per_m3",
        ]
        units = ["pascal", "mole / meter ** 3", "mole / meter ** 3"]
        if self._model.molar_mass_kg_per_mol[0] is not None:
            labels.extend(
                (
                    "vapor_mass_density_kg_per_m3",
                    "liquid_mass_density_kg_per_m3",
                )
            )
            units.extend(("kilogram / meter ** 3", "kilogram / meter ** 3"))
        return self._descriptor(
            active_parameters,
            identity,
            tuple(labels),
            tuple(units),
            (
                "absolute",
                "molar",
                "molar",
                *(("mass", "mass") if self._model.molar_mass_kg_per_mol[0] is not None else ()),
            ),
        )


@dataclass(frozen=True, slots=True, init=False)
class Mixture(_EosModel):
    """A strict multi-component EOS model with explicit composition."""

    def __init__(self, parameters: Parameters) -> None:
        if len(parameters.component_ids) < 2:
            raise _eos_error(
                "invalid-model",
                reason="mixture-requires-at-least-two-components",
                components=len(parameters.component_ids),
            )
        _EosModel.__init__(self, parameters)

    @overload
    def state(
        self,
        *,
        T: pint.Quantity[Any],
        rho: pint.Quantity[Any],
        P: None = None,
        x: Iterable[SupportsFloat],
        phase: None = None,
        density_anchor: None = None,
    ) -> State: ...

    @overload
    def state(
        self,
        *,
        T: pint.Quantity[Any],
        rho: None = None,
        P: pint.Quantity[Any],
        x: Iterable[SupportsFloat],
        phase: Phase,
        density_anchor: State | None = None,
    ) -> State: ...

    def state(
        self,
        *,
        T: object,
        rho: object | None = None,
        P: object | None = None,
        x: object,
        phase: str | None = None,
        density_anchor: object | None = None,
    ) -> State:
        return _EosModel._state(
            self,
            T=T,
            rho=rho,
            P=P,
            x=x,
            phase=phase,
            density_anchor=density_anchor,
        )

    def observable_block(
        self,
        active_parameters: ActiveParameterSet,
        *,
        identity: ObservableIdentity,
        T: object,
        x: object,
        rho: object | None = None,
        P: object | None = None,
        phase: Phase | None = None,
    ) -> ObservableBlock:
        if identity == "eos.fixed-state" and rho is not None and P is None and phase is None:
            return _EosModel._fixed_state_observable_block(
                self, active_parameters, T=T, rho=rho, x=x
            )
        if identity == "eos.fixed-pressure-fugacity" and rho is None and P is not None and phase:
            return _EosModel._fixed_pressure_fugacity_observable_block(
                self, active_parameters, T=T, P=P, x=x, phase=phase
            )
        raise _eos_error("invalid-observable-inputs", identity=identity)

    def reference(
        self,
        *,
        T: object,
        P: object,
        solvent_x: object,
        phase: Literal["liquid"],
        ion_component_ids: tuple[str, str],
        formula_stoichiometry: tuple[int, int],
        reference_convention: Literal["molality-infinite-dilution"],
    ) -> ElectrolyteReference:
        return _EosModel._reference(
            self,
            T=T,
            P=P,
            solvent_x=solvent_x,
            phase=phase,
            ion_component_ids=ion_component_ids,
            formula_stoichiometry=formula_stoichiometry,
            reference_convention=reference_convention,
        )


def evaluate_observable(
    model: Pure | Mixture,
    active_parameters: ActiveParameterSet,
    descriptor: ObservableDescriptor,
    state: Mapping[str, float],
    *,
    composition: tuple[float, ...] | None,
) -> ObservableBlock:
    """Evaluate an EOS-owned block from canonical source-state fields."""

    if descriptor.parameter_fingerprint != model.parameter_fingerprint:
        raise _eos_error("invalid-observable-descriptor", reason="different-parameters")
    if descriptor.active_parameter_identities != active_parameters.identities:
        raise _eos_error("invalid-active-parameter-set", reason="different-identities")
    family = next(
        (
            candidate
            for candidate in OBSERVABLE_FAMILIES.values()
            if candidate.descriptor_identity == descriptor.identity
            and ("pure" if isinstance(model, Pure) else "mixture") in candidate.model_kinds
        ),
        None,
    )
    if family is None or set(state) != set(family.state_fields):
        raise _eos_error(
            "invalid-observable-inputs",
            identity=descriptor.identity,
            fields=tuple(sorted(state)),
        )
    temperature = state["temperature_k"] * unit_registry.kelvin
    if descriptor.identity == "eos.pure-saturation":
        if not isinstance(model, Pure) or composition is not None:
            raise _eos_error("invalid-observable-inputs", identity=descriptor.identity)
        return model.observable_block(
            active_parameters, identity=descriptor.identity, T=temperature
        )
    density = state["molar_density_mol_per_m3"] * unit_registry.mole / unit_registry.meter**3
    if isinstance(model, Pure):
        if composition is not None:
            raise _eos_error("invalid-observable-inputs", identity=descriptor.identity)
        return model.observable_block(
            active_parameters,
            identity=descriptor.identity,
            T=temperature,
            rho=density,
        )
    if composition is None:
        raise _eos_error("invalid-observable-inputs", identity=descriptor.identity)
    return model.observable_block(
        active_parameters,
        identity=descriptor.identity,
        T=temperature,
        rho=density,
        x=composition,
    )


def _saturation_diagnostics(payload: Mapping[str, object]) -> SaturationDiagnostics:
    def optional_float(key: str) -> float | None:
        value = payload[key]
        if value is None:
            return None
        numeric = float(cast(float, value))
        return numeric if math.isfinite(numeric) else None

    def optional_bounds(key: str) -> tuple[float, float] | None:
        value = payload[key]
        if value is None:
            return None
        return cast(
            tuple[float, float],
            tuple(float(item) for item in cast(Iterable[SupportsFloat], value)),
        )

    candidate_available = bool(payload["candidate_available"])
    sensitivity_status = str(payload["sensitivity_status"])
    return SaturationDiagnostics(
        candidate_available=candidate_available,
        optimizer_converged=bool(payload["optimizer_converged"]),
        optimizer_status=str(payload["optimizer_status"]),
        vapor_density_bounds_mol_m3=optional_bounds("vapor_density_bounds_mol_m3"),
        liquid_density_bounds_mol_m3=optional_bounds("liquid_density_bounds_mol_m3"),
        feasibility_measure=optional_float("feasibility_measure"),
        feasibility_tolerance=optional_float("feasibility_tolerance"),
        physical_accepted=bool(payload["physical_accepted"]),
        physical_status=str(payload["physical_status"]),
        physical_status_count=(
            None
            if not candidate_available or payload["physical_status_count"] is None
            else int(cast(int, payload["physical_status_count"]))
        ),
        pressure_relative_residual=optional_float("pressure_relative_residual"),
        vapor_pressure_relative_residual=optional_float("vapor_pressure_relative_residual"),
        liquid_pressure_relative_residual=optional_float("liquid_pressure_relative_residual"),
        chemical_potential_absolute_residual=optional_float("chemical_potential_absolute_residual"),
        phase_density_distance=optional_float("phase_density_distance"),
        vapor_residual_density_derivative_m3_per_mol=optional_float(
            "vapor_residual_density_derivative_m3_per_mol"
        ),
        liquid_residual_density_derivative_m3_per_mol=optional_float(
            "liquid_residual_density_derivative_m3_per_mol"
        ),
        rejection_status=cast(
            Literal["none", "solve-failure", "root-swap", "coalescence", "branch-loss"],
            str(payload["rejection_status"]),
        ),
        exact_derivatives=bool(payload["exact_derivatives"]),
        sensitivity_status=sensitivity_status,
        sensitivity_failure_reason=str(payload["sensitivity_failure_reason"]),
        sensitivity_rank=(
            None
            if sensitivity_status == "not_requested" or payload["sensitivity_rank"] is None
            else int(cast(int, payload["sensitivity_rank"]))
        ),
        sensitivity_dimension=(
            None
            if sensitivity_status == "not_requested" or payload["sensitivity_dimension"] is None
            else int(cast(int, payload["sensitivity_dimension"]))
        ),
        sensitivity_condition_number=optional_float("sensitivity_condition_number"),
        sensitivity_backward_error=optional_float("sensitivity_backward_error"),
        globality_certificate=bool(payload["globality_certificate"]),
        failure_reason=str(payload["failure_reason"]),
    )


def _saturation_phase(payload: Mapping[str, object]) -> SaturationPhase:
    raw_chemical_potential = payload["chemical_potential_over_rt"]
    if type(raw_chemical_potential) is not float:
        raise RuntimeError("native pure saturation chemical potential must be a scalar float")
    return SaturationPhase(
        amount_mol=float(cast(float, payload["amount_mol"])),
        volume_m3=float(cast(float, payload["volume_m3"])),
        molar_density_mol_m3=float(cast(float, payload["molar_density_mol_m3"])),
        pressure_pa=float(cast(float, payload["pressure_pa"])),
        chemical_potential_over_rt=raw_chemical_potential,
    )


def _failed_saturation_diagnostics(reason: str) -> dict[str, object]:
    return {
        "candidate_available": False,
        "optimizer_converged": False,
        "optimizer_status": "native_exception",
        "vapor_density_bounds_mol_m3": None,
        "liquid_density_bounds_mol_m3": None,
        "feasibility_measure": None,
        "feasibility_tolerance": None,
        "physical_accepted": False,
        "physical_status": "not_adjudicated",
        "physical_status_count": None,
        "pressure_relative_residual": None,
        "vapor_pressure_relative_residual": None,
        "liquid_pressure_relative_residual": None,
        "chemical_potential_absolute_residual": None,
        "phase_density_distance": None,
        "vapor_residual_density_derivative_m3_per_mol": None,
        "liquid_residual_density_derivative_m3_per_mol": None,
        "rejection_status": "solve-failure",
        "exact_derivatives": False,
        "sensitivity_status": "unavailable",
        "sensitivity_failure_reason": "native_exception",
        "sensitivity_rank": None,
        "sensitivity_dimension": None,
        "sensitivity_condition_number": None,
        "sensitivity_backward_error": None,
        "globality_certificate": False,
        "failure_reason": reason,
    }


def _native_topology_fingerprint(model: _EosModel) -> str:
    resolved = model._model
    charges = resolved.charge_number or (0,) * len(resolved.component_ids)
    neutral_phase = (
        resolved.topology == "neutral"
        and len(resolved.component_ids) in (1, 2)
        and not any(charges)
        and not resolved.site_ids
        and not resolved.association_pairs
    )
    payload: dict[str, object]
    if neutral_phase:
        pure = len(resolved.component_ids) == 1
        payload = {
            "schema": (
                "epcsaft-native-neutral-pure-topology-v1"
                if pure
                else "epcsaft-native-neutral-binary-topology-v1"
            ),
            "component_ids": resolved.component_ids,
            "topology": resolved.topology,
            "charge_number": charges,
            "site_ids": resolved.site_ids,
            "site_multiplicity": resolved.site_multiplicity,
            "association_pair_count": len(resolved.association_pairs),
            "helmholtz_terms": (
                "ideal",
                "hard-chain",
                "dispersion",
            ),
        }
    else:
        payload = {
            "schema": "epcsaft-native-model-topology-v1",
            "component_ids": resolved.component_ids,
            "topology": resolved.topology,
            "charge_number": charges,
            "site_ids": resolved.site_ids,
            "site_multiplicity": resolved.site_multiplicity,
            "association_endpoints": tuple(
                (pair.endpoint_a, pair.endpoint_b) for pair in resolved.association_pairs
            ),
            "relative_permittivity_formulation": (resolved.relative_permittivity_formulation),
        }
    encoded = json.dumps(payload, separators=(",", ":"), sort_keys=True).encode("utf-8")
    return f"sha256:{hashlib.sha256(encoded).hexdigest()}"


def _native_identity_fingerprint(schema: str, payload: object) -> str:
    encoded = json.dumps(
        {"schema": schema, "payload": payload},
        separators=(",", ":"),
        sort_keys=True,
    ).encode("utf-8")
    return f"sha256:{hashlib.sha256(encoded).hexdigest()}"


def _state_magnitude(value: object, unit: str, field: str) -> float:
    try:
        quantity = _quantity(value, unit, field=field)
    except ParameterError as exc:
        details = dict(exc.context)
        if exc.code == "invalid-value":
            rejected = details.get("value")
            if rejected is not None:
                try:
                    rejected_number = float(rejected)
                except ValueError:
                    pass
                else:
                    if not math.isfinite(rejected_number):
                        raise _eos_error(
                            "invalid-state",
                            field=field,
                            reason="nonfinite",
                            value=rejected,
                        ) from exc
            context = {"field": field, "reason": exc.code}
            for key in ("value", "value_type"):
                if key in details:
                    context[key] = details[key]
            raise _eos_error("invalid-state", **context) from exc
        raise _eos_error(
            "invalid-state",
            field=field,
            reason=exc.code,
            expected_unit=details.get("expected", unit),
        ) from exc
    magnitude = float(quantity.magnitude)
    if not math.isfinite(magnitude):
        raise _eos_error(
            "invalid-state",
            field=field,
            reason="nonfinite",
            value=magnitude,
        )
    if magnitude <= 0.0:
        raise _eos_error("invalid-state", field=field, reason="nonpositive", value=magnitude)
    return magnitude


@dataclass(frozen=True, slots=True)
class Fugacity:
    """Immutable component fugacity projections in the mixture component order."""

    component_ids: tuple[str, ...]
    ln_coefficient: tuple[float, ...]
    coefficient: tuple[float, ...]
    value: tuple[pint.Quantity[float], ...]


def _basis(value: str) -> PropertyBasis:
    if value not in ("molar", "reduced"):
        raise _eos_error("invalid-property-basis", value=value)
    return value  # type: ignore[return-value]


@dataclass(frozen=True, slots=True)
class State:
    """Immutable properties plus the resolved state and parameter identity."""

    component_ids: tuple[str, ...]
    temperature: pint.Quantity[float]
    mole_fractions: tuple[float, ...]
    parameter_fingerprint: str
    topology_fingerprint: str
    component_order_fingerprint: str
    artifact_fingerprint: str
    requested_phase: Phase | None
    certified_branch: Mapping[str, object] | None
    state_identity: str
    hard_chain: float
    dispersion: float
    association: float
    debye_huckel: float
    born: float
    bulk_relative_permittivity: float | None
    _residual_helmholtz: float
    compressibility_factor: float
    pressure: pint.Quantity[float]
    residual_chemical_potential_over_rt: tuple[float, ...]
    _log_fugacity_coefficient: tuple[float, ...] | None
    _residual_entropy: pint.Quantity[float]
    _residual_internal_energy: pint.Quantity[float]
    _residual_enthalpy: pint.Quantity[float]
    _residual_isochoric_heat_capacity: pint.Quantity[float]
    _residual_isobaric_heat_capacity: pint.Quantity[float] | None
    molar_density: pint.Quantity[float]
    density_diagnostics: DensityDiagnostics | None
    fixed_pressure_composition_derivatives: Mapping[str, object] | None = None
    _activity: Activity | None = None

    def _rt(self) -> pint.Quantity[float]:
        return cast(
            pint.Quantity[float],
            _MOLAR_GAS_CONSTANT
            * self.temperature.to("kelvin").magnitude
            * unit_registry.joule
            / unit_registry.mole,
        )

    def ares(self, *, basis: PropertyBasis = "molar") -> pint.Quantity[float] | float:
        """Residual molar Helmholtz energy, or its dimensionless ``Ares/(nRT)``."""
        return (
            self._residual_helmholtz * self._rt()
            if _basis(basis) == "molar"
            else self._residual_helmholtz
        )

    def ures(self, *, basis: PropertyBasis = "molar") -> pint.Quantity[float] | float:
        """Residual molar internal energy, or its dimensionless ``Ures/(nRT)``."""
        return (
            self._residual_internal_energy
            if _basis(basis) == "molar"
            else (self._residual_internal_energy / self._rt()).to("").magnitude
        )

    def hres(self, *, basis: PropertyBasis = "molar") -> pint.Quantity[float] | float:
        """Residual molar enthalpy, or its dimensionless ``Hres/(nRT)``."""
        return (
            self._residual_enthalpy
            if _basis(basis) == "molar"
            else (self._residual_enthalpy / self._rt()).to("").magnitude
        )

    def gres(self, *, basis: PropertyBasis = "molar") -> pint.Quantity[float] | float:
        """Residual molar Gibbs energy from the exact chemical-potential projection."""
        reduced = math.fsum(
            fraction * chemical_potential
            for fraction, chemical_potential in zip(
                self.mole_fractions,
                self.residual_chemical_potential_over_rt,
                strict=True,
            )
        )
        return reduced * self._rt() if _basis(basis) == "molar" else reduced

    def sres(self, *, basis: PropertyBasis = "molar") -> pint.Quantity[float] | float:
        """Residual molar entropy, or its dimensionless ``Sres/(nR)``."""
        return (
            self._residual_entropy
            if _basis(basis) == "molar"
            else (
                self._residual_entropy
                / (
                    _MOLAR_GAS_CONSTANT
                    * unit_registry.joule
                    / unit_registry.mole
                    / unit_registry.kelvin
                )
            )
            .to("")
            .magnitude
        )

    def cvres(self, *, basis: PropertyBasis = "molar") -> pint.Quantity[float] | float:
        """Residual molar isochoric heat capacity, or its dimensionless value over R."""
        return (
            self._residual_isochoric_heat_capacity
            if _basis(basis) == "molar"
            else (
                self._residual_isochoric_heat_capacity
                / (
                    _MOLAR_GAS_CONSTANT
                    * unit_registry.joule
                    / unit_registry.mole
                    / unit_registry.kelvin
                )
            )
            .to("")
            .magnitude
        )

    def cpres(self, *, basis: PropertyBasis = "molar") -> pint.Quantity[float] | float | None:
        """Residual molar isobaric heat capacity, or its dimensionless value over R."""
        if self._residual_isobaric_heat_capacity is None:
            return None
        return (
            self._residual_isobaric_heat_capacity
            if _basis(basis) == "molar"
            else (
                self._residual_isobaric_heat_capacity
                / (
                    _MOLAR_GAS_CONSTANT
                    * unit_registry.joule
                    / unit_registry.mole
                    / unit_registry.kelvin
                )
            )
            .to("")
            .magnitude
        )

    @property
    def fugacity(self) -> Fugacity | None:
        """Component fugacity coefficients and fugacities, when defined."""
        if self._log_fugacity_coefficient is None:
            return None
        coefficients = tuple(math.exp(value) for value in self._log_fugacity_coefficient)
        values = tuple(
            fraction * coefficient * self.pressure
            for fraction, coefficient in zip(self.mole_fractions, coefficients, strict=True)
        )
        return Fugacity(
            component_ids=self.component_ids,
            ln_coefficient=self._log_fugacity_coefficient,
            coefficient=coefficients,
            value=values,
        )

    @property
    def activity(self) -> Activity | None:
        """Electrolyte activity projections, present only for reference-created states."""
        return self._activity


@dataclass(frozen=True, slots=True)
class DensityDiagnostics:
    bracket_lower: pint.Quantity[float]
    bracket_upper: pint.Quantity[float]
    pressure_residual: pint.Quantity[float]
    requested_pressure: pint.Quantity[float]
    pressure_density_derivative: pint.Quantity[float]
    iterations: int
    stable: bool
    branch: str
    stable_root_count: int | None
    selected_stable_root_index: int | None
    anchor_status: Literal["not_supplied", "accepted"]
    anchor_density: pint.Quantity[float] | None
    anchor_source_composition: tuple[float, ...] | None
    anchor_source_state_identity: str | None
    requested_phase: Phase | None
    certified_branch: Mapping[str, object] | None
    certificate_kind: str
    local_search_outcome: str
    local_root_count: int | None
    reachability_envelope: Mapping[str, object] | None


def _certified_branch(
    diagnostics: _native_core.NativeDensityDiagnostics | None,
    phase: str | None,
) -> Mapping[str, object] | None:
    if diagnostics is None or phase is None:
        return None
    return MappingProxyType(
        {
            "schema": "epcsaft-certified-density-branch-v1",
            "requested_phase": phase,
            "stable_root_index": diagnostics.certified_stable_root_index,
            "origin_molar_density_hex": diagnostics.certified_origin_molar_density.hex(),
        }
    )


def _state_identity(
    *,
    artifact_fingerprint: str,
    parameter_fingerprint: str,
    topology_fingerprint: str,
    component_order_fingerprint: str,
    component_ids: tuple[str, ...],
    temperature_k: float,
    pressure_pa: float,
    mole_fractions: tuple[float, ...],
    molar_density: float,
    requested_phase: str | None,
    certified_branch: Mapping[str, object] | None,
) -> str:
    encoded = json.dumps(
        {
            "schema": "epcsaft-state-identity-v1",
            "artifact_fingerprint": artifact_fingerprint,
            "parameter_fingerprint": parameter_fingerprint,
            "topology_fingerprint": topology_fingerprint,
            "component_order_fingerprint": component_order_fingerprint,
            "component_ids": component_ids,
            "temperature_k": temperature_k.hex(),
            "pressure_pa": pressure_pa.hex(),
            "mole_fractions": tuple(value.hex() for value in mole_fractions),
            "molar_density_mol_per_m3": molar_density.hex(),
            "requested_phase": requested_phase,
            "certified_branch": None if certified_branch is None else dict(certified_branch),
        },
        allow_nan=False,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("utf-8")
    return f"sha256:{hashlib.sha256(encoded).hexdigest()}"


def _validate_density_anchor(
    owner: _EosModel,
    candidate: object,
    temperature_k: float,
    pressure_pa: float,
    phase: str,
) -> State:
    if not isinstance(candidate, State):
        raise _eos_error("invalid-density-anchor", reason="expected-State")
    try:
        density = float(candidate.molar_density.to("mole / meter ** 3").magnitude)
        certified = cast(Mapping[str, object], candidate.certified_branch)
        origin_density = float.fromhex(str(certified["origin_molar_density_hex"]))
        diagnostics = candidate.density_diagnostics
        if diagnostics is None:
            raise ValueError("missing density diagnostics")
        requested_pressure = float(diagnostics.requested_pressure.to("pascal").magnitude)
        expected_identity = _state_identity(
            artifact_fingerprint=candidate.artifact_fingerprint,
            parameter_fingerprint=candidate.parameter_fingerprint,
            topology_fingerprint=candidate.topology_fingerprint,
            component_order_fingerprint=candidate.component_order_fingerprint,
            component_ids=candidate.component_ids,
            temperature_k=float(candidate.temperature.to("kelvin").magnitude),
            pressure_pa=requested_pressure,
            mole_fractions=candidate.mole_fractions,
            molar_density=density,
            requested_phase=candidate.requested_phase,
            certified_branch=certified,
        )
    except (AttributeError, KeyError, TypeError, ValueError):
        raise _eos_error("invalid-density-anchor", reason="missing-or-invalid-identity") from None
    if (
        not math.isfinite(density)
        or density <= 0.0
        or not math.isfinite(origin_density)
        or origin_density <= 0.0
        or candidate.state_identity != expected_identity
        or certified.get("schema") != "epcsaft-certified-density-branch-v1"
    ):
        raise _eos_error("invalid-density-anchor", reason="missing-or-invalid-identity")
    if (
        candidate.parameter_fingerprint != owner._parameter_fingerprint
        or candidate.component_ids != owner._component_ids
        or candidate.topology_fingerprint != owner._topology_fingerprint
        or candidate.component_order_fingerprint != owner._component_order_fingerprint
        or candidate.artifact_fingerprint != owner._artifact_fingerprint
    ):
        raise _eos_error("density-anchor-model-mismatch")
    anchor_temperature = float(candidate.temperature.to("kelvin").magnitude)
    if anchor_temperature != temperature_k or requested_pressure != pressure_pa:
        raise _eos_error("density-anchor-state-mismatch")
    if candidate.requested_phase != phase or certified.get("requested_phase") != phase:
        raise _eos_error("density-anchor-phase-mismatch")
    return candidate


def _fixed_pressure_composition_derivatives(
    owner: _EosModel,
    evaluator: _native_core.NativeModel,
    temperature_k: float,
    fractions: tuple[float, ...],
    molar_density: float,
) -> Mapping[str, object]:
    component_ids = owner._component_ids
    charge_numbers = owner._model.charge_number
    charged = any(charge != 0 for charge in charge_numbers)
    constraints = ("normalization", "electroneutrality") if charged else ("normalization",)
    artifact = _python_native_artifact_fingerprint()

    def unavailable(
        code: str,
        condition: float | None = None,
        diagnostic: str | None = None,
    ) -> Mapping[str, object]:
        return MappingProxyType(
            {
                "status": "non_evaluable",
                "failure": NonEvaluableTrial(cast(Any, code), diagnostic or code),
                "condition_measure": condition,
                "artifact_fingerprint": artifact,
            }
        )

    if any(fraction <= 0.0 for fraction in fractions):
        return unavailable("eos_domain_rejection")
    try:
        projected = evaluator.evaluate_fixed_pressure_composition_derivatives(
            temperature_k,
            fractions,
            1.0 / molar_density,
        )
    except (RuntimeError, ValueError) as exc:
        return unavailable("native_residual_failure", diagnostic=str(exc))

    condition = float(projected["condition_measure"])
    condition_value = condition if math.isfinite(condition) and condition > 0.0 else None
    failure_code = str(projected["failure_code"])
    if failure_code:
        return unavailable(failure_code, condition_value)
    coordinate_indices = tuple(
        int(index) for index in cast(Iterable[int], projected["coordinate_indices"])
    )
    dependent_indices = tuple(
        int(index) for index in cast(Iterable[int], projected["dependent_indices"])
    )
    coordinate_count = len(coordinate_indices)

    def rows(name: str) -> tuple[tuple[float, ...], ...]:
        flat = tuple(float(value) for value in cast(Iterable[float], projected[name]))
        return tuple(
            flat[component * coordinate_count : (component + 1) * coordinate_count]
            for component in range(len(fractions))
        )

    return MappingProxyType(
        {
            "status": "available",
            "failure": None,
            "constraint_identities": constraints,
            "coordinate_component_ids": tuple(component_ids[index] for index in coordinate_indices),
            "dependent_component_ids": tuple(component_ids[index] for index in dependent_indices),
            "log_composition_basis": tuple(
                tuple(value / fractions[row] for value in values)
                for row, values in enumerate(rows("composition_derivatives"))
            ),
            "chemical_potential_derivatives_over_rt": rows("chemical_derivatives"),
            "molar_density_derivatives_mol_per_m3": tuple(
                float(value) for value in cast(Iterable[float], projected["density_derivatives"])
            ),
            "condition_measure": condition_value,
            "artifact_fingerprint": artifact,
        }
    )


__all__ = (
    "OBSERVABLE_FAMILIES",
    "OBSERVABLE_TARGETS",
    "DensityDiagnostics",
    "EosError",
    "Fugacity",
    "Mixture",
    "ObservableBasis",
    "ObservableBlock",
    "ObservableBlockStatus",
    "ObservableDescriptor",
    "ObservableFamily",
    "ObservableIdentity",
    "ObservableRowStatus",
    "ObservableState",
    "ObservableTarget",
    "Phase",
    "PropertyBasis",
    "Pure",
    "SaturationBranchEvidence",
    "SaturationDiagnostics",
    "SaturationError",
    "SaturationPhase",
    "SaturationResult",
    "SaturationStart",
    "State",
    "evaluate_observable",
)
