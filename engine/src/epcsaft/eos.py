from __future__ import annotations

import hashlib
import json
import math
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import TYPE_CHECKING, Any, Literal, SupportsFloat, cast, overload

import pint

import epcsaft._core as _native_core

from ._resolve import _resolve_model, _ResolvedModel
from .parameters import ActiveParameterSet, Parameters
from .records import (
    ConstantCorrelation,
    ConstantPlusSumOfExponentialsCorrelation,
    SingleParameterRecord,
)
from .units import ParameterError, _error, _quantity, unit_registry

if TYPE_CHECKING:
    from .reference import Activity, ElectrolyteReference


type Phase = Literal["vapor", "liquid"]
type PropertyBasis = Literal["molar", "reduced"]
_MOLAR_GAS_CONSTANT = 8.31446261815324
_native_transport_factory = _native_core.native_transport


@lru_cache(maxsize=1)
def _python_native_artifact_fingerprint() -> str:
    return _file_fingerprint(Path(_native_core.__file__))


def _file_fingerprint(path: Path) -> str:
    return f"sha256:{hashlib.sha256(path.read_bytes()).hexdigest()}"


def _magnitude(
    record: SingleParameterRecord | ConstantCorrelation | ConstantPlusSumOfExponentialsCorrelation,
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
    """One pure phase on a locally certified saturation boundary."""

    amount_mol: float
    volume_m3: float
    molar_density_mol_m3: float
    pressure_pa: float
    chemical_potential_over_rt: float


@dataclass(frozen=True, slots=True)
class SaturationAttempt:
    role: str
    initial_guess: tuple[float, float, float]
    solver_converged: bool
    solver_status: str
    iterations: int
    constraint_violation: float
    callback_error: str


@dataclass(frozen=True, slots=True)
class SaturationDiagnostics:
    solver_converged: bool
    solver_status: str
    iterations: int
    attempts: int
    attempt_log: tuple[SaturationAttempt, ...]
    solver_lower_bounds: tuple[float, float, float]
    solver_upper_bounds: tuple[float, float, float]
    solver_constraint_violation: float
    numerical_converged: bool
    confirmation_solves: int
    confirmation_max_relative_difference: float
    physical_accepted: bool
    pressure_relative_residual: float
    chemical_potential_absolute_residual: float
    phase_density_distance: float
    exact_derivatives: bool
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


@dataclass(frozen=True, slots=True)
class ParameterCapability:
    observation_family: str
    parameter_identities: tuple[str, ...]
    parameter_values: tuple[float, ...]
    parameter_units: tuple[str, ...]
    model_domain: str
    temperature_domain_k: tuple[float | None, float | None]
    value_labels: tuple[str, ...]
    value_units: tuple[str, ...]
    value_shape: tuple[int, ...]
    derivative_shape: tuple[int, int]
    parameter_fingerprint: str
    active_parameter_fingerprint: str
    topology_fingerprint: str


@dataclass(frozen=True, slots=True)
class ParameterDerivativeResult:
    capability: ParameterCapability
    values: tuple[float, ...]
    derivatives: tuple[tuple[float, ...], ...]


def _direct_component_identity(identity: str, family: str) -> bool:
    parts = identity.split("/")
    return len(parts) == 3 and parts[0] == "component" and parts[2] == family


def _segment_diameter_parameter_identity(identity: str, family: str) -> bool:
    parts = identity.split("/")
    if len(parts) == 3:
        return family == "segment_diameter" and parts[0] == "component" and parts[2] == family
    if len(parts) == 5:
        return (
            family == "segment_diameter"
            and parts[0] == "component"
            and parts[2] == "segment_diameter"
            and parts[4] == "constant"
        )
    return (
        len(parts) == 6
        and parts[0] == "component"
        and parts[2] == "segment_diameter"
        and parts[4].startswith("term-")
        and (
            (family == "segment_diameter" and parts[5] == "amplitude")
            or (
                family == "exponential_temperature_coefficient"
                and parts[5] == "exponent_coefficient"
            )
        )
    )


def _association_identity(identity: str, family: str) -> bool:
    parts = identity.split("/")
    return (
        len(parts) == 6
        and parts[0] == "association"
        and parts[5] == family.removeprefix("association_")
    )


def _pair_identity(identity: str, family: str) -> bool:
    parts = identity.split("/")
    return len(parts) == 4 and parts[0] == "pair" and parts[3] == family


def _cross_association_identity(identity: str, family: str) -> bool:
    parts = identity.split("/")
    return (
        family == "k_hb_ij"
        and len(parts) == 7
        and parts[0] == "pair"
        and parts[3] == "association"
        and parts[6] == family
    )


def _require_supported_active_parameters(
    active_parameters: ActiveParameterSet,
    observation: str,
    admitted: Iterable[bool],
) -> None:
    unsupported = tuple(
        coordinate.identity
        for coordinate, supported in zip(active_parameters.coordinates, admitted, strict=True)
        if not supported
    )
    if unsupported:
        raise _eos_error(
            "unsupported-parameter-derivative",
            observation=observation,
            identities=unsupported,
        )


class SaturationError(RuntimeError):
    """Raised when no seed produces a confirmed physical saturation boundary."""

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
        "ionic-ssm-ds",
        "ionic-original-born",
    }:
        molar_mass = [0.0 if value is None else value for value in model.molar_mass_kg_per_mol]
        charge_number = list(model.charge_number)
        relative_permittivity = [
            0.0 if record is None else _magnitude(record, "dimensionless")
            for record in model.relative_permittivity
        ]
        if model.topology == "ionic-ssm-ds":
            dielectric_suppression = (
                0.0
                if model.ion_fraction_suppression_coefficient is None
                else model.ion_fraction_suppression_coefficient
            )
            born_diameter = [
                0.0 if value is None else value for value in model.born_diameter_angstrom
            ]
            born_diameter_sigma_scale = [0.0] * components
            solvation_factor = [
                _magnitude(record, "dimensionless")
                for record in model.solvation_factor
                if record is not None
            ]
            ionic_permittivity = model.ionic_region_relative_permittivity
            electrolyte_model = 1
            assert ionic_permittivity is not None
            assert len(solvation_factor) == components
        elif model.topology == "ionic-fully-dissociated-debye-huckel":
            dielectric_suppression = 0.0
            born_diameter = [0.0] * components
            born_diameter_sigma_scale = [0.0] * components
            solvation_factor = [0.0] * components
            ionic_permittivity = 0.0
            electrolyte_model = 3
        else:
            dielectric_suppression = 0.0
            born_diameter = [0.0] * components
            born_diameter_sigma_scale = [
                0.0 if charge == 0 else 1.0 for charge in charge_number
            ]
            solvation_factor = [1.0] * components
            ionic_permittivity = 1.0
            electrolyte_model = 2
    else:
        molar_mass = [0.0] * components
        charge_number = [0] * components
        relative_permittivity = [0.0] * components
        dielectric_suppression = 0.0
        born_diameter = [0.0] * components
        born_diameter_sigma_scale = [0.0] * components
        solvation_factor = [0.0] * components
        ionic_permittivity = 0.0
        electrolyte_model = 0

    return (
        list(model.segment_count),
        sigma_constant,
        term_component,
        term_amplitude,
        term_exponent,
        list(model.dispersion_energy_over_k),
        list(model.k_ij),
        list(model.l_ij),
        [component_position[component_id] for component_id, _ in model.site_ids],
        list(model.site_multiplicity),
        association_active,
        association_energy,
        association_volume,
        molar_mass,
        charge_number,
        relative_permittivity,
        dielectric_suppression,
        born_diameter,
        born_diameter_sigma_scale,
        solvation_factor,
        ionic_permittivity,
        electrolyte_model,
        list(model.dielectric_volume_m3_per_mol) or [0.0] * components,
        list(model.dielectric_temperature_k) or [0.0] * components,
        list(model.ion_specific_suppression_coefficient) or [0.0] * components,
        list(model.dipole_scaling_debye2) or [0.0] * components,
        list(model.polarizability_scaling_angstrom3) or [0.0] * components,
        list(model.correlation_integral_parameter) or [0.0] * components,
        relative_permittivity_formulation_codes[model.relative_permittivity_formulation],
        list(model.dipole_moment_debye) or [0.0] * components,
        list(model.quadrupole_moment_debye_angstrom) or [0.0] * components,
        0 if model.polar_formulation is None else 1,
    )


@dataclass(frozen=True, slots=True, init=False)
class _EosModel:
    _component_ids: tuple[str, ...]
    _parameter_fingerprint: str
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
        native = _native_core.NativeModel(*native_arguments).with_parameter_identities(
            parameters.component_ids,
            tuple(f"{component}/{site}" for component, site in model.site_ids),
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

    def _active_records(
        self, active_parameters: ActiveParameterSet
    ) -> tuple[_native_core.NativeParameterRecord, ...]:
        if not isinstance(active_parameters, ActiveParameterSet):
            raise _eos_error("invalid-active-parameter-set", expected="ActiveParameterSet")
        if active_parameters.parameters.fingerprint != self._parameter_fingerprint:
            raise _eos_error("invalid-active-parameter-set", reason="different-parameters")
        records = []
        for coordinate in active_parameters.coordinates:
            records.append(
                _native_core.NativeParameterRecord(coordinate.identity, coordinate.value)
            )
        return tuple(records)

    def _pair_coordinates(
        self, active_parameters: ActiveParameterSet
    ) -> tuple[_native_core.NativePairCoordinate, ...]:
        result = []
        for coordinate in active_parameters.coordinates:
            applicability = coordinate.applicability
            if coordinate.family not in {"k_ij", "l_ij", "k_hb_ij"} or applicability is None:
                raise _eos_error(
                    "unsupported-parameter-derivative",
                    observation="fixed-state-helmholtz",
                    identity=coordinate.identity,
                    reason="unresolved-pair-application",
                )
            roles = applicability.association_role_pair or ("", "")
            result.append(
                _native_core.NativePairCoordinate(
                    coordinate.identity,
                    coordinate.family,
                    applicability.component_indices[0],
                    applicability.component_indices[1],
                    roles[0],
                    roles[1],
                    applicability.association_edges,
                    applicability.association_energy_scales,
                    coordinate.value,
                )
            )
        return tuple(result)

    def _capability(
        self,
        active_parameters: ActiveParameterSet,
        observation_family: str,
        labels: tuple[str, ...],
        units: tuple[str, ...],
    ) -> ParameterCapability:
        return ParameterCapability(
            observation_family=observation_family,
            parameter_identities=active_parameters.identities,
            parameter_values=active_parameters.values,
            parameter_units=tuple(coordinate.unit for coordinate in active_parameters.coordinates),
            model_domain=self._model.topology,
            temperature_domain_k=(
                self._model.temperature_min_k,
                self._model.temperature_max_k,
            ),
            value_labels=labels,
            value_units=units,
            value_shape=(len(labels),),
            derivative_shape=(len(labels), len(active_parameters.coordinates)),
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
        polar_boundary = self._model.polar_formulation is not None and any(
            fraction <= 0.0
            and (
                self._model.dipole_moment_debye[component] != 0.0
                or self._model.quadrupole_moment_debye_angstrom[component] != 0.0
            )
            for component, fraction in enumerate(fractions)
        )
        if polar_boundary:
            raise _eos_error(
                "invalid-state",
                field="x",
                reason="polar-composition-boundary",
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

        if self._model.topology in {"ionic-ssm-ds", "ionic-original-born"}:
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

        if rho is not None:
            density_molar = _state_magnitude(rho, "mole / meter ** 3", "rho")
            try:
                native = self._native.evaluate(temperature_k, density_molar, fractions)
            except (RuntimeError, ValueError) as exc:
                if str(exc) == "polar-domain-error":
                    raise _eos_error("polar-domain-error") from exc
                raise _eos_error("evaluation-failed", reason=str(exc)) from exc
        else:
            pressure_pa = _state_magnitude(P, "pascal", "P")
            try:
                native = self._native.evaluate_pressure(
                    temperature_k, pressure_pa, fractions, phase or ""
                )
            except (RuntimeError, ValueError) as exc:
                if str(exc) == "polar-domain-error":
                    raise _eos_error("polar-domain-error") from exc
                raise _eos_error("evaluation-failed", reason=str(exc)) from exc
        diagnostics = native.density_diagnostics
        return State(
            component_ids=self._component_ids,
            temperature=temperature_k * unit_registry.kelvin,
            mole_fractions=fractions,
            parameter_fingerprint=self._parameter_fingerprint,
            hard_chain=native.hard_chain,
            dispersion=native.dispersion,
            association=native.association,
            ares_polar=native.ares_polar,
            debye_huckel=native.debye_huckel,
            born=native.born_ssm_ds,
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
                pressure_density_derivative=diagnostics.pressure_density_derivative
                * unit_registry.pascal
                * unit_registry.meter**3
                / unit_registry.mole,
                iterations=diagnostics.iterations,
                stable=diagnostics.stable,
                branch=diagnostics.branch,
            ),
        )

    def _phase_parameter_derivatives(
        self,
        active_parameters: ActiveParameterSet,
        *,
        T: object,
        rho: object,
        x: object,
    ) -> ParameterDerivativeResult:
        """Evaluate one fixed-density phase and exact selected-parameter derivatives."""

        associating = bool(self._model.association_pairs)
        pair_active = all(
            coordinate.family in {"k_ij", "l_ij", "k_hb_ij"}
            and coordinate.applicability is not None
            for coordinate in active_parameters.coordinates
        )
        if pair_active:
            admitted = (
                _pair_identity(coordinate.identity, coordinate.family)
                and coordinate.family in {"k_ij", "l_ij"}
                or _cross_association_identity(coordinate.identity, coordinate.family)
                for coordinate in active_parameters.coordinates
            )
        elif associating:
            admitted = (
                (
                    _direct_component_identity(coordinate.identity, coordinate.family)
                    and coordinate.family in {"segment_count", "dispersion_energy_over_k"}
                    or _segment_diameter_parameter_identity(coordinate.identity, coordinate.family)
                )
                or _association_identity(coordinate.identity, coordinate.family)
                and coordinate.family in {"association_energy_over_k", "association_volume"}
                for coordinate in active_parameters.coordinates
            )
        elif len(self._component_ids) == 1:
            admitted = (
                _direct_component_identity(coordinate.identity, coordinate.family)
                and coordinate.family in {"segment_count", "dispersion_energy_over_k"}
                or _segment_diameter_parameter_identity(coordinate.identity, coordinate.family)
                for coordinate in active_parameters.coordinates
            )
        else:
            admitted = (False for _ in active_parameters.coordinates)
        _require_supported_active_parameters(active_parameters, "fixed-state-helmholtz", admitted)
        records = self._active_records(active_parameters)
        state = self._state(T=T, rho=rho, x=x)
        temperature_k = float(state.temperature.to("kelvin").magnitude)
        density = float(state.molar_density.to("mole / meter ** 3").magnitude)
        amounts = state.mole_fractions
        volume = 1.0 / density
        parameter_count = len(records)
        component_count = len(amounts)
        parameter_coordinates: tuple[int, ...]
        block: Any
        try:
            if pair_active:
                block = self._native.evaluate_mixture_phase_pair_coordinates(
                    temperature_k,
                    amounts,
                    volume,
                    self._pair_coordinates(active_parameters),
                )
                width = component_count + 1 + parameter_count
                parameter_coordinates = tuple(range(component_count + 1, width))
            elif self._model.association_pairs:
                block = (
                    self._native.evaluate_fixed_topology_association_active_parameters(
                        temperature_k, 1.0, volume, records
                    )
                    if component_count == 1
                    else self._native.evaluate_fixed_topology_associating_mixture_active_parameters(
                        temperature_k, amounts, volume, records
                    )
                )
                width = component_count + 1 + parameter_count
                parameter_coordinates = tuple(range(component_count + 1, width))
            elif component_count == 1:
                block = self._native.evaluate_fixed_topology_association_active_parameters(
                    temperature_k, 1.0, volume, records
                )
                width = 2 + parameter_count
                parameter_coordinates = tuple(range(2, width))
            else:
                raise _eos_error(
                    "unsupported-parameter-derivative",
                    observation="fixed-state-helmholtz",
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
        capability = self._capability(
            active_parameters,
            "fixed-state-helmholtz",
            labels,
            ("1", "pascal", *("1" for _ in self._component_ids)),
        )
        return ParameterDerivativeResult(
            capability=capability,
            values=(
                float(block.helmholtz_over_rt_reference_amount),
                float(block.pressure_pa),
                *(gradient[component] for component in range(component_count)),
            ),
            derivatives=(
                helmholtz_derivatives,
                pressure_derivatives,
                *chemical_derivatives,
            ),
        )

    def _fugacity_parameter_derivatives(
        self,
        active_parameters: ActiveParameterSet,
        *,
        T: object,
        P: object,
        x: object,
        phase: Phase,
    ) -> ParameterDerivativeResult:
        """Evaluate fixed-pressure log fugacity and exact implicit derivatives."""

        _require_supported_active_parameters(
            active_parameters,
            "fixed-pressure-fugacity",
            (
                _direct_component_identity(coordinate.identity, coordinate.family)
                and coordinate.family in {"segment_diameter", "dispersion_energy_over_k"}
                for coordinate in active_parameters.coordinates
            ),
        )
        state = self._state(T=T, P=P, x=x, phase=phase)
        temperature_k = float(state.temperature.to("kelvin").magnitude)
        pressure_pa = float(state.pressure.to("pascal").magnitude)
        records = self._active_records(active_parameters)
        try:
            block = self._native.evaluate_fugacity_derivatives_active(
                temperature_k,
                pressure_pa,
                state.mole_fractions,
                phase,
                records,
            )
        except (RuntimeError, ValueError) as exc:
            raise _eos_error(
                "unsupported-parameter-derivative",
                observation="fixed-pressure-fugacity",
                reason=str(exc),
            ) from exc
        parameter_count = len(records)
        values = tuple(float(value) for value in block.log_fugacity_coefficient)
        flat = tuple(float(value) for value in block.parameter_derivatives)
        derivatives = tuple(
            flat[row * parameter_count : (row + 1) * parameter_count] for row in range(len(values))
        )
        capability = self._capability(
            active_parameters,
            "fixed-pressure-fugacity",
            tuple(f"log_fugacity_coefficient/{component}" for component in self._component_ids),
            tuple("1" for _ in self._component_ids),
        )
        return ParameterDerivativeResult(capability, values, derivatives)

    def parameter_capability(
        self,
        active_parameters: ActiveParameterSet,
        *,
        observation_family: str,
    ) -> ParameterCapability:
        if observation_family != "fixed_state_pressure":
            raise _eos_error(
                "unsupported-observation-family", observation_family=observation_family
            )
        self._active_records(active_parameters)
        associating = bool(self._model.association_pairs)
        if associating:
            admitted = (
                (
                    _direct_component_identity(coordinate.identity, coordinate.family)
                    and coordinate.family in {"segment_count", "dispersion_energy_over_k"}
                    or _segment_diameter_parameter_identity(coordinate.identity, coordinate.family)
                )
                or _association_identity(coordinate.identity, coordinate.family)
                and coordinate.family in {"association_energy_over_k", "association_volume"}
                for coordinate in active_parameters.coordinates
            )
        elif len(self._component_ids) == 1:
            admitted = (
                _direct_component_identity(coordinate.identity, coordinate.family)
                and coordinate.family in {"segment_count", "dispersion_energy_over_k"}
                or _segment_diameter_parameter_identity(coordinate.identity, coordinate.family)
                for coordinate in active_parameters.coordinates
            )
        else:
            if (
                len(self._component_ids) != 2
                or len(active_parameters.coordinates) != 1
                or any(self._model.charge_number)
            ):
                raise _eos_error(
                    "unsupported-parameter-derivative",
                    observation=observation_family,
                    reason="binary-pair-coordinate-required",
                )
            admitted = (
                _pair_identity(coordinate.identity, coordinate.family)
                and coordinate.family in {"k_ij", "l_ij"}
                for coordinate in active_parameters.coordinates
            )
        _require_supported_active_parameters(active_parameters, observation_family, admitted)
        return self._capability(
            active_parameters,
            observation_family,
            ("pressure_pa",),
            ("pascal",),
        )

    def _reference(
        self,
        *,
        T: object,
        P: object,
        solvent_x: object,
        phase: Literal["liquid"],
        ion_component_ids: tuple[str, str] | None = None,
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
        T: pint.Quantity[Any],
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

    def phase_parameter_derivatives(
        self,
        active_parameters: ActiveParameterSet,
        *,
        T: object,
        rho: object,
    ) -> ParameterDerivativeResult:
        return _EosModel._phase_parameter_derivatives(
            self, active_parameters, T=T, rho=rho, x=(1.0,)
        )

    def fugacity_parameter_derivatives(
        self,
        active_parameters: ActiveParameterSet,
        *,
        T: object,
        P: object,
        phase: Phase,
    ) -> ParameterDerivativeResult:
        return _EosModel._fugacity_parameter_derivatives(
            self, active_parameters, T=T, P=P, x=(1.0,), phase=phase
        )

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
        if not bool(native["accepted"]):
            reason = diagnostics.failure_reason or "local saturation boundary was rejected"
            raise SaturationError(reason, diagnostics)
        return SaturationResult(
            temperature_k=temperature_k,
            saturation_pressure_pa=float(cast(float, native["saturation_pressure_pa"])),
            parameter_fingerprint=self._parameter_fingerprint,
            vapor=_saturation_phase(cast(Mapping[str, object], native["vapor"])),
            liquid=_saturation_phase(cast(Mapping[str, object], native["liquid"])),
            diagnostics=diagnostics,
        )

    def saturation(self, *, T: pint.Quantity[Any]) -> SaturationResult:
        """Solve a local pure-component vapor-liquid saturation boundary."""

        temperature_k = self._saturation_temperature(T)
        try:
            liquid_density_upper = self._native.molar_density_upper_bound(temperature_k, (1.0,))
            native = cast(
                Mapping[str, object],
                _native_core._solve_saturation(
                    self._native,
                    temperature_k,
                    liquid_density_upper,
                ),
            )
        except (RuntimeError, ValueError) as exc:
            failure = _saturation_diagnostics(_failed_saturation_diagnostics(str(exc)))
            raise SaturationError(str(exc), failure) from exc
        return self._saturation_result(native, temperature_k)

    def saturation_parameter_derivatives(
        self,
        active_parameters: ActiveParameterSet,
        *,
        T: pint.Quantity[Any],
    ) -> ParameterDerivativeResult:
        """Solve saturation and return exact implicit selected-parameter derivatives."""

        capability = self.parameter_capability(
            active_parameters, observation_family="pure_saturation"
        )
        temperature_k = self._saturation_temperature(T)
        records = self._active_records(active_parameters)
        try:
            liquid_density_upper = self._native.molar_density_upper_bound(temperature_k, (1.0,))
            native = cast(
                Mapping[str, object],
                _native_core._solve_saturation_active(
                    self._native,
                    temperature_k,
                    liquid_density_upper,
                    records,
                ),
            )
        except (RuntimeError, ValueError) as exc:
            failure = _saturation_diagnostics(_failed_saturation_diagnostics(str(exc)))
            raise SaturationError(str(exc), failure) from exc
        result = self._saturation_result(native, temperature_k)
        parameter_count = len(records)
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
            if not all(math.isfinite(value) for value in mass_values) or not all(
                math.isfinite(value) for row in mass_derivatives for value in row
            ):
                raise _eos_error(
                    "nonfinite-parameter-derivative",
                    observation_family=capability.observation_family,
                    target="mass-density",
                )
            values = (
                *values,
                *mass_values,
            )
            derivatives = (
                *derivatives,
                *mass_derivatives,
            )
        return ParameterDerivativeResult(
            capability,
            values,
            derivatives,
        )

    def parameter_capability(
        self,
        active_parameters: ActiveParameterSet,
        *,
        observation_family: str,
    ) -> ParameterCapability:
        """Return one exact EOS-owned derivative capability without evaluating a state."""

        if observation_family == "fixed_state_pressure":
            return _EosModel.parameter_capability(
                self, active_parameters, observation_family=observation_family
            )
        if observation_family != "pure_saturation":
            raise _eos_error(
                "unsupported-observation-family", observation_family=observation_family
            )
        self._active_records(active_parameters)
        associating = bool(self._model.association_pairs)
        _require_supported_active_parameters(
            active_parameters,
            observation_family,
            (
                (
                    _direct_component_identity(coordinate.identity, coordinate.family)
                    and coordinate.family in {"segment_count", "dispersion_energy_over_k"}
                    or _segment_diameter_parameter_identity(coordinate.identity, coordinate.family)
                )
                or associating
                and _association_identity(coordinate.identity, coordinate.family)
                and coordinate.family in {"association_energy_over_k", "association_volume"}
                for coordinate in active_parameters.coordinates
            ),
        )
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
        return self._capability(
            active_parameters,
            observation_family,
            tuple(labels),
            tuple(units),
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
    ) -> State: ...

    @overload
    def state(
        self,
        *,
        T: pint.Quantity[Any],
        rho: None = None,
        P: pint.Quantity[Any],
        x: Iterable[SupportsFloat],
        phase: Phase | None = None,
    ) -> State: ...

    def state(
        self,
        *,
        T: object,
        rho: object | None = None,
        P: object | None = None,
        x: object,
        phase: str | None = None,
    ) -> State:
        return _EosModel._state(self, T=T, rho=rho, P=P, x=x, phase=phase)

    def phase_parameter_derivatives(
        self,
        active_parameters: ActiveParameterSet,
        *,
        T: object,
        rho: object,
        x: object,
    ) -> ParameterDerivativeResult:
        return _EosModel._phase_parameter_derivatives(self, active_parameters, T=T, rho=rho, x=x)

    def fugacity_parameter_derivatives(
        self,
        active_parameters: ActiveParameterSet,
        *,
        T: object,
        P: object,
        x: object,
        phase: Phase,
    ) -> ParameterDerivativeResult:
        return _EosModel._fugacity_parameter_derivatives(
            self, active_parameters, T=T, P=P, x=x, phase=phase
        )

    def ion_solvation_parameter_derivatives(
        self,
        active_parameters: ActiveParameterSet,
        *,
        reference: ElectrolyteReference,
    ) -> ParameterDerivativeResult:
        """Evaluate selected single-ion solvation energies and exact derivatives."""

        from .reference import ElectrolyteReference as Reference

        records = self._active_records(active_parameters)
        if (
            not isinstance(reference, Reference)
            or reference.component_ids != self._component_ids
            or reference._mixture.parameter_fingerprint != self._parameter_fingerprint
        ):
            raise _eos_error(
                "unsupported-parameter-derivative",
                observation="ion-solvation-gibbs",
                reason="reference-model-mismatch",
            )
        _require_supported_active_parameters(
            active_parameters,
            "ion-solvation-gibbs",
            (
                _direct_component_identity(coordinate.identity, coordinate.family)
                and coordinate.family
                in {
                    "born_diameter",
                    "solvation_factor",
                    "segment_diameter",
                    "dispersion_energy_over_k",
                }
                or coordinate.family == "ionic_region_relative_permittivity"
                and coordinate.identity
                == "model/relative_permittivity/ionic_region_relative_permittivity"
                or coordinate.family == "relative_permittivity"
                and len(coordinate.identity.split("/")) == 5
                and coordinate.identity.split("/")[0] == "component"
                and coordinate.identity.split("/")[1] in self._component_ids
                and coordinate.identity.split("/")[2:] == [
                    "relative_permittivity",
                    "constant",
                    "constant",
                ]
                for coordinate in active_parameters.coordinates
            ),
        )
        try:
            block = self._native.evaluate_ion_solvation_derivatives(
                float(reference.temperature.to("kelvin").magnitude),
                float(reference.pressure.to("pascal").magnitude),
                reference.solvent_mole_fractions,
                reference._ion_indices,
                reference.formula_stoichiometry,
                reference.phase,
                records,
            )
        except (RuntimeError, ValueError) as exc:
            raise _eos_error(
                "unsupported-parameter-derivative",
                observation="ion-solvation-gibbs",
                reason=str(exc),
            ) from exc
        capability = self._capability(
            active_parameters,
            "ion-solvation-gibbs",
            tuple(
                f"solvation_gibbs_j_per_mol/{component}"
                for component in reference.ion_component_ids
            ),
            ("joule / mole", "joule / mole"),
        )
        values = tuple(float(value) for value in block.solvation_gibbs_j_per_mol)
        flat = tuple(float(value) for value in block.parameter_derivatives_j_per_mol)
        parameter_count = len(records)
        return ParameterDerivativeResult(
            capability,
            values,
            tuple(
                flat[row * parameter_count : (row + 1) * parameter_count]
                for row in range(len(values))
            ),
        )

    def mean_ionic_activity_parameter_derivatives(
        self,
        active_parameters: ActiveParameterSet,
        *,
        T: object,
        P: object,
        molality: object,
    ) -> ParameterDerivativeResult:
        """Evaluate one admitted solvation-factor MIAC coordinate exactly."""

        records = self._active_records(active_parameters)
        if len(records) != 1 or active_parameters.coordinates[0].family != "solvation_factor":
            raise _eos_error(
                "unsupported-parameter-derivative",
                observation="mean-ionic-activity",
            )
        coordinate = active_parameters.coordinates[0]
        _require_supported_active_parameters(
            active_parameters,
            "mean-ionic-activity",
            (
                _direct_component_identity(coordinate.identity, coordinate.family)
                and coordinate.identity.split("/")[1] == self._component_ids[0],
            ),
        )
        temperature_k = _state_magnitude(T, "kelvin", "T")
        pressure_pa = _state_magnitude(P, "pascal", "P")
        molality_value = _state_magnitude(molality, "mole / kilogram", "molality")
        try:
            block = self._native.evaluate_aqueous_miac_solvation_factor(
                temperature_k,
                pressure_pa,
                molality_value,
                records[0].value,
            )
        except (RuntimeError, ValueError) as exc:
            raise _eos_error(
                "unsupported-parameter-derivative",
                observation="mean-ionic-activity",
                reason=str(exc),
            ) from exc
        capability = self._capability(
            active_parameters,
            "mean-ionic-activity",
            ("log_mean_ionic_activity_coefficient_molality",),
            ("1",),
        )
        return ParameterDerivativeResult(
            capability,
            (float(block.log_mean_ionic_activity_coefficient_molality),),
            ((float(block.derivative),),),
        )

    def relative_permittivity_parameter_derivatives(
        self,
        active_parameters: ActiveParameterSet,
        *,
        ion_mole_fraction: SupportsFloat,
    ) -> ParameterDerivativeResult:
        """Evaluate ion-fraction permittivity and one exact coefficient derivative."""

        records = self._active_records(active_parameters)
        if (
            len(records) != 1
            or active_parameters.coordinates[0].family != "ion_fraction_suppression_coefficient"
        ):
            raise _eos_error(
                "unsupported-parameter-derivative",
                observation="relative-permittivity",
            )
        _require_supported_active_parameters(
            active_parameters,
            "relative-permittivity",
            (
                active_parameters.coordinates[0].identity
                == "model/relative_permittivity/ion_fraction_suppression_coefficient",
            ),
        )
        fraction = float(ion_mole_fraction)
        try:
            block = self._native.evaluate_ion_fraction_suppression(fraction, records[0].value)
        except (RuntimeError, ValueError) as exc:
            raise _eos_error(
                "unsupported-parameter-derivative",
                observation="relative-permittivity",
                reason=str(exc),
            ) from exc
        capability = self._capability(
            active_parameters,
            "relative-permittivity",
            ("bulk_relative_permittivity",),
            ("1",),
        )
        return ParameterDerivativeResult(
            capability,
            (float(block.bulk_relative_permittivity),),
            ((float(block.derivative),),),
        )

    def reference(
        self,
        *,
        T: object,
        P: object,
        solvent_x: object,
        phase: Literal["liquid"],
        ion_component_ids: tuple[str, str] | None = None,
    ) -> ElectrolyteReference:
        return _EosModel._reference(
            self,
            T=T,
            P=P,
            solvent_x=solvent_x,
            phase=phase,
            ion_component_ids=ion_component_ids,
        )


def _saturation_triple(payload: object) -> tuple[float, float, float]:
    values = tuple(float(value) for value in cast(Iterable[SupportsFloat], payload))
    if len(values) != 3:
        raise RuntimeError("native saturation diagnostic triple has the wrong size")
    return values


def _saturation_attempt(payload: Mapping[str, object]) -> SaturationAttempt:
    return SaturationAttempt(
        role=str(payload["role"]),
        initial_guess=_saturation_triple(payload["initial_guess"]),
        solver_converged=bool(payload["solver_converged"]),
        solver_status=str(payload["solver_status"]),
        iterations=int(cast(int, payload["iterations"])),
        constraint_violation=float(cast(float, payload["constraint_violation"])),
        callback_error=str(payload["callback_error"]),
    )


def _saturation_diagnostics(payload: Mapping[str, object]) -> SaturationDiagnostics:
    attempts = cast(Iterable[Mapping[str, object]], payload["attempt_log"])
    return SaturationDiagnostics(
        solver_converged=bool(payload["solver_converged"]),
        solver_status=str(payload["solver_status"]),
        iterations=int(cast(int, payload["iterations"])),
        attempts=int(cast(int, payload["attempts"])),
        attempt_log=tuple(_saturation_attempt(item) for item in attempts),
        solver_lower_bounds=_saturation_triple(payload["solver_lower_bounds"]),
        solver_upper_bounds=_saturation_triple(payload["solver_upper_bounds"]),
        solver_constraint_violation=float(cast(float, payload["solver_constraint_violation"])),
        numerical_converged=bool(payload["numerical_converged"]),
        confirmation_solves=int(cast(int, payload["confirmation_solves"])),
        confirmation_max_relative_difference=float(
            cast(float, payload["confirmation_max_relative_difference"])
        ),
        physical_accepted=bool(payload["physical_accepted"]),
        pressure_relative_residual=float(cast(float, payload["pressure_relative_residual"])),
        chemical_potential_absolute_residual=float(
            cast(float, payload["chemical_potential_absolute_residual"])
        ),
        phase_density_distance=float(cast(float, payload["phase_density_distance"])),
        exact_derivatives=bool(payload["exact_derivatives"]),
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
        "solver_converged": False,
        "solver_status": "native_exception",
        "iterations": 0,
        "attempts": 0,
        "attempt_log": (),
        "solver_lower_bounds": (math.nan, math.nan, math.nan),
        "solver_upper_bounds": (math.nan, math.nan, math.nan),
        "solver_constraint_violation": math.inf,
        "numerical_converged": False,
        "confirmation_solves": 0,
        "confirmation_max_relative_difference": math.inf,
        "physical_accepted": False,
        "pressure_relative_residual": math.inf,
        "chemical_potential_absolute_residual": math.inf,
        "phase_density_distance": 0.0,
        "exact_derivatives": False,
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
            "polar_formulation": resolved.polar_formulation,
            "dipole_presence": tuple(value > 0.0 for value in resolved.dipole_moment_debye),
            "quadrupole_presence": tuple(
                value > 0.0 for value in resolved.quadrupole_moment_debye_angstrom
            ),
            "helmholtz_terms": (
                "ideal",
                "hard-chain",
                "dispersion",
                *(("polar",) if resolved.polar_formulation is not None else ()),
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
            "polar_formulation": resolved.polar_formulation,
            "dipole_presence": tuple(value > 0.0 for value in resolved.dipole_moment_debye),
            "quadrupole_presence": tuple(
                value > 0.0 for value in resolved.quadrupole_moment_debye_angstrom
            ),
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


def _native_transport_configuration(
    model: _EosModel,
    artifact_fingerprint: str,
) -> dict[str, object]:
    if type(model) not in (Pure, Mixture):
        raise _eos_error("unsupported-native-transport", reason="expected-exact-eos-model")
    fingerprint = model._parameter_fingerprint
    resolved = model._model
    charges = tuple(resolved.charge_number)
    minimum = math.nan if resolved.temperature_min_k is None else resolved.temperature_min_k
    maximum = math.nan if resolved.temperature_max_k is None else resolved.temperature_max_k
    pressure_minimum = math.nan if resolved.pressure_min_pa is None else resolved.pressure_min_pa
    pressure_maximum = math.nan if resolved.pressure_max_pa is None else resolved.pressure_max_pa
    maximum_ion_fraction = (
        math.nan if resolved.ion_mole_fraction_max is None else resolved.ion_mole_fraction_max
    )
    return {
        "parameter_fingerprint": fingerprint,
        "temperature_min_k": minimum,
        "temperature_max_k": maximum,
        "component_ids": model._component_ids,
        "ion_mole_fraction_max": maximum_ion_fraction,
        "pressure_min_pa": pressure_minimum,
        "pressure_max_pa": pressure_maximum,
        "topology_fingerprint": _native_topology_fingerprint(model),
        "component_order_fingerprint": _native_identity_fingerprint(
            "epcsaft-native-component-order-v1",
            {
                "component_ids": model._component_ids,
                "charges": charges,
            },
        ),
        "domain_fingerprint": _native_identity_fingerprint(
            "epcsaft-native-source-domain-v1",
            {
                "temperature_min_k": resolved.temperature_min_k,
                "temperature_max_k": resolved.temperature_max_k,
                "pressure_min_pa": resolved.pressure_min_pa,
                "pressure_max_pa": resolved.pressure_max_pa,
                "ion_mole_fraction_max": resolved.ion_mole_fraction_max,
            },
        ),
        "artifact_fingerprint": artifact_fingerprint,
    }


def _native_transport(model: _EosModel) -> object:
    """Return the model-bound capsule used only by Engine's native modules."""
    configuration = _native_transport_configuration(model, _python_native_artifact_fingerprint())
    try:
        return _native_transport_factory(model._native, **configuration)
    except (RuntimeError, ValueError) as exc:
        raise _eos_error("unsupported-native-transport", reason=str(exc)) from exc


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
    hard_chain: float
    dispersion: float
    association: float
    ares_polar: float
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
    pressure_density_derivative: pint.Quantity[float]
    iterations: int
    stable: bool
    branch: str


__all__ = (
    "DensityDiagnostics",
    "EosError",
    "Fugacity",
    "Mixture",
    "Phase",
    "PropertyBasis",
    "Pure",
    "SaturationAttempt",
    "SaturationDiagnostics",
    "SaturationError",
    "SaturationPhase",
    "SaturationResult",
    "State",
)
