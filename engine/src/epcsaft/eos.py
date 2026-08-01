from __future__ import annotations

import hashlib
import json
import math
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from functools import lru_cache
from os import PathLike
from pathlib import Path
from typing import TYPE_CHECKING, Any, Literal, SupportsFloat, TypeAlias, overload

import pint

from . import _core as _native_core
from ._resolve import _ResolvedModel, _resolve_model
from .parameters import Parameters
from .records import (
    ConstantCorrelation,
    ConstantPlusSumOfExponentialsCorrelation,
    SingleParameterRecord,
)
from .units import ParameterError, _error, _quantity, unit_registry

if TYPE_CHECKING:
    from .reference import Activity, ElectrolyteReference


Phase: TypeAlias = Literal["vapor", "liquid"]
PropertyBasis: TypeAlias = Literal["molar", "reduced"]
_MOLAR_GAS_CONSTANT = 8.31446261815324
_NativeModel = _native_core.NativeModel
_native_sdk_v1 = _native_core.native_sdk_v1


@lru_cache(maxsize=1)
def _python_native_artifact_fingerprint() -> str:
    return _file_fingerprint(Path(_native_core.__file__))


@lru_cache(maxsize=1)
def _static_native_artifact_fingerprint() -> str:
    return _file_fingerprint(Path(_native_core.__file__).parent / "lib" / "libepcsaft_native_sdk.a")


def _file_fingerprint(path: Path) -> str:
    return f"sha256:{hashlib.sha256(path.read_bytes()).hexdigest()}"


def _magnitude(
    record: SingleParameterRecord | ConstantCorrelation,
    unit: str,
) -> float:
    value = record.value if isinstance(record, SingleParameterRecord) else record.constant
    if unit != "dimensionless":
        value = value.to(unit).magnitude
    return float(value)


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
            sigma_constant.append(float(record.constant.to("angstrom").magnitude))
            if isinstance(record, ConstantPlusSumOfExponentialsCorrelation):
                for term in record.terms:
                    term_component.append(component)
                    term_amplitude.append(float(term.amplitude.to("angstrom").magnitude))
                    term_exponent.append(
                        float(term.exponent_coefficient.to("1 / kelvin").magnitude)
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
            solvation_factor = [0.0] * components
            ionic_permittivity = 0.0
            electrolyte_model = 3
        else:
            dielectric_suppression = 0.0
            born_diameter = [0.0] * components
            solvation_factor = [0.0] * components
            ionic_permittivity = 0.0
            electrolyte_model = 2
    else:
        molar_mass = [0.0] * components
        charge_number = [0] * components
        relative_permittivity = [0.0] * components
        dielectric_suppression = 0.0
        born_diameter = [0.0] * components
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
    )


@dataclass(frozen=True, slots=True, init=False)
class Mixture:
    _component_ids: tuple[str, ...]
    _parameter_fingerprint: str
    _model: _ResolvedModel
    _native: _NativeModel

    def __init__(self, parameters: Parameters) -> None:
        model = _resolve_model(parameters)
        native_arguments = _native_arguments(model)
        object.__setattr__(self, "_component_ids", parameters.component_ids)
        object.__setattr__(self, "_parameter_fingerprint", parameters.fingerprint)
        object.__setattr__(self, "_model", model)
        object.__setattr__(self, "_native", _NativeModel(*native_arguments))

    @classmethod
    def from_dictionary(cls, mapping: Mapping[str, object]) -> Mixture:
        """Construct a mixture from the strict fixed-unit parameter dictionary."""
        return cls(Parameters.from_dictionary(mapping))

    @property
    def component_ids(self) -> tuple[str, ...]:
        return self._component_ids

    @property
    def relative_permittivity_formulation(self) -> str:
        """Return the generic relative-permittivity formulation used by this state."""
        return self._model.relative_permittivity_formulation

    @property
    def parameter_fingerprint(self) -> str:
        return self._parameter_fingerprint

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
            fractions = tuple(float(value) for value in x)  # type: ignore[union-attr]
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
                raise _eos_error("evaluation-failed", reason=str(exc)) from exc
        else:
            pressure_pa = _state_magnitude(P, "pascal", "P")
            try:
                native = self._native.evaluate_pressure(
                    temperature_k, pressure_pa, fractions, phase or ""
                )
            except (RuntimeError, ValueError) as exc:
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

    def reference(
        self,
        *,
        T: object,
        P: object,
        solvent_x: object,
        phase: Literal["liquid"],
    ) -> ElectrolyteReference:
        """Build one liquid electrolyte reference on the selected component order.

        The model must contain one source-complete cation/anion formula unit and
        solvent-only normalized reference fractions. Failures raise ``EosError``.
        """
        from .reference import build_reference_state

        return build_reference_state(
            self,
            temperature=T,
            pressure=P,
            solvent_mole_fractions=solvent_x,
            phase=phase,
        )


def _native_topology_fingerprint(model: Mixture) -> str:
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
            "helmholtz_terms": ("ideal", "hard-chain", "dispersion"),
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


def _native_sdk_configuration(
    model: Mixture,
    artifact_fingerprint: str,
) -> dict[str, object]:
    if type(model) is not Mixture:
        raise _eos_error("unsupported-native-sdk", reason="expected-exact-Mixture-model")
    fingerprint = model._parameter_fingerprint
    resolved = model._model
    charges = tuple(resolved.charge_number)
    corrected_triplet = (
        resolved.topology == "ionic-ssm-ds"
        and len(model._component_ids) == 3
        and charges[0] == 0
        and charges[1] != 0
        and charges[2] == -charges[1]
    )
    ordered_electrolyte_triplet = corrected_triplet and charges == (0, 1, -1)
    fixed_reference_conditions = (
        resolved.temperature_min_k == resolved.temperature_max_k == 298.15
        and resolved.pressure_min_pa == resolved.pressure_max_pa == 100000.0
    )
    complete_reference_domain = all(
        value is not None
        for value in (
            resolved.temperature_min_k,
            resolved.temperature_max_k,
            resolved.pressure_min_pa,
            resolved.pressure_max_pa,
        )
    )
    phase_admitted = True
    active_born_admitted = corrected_triplet and fixed_reference_conditions
    aqueous_kij_admitted = ordered_electrolyte_triplet and fixed_reference_conditions
    ion_fraction_suppression_admitted = (
        ordered_electrolyte_triplet
        and fixed_reference_conditions
        and resolved.relative_permittivity_formulation == "ion-fraction-suppression"
    )
    ion_solvation_kij_admitted = ordered_electrolyte_triplet and fixed_reference_conditions
    aqueous_solvation_factor_admitted = (
        ordered_electrolyte_triplet
        and fixed_reference_conditions
        and resolved.ion_mole_fraction_max is not None
    )
    neutral_reference_admitted = (
        complete_reference_domain
        and len(charges) == len(model._component_ids)
        and sum(charge == 0 for charge in charges) == 1
        and any(charge > 0 for charge in charges)
        and any(charge < 0 for charge in charges)
    )
    minimum = math.nan if resolved.temperature_min_k is None else resolved.temperature_min_k
    maximum = math.nan if resolved.temperature_max_k is None else resolved.temperature_max_k
    pressure_minimum = math.nan if resolved.pressure_min_pa is None else resolved.pressure_min_pa
    pressure_maximum = math.nan if resolved.pressure_max_pa is None else resolved.pressure_max_pa
    maximum_ion_fraction = (
        math.nan if resolved.ion_mole_fraction_max is None else resolved.ion_mole_fraction_max
    )
    unknown_applicability = all(
        value is None
        for value in (
            resolved.temperature_min_k,
            resolved.temperature_max_k,
            resolved.pressure_min_pa,
            resolved.pressure_max_pa,
            resolved.ion_mole_fraction_max,
        )
    )
    return {
        "parameter_fingerprint": fingerprint,
        "temperature_min_k": minimum,
        "temperature_max_k": maximum,
        "component_ids": model._component_ids,
        "site_ids": tuple(site_id for _, site_id in model._model.site_ids),
        "ion_mole_fraction_max": maximum_ion_fraction,
        "phase_callbacks_admitted": phase_admitted,
        "active_born_admitted": active_born_admitted,
        "aqueous_kij_admitted": aqueous_kij_admitted,
        "ion_fraction_suppression_admitted": ion_fraction_suppression_admitted,
        "ion_solvation_kij_admitted": ion_solvation_kij_admitted,
        "pressure_min_pa": pressure_minimum,
        "pressure_max_pa": pressure_maximum,
        "neutral_reference_admitted": neutral_reference_admitted,
        "unknown_applicability_admitted": unknown_applicability,
        "aqueous_solvation_factor_admitted": aqueous_solvation_factor_admitted,
        "topology_fingerprint": _native_topology_fingerprint(model),
        "artifact_fingerprint": artifact_fingerprint,
    }


def native_sdk(model: Mixture) -> object:
    """Return the model-bound version-1 native provider capsule."""
    configuration = _native_sdk_configuration(model, _python_native_artifact_fingerprint())
    try:
        return _native_sdk_v1(model._native, **configuration)
    except (RuntimeError, ValueError) as exc:
        raise _eos_error("unsupported-native-sdk", reason=str(exc)) from exc


def export_native_model(model: Mixture, path: str | PathLike[str]) -> None:
    """Write one immutable resolved-model artifact for the installed native SDK."""
    configuration = _native_sdk_configuration(model, _static_native_artifact_fingerprint())
    names = (
        "segment_count",
        "sigma_constant_angstrom",
        "sigma_term_component",
        "sigma_term_amplitude_angstrom",
        "sigma_term_exponent_per_k",
        "dispersion_energy_over_k_k",
        "k_ij",
        "l_ij",
        "site_component",
        "site_multiplicity",
        "association_active",
        "association_energy_over_k_k",
        "association_volume",
        "molar_mass_kg_per_mol",
        "charge_number",
        "relative_permittivity",
        "ion_fraction_suppression_coefficient",
        "born_diameter_angstrom",
        "solvation_factor",
        "ionic_region_relative_permittivity",
        "electrolyte_model",
        "dielectric_volume_m3_per_mol",
        "dielectric_temperature_k",
        "ion_specific_suppression_coefficient",
        "dipole_scaling_debye2",
        "polarizability_scaling_angstrom3",
        "correlation_integral_parameter",
        "relative_permittivity_formulation",
    )
    native_parameters = dict(zip(names, _native_arguments(model._model), strict=True))
    for key in (
        "temperature_min_k",
        "temperature_max_k",
        "ion_mole_fraction_max",
        "pressure_min_pa",
        "pressure_max_pa",
    ):
        if math.isnan(float(configuration[key])):
            configuration[key] = None
    payload = {
        "schema": "epcsaft-resolved-native-model-v1",
        "native_parameters": native_parameters,
        "sdk": configuration,
    }
    target = Path(path)
    with target.open("x", encoding="utf-8") as stream:
        json.dump(payload, stream, allow_nan=False, separators=(",", ":"), sort_keys=True)
        stream.write("\n")


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
        return (
            _MOLAR_GAS_CONSTANT
            * self.temperature.to("kelvin").magnitude
            * unit_registry.joule
            / unit_registry.mole
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
    "State",
    "export_native_model",
    "native_sdk",
)
