from __future__ import annotations

from dataclasses import dataclass
import math
from typing import TYPE_CHECKING, Iterable

import pint

from ._core import NativeModel as _NativeModel
from ._core import native_sdk_v1 as _native_sdk_v1
from ._resolve import _ResolvedModel, _resolve_model
from .parameters import ParameterSet
from .records import (
    ConstantCorrelation,
    ConstantPlusSumOfExponentialsCorrelation,
    SingleParameterRecord,
)
from .units import ParameterError, _error, _quantity, unit_registry

if TYPE_CHECKING:
    from .reference import ElectrolyteReference


_NATIVE_SDK_V1_PARAMETER_FINGERPRINTS = frozenset(
    {
        "sha256:5f836aa84935df70be2e5cffae51b178a7b797c2cee036e9ff47d8097ca94bbf",
        "sha256:288fbcaa1304881c16f64c3a784eeed19b75c58cca4558f92a21268e5e91258a",
        "sha256:9bfbc8d7789e51609945e61dbdf7a020decc8f9e31b408b0977724c7cb3e1551",
        "sha256:307fcb28d535b94782f3e3caf4012c0c8c0dc87ee4239d6c316de56553543286",
        "sha256:cb0c63530f96d51d5a24f48a2d02ae6567d163abdc373099ce096f7ff0a896f5",
        "sha256:7c637771bc9f717b8f47b44bb2a61044c3fe83084dca7c3c16102fba0989912d",
        "sha256:7551f1eee5903b66061cf7520f3bb59b169896ce372f3df3d48aa7ec778c39d4",
    }
)

_NATIVE_SDK_V1_ACTIVE_BORN_FINGERPRINTS = frozenset(
    {
        "sha256:1bb528ebe8f5612757e148608fc55821f9fb03737dbcec6d0bc4fffd0f4cbc4c",
        "sha256:7c637771bc9f717b8f47b44bb2a61044c3fe83084dca7c3c16102fba0989912d",
        "sha256:d29cef0c0f63034436d547d0aafa57934effe06783c8dffd89c94fa85e6940f4",
        "sha256:7551f1eee5903b66061cf7520f3bb59b169896ce372f3df3d48aa7ec778c39d4",
        "sha256:70ae04599dfa8338175e793bac6b9e4dfad37a9b96a568b5484dc87f104ef1a9",
    }
)

_NATIVE_SDK_V1_PARAMETER_BUNDLE_SELECTIONS = {
    "sha256:8ea2d7acbc179b82348fa938d081fe8d56f274bfe1f2167d10a507b7711aa1c2": frozenset(
        {
            "water",
            "hydronium-cation",
            "hydroxide-anion",
        }
    ),
    "sha256:03f9b1d88843a2d53da3420d96987ec4a2305b2814cc4946cdb25d1db67f5f72": frozenset(
        {
            "water",
            "ethanol",
            "isobutanol",
            "sodium-cation",
            "chloride-anion",
        }
    ),
    "sha256:b7a54dab41c61a71803f1d63bf25396188a2a1edc55172d5beddb99c3368c9fd": frozenset(
        {
            "water",
            "1-butanol",
            "sodium-cation",
            "potassium-cation",
            "chloride-anion",
        }
    ),
}


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
    if model.topology in {
        "ionic-held-2008-debye-huckel",
        "ionic-ssm-ds",
        "ionic-ascani-original-born",
    }:
        molar_mass = [0.0 if value is None else value for value in model.molar_mass_kg_per_mol]
        charge_number = list(model.charge_number)
        relative_permittivity = [
            0.0 if record is None else _magnitude(record, "dimensionless")
            for record in model.relative_permittivity
        ]
        if model.topology == "ionic-ssm-ds":
            dielectric_suppression = model.dielectric_ion_suppression_coefficient
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
            assert dielectric_suppression is not None
            assert ionic_permittivity is not None
            assert len(solvation_factor) == components
        elif model.topology == "ionic-held-2008-debye-huckel":
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
    )


@dataclass(frozen=True, slots=True, init=False)
class EPCSAFT:
    _component_ids: tuple[str, ...]
    _parameter_fingerprint: str
    _parameter_bundle_fingerprint: str
    _model: _ResolvedModel
    _native: _NativeModel

    def __init__(self, parameters: ParameterSet) -> None:
        model = _resolve_model(parameters)
        native_arguments = _native_arguments(model)
        object.__setattr__(self, "_component_ids", parameters.component_ids)
        object.__setattr__(self, "_parameter_fingerprint", parameters.fingerprint)
        object.__setattr__(self, "_parameter_bundle_fingerprint", parameters.bundle_fingerprint)
        object.__setattr__(self, "_model", model)
        object.__setattr__(self, "_native", _NativeModel(*native_arguments))

    @property
    def component_ids(self) -> tuple[str, ...]:
        return self._component_ids

    @property
    def parameter_fingerprint(self) -> str:
        return self._parameter_fingerprint

    def evaluate(
        self,
        *,
        temperature: object,
        molar_density: object | None = None,
        pressure: object | None = None,
        mole_fractions: object,
        phase: str | None = None,
    ) -> EosResult:
        if (molar_density is None) == (pressure is None):
            raise _eos_error("invalid-state", reason="specify-exactly-one-of-density-or-pressure")
        if phase not in (None, "vapor", "liquid"):
            raise _eos_error("invalid-state", field="phase")
        if molar_density is not None and phase is not None:
            raise _eos_error("invalid-state", field="phase")
        temperature_k = _state_magnitude(temperature, "kelvin", "temperature")
        try:
            fractions = tuple(float(value) for value in mole_fractions)  # type: ignore[union-attr]
        except (TypeError, ValueError) as exc:
            raise _eos_error("invalid-state", field="mole_fractions") from exc
        if (
            len(fractions) != len(self._component_ids)
            or any(not math.isfinite(value) or value < 0.0 for value in fractions)
            or abs(math.fsum(fractions) - 1.0) > 1.0e-12
        ):
            raise _eos_error("invalid-state", field="mole_fractions")

        if (
            self._model.temperature_min_k is not None
            and temperature_k < self._model.temperature_min_k
        ) or (
            self._model.temperature_max_k is not None
            and temperature_k > self._model.temperature_max_k
        ):
            raise _eos_error("invalid-state", field="temperature")

        if self._model.topology in {"ionic-ssm-ds", "ionic-ascani-original-born"}:
            ion_fraction = math.fsum(
                fraction
                for fraction, charge in zip(fractions, self._model.charge_number, strict=True)
                if charge != 0
            )
            if (
                self._model.ion_mole_fraction_max is not None
                and ion_fraction > self._model.ion_mole_fraction_max + 1.0e-12
            ):
                raise _eos_error("invalid-state", field="mole_fractions")
            charge_balance = math.fsum(
                fraction * charge
                for fraction, charge in zip(fractions, self._model.charge_number, strict=True)
            )
            if abs(charge_balance) > 1.0e-12:
                raise _eos_error("invalid-state", field="electroneutrality")

        try:
            if molar_density is not None:
                density_molar = _state_magnitude(
                    molar_density, "mole / meter ** 3", "molar_density"
                )
                native = self._native.evaluate(temperature_k, density_molar, fractions)
            else:
                pressure_pa = _state_magnitude(pressure, "pascal", "pressure")
                native = self._native.evaluate_pressure(
                    temperature_k, pressure_pa, fractions, phase or ""
                )
        except (RuntimeError, ValueError) as exc:
            raise _eos_error("evaluation-failed", reason=str(exc)) from exc
        diagnostics = native.density_diagnostics
        return EosResult(
            hard_chain=native.hard_chain,
            dispersion=native.dispersion,
            association=native.association,
            debye_huckel=native.debye_huckel,
            born_ssm_ds=native.born_ssm_ds,
            residual_helmholtz=native.residual_helmholtz,
            compressibility_factor=native.compressibility_factor,
            pressure=native.pressure * unit_registry.pascal,
            residual_chemical_potential_over_rt=tuple(native.residual_chemical_potential_over_rt),
            log_fugacity_coefficient=None
            if native.log_fugacity_coefficient is None
            else tuple(native.log_fugacity_coefficient),
            residual_entropy=native.residual_entropy
            * unit_registry.joule
            / unit_registry.mole
            / unit_registry.kelvin,
            residual_internal_energy=native.residual_internal_energy
            * unit_registry.joule
            / unit_registry.mole,
            residual_enthalpy=native.residual_enthalpy * unit_registry.joule / unit_registry.mole,
            residual_isochoric_heat_capacity=native.residual_isochoric_heat_capacity
            * unit_registry.joule
            / unit_registry.mole
            / unit_registry.kelvin,
            residual_isobaric_heat_capacity=None
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

    def reference_state(
        self,
        *,
        temperature: object,
        pressure: object,
        solvent_mole_fractions: object,
        phase: str,
    ) -> ElectrolyteReference:
        from .reference import build_reference_state

        return build_reference_state(
            self,
            temperature=temperature,
            pressure=pressure,
            solvent_mole_fractions=solvent_mole_fractions,
            phase=phase,
        )


def native_sdk(model: EPCSAFT) -> object:
    """Return the model-bound version-1 native provider capsule."""
    if type(model) is not EPCSAFT:
        raise _eos_error("unsupported-native-sdk", reason="expected-exact-EPCSAFT-model")
    fingerprint = model._parameter_fingerprint
    bundle_selection = _NATIVE_SDK_V1_PARAMETER_BUNDLE_SELECTIONS.get(
        model._parameter_bundle_fingerprint
    )
    phase_admitted = fingerprint in _NATIVE_SDK_V1_PARAMETER_FINGERPRINTS or (
        bundle_selection is not None
        and len(model._component_ids) == len(bundle_selection)
        and frozenset(model._component_ids) == bundle_selection
    )
    active_born_admitted = fingerprint in _NATIVE_SDK_V1_ACTIVE_BORN_FINGERPRINTS
    if not phase_admitted and not active_born_admitted:
        raise _eos_error(
            "unsupported-native-sdk",
            reason="unadmitted-parameter-fingerprint",
        )
    minimum = math.nan if model._model.temperature_min_k is None else model._model.temperature_min_k
    maximum = math.nan if model._model.temperature_max_k is None else model._model.temperature_max_k
    minimum_pressure = (
        math.nan if model._model.pressure_min_pa is None else model._model.pressure_min_pa
    )
    maximum_pressure = (
        math.nan if model._model.pressure_max_pa is None else model._model.pressure_max_pa
    )
    maximum_ion_fraction = (
        math.nan
        if model._model.ion_mole_fraction_max is None
        else model._model.ion_mole_fraction_max
    )
    try:
        return _native_sdk_v1(
            model._native,
            fingerprint,
            minimum,
            maximum,
            minimum_pressure,
            maximum_pressure,
            model._component_ids,
            maximum_ion_fraction,
            phase_admitted,
            active_born_admitted,
        )
    except (RuntimeError, ValueError) as exc:
        raise _eos_error("unsupported-native-sdk", reason=str(exc)) from exc


def _state_magnitude(value: object, unit: str, field: str) -> float:
    try:
        quantity = _quantity(value, unit, field=field)
    except ParameterError as exc:
        raise _eos_error("invalid-state", field=field) from exc
    magnitude = float(quantity.magnitude)
    if not math.isfinite(magnitude) or magnitude <= 0.0:
        raise _eos_error("invalid-state", field=field)
    return magnitude


@dataclass(frozen=True, slots=True)
class EosResult:
    hard_chain: float
    dispersion: float
    association: float
    debye_huckel: float
    born_ssm_ds: float
    residual_helmholtz: float
    compressibility_factor: float
    pressure: pint.Quantity
    residual_chemical_potential_over_rt: tuple[float, ...]
    log_fugacity_coefficient: tuple[float, ...] | None
    residual_entropy: pint.Quantity
    residual_internal_energy: pint.Quantity
    residual_enthalpy: pint.Quantity
    residual_isochoric_heat_capacity: pint.Quantity
    residual_isobaric_heat_capacity: pint.Quantity | None
    molar_density: pint.Quantity
    density_diagnostics: DensityDiagnostics | None


@dataclass(frozen=True, slots=True)
class DensityDiagnostics:
    bracket_lower: pint.Quantity
    bracket_upper: pint.Quantity
    pressure_residual: pint.Quantity
    pressure_density_derivative: pint.Quantity
    iterations: int
    stable: bool
    branch: str


__all__ = ("EPCSAFT", "EosError", "EosResult", "DensityDiagnostics", "native_sdk")
