from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Iterable

import pint

from ._core import NativeModel as _NativeModel
from ._resolve import _ResolvedModel, _resolve_model
from .parameters import ParameterSet
from .records import (
    ConstantCorrelation,
    ConstantPlusSumOfExponentialsCorrelation,
    SingleParameterRecord,
)
from .units import ParameterError, _error, _quantity, unit_registry


def _magnitude(
    record: SingleParameterRecord | ConstantCorrelation,
    unit: str,
) -> float:
    value = (
        record.value
        if isinstance(record, SingleParameterRecord)
        else record.constant
    )
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
    if model.topology == "ionic-ssm-ds":
        molar_mass = [
            0.0 if value is None else value
            for value in model.molar_mass_kg_per_mol
        ]
        charge_number = list(model.charge_number)
        relative_permittivity = [
            0.0 if record is None else _magnitude(record, "dimensionless")
            for record in model.relative_permittivity
        ]
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
        assert dielectric_suppression is not None
        assert ionic_permittivity is not None
        assert len(solvation_factor) == components
    else:
        molar_mass = [0.0] * components
        charge_number = [0] * components
        relative_permittivity = [0.0] * components
        dielectric_suppression = 0.0
        born_diameter = [0.0] * components
        solvation_factor = [0.0] * components
        ionic_permittivity = 0.0

    return (
        list(model.segment_count),
        sigma_constant,
        term_component,
        term_amplitude,
        term_exponent,
        list(model.dispersion_energy_over_k),
        list(model.k_ij),
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
    )


@dataclass(frozen=True, slots=True, init=False)
class EPCSAFT:
    _component_ids: tuple[str, ...]
    _parameter_fingerprint: str
    _model: _ResolvedModel
    _native: _NativeModel

    def __init__(self, parameters: ParameterSet) -> None:
        model = _resolve_model(parameters)
        native_arguments = _native_arguments(model)
        object.__setattr__(self, "_component_ids", parameters.component_ids)
        object.__setattr__(self, "_parameter_fingerprint", parameters.fingerprint)
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
        molar_density: object,
        mole_fractions: object,
    ) -> EosResult:
        temperature_k = _state_magnitude(temperature, "kelvin", "temperature")
        density_molar = _state_magnitude(
            molar_density, "mole / meter ** 3", "molar_density"
        )
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

        if self._model.topology == "ionic-ssm-ds":
            if (
                self._model.temperature_min_k is not None
                and temperature_k < self._model.temperature_min_k
            ) or (
                self._model.temperature_max_k is not None
                and temperature_k > self._model.temperature_max_k
            ):
                raise _eos_error("invalid-state", field="temperature")
            ion_fraction = math.fsum(
                fraction
                for fraction, charge in zip(
                    fractions, self._model.charge_number, strict=True
                )
                if charge != 0
            )
            if (
                self._model.ion_mole_fraction_max is not None
                and ion_fraction > self._model.ion_mole_fraction_max + 1.0e-12
            ):
                raise _eos_error("invalid-state", field="mole_fractions")
            charge_balance = math.fsum(
                fraction * charge
                for fraction, charge in zip(
                    fractions, self._model.charge_number, strict=True
                )
            )
            if abs(charge_balance) > 1.0e-12:
                raise _eos_error("invalid-state", field="electroneutrality")

        try:
            native = self._native.evaluate(temperature_k, density_molar, fractions)
        except (RuntimeError, ValueError) as exc:
            raise _eos_error("evaluation-failed", reason=str(exc)) from exc
        return EosResult(
            hard_chain=native.hard_chain,
            dispersion=native.dispersion,
            association=native.association,
            debye_huckel=native.debye_huckel,
            born_ssm_ds=native.born_ssm_ds,
            residual_helmholtz=native.residual_helmholtz,
            compressibility_factor=native.compressibility_factor,
            pressure=native.pressure * unit_registry.pascal,
            residual_chemical_potential_over_rt=tuple(
                native.residual_chemical_potential_over_rt
            ),
            log_fugacity_coefficient=None
            if native.log_fugacity_coefficient is None
            else tuple(native.log_fugacity_coefficient),
            residual_entropy=native.residual_entropy
            * unit_registry.joule / unit_registry.mole / unit_registry.kelvin,
            residual_internal_energy=native.residual_internal_energy
            * unit_registry.joule / unit_registry.mole,
            residual_enthalpy=native.residual_enthalpy
            * unit_registry.joule / unit_registry.mole,
            residual_isochoric_heat_capacity=
            native.residual_isochoric_heat_capacity
            * unit_registry.joule / unit_registry.mole / unit_registry.kelvin,
            residual_isobaric_heat_capacity=None
            if native.residual_isobaric_heat_capacity is None
            else native.residual_isobaric_heat_capacity
            * unit_registry.joule / unit_registry.mole / unit_registry.kelvin,
        )


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


__all__ = ("EPCSAFT", "EosError", "EosResult")
