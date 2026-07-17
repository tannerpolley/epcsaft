from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Iterable

import pint

from ._core import NativeModel as _NativeModel
from ._resolve import _ResolvedModel, _resolve_model
from .parameters import ParameterSet
from .records import SingleParameterRecord
from .units import ParameterError, _error, _quantity, unit_registry


def _magnitude(record: SingleParameterRecord, unit: str) -> float:
    value = record.value
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


def _native_arguments(model: _ResolvedModel) -> tuple[list[float], ...]:
    if model.topology != "neutral" or not all(
        isinstance(record, SingleParameterRecord) for record in model.segment_diameter
    ):
        raise _error("invalid-value", topology=model.topology, reason="native-not-admitted")
    return (
        list(model.segment_count),
        [_magnitude(record, "angstrom") for record in model.segment_diameter],
        list(model.dispersion_energy_over_k),
        list(model.k_ij),
    )


@dataclass(frozen=True, slots=True, init=False)
class EPCSAFT:
    _component_ids: tuple[str, ...]
    _parameter_fingerprint: str
    _model: _ResolvedModel
    _native: _NativeModel

    def __init__(self, parameters: ParameterSet) -> None:
        model = _resolve_model(parameters)
        m, sigma, epsilon_over_k, kij = _native_arguments(model)
        object.__setattr__(self, "_component_ids", parameters.component_ids)
        object.__setattr__(self, "_parameter_fingerprint", parameters.fingerprint)
        object.__setattr__(self, "_model", model)
        object.__setattr__(self, "_native", _NativeModel(m, sigma, epsilon_over_k, kij))

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

        try:
            native = self._native.evaluate(temperature_k, density_molar, fractions)
        except (RuntimeError, ValueError) as exc:
            raise _eos_error("evaluation-failed", reason=str(exc)) from exc
        return EosResult(
            hard_chain=native.hard_chain,
            dispersion=native.dispersion,
            residual_helmholtz=native.residual_helmholtz,
            compressibility_factor=native.compressibility_factor,
            pressure=native.pressure * unit_registry.pascal,
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
    residual_helmholtz: float
    compressibility_factor: float
    pressure: pint.Quantity


__all__ = ("EPCSAFT", "EosError", "EosResult")
