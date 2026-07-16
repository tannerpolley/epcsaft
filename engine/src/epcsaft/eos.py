from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
import math
from typing import Final

import pint

from ._core import NativeModel as _NativeModel
from .parameters import ParameterSet
from .records import PairParameterRecord, SingleParameterRecord
from .units import _error, _quantity, unit_registry


_REQUIRED_SINGLE: Final = {
    "segment_count": "dimensionless",
    "segment_diameter": "angstrom",
    "dispersion_energy_over_k": "kelvin",
}
_ALLOWED_SINGLE = frozenset({"molar_mass", *_REQUIRED_SINGLE})
_ALLOWED_PAIR = frozenset({"k_ij"})


def _reject(record: object) -> None:
    raise _error(
        "invalid-value",
        record_id=getattr(record, "record_id", type(record).__name__),
        family=getattr(record, "family", "site"),
    )


def _magnitude(record: SingleParameterRecord, unit: str) -> float:
    value = record.value
    if unit != "dimensionless":
        value = value.to(unit).magnitude
    return float(value)


def _resolve(parameters: ParameterSet) -> tuple[list[float], ...]:
    if not isinstance(parameters, ParameterSet):
        raise _error("invalid-value", field="parameters", expected="ParameterSet")

    for record in parameters.records:
        if isinstance(record, SingleParameterRecord):
            if record.family not in _ALLOWED_SINGLE:
                _reject(record)
        elif isinstance(record, PairParameterRecord):
            if record.family not in _ALLOWED_PAIR:
                _reject(record)
        else:
            _reject(record)

    singles = {
        (record.component_id, record.family): record
        for record in parameters.records
        if isinstance(record, SingleParameterRecord)
    }
    resolved: list[list[float]] = []
    for family, unit in _REQUIRED_SINGLE.items():
        values: list[float] = []
        for component_id in parameters.component_ids:
            record = singles.get((component_id, family))
            if record is None:
                raise _error(
                    "missing-record", component_id=component_id, family=family
                )
            values.append(_magnitude(record, unit))
        resolved.append(values)

    size = len(parameters.component_ids)
    kij = [0.0] * (size * size)
    pairs = {
        (record.component_id_a, record.component_id_b): record
        for record in parameters.records
        if isinstance(record, PairParameterRecord)
    }
    positions = {
        component_id: index
        for index, component_id in enumerate(parameters.component_ids)
    }
    for component_a, component_b in combinations(parameters.component_ids, 2):
        key = tuple(sorted((component_a, component_b)))
        record = pairs.get(key)
        if record is None:
            raise _error(
                "missing-record",
                component_id_a=key[0],
                component_id_b=key[1],
                family="k_ij",
            )
        i, j = positions[component_a], positions[component_b]
        kij[i * size + j] = kij[j * size + i] = float(record.value)
    return (*resolved, kij)


@dataclass(frozen=True, slots=True, init=False)
class EPCSAFT:
    _component_ids: tuple[str, ...]
    _parameter_fingerprint: str
    _native: _NativeModel

    def __init__(self, parameters: ParameterSet) -> None:
        m, sigma, epsilon_over_k, kij = _resolve(parameters)
        object.__setattr__(self, "_component_ids", parameters.component_ids)
        object.__setattr__(self, "_parameter_fingerprint", parameters.fingerprint)
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
            raise ValueError("invalid mole_fractions") from exc
        if (
            len(fractions) != len(self._component_ids)
            or any(not math.isfinite(value) or value < 0.0 for value in fractions)
            or abs(math.fsum(fractions) - 1.0) > 1.0e-12
        ):
            raise ValueError("invalid mole_fractions")

        native = self._native.evaluate(temperature_k, density_molar, fractions)
        return EosResult(
            hard_chain=native.hard_chain,
            dispersion=native.dispersion,
            residual_helmholtz=native.residual_helmholtz,
            compressibility_factor=native.compressibility_factor,
            pressure=native.pressure * unit_registry.pascal,
        )


def _state_magnitude(value: object, unit: str, field: str) -> float:
    quantity = _quantity(value, unit, field=field)
    magnitude = float(quantity.magnitude)
    if not math.isfinite(magnitude) or magnitude <= 0.0:
        raise ValueError(f"{field} must be positive and finite")
    return magnitude


@dataclass(frozen=True, slots=True)
class EosResult:
    hard_chain: float
    dispersion: float
    residual_helmholtz: float
    compressibility_factor: float
    pressure: pint.Quantity


__all__ = ("EPCSAFT", "EosResult")
