from __future__ import annotations

from collections.abc import Iterable
from decimal import Decimal, InvalidOperation
from typing import Any

import pint

unit_registry: pint.UnitRegistry[Any] = pint.UnitRegistry()


class ParameterError(ValueError):
    """One stable public error for invalid parameter input."""

    __slots__ = ("_code", "_context")

    def __init__(
        self,
        code: str,
        context: Iterable[tuple[str, str]] = (),
    ) -> None:
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


def _error(code: str, **context: object) -> ParameterError:
    return ParameterError(code, tuple((key, str(value)) for key, value in context.items()))


def _decimal_number(value: object, *, field: str) -> Decimal:
    if isinstance(value, bool) or not isinstance(value, (int, float, Decimal)):
        raise _error("invalid-value", field=field, value_type=type(value).__name__)
    try:
        result = Decimal(str(value))
    except (InvalidOperation, ValueError) as exc:
        raise _error("invalid-value", field=field, value=value) from exc
    if not result.is_finite():
        raise _error("invalid-value", field=field, value=value)
    return Decimal(0) if result.is_zero() else result


def _exact_integer(value: object, *, field: str, positive: bool = False) -> int:
    if type(value) is not int or (positive and value <= 0):
        raise _error("invalid-value", field=field, value=value)
    return value


def _quantity(value: object, target_unit: str, *, field: str) -> pint.Quantity[Any]:
    if not isinstance(value, pint.Quantity):
        raise _error("invalid-unit", field=field, expected=target_unit)
    try:
        magnitude = _decimal_number(value.magnitude, field=field)
        declared_unit = str(value.units)
        quantity = unit_registry.Quantity(magnitude, declared_unit).to(target_unit)
    except ParameterError:
        raise
    except (AttributeError, TypeError, ValueError, pint.errors.PintError) as exc:
        raise _error("invalid-unit", field=field, expected=target_unit) from exc
    return quantity


def _declared_quantity(value: object, *, field: str) -> pint.Quantity[Any]:
    if not isinstance(value, pint.Quantity):
        raise _error("invalid-unit", field=field, expected="quantity")
    try:
        magnitude = _decimal_number(value.magnitude, field=field)
        declared_unit = str(value.units)
        quantity = unit_registry.Quantity(magnitude, declared_unit)
    except ParameterError:
        raise
    except (AttributeError, TypeError, ValueError, pint.errors.PintError) as exc:
        raise _error("invalid-unit", field=field, expected="quantity") from exc
    return quantity


__all__ = ("unit_registry",)
