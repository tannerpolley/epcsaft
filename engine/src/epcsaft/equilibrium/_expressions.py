"""Small, declarative algebra for derived equilibrium output rows.

The expression layer deliberately operates only on rows that an equilibrium
owner has already evaluated.  It does not evaluate thermodynamics, discover
species, or contain source-specific observation logic.  Its Jacobian is the
exact chain-rule Jacobian with respect to the ordered coordinates carried by
the input rows.
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from typing import Final, Protocol

_DEFAULT_DENOMINATOR_FLOOR: Final[float] = 1.0e-14


class ExpressionError(ValueError):
    """Base error for invalid expression definitions or evaluations."""


class ExpressionEvaluationError(ExpressionError):
    """Raised when rows or an evaluated expression are not scientifically valid."""


class EvaluatedRowLike(Protocol):
    """The already-evaluated row surface consumed by this module."""

    identity: str
    status: str
    value: float | None
    jacobian: Sequence[float] | None


@dataclass(frozen=True, slots=True)
class EvaluatedRow:
    """A compact row value and exact Jacobian supplied by an owner module.

    ``status`` is intentionally retained even though successful rows are the
    only rows this algebra can consume.  A failed or pending owner row can be
    passed through and is rejected with a structured expression error rather
    than being silently turned into a number.
    """

    identity: str
    value: float | None
    jacobian: tuple[float, ...] | None
    status: str = "evaluated"


@dataclass(frozen=True, slots=True)
class Constant:
    """A finite scalar constant."""

    value: float

    def __post_init__(self) -> None:
        object.__setattr__(self, "value", _finite_scalar(self.value, "constant"))


@dataclass(frozen=True, slots=True)
class RowReference:
    """Reference one owner-produced row by its stable identity."""

    identity: str

    def __post_init__(self) -> None:
        if not isinstance(self.identity, str) or not self.identity.strip():
            raise ExpressionError("row reference identity must be a non-empty string")


@dataclass(frozen=True, slots=True)
class LinearSum:
    """A finite linear combination of expression nodes plus a constant."""

    terms: tuple[tuple[float, Expression], ...] = ()
    constant: float = 0.0

    def __post_init__(self) -> None:
        normalized: list[tuple[float, Expression]] = []
        for term in self.terms:
            try:
                coefficient, expression = term
            except (TypeError, ValueError) as error:
                raise ExpressionError(
                    "linear-sum terms must be (coefficient, expression) pairs"
                ) from error
            normalized.append(
                (_finite_scalar(coefficient, "linear-sum coefficient"), _ensure_expression(expression))
            )
        object.__setattr__(self, "terms", tuple(normalized))
        object.__setattr__(self, "constant", _finite_scalar(self.constant, "linear-sum constant"))


@dataclass(frozen=True, slots=True)
class Product:
    """Product of two expression nodes."""

    left: Expression
    right: Expression

    def __post_init__(self) -> None:
        _ensure_expression(self.left)
        _ensure_expression(self.right)


@dataclass(frozen=True, slots=True)
class Ratio:
    """Ratio of two expression nodes with an explicit near-zero guard."""

    numerator: Expression
    denominator: Expression
    denominator_floor: float = _DEFAULT_DENOMINATOR_FLOOR

    def __post_init__(self) -> None:
        _ensure_expression(self.numerator)
        _ensure_expression(self.denominator)
        floor = _finite_scalar(self.denominator_floor, "denominator floor")
        if floor <= 0.0:
            raise ExpressionError("denominator floor must be positive")
        object.__setattr__(self, "denominator_floor", floor)


type Expression = Constant | RowReference | LinearSum | Product | Ratio


@dataclass(frozen=True, slots=True)
class EvaluatedExpression:
    """Value and exact total Jacobian returned by :func:`evaluate_expression`."""

    value: float
    jacobian: tuple[float, ...]

    def __post_init__(self) -> None:
        value = _finite_scalar(self.value, "expression value")
        jacobian = tuple(_finite_scalar(item, "expression Jacobian entry") for item in self.jacobian)
        object.__setattr__(self, "value", value)
        object.__setattr__(self, "jacobian", jacobian)


@dataclass(frozen=True, slots=True)
class _Dual:
    value: float
    jacobian: tuple[float, ...]


def evaluate_expression(
    expression: Expression,
    rows: Iterable[EvaluatedRowLike],
) -> EvaluatedExpression:
    """Evaluate a strict expression tree over owner-produced rows.

    Rows establish one common ordered Jacobian coordinate basis.  A constant
    expression can be evaluated with no rows and consequently has an empty
    Jacobian.  Missing references, duplicate row identities, failed rows,
    inconsistent derivative dimensions, non-finite inputs/results, and
    near-zero ratio denominators all fail closed.
    """

    expression = _ensure_expression(expression)
    indexed, dimension = _index_rows(rows)
    result = _evaluate(expression, indexed, dimension)
    if not math.isfinite(result.value) or any(not math.isfinite(item) for item in result.jacobian):
        raise ExpressionEvaluationError("expression result is non-finite")
    return EvaluatedExpression(result.value, result.jacobian)


def _finite_scalar(value: float, label: str) -> float:
    try:
        converted = float(value)
    except (TypeError, ValueError, OverflowError) as error:
        raise ExpressionError(f"{label} must be a finite real number") from error
    if not math.isfinite(converted):
        raise ExpressionError(f"{label} must be finite")
    return converted


def _ensure_expression(expression: object) -> Expression:
    if not isinstance(expression, (Constant, RowReference, LinearSum, Product, Ratio)):
        raise ExpressionError("expression contains an unsupported node")
    return expression


def _index_rows(rows: Iterable[EvaluatedRowLike]) -> tuple[dict[str, _Dual], int]:
    indexed: dict[str, _Dual] = {}
    dimension: int | None = None
    try:
        iterator = iter(rows)
    except TypeError as error:
        raise ExpressionEvaluationError("rows must be an iterable of evaluated rows") from error

    for row in iterator:
        try:
            identity = row.identity
            status = row.status
            value = row.value
            jacobian = row.jacobian
        except AttributeError as error:
            raise ExpressionEvaluationError(
                "each row must expose identity, status, value, and jacobian"
            ) from error
        if not isinstance(identity, str) or not identity.strip():
            raise ExpressionEvaluationError("row identity must be a non-empty string")
        if identity in indexed:
            raise ExpressionEvaluationError(f"duplicate row identity: {identity!r}")
        if status != "evaluated" or value is None or jacobian is None:
            raise ExpressionEvaluationError(f"row {identity!r} is not evaluable")
        try:
            row_value = float(value)
            row_jacobian = tuple(float(item) for item in jacobian)
        except (TypeError, ValueError, OverflowError) as error:
            raise ExpressionEvaluationError(f"row {identity!r} is not numeric") from error
        if not math.isfinite(row_value) or any(not math.isfinite(item) for item in row_jacobian):
            raise ExpressionEvaluationError(f"row {identity!r} contains a non-finite value")
        row_dimension = len(row_jacobian)
        if dimension is None:
            dimension = row_dimension
        elif row_dimension != dimension:
            raise ExpressionEvaluationError("row Jacobians have inconsistent dimensions")
        indexed[identity] = _Dual(row_value, row_jacobian)
    return indexed, 0 if dimension is None else dimension


def _zeros(dimension: int) -> tuple[float, ...]:
    return (0.0,) * dimension


def _evaluate(expression: Expression, rows: dict[str, _Dual], dimension: int) -> _Dual:
    if isinstance(expression, Constant):
        return _Dual(expression.value, _zeros(dimension))
    if isinstance(expression, RowReference):
        try:
            return rows[expression.identity]
        except KeyError as error:
            raise ExpressionEvaluationError(
                f"missing evaluated row: {expression.identity!r}"
            ) from error
    if isinstance(expression, LinearSum):
        value = expression.constant
        linear_jacobian = list(_zeros(dimension))
        for coefficient, child in expression.terms:
            evaluated = _evaluate(child, rows, dimension)
            value += coefficient * evaluated.value
            for index, derivative in enumerate(evaluated.jacobian):
                linear_jacobian[index] += coefficient * derivative
        return _checked_dual(value, linear_jacobian)
    if isinstance(expression, Product):
        left = _evaluate(expression.left, rows, dimension)
        right = _evaluate(expression.right, rows, dimension)
        value = left.value * right.value
        product_jacobian = tuple(
            left_derivative * right.value + right_derivative * left.value
            for left_derivative, right_derivative in zip(left.jacobian, right.jacobian, strict=True)
        )
        return _checked_dual(value, product_jacobian)
    numerator = _evaluate(expression.numerator, rows, dimension)
    denominator = _evaluate(expression.denominator, rows, dimension)
    if abs(denominator.value) <= expression.denominator_floor:
        raise ExpressionEvaluationError("ratio denominator is zero or near zero")
    denominator_squared = denominator.value * denominator.value
    value = numerator.value / denominator.value
    ratio_jacobian = tuple(
        (numerator_derivative * denominator.value - denominator_derivative * numerator.value)
        / denominator_squared
        for numerator_derivative, denominator_derivative in zip(
            numerator.jacobian, denominator.jacobian, strict=True
        )
    )
    return _checked_dual(value, ratio_jacobian)


def _checked_dual(value: float, jacobian: Iterable[float]) -> _Dual:
    normalized = tuple(jacobian)
    if not math.isfinite(value) or any(not math.isfinite(item) for item in normalized):
        raise ExpressionEvaluationError("expression intermediate is non-finite")
    return _Dual(value, normalized)


__all__ = [
    "Constant",
    "EvaluatedExpression",
    "EvaluatedRow",
    "EvaluatedRowLike",
    "Expression",
    "ExpressionError",
    "ExpressionEvaluationError",
    "LinearSum",
    "Product",
    "Ratio",
    "RowReference",
    "evaluate_expression",
]
