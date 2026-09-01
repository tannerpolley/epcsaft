"""Owner-neutral Ceres regression entry points.

The public layer only assembles immutable request data.  Thermodynamic and
equilibrium evaluation remains in the native owner evaluators, so a fit never
invokes a Python callback while Ceres is running.
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from types import MappingProxyType
from typing import Any, cast

from epcsaft import _core


def _freeze(value: object) -> object:
    if isinstance(value, Mapping):
        return MappingProxyType({str(key): _freeze(item) for key, item in value.items()})
    if isinstance(value, (list, tuple)):
        return tuple(_freeze(item) for item in value)
    return value


def _vector(value: Iterable[Any], name: str) -> tuple[float, ...]:
    try:
        result = tuple(float(item) for item in value)
    except (TypeError, ValueError) as error:
        raise TypeError(f"{name} must be a numeric sequence") from error
    if any(not math.isfinite(item) for item in result):
        raise ValueError(f"{name} must be finite")
    return result


def _bound_vector(value: Iterable[Any], name: str) -> tuple[float, ...]:
    try:
        result = tuple(float(item) for item in value)
    except (TypeError, ValueError) as error:
        raise TypeError(f"{name} must be a numeric sequence") from error
    if any(math.isnan(item) for item in result):
        raise ValueError(f"{name} must not contain NaN")
    return result


def _native_evaluators(evaluators: Iterable[CompiledEvaluator]) -> tuple[_core.CompiledOwner, ...]:
    result = tuple(evaluators)
    if not result or any(not isinstance(item, CompiledEvaluator) for item in result):
        raise TypeError("evaluators must contain one or more compiled evaluators")
    return tuple(cast(_core.CompiledOwner, item._native) for item in result)


@dataclass(frozen=True, slots=True)
class CompiledEvaluator:
    """One immutable opaque native owner and its audit sidecar."""

    _native: Any
    sidecar: Mapping[str, object]

    @classmethod
    def from_native(cls, native: Any) -> CompiledEvaluator:
        sidecar = cast(Mapping[str, object], native.sidecar)
        return cls(native, cast(Mapping[str, object], _freeze(sidecar)))

    @property
    def parameter_count(self) -> int:
        return int(self._native.parameter_count)

    @property
    def output_count(self) -> int:
        return int(self._native.output_count)

    @property
    def parameter_contract_fingerprint(self) -> str:
        return str(self._native.parameter_contract_fingerprint)

    @property
    def evaluator_fingerprint(self) -> str:
        return str(self._native.evaluator_fingerprint)


@dataclass(frozen=True, slots=True)
class LeastSquaresResult:
    """Fit values plus the complete native diagnostic payload."""

    status: str
    termination: str
    message: str
    usable: bool
    parameter_contract_fingerprint: str
    residual_scale_metadata: str
    initial_parameters: tuple[float, ...]
    lower_bounds: tuple[float, ...]
    upper_bounds: tuple[float, ...]
    affine_origins: tuple[float, ...]
    affine_scales: tuple[float, ...]
    parameters: tuple[float, ...]
    observations: tuple[float, ...]
    predictions: tuple[float, ...]
    raw_residuals: tuple[float, ...]
    residuals: tuple[float, ...]
    residual_scales: tuple[float, ...]
    physical_jacobian: tuple[tuple[float, ...], ...]
    optimizer_jacobian: tuple[tuple[float, ...], ...]
    evaluator_sidecars: tuple[Mapping[str, object], ...]
    evaluator_rows: tuple[Mapping[str, object], ...]
    failures: tuple[Mapping[str, object], ...]
    requested_options: Mapping[str, object]
    resolved_options: Mapping[str, object]
    counters: Mapping[str, object]
    covariance: Mapping[str, object]

    @property
    def success(self) -> bool:
        return self.usable

    @property
    def x(self) -> tuple[float, ...]:
        """SciPy-compatible fitted coordinate spelling."""
        return self.parameters

    @property
    def cost(self) -> float:
        return 0.5 * math.fsum(value * value for value in self.residuals)


def _matrix(values: Sequence[Any], rows: int, columns: int) -> tuple[tuple[float, ...], ...]:
    flat = tuple(float(value) for value in values)
    if not flat:
        return ()
    if len(flat) != rows * columns:
        raise RuntimeError("native Jacobian dimensions are inconsistent")
    return tuple(
        flat[row * columns : (row + 1) * columns] for row in range(rows)
    )


def _result(
    native: Mapping[str, object],
    observations: tuple[float, ...],
    residual_scales: tuple[float, ...],
    parameter_count: int,
    initial_parameters: tuple[float, ...],
    lower_bounds: tuple[float, ...],
    upper_bounds: tuple[float, ...],
    affine_origins: tuple[float, ...],
    affine_scales: tuple[float, ...],
    evaluator_sidecars: tuple[Mapping[str, object], ...],
) -> LeastSquaresResult:
    predictions = tuple(float(value) for value in cast(Sequence[Any], native["predictions"]))
    frozen_rows = tuple(
        cast(Mapping[str, object], _freeze(row))
        for row in cast(Sequence[Any], native["evaluator_rows"])
    )
    failures = tuple(
        cast(Mapping[str, object], _freeze(row))
        for row in cast(Sequence[Any], native["failures"])
    )
    requested = cast(Mapping[str, object], _freeze(native["requested_options"]))
    resolved = cast(Mapping[str, object], _freeze(native["resolved_options"]))
    counters = cast(Mapping[str, object], _freeze(native["counters"]))
    covariance = cast(Mapping[str, object], _freeze(native["covariance"]))
    raw_residuals = tuple(prediction - observation for prediction, observation in zip(
        predictions, observations, strict=True
    )) if predictions else ()
    residuals = tuple(float(value) for value in cast(Sequence[Any], native["residuals"]))
    physical_jacobian = _matrix(
        cast(Sequence[object], native["physical_jacobian"]), len(predictions), parameter_count
    )
    optimizer_jacobian = _matrix(
        cast(Sequence[object], native["optimizer_jacobian"]), len(predictions), parameter_count
    )
    return LeastSquaresResult(
        str(native["status"]),
        str(native["termination"]),
        str(native["message"]),
        bool(native["usable"]),
        str(native["parameter_contract_fingerprint"]),
        str(native["residual_scale_metadata"]),
        initial_parameters,
        lower_bounds,
        upper_bounds,
        affine_origins,
        affine_scales,
        tuple(float(value) for value in cast(Sequence[Any], native["parameters"])),
        observations,
        predictions,
        raw_residuals,
        residuals,
        residual_scales,
        physical_jacobian,
        optimizer_jacobian,
        evaluator_sidecars,
        frozen_rows,
        failures,
        requested,
        resolved,
        counters,
        covariance,
    )


def least_squares(
    evaluators: Iterable[CompiledEvaluator],
    observations: Iterable[Any],
    initial_parameters: Iterable[Any],
    *,
    bounds: tuple[Iterable[Any], Iterable[Any]] | None = None,
    affine_origin: Iterable[Any] | None = None,
    affine_scale: Iterable[Any] | None = None,
    residual_scales: Iterable[Any] | None = None,
    residual_scale_metadata: str = "caller-defined discrepancy weights",
    function_tolerance: float = 1.0e-12,
    gradient_tolerance: float = 1.0e-12,
    parameter_tolerance: float = 1.0e-12,
    maximum_iterations: int = 1000,
    maximum_elapsed_time_seconds: float = 300.0,
    maximum_uncached_owner_calls: int = 10000,
) -> LeastSquaresResult:
    """Solve ordered owner observations with native Ceres and exact Jacobians."""

    compiled = tuple(evaluators)
    native_evaluators = _native_evaluators(compiled)
    initial = _vector(initial_parameters, "initial_parameters")
    if not initial:
        raise ValueError("initial_parameters must not be empty")
    parameter_count = len(initial)
    if any(item.parameter_count != parameter_count for item in compiled):
        raise ValueError("all evaluators must use the initial parameter dimension")
    y = _vector(observations, "observations")
    if not y:
        raise ValueError("observations must not be empty")
    if bounds is not None:
        lower = _bound_vector(bounds[0], "lower_bounds")
        upper = _bound_vector(bounds[1], "upper_bounds")
    else:
        lower = tuple(-math.inf for _ in initial)
        upper = tuple(math.inf for _ in initial)
    if len(lower) != parameter_count or len(upper) != parameter_count:
        raise ValueError("bounds must match initial_parameters")
    if any(math.isnan(value) for value in (*lower, *upper)):
        raise ValueError("bounds must not contain NaN")
    origins = (
        _vector(affine_origin, "affine_origin")
        if affine_origin is not None
        else tuple(0.0 for _ in initial)
    )
    scales = (
        _vector(affine_scale, "affine_scale")
        if affine_scale is not None
        else tuple(1.0 for _ in initial)
    )
    if len(origins) != parameter_count or len(scales) != parameter_count:
        raise ValueError("affine coordinates must match initial_parameters")
    if any(scale <= 0.0 for scale in scales):
        raise ValueError("affine_scale must be positive")
    scale_values = (
        _vector(residual_scales, "residual_scales")
        if residual_scales is not None
        else ()
    )
    if scale_values and len(scale_values) not in {1, len(y)}:
        raise ValueError("residual_scales must be scalar or match observations")
    expanded_scales = (
        tuple(scale_values[0] for _ in y) if len(scale_values) == 1 else scale_values
    ) or tuple(1.0 for _ in y)
    if any(scale <= 0.0 for scale in expanded_scales):
        raise ValueError("residual_scales must be positive")
    if not isinstance(residual_scale_metadata, str) or not residual_scale_metadata:
        raise ValueError("residual_scale_metadata must be non-empty text")
    controls = {
        "function_tolerance": float(function_tolerance),
        "gradient_tolerance": float(gradient_tolerance),
        "parameter_tolerance": float(parameter_tolerance),
        "maximum_iterations": int(maximum_iterations),
        "maximum_elapsed_time_seconds": float(maximum_elapsed_time_seconds),
        "maximum_uncached_owner_calls": int(maximum_uncached_owner_calls),
    }
    native_result = _core._solve_least_squares(
        native_evaluators,
        y,
        scale_values,
        residual_scale_metadata,
        initial,
        lower,
        upper,
        origins,
        scales,
        controls,
    )
    return _result(
        native_result,
        y,
        expanded_scales,
        parameter_count,
        initial,
        lower,
        upper,
        origins,
        scales,
        tuple(item.sidecar for item in compiled),
    )


def curve_fit(
    evaluators: Iterable[CompiledEvaluator],
    observations: Iterable[object],
    initial_parameters: Iterable[object],
    **kwargs: Any,
) -> LeastSquaresResult:
    """Small curve-fitting spelling over the same owner-neutral least-squares core."""

    return least_squares(evaluators, observations, initial_parameters, **kwargs)
