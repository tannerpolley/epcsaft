from __future__ import annotations

import hashlib
import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from importlib import metadata
from os import PathLike
from pathlib import Path
from typing import Any, SupportsFloat, cast

import pint

from .. import _core
from ..eos import (
    OBSERVABLE_TARGETS,
    EosError,
    Mixture,
    ObservableDescriptor,
    Pure,
    SaturationError,
    _python_native_artifact_fingerprint,
    evaluate_observable,
)
from ..parameters import ActiveParameterSet, Parameters
from ._model import ParameterFit, TableRow
from ._resolve import ResolvedFit, ResolvedObservation
from ._result import FitResult

_MAXIMUM_ITERATIONS = 100
_FUNCTION_TOLERANCE = 1.0e-12
_GRADIENT_TOLERANCE = 1.0e-12
_PARAMETER_TOLERANCE = 1.0e-12


def _root(base_path: str | PathLike[str] | None, source: Path | None) -> Path:
    candidate = Path(base_path) if base_path is not None else source or Path.cwd()
    return candidate if candidate.is_dir() else candidate.parent


def _load(
    value: ParameterFit | Mapping[str, object] | str | PathLike[str],
) -> tuple[ParameterFit, Path | None]:
    if isinstance(value, ParameterFit):
        return value, None
    if isinstance(value, Mapping):
        return ParameterFit.from_mapping(value), None
    source = Path(value)
    return ParameterFit.from_json(source), source


def _parameter_input(resolved: ResolvedFit, identity: str, unit: str) -> float:
    specification = next(
        item for item in resolved.parameters.parameter_specs if item.identity == identity
    )
    value = specification.value
    return (
        float(value.to(unit).magnitude)
        if isinstance(value, pint.Quantity)
        else float(cast(SupportsFloat, value))
    )


@dataclass(frozen=True, slots=True)
class _OwnerRequest:
    descriptor: ObservableDescriptor
    state: tuple[tuple[str, float], ...]
    composition: tuple[float, ...] | None
    row_identities: tuple[str, ...]


def _rows(
    observation: ResolvedObservation,
) -> tuple[tuple[TableRow, float], ...]:
    if observation.rows is not None:
        return tuple((row, observation.specification.multiplier) for row in observation.rows)
    if observation.quadrature is None:
        raise RuntimeError("resolved observation has neither rows nor quadrature")
    rule = observation.quadrature.rules[-1]
    target = observation.specification.targets[0]
    return tuple(
        (
            TableRow((("temperature_k", temperature), (target, observed))),
            observation.specification.multiplier * math.sqrt(weight),
        )
        for temperature, weight, observed in zip(
            rule.nodes, rule.weights, rule.observed_values, strict=True
        )
    )


def _owner_requests(
    resolved: ResolvedFit,
) -> tuple[tuple[_OwnerRequest, ...], list[dict[str, object]]]:
    requests: list[_OwnerRequest] = []
    residual_blocks: list[dict[str, object]] = []
    for observation in resolved.observations:
        specification = observation.specification
        state_fields = tuple(name for name, _ in specification.source.domain.bounds)
        residuals = dict(specification.residuals)
        row_identities = tuple(
            OBSERVABLE_TARGETS[target].row_identity for target in specification.targets
        )
        for row, multiplier in _rows(observation):
            requests.append(
                _OwnerRequest(
                    observation.descriptor,
                    tuple((name, row.value(name)) for name in state_fields),
                    specification.composition,
                    row_identities,
                )
            )
            residual_blocks.append(
                {
                    "multiplier": multiplier,
                    "targets": [
                        {
                            "observed": row.value(target),
                            "residual": residuals[target],
                            "scale": row.value(target),
                        }
                        for target in specification.targets
                    ],
                }
            )
    return tuple(requests), residual_blocks


class _OwnerEvaluator:
    __slots__ = ("_model", "_parameters", "_requests", "_selected_identities")

    def __init__(
        self,
        model: Pure | Mixture,
        resolved: ResolvedFit,
        requests: Sequence[_OwnerRequest],
    ) -> None:
        self._model = model
        self._parameters = resolved.parameters
        self._selected_identities = resolved.selected_identities
        self._requests = tuple(requests)

    @staticmethod
    def _rejected(code: str, reason: str) -> dict[str, object]:
        return {
            "status": "rejected",
            "failure": {"code": code, "reason": reason},
        }

    def __call__(self, index: int, values: Sequence[float]) -> dict[str, object]:
        request = self._requests[index]
        try:
            active_parameters = ActiveParameterSet(
                self._parameters,
                self._selected_identities,
                values=values,
            )
            block = evaluate_observable(
                self._model,
                active_parameters,
                request.descriptor,
                dict(request.state),
                composition=request.composition,
            )
        except EosError as exc:
            return self._rejected(exc.code, str(exc))
        except SaturationError as exc:
            return self._rejected(f"saturation-{exc.diagnostics.rejection_status}", str(exc))
        except (TypeError, ValueError) as exc:
            return self._rejected("invalid-active-parameter-value", str(exc))
        if block.status.status != "accepted":
            code = (
                f"saturation-{block.saturation_diagnostics.rejection_status}"
                if block.saturation_diagnostics is not None
                else f"observable-{block.status.scope}-rejected"
            )
            return self._rejected(code, block.status.reason)
        if block.descriptor.active_parameter_identities != self._selected_identities:
            return self._rejected(
                "owner-parameter-order-mismatch",
                "owner block did not preserve the selected parameter order",
            )
        rows = {
            identity: (value, jacobian, status)
            for identity, value, jacobian, status in zip(
                block.descriptor.row_identities,
                block.values,
                block.jacobian,
                block.row_statuses,
                strict=True,
            )
        }
        selected = tuple(rows.get(identity) for identity in request.row_identities)
        if any(
            row is None or row[0] is None or row[1] is None or row[2].status != "accepted"
            for row in selected
        ):
            return self._rejected(
                "owner-row-rejected",
                "owner block did not return every selected row",
            )
        accepted = cast(tuple[tuple[float, tuple[float, ...], object], ...], selected)
        return {
            "status": "accepted",
            "predictions": [row[0] for row in accepted],
            "jacobian": [value for row in accepted for value in row[1]],
        }


def _package_fingerprint() -> str:
    record = metadata.distribution("epcsaft").read_text("RECORD")
    if record is None or not record.strip():
        raise RuntimeError("installed epcsaft distribution has no RECORD identity")
    return f"sha256:{hashlib.sha256(record.encode()).hexdigest()}"


def fit(
    value: ParameterFit | Mapping[str, object] | str | PathLike[str],
    *,
    base_path: str | PathLike[str] | None = None,
) -> FitResult:
    specification, source = _load(value)
    root = _root(base_path, source)
    parameter_path = Path(specification.parameters.path)
    if not parameter_path.is_absolute():
        parameter_path = root / parameter_path
    parameters = Parameters.from_json(parameter_path)
    resolved = specification.resolve(parameters, base_path=root)
    model: Pure | Mixture = (
        Pure(resolved.parameters)
        if specification.model.kind == "pure"
        else Mixture(resolved.parameters)
    )
    coordinate_indices = {
        identity: index
        for index, coordinate in enumerate(resolved.coordinates)
        for identity in coordinate.tied_identities
    }
    requests, residual_blocks = _owner_requests(resolved)
    native = _core._fit_resolved(
        _OwnerEvaluator(model, resolved, requests),
        len(resolved.selected_identities),
        [coordinate_indices[identity] for identity in resolved.selected_identities],
        [
            {
                "scale": coordinate.scale,
                "lower": coordinate.bounds[0],
                "upper": coordinate.bounds[1],
                "input": _parameter_input(resolved, coordinate.identity, coordinate.unit),
            }
            for coordinate in resolved.coordinates
        ],
        residual_blocks,
        {
            "start_count": specification.solver.multistart.count,
            "include_input": specification.solver.multistart.include_input,
            "workers": 1,
            "maximum_iterations": _MAXIMUM_ITERATIONS,
            "function_tolerance": _FUNCTION_TOLERANCE,
            "gradient_tolerance": _GRADIENT_TOLERANCE,
            "parameter_tolerance": _PARAMETER_TOLERANCE,
        },
    )
    return FitResult.from_native(
        resolved,
        cast(Mapping[str, Any], native),
        artifact_fingerprint=_python_native_artifact_fingerprint(),
        package_fingerprint=_package_fingerprint(),
    )


__all__ = ("fit",)
