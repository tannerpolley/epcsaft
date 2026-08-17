from __future__ import annotations

import hashlib
import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from importlib import metadata
from os import PathLike
from pathlib import Path
from typing import Any, cast

from .. import _core
from .._trials import NON_EVALUABLE_TRIAL_CODES, NonEvaluableTrialCode
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
from ._result import FitResult, FitSpecificationReceipt

_MAXIMUM_ITERATIONS = 25
_FUNCTION_TOLERANCE = 1.0e-6
_GRADIENT_TOLERANCE = 1.0e-6
_PARAMETER_TOLERANCE = 1.0e-6


class FitPreflightError(ValueError):
    """A typed non-evaluable exact start encountered before native optimization."""

    def __init__(
        self,
        code: NonEvaluableTrialCode,
        start_identity: str,
        batch_identity: str,
        diagnostic: str,
    ) -> None:
        if code not in NON_EVALUABLE_TRIAL_CODES:
            raise ValueError("preflight failure code is outside the public vocabulary")
        self.code = code
        self.start_identity = start_identity
        self.batch_identity = batch_identity
        self.diagnostic = diagnostic
        super().__init__(
            f"{code} during derivative preflight for {start_identity}/{batch_identity}: "
            f"{diagnostic}"
        )


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

    def __call__(
        self,
        _start_index: int,
        index: int,
        values: Sequence[float],
        _receipt_requested: bool,
    ) -> dict[str, object]:
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
            return self._rejected("provider_domain_rejection", f"{exc.code}: {exc}")
        except SaturationError as exc:
            return self._rejected(
                "physical_rejection",
                f"{exc.diagnostics.rejection_status}: {exc}",
            )
        except (TypeError, ValueError) as exc:
            return self._rejected("parameter_topology_mismatch", str(exc))
        if block.status.status != "accepted":
            return self._rejected(
                "physical_rejection"
                if block.status.scope in {"shared-solve", "branch"}
                else "observation_out_of_support",
                block.status.reason,
            )
        if block.descriptor.active_parameter_identities != self._selected_identities:
            return self._rejected(
                "owner_contract_failure",
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
                "owner_contract_failure",
                "owner block did not return every selected row",
            )
        accepted = cast(tuple[tuple[float, tuple[float, ...], object], ...], selected)
        return {
            "status": "accepted",
            "predictions": [row[0] for row in accepted],
            "jacobian": [value for row in accepted for value in row[1]],
        }


def _preflight(
    owner: _OwnerEvaluator,
    resolved: ResolvedFit,
    requests: Sequence[_OwnerRequest],
) -> None:
    coordinate_by_parameter = {
        identity: index
        for index, coordinate in enumerate(resolved.coordinates)
        for identity in coordinate.tied_identities
    }
    for start in resolved.fit.solver.multistart.starts:
        physical = start.values
        active = tuple(
            physical[coordinate_by_parameter[identity]] for identity in resolved.selected_identities
        )
        for index, request in enumerate(requests):
            response = owner(0, index, active, False)
            if response.get("status") != "accepted":
                failure = response.get("failure")
                failure_mapping = failure if isinstance(failure, Mapping) else {}
                raw_code = str(failure_mapping.get("code", "owner_contract_failure"))
                code = (
                    raw_code if raw_code in NON_EVALUABLE_TRIAL_CODES else "owner_contract_failure"
                )
                raise FitPreflightError(
                    cast(NonEvaluableTrialCode, code),
                    start.identity,
                    f"observation-{index}",
                    str(failure_mapping.get("reason", failure)),
                )
            predictions = response.get("predictions")
            jacobian = response.get("jacobian")
            if (
                not isinstance(predictions, Sequence)
                or isinstance(predictions, (str, bytes))
                or not isinstance(jacobian, Sequence)
                or isinstance(jacobian, (str, bytes))
                or len(predictions) != len(request.row_identities)
                or len(jacobian) != len(request.row_identities) * len(active)
                or any(not math.isfinite(float(value)) for value in predictions)
                or any(not math.isfinite(float(value)) for value in jacobian)
            ):
                raise FitPreflightError(
                    "owner_contract_failure",
                    start.identity,
                    f"observation-{index}",
                    "owner returned an incomplete value/Jacobian block",
                )


def _package_fingerprint() -> str:
    record = metadata.distribution("epcsaft").read_text("RECORD")
    if record is None or not record.strip():
        raise RuntimeError("installed epcsaft distribution has no RECORD identity")
    return f"sha256:{hashlib.sha256(record.encode()).hexdigest()}"


def fit(
    value: object,
    *,
    base_path: str | PathLike[str] | None = None,
) -> FitResult:
    if not isinstance(value, (ParameterFit, Mapping, str, PathLike)):
        raise TypeError("fit requires a ParameterFit, mapping, or path")
    specification, source = _load(value)
    root = _root(base_path, source)
    parameter_path = Path(specification.parameters.path)
    if not parameter_path.is_absolute():
        parameter_path = root / parameter_path
    parameters = Parameters.from_json(parameter_path)
    resolved = specification.resolve(parameters, base_path=root)
    if not resolved.observations:
        from ._equilibrium_fit import _fit_equilibrium_parameter_fit

        return _fit_equilibrium_parameter_fit(resolved)
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
    owner = _OwnerEvaluator(model, resolved, requests)
    _preflight(owner, resolved, requests)
    native = _core._fit_resolved(
        owner,
        len(resolved.selected_identities),
        [coordinate_indices[identity] for identity in resolved.selected_identities],
        [
            {
                "origin": coordinate.origin,
                "scale": coordinate.scale,
                "lower_affine": coordinate.affine_bounds[0],
                "upper_affine": coordinate.affine_bounds[1],
            }
            for coordinate in resolved.coordinates
        ],
        residual_blocks,
        {
            "starts": [
                {
                    "identity": start.identity,
                    "affine": [
                        (value - coordinate.origin) / coordinate.scale
                        for value, coordinate in zip(
                            start.values,
                            resolved.coordinates,
                            strict=True,
                        )
                    ],
                }
                for start in specification.solver.multistart.starts
            ],
            "affine_bound_margin": specification.solver.multistart.affine_bound_margin,
            "workers": 1,
            "maximum_iterations": specification.solver.controls.maximum_iterations,
            "function_tolerance": specification.solver.controls.function_tolerance,
            "gradient_tolerance": specification.solver.controls.gradient_tolerance,
            "parameter_tolerance": specification.solver.controls.parameter_tolerance,
        },
        tuple(sorted(NON_EVALUABLE_TRIAL_CODES)),
    )
    return FitResult.from_native(
        FitSpecificationReceipt.from_resolved_fit(resolved),
        cast(Mapping[str, Any], native),
        artifact_fingerprint=_python_native_artifact_fingerprint(),
        package_fingerprint=_package_fingerprint(),
    )


__all__ = ("FitPreflightError", "fit")
