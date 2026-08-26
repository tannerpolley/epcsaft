from __future__ import annotations

import hashlib
import math
from collections.abc import Callable, Mapping, Sequence
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
from ._model import EquilibriumFitObservation, ParameterFit, TableRow
from ._resolve import ResolvedFit, ResolvedObservation
from ._result import FitResult, FitRowReceipt, FitSpecificationReceipt


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
) -> dict[tuple[str, int, tuple[str, ...], tuple[float, ...]], dict[str, object]]:
    exact_start_blocks: dict[
        tuple[str, int, tuple[str, ...], tuple[float, ...]], dict[str, object]
    ] = {}
    coordinate_by_parameter = {
        identity: index
        for index, coordinate in enumerate(resolved.coordinates)
        for identity in coordinate.tied_identities
    }
    for start_index, start in enumerate(resolved.fit.solver.multistart.starts):
        physical = start.values
        active = tuple(
            physical[coordinate_by_parameter[identity]] for identity in resolved.selected_identities
        )
        for index, request in enumerate(requests):
            response = owner(start_index, index, active, False)
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
            exact_start_blocks[(start.identity, index, resolved.selected_identities, active)] = (
                response
            )
    return exact_start_blocks


def _package_fingerprint() -> str:
    record = metadata.distribution("epcsaft").read_text("RECORD")
    if record is None or not record.strip():
        raise RuntimeError("installed epcsaft distribution has no RECORD identity")
    return f"sha256:{hashlib.sha256(record.encode()).hexdigest()}"


def _fit_owner_batches(
    resolved: ResolvedFit,
    owner: Callable[[int, int, Sequence[float], bool], Mapping[str, object]],
    residual_blocks: Sequence[Mapping[str, object]],
    specification: FitSpecificationReceipt,
    full_space: tuple[
        Callable[[int, Sequence[float]], Mapping[str, object]],
        Sequence[Mapping[str, object]],
    ]
    | None = None,
) -> FitResult:
    coordinate_indices = {
        identity: index
        for index, coordinate in enumerate(resolved.coordinates)
        for identity in coordinate.tied_identities
    }
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
                for start in resolved.fit.solver.multistart.starts
            ],
            "affine_bound_margin": resolved.fit.solver.multistart.affine_bound_margin,
            "maximum_iterations": resolved.fit.solver.controls.maximum_iterations,
            "function_tolerance": resolved.fit.solver.controls.function_tolerance,
            "gradient_tolerance": resolved.fit.solver.controls.gradient_tolerance,
            "parameter_tolerance": resolved.fit.solver.controls.parameter_tolerance,
        },
        tuple(sorted(NON_EVALUABLE_TRIAL_CODES)),
        None if full_space is None else full_space[0],
        () if full_space is None else full_space[1],
    )
    return FitResult.from_native(
        specification,
        cast(Mapping[str, Any], native),
        artifact_fingerprint=_python_native_artifact_fingerprint(),
        package_fingerprint=_package_fingerprint(),
    )


def _receipt_blocks(
    receipt: FitSpecificationReceipt,
) -> tuple[tuple[FitRowReceipt, ...], ...]:
    offset = 0
    blocks = []
    for size in receipt.block_sizes:
        blocks.append(receipt.rows[offset : offset + size])
        offset += size
    return tuple(blocks)


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
    model: Pure | Mixture = (
        Pure(resolved.parameters)
        if specification.model.kind == "pure"
        else Mixture(resolved.parameters)
    )
    requests, residual_blocks = _owner_requests(resolved)
    direct_owner = _OwnerEvaluator(model, resolved, requests)
    exact_start_blocks = _preflight(direct_owner, resolved, requests)

    def evaluate_owner(
        start_index: int,
        index: int,
        values: Sequence[float],
        receipt_requested: bool,
    ) -> dict[str, object]:
        key = (
            specification.solver.multistart.starts[start_index].identity,
            index,
            resolved.selected_identities,
            tuple(values),
        )
        if not receipt_requested and (block := exact_start_blocks.get(key)) is not None:
            return block
        return direct_owner(start_index, index, values, receipt_requested)

    direct_receipt = FitSpecificationReceipt.from_resolved_fit(resolved)
    direct_receipt_blocks = _receipt_blocks(direct_receipt)
    owner_entries: list[
        tuple[Callable[[int, int, Sequence[float], bool], Mapping[str, object]], int]
    ] = []
    ordered_residual_blocks: list[Mapping[str, object]] = []
    ordered_receipt_rows: list[FitRowReceipt] = []
    ordered_block_sizes: list[int] = []
    equilibrium_all_batches: Sequence[Any] = ()
    equilibrium_owner = None
    equilibrium_batches: Sequence[Any] = ()
    equilibrium_receipt = None
    equilibrium_residual_blocks: tuple[Mapping[str, object], ...] = ()
    full_space = None
    if any(
        isinstance(observation, EquilibriumFitObservation)
        for observation in resolved.fit.observations
    ):
        from ._equilibrium_owner import _equilibrium_owner_batches, _full_space_batches

        (
            equilibrium_all_batches,
            equilibrium_batches,
            equilibrium_owner,
            equilibrium_residual_blocks,
            equilibrium_receipt,
        ) = _equilibrium_owner_batches(resolved)
        if len(equilibrium_batches) == len(resolved.fit.observations):
            full_space = _full_space_batches(resolved, equilibrium_batches)
    equilibrium_receipt_blocks = (
        () if equilibrium_receipt is None else _receipt_blocks(equilibrium_receipt)
    )
    equilibrium_index_by_identity = {
        batch.identity: index for index, batch in enumerate(equilibrium_batches)
    }
    direct_observation_index = 0
    direct_batch_index = 0
    equilibrium_batch_index = 0
    for observation in resolved.fit.observations:
        if isinstance(observation, EquilibriumFitObservation):
            assert equilibrium_owner is not None
            batch = equilibrium_all_batches[equilibrium_batch_index]
            equilibrium_batch_index += 1
            active_index = equilibrium_index_by_identity.get(batch.identity)
            if active_index is None:
                continue
            owner_entries.append((equilibrium_owner, active_index))
            ordered_residual_blocks.append(equilibrium_residual_blocks[active_index])
            rows = equilibrium_receipt_blocks[active_index]
            ordered_receipt_rows.extend(rows)
            ordered_block_sizes.append(len(rows))
            continue
        row_count = len(_rows(resolved.observations[direct_observation_index]))
        direct_observation_index += 1
        for _ in range(row_count):
            owner_entries.append((evaluate_owner, direct_batch_index))
            ordered_residual_blocks.append(residual_blocks[direct_batch_index])
            rows = direct_receipt_blocks[direct_batch_index]
            ordered_receipt_rows.extend(rows)
            ordered_block_sizes.append(len(rows))
            direct_batch_index += 1

    def evaluate_ordered_owner(
        start_index: int,
        index: int,
        values: Sequence[float],
        receipt_requested: bool,
    ) -> Mapping[str, object]:
        owner, local_index = owner_entries[index]
        return owner(start_index, local_index, values, receipt_requested)

    combined_receipt = FitSpecificationReceipt.create(
        resolved.to_mapping(),
        parameters=resolved.parameters,
        coordinates=direct_receipt.coordinates,
        rows=ordered_receipt_rows,
        block_sizes=ordered_block_sizes,
        input_row_count=direct_receipt.input_row_count
        + (0 if equilibrium_receipt is None else equilibrium_receipt.input_row_count),
        dropped_row_count=direct_receipt.dropped_row_count
        + (0 if equilibrium_receipt is None else equilibrium_receipt.dropped_row_count),
    )

    result = _fit_owner_batches(
        resolved,
        evaluate_ordered_owner,
        ordered_residual_blocks,
        combined_receipt,
        full_space,
    )
    if equilibrium_all_batches:
        from ._equilibrium_owner import _validated_owner_receipts

        for start in result.starts:
            _validated_owner_receipts(resolved, equilibrium_all_batches, start)
    return result


__all__ = ("FitPreflightError", "fit")
