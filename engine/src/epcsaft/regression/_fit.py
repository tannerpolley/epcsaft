from __future__ import annotations

import hashlib
import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from importlib import metadata
from os import PathLike
from pathlib import Path
from typing import TYPE_CHECKING, Any, cast

from .. import _core
from .._trials import NON_EVALUABLE_TRIAL_CODES, NonEvaluableTrialCode
from ..eos import (
    OBSERVABLE_TARGETS,
    Mixture,
    ObservableDescriptor,
    Pure,
    _python_native_artifact_fingerprint,
)
from ..parameters import ActiveParameterSet, Parameters
from ._model import EquilibriumFitObservation, ParameterFit, TableRow
from ._resolve import ResolvedFit, ResolvedObservation
from ._result import FitResult, FitRowReceipt, FitSpecificationReceipt

if TYPE_CHECKING:
    from ._equilibrium_owner import _ReactiveBatch


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
                            "scale": 1.0
                            if residuals[target] == "log_ratio"
                            else row.value(target),
                        }
                        for target in specification.targets
                    ],
                }
            )
    return tuple(requests), residual_blocks


def _direct_experiments(
    model: Pure | Mixture,
    resolved: ResolvedFit,
    requests: Sequence[_OwnerRequest],
) -> tuple[Mapping[str, object], ...]:
    active = ActiveParameterSet(resolved.parameters, resolved.selected_identities)
    native_model = model._native_model(active)
    return tuple(
        {
            "kind": "direct_eos",
            "model": native_model,
            "descriptor_identity": request.descriptor.identity,
            "state": dict(request.state),
            "composition": request.composition,
            "selected_rows": [
                request.descriptor.row_identities.index(identity)
                for identity in request.row_identities
            ],
        }
        for request in requests
    )


def _package_fingerprint() -> str:
    record = metadata.distribution("epcsaft").read_text("RECORD")
    if record is None or not record.strip():
        raise RuntimeError("installed epcsaft distribution has no RECORD identity")
    return f"sha256:{hashlib.sha256(record.encode()).hexdigest()}"


def _native_owner_diagnostic(batch: _ReactiveBatch, result: object) -> Mapping[str, object]:
    native = cast(Any, result).equilibrium

    def finite(value: Any) -> float | None:
        number = float(value)
        return number if math.isfinite(number) else None

    return {
        "batch_identity": batch.identity,
        "accepted": bool(native.accepted),
        "solver_status": str(native.solver_status),
        "numerical_status": str(native.numerical_status),
        "physical_status": str(native.physical_status),
        "eos_domain_status": str(native.eos_domain_status),
        "failure_code": str(cast(Any, result).failure_code),
        "failure_reason": str(native.failure_reason),
        "local_variables": [finite(value) for value in native.local_variables],
        "local_lower_bounds": [finite(value) for value in native.local_lower_bounds],
        "local_upper_bounds": [finite(value) for value in native.local_upper_bounds],
        "resolved_policy": [
            {"name": str(name), "value": finite(value)}
            for name, value in native.resolved_policy
        ],
        "certificate": {
            name: finite(getattr(native, name))
            for name in (
                "balance_inf_norm",
                "pressure_relative_inf_norm",
                "reaction_affinity_inf_norm",
                "transfer_chemical_potential_inf_norm",
                "kkt_stationarity_inf_norm",
                "phase_distinctness",
                "minimum_amount_mol",
                "minimum_bound_margin",
                "reduced_hessian_min_eigenvalue",
                "reduced_hessian_scale",
                "reduced_hessian_threshold",
            )
        },
    }


def _fit_owner_batches(
    resolved: ResolvedFit,
    residual_blocks: Sequence[Mapping[str, object]],
    experiments: Sequence[Mapping[str, object]],
    owner_batches: Sequence[_ReactiveBatch | None],
    specification: FitSpecificationReceipt,
) -> FitResult:
    coordinate_indices = {
        identity: index
        for index, coordinate in enumerate(resolved.coordinates)
        for identity in coordinate.tied_identities
    }
    native = _core._fit_resolved(
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
            "optimizer_tolerance": resolved.fit.solver.controls.optimizer_tolerance,
        },
        tuple(sorted(NON_EVALUABLE_TRIAL_CODES)),
        experiments,
        resolved.fit.solver.evaluation_strategy,
    )
    if (failure := native.get("preflight_failure")) is not None:
        record = cast(Mapping[str, object], failure)
        raw_code = str(record["code"])
        code = raw_code if raw_code in NON_EVALUABLE_TRIAL_CODES else "owner_contract_failure"
        raise FitPreflightError(
            cast(NonEvaluableTrialCode, code),
            str(record["start_identity"]),
            str(record["batch_identity"]),
            str(record["diagnostic"]),
        )
    if any(batch is not None for batch in owner_batches):
        from ._equilibrium_owner import _serialize_native_owner_result

        coordinate_by_parameter = {
            identity: index
            for index, coordinate in enumerate(resolved.coordinates)
            for identity in coordinate.tied_identities
        }
        for start in cast(list[dict[str, Any]], native["starts"]):
            owner_results = cast(Sequence[object | None], start.pop("owner_results"))
            start["owner_diagnostics"] = [
                _native_owner_diagnostic(batch, owner_result)
                for batch, owner_result in zip(owner_batches, owner_results, strict=True)
                if batch is not None and owner_result is not None
            ]
            if not bool(start["usable"]):
                continue
            physical = cast(Sequence[float], start["final_physical"])
            active = ActiveParameterSet(
                resolved.parameters,
                resolved.selected_identities,
                values=tuple(
                    physical[coordinate_by_parameter[identity]]
                    for identity in resolved.selected_identities
                ),
            )
            start["owner_receipts_json"] = [
                _serialize_native_owner_result(batch, active, owner_result)
                for batch, owner_result in zip(owner_batches, owner_results, strict=True)
                if batch is not None and owner_result is not None
            ]
    else:
        for start in cast(list[dict[str, Any]], native["starts"]):
            start.pop("owner_results", None)
            start["owner_diagnostics"] = []
    return FitResult.from_native(
        specification,
        native,
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
    direct_experiments = _direct_experiments(model, resolved, requests)

    direct_receipt = FitSpecificationReceipt.from_resolved_fit(resolved)
    direct_receipt_blocks = _receipt_blocks(direct_receipt)
    ordered_experiments: list[Mapping[str, object]] = []
    ordered_owner_batches: list[_ReactiveBatch | None] = []
    ordered_residual_blocks: list[Mapping[str, object]] = []
    ordered_receipt_rows: list[FitRowReceipt] = []
    ordered_block_sizes: list[int] = []
    equilibrium_all_batches: Sequence[_ReactiveBatch] = ()
    equilibrium_experiments: Sequence[Mapping[str, object]] = ()
    equilibrium_batches: Sequence[_ReactiveBatch] = ()
    equilibrium_receipt = None
    equilibrium_residual_blocks: tuple[Mapping[str, object], ...] = ()
    if any(
        isinstance(observation, EquilibriumFitObservation)
        for observation in resolved.fit.observations
    ):
        from ._equilibrium_owner import _equilibrium_owner_batches

        (
            equilibrium_all_batches,
            equilibrium_batches,
            equilibrium_experiments,
            equilibrium_residual_blocks,
            equilibrium_receipt,
        ) = _equilibrium_owner_batches(resolved)
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
            batch = equilibrium_all_batches[equilibrium_batch_index]
            equilibrium_batch_index += 1
            active_index = equilibrium_index_by_identity.get(batch.identity)
            if active_index is None:
                continue
            ordered_experiments.append(equilibrium_experiments[active_index])
            ordered_owner_batches.append(equilibrium_batches[active_index])
            ordered_residual_blocks.append(equilibrium_residual_blocks[active_index])
            rows = equilibrium_receipt_blocks[active_index]
            ordered_receipt_rows.extend(rows)
            ordered_block_sizes.append(len(rows))
            continue
        row_count = len(_rows(resolved.observations[direct_observation_index]))
        direct_observation_index += 1
        for _ in range(row_count):
            ordered_experiments.append(direct_experiments[direct_batch_index])
            ordered_owner_batches.append(None)
            ordered_residual_blocks.append(residual_blocks[direct_batch_index])
            rows = direct_receipt_blocks[direct_batch_index]
            ordered_receipt_rows.extend(rows)
            ordered_block_sizes.append(len(rows))
            direct_batch_index += 1

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
        ordered_residual_blocks,
        ordered_experiments,
        ordered_owner_batches,
        combined_receipt,
    )
    if equilibrium_all_batches:
        from ._equilibrium_owner import _validated_owner_receipts

        for start in result.starts:
            _validated_owner_receipts(resolved, equilibrium_all_batches, start)
    return result


__all__ = ("FitPreflightError", "fit")
