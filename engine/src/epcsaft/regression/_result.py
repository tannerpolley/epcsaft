from __future__ import annotations

import hashlib
import json
import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from os import PathLike
from pathlib import Path
from typing import Any

from ..parameters import Parameters
from ..units import unit_registry
from ._resolve import ResolvedFit


def _canonical(value: object) -> bytes:
    return json.dumps(
        value,
        ensure_ascii=False,
        allow_nan=False,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("utf-8")


def _finite_or_none(value: float) -> float | None:
    return value if math.isfinite(value) else None


def _finite_values(values: Sequence[float]) -> list[float | None]:
    return [_finite_or_none(value) for value in values]


@dataclass(frozen=True, slots=True)
class FitStart:
    index: int
    initial_physical: tuple[float, ...]
    final_physical: tuple[float, ...]
    final_affine: tuple[float, ...]
    residuals: tuple[float, ...]
    physical_jacobian: tuple[float, ...]
    affine_jacobian: tuple[float, ...]
    singular_values: tuple[float, ...]
    rank: int
    condition_number: float
    active_bounds: tuple[int, ...]
    solver_status: str
    failure_reason: str
    iterations: int
    non_evaluable_count: int
    cost: float
    usable: bool
    basin: int | None

    @classmethod
    def from_native(cls, value: Mapping[str, Any]) -> FitStart:
        diagnostics = value["diagnostics"]
        if not isinstance(diagnostics, Mapping):
            raise TypeError("native fit diagnostics must be a mapping")
        return cls(
            index=int(value["index"]),
            initial_physical=tuple(float(item) for item in value["initial_physical"]),
            final_physical=tuple(float(item) for item in value["final_physical"]),
            final_affine=tuple(float(item) for item in value["final_affine"]),
            residuals=tuple(float(item) for item in value["residuals"]),
            physical_jacobian=tuple(float(item) for item in value["physical_jacobian"]),
            affine_jacobian=tuple(float(item) for item in value["affine_jacobian"]),
            singular_values=tuple(float(item) for item in diagnostics["singular_values"]),
            rank=int(diagnostics["rank"]),
            condition_number=float(diagnostics["condition_number"]),
            active_bounds=tuple(int(item) for item in value["active_bounds"]),
            solver_status=str(value["solver_status"]),
            failure_reason=str(value["failure_reason"]),
            iterations=int(value["iterations"]),
            non_evaluable_count=int(value["non_evaluable_count"]),
            cost=float(value["cost"]),
            usable=bool(value["usable"]),
            basin=None if value["basin"] is None else int(value["basin"]),
        )

    def to_mapping(self, block_sizes: Sequence[int]) -> dict[str, object]:
        contributions: list[float | None] = []
        offset = 0
        for size in block_sizes:
            block = self.residuals[offset : offset + size]
            contributions.append(
                0.5 * sum(value * value for value in block)
                if all(math.isfinite(value) for value in block)
                else None
            )
            offset += size
        return {
            "index": self.index,
            "initial_physical": _finite_values(self.initial_physical),
            "final_physical": _finite_values(self.final_physical),
            "final_affine": _finite_values(self.final_affine),
            "residuals": _finite_values(self.residuals),
            "block_contributions": contributions,
            "physical_jacobian": _finite_values(self.physical_jacobian),
            "affine_jacobian": _finite_values(self.affine_jacobian),
            "singular_values": _finite_values(self.singular_values),
            "rank": self.rank,
            "condition_number": _finite_or_none(self.condition_number),
            "active_bounds": list(self.active_bounds),
            "solver_status": self.solver_status,
            "failure_reason": self.failure_reason,
            "iterations": self.iterations,
            "non_evaluable_count": self.non_evaluable_count,
            "cost": _finite_or_none(self.cost),
            "usable": self.usable,
            "basin": self.basin,
        }


@dataclass(frozen=True, slots=True, init=False)
class FitResult:
    resolved_fit: ResolvedFit
    starts: tuple[FitStart, ...]
    artifact_fingerprint: str
    package_fingerprint: str
    digest: str

    def __init__(
        self,
        resolved_fit: ResolvedFit,
        starts: Sequence[FitStart],
        *,
        artifact_fingerprint: str,
        package_fingerprint: str,
    ) -> None:
        if not isinstance(resolved_fit, ResolvedFit):
            raise TypeError("resolved_fit must be a ResolvedFit")
        frozen_starts = tuple(starts)
        if not frozen_starts or any(not isinstance(item, FitStart) for item in frozen_starts):
            raise ValueError("starts must contain every native start")
        if tuple(item.index for item in frozen_starts) != tuple(range(len(frozen_starts))):
            raise ValueError("starts must retain deterministic native order")
        object.__setattr__(self, "resolved_fit", resolved_fit)
        object.__setattr__(self, "starts", frozen_starts)
        object.__setattr__(self, "artifact_fingerprint", artifact_fingerprint)
        object.__setattr__(self, "package_fingerprint", package_fingerprint)
        object.__setattr__(self, "digest", "")
        object.__setattr__(
            self,
            "digest",
            f"sha256:{hashlib.sha256(_canonical(self._payload())).hexdigest()}",
        )

    @classmethod
    def from_native(
        cls,
        resolved_fit: ResolvedFit,
        native: Mapping[str, Any],
        *,
        artifact_fingerprint: str,
        package_fingerprint: str,
    ) -> FitResult:
        starts = native.get("starts")
        if not isinstance(starts, Sequence) or isinstance(starts, (str, bytes)):
            raise TypeError("native fit starts must be a sequence")
        return cls(
            resolved_fit,
            tuple(FitStart.from_native(item) for item in starts),
            artifact_fingerprint=artifact_fingerprint,
            package_fingerprint=package_fingerprint,
        )

    @property
    def accepted_start(self) -> FitStart | None:
        usable = tuple(start for start in self.starts if start.usable)
        return min(usable, key=lambda start: (start.cost, start.index)) if usable else None

    def _observation_layout(self) -> tuple[tuple[str, ...], tuple[int, ...]]:
        identities: list[str] = []
        block_sizes: list[int] = []
        for observation_index, observation in enumerate(self.resolved_fit.observations):
            states = (
                len(observation.rows)
                if observation.rows is not None
                else len(observation.quadrature.rules[-1].nodes)
                if observation.quadrature is not None
                else 0
            )
            for state_index in range(states):
                block_sizes.append(len(observation.specification.targets))
                identities.extend(
                    f"observation/{observation_index}/state/{state_index}/{target}"
                    for target in observation.specification.targets
                )
        return tuple(identities), tuple(block_sizes)

    def fitted_values(self) -> tuple[tuple[str, float, str], ...]:
        accepted = self.accepted_start
        if accepted is None:
            return ()
        values: list[tuple[str, float, str]] = []
        for coordinate, value in zip(
            self.resolved_fit.coordinates, accepted.final_physical, strict=True
        ):
            values.extend(
                (identity, value, coordinate.unit) for identity in coordinate.tied_identities
            )
        return tuple(values)

    def _payload(self) -> dict[str, object]:
        row_identities, block_sizes = self._observation_layout()
        accepted = self.accepted_start
        return {
            "schema": "epcsaft.fit-result",
            "schema_version": 1,
            "resolved_fit_digest": self.resolved_fit.digest,
            "resolved_fit": self.resolved_fit.to_mapping(),
            "input_parameters": self.resolved_fit.parameters.to_mapping(),
            "artifact_fingerprint": self.artifact_fingerprint,
            "package_fingerprint": self.package_fingerprint,
            "row_identities": list(row_identities),
            "fitted_values": [
                {"identity": identity, "value": _finite_or_none(value), "unit": unit}
                for identity, value, unit in self.fitted_values()
            ],
            "accepted_start": None if accepted is None else accepted.index,
            "starts": [start.to_mapping(block_sizes) for start in self.starts],
        }

    def to_mapping(self) -> dict[str, object]:
        return {**self._payload(), "digest": self.digest}

    def to_json_bytes(self) -> bytes:
        return _canonical(self.to_mapping())

    def to_json(self, path: str | PathLike[str]) -> Path:
        target = Path(path)
        target.write_bytes(self.to_json_bytes() + b"\n")
        return target

    def to_parameters(self) -> Parameters:
        values = {
            identity: unit_registry.Quantity(value, unit)
            for identity, value, unit in self.fitted_values()
        }
        if not values:
            raise RuntimeError("fit has no usable start")
        suffix = self.digest.removeprefix("sha256:")[:16]
        return self.resolved_fit.parameters.with_values(
            values,
            source_id=f"fit-result-{suffix}",
            locator=f"FitResult {self.digest}",
        )


__all__ = ("FitResult", "FitStart")
