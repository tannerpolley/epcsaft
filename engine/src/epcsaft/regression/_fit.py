from __future__ import annotations

import hashlib
import math
from collections.abc import Mapping
from importlib import metadata
from os import PathLike
from pathlib import Path
from typing import Any, SupportsFloat, cast

import pint

from .. import _core
from ..eos import Mixture, Pure, _python_native_artifact_fingerprint
from ..parameters import Parameters
from ._model import SCHEMA_REGISTRY, ParameterFit, TableRow
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


def _target(
    model: Pure | Mixture,
    observation: ResolvedObservation,
    row: TableRow,
    target: str,
) -> dict[str, object]:
    label = SCHEMA_REGISTRY.targets[target].capability_label
    factor = 1.0
    if label.endswith("_mass_density_kg_per_m3"):
        molar_masses = model._model.molar_mass_kg_per_mol
        if len(molar_masses) != 1 or molar_masses[0] is None:
            raise ValueError("pure mass-density fitting requires one declared molar mass")
        factor = molar_masses[0]
        label = label.replace("_mass_density_kg_per_m3", "_molar_density_mol_per_m3")
    observed = row.value(target)
    return {
        "label": label,
        "observed": observed,
        "residual": dict(observation.specification.residuals)[target],
        "scale": observed,
        "value_factor": factor,
    }


def _observation(
    model: Pure | Mixture,
    resolved: ResolvedObservation,
    row: TableRow,
    multiplier: float,
) -> dict[str, object]:
    specification = resolved.specification
    temperature = row.value("temperature_k")
    result: dict[str, object] = {
        "family": specification.family,
        "temperature_k": temperature,
        "multiplier": multiplier,
        "targets": [_target(model, resolved, row, target) for target in specification.targets],
    }
    if specification.family == "pure_saturation":
        result["liquid_density_upper_mol_m3"] = model._native.molar_density_upper_bound(
            temperature, [1.0]
        )
    else:
        result["molar_density_mol_m3"] = row.value("molar_density_mol_per_m3")
        result["composition"] = specification.composition
    return result


def _observations(
    model: Pure | Mixture, resolved: ResolvedFit
) -> list[dict[str, object]]:
    result: list[dict[str, object]] = []
    for observation in resolved.observations:
        if observation.rows is not None:
            result.extend(
                _observation(model, observation, row, observation.specification.multiplier)
                for row in observation.rows
            )
            continue
        if observation.quadrature is None:
            raise RuntimeError("resolved observation has neither rows nor quadrature")
        rule = observation.quadrature.rules[-1]
        target = observation.specification.targets[0]
        for temperature, weight, observed in zip(
            rule.nodes, rule.weights, rule.observed_values, strict=True
        ):
            row = TableRow((("temperature_k", temperature), (target, observed)))
            result.append(
                _observation(
                    model,
                    observation,
                    row,
                    observation.specification.multiplier * math.sqrt(weight),
                )
            )
    return result


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
    native = _core._fit_resolved(
        model._native,
        model._active_coordinates(resolved.active_parameters),
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
        _observations(model, resolved),
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
