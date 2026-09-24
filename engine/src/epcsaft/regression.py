"""Thin public fitting exposure (#49): names map to native indices; the native owners do the science."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence

from . import _core
from ._core import (
    Coordinate, CoordinateSet, FitControls, FitCovariance, FitResult, Observation, ObservationKind,
    ObservationStatus, Prediction, PredictionGradient, ResidualForm, TrialFailure,
)
from .eos import Mixture, Pure, ThermochemistryRecord, _branch
from .equilibrium import Discover, Problem
from .parameters import ParameterError, Parameters

_FAMILIES = {
    "segment_count": _core.ParameterFamily.SegmentCount,
    "sigma_constant": _core.ParameterFamily.SigmaConstant,
    "dispersion_energy_over_k": _core.ParameterFamily.DispersionEnergyOverK,
    "k_ij": _core.ParameterFamily.Kij,
    "l_ij": _core.ParameterFamily.Lij,
    "association_energy_over_k": _core.ParameterFamily.AssociationEnergyOverK,
    "association_volume": _core.ParameterFamily.AssociationVolume,
    "k_ij_reciprocal_temperature_slope": _core.ParameterFamily.KijReciprocalTemperatureSlope,
}
_PAIRS = {"k_ij", "l_ij", "k_ij_reciprocal_temperature_slope"}
_SITE_PAIRS = {"association_energy_over_k", "association_volume"}
_FORMS = {"difference": ResidualForm.Difference, "log_ratio": ResidualForm.LogRatio}
_KINDS = {
    "phase_density": ObservationKind.PhaseDensity,
    "pressure": ObservationKind.EquilibriumPressure,
    "temperature": ObservationKind.EquilibriumTemperature,
    "mole_fraction": ObservationKind.PhaseMoleFraction,
    "phase_fraction": ObservationKind.PhaseFraction,
    "partial_pressure": ObservationKind.PartialPressure,
}
__all__ = ("Coordinate", "CoordinateSet", "FitControls", "FitCovariance", "FitResult", "Observation", "ObservationKind",
           "ObservationStatus", "Prediction", "PredictionGradient", "ResidualForm", "TrialFailure", "coordinate", "density",
           "differentiate", "fit", "fitted_parameters", "observation", "predict")


def _index(ids: Sequence[object], name: object, what: str) -> int:
    try:
        return ids.index(name)
    except ValueError as error:
        raise ParameterError(f"unknown {what} {name}") from error


def _slot(parameters: Parameters, family: str, target: object) -> list[int]:
    """Native slot indices of one target: a component, a mirrored component pair or a site pair."""
    if family not in _PAIRS and family not in _SITE_PAIRS:
        return [_index(parameters.component_ids, target, "component")]
    if not isinstance(target, tuple) or len(target) != 2:
        raise ParameterError(f"{family} needs a pair target")
    ids: Sequence[object] = parameters.site_ids if family in _SITE_PAIRS else parameters.component_ids
    i, j = sorted(_index(ids, item, "site" if family in _SITE_PAIRS else "component") for item in target)
    if i == j and family in _PAIRS:
        raise ParameterError(f"{family} needs two distinct components")
    return [i * len(ids) + i] if i == j else [i * len(ids) + j, j * len(ids) + i]


def coordinate(parameters: Parameters, family: str, *targets: object, origin: float, scale: float,
               bounds: tuple[float, float] | None = None) -> Coordinate:
    """One optimizer coordinate u, theta = origin + scale u in the family's native unit (1, angstrom, K),
    driving every tied target: component ids, (a, b) component pairs, or ((component, site), ...) site pairs."""
    if family not in _FAMILIES:
        raise ParameterError(f"unknown parameter family {family}")
    if not targets:
        raise ParameterError("a coordinate needs at least one target")
    result = Coordinate()
    result.family = _FAMILIES[family]
    result.slots = [_slot(parameters, family, target) for target in targets]
    result.origin, result.scale = origin, scale
    result.lower, result.upper = (-math.inf, math.inf) if bounds is None else bounds
    return result


def _finish(result: Observation, observed: float, form: str, scale: float | None, mass_basis: bool) -> Observation:
    if form not in _FORMS:
        raise ParameterError(f"unknown residual form {form}")
    if (scale is None) != (form == "log_ratio"):
        raise ParameterError("a difference residual needs its scale; a log_ratio residual has none")
    result.observed, result.form, result.mass_basis = observed, _FORMS[form], mass_basis
    result.scale = 1.0 if scale is None else scale
    return result


def density(parameters: Parameters, T: float, P: float, x: Mapping[str, float], *, phase: str, observed: float,
            form: str, mass_basis: bool, scale: float | None = None) -> Observation:
    """A pressure-closed single-phase density (mol/m^3, or kg/m^3 with mass_basis) on the declared branch."""
    unknown = next((name for name in x if name not in parameters.component_ids), None)
    if unknown is not None:
        raise ParameterError(f"unknown species {unknown}")
    result = Observation()
    result.kind, result.temperature, result.pressure = ObservationKind.Density, T, P
    result.composition = [x.get(name, 0.0) for name in parameters.component_ids]
    result.branch = _branch(phase)
    return _finish(result, observed, form, scale, mass_basis)


def observation(parameters: Parameters, kind: str, problem: Problem, *, observed: float, form: str, mass_basis: bool,
                phase: str | None = None, species: Sequence[str] | None = None, scale: float | None = None) -> Observation:
    """One solved-state observation of a declared problem: pressure, temperature, phase_density,
    mole_fraction (summed over `species` in `phase`), partial_pressure (P times that sum, Pa) or
    phase_fraction (mole basis)."""
    if kind not in _KINDS:
        raise ParameterError(f"unknown observation kind {kind}")
    if isinstance(problem.phases, Discover):
        raise ParameterError("a regression observation needs declared phases")
    result = Observation()
    result.kind = _KINDS[kind]
    if phase is None and kind not in ("pressure", "temperature"):
        raise ParameterError(f"{kind} needs its phase")
    result.problem = problem._native(parameters)
    result.phase = 0 if phase is None else _index([item.name for item in problem.phases], phase, "phase")
    result.species = [_index(parameters.component_ids, name, "species") for name in species or ()]
    return _finish(result, observed, form, scale, mass_basis)


def predict(model: Mixture | Pure, item: Observation) -> Prediction:
    """The observation's value and residual from the model's shared owners, with its typed status."""
    return _core.predict(model._model, item)


def differentiate(model: Mixture | Pure, item: Observation) -> PredictionGradient:
    """d(value)/d(theta) per active parameter of the model, in its active order, at a fresh prediction."""
    return _core.differentiate(model._model, item)


def fit(parameters: Parameters, coordinates: Sequence[Coordinate], observations: Sequence[Observation], *,
        weights: Sequence[float] | None = None, controls: FitControls | None = None,
        thermochemistry: ThermochemistryRecord | None = None) -> FitResult:
    """Bounded least squares of sqrt(w_i) r_i over the optimizer coordinates (ALG-FIT)."""
    return _core.fit(parameters._native, thermochemistry, list(coordinates), list(observations), list(weights or ()),
                     FitControls() if controls is None else controls)


def fitted_parameters(start: Parameters, result: FitResult) -> Parameters:
    """The fitted native document under the start's component and site names, for a new Mixture or Pure."""
    return Parameters(result.parameters, start.component_ids, start.site_ids)
