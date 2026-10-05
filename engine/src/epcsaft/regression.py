"""Thin public fitting exposure (#49): names map to native indices; the native owners do the science."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from typing import overload

from . import _core
from ._core import (
    MinimaxRequest, MinimaxResult, MinimaxRow, ObservationEvaluation, HydrationHold, HoldResult, EvaluatedCoordinates, CheckMetric, CheckGroup, CheckResult,
    Coordinate, CoordinateSet, FitControls, FitCovariance, FitResult, Observation, ObservationKind,
    ObservationStatus, Prediction, PredictionGradient, ResidualForm, TrialFailure,
)
from .eos import Mixture, Pure, ThermochemistryRecord, _branch
from .equilibrium import Discover, Problem
from .parameters import ParameterError, Parameters

_FAMILIES = {
    "segment_count": _core.ParameterFamily.SegmentCount,
    "packing_diameter_angstrom": _core.ParameterFamily.PackingDiameter,
    "debye_huckel_diameter_angstrom": _core.ParameterFamily.DebyeHuckelDiameter,
    "born_diameter_angstrom": _core.ParameterFamily.BornDiameter,
    "ion_suppression_coefficient": _core.ParameterFamily.IonSuppressionCoefficient,
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
    "amount_ratio": ObservationKind.PhaseAmountRatio,
}
__all__ = ("MinimaxRequest", "MinimaxResult", "MinimaxRow", "ObservationEvaluation", "Coordinate", "CoordinateSet", "FitControls", "FitCovariance", "FitResult", "Observation", "ObservationKind",
           "ObservationStatus", "Prediction", "PredictionGradient", "ResidualForm", "TrialFailure", "addition_enthalpy",
           "coordinate", "density", "isobaric_heat_capacity", "differentiate", "fit", "fitted_parameters", "observation",
           "predict", "reaction_coordinate", "HydrationHold", "HoldResult", "EvaluatedCoordinates", "CheckMetric",
           "CheckGroup", "CheckResult", "electrolyte", "hydration_hold", "evaluate_coordinates", "evaluate_checks")


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


def reaction_coordinate(name: str, coefficient: str, *, origin: float, scale: float,
                        bounds: tuple[float, float] | None = None) -> Coordinate:
    """Fit a/ln K (1), b (K), c (1) or d (1/K), tied on the named source-basis reaction."""
    if coefficient not in ("a", "b", "c", "d"):
        raise ParameterError(f"unknown reaction coefficient {coefficient}")
    target = _core.ReactionCoordinate()
    target.name, target.coefficient = name, getattr(_core.ReactionCoefficient, coefficient.upper())
    result = Coordinate()
    result.reaction, result.origin, result.scale = target, origin, scale
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
    return _pressure_closed(parameters, ObservationKind.Density, T, P, x, phase, observed, form, scale, mass_basis)


def isobaric_heat_capacity(parameters: Parameters, T: float, P: float, x: Mapping[str, float], *, phase: str,
                          observed: float, form: str, scale: float | None = None) -> Observation:
    """Total molar c_p (J/(mol K)) at fixed T, P, x and branch; thermochemistry is a fixed model input."""
    return _pressure_closed(parameters, ObservationKind.IsobaricHeatCapacity, T, P, x, phase, observed, form, scale, False)


def _pressure_closed(parameters: Parameters, kind: ObservationKind, T: float, P: float, x: Mapping[str, float],
                     phase: str, observed: float, form: str, scale: float | None, mass_basis: bool) -> Observation:
    for name in x:
        _index(parameters.component_ids, name, "species")
    result = Observation()
    result.kind, result.temperature, result.pressure = kind, T, P
    result.composition = [x.get(name, 0.0) for name in parameters.component_ids]
    result.branch = _branch(phase)
    return _finish(result, observed, form, scale, mass_basis)


def observation(parameters: Parameters, kind: str, problem: Problem, *, observed: float, form: str, mass_basis: bool,
                phase: str | None = None, species: Sequence[str] | None = None, scale: float | None = None,
                denominator_phase: str | None = None, denominator_species: Sequence[str] | None = None) -> Observation:
    """One solved-state observation of a declared problem: pressure, temperature, phase_density,
    mole_fraction (summed over `species` in `phase`), partial_pressure (P times that sum, Pa) or
    phase_fraction (mole basis), or amount_ratio (phase-component amounts)."""
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
    if kind == "amount_ratio":
        result.denominator_phase = _index([item.name for item in problem.phases], denominator_phase, "denominator phase")
        result.denominator_species = [_index(parameters.component_ids, name, "species") for name in denominator_species or ()]
    return _finish(result, observed, form, scale, mass_basis)


def addition_enthalpy(parameters: Parameters, initial: Problem, final: Problem, feed: Mapping[str, float], amount: float,
                      T: float, P: float, *, phase: str, observed: float, form: str,
                      scale: float | None = None) -> Observation:
    """Signed total enthalpy increment per mole of added feed, less its enthalpy at (T, P, feed, phase)."""
    if isinstance(initial.phases, Discover) or isinstance(final.phases, Discover):
        raise ParameterError("an addition-enthalpy observation needs declared phases at both endpoints")
    result = _pressure_closed(parameters, ObservationKind.AdditionEnthalpy, T, P, feed, phase,
                              observed, form, scale, False)
    result.problem, result.final_problem = initial._native(parameters), final._native(parameters)
    result.feed_amount = amount
    return result


def predict(model: Mixture | Pure, item: Observation) -> Prediction:
    """The observation's value and residual from the model's shared owners, with its typed status."""
    return _core.predict(model._model, item)


def differentiate(model: Mixture | Pure, item: Observation) -> PredictionGradient:
    """d(value)/d(theta) per active parameter of the model, in its active order, at a fresh prediction."""
    return _core.differentiate(model._model, item)


@overload
def fit(parameters: Parameters, coordinates: Sequence[Coordinate], observations: Sequence[Observation], *,
        weights: Sequence[float] | None = None, controls: FitControls | None = None,
        thermochemistry: ThermochemistryRecord | None = None, holds: Sequence[HydrationHold] | None = None,
        checks: Sequence[CheckGroup] | None = None, baseline: Parameters | None = None, singular_vectors: bool | None = None, minimax: MinimaxRequest | None = None) -> FitResult: ...


@overload
def fit(parameters: Parameters, coordinates: Sequence[Coordinate], observations: Sequence[Observation], *,
        weights: Sequence[float] | None = None, controls: FitControls | None = None,
        thermochemistry: ThermochemistryRecord | None = None, holds: Sequence[HydrationHold] | None = None,
        checks: Sequence[CheckGroup] | None = None, baseline: Parameters | None = None, singular_vectors: bool | None = None,
        starts: Sequence[Sequence[float]], shared_iterations: int | None = None,
        shared_seconds: float | None = None, minimax: MinimaxRequest | None = None) -> list[FitResult]: ...


def fit(parameters: Parameters, coordinates: Sequence[Coordinate], observations: Sequence[Observation], *,
        weights: Sequence[float] | None = None, controls: FitControls | None = None,
        thermochemistry: ThermochemistryRecord | None = None, holds: Sequence[HydrationHold] | None = None,
        checks: Sequence[CheckGroup] | None = None, baseline: Parameters | None = None, singular_vectors: bool | None = None,
        starts: Sequence[Sequence[float]] | None = None, shared_iterations: int | None = None,
        shared_seconds: float | None = None, minimax: MinimaxRequest | None = None) -> FitResult | list[FitResult]:
    """Native least squares, or declared guarded minimax, over the held chart (ALG-FIT)."""
    return _core.fit(parameters._native, thermochemistry, list(coordinates), list(observations), list(weights or ()),
                     FitControls() if controls is None else controls, list(holds or ()), list(checks or ()),
                     None if baseline is None else baseline._native, bool(singular_vectors), None if starts is None else [list(start) for start in starts], shared_iterations, shared_seconds, minimax)


def fitted_parameters(start: Parameters, result: FitResult) -> Parameters:
    """The fitted native document under the start's component and site names, for a new Mixture or Pure."""
    return Parameters(result.parameters, start.component_ids, start.site_ids)


def electrolyte(parameters: Parameters, kind: str, T: float, P: float, x: Mapping[str, float], *,
                ions: tuple[str, str], stoichiometry: tuple[int, int], molality: float,
                species: str | None = None, reference: Mapping[str, float] | None = None,
                observed: float, form: str, scale: float | None = None) -> Observation:
    """Molality activity, osmotic coefficient, single-ion energy (J/mol), or dimensionless permittivity."""
    kinds = {"mean_ionic_activity_coefficient": ObservationKind.MeanIonicActivityCoefficient,
             "osmotic_coefficient": ObservationKind.OsmoticCoefficient, "solvation_gibbs_energy": ObservationKind.SolvationGibbsEnergy,
             "transfer_gibbs_energy": ObservationKind.TransferGibbsEnergy, "bulk_relative_permittivity": ObservationKind.BulkRelativePermittivity}
    if kind not in kinds:
        raise ParameterError(f"unknown electrolyte observation {kind}")
    result = _pressure_closed(parameters, kinds[kind], T, P, x, "liquid", observed, form, scale, False)
    result.ion_indices = (_index(parameters.component_ids, ions[0], "ion"), _index(parameters.component_ids, ions[1], "ion"))
    result.stoichiometry, result.molality = stoichiometry, molality
    result.species = [] if species is None else [_index(parameters.component_ids, species, "ion")]
    if reference is not None:
        for name in reference:
            _index(parameters.component_ids, name, "reference component")
        result.reference_composition = [reference.get(name, 0.0) for name in parameters.component_ids]
    return result


def hydration_hold(reference: Observation, *, bracket: tuple[float, float], energy_tolerance: float,
                   diameter_tolerance: float, maximum_iterations: int) -> HydrationHold:
    """Declare H = RT ln(phi_i at infinite dilution) - target, with a dependent Born diameter."""
    result = HydrationHold()
    result.reference, result.lower, result.upper = reference, *bracket
    result.energy_tolerance, result.diameter_tolerance, result.maximum_iterations = energy_tolerance, diameter_tolerance, maximum_iterations
    return result


def evaluate_coordinates(parameters: Parameters, coordinates: Sequence[Coordinate], optimizer: Sequence[float],
                         holds: Sequence[HydrationHold], thermochemistry: ThermochemistryRecord | None = None) -> EvaluatedCoordinates:
    """Evaluate the native physical chart and exact D at one declared optimizer point."""
    chart = CoordinateSet()
    chart.coordinates = list(coordinates)
    return _core.evaluate_coordinates(parameters._native, chart, list(optimizer), list(holds), thermochemistry)


def evaluate_checks(parameters: Parameters, observations: Sequence[Observation], checks: Sequence[CheckGroup],
                    baseline: Parameters | None = None) -> list[CheckResult]:
    """Native unweighted final metrics on stored parameters, without re-solving holds."""
    return _core.evaluate_checks(parameters._native, list(observations), list(checks), None if baseline is None else baseline._native)
