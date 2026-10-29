"""Declare parameter estimation and read native predictions and exact derivatives.

Names map to native component/site indices. All property evaluations,
equilibrium solves, coordinate transforms, hydration constraints and fitting
iterations run in the native owners (ALG-FIT). Temperatures use K, pressures
Pa and feed amounts mol. Observation units are stated by each constructor.
Difference residuals are ``(prediction-observed)/scale``; log-ratio residuals
are natural ``ln(prediction/observed)`` and accept no supplied scale.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from typing import Any, overload

from . import _core
from ._core import (
    MinimaxRequest, MinimaxResult, MinimaxRow, ObservationEvaluation, HydrationHold, HoldResult, EvaluatedCoordinates, CheckMetric, CheckGroup, CheckResult,
    CoordinateTarget, Coordinate, CoordinateSet, FitControls, FitCovariance, FitResult, Observation, ObservationKind,
    ObservationStatus, Prediction, PredictionGradient, ResidualForm, TrialFailure, reaction_log10_reciprocal,
    MeasuredInputKind, MeasuredInput, MeasurementGroup, MeasurementUncertainty, InputReconciliation,
)
from .eos import Mixture, Pure, ThermochemistryRecord, _branch
from .equilibrium import Discover, Problem
from .parameters import ParameterError, Parameters

_FAMILIES = {
    "born_c_dielectric": _core.ParameterFamily.BornCDielectric,
    "ionic_region_relative_permittivity": _core.ParameterFamily.IonicRegionRelativePermittivity,
    "segment_count": _core.ParameterFamily.SegmentCount,
    "packing_diameter": _core.ParameterFamily.PackingDiameter,
    "debye_huckel_diameter": _core.ParameterFamily.DebyeHuckelDiameter,
    "born_diameter": _core.ParameterFamily.BornDiameter,
    "ion_suppression_coefficient": _core.ParameterFamily.IonSuppressionCoefficient,
    "segment_diameter": _core.ParameterFamily.SigmaConstant,
    "dispersion_energy_over_k": _core.ParameterFamily.DispersionEnergyOverK,
    "k_ij": _core.ParameterFamily.Kij,
    "l_ij": _core.ParameterFamily.Lij,
    "association_energy_over_k": _core.ParameterFamily.AssociationEnergyOverK,
    "association_volume": _core.ParameterFamily.AssociationVolume,
    "k_ij_reciprocal_temperature_slope": _core.ParameterFamily.KijReciprocalTemperatureSlope,
    "dielectric_volume": _core.ParameterFamily.DielectricVolume,  # m3/mol; an L/mol coordinate takes coefficient 1e-3
    "dielectric_temperature": _core.ParameterFamily.DielectricTemperature,
    "temperature_law_coefficient": _core.ParameterFamily.TemperatureLawCoefficient,  # target: the coefficient identity
    "cation_pair_mixing": _core.ParameterFamily.CationPairMixing,
    "cation_triple_mixing": _core.ParameterFamily.CationTripleMixing,
}
_PAIRS = {"k_ij", "l_ij", "k_ij_reciprocal_temperature_slope"}
_SITE_PAIRS = {"association_energy_over_k", "association_volume"}
_FORMS = {"difference": ResidualForm.Difference, "log_ratio": ResidualForm.LogRatio}
_KINDS = {
    "equilibrium_isobaric_heat_capacity": ObservationKind.EquilibriumIsobaricHeatCapacity,
    "phase_density": ObservationKind.PhaseDensity,
    "pressure": ObservationKind.EquilibriumPressure,
    "temperature": ObservationKind.EquilibriumTemperature,
    "mole_fraction": ObservationKind.PhaseMoleFraction,
    "phase_fraction": ObservationKind.PhaseFraction,
    "partial_pressure": ObservationKind.PartialPressure,
    "amount_ratio": ObservationKind.PhaseAmountRatio,
}
__all__ = ("MinimaxRequest", "MinimaxResult", "MinimaxRow", "ObservationEvaluation", "CoordinateTarget", "Coordinate", "CoordinateSet", "FitControls", "FitCovariance", "FitResult", "Observation", "ObservationKind",
           "ObservationStatus", "Prediction", "PredictionGradient", "ResidualForm", "TrialFailure", "addition_enthalpy", "saturation_capped_addition_enthalpy",
           "active_parameter", "coordinate", "density", "isobaric_heat_capacity", "differentiate", "fit", "fitted_parameters", "observation",
           "predict", "reaction_coordinate", "reaction_enthalpy_shift", "reaction_log10_reciprocal", "HydrationHold", "HoldResult", "EvaluatedCoordinates", "CheckMetric",
           "CheckGroup", "CheckResult", "electrolyte", "mixture_osmotic_coefficient", "hydration_hold", "evaluate_coordinates", "evaluate_checks",
           "MeasuredInputKind", "MeasuredInput", "MeasurementGroup", "MeasurementUncertainty", "InputReconciliation",
           "measured_input", "measurement_group")


def measured_input(parameters: Parameters, kind: str, *, sigma: float, bounds: tuple[float, float],
                   unit: str, species: str | None = None) -> MeasuredInput:
    """Declare one uncertain physical input for native errors-in-variables estimation.

    Parameters
    ----------
    parameters : epcsaft.parameters.Parameters
        Ordered names used to resolve a mole-fraction species.
    kind : {'temperature', 'pressure', 'mole_fraction', 'molality'}
        Measured input type; respective units are K, Pa, dimensionless and
        mol/kg solvent.
    sigma : float
        Measurement standard deviation in the input's physical unit. Native
        fitting admits finite nonnegative values; zero fixes this input and
        removes its fitted correction.
    bounds : tuple of float
        Lower and upper true-input bounds in the same unit.
    unit : str
        Exact native label: 'K', 'Pa', 'mol/mol' or 'mol/kg' for temperature,
        pressure, mole fraction or molality, respectively.
    species : str, optional
        Required only for ``mole_fraction``; identifies its component.

    Returns
    -------
    input : MeasuredInput
        Native input-error declaration; no reconciliation has been solved.

    Raises
    ------
    ParameterError
        Unknown kind/species, or species supplied for a non-composition input
        or omitted for a mole-fraction input.

    Notes
    -----
    Native fitting checks finite bounds, the physical domain and the measured
    value's inclusion in its bounds. MeasurementUncertainty may declare a joint
    covariance or correlation; without one, the measurement errors are treated
    independently (ALG-FIT).
    """
    kinds = {"temperature": MeasuredInputKind.Temperature, "pressure": MeasuredInputKind.Pressure,
             "mole_fraction": MeasuredInputKind.MoleFraction, "molality": MeasuredInputKind.Molality}
    if kind not in kinds or (kind == "mole_fraction") != (species is not None):
        raise ParameterError("a measured input needs a known kind and only mole_fraction needs a species")
    result = MeasuredInput()
    result.kind, result.sigma, result.lower, result.upper, result.unit = kinds[kind], sigma, *bounds, unit
    if species is not None:
        result.species = _index(parameters.component_ids, species, "measured component")
    return result


def measurement_group(parameters: Parameters, rows: Sequence[int], inputs: Sequence[MeasuredInput], *,
                      dependent_species: Sequence[str] | None = None, feed_total_mol: float | None = None) -> MeasurementGroup:
    """Group observation rows sharing the same uncertain input measurement.

    Parameters
    ----------
    parameters : epcsaft.parameters.Parameters
        Component names used for dependent-species indices.
    rows : sequence of int
        Zero-based observation-row indices sharing the inputs.
    inputs : sequence of MeasuredInput
        Measurement declarations with K, Pa, dimensionless mole-fraction or
        mol/kg solvent units as appropriate.
    dependent_species : sequence of str, optional
        Explicit dependent feed species for composition closure.
    feed_total_mol : float, optional
        Fixed total feed amount in mol when required by the measurement chart.

    Returns
    -------
    group : MeasurementGroup
        Native row/input group; independent input errors are not fitted here.

    Raises
    ------
    ParameterError
        A dependent-species identifier is unknown.

    Notes
    -----
    Native fitting validates row membership, closure and shared-measurement
    requirements before solving the input-error objective (ALG-FIT).
    """
    result = MeasurementGroup()
    result.rows, result.inputs = list(rows), list(inputs)
    if dependent_species is not None:
        result.dependent_species = [_index(parameters.component_ids, name, "dependent component") for name in dependent_species]
    if feed_total_mol is not None:
        result.feed_total_mol = feed_total_mol
    return result


def _index(ids: Sequence[object], name: object, what: str) -> int:
    try:
        return ids.index(name)
    except ValueError as error:
        raise ParameterError(f"unknown {what} {name}; expected one of {', '.join(map(str, ids))}") from error


def _slot(parameters: Parameters, family: str, target: object) -> list[int]:
    """Native slot indices of one target: a component, a mirrored component pair or a site pair."""
    if family not in _FAMILIES:
        raise ParameterError(f"unknown parameter family {family}; expected one of {', '.join(_FAMILIES)}")
    if family in ("born_c_dielectric", "ionic_region_relative_permittivity") or family.startswith("cation_"):
        if target is not None:
            raise ParameterError(f"{family} is global and takes no component target")
        return [0]
    if family == "temperature_law_coefficient":
        return [_index(parameters._native.temperature_law_coefficient_identities, target, "temperature-law coefficient")]
    if family not in _PAIRS and family not in _SITE_PAIRS:
        return [_index(parameters.component_ids, target, "component")]
    if not isinstance(target, tuple) or len(target) != 2:
        raise ParameterError(f"{family} needs a pair target")
    ids: Sequence[object] = parameters.site_ids if family in _SITE_PAIRS else parameters.component_ids
    i, j = sorted(_index(ids, item, "site" if family in _SITE_PAIRS else "component") for item in target)
    if i == j and family in _PAIRS:
        raise ParameterError(f"{family} needs two distinct components")
    return [i * len(ids) + i] if i == j else [i * len(ids) + j, j * len(ids) + i]


def active_parameter(parameters: Parameters, family: str, target: object) -> _core.ActiveParameter:
    """Resolve a named physical parameter into exact native derivative slots.

    Parameters
    ----------
    parameters : epcsaft.parameters.Parameters
        Component/site order and temperature-law coefficient identities.
    family : str
        Native parameter-family name, e.g. ``segment_diameter`` or ``k_ij``.
    target : object
        Component identifier, unordered component pair, association-site pair,
        temperature-law coefficient identity, or ``None`` for a global scalar.

    Returns
    -------
    parameter : ActiveParameter
        One physical derivative coordinate owning the selected native slots.
        Native units include Å for diameters, K for energy/k, m³/mol for
        dielectric volume and dimensionless k_ij/l_ij (ALG-PARAMETERS).

    Raises
    ------
    ParameterError
        Family/target is unknown, a global target is supplied, or an unlike
        component pair has equal endpoints.

    Notes
    -----
    Unlike pairs own both mirrored entries. Same-site association may own its
    one diagonal slot. No logarithmic coordinate transform is implied.
    """
    slots = _slot(parameters, family, target)
    return _core.ActiveParameter(_FAMILIES[family], slots)


def coordinate(parameters: Parameters, targets: Sequence[tuple[str, object, float]], *, origin: float, scale: float,
               bounds: tuple[float, float] | None = None) -> Coordinate:
    """Declare an affine optimizer coordinate tied to named physical parameters.

    Parameters
    ----------
    parameters : epcsaft.parameters.Parameters
        Names and native slots to resolve.
    targets : sequence of tuple
        (Family, component/pair/coefficient identity or None, coefficient).
        A target's value is coefficient times physical coordinate q; the
        coefficient converts q units into the native target-family unit.
    origin : float
        Physical q origin, in the caller-declared coordinate unit.
    scale : float
        Physical q increment per dimensionless optimizer increment.
    bounds : tuple of float, optional
        Inclusive q bounds in the same physical unit; absent means unbounded.

    Returns
    -------
    coordinate : Coordinate
        Declaration with ``q = origin + scale*v`` for dimensionless optimizer v.
        Native target values are coefficient*q.

    Raises
    ------
    ParameterError
        A family or target cannot be resolved by ``active_parameter`` rules.

    Notes
    -----
    Native chart admission checks scales, bounds and conflicting targets.
    For a source dielectric-volume coordinate in L/mol, coefficient 1e-3
    converts to native m³/mol (ALG-PARAMETERS and ALG-FIT).
    """
    result = Coordinate()
    rows = []
    for family, name, coefficient in targets:
        target = CoordinateTarget()
        target.indices, target.family, target.coefficient = _slot(parameters, family, name), _FAMILIES[family], coefficient
        rows.append(target)
    result.targets = rows
    result.origin, result.scale = origin, scale
    result.lower, result.upper = (-math.inf, math.inf) if bounds is None else bounds
    return result


def reaction_coordinate(name: str, coefficient: str, *, origin: float, scale: float,
                        bounds: tuple[float, float] | None = None) -> Coordinate:
    """Declare a fitted coefficient of a named source-basis reaction law.

    Parameters
    ----------
    name : str
        Reaction name shared by the fitted observation declarations.
    coefficient : {'a', 'b', 'c', 'd'}
        Coefficient of natural ln K(T) = a + b/T + c*ln(T/T_ref) + d*T.
        Units are dimensionless, K, dimensionless and 1/K respectively.
    origin : float
        Physical coefficient origin in its stated unit.
    scale : float
        Coefficient increment per dimensionless optimizer increment.
    bounds : tuple of float, optional
        Inclusive physical coefficient bounds; default is unbounded.

    Returns
    -------
    coordinate : Coordinate
        Named reaction coordinate with the affine optimizer chart.

    Raises
    ------
    ParameterError
        Coefficient is not a, b, c or d.

    Notes
    -----
    For a scalar reaction K, the a coordinate fits natural ln K. Native
    admission checks named-reaction consistency and source reference (ALG-FIT).

    Examples
    --------
    >>> from epcsaft.regression import reaction_coordinate
    >>> coordinate = reaction_coordinate('reaction', 'b', origin=100.0, scale=10.0)
    >>> coordinate.origin, coordinate.scale
    (100.0, 10.0)
    """
    if coefficient not in ("a", "b", "c", "d"):
        raise ParameterError(f"unknown reaction coefficient {coefficient}; expected a, b, c or d")
    target = _core.ReactionCoordinate()
    target.name, target.coefficient = name, getattr(_core.ReactionCoefficient, coefficient.upper())
    result = Coordinate()
    result.reaction, result.origin, result.scale = target, origin, scale
    result.lower, result.upper = (-math.inf, math.inf) if bounds is None else bounds
    return result


def reaction_enthalpy_shift(name: str, *, pivot_temperature: float, origin: float, scale: float,
                            bounds: tuple[float, float]) -> Coordinate:
    """Declare a reaction-enthalpy shift preserving ln K at a fixed pivot.

    Parameters
    ----------
    name : str
        Named source-basis reaction to shift.
    pivot_temperature : float
        Required pivot temperature in K at which supplied natural ln K is kept.
    origin : float
        Enthalpy-shift origin in J/mol of reaction as written.
    scale : float
        J/mol shift per dimensionless optimizer increment.
    bounds : tuple of float
        Inclusive enthalpy-shift bounds in J/mol.

    Returns
    -------
    coordinate : Coordinate
        Native EnthalpyShift declaration with the pivot recorded in K.

    Notes
    -----
    The native reaction-coordinate owner performs the coupled coefficient
    change; this function does not evaluate or refit a reaction law (ALG-FIT).

    Examples
    --------
    >>> from epcsaft.regression import reaction_enthalpy_shift
    >>> shift = reaction_enthalpy_shift('reaction', pivot_temperature=300.0,
    ...                                origin=0.0, scale=1000.0, bounds=(-2000.0, 2000.0))
    >>> shift.reaction.pivot_temperature
    300.0
    """
    result = reaction_coordinate(name, "a", origin=origin, scale=scale, bounds=bounds)
    assert result.reaction is not None
    result.reaction.coefficient = _core.ReactionCoefficient.EnthalpyShift
    result.reaction.pivot_temperature = pivot_temperature
    return result


def _finish(result: Observation, observed: float, form: str, scale: float | None, mass_basis: bool) -> Observation:
    if form not in _FORMS:
        raise ParameterError(f"unknown residual form {form}; expected one of {', '.join(_FORMS)}")
    if (scale is None) != (form == "log_ratio"):
        raise ParameterError("a difference residual needs its scale; a log_ratio residual has none")
    result.observed, result.form, result.mass_basis = observed, _FORMS[form], mass_basis
    result.scale = 1.0 if scale is None else scale
    return result


def density(parameters: Parameters, T: float, P: float, x: Mapping[str, float], *, phase: str, observed: float,
            form: str, mass_basis: bool, scale: float | None = None) -> Observation:
    """Declare a pressure-closed single-phase density observation.

    Molar density in mol/m³, or mass density in kg/m³ with ``mass_basis``.

    Parameters
    ----------
    parameters : epcsaft.parameters.Parameters
        Component order used to resolve named compositions.
    T : float
        Absolute temperature in K.
    P : float
        Pressure in Pa.
    x : mapping of str to float
        Dimensionless mole fractions; missing model species are zero.
    phase : {'liquid', 'vapor'}
        Required pressure-density branch.
    observed : float
        Measured value in the units stated above.
    form : {'difference', 'log_ratio'}
        Difference divided by scale, or natural ln(prediction/observed).
    scale : float, optional
        Required positive scale in observed units for a difference residual;
        must be absent for a log-ratio residual.
    mass_basis : bool
        False for molar density, true for mass density.

    Returns
    -------
    observation : Observation
        Native physical observation declaration; prediction is deferred.

    Raises
    ------
    ParameterError
        A species/phase is unknown, residual form is unsupported or scale
        presence conflicts with the residual form.

    Notes
    -----
    Native prediction owns density closure, physical admission and typed
    unavailability. Total heat capacity requires fixed ideal thermochemistry;
    this declaration supplies no missing thermal reference (ALG-FIT).
    """
    return _pressure_closed(parameters, ObservationKind.Density, T, P, x, phase, observed, form, scale, mass_basis)


def isobaric_heat_capacity(parameters: Parameters, T: float, P: float, x: Mapping[str, float], *, phase: str,
                          observed: float, form: str, scale: float | None = None) -> Observation:
    """Declare a fixed-composition total molar heat-capacity observation.

    Total molar isobaric heat capacity in J/(mol K).

    Parameters
    ----------
    parameters : epcsaft.parameters.Parameters
        Component order used to resolve named compositions.
    T : float
        Absolute temperature in K.
    P : float
        Pressure in Pa.
    x : mapping of str to float
        Dimensionless mole fractions; missing model species are zero.
    phase : {'liquid', 'vapor'}
        Required pressure-density branch.
    observed : float
        Measured value in the units stated above.
    form : {'difference', 'log_ratio'}
        Difference divided by scale, or natural ln(prediction/observed).
    scale : float, optional
        Required positive scale in observed units for a difference residual;
        must be absent for a log-ratio residual.

    Returns
    -------
    observation : Observation
        Native physical observation declaration; prediction is deferred.

    Raises
    ------
    ParameterError
        A species/phase is unknown, residual form is unsupported or scale
        presence conflicts with the residual form.

    Notes
    -----
    Native prediction owns density closure, physical admission and typed
    unavailability. Total heat capacity requires fixed ideal thermochemistry;
    this declaration supplies no missing thermal reference (ALG-FIT).
    """
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


def observation(parameters: Parameters, kind: str, problem: Problem, *, observed: float, form: str, mass_basis: bool | None = None,
                phase: str | None = None, species: Sequence[str] | None = None, scale: float | None = None,
                denominator_phase: str | None = None, denominator_species: Sequence[str] | None = None) -> Observation:
    """Declare an observable of a solved equilibrium problem.

    Parameters
    ----------
    parameters : epcsaft.parameters.Parameters
        Model component order for species and problem conversion.
    kind : str
        ``pressure`` (Pa), ``temperature`` (K), ``phase_density`` (mol/m³ or
        kg/m³), ``mole_fraction`` (dimensionless sum over selected species),
        ``partial_pressure`` (Pa, P times that sum), ``phase_fraction``
        (dimensionless molar fraction), ``amount_ratio`` (dimensionless ratio
        of selected phase-component mol), or
        ``equilibrium_isobaric_heat_capacity`` (J/(mol-feed K), driven T/P).
    problem : Problem
        Declared phases and physical constraints; discovery is refused.
    observed : float
        Measured quantity in the unit selected by ``kind``.
    form : {'difference', 'log_ratio'}
        Difference divided by scale or natural ln(prediction/observed).
    mass_basis : bool, optional
        Required for phase density: true selects kg/m³, false mol/m³.
    phase : str, optional
        Required phase name except for pressure, temperature and equilibrium
        heat capacity.
    species : sequence of str, optional
        Selected numerator components for composition, pressure or amount rows.
    scale : float, optional
        Required positive scale in observed units for difference residuals;
        absent for log-ratio residuals.
    denominator_phase : str, optional
        Required denominator phase for an amount ratio.
    denominator_species : sequence of str, optional
        Selected denominator components for an amount ratio.

    Returns
    -------
    observation : Observation
        Native problem and requested quantity with resolved component indices.

    Raises
    ------
    ParameterError
        Unknown kind/phase/species, discovery, missing required phase/density
        basis, or inconsistent residual/scale declarations.

    Notes
    -----
    The native owner checks structural observability and availability after
    solving the declared problem. No equilibrium solve occurs here (ALG-FIT).
    """
    if kind not in _KINDS:
        raise ParameterError(f"unknown observation kind {kind}; expected one of {', '.join(_KINDS)}")
    if isinstance(problem.phases, Discover):
        raise ParameterError("a regression observation needs declared phases")
    result = Observation()
    result.kind = _KINDS[kind]
    if phase is None and kind not in ("pressure", "temperature", "equilibrium_isobaric_heat_capacity"):
        raise ParameterError(f"{kind} needs its phase")
    if mass_basis is None and kind == "phase_density":  # only a density has a basis
        raise ParameterError("phase_density needs mass_basis: True for kg/m^3, False for mol/m^3")
    result.problem = problem._native(parameters)
    result.phase = 0 if phase is None else _index([item.name for item in problem.phases], phase, "phase")
    result.species = [_index(parameters.component_ids, name, "species") for name in species or ()]
    if kind == "amount_ratio":
        result.denominator_phase = _index([item.name for item in problem.phases], denominator_phase, "denominator phase")
        result.denominator_species = [_index(parameters.component_ids, name, "species") for name in denominator_species or ()]
    return _finish(result, observed, form, scale, False if mass_basis is None else mass_basis)


def saturation_capped_addition_enthalpy(parameters: Parameters, initial: Problem, bubble: Problem,
        feed: Mapping[str, float], *, basis_amount: float, loading_cap: float,
        loading_bounds: tuple[float, float], loading_grid: Sequence[float], T: float, P: float,
        phase: str, observed: float, form: str, scale: float | None = None) -> Observation:
    """Declare addition heat capped by the native model's saturation loading.

    Parameters
    ----------
    parameters : Parameters
        Component and fitted-parameter declaration shared by the states.
    initial, bubble : Problem
        Zero-addition state at fixed T/P and declared free-pressure bubble state.
    feed : mapping of str to float
        Neutral pure added feed, with mole fraction one.
    basis_amount : float
        Fixed analytical solvent basis in mol.
    loading_cap : float
        Fixed source loading cap in mol added feed per mol analytical basis.
    loading_bounds : tuple of float
        Positive lower and upper loading bounds in mol/mol.
    loading_grid : sequence of float
        Positive increasing native search nodes, including both bounds.
    T, P : float
        Common temperature in K and physical pressure in Pa.
    phase : str
        Declared feed density branch, liquid or vapor.
    observed : float
        Signed measured heat in J/mol added feed; exothermic heat is negative.
    form : str
        Difference residual with its fixed source scale.
    scale : float, optional
        Fixed positive difference scale in J/mol added feed.

    Returns
    -------
    observation : Observation
        Native declaration; no solve or thermodynamic calculation is done here.

    Raises
    ------
    ParameterError
        The feed, state, grid, branch or residual declaration is inadmissible.

    Notes
    -----
    The grid and successive native seeds define a finite search. Unsampled
    tangencies, paired crossings or multiple roots within a cell, turns between
    nodes and undiscovered branches can escape detection. Availability means
    unique within the declared search. At an exact cap switch the classical
    gradient is unavailable. Native results retain search and action causes.
    """
    result = addition_enthalpy(parameters, initial, bubble, feed, basis_amount, T, P,
                              phase=phase, observed=observed, form=form, scale=scale)
    result.kind = ObservationKind.SaturationCappedAdditionEnthalpy
    result.basis_amount, result.loading_cap = basis_amount, loading_cap
    result.loading_bounds, result.loading_grid = loading_bounds, list(loading_grid)
    return result


def addition_enthalpy(parameters: Parameters, initial: Problem, final: Problem, feed: Mapping[str, float], amount: float,
                      T: float, P: float, *, phase: str, observed: float, form: str,
                      scale: float | None = None) -> Observation:
    """Declare signed enthalpy of addition between two equilibrium endpoints.

    Parameters
    ----------
    parameters : epcsaft.parameters.Parameters
        Component order used by endpoints and added-feed composition.
    initial, final : Problem
        Declared endpoint problems with total enthalpies in J; discovery is
        refused at either endpoint.
    feed : mapping of str to float
        Dimensionless mole fractions of the added feed.
    amount : float
        Positive total added-feed amount in mol.
    T : float
        Temperature of the added feed in K.
    P : float
        Pressure of the added feed in Pa.
    phase : {'liquid', 'vapor'}
        Branch for the added-feed enthalpy.
    observed : float
        Signed molar addition enthalpy in J/mol of added feed.
    form : {'difference', 'log_ratio'}
        Difference divided by scale or natural ln(prediction/observed).
    scale : float, optional
        Difference scale in J/mol; absent for log-ratio residuals.

    Returns
    -------
    observation : Observation
        Native two-endpoint declaration; no endpoint solve is done here.

    Raises
    ------
    ParameterError
        Discovery is requested, a species/branch is unknown, or residual
        form/scale presence is inconsistent.

    Notes
    -----
    Native prediction computes ``(H_final-H_initial)/amount-h_added_feed``
    with total endpoint H in J and added-feed molar h in J/mol (ALG-FIT).
    """
    if isinstance(initial.phases, Discover) or isinstance(final.phases, Discover):
        raise ParameterError("an addition-enthalpy observation needs declared phases at both endpoints")
    result = _pressure_closed(parameters, ObservationKind.AdditionEnthalpy, T, P, feed, phase,
                              observed, form, scale, False)
    result.problem, result.final_problem = initial._native(parameters), final._native(parameters)
    result.feed_amount = amount
    return result


def predict(model: Mixture | Pure, item: Observation) -> Prediction:
    """Evaluate one observation using the model's shared native calculations.

    Parameters
    ----------
    model : Mixture or Pure
        Model with any required thermochemistry and active-parameter declaration.
    item : Observation
        Physical observation, including its observed value, units and residual.

    Returns
    -------
    prediction : Prediction
        Value in observation units and dimensionless residual, with typed
        unavailability, state/equilibrium/reference diagnostics as applicable.

    Raises
    ------
    ParameterError
        Native physical declaration admission fails.

    Notes
    -----
    Inspect the returned prediction cause before consuming its value. Failed
    trials are not observations of zero (ALG-FIT).
    """
    return _core.predict(model._model, item)


def differentiate(model: Mixture | Pure | Parameters, item: Observation | Sequence[Coordinate],
                  observations: Sequence[Observation] | None = None, *, optimizer: Sequence[float] | None = None,
                  thermochemistry: ThermochemistryRecord | None = None,
                  holds: Sequence[HydrationHold] | None = None,
                  solver_options: Mapping[str, float] | None = None) -> Any:
    """Evaluate exact native prediction or coordinate-aware residual derivatives.

    Parameters
    ----------
    model : Mixture, Pure or epcsaft.parameters.Parameters
        Fluid model for one-observation differentiation, or a parameter set
        for coordinate-aware evaluation.
    item : Observation or sequence of Coordinate
        One observation for a model; coordinates for a parameter set.
    observations : sequence of Observation, optional
        Required with a parameter set, in row order.
    optimizer : sequence of float, optional
        Required dimensionless optimizer point for coordinate-aware evaluation.
    thermochemistry : ThermochemistryRecord, optional
        Fixed ideal-gas molar thermochemistry for coordinate evaluation;
        temperature intervals K and standard pressure Pa.
    holds : sequence of HydrationHold, optional
        Dependent Born-diameter constraints for coordinate evaluation.
    solver_options : mapping of str to float, optional
        Native tolerances and work budgets for coordinate evaluation.

    Returns
    -------
    evaluation : PredictionGradient or ObservationEvaluation
        For a model, prediction derivatives in active-parameter order, with
        observation units per native parameter unit; inspect status/message for
        unavailable derivatives. For a parameter set, gradient is row-major in
        observation order, with physical columns ordered as independent EOS
        targets, held Born diameters, then reaction coordinates. Optimizer columns
        follow coordinate order; optimizer_gradient contains prediction columns
        in observation units, and residual_optimizer_gradient contains
        dimensionless residual columns at unit observation weights. This path
        returns only a completed evaluation; unavailable rows raise ParameterError.

    Raises
    ------
    ParameterError
        Parameters lack coordinates/observations/optimizer, a model receives
        multiple observations, native physical admission fails, or a requested
        coordinate-aware prediction, chart, hold or derivative is unavailable.

    Notes
    -----
    Thermochemistry, holds and solver_options are consumed only by the parameter
    set path. Derivatives follow the same native equations and solved branches
    as prediction (ALG-FIT and ALG-ACTIONS).
    """
    if isinstance(model, Parameters):
        if isinstance(item, Observation) or observations is None or optimizer is None:
            raise ParameterError("coordinate differentiation needs coordinates and observations")
        return _core.evaluate_observations(model._native, list(item), list(optimizer), list(holds or ()),
                                           list(observations), thermochemistry, solver_options=solver_options)
    if not isinstance(item, Observation):
        raise ParameterError("model differentiation needs one observation")
    return _core.differentiate(model._model, item)


@overload
def fit(parameters: Parameters, coordinates: Sequence[Coordinate], observations: Sequence[Observation], *,
        weights: Sequence[float] | None = None, controls: FitControls,
        thermochemistry: ThermochemistryRecord | None = None, holds: Sequence[HydrationHold] | None = None,
        checks: Sequence[CheckGroup] | None = None, baseline: Parameters | None = None, singular_vectors: bool | None = None,
        minimax: MinimaxRequest | None = None, uncertainty: MeasurementUncertainty | None = None) -> FitResult: ...


@overload
def fit(parameters: Parameters, coordinates: Sequence[Coordinate], observations: Sequence[Observation], *,
        weights: Sequence[float] | None = None, controls: FitControls,
        thermochemistry: ThermochemistryRecord | None = None, holds: Sequence[HydrationHold] | None = None,
        checks: Sequence[CheckGroup] | None = None, baseline: Parameters | None = None, singular_vectors: bool | None = None,
        starts: Sequence[Sequence[float]], shared_iterations: int | None = None,
        shared_seconds: float | None = None, minimax: MinimaxRequest | None = None,
        uncertainty: MeasurementUncertainty | None = None) -> list[FitResult]: ...


def fit(parameters: Parameters, coordinates: Sequence[Coordinate], observations: Sequence[Observation], *,
        weights: Sequence[float] | None = None, controls: FitControls | None = None,
        thermochemistry: ThermochemistryRecord | None = None, holds: Sequence[HydrationHold] | None = None,
        checks: Sequence[CheckGroup] | None = None, baseline: Parameters | None = None, singular_vectors: bool | None = None,
        starts: Sequence[Sequence[float]] | None = None, shared_iterations: int | None = None,
        shared_seconds: float | None = None, minimax: MinimaxRequest | None = None,
        uncertainty: MeasurementUncertainty | None = None) -> FitResult | list[FitResult]:
    """Fit physical parameters with native least squares or declared guarded minimax.

    Parameters
    ----------
    parameters : epcsaft.parameters.Parameters
        Starting physical parameters in native units.
    coordinates : sequence of Coordinate
        Affine physical-to-dimensionless-optimizer chart and bounds.
    observations : sequence of Observation
        Ordered physical observations and residual definitions.
    weights : sequence of float, optional
        Nonnegative dimensionless least-squares weights; residuals are multiplied
        by their square roots. Omitted/empty means unit weights.
    controls : FitControls
        Required iteration and wall-time (s) budgets and stopping controls.
    thermochemistry : ThermochemistryRecord, optional
        Fixed molar ideal-gas thermal record, interval bounds K and standard Pa.
    holds : sequence of HydrationHold, optional
        Dependent physical constraints solved natively within the chart.
    checks : sequence of CheckGroup, optional
        Declared final physical metrics and limits.
    baseline : epcsaft.parameters.Parameters, optional
        Comparison parameter set for final check groups.
    singular_vectors : bool, optional
        Request native singular vectors; default false.
    starts : sequence of sequence of float, optional
        Dimensionless optimizer starts, attempted in supplied order while shared
        budgets remain. Reaction-enthalpy-shift coordinates are refused here.
    shared_iterations : int, optional
        Positive shared iteration budget; defaults to controls.maximum_iterations.
    shared_seconds : float, optional
        Positive shared wall-time budget in s; defaults to
        controls.maximum_elapsed_time_seconds.
    minimax : MinimaxRequest, optional
        Explicit guarded minimax request instead of ordinary least squares.
    uncertainty : MeasurementUncertainty, optional
        Input/output measurement-error declarations for native reconciliation.

    Returns
    -------
    result : FitResult or list of FitResult
        Fitted native parameters, physical/optimizer coordinates, costs,
        diagnostics and any admitted local covariance. Supplied starts return
        one entry per attempted start, in input order; exhausted shared budgets
        leave a shorter list. Inspect stopping and usability before use.

    Raises
    ------
    ParameterError
        The required controls argument is omitted, solver options or shared
        budgets are invalid, or starts are combined with a reaction-enthalpy shift.

    Notes
    -----
    Native invalid fit declarations can return an unusable result with their
    message; unsupported uncertainty with minimax also returns a refusal.
    All fitting iterations, derivative columns, constraint solves and trial
    failures belong to ALG-FIT. A fitted parameter set is calibration evidence,
    not independent physical validation. Local covariance requires the native
    rank and error-model assumptions; unavailable covariance is not zero.
    """
    if controls is None:
        raise ParameterError("fit requires controls: the iteration and time budgets are scientific inputs (ALG-FIT)")
    return _core.fit(parameters._native, thermochemistry, list(coordinates), list(observations), list(weights or ()),
                     controls, list(holds or ()), list(checks or ()),
                     None if baseline is None else baseline._native, bool(singular_vectors), None if starts is None else [list(start) for start in starts], shared_iterations, shared_seconds, minimax, uncertainty)


def fitted_parameters(start: Parameters, result: FitResult) -> Parameters:
    """Associate fitted native values with the starting component and site names.

    Parameters
    ----------
    start : epcsaft.parameters.Parameters
        Source names and, when available, retained unit-bearing document.
    result : FitResult
        Native fit result whose ``parameters`` are to be retained.

    Returns
    -------
    parameters : epcsaft.parameters.Parameters
        Fitted native values with the starting names. A representable source
        document is updated for export; otherwise its precise refusal is stored.

    Notes
    -----
    This function does not require a usable fit result; inspect fit diagnostics
    first. A source document that cannot express the fitted values does not
    invalidate the native parameter set. Subsequent ASL export raises the stored
    cause. No absent source document is invented.
    """
    fitted = Parameters(result.parameters, start.component_ids, start.site_ids)
    if isinstance(start._document, _core.ParameterDocument):
        try:
            fitted._document = _core.fitted_document(start._document, result.parameters)
        except ValueError as refusal:  # the fit stays usable; the ASL record builders raise this cause
            fitted._document = f"the start document cannot represent the fit: {refusal}"
    return fitted


def electrolyte(parameters: Parameters, kind: str, T: float, P: float, x: Mapping[str, float], *,
                ions: tuple[str, str], stoichiometry: tuple[int, int], molality: float,
                species: str | None = None, reference: Mapping[str, float] | None = None,
                observed: float, form: str, scale: float | None = None) -> Observation:
    """Declare a liquid electrolyte or single-ion-energy observation.

    Parameters
    ----------
    parameters : epcsaft.parameters.Parameters
        Names and order for compositions and ions.
    kind : str
        ``mean_ionic_activity_coefficient`` (dimensionless molality scale),
        ``osmotic_coefficient`` (dimensionless), ``solvation_gibbs_energy``
        or ``transfer_gibbs_energy`` (J/mol of the selected ion), or
        ``bulk_relative_permittivity`` (dimensionless).
    T : float
        Temperature in K.
    P : float
        Pressure in Pa.
    x : mapping of str to float
        Dimensionless mole fractions: full composition for bulk permittivity,
        solvent-only composition for the ionic reference/activity path;
        missing components are zero.
    ions : tuple of str
        Two formula-ion identifiers, in declared ion order.
    stoichiometry : tuple of int
        Positive dimensionless formula counts in ion order.
    molality : float
        Formula molality in mol/kg solvent.
    species : str, optional
        Selected ion for a single-ion energy observation.
    reference : mapping of str to float, optional
        Dimensionless second solvent-reference composition for transfer energy;
        for bulk permittivity, a second full composition requests the relative
        permittivity at x divided by the relative permittivity at reference.
    observed : float
        Value in the unit selected by kind.
    form : {'difference', 'log_ratio'}
        Difference divided by scale or natural ln(prediction/observed).
    scale : float, optional
        Required scale in observed units for difference residuals; absent for
        log-ratio residuals.

    Returns
    -------
    observation : Observation
        Liquid-branch native declaration with named indices resolved.

    Raises
    ------
    ParameterError
        Kind/species/reference names are unknown, or residual/scale declarations
        conflict. Native physical admission is deferred to prediction/fitting.

    Notes
    -----
    Single-ion energies use RT times natural ln(phi) at infinite dilution or
    its transfer difference. They retain the declared single-ion convention
    (ALG-REFERENCE and ALG-FIT).
    """
    kinds = {"mean_ionic_activity_coefficient": ObservationKind.MeanIonicActivityCoefficient,
             "osmotic_coefficient": ObservationKind.OsmoticCoefficient, "solvation_gibbs_energy": ObservationKind.SolvationGibbsEnergy,
             "transfer_gibbs_energy": ObservationKind.TransferGibbsEnergy, "bulk_relative_permittivity": ObservationKind.BulkRelativePermittivity}
    if kind not in kinds:
        raise ParameterError(f"unknown electrolyte observation {kind}; expected one of {', '.join(kinds)}")
    result = _pressure_closed(parameters, kinds[kind], T, P, x, "liquid", observed, form, scale, False)
    result.ion_indices = (_index(parameters.component_ids, ions[0], "ion"), _index(parameters.component_ids, ions[1], "ion"))
    result.stoichiometry, result.molality = stoichiometry, molality
    result.species = [] if species is None else [_index(parameters.component_ids, species, "ion")]
    if reference is not None:
        for name in reference:
            _index(parameters.component_ids, name, "reference component")
        result.reference_composition = [reference.get(name, 0.0) for name in parameters.component_ids]
    return result


def mixture_osmotic_coefficient(parameters: Parameters, T: float, P: float, x: Mapping[str, float], *,
                                observed: float, scale: float) -> Observation:
    """Declare the practical osmotic coefficient of a full liquid composition.

    Parameters
    ----------
    parameters : epcsaft.parameters.Parameters
        Component order and solvent molar mass in kg/mol.
    T : float
        Temperature in K.
    P : float
        Pressure in Pa.
    x : mapping of str to float
        Full electroneutral dimensionless mole fractions with one solvent;
        omitted species are zero.
    observed : float
        Dimensionless practical osmotic coefficient.
    scale : float
        Positive dimensionless difference-residual scale.

    Returns
    -------
    observation : Observation
        Liquid observation with residual ``(Phi-observed)/scale``.

    Raises
    ------
    ParameterError
        Unknown species; native prediction checks the physical declaration.

    Notes
    -----
    Native evaluation uses ``Phi = -ln(a_w)/(M_w*sum_j m_j)`` with natural
    logarithm, solvent molar mass M_w in kg/mol and individual solute-species
    molalities m_j in mol/kg solvent. This is a full-composition observation,
    not a single formula-ion dilution path (ALG-REFERENCE).
    """
    return _pressure_closed(parameters, ObservationKind.MixtureOsmoticCoefficient, T, P, x, "liquid", observed, "difference", scale, False)


def hydration_hold(reference: Observation, *, bracket: tuple[float, float], energy_tolerance: float,
                   diameter_tolerance: float, maximum_iterations: int) -> HydrationHold:
    """Declare an infinite-dilution energy constraint with dependent Born diameter.

    Parameters
    ----------
    reference : Observation
        Single-ion reference energy target in J/mol, at declared T (K), P (Pa)
        and solvent composition.
    bracket : tuple of float
        Lower and upper dependent Born diameters in native Å.
    energy_tolerance : float
        Positive energy-residual tolerance in J/mol.
    diameter_tolerance : float
        Positive diameter-bracket tolerance in Å.
    maximum_iterations : int
        Native constraint-solve iteration budget.

    Returns
    -------
    hold : HydrationHold
        Native declaration for ``H = RT*ln(phi_i at infinite dilution)-target``;
        H has units J/mol. No constraint solve occurs here.

    Notes
    -----
    Native coordinate evaluation owns the bracket solve and implicit exact
    chart derivative. Invalid or unresolved holds retain their causes (ALG-FIT).
    """
    result = HydrationHold()
    result.reference, result.lower, result.upper = reference, *bracket
    result.energy_tolerance, result.diameter_tolerance, result.maximum_iterations = energy_tolerance, diameter_tolerance, maximum_iterations
    return result


def evaluate_coordinates(parameters: Parameters, coordinates: Sequence[Coordinate], optimizer: Sequence[float],
                         holds: Sequence[HydrationHold], thermochemistry: ThermochemistryRecord | None = None,
                         solver_options: Mapping[str, float] | None = None) -> EvaluatedCoordinates:
    """Evaluate the physical parameter chart and its exact native derivative.

    Parameters
    ----------
    parameters : epcsaft.parameters.Parameters
        Starting native physical parameter set.
    coordinates : sequence of Coordinate
        Ordered affine coordinates and physical bounds.
    optimizer : sequence of float
        Dimensionless optimizer point in coordinate order.
    holds : sequence of HydrationHold
        Dependent Born-diameter constraints; empty means none.
    thermochemistry : ThermochemistryRecord, optional
        Fixed ideal-gas molar thermochemistry with intervals in K and standard Pa.
    solver_options : mapping of str to float, optional
        Native solver tolerances and work limits.

    Returns
    -------
    evaluation : EvaluatedCoordinates
        Completed physical chart values in native parameter units and exact D
        columns per dimensionless optimizer increment, including dependent holds.
        No unavailable chart or derivative record is returned.

    Raises
    ------
    ParameterError
        The chart, optimizer dimensions, bounds or hold declarations are invalid,
        or a hold reference/action is unavailable, its root is not bracketed,
        its derivative is singular, or its iteration budget is exhausted (ALG-FIT).
    """
    chart = CoordinateSet()
    chart.coordinates = list(coordinates)
    return _core.evaluate_coordinates(parameters._native, chart, list(optimizer), list(holds), thermochemistry, solver_options)


def evaluate_checks(parameters: Parameters, observations: Sequence[Observation], checks: Sequence[CheckGroup],
                    baseline: Parameters | None = None, solver_options: Mapping[str, float] | None = None) -> list[CheckResult]:
    """Evaluate unweighted final physical metrics on stored parameters.

    Parameters
    ----------
    parameters : epcsaft.parameters.Parameters
        Physical parameter values to check.
    observations : sequence of Observation
        Observation declarations used by each check group's row indices.
    checks : sequence of CheckGroup
        Native metric groups and admitted limits; units follow each named
        physical metric or its dimensionless relative/log definition.
    baseline : epcsaft.parameters.Parameters, optional
        Comparison physical parameter set for baseline checks.
    solver_options : mapping of str to float, optional
        Native solver tolerances and work budgets.

    Returns
    -------
    results : list of CheckResult
        Unweighted metrics and check diagnostics in check-group order, retaining
        unavailable observations and their causes.

    Raises
    ------
    ParameterError
        Native observations, groups or baseline requirements are inadmissible.

    Notes
    -----
    No hydration holds are re-solved. These metrics evaluate stored physical
    parameters, separately from weighted fit objectives (ALG-FIT).
    """
    return _core.evaluate_checks(parameters._native, list(observations), list(checks), None if baseline is None else baseline._native, solver_options)
