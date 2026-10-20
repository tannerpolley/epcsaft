"""Construct native Helmholtz fluid models and request state properties.

Temperatures use K, pressure Pa, molar density mol/m³ and mole fractions a
dimensionless molar basis. Total caloric state properties are molar and
require explicit ideal-gas thermochemistry; residual properties retain their
native reference. Exact physical directional actions follow ALG-ACTIONS.
Parameter names and array order come from ``Parameters``.
"""

from __future__ import annotations

import json
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import TYPE_CHECKING, Any, cast

from . import _core
from .equilibrium.problem import _vector
from .parameters import Parameters
from .reference import Reference

if TYPE_CHECKING:
    from .equilibrium import Problem

State = _core.State
IdealPolynomial = _core.IdealPolynomial
IdealShomate = _core.IdealShomate
IdealNASA7 = _core.IdealNASA7
IdealNASA9 = _core.IdealNASA9
IdealDIPPR107 = _core.IdealDIPPR107
IdealCorrelation = _core.IdealCorrelation
IdealInterval = _core.IdealInterval
ThermochemistryRecord = _core.ThermochemistryRecord
ThermochemistryStatus = _core.ThermochemistryStatus
DensityBranch = _core.DensityBranch
DensityDiagnostics = _core.DensityDiagnostics
ParameterFamily = _core.ParameterFamily
ActiveParameter = _core.ActiveParameter
GAS_CONSTANT_J_PER_MOL_K = _core.GAS_CONSTANT_J_PER_MOL_K
PropertyObservable = _core.PropertyObservable
PropertyActionStatus = _core.PropertyActionStatus
DensityActionDirection = _core.DensityActionDirection
PressureActionDirection = _core.PressureActionDirection
ParameterActionDirection = _core.ParameterActionDirection
PropertyActionDiagnostics = _core.PropertyActionDiagnostics
PropertyActionResult = _core.PropertyActionResult
PropertyActionBatch = _core.PropertyActionBatch
ParameterError = _core.ParameterError
EosError = _core.EosError
DensityStructureError = _core.DensityStructureError
property_asl_record = _core.property_asl_record


class ThermochemistryPacket:
    """Load and prepare source NASA9 ideal-gas thermochemistry once.

    Attributes
    ----------
    species_ids : tuple of str
        Source species identifiers in packet order; no physical unit.

    Notes
    -----
    Input temperatures are K and standard pressure is Pa. Native preparation
    rescales the nine NASA9 coefficients by the source gas constant divided by
    the Engine gas constant (both J/(mol K)). Intervals are lower-inclusive and
    upper-exclusive except for the last interval, whose upper end is included.
    ``select`` only copies and reorders the prepared records (ALG-CALORICS).
    """

    __slots__ = ("species_ids", "_record")

    def __init__(self, path: str | Path) -> None:
        """Read a NASA9 record file and prepare it in the Engine's gas-constant basis.

        Parameters
        ----------
        path : str or pathlib.Path
            JSON file with ``nasa-glenn-nasa9-records`` version 1 and coefficient
            order ``a1 a2 a3 a4 a5 a6 a7 b1 b2``. Interval bounds use K, standard
            pressure Pa, and the source gas constant J/(mol K).

        Raises
        ------
        ParameterError
            The form/version/order is unsupported, species identifiers repeat,
            or native thermochemistry admission refuses the records.
        OSError
            The file cannot be read.
        json.JSONDecodeError
            The file is not valid JSON. Missing required keys raise ``KeyError``.
        """
        document = json.loads(Path(path).read_text())
        form = (document.get("schema"), document.get("schema_version"), document.get("coefficient_order"))
        if form != ("nasa-glenn-nasa9-records", 1, "a1 a2 a3 a4 a5 a6 a7 b1 b2".split()):
            raise ParameterError("unsupported thermochemistry packet form or version")
        p0, species = document["standard_pressure_Pa"], document["species"]
        self.species_ids: tuple[str, ...] = tuple(item["species_id"] for item in species)
        if len(set(self.species_ids)) != len(self.species_ids):
            raise ParameterError("duplicate thermochemistry species")

        def interval(i: Any, last: bool) -> IdealInterval:
            return IdealInterval(i["lower_K"], i["upper_K"], True, last, IdealCorrelation(IdealNASA9(i["coefficients"]), p0))

        ordered = [sorted(item["intervals"], key=lambda record: record["lower_K"]) for item in species]
        record = ThermochemistryRecord(p0, [[interval(item, item is items[-1]) for item in items] for items in ordered])
        self._record = _core.prepare_thermochemistry(record, document["gas_constant_J_per_mol_K"])

    def select(self, species_ids: Sequence[str]) -> ThermochemistryRecord:
        """Copy prepared thermochemistry into the requested species order.

        Parameters
        ----------
        species_ids : sequence of str
            Identifiers to select, in the consuming model's component order.

        Returns
        -------
        record : ThermochemistryRecord
            Prepared NASA9 intervals, retaining K bounds and standard pressure in Pa.

        Raises
        ------
        ParameterError
            Any identifier is absent from the packet.
        """
        if missing := [item for item in species_ids if item not in self.species_ids]:
            raise ParameterError(f"thermochemistry packet has no species {missing[0]}")
        prepared = self._record.component_intervals
        return ThermochemistryRecord(self._record.standard_pressure, [prepared[self.species_ids.index(i)] for i in species_ids])


def _branch(phase: str | None) -> _core.DensityBranch:
    if phase is None:
        return _core.DensityBranch.Unspecified
    if phase == "vapor":
        return _core.DensityBranch.Vapor
    if phase == "liquid":
        return _core.DensityBranch.Liquid
    raise ParameterError(f"unknown phase {phase}; expected liquid, vapor or None")


class Mixture:
    """Evaluate a fluid mixture with an explicit mole-fraction composition.

    Attributes
    ----------
    parameters : epcsaft.parameters.Parameters
        Parameter declaration used to construct the model.
    component_ids : tuple of str
        Component order for all compositions and component result arrays.
    site_ids : tuple of tuple of str
        Association-site order as (component identifier, site identifier).

    Notes
    -----
    At least two components are required for ``Mixture`` itself. ``Pure``
    uses the same calculations at dimensionless mole fraction ``x = (1,)``.
    Native construction copies parameter, thermochemistry and solver inputs;
    subsequent caller edits do not redefine the constructed native model.
    An explicit x sequence needs one finite nonnegative entry per component and
    a sum of one within native tolerance; an empty sequence is refused. State
    and property-action requests also require electroneutrality.
    """

    __slots__ = ("_model", "component_ids", "parameters", "site_ids")

    parameters: Parameters
    component_ids: tuple[str, ...]
    site_ids: tuple[tuple[str, str], ...]
    _model: _core.HelmholtzModel

    def __init__(
        self,
        parameters: Parameters,
        *,
        thermochemistry: ThermochemistryRecord | None = None,
        active_parameters: Sequence[ActiveParameter] | None = None,
        solver_options: Mapping[str, float] | None = None,
    ) -> None:
        """Construct the native Helmholtz model for the declared components.

        Parameters
        ----------
        parameters : epcsaft.parameters.Parameters
            Ordered physical parameter declaration; at least two components.
        thermochemistry : ThermochemistryRecord, optional
            Ideal-gas molar thermochemistry in component order, with temperature
            intervals in K and reference pressure in Pa. Required for total calorics.
        active_parameters : sequence of ActiveParameter, optional
            Native parameter coordinates whose exact derivatives are requested.
            Each coordinate retains the physical unit of its parameter family.
        solver_options : mapping of str to float, optional
            Native solver tolerances and work limits; each named option retains
            its native unit (ALG-TERMINATION).

        Raises
        ------
        ParameterError
            Fewer than two components, or native input admission fails.
        """
        if len(parameters.component_ids) < 2 and type(self) is Mixture:
            raise ParameterError("Mixture requires at least two components")
        self.parameters = parameters
        self.component_ids = parameters.component_ids
        self.site_ids = parameters.site_ids
        self._model = _core.HelmholtzModel(parameters._native, active_parameters or (), thermochemistry, solver_options)

    def state(
        self,
        T: float,
        *,
        rho: float | None = None,
        P: float | None = None,
        x: Sequence[float] | Mapping[str, float] | None = None,
        phase: str | None = None,
        anchor: State | None = None,
        derived_properties: bool | None = None,
    ) -> State:
        """Evaluate a state at fixed molar density or solve density at fixed pressure.

        Parameters
        ----------
        T : float
            Absolute temperature in K.
        rho : float, optional
            Molar density in mol/m³; specify exactly one of ``rho`` and ``P``.
        P : float, optional
            Pressure in Pa for a native density solve.
        x : sequence of float or mapping of str to float, optional
            Dimensionless mole fractions in component order. Unnamed mapping
            components are zero; a mixture requires an explicit composition.
        phase : {'liquid', 'vapor', None}, optional
            Density branch for a pressure solve; ``None`` leaves it unspecified.
        anchor : State, optional
            Previous density (mol/m³) used as a pressure-solve anchor.
        derived_properties : bool, optional
            Request derived properties; ``None`` is passed as false.

        Returns
        -------
        state : State
            Native state: pressure in Pa, density in mol/m³, dimensionless residual
            Helmholtz energy divided by RT, and available molar calorics in J/mol
            or J/(mol K). Optional properties retain their availability diagnostics.

        Raises
        ------
        ParameterError
            Both/neither density and pressure are supplied, phase is unknown,
            composition is missing or names/physical inputs are inadmissible.
        DensityStructureError
            The requested pressure-density branch cannot be established.
        EosError
            Native thermodynamic evaluation fails.

        Notes
        -----
        For a density input, branch and anchor do not select a root. Density solves
        and their finite-search limits follow ALG-DENSITY. State availability does
        not establish global phase stability.
        """
        if (rho is None) == (P is None):
            raise ParameterError("exactly one of rho or P is required")
        branch = _branch(phase)
        if rho is not None:
            return self._model.state(T, rho, self._composition(x), bool(derived_properties))
        assert P is not None
        return self._model.state_at_pressure(
            T, P, self._composition(x), branch, anchor.molar_density if anchor else None, bool(derived_properties)
        )

    def reference(
        self,
        T: float,
        P: float,
        *,
        solvent_fractions: Sequence[float],
        ions: tuple[str, str] | None = None,
        stoichiometry: tuple[int, int] | None = None,
    ) -> Reference:
        """Solve a solvent or electrolyte infinite-dilution reference at fixed T and P.

        Parameters
        ----------
        T : float
            Reference temperature in K.
        P : float
            Reference pressure in Pa.
        solvent_fractions : sequence of float
            Dimensionless solvent-only mole fractions in component order.
        ions : tuple of str, optional
            Two ion identifiers in the order used for returned ionic quantities.
        stoichiometry : tuple of int, optional
            Positive formula stoichiometric counts in ion order; required exactly
            when ``ions`` is supplied.

        Returns
        -------
        reference : Reference
            Immutable physical inputs and native reference fugacity coefficients.

        Raises
        ------
        ParameterError
            Only one of ions/stoichiometry is given, identifiers are unknown,
            ionic reference/composition inputs are invalid, or native model
            input admission, such as temperature or pressure, fails.
        EosError
            An ionic reference state or dilution limit cannot be evaluated.

        Notes
        -----
        The solvent-only path supports solvent activity and neutral Henry constants.
        An inadmissible solvent composition or numerical density failure returns
        a Reference with pure_solvent_log_fugacity set to None; its
        solvent_activity result retains the native status and cause.
        The ionic path follows the native charge-neutral dilution ladder and
        convergence checks (ALG-REFERENCE). Input sequences are copied.
        """
        if (ions is None) != (stoichiometry is None):
            raise ParameterError("reference requires both ions and stoichiometry or neither")
        if ions is None:
            fractions = tuple(solvent_fractions)
            native = _core.solvent_reference(self._model, T, P, fractions)
            return Reference(self._model, self.component_ids, T, P, fractions, None, None, native)
        assert stoichiometry is not None
        try:
            ion_indices = (self.component_ids.index(ions[0]), self.component_ids.index(ions[1]))
        except (ValueError, IndexError) as error:
            raise ParameterError(f"unknown ion component in {ions}") from error
        # The reference keeps its own copies: later caller edits cannot mix states with the native fugacities.
        fractions, formula = tuple(solvent_fractions), cast(tuple[int, int], tuple(stoichiometry))
        native = _core.infinite_dilution_reference(self._model, T, P, fractions, ion_indices, formula)
        return Reference(self._model, self.component_ids, T, P, fractions, ion_indices, formula, native)

    def density_actions(
        self,
        T: float,
        rho: float,
        *,
        x: Sequence[float] | None = None,
        observables: Sequence[PropertyObservable],
        directions: Sequence[DensityActionDirection],
    ) -> PropertyActionBatch:
        """Evaluate exact first through third directional property actions.

        Parameters
        ----------
        T : float
            Temperature in K.
        rho : float
            Molar density in mol/m³.
        x : sequence of float, optional
            Dimensionless mole fractions in component order; required for a mixture.
        observables : sequence of PropertyObservable
            Requested native properties in output order.
        directions : sequence of DensityActionDirection
            One to three physical directions: temperature in K, rho in mol/m³,
            composition dimensionless and parameters in native family units.
            Each composition increment must have one entry per component; an
            empty increment is refused. Supply an explicit zero vector when
            composition is fixed.

        Returns
        -------
        batch : PropertyActionBatch
            Values and directional contractions in each observable's unit:
            molar energy in J/mol and entropy/heat capacity in J/(mol K).
            LogFugacityCoefficient requires pressure_actions. Unavailable entries
            retain their native causes.

        Raises
        ------
        ParameterError
            Inputs, branch, composition or direction dimensions are inadmissible.
        EosError
            The native state evaluation fails.

        Notes
        -----
        Composition directions must preserve normalization and charge balance.
        Directions are physical increments, not dimensionless coordinate labels;
        contracted derivatives retain the observable unit. Higher mixed actions
        use the supplied directions together, rather than separate first actions
        (ALG-ACTIONS).
        """
        return self._model.density_actions(T, rho, self._composition(x), observables, directions)

    def pressure_actions(
        self,
        T: float,
        P: float,
        *,
        x: Sequence[float] | None = None,
        phase: str,
        observables: Sequence[PropertyObservable],
        directions: Sequence[PressureActionDirection],
        anchor: State | None = None,
    ) -> PropertyActionBatch:
        """Evaluate exact first through third directional property actions.

        Parameters
        ----------
        T : float
            Temperature in K.
        P : float
            Pressure in Pa.
        x : sequence of float, optional
            Dimensionless mole fractions in component order; required for a mixture.
        observables : sequence of PropertyObservable
            Requested native properties in output order.
        directions : sequence of PressureActionDirection
            One to three physical directions: temperature in K, P in Pa,
            composition dimensionless and parameters in native family units.
            Each composition increment must have one entry per component; an
            empty increment is refused. Supply an explicit zero vector when
            composition is fixed.
        phase : {'liquid', 'vapor'}
            Branch retained by the pressure-closed derivative calculation.
        anchor : State, optional
            Density anchor in mol/m³ for that branch.

        Returns
        -------
        batch : PropertyActionBatch
            Values and directional contractions in each observable's unit:
            molar energy in J/mol, entropy/heat capacity in J/(mol K), and
            dimensionless natural ln(phi), with one fugacity result per component
            in component order. Unavailable entries retain their native causes.

        Raises
        ------
        ParameterError
            Inputs, branch, composition or direction dimensions are inadmissible.
        EosError
            The native state evaluation fails.

        Notes
        -----
        Composition directions must preserve normalization and charge balance.
        Directions are physical increments, not dimensionless coordinate labels;
        contracted derivatives retain the observable unit. Higher mixed actions
        use the supplied directions together, rather than separate first actions
        (ALG-ACTIONS).
        """
        return self._model.pressure_actions(
            T, P, self._composition(x), _branch(phase), observables, directions,
            anchor.molar_density if anchor else None,
        )

    def density_diagnostics(
        self, T: float, P: float, *, x: Sequence[float] | None = None, nodes: Sequence[float],
        density_tolerance: float, log_tolerance: float,
        max_evaluations: int, max_seconds: float, adaptive: bool,
        finite_upper: float, evaluation_order: _core.DerivativeOrder,
        finite_anchor: float | None = None,
    ) -> DensityDiagnostics:
        """Inspect pressure roots and mechanical turning points on a declared grid.

        Parameters
        ----------
        T : float
            Temperature in K.
        P : float
            Target pressure in Pa.
        x : sequence of float, optional
            Dimensionless mole fractions in component order; required for a mixture.
        nodes : sequence of float
            At least two increasing positive densities in mol/m³, below the
            geometric density cap.
        density_tolerance : float
            Positive absolute density-bracket width in mol/m³.
        log_tolerance : float
            Positive dimensionless bracket width in natural log density ratio.
        max_evaluations : int
            Positive maximum number of native pressure evaluations.
        max_seconds : float
            Positive wall-time budget in s.
        adaptive : bool
            Enable native adaptive grid investigation.
        finite_upper : float
            Finite-domain upper density in mol/m³ when an anchor is supplied;
            otherwise unused.
        evaluation_order : DerivativeOrder
            Gradient, PressureSlope or Hessian for finite-domain evaluation.
        finite_anchor : float, optional
            Positive finite-domain anchor density in mol/m³, below the geometric cap
            and ``finite_upper``.

        Returns
        -------
        diagnostics : DensityDiagnostics
            Sampled pressure in Pa, densities and brackets in mol/m³, mechanical
            signs, work counts and native stop causes. Root residuals use Pa;
            spinodal residuals are pressure slopes divided by RT, dimensionless.

        Raises
        ------
        ParameterError
            The grid, budgets, composition, domain or derivative order is invalid.
        EosError
            Thermodynamic evaluation cannot be prepared.

        Notes
        -----
        Bracket refinement requires both absolute and natural-log width criteria.
        Partial results retain budget/domain causes; a finite scan is not a proof
        that no additional root exists (ALG-DENSITY).
        """
        return self._model.density_diagnostics(
            T, P, self._composition(x), nodes, density_tolerance, log_tolerance, max_evaluations,
            max_seconds, adaptive, finite_anchor, finite_upper, evaluation_order,
        )

    def molar_density_upper_bound(self, T: float, *, x: Sequence[float] | None = None) -> float:
        """Return the geometric packing limit for the declared composition.

        Parameters
        ----------
        T : float
            Temperature in K.
        x : sequence of float, optional
            Dimensionless mole fractions in component order; required for a mixture.

        Returns
        -------
        rho_cap : float
            Upper geometric density in mol/m³; not a stable-liquid density.

        Raises
        ------
        ParameterError
            Composition or temperature is inadmissible.
        """
        return self._model.molar_density_upper_bound(T, self._composition(x))

    def _composition(self, x: Sequence[float] | Mapping[str, float] | None) -> Sequence[float]:
        if x is not None:  # names absent from a mapping are zero, as in feeds
            return _vector(x, self.component_ids) if isinstance(x, Mapping) else x
        if len(self.component_ids) != 1:
            raise ParameterError("a mixture requires its composition x")
        return (1.0,)

    def counters(self) -> _core.SweepCounters:
        """Read native work counts accumulated by this model.

        Returns
        -------
        counters : SweepCounters
            Dimensionless counts of native evaluations and solver work since reset.
        """
        return self._model.counters()

    def reset_counters(self) -> None:
        """Reset this model's accumulated native work counts to zero.

        Returns
        -------
        None
            Subsequent ``counters`` reports work accumulated after this call.
        """
        self._model.reset_counters()


class Pure(Mixture):
    """Evaluate one component at mole fraction one, including saturation declarations.

    Notes
    -----
    Exactly one component is required. Component names, sites and parameters
    are the inherited ``Mixture`` attributes. State, reference and work-count
    methods follow ``Mixture``; no explicit composition is needed for a state.
    Saturation returns a declaration for the native equilibrium driver.
    """

    __slots__ = ()

    def __init__(
        self,
        parameters: Parameters,
        *,
        thermochemistry: ThermochemistryRecord | None = None,
        active_parameters: Sequence[ActiveParameter] | None = None,
        solver_options: Mapping[str, float] | None = None,
    ) -> None:
        """Construct the one-component native Helmholtz model.

        Parameters
        ----------
        parameters : epcsaft.parameters.Parameters
            Exactly one component, with source-declared physical parameter units.
        thermochemistry : ThermochemistryRecord, optional
            One component's ideal-gas molar thermochemistry; interval bounds in K
            and standard pressure in Pa. Required for total caloric properties.
        active_parameters : sequence of ActiveParameter, optional
            Exact derivative coordinates in native parameter-family units.
        solver_options : mapping of str to float, optional
            Native tolerances and work limits, as for ``Mixture``.

        Raises
        ------
        ParameterError
            Component count differs from one, or native model admission fails.
        """
        if len(parameters.component_ids) != 1:
            raise ParameterError("Pure requires exactly one component")
        super().__init__(parameters, thermochemistry=thermochemistry, active_parameters=active_parameters, solver_options=solver_options)

    def density_actions(
        self,
        T: float,
        rho: float,
        *,
        x: Sequence[float] | None = None,
        observables: Sequence[PropertyObservable],
        directions: Sequence[DensityActionDirection],
    ) -> PropertyActionBatch:
        """Evaluate exact first through third directional property actions.

        Parameters
        ----------
        T : float
            Temperature in K.
        rho : float
            Molar density in mol/m³.
        x : sequence of float, optional
            Dimensionless mole fractions in component order; defaults to (1,).
        observables : sequence of PropertyObservable
            Requested native properties in output order.
        directions : sequence of DensityActionDirection
            Physical increments: temperature in K, rho in mol/m³,
            composition dimensionless and parameters in their native family units.

        Returns
        -------
        batch : PropertyActionBatch
            Values and directional contractions in each observable's unit:
            molar energy in J/mol and entropy/heat capacity in J/(mol K).
            LogFugacityCoefficient requires pressure_actions. Unavailable entries
            retain their native causes.

        Raises
        ------
        ParameterError
            Inputs, branch, composition or direction dimensions are inadmissible.
        EosError
            The native state evaluation fails.

        Notes
        -----
        Composition directions must preserve normalization and charge balance.
        Directions are physical increments, not dimensionless coordinate labels;
        contracted derivatives retain the observable unit. Higher mixed actions
        use the supplied directions together, rather than separate first actions
        (ALG-ACTIONS). An omitted composition increment is replaced by zero.
        """
        native = tuple(DensityActionDirection(item.temperature, item.molar_density, item.composition or (0.0,), item.parameter)
                       for item in directions)
        return super().density_actions(T, rho, x=x, observables=observables, directions=native)

    def pressure_actions(
        self,
        T: float,
        P: float,
        *,
        x: Sequence[float] | None = None,
        phase: str,
        observables: Sequence[PropertyObservable],
        directions: Sequence[PressureActionDirection],
        anchor: State | None = None,
    ) -> PropertyActionBatch:
        """Evaluate exact first through third directional property actions.

        Parameters
        ----------
        T : float
            Temperature in K.
        P : float
            Pressure in Pa.
        x : sequence of float, optional
            Dimensionless mole fractions in component order; defaults to (1,).
        observables : sequence of PropertyObservable
            Requested native properties in output order.
        directions : sequence of PressureActionDirection
            Physical increments: temperature in K, P in Pa,
            composition dimensionless and parameters in their native family units.
        phase : {'liquid', 'vapor'}
            Branch retained by the pressure-closed derivative calculation.
        anchor : State, optional
            Density anchor in mol/m³ for that branch.

        Returns
        -------
        batch : PropertyActionBatch
            Values and directional contractions in each observable's unit:
            molar energy in J/mol, entropy/heat capacity in J/(mol K), and
            dimensionless natural ln(phi), with one fugacity result per component
            in component order. Unavailable entries retain their native causes.

        Raises
        ------
        ParameterError
            Inputs, branch, composition or direction dimensions are inadmissible.
        EosError
            The native state evaluation fails.

        Notes
        -----
        Composition directions must preserve normalization and charge balance.
        Directions are physical increments, not dimensionless coordinate labels;
        contracted derivatives retain the observable unit. Higher mixed actions
        use the supplied directions together, rather than separate first actions
        (ALG-ACTIONS). An omitted composition increment is replaced by zero.
        """
        native = tuple(PressureActionDirection(item.temperature, item.pressure, item.composition or (0.0,), item.parameter)
                       for item in directions)
        return super().pressure_actions(T, P, x=x, phase=phase, observables=observables, directions=native, anchor=anchor)

    def saturation(self, T: float) -> Problem:
        """Declare isothermal liquid-vapor saturation without solving it.

        Parameters
        ----------
        T : float
            Fixed temperature in K.

        Returns
        -------
        problem : Problem
            One mol of this component, liquid ``L`` and incipient vapor ``V``
            pinned at zero mol, with pressure free in Pa.

        Notes
        -----
        Call ``solve_equilibrium`` to solve the returned problem. This declaration
        alone establishes neither existence nor stability of a saturation state.
        """
        from .equilibrium import Amounts, Free, Phase, Pinned, Problem

        return Problem(
            phases=[Phase("L", kind="liquid"), Phase("V", kind="vapor", amount=Pinned(0.0))],
            T=T,
            P=Free(),
            feed=Amounts({self.component_ids[0]: 1.0}),
        )
