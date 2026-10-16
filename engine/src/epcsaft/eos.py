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
    """NASA9 species records loaded once: native preparation rescales the nine coefficients by
    R_source/R_Engine; select only copies and reorders. Intervals are [lower, upper), the last closed."""

    __slots__ = ("species_ids", "_record")

    def __init__(self, path: str | Path) -> None:
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
    """A mixture states its composition x; Pure, with one component, has x = (1,)."""

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
        """Thermal intervals follow component_ids order; the native model copies every input."""
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
    ) -> State:
        if (rho is None) == (P is None):
            raise ParameterError("exactly one of rho or P is required")
        branch = _branch(phase)
        if rho is not None:
            return self._model.state(T, rho, self._composition(x))
        assert P is not None
        return self._model.state_at_pressure(T, P, self._composition(x), branch, anchor.molar_density if anchor else None)

    def reference(
        self,
        T: float,
        P: float,
        *,
        solvent_fractions: Sequence[float],
        ions: tuple[str, str] | None = None,
        stoichiometry: tuple[int, int] | None = None,
    ) -> Reference:
        if (ions is None) != (stoichiometry is None):
            raise ParameterError("reference requires both ions and stoichiometry or neither")
        if ions is None:
            fractions = tuple(solvent_fractions)
            native = _core.solvent_reference(self._model, T, P, fractions)
            return Reference(self._model, T, P, fractions, None, None, native)
        assert stoichiometry is not None
        try:
            ion_indices = (self.component_ids.index(ions[0]), self.component_ids.index(ions[1]))
        except (ValueError, IndexError) as error:
            raise ParameterError(f"unknown ion component in {ions}") from error
        # The reference keeps its own copies: later caller edits cannot mix states with the native fugacities.
        fractions, formula = tuple(solvent_fractions), cast(tuple[int, int], tuple(stoichiometry))
        native = _core.infinite_dilution_reference(self._model, T, P, fractions, ion_indices, formula)
        return Reference(self._model, T, P, fractions, ion_indices, formula, native)

    def density_actions(
        self,
        T: float,
        rho: float,
        *,
        x: Sequence[float] | None = None,
        observables: Sequence[PropertyObservable],
        directions: Sequence[DensityActionDirection],
    ) -> PropertyActionBatch:
        """Evaluate order 1-3 SI-coordinate actions in the observable's unit."""
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
        """Evaluate order 1-3 SI-coordinate actions in the observable's unit."""
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
        return self._model.density_diagnostics(
            T, P, self._composition(x), nodes, density_tolerance, log_tolerance, max_evaluations,
            max_seconds, adaptive, finite_anchor, finite_upper, evaluation_order,
        )

    def molar_density_upper_bound(self, T: float, *, x: Sequence[float] | None = None) -> float:
        return self._model.molar_density_upper_bound(T, self._composition(x))

    def _composition(self, x: Sequence[float] | Mapping[str, float] | None) -> Sequence[float]:
        if x is not None:  # names absent from a mapping are zero, as in feeds
            return _vector(x, self.component_ids) if isinstance(x, Mapping) else x
        if len(self.component_ids) != 1:
            raise ParameterError("a mixture requires its composition x")
        return (1.0,)

    def counters(self) -> _core.SweepCounters:
        return self._model.counters()

    def reset_counters(self) -> None:
        self._model.reset_counters()


class Pure(Mixture):
    """One component: the Mixture methods at x = (1,), plus saturation."""

    __slots__ = ()

    def __init__(
        self,
        parameters: Parameters,
        *,
        thermochemistry: ThermochemistryRecord | None = None,
        active_parameters: Sequence[ActiveParameter] | None = None,
        solver_options: Mapping[str, float] | None = None,
    ) -> None:
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
        """A direction may omit its composition change, which is zero for one component."""
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
        """A direction may omit its composition change, which is zero for one component."""
        native = tuple(PressureActionDirection(item.temperature, item.pressure, item.composition or (0.0,), item.parameter)
                       for item in directions)
        return super().pressure_actions(T, P, x=x, phase=phase, observables=observables, directions=native, anchor=anchor)

    def saturation(self, T: float) -> Problem:
        from .equilibrium import Amounts, Free, Phase, Pinned, Problem

        return Problem(
            phases=[Phase("L", kind="liquid"), Phase("V", kind="vapor", amount=Pinned(0.0))],
            T=T,
            P=Free(),
            feed=Amounts({self.component_ids[0]: 1.0}),
        )
