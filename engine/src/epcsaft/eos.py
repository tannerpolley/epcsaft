from __future__ import annotations

from collections.abc import Sequence
from typing import TYPE_CHECKING, cast

from . import _core
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


def _branch(phase: str | None) -> _core.DensityBranch:
    if phase is None:
        return _core.DensityBranch.Unspecified
    if phase == "vapor":
        return _core.DensityBranch.Vapor
    if phase == "liquid":
        return _core.DensityBranch.Liquid
    raise ParameterError(f"unknown phase {phase}")


class Mixture:
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
    ) -> None:
        """Thermal intervals follow component_ids order; the native model copies every input."""
        if len(parameters.component_ids) < 2:
            raise ParameterError("Mixture requires at least two components")
        self.parameters = parameters
        self.component_ids = parameters.component_ids
        self.site_ids = parameters.site_ids
        self._model = _core.HelmholtzModel(parameters._native, active_parameters or (), thermochemistry)

    def state(
        self,
        T: float,
        *,
        rho: float | None = None,
        P: float | None = None,
        x: Sequence[float],
        phase: str | None = None,
        anchor: State | None = None,
    ) -> State:
        if (rho is None) == (P is None):
            raise ParameterError("exactly one of rho or P is required")
        branch = _branch(phase)
        if rho is not None:
            return self._model.state(T, rho, x)
        assert P is not None
        return self._model.state_at_pressure(T, P, x, branch, anchor.molar_density if anchor else None)

    def reference(
        self,
        T: float,
        P: float,
        *,
        solvent_fractions: Sequence[float],
        ions: tuple[str, str],
        stoichiometry: tuple[int, int],
    ) -> Reference:
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
        x: Sequence[float],
        observables: Sequence[PropertyObservable],
        directions: Sequence[DensityActionDirection],
    ) -> PropertyActionBatch:
        """Evaluate order 1-3 SI-coordinate actions in the observable's unit."""
        return self._model.density_actions(T, rho, x, observables, directions)

    def pressure_actions(
        self,
        T: float,
        P: float,
        *,
        x: Sequence[float],
        phase: str,
        observables: Sequence[PropertyObservable],
        directions: Sequence[PressureActionDirection],
        anchor: State | None = None,
    ) -> PropertyActionBatch:
        """Evaluate order 1-3 SI-coordinate actions in the observable's unit."""
        return self._model.pressure_actions(
            T, P, x, _branch(phase), observables, directions,
            anchor.molar_density if anchor else None,
        )

    def density_diagnostics(
        self, T: float, P: float, *, x: Sequence[float], nodes: Sequence[float],
        density_tolerance: float = 1e-9, log_tolerance: float = 1e-14,
        max_evaluations: int = 100000, max_seconds: float = 60.0, adaptive: bool = True,
        finite_anchor: float | None = None, finite_upper: float = 0.0,
        evaluation_order: _core.DerivativeOrder = _core.DerivativeOrder.PressureSlope,
    ) -> DensityDiagnostics:
        return self._model.density_diagnostics(
            T, P, x, nodes, density_tolerance, log_tolerance, max_evaluations,
            max_seconds, adaptive, finite_anchor, finite_upper, evaluation_order,
        )

    def molar_density_upper_bound(self, T: float, *, x: Sequence[float]) -> float:
        return self._model.molar_density_upper_bound(T, x)

    def counters(self) -> _core.SweepCounters:
        return self._model.counters()

    def reset_counters(self) -> None:
        self._model.reset_counters()


class Pure:
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
    ) -> None:
        if len(parameters.component_ids) != 1:
            raise ParameterError("Pure requires exactly one component")
        self.parameters = parameters
        self.component_ids = parameters.component_ids
        self.site_ids = parameters.site_ids
        self._model = _core.HelmholtzModel(parameters._native, active_parameters or (), thermochemistry)

    def state(
        self,
        T: float,
        *,
        rho: float | None = None,
        P: float | None = None,
        phase: str | None = None,
    ) -> State:
        if (rho is None) == (P is None):
            raise ParameterError("exactly one of rho or P is required")
        branch = _branch(phase)
        if rho is not None:
            return self._model.state(T, rho, (1.0,))
        assert P is not None
        return self._model.state_at_pressure(T, P, (1.0,), branch, None)

    def density_actions(
        self,
        T: float,
        rho: float,
        *,
        observables: Sequence[PropertyObservable],
        directions: Sequence[DensityActionDirection],
    ) -> PropertyActionBatch:
        """Evaluate order 1-3 SI-coordinate actions in the observable's unit."""
        native = tuple(
            DensityActionDirection(
                item.temperature, item.molar_density,
                item.composition or (0.0,), item.parameter,
            )
            for item in directions
        )
        return self._model.density_actions(T, rho, (1.0,), observables, native)

    def pressure_actions(
        self,
        T: float,
        P: float,
        *,
        phase: str,
        observables: Sequence[PropertyObservable],
        directions: Sequence[PressureActionDirection],
    ) -> PropertyActionBatch:
        """Evaluate order 1-3 SI-coordinate actions in the observable's unit."""
        native = tuple(
            PressureActionDirection(
                item.temperature, item.pressure,
                item.composition or (0.0,), item.parameter,
            )
            for item in directions
        )
        return self._model.pressure_actions(
            T, P, (1.0,), _branch(phase), observables, native, None,
        )

    def density_diagnostics(
        self, T: float, P: float, *, nodes: Sequence[float],
        density_tolerance: float = 1e-9, log_tolerance: float = 1e-14,
        max_evaluations: int = 100000, max_seconds: float = 60.0, adaptive: bool = True,
        finite_anchor: float | None = None, finite_upper: float = 0.0,
        evaluation_order: _core.DerivativeOrder = _core.DerivativeOrder.PressureSlope,
    ) -> DensityDiagnostics:
        return self._model.density_diagnostics(
            T, P, (1.0,), nodes, density_tolerance, log_tolerance, max_evaluations,
            max_seconds, adaptive, finite_anchor, finite_upper, evaluation_order,
        )

    def molar_density_upper_bound(self, T: float) -> float:
        return self._model.molar_density_upper_bound(T, (1.0,))

    def saturation(self, T: float) -> Problem:
        from .equilibrium import Amounts, Free, Phase, Pinned, Problem

        return Problem(
            phases=[Phase("L", kind="liquid"), Phase("V", kind="vapor", amount=Pinned(0.0))],
            T=T,
            P=Free(),
            feed=Amounts({self.component_ids[0]: 1.0}),
        )
