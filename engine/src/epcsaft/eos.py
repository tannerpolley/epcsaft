from __future__ import annotations

from collections.abc import Sequence
from typing import TYPE_CHECKING

from . import _core
from .parameters import Parameters
from .reference import Reference

if TYPE_CHECKING:
    from .equilibrium import Problem

State = _core.State
PropertyObservable = _core.PropertyObservable
PropertyActionStatus = _core.PropertyActionStatus
DensityActionDirection = _core.DensityActionDirection
PressureActionDirection = _core.PressureActionDirection
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

    def __init__(self, parameters: Parameters) -> None:
        if len(parameters.component_ids) < 2:
            raise ParameterError("Mixture requires at least two components")
        self.parameters = parameters
        self.component_ids = parameters.component_ids
        self.site_ids = parameters.site_ids
        self._model = _core.HelmholtzModel(parameters._native)

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
        native = _core.infinite_dilution_reference(
            self._model, T, P, solvent_fractions, ion_indices, stoichiometry
        )
        return Reference(self._model, T, P, solvent_fractions, ion_indices, stoichiometry, native)

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

    def __init__(self, parameters: Parameters) -> None:
        if len(parameters.component_ids) != 1:
            raise ParameterError("Pure requires exactly one component")
        self.parameters = parameters
        self.component_ids = parameters.component_ids
        self.site_ids = parameters.site_ids
        self._model = _core.HelmholtzModel(parameters._native)

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
                item.composition or (0.0,),
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
                item.composition or (0.0,),
            )
            for item in directions
        )
        return self._model.pressure_actions(
            T, P, (1.0,), _branch(phase), observables, native, None,
        )

    def saturation(self, T: float) -> Problem:
        from .equilibrium import Amounts, Free, Phase, Pinned, Problem

        return Problem(
            phases=[Phase("L"), Phase("V", amount=Pinned(0.0))],
            T=T,
            P=Free(),
            feed=Amounts({self.component_ids[0]: 1.0}),
        )
