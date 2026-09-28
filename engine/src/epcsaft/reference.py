from __future__ import annotations

from dataclasses import dataclass

from . import _core

Activity = _core.ElectrolyteActivity


@dataclass(frozen=True, slots=True)
class Reference:
    """Immutable snapshot of the physical inputs at which the native reference fugacities were solved."""

    _model: _core.HelmholtzModel
    T: float
    P: float
    solvent_fractions: tuple[float, ...]
    ion_indices: tuple[int, int]
    stoichiometry: tuple[int, int]
    _native: _core.InfiniteDilutionReference

    @property
    def log_fugacity(self) -> tuple[float, float]:
        return self._native.log_fugacity

    @property
    def molality(self) -> float:
        return self._native.molality

    @property
    def convergence_error(self) -> float:
        return self._native.convergence_error

    @property
    def pure_solvent_log_fugacity(self) -> float | None:
        return self._native.pure_solvent_log_fugacity

    def activity(self, molality: float) -> Activity:
        return _core.electrolyte_activity(
            self._model, self.T, self.P, self.solvent_fractions, self.ion_indices, self.stoichiometry, molality, self._native
        )
