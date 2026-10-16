from __future__ import annotations

from dataclasses import dataclass
from collections.abc import Sequence

from . import _core

Activity = _core.ElectrolyteActivity
SolventActivity = _core.SolventActivity


@dataclass(frozen=True, slots=True)
class Reference:
    """Immutable snapshot of the physical inputs at which the native reference fugacities were solved."""

    _model: _core.HelmholtzModel
    T: float
    P: float
    solvent_fractions: tuple[float, ...]
    _ion_indices: tuple[int, int] | None
    _stoichiometry: tuple[int, int] | None
    _native: _core.InfiniteDilutionReference

    def _require_ions(self) -> None:
        if self._ion_indices is None:
            raise _core.ParameterError("IonicReferenceUnavailable: solvent-only reference")

    @property
    def ion_indices(self) -> tuple[int, int]:
        self._require_ions()
        assert self._ion_indices is not None
        return self._ion_indices

    @property
    def stoichiometry(self) -> tuple[int, int]:
        self._require_ions()
        assert self._stoichiometry is not None
        return self._stoichiometry

    @property
    def log_fugacity(self) -> tuple[float, float]:
        self._require_ions()
        return self._native.log_fugacity

    @property
    def molality(self) -> float:
        self._require_ions()
        return self._native.molality

    @property
    def convergence_error(self) -> float:
        self._require_ions()
        return self._native.convergence_error

    @property
    def pure_solvent_log_fugacity(self) -> float | None:
        return self._native.pure_solvent_log_fugacity

    def activity(self, molality: float) -> Activity:
        return _core.electrolyte_activity(
            self._model, self.T, self.P, self.solvent_fractions, self.ion_indices, self.stoichiometry, molality, self._native
        )

    def solvent_activity(self, *, x: Sequence[float]) -> SolventActivity:
        return _core.solvent_activity(self._model, self.T, self.P, x, self._native)
