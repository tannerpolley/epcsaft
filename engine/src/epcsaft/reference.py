from __future__ import annotations

from collections.abc import Sequence

from . import _core

Activity = _core.ElectrolyteActivity


class Reference:
    __slots__ = ("_model", "T", "P", "solvent_fractions", "ion_indices", "stoichiometry", "_native", "log_fugacity", "molality", "convergence_error", "pure_solvent_log_fugacity")

    def __init__(
        self,
        model: _core.HelmholtzModel,
        T: float,
        P: float,
        solvent_fractions: Sequence[float],
        ion_indices: tuple[int, int],
        stoichiometry: tuple[int, int],
        native: _core.InfiniteDilutionReference,
    ) -> None:
        self._model = model
        self.T = T
        self.P = P
        self.solvent_fractions = solvent_fractions
        self.ion_indices = ion_indices
        self.stoichiometry = stoichiometry
        self._native = native
        self.log_fugacity = native.log_fugacity
        self.molality = native.molality
        self.convergence_error = native.convergence_error
        self.pure_solvent_log_fugacity = native.pure_solvent_log_fugacity

    def activity(self, molality: float) -> Activity:
        return _core.electrolyte_activity(
            self._model,
            self.T,
            self.P,
            self.solvent_fractions,
            self.ion_indices,
            self.stoichiometry,
            molality,
            self._native,
        )
