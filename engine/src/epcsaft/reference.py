"""Read solvent and infinite-dilution references solved by the native model.

``Mixture.reference`` constructs these immutable physical-input snapshots
at temperature in K and pressure in Pa. Log fugacity fields are natural
logarithms of dimensionless fugacity coefficients, not logarithms of Pa.
Activity coefficients use the declared molality or solvent reference.
"""

from __future__ import annotations

from dataclasses import dataclass
from collections.abc import Sequence

from . import _core

Activity = _core.ElectrolyteActivity
SolventActivity = _core.SolventActivity


@dataclass(frozen=True, slots=True)
class Reference:
    """Retain the physical inputs and native reference fugacity solution.

    Attributes
    ----------
    component_ids : tuple of str
        Component order for all compositions and returned component quantities.
    T : float
        Reference temperature in K.
    P : float
        Reference pressure in Pa.
    solvent_fractions : tuple of float
        Dimensionless solvent-only mole fractions in component order.
    ion_indices : tuple of int
        Zero-based component positions in declared ion order; solvent-only
        references refuse access with IonicReferenceUnavailable.
    stoichiometry : tuple of int
        Dimensionless positive formula counts in ion order; solvent-only
        references refuse access with IonicReferenceUnavailable.
    log_fugacity : tuple of float
        Dimensionless natural-log ionic fugacity-coefficient limits in ion
        order; solvent-only references refuse access.
    molality : float
        Zero formula molality in mol/kg solvent at the analytic limit;
        solvent-only references refuse access.
    convergence_error : float or None
        None: the analytic limit has no terminal sequence comparison.
        Solvent density and association diagnostics describe its inner calculations.
    pure_solvent_log_fugacity : float or None
        Dimensionless natural-log pure-solvent fugacity coefficient at the
        stored T and P. None retains native unavailability, not zero.

    Notes
    -----
    Created by ``Mixture.reference`` with copied physical inputs and native
    reference results. Ionic properties/methods refuse solvent-only references.
    The native model and reference records are internal constructor inputs;
    use the model's reference method to obtain a consistent snapshot.
    ALG-REFERENCE defines the dilution limit and its convergence diagnostic.
    """

    _model: _core.HelmholtzModel
    component_ids: tuple[str, ...]
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
        """Read ``ion_indices`` from the stored reference.

        Returns
        -------
        value : tuple of int
            Field documented with its unit, basis and availability in the
            ``Reference`` Attributes section.

        Raises
        ------
        ParameterError
            ``IonicReferenceUnavailable``: this reference is solvent-only.
        """
        self._require_ions()
        assert self._ion_indices is not None
        return self._ion_indices

    @property
    def stoichiometry(self) -> tuple[int, int]:
        """Read ``stoichiometry`` from the stored reference.

        Returns
        -------
        value : tuple of int
            Field documented with its unit, basis and availability in the
            ``Reference`` Attributes section.

        Raises
        ------
        ParameterError
            ``IonicReferenceUnavailable``: this reference is solvent-only.
        """
        self._require_ions()
        assert self._stoichiometry is not None
        return self._stoichiometry

    @property
    def log_fugacity(self) -> tuple[float, float]:
        """Read ``log_fugacity`` from the stored reference.

        Returns
        -------
        value : tuple of float
            Field documented with its unit, basis and availability in the
            ``Reference`` Attributes section.

        Raises
        ------
        ParameterError
            ``IonicReferenceUnavailable``: this reference is solvent-only.
        """
        self._require_ions()
        return self._native.log_fugacity

    @property
    def molality(self) -> float:
        """Read ``molality`` from the stored reference.

        Returns
        -------
        value : float
            Field documented with its unit, basis and availability in the
            ``Reference`` Attributes section.

        Raises
        ------
        ParameterError
            ``IonicReferenceUnavailable``: this reference is solvent-only.
        """
        self._require_ions()
        return self._native.molality

    @property
    def convergence_error(self) -> float | None:
        """Read ``convergence_error`` from the stored reference.

        Returns
        -------
        value : float
            Field documented with its unit, basis and availability in the
            ``Reference`` Attributes section.

        Raises
        ------
        ParameterError
            ``IonicReferenceUnavailable``: this reference is solvent-only.
        """
        self._require_ions()
        return self._native.convergence_error

    @property
    def pure_solvent_log_fugacity(self) -> float | None:
        """Read ``pure_solvent_log_fugacity`` from the stored reference.

        Returns
        -------
        value : float or None
            Field documented with its unit, basis and availability in the
            ``Reference`` Attributes section.
        """
        return self._native.pure_solvent_log_fugacity

    def activity(self, molality: float) -> Activity:
        """Evaluate electrolyte activities along this reference's formula dilution path.

        Parameters
        ----------
        molality : float
            Finite, strictly positive formula molality in mol/kg solvent at the
            stored T (K) and P (Pa). Zero molality raises ParameterError.

        Returns
        -------
        activity : Activity
            Dimensionless molality-scale mean ionic coefficient and osmotic
            coefficient, single-ion quantities and native availability diagnostics.

        Raises
        ------
        ParameterError
            The reference is solvent-only, or formula/composition inputs are invalid.
        EosError
            A liquid state or reference evaluation fails.

        Notes
        -----
        Stoichiometric counts convert formula molality to individual-ion molalities.
        Single-ion quantities retain the declared native reference convention;
        only neutral combinations are convention-invariant (ALG-REFERENCE).
        """
        return _core.electrolyte_activity(
            self._model, self.T, self.P, self.solvent_fractions, self.ion_indices, self.stoichiometry, molality, self._native
        )

    def solvent_activity(self, *, x: Sequence[float]) -> SolventActivity:
        """Evaluate solvent activity for an explicit liquid composition.

        Parameters
        ----------
        x : sequence of float
            Dimensionless mole fractions in component order, at stored T and P.

        Returns
        -------
        activity : SolventActivity
            Dimensionless solvent fugacity relative to the pure-liquid reference
            and the practical osmotic coefficient when available. An unavailable
            pure reference, zero solvent fraction or failed liquid-density solve
            returns status/cause with unavailable values. Osmotic unavailability
            has its own osmotic_cause and can coexist with available activity.

        Raises
        ------
        ParameterError
            With an available pure reference, the composition has the wrong
            length, nonfinite or negative entries, or violates normalization
            or charge balance.

        Notes
        -----
        The native definition is ``a_w = f_w / f_w_pure`` at the same K and Pa;
        its logarithm is natural. It does not replace ionic mean activities.
        """
        return _core.solvent_activity(self._model, self.T, self.P, x, self._native)

    def henry_constants(self, *, solutes: Sequence[str]) -> list[_core.HenryConstant]:
        """Evaluate neutral-solute fugacity-based Henry constants in a pure solvent.

        Parameters
        ----------
        solutes : sequence of str
            Nonempty, distinct neutral component identifiers excluding the solvent.
            Results are sorted by component index, regardless of this input order.

        Returns
        -------
        constants : list of HenryConstant
            H_i^px = lim(f_i/x_i) in Pa, sorted by component_index, with
            dimensionless x_i and fugacity f_i in Pa. An ion-containing model or
            failed reference evaluation returns unavailable entries with status
            and message; inspect each henry_constant_px before use.

        Raises
        ------
        ParameterError
            This is an ionic Reference, a solute name is unknown, the solvent
            composition is not pure, or the solute set is empty, repeated or
            contains the solvent.
        """
        if self._ion_indices is not None:
            raise _core.ParameterError("NonNeutralDeclaration: Henry constants need a solvent-only reference")
        try:
            indices = [self.component_ids.index(solute) for solute in solutes]
        except ValueError as error:
            raise _core.ParameterError(f"unknown solute component in {list(solutes)}") from error
        return _core.henry_constants(self._model, self.T, self.P, self.solvent_fractions, indices)
