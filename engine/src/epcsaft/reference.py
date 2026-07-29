from __future__ import annotations

from dataclasses import dataclass, replace
import math
from math import gcd
from typing import TYPE_CHECKING, Literal

import pint

from .eos import _eos_error, _state_magnitude
from .units import unit_registry

if TYPE_CHECKING:
    from .eos import Mixture, State


_MOLAR_GAS_CONSTANT = 8.31446261815324


# EqExclude: transport and buffer mechanics
def _log_fugacity(
    mixture: Mixture,
    *,
    temperature: pint.Quantity,
    pressure: pint.Quantity,
    phase: str,
    composition: tuple[float, ...],
    required_indices: tuple[int, ...],
) -> tuple[float, ...]:
    result = mixture.state(
        T=temperature,
        P=pressure,
        x=composition,
        phase=phase,
    )
    fugacity = result.fugacity
    if fugacity is None:
        raise _eos_error("reference-state-failed", reason="fugacity-unavailable")
    values = fugacity.ln_coefficient
    if any(not math.isfinite(values[index]) for index in required_indices):
        raise _eos_error("reference-state-failed", reason="nonfinite-fugacity")
    return values


@dataclass(frozen=True, slots=True)
class Activity:
    """Immutable electrolyte activity result in the complete model component order."""

    component_ids: tuple[str, ...]
    formula_unit_molality: pint.Quantity[float]
    mole_fractions: tuple[float, ...]
    ion_component_ids: tuple[str, str]
    single_ion: tuple[float, float]
    mean_mole_fraction: float
    mean_molality: float
    osmotic_coefficient: float | None


@dataclass(frozen=True, slots=True)
class TransferFreeEnergy:
    ion_component_ids: tuple[str, str]
    values: tuple[pint.Quantity[float], pint.Quantity[float]]


@dataclass(frozen=True, slots=True)
class ElectrolyteReference:
    component_ids: tuple[str, ...]
    ion_component_ids: tuple[str, str]
    formula_stoichiometry: tuple[int, int]
    solvent_mole_fractions: tuple[float, ...]
    temperature: pint.Quantity[float]
    pressure: pint.Quantity[float]
    phase: str
    infinite_dilution_log_fugacity_coefficient: tuple[float, float]
    reference_molality: pint.Quantity[float]
    convergence_error: float
    _mixture: Mixture
    _ion_indices: tuple[int, int]
    _molar_masses: tuple[float | None, ...]
    _pure_water_log_fugacity: float | None

    @property
    def solvation_free_energy(
        self,
    ) -> tuple[pint.Quantity[float], pint.Quantity[float]]:
        # EqID: infinite_dilution_solvation_energy
        scale = (
            _MOLAR_GAS_CONSTANT
            * self.temperature.to("kelvin").magnitude
            * unit_registry.joule
            / unit_registry.mole
        )
        return tuple(  # type: ignore[return-value]
            value * scale for value in self.infinite_dilution_log_fugacity_coefficient
        )

    def state(self, *, molality: object) -> State:
        """Evaluate an electrolyte state on this explicit formula-unit reference."""
        molality = _state_magnitude(
            molality,
            "mole / kilogram",
            "molality",
        )
        composition = tuple(
            self._mixture._native.reference_composition(
                self.solvent_mole_fractions,
                self._ion_indices,
                self.formula_stoichiometry,
                molality,
            )
        )
        water_index = (
            self.component_ids.index("water") if self._pure_water_log_fugacity is not None else None
        )
        finite_log_fugacity = _log_fugacity(
            self._mixture,
            temperature=self.temperature,
            pressure=self.pressure,
            phase=self.phase,
            composition=composition,
            required_indices=self._ion_indices
            + ((water_index,) if water_index is not None else ()),
        )
        log_activity = tuple(
            finite_log_fugacity[index] - reference
            for index, reference in zip(
                self._ion_indices,
                self.infinite_dilution_log_fugacity_coefficient,
                strict=True,
            )
        )
        # EqID: single_ion_activity_coefficient
        single_ion = tuple(math.exp(value) for value in log_activity)
        coefficient_sum = sum(self.formula_stoichiometry)
        # EqID: mean_ionic_activity_coefficient
        log_mean = (
            math.fsum(
                coefficient * value
                for coefficient, value in zip(
                    self.formula_stoichiometry,
                    log_activity,
                    strict=True,
                )
            )
            / coefficient_sum
        )
        mean_mole_fraction = math.exp(log_mean)
        solvent_molar_mass = math.fsum(
            fraction * (self._molar_masses[index] or 0.0)
            for index, fraction in enumerate(self.solvent_mole_fractions)
        )
        # EqID: molality_activity_conversion
        mean_molality = mean_mole_fraction / (1.0 + solvent_molar_mass * molality * coefficient_sum)
        osmotic_coefficient = None
        if self._pure_water_log_fugacity is not None and water_index is not None:
            log_solvent_activity = (
                math.log(composition[water_index])
                + finite_log_fugacity[water_index]
                - self._pure_water_log_fugacity
            )
            # EqID: practical_osmotic_coefficient
            osmotic_coefficient = -log_solvent_activity / (
                coefficient_sum * molality * solvent_molar_mass
            )
            if not math.isfinite(osmotic_coefficient):
                raise _eos_error("reference-state-failed", reason="nonfinite-osmotic-coefficient")
        activity = Activity(
            component_ids=self.component_ids,
            formula_unit_molality=molality * unit_registry.mole / unit_registry.kilogram,
            mole_fractions=composition,
            ion_component_ids=self.ion_component_ids,
            single_ion=single_ion,  # type: ignore[arg-type]
            mean_mole_fraction=mean_mole_fraction,
            mean_molality=mean_molality,
            osmotic_coefficient=osmotic_coefficient,
        )
        state = self._mixture.state(
            T=self.temperature,
            P=self.pressure,
            x=composition,
            phase=self.phase,  # type: ignore[arg-type]
        )
        return replace(state, _activity=activity)

    def transfer_to(self, target: ElectrolyteReference) -> TransferFreeEnergy:
        if not isinstance(target, ElectrolyteReference):
            raise _eos_error("invalid-reference-state", field="target")
        if (
            self.ion_component_ids != target.ion_component_ids
            or self.formula_stoichiometry != target.formula_stoichiometry
            or self.phase != target.phase
            or self.temperature.to("kelvin").magnitude != target.temperature.to("kelvin").magnitude
            or self.pressure.to("pascal").magnitude != target.pressure.to("pascal").magnitude
        ):
            raise _eos_error("invalid-reference-state", reason="incompatible-transfer")
        scale = (
            _MOLAR_GAS_CONSTANT
            * self.temperature.to("kelvin").magnitude
            * unit_registry.joule
            / unit_registry.mole
        )
        # EqID: transfer_free_energy
        values = tuple(
            (target_value - source_value) * scale
            for source_value, target_value in zip(
                self.infinite_dilution_log_fugacity_coefficient,
                target.infinite_dilution_log_fugacity_coefficient,
                strict=True,
            )
        )
        return TransferFreeEnergy(
            ion_component_ids=self.ion_component_ids,
            values=values,  # type: ignore[arg-type]
        )


def build_reference_state(
    eos: Mixture,
    *,
    temperature: object,
    pressure: object,
    solvent_mole_fractions: object,
    phase: Literal["liquid"],
) -> ElectrolyteReference:
    temperature_k = _state_magnitude(temperature, "kelvin", "temperature")
    pressure_pa = _state_magnitude(pressure, "pascal", "pressure")
    if phase != "liquid":
        raise _eos_error("invalid-reference-state", field="phase")
    try:
        solvent_fractions = tuple(float(value) for value in solvent_mole_fractions)  # type: ignore[union-attr]
    except (TypeError, ValueError) as exc:
        raise _eos_error("invalid-reference-state", field="solvent_mole_fractions") from exc

    model = eos._model
    charges = model.charge_number
    if model.topology not in {
        "ionic-ssm-ds",
        "ionic-fully-dissociated-debye-huckel",
    } or len(charges) != len(eos.component_ids):
        raise _eos_error("invalid-reference-state", reason="electrolyte-required")
    cations = tuple(index for index, charge in enumerate(charges) if charge > 0)
    anions = tuple(index for index, charge in enumerate(charges) if charge < 0)
    if len(cations) != 1 or len(anions) != 1:
        raise _eos_error("invalid-reference-state", reason="one-formula-unit-required")
    if (
        len(solvent_fractions) != len(eos.component_ids)
        or any(not math.isfinite(value) or value < 0.0 for value in solvent_fractions)
        or abs(math.fsum(solvent_fractions) - 1.0) > 1.0e-12
        or any(solvent_fractions[index] != 0.0 for index in (*cations, *anions))
    ):
        raise _eos_error("invalid-reference-state", field="solvent_mole_fractions")
    solvent_indices = tuple(index for index, charge in enumerate(charges) if charge == 0)
    if not solvent_indices or any(
        solvent_fractions[index] > 0.0 and (model.molar_mass_kg_per_mol[index] or 0.0) <= 0.0
        for index in solvent_indices
    ):
        raise _eos_error("invalid-reference-state", reason="solvent-molar-mass-required")

    cation, anion = cations[0], anions[0]
    divisor = gcd(abs(charges[cation]), abs(charges[anion]))
    # EqID: charge_balanced_formula_stoichiometry
    stoichiometry = (
        abs(charges[anion]) // divisor,
        abs(charges[cation]) // divisor,
    )
    ion_indices = (cation, anion)
    temperature_quantity = temperature_k * unit_registry.kelvin
    pressure_quantity = pressure_pa * unit_registry.pascal
    try:
        sequence = eos._native.electrolyte_reference_sequence(
            temperature_k,
            pressure_pa,
            solvent_fractions,
            ion_indices,
            stoichiometry,
            phase,
        )
    except (RuntimeError, ValueError) as exc:
        raise _eos_error("reference-state-failed", reason=str(exc)) from exc
    current = tuple(sequence.log_fugacity)
    convergence_error = sequence.convergence_error

    active_solvents = tuple(index for index in solvent_indices if solvent_fractions[index] > 0.0)
    pure_water_log_fugacity = None
    if len(active_solvents) == 1 and eos.component_ids[active_solvents[0]] == "water":
        pure_water_log_fugacity = sequence.pure_solvent_log_fugacity_coefficient
        if pure_water_log_fugacity is None or not math.isfinite(pure_water_log_fugacity):
            raise _eos_error("reference-state-failed", reason="pure-solvent-fugacity-unavailable")

    return ElectrolyteReference(
        component_ids=eos.component_ids,
        ion_component_ids=(eos.component_ids[cation], eos.component_ids[anion]),
        formula_stoichiometry=stoichiometry,
        solvent_mole_fractions=solvent_fractions,
        temperature=temperature_quantity,
        pressure=pressure_quantity,
        phase=phase,
        infinite_dilution_log_fugacity_coefficient=current,
        reference_molality=sequence.reference_molality
        * unit_registry.mole
        / unit_registry.kilogram,
        convergence_error=convergence_error,
        _mixture=eos,
        _ion_indices=ion_indices,
        _molar_masses=model.molar_mass_kg_per_mol,
        _pure_water_log_fugacity=pure_water_log_fugacity,
    )


__all__ = ("Activity", "ElectrolyteReference", "TransferFreeEnergy")
