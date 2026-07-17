from __future__ import annotations

from dataclasses import dataclass
import math
from math import gcd
from typing import TYPE_CHECKING

import pint

from .eos import _eos_error, _state_magnitude
from .units import unit_registry

if TYPE_CHECKING:
    from .eos import EPCSAFT


_REFERENCE_MOLALITIES = (1.0e-6, 1.0e-8, 1.0e-10, 1.0e-12)
_REFERENCE_CONVERGENCE_TOLERANCE = 5.0e-5
_MOLAR_GAS_CONSTANT = 8.31446261815324


def _formula_composition(
    solvent_fractions: tuple[float, ...],
    molar_masses: tuple[float | None, ...],
    ion_indices: tuple[int, int],
    stoichiometry: tuple[int, int],
    molality: float,
) -> tuple[float, ...]:
    solvent_molar_mass = math.fsum(
        fraction * (molar_masses[index] or 0.0)
        for index, fraction in enumerate(solvent_fractions)
    )
    solvent_moles = 1.0 / solvent_molar_mass
    total_moles = solvent_moles + molality * sum(stoichiometry)
    composition = [fraction * solvent_moles / total_moles for fraction in solvent_fractions]
    for index, coefficient in zip(ion_indices, stoichiometry, strict=True):
        composition[index] = coefficient * molality / total_moles
    return tuple(composition)


def _ion_log_fugacity(
    eos: EPCSAFT,
    *,
    temperature: pint.Quantity,
    pressure: pint.Quantity,
    phase: str,
    composition: tuple[float, ...],
    ion_indices: tuple[int, int],
) -> tuple[float, float]:
    result = eos.evaluate(
        temperature=temperature,
        pressure=pressure,
        mole_fractions=composition,
        phase=phase,
    )
    values = result.log_fugacity_coefficient
    if values is None:
        raise _eos_error("reference-state-failed", reason="fugacity-unavailable")
    selected = tuple(values[index] for index in ion_indices)
    if any(not math.isfinite(value) for value in selected):
        raise _eos_error("reference-state-failed", reason="nonfinite-fugacity")
    return selected  # type: ignore[return-value]


@dataclass(frozen=True, slots=True)
class IonicActivity:
    formula_unit_molality: pint.Quantity
    mole_fractions: tuple[float, ...]
    ion_component_ids: tuple[str, str]
    single_ion_activity_coefficients: tuple[float, float]
    mean_ionic_activity_coefficient_mole_fraction: float
    mean_ionic_activity_coefficient_molality: float


@dataclass(frozen=True, slots=True)
class TransferFreeEnergy:
    ion_component_ids: tuple[str, str]
    values: tuple[pint.Quantity, pint.Quantity]


@dataclass(frozen=True, slots=True)
class ElectrolyteReference:
    component_ids: tuple[str, ...]
    ion_component_ids: tuple[str, str]
    formula_stoichiometry: tuple[int, int]
    solvent_mole_fractions: tuple[float, ...]
    temperature: pint.Quantity
    pressure: pint.Quantity
    phase: str
    infinite_dilution_log_fugacity: tuple[float, float]
    reference_molality: pint.Quantity
    convergence_error: float
    _eos: EPCSAFT
    _ion_indices: tuple[int, int]
    _molar_masses: tuple[float | None, ...]

    @property
    def solvation_free_energy(self) -> tuple[pint.Quantity, pint.Quantity]:
        scale = (
            _MOLAR_GAS_CONSTANT
            * self.temperature.to("kelvin").magnitude
            * unit_registry.joule
            / unit_registry.mole
        )
        return tuple(  # type: ignore[return-value]
            value * scale for value in self.infinite_dilution_log_fugacity
        )

    def at_molality(self, formula_unit_molality: object) -> IonicActivity:
        molality = _state_magnitude(
            formula_unit_molality,
            "mole / kilogram",
            "formula_unit_molality",
        )
        composition = _formula_composition(
            self.solvent_mole_fractions,
            self._molar_masses,
            self._ion_indices,
            self.formula_stoichiometry,
            molality,
        )
        finite_log_fugacity = _ion_log_fugacity(
            self._eos,
            temperature=self.temperature,
            pressure=self.pressure,
            phase=self.phase,
            composition=composition,
            ion_indices=self._ion_indices,
        )
        log_activity = tuple(
            value - reference
            for value, reference in zip(
                finite_log_fugacity,
                self.infinite_dilution_log_fugacity,
                strict=True,
            )
        )
        single_ion = tuple(math.exp(value) for value in log_activity)
        coefficient_sum = sum(self.formula_stoichiometry)
        log_mean = math.fsum(
            coefficient * value
            for coefficient, value in zip(
                self.formula_stoichiometry,
                log_activity,
                strict=True,
            )
        ) / coefficient_sum
        mean_mole_fraction = math.exp(log_mean)
        solvent_molar_mass = math.fsum(
            fraction * (self._molar_masses[index] or 0.0)
            for index, fraction in enumerate(self.solvent_mole_fractions)
        )
        mean_molality = mean_mole_fraction / (
            1.0 + solvent_molar_mass * molality * coefficient_sum
        )
        return IonicActivity(
            formula_unit_molality=molality * unit_registry.mole / unit_registry.kilogram,
            mole_fractions=composition,
            ion_component_ids=self.ion_component_ids,
            single_ion_activity_coefficients=single_ion,  # type: ignore[arg-type]
            mean_ionic_activity_coefficient_mole_fraction=mean_mole_fraction,
            mean_ionic_activity_coefficient_molality=mean_molality,
        )

    def transfer_to(self, target: ElectrolyteReference) -> TransferFreeEnergy:
        if not isinstance(target, ElectrolyteReference):
            raise _eos_error("invalid-reference-state", field="target")
        if (
            self.ion_component_ids != target.ion_component_ids
            or self.formula_stoichiometry != target.formula_stoichiometry
            or self.phase != target.phase
            or self.temperature.to("kelvin").magnitude
            != target.temperature.to("kelvin").magnitude
            or self.pressure.to("pascal").magnitude
            != target.pressure.to("pascal").magnitude
        ):
            raise _eos_error("invalid-reference-state", reason="incompatible-transfer")
        scale = (
            _MOLAR_GAS_CONSTANT
            * self.temperature.to("kelvin").magnitude
            * unit_registry.joule
            / unit_registry.mole
        )
        values = tuple(
            (target_value - source_value) * scale
            for source_value, target_value in zip(
                self.infinite_dilution_log_fugacity,
                target.infinite_dilution_log_fugacity,
                strict=True,
            )
        )
        return TransferFreeEnergy(
            ion_component_ids=self.ion_component_ids,
            values=values,  # type: ignore[arg-type]
        )


def build_reference_state(
    eos: EPCSAFT,
    *,
    temperature: object,
    pressure: object,
    solvent_mole_fractions: object,
    phase: str,
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
    if model.topology != "ionic-ssm-ds" or len(charges) != len(eos.component_ids):
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
        solvent_fractions[index] > 0.0
        and (model.molar_mass_kg_per_mol[index] or 0.0) <= 0.0
        for index in solvent_indices
    ):
        raise _eos_error("invalid-reference-state", reason="solvent-molar-mass-required")

    cation, anion = cations[0], anions[0]
    divisor = gcd(abs(charges[cation]), abs(charges[anion]))
    stoichiometry = (
        abs(charges[anion]) // divisor,
        abs(charges[cation]) // divisor,
    )
    ion_indices = (cation, anion)
    temperature_quantity = temperature_k * unit_registry.kelvin
    pressure_quantity = pressure_pa * unit_registry.pascal
    previous: tuple[float, float] | None = None
    current: tuple[float, float] | None = None
    convergence_error = math.inf
    for molality in _REFERENCE_MOLALITIES:
        composition = _formula_composition(
            solvent_fractions,
            model.molar_mass_kg_per_mol,
            ion_indices,
            stoichiometry,
            molality,
        )
        current = _ion_log_fugacity(
            eos,
            temperature=temperature_quantity,
            pressure=pressure_quantity,
            phase=phase,
            composition=composition,
            ion_indices=ion_indices,
        )
        if previous is not None:
            convergence_error = max(
                abs(value - old) for value, old in zip(current, previous, strict=True)
            )
        previous = current
    assert current is not None
    if convergence_error > _REFERENCE_CONVERGENCE_TOLERANCE:
        raise _eos_error(
            "reference-state-failed",
            reason="limit-not-converged",
            error=convergence_error,
        )

    return ElectrolyteReference(
        component_ids=eos.component_ids,
        ion_component_ids=(eos.component_ids[cation], eos.component_ids[anion]),
        formula_stoichiometry=stoichiometry,
        solvent_mole_fractions=solvent_fractions,
        temperature=temperature_quantity,
        pressure=pressure_quantity,
        phase=phase,
        infinite_dilution_log_fugacity=current,
        reference_molality=_REFERENCE_MOLALITIES[-1]
        * unit_registry.mole
        / unit_registry.kilogram,
        convergence_error=convergence_error,
        _eos=eos,
        _ion_indices=ion_indices,
        _molar_masses=model.molar_mass_kg_per_mol,
    )


__all__ = ("ElectrolyteReference", "IonicActivity", "TransferFreeEnergy")
