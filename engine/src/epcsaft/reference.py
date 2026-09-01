from __future__ import annotations

import hashlib
import json
import math
from collections.abc import Iterable
from dataclasses import dataclass, replace
from typing import TYPE_CHECKING, Any, Literal, SupportsFloat, cast

import pint

import epcsaft._core as _native_core

from .eos import (
    EosError,
    _eos_error,
    _state_magnitude,
)
from .parameters import ActiveParameterSet
from .units import unit_registry

if TYPE_CHECKING:
    from .eos import Mixture, State


_MOLAR_GAS_CONSTANT = 8.31446261815324
_SOURCE_REFERENCE_DERIVATIVES = {
    "temperature": 1,
    "pressure": 2,
    "composition": 4,
    "parameters": 8,
}
_SOURCE_REFERENCE_ERRORS = {
    1: "unsupported-source-reference-transfer",
    2: "invalid-source-composition",
    3: "component-order-mismatch",
    4: "missing-parameter-or-topology",
    5: "eos-domain-rejection",
    6: "reference-limit-convergence-failure",
    7: "derivative-unavailable",
    8: "nonfinite-source-reference-state",
    9: "source-reference-transport-identity-mismatch",
    10: "source-reference-evaluation-failure",
}


# EqExclude: transport and buffer mechanics
def _log_fugacity(
    mixture: Mixture,
    *,
    temperature: pint.Quantity[Any],
    pressure: pint.Quantity[Any],
    phase: Literal["vapor", "liquid"],
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
    return tuple(values)


@dataclass(frozen=True, slots=True)
class Activity:
    """Immutable electrolyte activity result in the complete model component order."""

    component_ids: tuple[str, ...]
    formula_unit_molality: pint.Quantity[float]
    mole_fractions: tuple[float, ...]
    ion_component_ids: tuple[str, str]
    reference_convention: Literal["molality-infinite-dilution"]
    single_ion: tuple[float, float]
    mean_mole_fraction: float
    mean_molality: float
    osmotic_coefficient: float | None


@dataclass(frozen=True, slots=True)
class TransferFreeEnergy:
    ion_component_ids: tuple[str, str]
    reference_convention: Literal["molality-infinite-dilution"]
    values: tuple[pint.Quantity[float], pint.Quantity[float]]


@dataclass(frozen=True, slots=True)
class SourceReferenceIonPair:
    """One caller-declared electroneutral formula-unit basis row."""

    component_ids: tuple[str, str]
    stoichiometry: tuple[int, int]

    def __post_init__(self) -> None:
        object.__setattr__(self, "component_ids", tuple(self.component_ids))
        object.__setattr__(self, "stoichiometry", tuple(self.stoichiometry))
        if (
            len(self.component_ids) != 2
            or len(self.stoichiometry) != 2
            or any(type(value) is not int for value in self.stoichiometry)
            or any(value <= 0 or value > 2**31 - 1 for value in self.stoichiometry)
        ):
            raise ValueError(
                "source-reference stoichiometry must contain two positive int32 values"
            )


@dataclass(frozen=True, slots=True)
class SourceReferencePureComponent:
    """One neutral component evaluated on its declared pure-phase branch."""

    component_id: str
    phase: Literal["liquid", "vapor"]

    def __post_init__(self) -> None:
        if not self.component_id or self.phase not in ("liquid", "vapor"):
            raise ValueError("pure-component source references require an identity and phase")


@dataclass(frozen=True, slots=True)
class SourceReferenceDeclaration:
    """Exact caller-owned source standard state in installed component order."""

    reference_state_id: str
    component_ids: tuple[str, ...]
    solvent_mole_fractions: tuple[float, ...]
    ion_pairs: tuple[SourceReferenceIonPair, ...]
    phase: Literal["liquid"]
    reference_convention: Literal["molality-infinite-dilution", "mole-fraction-mixed-reference"]
    activity_convention_id: str
    standard_molality_mol_per_kg: float
    reference_pressure_pa: float
    required_derivatives: tuple[str, ...]
    pure_components: tuple[SourceReferencePureComponent, ...] = ()

    def __post_init__(self) -> None:
        object.__setattr__(self, "component_ids", tuple(self.component_ids))
        object.__setattr__(self, "solvent_mole_fractions", tuple(self.solvent_mole_fractions))
        object.__setattr__(self, "ion_pairs", tuple(self.ion_pairs))
        object.__setattr__(self, "required_derivatives", tuple(self.required_derivatives))
        object.__setattr__(self, "pure_components", tuple(self.pure_components))


@dataclass(frozen=True, slots=True)
class SourceReferenceTransfer:
    """Charge-neutral source-to-EOS transfer and immutable EOS receipt."""

    status: Literal["ok"]
    domain_status: Literal["admitted-domain"]
    convergence_status: Literal["converged-limit"]
    reference_state_id: str
    activity_convention_id: str
    component_ids: tuple[str, ...]
    neutral_basis: tuple[tuple[float, ...], ...]
    log_fugacity_contractions: tuple[float, ...]
    activity_scale_log_contractions: tuple[float, ...]
    transfer_log_contractions: tuple[float, ...]
    source_composition: tuple[float, ...]
    derivative_availability: tuple[str, ...]
    pressure_derivatives_per_pa: tuple[float, ...]
    active_parameter_identities: tuple[str, ...]
    active_parameter_fingerprint: str | None
    parameter_derivatives: tuple[tuple[float, ...], ...]
    temperature_k: float
    pressure_pa: float
    source_reference_pressure_pa: float
    standard_molality_mol_per_kg: float
    solvent_molar_mass_kg_per_mol: float
    reference_limit_molality_mol_per_kg: float
    reference_convergence_error: float
    source_temperature_interval_k: tuple[float, float]
    source_pressure_interval_pa: tuple[float, float]
    candidate_temperature_interval_k: tuple[float, float]
    candidate_pressure_interval_pa: tuple[float, float]
    candidate_qualification: str
    candidate_source_fingerprint: str
    parameter_fingerprint: str
    topology_fingerprint: str
    component_order_fingerprint: str
    reference_state_fingerprint: str
    domain_fingerprint: str
    artifact_fingerprint: str
    helmholtz_basis_id: str
    eos_evaluation_count: int


def source_reference_transfer(
    model: Mixture,
    declaration: SourceReferenceDeclaration,
    *,
    T: object,
    P: object,
    active_parameters: ActiveParameterSet | None = None,
) -> SourceReferenceTransfer:
    """Evaluate a declared-solvent source reference through the installed transport."""
    if not isinstance(declaration, SourceReferenceDeclaration):
        raise _eos_error("unsupported-source-reference-transfer", field="declaration")
    temperature_k = _state_magnitude(T, "kelvin", "temperature")
    pressure_pa = _state_magnitude(P, "pascal", "pressure")
    component_ids = tuple(declaration.component_ids)
    try:
        solvent_mole_fractions = tuple(float(value) for value in declaration.solvent_mole_fractions)
        standard_molality_mol_per_kg = float(declaration.standard_molality_mol_per_kg)
        reference_pressure_pa = float(declaration.reference_pressure_pa)
    except (OverflowError, TypeError, ValueError) as exc:
        raise _eos_error("invalid-source-composition", reason=str(exc)) from exc
    if any(not math.isfinite(value) for value in solvent_mole_fractions) or any(
        not math.isfinite(value) for value in (standard_molality_mol_per_kg, reference_pressure_pa)
    ):
        raise _eos_error("nonfinite-source-reference-state")
    required_derivatives = tuple(declaration.required_derivatives)
    if active_parameters is not None and not isinstance(active_parameters, ActiveParameterSet):
        raise _eos_error("derivative-unavailable", derivative="parameters")
    if component_ids != tuple(model.component_ids):
        raise _eos_error("component-order-mismatch")
    invalid_derivative = next(
        (
            derivative
            for derivative in required_derivatives
            if derivative not in _SOURCE_REFERENCE_DERIVATIVES
        ),
        None,
    )
    if invalid_derivative is not None or len(set(required_derivatives)) != len(
        required_derivatives
    ):
        raise _eos_error("derivative-unavailable", derivative=invalid_derivative or "duplicate")
    if "parameters" in required_derivatives and active_parameters is None:
        raise _eos_error("derivative-unavailable", derivative="parameters")
    try:
        derivative_mask = sum(
            _SOURCE_REFERENCE_DERIVATIVES[derivative] for derivative in required_derivatives
        )
        payload = {
            "schema": "epcsaft-source-reference-declaration-v1",
            "reference_state_id": declaration.reference_state_id,
            "component_ids": component_ids,
            "solvent_mole_fractions": solvent_mole_fractions,
            "ion_pairs": tuple(
                (pair.component_ids, pair.stoichiometry) for pair in declaration.ion_pairs
            ),
            "phase": declaration.phase,
            "reference_convention": declaration.reference_convention,
            "activity_convention_id": declaration.activity_convention_id,
            "standard_molality_mol_per_kg": standard_molality_mol_per_kg,
            "reference_pressure_pa": reference_pressure_pa,
            "required_derivatives": required_derivatives,
            "pure_components": tuple(
                (component.component_id, component.phase)
                for component in declaration.pure_components
            ),
        }
        reference_fingerprint = (
            "sha256:"
            + hashlib.sha256(
                json.dumps(
                    payload,
                    allow_nan=False,
                    separators=(",", ":"),
                    sort_keys=True,
                ).encode("utf-8")
            ).hexdigest()
        )
        value = _native_core.source_reference_transfer(
            model._native_model(active_parameters),
            temperature_k,
            pressure_pa,
            component_ids,
            solvent_mole_fractions,
            tuple(
                tuple(component_ids.index(component_id) for component_id in pair.component_ids)
                for pair in declaration.ion_pairs
            ),
            tuple(pair.stoichiometry for pair in declaration.ion_pairs),
            tuple(
                component_ids.index(component.component_id)
                for component in declaration.pure_components
            ),
            tuple(component.phase for component in declaration.pure_components),
            declaration.phase,
            declaration.reference_convention,
            declaration.reference_state_id,
            declaration.activity_convention_id,
            reference_fingerprint,
            standard_molality_mol_per_kg,
            reference_pressure_pa,
            derivative_mask,
        )
    except EosError:
        raise
    except (OverflowError, TypeError, ValueError) as exc:
        if str(exc) == "unsupported-source-reference-transfer":
            raise _eos_error("unsupported-source-reference-transfer", reason=str(exc)) from exc
        raise _eos_error("source-reference-evaluation-failure", reason=str(exc)) from exc
    status = int(value["status"])
    if status != 0:
        raise _eos_error(
            _SOURCE_REFERENCE_ERRORS.get(status, "source-reference-evaluation-failure"),
            reason=value.get("error", "unknown failure"),
        )
    component_count = int(value["component_count"])
    row_count = int(value["neutral_basis_row_count"])
    flat_basis = tuple(float(item) for item in value["neutral_basis"])
    derivative_mask = int(value["derivative_availability"])
    derivative_availability = tuple(
        derivative
        for derivative in required_derivatives
        if derivative_mask & _SOURCE_REFERENCE_DERIVATIVES[derivative]
    )
    if derivative_availability != required_derivatives:
        raise _eos_error("derivative-unavailable", derivative="incomplete")
    active_parameter_identities = () if active_parameters is None else active_parameters.identities
    parameter_count = len(active_parameter_identities)
    flat_parameter_derivatives = tuple(float(item) for item in value["parameter_derivatives"])
    pressure_derivatives = tuple(float(item) for item in value["pressure_derivatives_per_pa"])
    expected_pressure_count = row_count if "pressure" in required_derivatives else 0
    expected_parameter_count = (
        row_count * parameter_count if "parameters" in required_derivatives else 0
    )
    if len(pressure_derivatives) != expected_pressure_count or (
        len(flat_parameter_derivatives) != expected_parameter_count
    ):
        raise _eos_error("source-reference-evaluation-failure", reason="parameter dimensions")
    return SourceReferenceTransfer(
        status="ok",
        domain_status="admitted-domain",
        convergence_status="converged-limit",
        reference_state_id=declaration.reference_state_id,
        activity_convention_id=declaration.activity_convention_id,
        component_ids=component_ids,
        neutral_basis=tuple(
            flat_basis[offset : offset + component_count]
            for offset in range(0, row_count * component_count, component_count)
        ),
        log_fugacity_contractions=tuple(value["log_fugacity_contractions"]),
        activity_scale_log_contractions=tuple(value["activity_scale_log_contractions"]),
        transfer_log_contractions=tuple(value["transfer_log_contractions"]),
        source_composition=tuple(value["source_composition"]),
        derivative_availability=derivative_availability,
        pressure_derivatives_per_pa=pressure_derivatives,
        active_parameter_identities=active_parameter_identities,
        active_parameter_fingerprint=(
            None if active_parameters is None else active_parameters.fingerprint
        ),
        parameter_derivatives=(
            tuple(
                flat_parameter_derivatives[offset : offset + parameter_count]
                for offset in range(0, row_count * parameter_count, parameter_count)
            )
            if "parameters" in required_derivatives
            else ()
        ),
        temperature_k=float(value["temperature_k"]),
        pressure_pa=float(value["pressure_pa"]),
        source_reference_pressure_pa=float(value["source_reference_pressure_pa"]),
        standard_molality_mol_per_kg=float(value["standard_molality_mol_per_kg"]),
        solvent_molar_mass_kg_per_mol=float(value["solvent_molar_mass_kg_per_mol"]),
        reference_limit_molality_mol_per_kg=float(value["reference_limit_molality_mol_per_kg"]),
        reference_convergence_error=float(value["reference_convergence_error"]),
        source_temperature_interval_k=(
            float(value["source_temperature_min_k"]),
            float(value["source_temperature_max_k"]),
        ),
        source_pressure_interval_pa=(
            float(value["source_pressure_min_pa"]),
            float(value["source_pressure_max_pa"]),
        ),
        candidate_temperature_interval_k=(
            float(value["candidate_temperature_min_k"]),
            float(value["candidate_temperature_max_k"]),
        ),
        candidate_pressure_interval_pa=(
            float(value["candidate_pressure_min_pa"]),
            float(value["candidate_pressure_max_pa"]),
        ),
        candidate_qualification=str(value["candidate_qualification"]),
        candidate_source_fingerprint=str(value["candidate_source_fingerprint"]),
        parameter_fingerprint=str(value["parameter_fingerprint"]),
        topology_fingerprint=str(value["topology_fingerprint"]),
        component_order_fingerprint=str(value["component_order_fingerprint"]),
        reference_state_fingerprint=str(value["reference_state_fingerprint"]),
        domain_fingerprint=str(value["domain_fingerprint"]),
        artifact_fingerprint=str(value["artifact_fingerprint"]),
        helmholtz_basis_id=str(value["helmholtz_basis_id"]),
        eos_evaluation_count=int(value["eos_evaluation_count"]),
    )


@dataclass(frozen=True, slots=True)
class ElectrolyteReference:
    component_ids: tuple[str, ...]
    ion_component_ids: tuple[str, str]
    formula_stoichiometry: tuple[int, int]
    reference_convention: Literal["molality-infinite-dilution"]
    solvent_mole_fractions: tuple[float, ...]
    temperature: pint.Quantity[float]
    pressure: pint.Quantity[float]
    phase: Literal["liquid"]
    infinite_dilution_log_fugacity_coefficient: tuple[float, float]
    reference_molality: pint.Quantity[float]
    convergence_error: float
    _mixture: Mixture
    _ion_indices: tuple[int, int]
    _molar_masses: tuple[float | None, ...]
    _osmotic_solvent_index: int | None
    _pure_solvent_log_fugacity: float | None

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
        return (
            self.infinite_dilution_log_fugacity_coefficient[0] * scale,
            self.infinite_dilution_log_fugacity_coefficient[1] * scale,
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
        osmotic_solvent_index = self._osmotic_solvent_index
        finite_log_fugacity = _log_fugacity(
            self._mixture,
            temperature=self.temperature,
            pressure=self.pressure,
            phase=self.phase,
            composition=composition,
            required_indices=self._ion_indices
            + ((osmotic_solvent_index,) if osmotic_solvent_index is not None else ()),
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
        single_ion = (math.exp(log_activity[0]), math.exp(log_activity[1]))
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
        if self._pure_solvent_log_fugacity is not None and osmotic_solvent_index is not None:
            log_solvent_activity = (
                math.log(composition[osmotic_solvent_index])
                + finite_log_fugacity[osmotic_solvent_index]
                - self._pure_solvent_log_fugacity
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
            reference_convention=self.reference_convention,
            single_ion=single_ion,
            mean_mole_fraction=mean_mole_fraction,
            mean_molality=mean_molality,
            osmotic_coefficient=osmotic_coefficient,
        )
        state = self._mixture.state(
            T=self.temperature,
            P=self.pressure,
            x=composition,
            phase=self.phase,
        )
        return replace(state, _activity=activity)

    def transfer_to(self, target: ElectrolyteReference) -> TransferFreeEnergy:
        if not isinstance(target, ElectrolyteReference):
            raise _eos_error("invalid-reference-state", field="target")
        if (
            self.ion_component_ids != target.ion_component_ids
            or self.formula_stoichiometry != target.formula_stoichiometry
            or self.reference_convention != target.reference_convention
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
        values = (
            (
                target.infinite_dilution_log_fugacity_coefficient[0]
                - self.infinite_dilution_log_fugacity_coefficient[0]
            )
            * scale,
            (
                target.infinite_dilution_log_fugacity_coefficient[1]
                - self.infinite_dilution_log_fugacity_coefficient[1]
            )
            * scale,
        )
        return TransferFreeEnergy(
            ion_component_ids=self.ion_component_ids,
            reference_convention=self.reference_convention,
            values=values,
        )


def build_reference_state(
    eos: Mixture,
    *,
    temperature: object,
    pressure: object,
    solvent_mole_fractions: object,
    phase: Literal["liquid"],
    ion_component_ids: tuple[str, str],
    formula_stoichiometry: tuple[int, int],
    reference_convention: Literal["molality-infinite-dilution"],
) -> ElectrolyteReference:
    temperature_k = _state_magnitude(temperature, "kelvin", "temperature")
    pressure_pa = _state_magnitude(pressure, "pascal", "pressure")
    if phase != "liquid":
        raise _eos_error("invalid-reference-state", field="phase")
    if reference_convention != "molality-infinite-dilution":
        raise _eos_error("invalid-reference-state", field="reference_convention")
    try:
        solvent_fractions = tuple(
            float(cast(SupportsFloat, value))
            for value in cast(Iterable[object], solvent_mole_fractions)
        )
    except (TypeError, ValueError) as exc:
        raise _eos_error("invalid-reference-state", field="solvent_mole_fractions") from exc

    model = eos._model
    charges = model.charge_number
    if model.topology not in {
        "ionic-born",
        "ionic-fully-dissociated-debye-huckel",
    } or len(charges) != len(eos.component_ids):
        raise _eos_error("invalid-reference-state", reason="electrolyte-required")
    cations = tuple(index for index, charge in enumerate(charges) if charge > 0)
    anions = tuple(index for index, charge in enumerate(charges) if charge < 0)
    if (
        len(ion_component_ids) != 2
        or ion_component_ids[0] not in eos.component_ids
        or ion_component_ids[1] not in eos.component_ids
    ):
        raise _eos_error("invalid-reference-state", field="ion_component_ids")
    cation = eos.component_ids.index(ion_component_ids[0])
    anion = eos.component_ids.index(ion_component_ids[1])
    if charges[cation] <= 0 or charges[anion] >= 0:
        raise _eos_error("invalid-reference-state", field="ion_component_ids")
    if (
        len(formula_stoichiometry) != 2
        or any(type(value) is not int or value <= 0 for value in formula_stoichiometry)
        or formula_stoichiometry[0] * charges[cation] + formula_stoichiometry[1] * charges[anion]
        != 0
    ):
        raise _eos_error("invalid-reference-state", field="formula_stoichiometry")
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

    stoichiometry = formula_stoichiometry
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
    if len(current) != 2:
        raise _eos_error("reference-state-failed", reason="invalid-ion-reference-count")
    ion_reference = (current[0], current[1])
    convergence_error = sequence.convergence_error

    active_solvents = tuple(index for index in solvent_indices if solvent_fractions[index] > 0.0)
    osmotic_solvent_index = active_solvents[0] if len(active_solvents) == 1 else None
    pure_solvent_log_fugacity = None
    if osmotic_solvent_index is not None:
        pure_solvent_log_fugacity = sequence.pure_solvent_log_fugacity_coefficient
        if pure_solvent_log_fugacity is None or not math.isfinite(pure_solvent_log_fugacity):
            raise _eos_error("reference-state-failed", reason="pure-solvent-fugacity-unavailable")

    return ElectrolyteReference(
        component_ids=eos.component_ids,
        ion_component_ids=(eos.component_ids[cation], eos.component_ids[anion]),
        formula_stoichiometry=stoichiometry,
        reference_convention=reference_convention,
        solvent_mole_fractions=solvent_fractions,
        temperature=temperature_quantity,
        pressure=pressure_quantity,
        phase=phase,
        infinite_dilution_log_fugacity_coefficient=ion_reference,
        reference_molality=sequence.reference_molality
        * unit_registry.mole
        / unit_registry.kilogram,
        convergence_error=convergence_error,
        _mixture=eos,
        _ion_indices=ion_indices,
        _molar_masses=model.molar_mass_kg_per_mol,
        _osmotic_solvent_index=osmotic_solvent_index,
        _pure_solvent_log_fugacity=pure_solvent_log_fugacity,
    )


__all__ = (
    "Activity",
    "ElectrolyteReference",
    "SourceReferenceDeclaration",
    "SourceReferenceIonPair",
    "SourceReferencePureComponent",
    "SourceReferenceTransfer",
    "TransferFreeEnergy",
    "source_reference_transfer",
)
