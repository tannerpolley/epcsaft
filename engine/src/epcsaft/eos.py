from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import math
from os import PathLike
from pathlib import Path
from typing import Any, TYPE_CHECKING, Iterable, Literal, SupportsFloat, TypeAlias, overload

import pint

from ._core import NativeModel as _NativeModel
from ._core import native_sdk_v1 as _native_sdk_v1
from ._resolve import _ResolvedModel, _resolve_model
from .parameters import BundlePurpose, ParameterBundle, ParameterSet
from .records import (
    ConstantCorrelation,
    ConstantPlusSumOfExponentialsCorrelation,
    SingleParameterRecord,
)
from .units import ParameterError, _error, _quantity, unit_registry

if TYPE_CHECKING:
    from .reference import ElectrolyteReference


Phase: TypeAlias = Literal["vapor", "liquid"]


_NATIVE_SDK_V1_PARAMETER_FINGERPRINTS = frozenset(
    {
        "sha256:5f836aa84935df70be2e5cffae51b178a7b797c2cee036e9ff47d8097ca94bbf",
        "sha256:288fbcaa1304881c16f64c3a784eeed19b75c58cca4558f92a21268e5e91258a",
        "sha256:9bfbc8d7789e51609945e61dbdf7a020decc8f9e31b408b0977724c7cb3e1551",
        "sha256:307fcb28d535b94782f3e3caf4012c0c8c0dc87ee4239d6c316de56553543286",
        "sha256:cb0c63530f96d51d5a24f48a2d02ae6567d163abdc373099ce096f7ff0a896f5",
        "sha256:7c637771bc9f717b8f47b44bb2a61044c3fe83084dca7c3c16102fba0989912d",
        "sha256:7551f1eee5903b66061cf7520f3bb59b169896ce372f3df3d48aa7ec778c39d4",
    }
)

_NATIVE_SDK_V1_ACTIVE_BORN_FINGERPRINTS = frozenset(
    {
        "sha256:1bb528ebe8f5612757e148608fc55821f9fb03737dbcec6d0bc4fffd0f4cbc4c",
        "sha256:7c637771bc9f717b8f47b44bb2a61044c3fe83084dca7c3c16102fba0989912d",
        "sha256:d29cef0c0f63034436d547d0aafa57934effe06783c8dffd89c94fa85e6940f4",
        "sha256:7551f1eee5903b66061cf7520f3bb59b169896ce372f3df3d48aa7ec778c39d4",
        "sha256:70ae04599dfa8338175e793bac6b9e4dfad37a9b96a568b5484dc87f104ef1a9",
    }
)

_NATIVE_SDK_V1_AQUEOUS_ACTIVE_KIJ_FINGERPRINTS = frozenset(
    {
        "sha256:1bb528ebe8f5612757e148608fc55821f9fb03737dbcec6d0bc4fffd0f4cbc4c",
        "sha256:7c637771bc9f717b8f47b44bb2a61044c3fe83084dca7c3c16102fba0989912d",
        "sha256:d29cef0c0f63034436d547d0aafa57934effe06783c8dffd89c94fa85e6940f4",
        "sha256:1258d8f76ddb4300d4983e70a2efa107b5cdc468748fe8760a04f82c4de38668",
        "sha256:b46a48ee5025e15a2d063b6f3ca4de54e29c33d429b94940913a9d225c78456b",
        "sha256:1d498278767a568e2120286f27f98a0b79268c92a68805d85ca30a21f14f3284",
    }
)

_NATIVE_SDK_V1_AQUEOUS_ACTIVE_SOLVATION_FACTOR_FINGERPRINTS = frozenset(
    {
        "sha256:b46a48ee5025e15a2d063b6f3ca4de54e29c33d429b94940913a9d225c78456b",
    }
)

_FIGIEL_STAGED_AQUEOUS_SELECTIONS = frozenset(
    ("water", cation, anion)
    for cation in ("lithium-cation", "sodium-cation", "potassium-cation")
    for anion in ("chloride-anion", "bromide-anion")
)

_FIGIEL_STAGED_ACTIVE_BORN_SELECTIONS = frozenset(
    {
        ("water", "lithium-cation", "chloride-anion"),
        ("water", "sodium-cation", "chloride-anion"),
        ("water", "potassium-cation", "chloride-anion"),
        ("water", "chloride-anion", "sodium-cation"),
        ("water", "bromide-anion", "sodium-cation"),
    }
)

_NATIVE_SDK_V1_NEUTRAL_REFERENCE_SELECTIONS = {
    "sha256:9f017007adbea2c943dc91c4a8df5271a5400392ae89796b74959aac913b920c": frozenset(
        frozenset(("water", cation, anion))
        for cation in ("lithium-cation", "sodium-cation", "potassium-cation")
        for anion in ("chloride-anion", "bromide-anion")
    ),
    "sha256:8ea2d7acbc179b82348fa938d081fe8d56f274bfe1f2167d10a507b7711aa1c2": (
        frozenset(
            {
                frozenset(
                    {
                        "water",
                        "hydronium-cation",
                        "hydroxide-anion",
                    }
                )
            }
        )
    ),
}

_NATIVE_SDK_V1_PARAMETER_BUNDLE_SELECTIONS = {
    "sha256:8ea2d7acbc179b82348fa938d081fe8d56f274bfe1f2167d10a507b7711aa1c2": frozenset(
        {
            "water",
            "hydronium-cation",
            "hydroxide-anion",
        }
    ),
    "sha256:03f9b1d88843a2d53da3420d96987ec4a2305b2814cc4946cdb25d1db67f5f72": frozenset(
        {
            "water",
            "ethanol",
            "isobutanol",
            "sodium-cation",
            "chloride-anion",
        }
    ),
    "sha256:b7a54dab41c61a71803f1d63bf25396188a2a1edc55172d5beddb99c3368c9fd": frozenset(
        {
            "water",
            "1-butanol",
            "sodium-cation",
            "potassium-cation",
            "chloride-anion",
        }
    ),
}


def _magnitude(
    record: SingleParameterRecord | ConstantCorrelation,
    unit: str,
) -> float:
    value = record.value if isinstance(record, SingleParameterRecord) else record.constant
    if unit != "dimensionless":
        value = value.to(unit).magnitude
    return float(value)


class EosError(ValueError):
    """One stable public error for invalid states and EOS solver failures."""

    __slots__ = ("_code", "_context")

    def __init__(self, code: str, context: Iterable[tuple[str, object]] = ()) -> None:
        self._code = code
        self._context = tuple((str(key), str(value)) for key, value in context)
        details = ", ".join(f"{key}={value}" for key, value in self._context)
        super().__init__(f"{code} ({details})" if details else code)

    @property
    def code(self) -> str:
        return self._code

    @property
    def context(self) -> tuple[tuple[str, str], ...]:
        return self._context


def _eos_error(code: str, **context: object) -> EosError:
    return EosError(code, context.items())


def _native_arguments(
    model: _ResolvedModel,
) -> tuple[object, ...]:
    sigma_constant: list[float] = []
    term_component: list[int] = []
    term_amplitude: list[float] = []
    term_exponent: list[float] = []
    for component, record in enumerate(model.segment_diameter):
        if isinstance(record, SingleParameterRecord):
            sigma_constant.append(_magnitude(record, "angstrom"))
        else:
            sigma_constant.append(float(record.constant.to("angstrom").magnitude))
            if isinstance(record, ConstantPlusSumOfExponentialsCorrelation):
                for term in record.terms:
                    term_component.append(component)
                    term_amplitude.append(float(term.amplitude.to("angstrom").magnitude))
                    term_exponent.append(
                        float(term.exponent_coefficient.to("1 / kelvin").magnitude)
                    )
            elif not isinstance(record, ConstantCorrelation):
                raise _error(
                    "invalid-value",
                    family="segment_diameter",
                    reason="unsupported-correlation",
                )

    component_position = {
        component_id: index for index, component_id in enumerate(model.component_ids)
    }
    site_position = {site_id: index for index, site_id in enumerate(model.site_ids)}
    sites = len(model.site_ids)
    association_active = [0] * (sites * sites)
    association_energy = [0.0] * (sites * sites)
    association_volume = [0.0] * (sites * sites)
    for pair in model.association_pairs:
        a = site_position[pair.endpoint_a]
        b = site_position[pair.endpoint_b]
        for index in (a * sites + b, b * sites + a):
            association_active[index] = 1
            association_energy[index] = pair.energy_over_k
            association_volume[index] = pair.volume

    components = len(model.component_ids)
    permittivity_model_codes = {
        "none": 0,
        "figiel-composition": 1,
        "ascani-mole-fraction": 2,
        "schreckenberg-density": 3,
        "zuber-ion-specific": 4,
        "rueben-perturbation": 5,
        "held-water": 6,
    }
    if model.topology in {
        "ionic-held-2008-debye-huckel",
        "ionic-ssm-ds",
        "ionic-ascani-original-born",
    }:
        molar_mass = [0.0 if value is None else value for value in model.molar_mass_kg_per_mol]
        charge_number = list(model.charge_number)
        relative_permittivity = [
            0.0 if record is None else _magnitude(record, "dimensionless")
            for record in model.relative_permittivity
        ]
        if model.topology == "ionic-ssm-ds":
            dielectric_suppression = (
                0.0
                if model.dielectric_ion_suppression_coefficient is None
                else model.dielectric_ion_suppression_coefficient
            )
            born_diameter = [
                0.0 if value is None else value for value in model.born_diameter_angstrom
            ]
            solvation_factor = [
                _magnitude(record, "dimensionless")
                for record in model.solvation_factor
                if record is not None
            ]
            ionic_permittivity = model.ionic_region_relative_permittivity
            electrolyte_model = 1
            assert ionic_permittivity is not None
            assert len(solvation_factor) == components
        elif model.topology == "ionic-held-2008-debye-huckel":
            dielectric_suppression = 0.0
            born_diameter = [0.0] * components
            solvation_factor = [0.0] * components
            ionic_permittivity = 0.0
            electrolyte_model = 3
        else:
            dielectric_suppression = 0.0
            born_diameter = [0.0] * components
            solvation_factor = [0.0] * components
            ionic_permittivity = 0.0
            electrolyte_model = 2
    else:
        molar_mass = [0.0] * components
        charge_number = [0] * components
        relative_permittivity = [0.0] * components
        dielectric_suppression = 0.0
        born_diameter = [0.0] * components
        solvation_factor = [0.0] * components
        ionic_permittivity = 0.0
        electrolyte_model = 0

    return (
        list(model.segment_count),
        sigma_constant,
        term_component,
        term_amplitude,
        term_exponent,
        list(model.dispersion_energy_over_k),
        list(model.k_ij),
        list(model.l_ij),
        [component_position[component_id] for component_id, _ in model.site_ids],
        list(model.site_multiplicity),
        association_active,
        association_energy,
        association_volume,
        molar_mass,
        charge_number,
        relative_permittivity,
        dielectric_suppression,
        born_diameter,
        solvation_factor,
        ionic_permittivity,
        electrolyte_model,
        list(model.schreckenberg_dielectric_volume_m3_per_mol)
        or [0.0] * components,
        list(model.schreckenberg_dielectric_temperature_k)
        or [0.0] * components,
        list(model.zuber_ion_suppression_coefficient) or [0.0] * components,
        list(model.rueben_dipole_scaling_debye2) or [0.0] * components,
        list(model.rueben_polarizability_scaling_angstrom3)
        or [0.0] * components,
        list(model.rueben_correlation_integral_parameter)
        or [0.0] * components,
        permittivity_model_codes[model.permittivity_model],
    )


def _figiel_staged_trial_bundle_admitted(
    parameters: ParameterSet, model: _ResolvedModel
) -> bool:
    admitted_component_order = (
        model.component_ids in _FIGIEL_STAGED_AQUEOUS_SELECTIONS
        and model.charge_number == (0, 1, -1)
    ) or (
        model.component_ids in _FIGIEL_STAGED_ACTIVE_BORN_SELECTIONS
        and model.charge_number == (0, -1, 1)
    )
    if (
        parameters.bundle_purpose != "user-provided"
        or model.topology != "ionic-ssm-ds"
        or not admitted_component_order
        or model.temperature_min_k != 298.15
        or model.temperature_max_k != 298.15
        or model.pressure_min_pa != 100000.0
        or model.pressure_max_pa != 100000.0
        or model.ion_mole_fraction_max != 0.38
    ):
        return False
    canonical = ParameterBundle.from_catalog(
        "figiel-2025-reference-electrolytes", version=1
    ).select(model.component_ids)
    if (
        parameters._source_domain_fingerprint
        != canonical._source_domain_fingerprint
    ):
        return False
    return all(
        (
            record.source_id == "figiel-yu-held-2025"
            and record.domain_id == "figiel-298-k-1-bar"
        )
        or (
            record.record_id == "water-molar-mass"
            and record.source_id == "nist-srd-69-water"
            and record.domain_id == "nist-water-identity"
        )
        for record in parameters.records
    )


@dataclass(frozen=True, slots=True, init=False)
class EPCSAFT:
    _component_ids: tuple[str, ...]
    _parameter_fingerprint: str
    _parameter_bundle_purpose: BundlePurpose
    _parameter_bundle_fingerprint: str
    _figiel_staged_trial_bundle_admitted: bool
    _model: _ResolvedModel
    _native: _NativeModel

    def __init__(self, parameters: ParameterSet) -> None:
        model = _resolve_model(parameters)
        native_arguments = _native_arguments(model)
        object.__setattr__(self, "_component_ids", parameters.component_ids)
        object.__setattr__(self, "_parameter_fingerprint", parameters.fingerprint)
        object.__setattr__(self, "_parameter_bundle_purpose", parameters.bundle_purpose)
        object.__setattr__(self, "_parameter_bundle_fingerprint", parameters.bundle_fingerprint)
        object.__setattr__(
            self,
            "_figiel_staged_trial_bundle_admitted",
            _figiel_staged_trial_bundle_admitted(parameters, model),
        )
        object.__setattr__(self, "_model", model)
        object.__setattr__(self, "_native", _NativeModel(*native_arguments))

    @property
    def component_ids(self) -> tuple[str, ...]:
        return self._component_ids

    @property
    def permittivity_model(self) -> str:
        """Return the source-named bulk-permittivity model used by this state."""
        return self._model.permittivity_model

    @property
    def parameter_fingerprint(self) -> str:
        return self._parameter_fingerprint

    @overload
    def evaluate(
        self,
        *,
        temperature: pint.Quantity[Any],
        molar_density: pint.Quantity[Any],
        pressure: None = None,
        mole_fractions: Iterable[SupportsFloat],
        phase: None = None,
    ) -> EosResult: ...

    @overload
    def evaluate(
        self,
        *,
        temperature: pint.Quantity[Any],
        molar_density: None = None,
        pressure: pint.Quantity[Any],
        mole_fractions: Iterable[SupportsFloat],
        phase: Phase | None = None,
    ) -> EosResult: ...

    def evaluate(
        self,
        *,
        temperature: object,
        molar_density: object | None = None,
        pressure: object | None = None,
        mole_fractions: object,
        phase: str | None = None,
    ) -> EosResult:
        """Evaluate one ordered composition at explicit density or specified pressure.

        Dimensioned inputs must be Pint quantities. Supply exactly one of
        ``molar_density`` and ``pressure``; ``phase`` is only valid for pressure
        closure. Invalid states and failed density/EOS evaluations raise ``EosError``.
        """
        if (molar_density is None) == (pressure is None):
            raise _eos_error("invalid-state", reason="specify-exactly-one-of-density-or-pressure")
        if phase not in (None, "vapor", "liquid"):
            raise _eos_error("invalid-state", field="phase")
        if molar_density is not None and phase is not None:
            raise _eos_error("invalid-state", field="phase")
        temperature_k = _state_magnitude(temperature, "kelvin", "temperature")
        try:
            fractions = tuple(float(value) for value in mole_fractions)  # type: ignore[union-attr]
        except (TypeError, ValueError) as exc:
            raise _eos_error(
                "invalid-state", field="mole_fractions", reason="not-numeric-sequence"
            ) from exc
        if len(fractions) != len(self._component_ids):
            raise _eos_error(
                "invalid-state",
                field="mole_fractions",
                reason="component-count",
                value=len(fractions),
                expected=len(self._component_ids),
            )
        if any(not math.isfinite(value) for value in fractions):
            raise _eos_error("invalid-state", field="mole_fractions", reason="nonfinite")
        if any(value < 0.0 for value in fractions):
            raise _eos_error("invalid-state", field="mole_fractions", reason="negative")
        fraction_sum = math.fsum(fractions)
        if abs(fraction_sum - 1.0) > 1.0e-12:
            raise _eos_error(
                "invalid-state",
                field="mole_fractions",
                reason="not-normalized",
                value=fraction_sum,
                expected=1.0,
                tolerance=1.0e-12,
            )

        if (
            self._model.temperature_min_k is not None
            and temperature_k < self._model.temperature_min_k
        ) or (
            self._model.temperature_max_k is not None
            and temperature_k > self._model.temperature_max_k
        ):
            context: dict[str, object] = {
                "field": "temperature",
                "reason": "outside-domain",
                "value": temperature_k,
            }
            if self._model.temperature_min_k is not None:
                context["minimum"] = self._model.temperature_min_k
            if self._model.temperature_max_k is not None:
                context["maximum"] = self._model.temperature_max_k
            context["unit"] = "kelvin"
            raise _eos_error("invalid-state", **context)

        if self._model.topology in {"ionic-ssm-ds", "ionic-ascani-original-born"}:
            ion_fraction = math.fsum(
                fraction
                for fraction, charge in zip(fractions, self._model.charge_number, strict=True)
                if charge != 0
            )
            if (
                self._model.ion_mole_fraction_max is not None
                and ion_fraction > self._model.ion_mole_fraction_max + 1.0e-12
            ):
                raise _eos_error(
                    "invalid-state",
                    field="mole_fractions",
                    reason="ion-fraction-limit",
                    value=ion_fraction,
                    maximum=self._model.ion_mole_fraction_max,
                )
            charge_balance = math.fsum(
                fraction * charge
                for fraction, charge in zip(fractions, self._model.charge_number, strict=True)
            )
            if abs(charge_balance) > 1.0e-12:
                raise _eos_error(
                    "invalid-state",
                    field="electroneutrality",
                    reason="charge-balance",
                    value=charge_balance,
                    tolerance=1.0e-12,
                )

        if molar_density is not None:
            density_molar = _state_magnitude(
                molar_density, "mole / meter ** 3", "molar_density"
            )
            try:
                native = self._native.evaluate(temperature_k, density_molar, fractions)
            except (RuntimeError, ValueError) as exc:
                raise _eos_error("evaluation-failed", reason=str(exc)) from exc
        else:
            pressure_pa = _state_magnitude(pressure, "pascal", "pressure")
            try:
                native = self._native.evaluate_pressure(
                    temperature_k, pressure_pa, fractions, phase or ""
                )
            except (RuntimeError, ValueError) as exc:
                raise _eos_error("evaluation-failed", reason=str(exc)) from exc
        diagnostics = native.density_diagnostics
        return EosResult(
            component_ids=self._component_ids,
            temperature=temperature_k * unit_registry.kelvin,
            mole_fractions=fractions,
            parameter_fingerprint=self._parameter_fingerprint,
            hard_chain=native.hard_chain,
            dispersion=native.dispersion,
            association=native.association,
            debye_huckel=native.debye_huckel,
            born=native.born_ssm_ds,
            bulk_relative_permittivity=None
            if self._model.permittivity_model == "none"
            else native.bulk_relative_permittivity,
            residual_helmholtz=native.residual_helmholtz,
            compressibility_factor=native.compressibility_factor,
            pressure=native.pressure * unit_registry.pascal,
            residual_chemical_potential_over_rt=tuple(native.residual_chemical_potential_over_rt),
            log_fugacity_coefficient=None
            if native.log_fugacity_coefficient is None
            else tuple(native.log_fugacity_coefficient),
            residual_entropy=native.residual_entropy
            * unit_registry.joule
            / unit_registry.mole
            / unit_registry.kelvin,
            residual_internal_energy=native.residual_internal_energy
            * unit_registry.joule
            / unit_registry.mole,
            residual_enthalpy=native.residual_enthalpy * unit_registry.joule / unit_registry.mole,
            residual_isochoric_heat_capacity=native.residual_isochoric_heat_capacity
            * unit_registry.joule
            / unit_registry.mole
            / unit_registry.kelvin,
            residual_isobaric_heat_capacity=None
            if native.residual_isobaric_heat_capacity is None
            else native.residual_isobaric_heat_capacity
            * unit_registry.joule
            / unit_registry.mole
            / unit_registry.kelvin,
            molar_density=native.molar_density * unit_registry.mole / unit_registry.meter**3,
            density_diagnostics=None
            if diagnostics is None
            else DensityDiagnostics(
                bracket_lower=diagnostics.bracket_lower
                * unit_registry.mole
                / unit_registry.meter**3,
                bracket_upper=diagnostics.bracket_upper
                * unit_registry.mole
                / unit_registry.meter**3,
                pressure_residual=diagnostics.pressure_residual * unit_registry.pascal,
                pressure_density_derivative=diagnostics.pressure_density_derivative
                * unit_registry.pascal
                * unit_registry.meter**3
                / unit_registry.mole,
                iterations=diagnostics.iterations,
                stable=diagnostics.stable,
                branch=diagnostics.branch,
            ),
        )

    def reference_state(
        self,
        *,
        temperature: object,
        pressure: object,
        solvent_mole_fractions: object,
        phase: Literal["liquid"],
    ) -> ElectrolyteReference:
        """Build one liquid electrolyte reference on the selected component order.

        The model must contain one source-complete cation/anion formula unit and
        solvent-only normalized reference fractions. Failures raise ``EosError``.
        """
        from .reference import build_reference_state

        return build_reference_state(
            self,
            temperature=temperature,
            pressure=pressure,
            solvent_mole_fractions=solvent_mole_fractions,
            phase=phase,
        )


def _native_topology_fingerprint(model: EPCSAFT) -> str:
    resolved = model._model
    charges = resolved.charge_number or (0,) * len(resolved.component_ids)
    neutral_phase = (
        resolved.topology == "neutral"
        and len(resolved.component_ids) in (1, 2)
        and not any(charges)
        and not resolved.site_ids
        and not resolved.association_pairs
    )
    if neutral_phase:
        pure = len(resolved.component_ids) == 1
        payload = {
            "schema": (
                "epcsaft-native-neutral-pure-topology-v1"
                if pure
                else "epcsaft-native-neutral-binary-topology-v1"
            ),
            "component_ids": resolved.component_ids,
            "topology": resolved.topology,
            "charge_number": charges,
            "site_ids": resolved.site_ids,
            "site_multiplicity": resolved.site_multiplicity,
            "association_pair_count": len(resolved.association_pairs),
            "helmholtz_terms": ("ideal", "hard-chain", "dispersion"),
        }
    else:
        payload = {
            "schema": "epcsaft-native-model-topology-v1",
            "component_ids": resolved.component_ids,
            "topology": resolved.topology,
            "charge_number": charges,
            "site_ids": resolved.site_ids,
            "site_multiplicity": resolved.site_multiplicity,
            "association_endpoints": tuple(
                (pair.endpoint_a, pair.endpoint_b)
                for pair in resolved.association_pairs
            ),
            "permittivity_model": resolved.permittivity_model,
        }
    encoded = json.dumps(
        payload, separators=(",", ":"), sort_keys=True
    ).encode("utf-8")
    return f"sha256:{hashlib.sha256(encoded).hexdigest()}"


def _native_sdk_configuration(model: EPCSAFT) -> dict[str, object]:
    if type(model) is not EPCSAFT:
        raise _eos_error("unsupported-native-sdk", reason="expected-exact-EPCSAFT-model")
    fingerprint = model._parameter_fingerprint
    bundle_selection = _NATIVE_SDK_V1_PARAMETER_BUNDLE_SELECTIONS.get(
        model._parameter_bundle_fingerprint
    )
    phase_admitted = (
        model._parameter_bundle_purpose == "user-provided"
        or fingerprint in _NATIVE_SDK_V1_PARAMETER_FINGERPRINTS
        or (
            bundle_selection is not None
            and len(model._component_ids) == len(bundle_selection)
            and frozenset(model._component_ids) == bundle_selection
        )
    )
    staged_trial_admitted = model._figiel_staged_trial_bundle_admitted
    active_born_admitted = (
        fingerprint in _NATIVE_SDK_V1_ACTIVE_BORN_FINGERPRINTS
        or (
            staged_trial_admitted
            and model._component_ids in _FIGIEL_STAGED_ACTIVE_BORN_SELECTIONS
        )
    )
    aqueous_kij_admitted = (
        fingerprint in _NATIVE_SDK_V1_AQUEOUS_ACTIVE_KIJ_FINGERPRINTS
        or (
            staged_trial_admitted
            and model._component_ids in _FIGIEL_STAGED_AQUEOUS_SELECTIONS
        )
    )
    figiel_dielectric_admitted = (
        len(model._component_ids) == 3
        and model._component_ids[0] in ("water", "methanol")
        and tuple(model._model.charge_number) == (0, 1, -1)
        and model._model.permittivity_model == "figiel-composition"
    )
    aqueous_solvation_factor_admitted = (
        fingerprint
        in _NATIVE_SDK_V1_AQUEOUS_ACTIVE_SOLVATION_FACTOR_FINGERPRINTS
        or (
            staged_trial_admitted
            and model._component_ids
            == ("water", "sodium-cation", "bromide-anion")
        )
    )
    catalog_neutral_reference_admitted = frozenset(model._component_ids) in (
        _NATIVE_SDK_V1_NEUTRAL_REFERENCE_SELECTIONS.get(
            model._parameter_bundle_fingerprint, frozenset()
        )
    )
    neutral_reference_admitted = catalog_neutral_reference_admitted or (
        phase_admitted
        and len(model._model.charge_number) == len(model._component_ids)
        and sum(
            component_id == "water" and charge == 0
            for component_id, charge in zip(
                model._component_ids, model._model.charge_number, strict=True
            )
        )
        == 1
        and any(charge > 0 for charge in model._model.charge_number)
        and any(charge < 0 for charge in model._model.charge_number)
    )
    if (
        not phase_admitted
        and not active_born_admitted
        and not aqueous_kij_admitted
        and not figiel_dielectric_admitted
        and not aqueous_solvation_factor_admitted
        and not neutral_reference_admitted
    ):
        raise _eos_error(
            "unsupported-native-sdk",
            reason="unadmitted-parameter-fingerprint",
        )
    minimum = math.nan if model._model.temperature_min_k is None else model._model.temperature_min_k
    maximum = math.nan if model._model.temperature_max_k is None else model._model.temperature_max_k
    pressure_minimum = (
        math.nan if model._model.pressure_min_pa is None else model._model.pressure_min_pa
    )
    pressure_maximum = (
        math.nan if model._model.pressure_max_pa is None else model._model.pressure_max_pa
    )
    maximum_ion_fraction = (
        math.nan
        if model._model.ion_mole_fraction_max is None
        else model._model.ion_mole_fraction_max
    )
    return {
        "parameter_fingerprint": fingerprint,
        "temperature_min_k": minimum,
        "temperature_max_k": maximum,
        "component_ids": model._component_ids,
        "ion_mole_fraction_max": maximum_ion_fraction,
        "phase_callbacks_admitted": phase_admitted,
        "active_born_admitted": active_born_admitted,
        "aqueous_kij_admitted": aqueous_kij_admitted,
        "figiel_dielectric_admitted": figiel_dielectric_admitted,
        "pressure_min_pa": pressure_minimum,
        "pressure_max_pa": pressure_maximum,
        "neutral_reference_admitted": neutral_reference_admitted,
        "unknown_applicability_admitted": model._parameter_bundle_purpose == "user-provided",
        "aqueous_solvation_factor_admitted": aqueous_solvation_factor_admitted,
        "staged_trial_bundle_admitted": staged_trial_admitted,
        "topology_fingerprint": _native_topology_fingerprint(model),
    }


def native_sdk(model: EPCSAFT) -> object:
    """Return the model-bound version-1 native provider capsule."""
    configuration = _native_sdk_configuration(model)
    try:
        return _native_sdk_v1(model._native, **configuration)
    except (RuntimeError, ValueError) as exc:
        raise _eos_error("unsupported-native-sdk", reason=str(exc)) from exc


def export_native_model(model: EPCSAFT, path: str | PathLike[str]) -> None:
    """Write one immutable resolved-model artifact for the installed native SDK."""
    configuration = _native_sdk_configuration(model)
    names = (
        "segment_count",
        "sigma_constant_angstrom",
        "sigma_term_component",
        "sigma_term_amplitude_angstrom",
        "sigma_term_exponent_per_k",
        "dispersion_energy_over_k_k",
        "k_ij",
        "l_ij",
        "site_component",
        "site_multiplicity",
        "association_active",
        "association_energy_over_k_k",
        "association_volume",
        "molar_mass_kg_per_mol",
        "charge_number",
        "relative_permittivity",
        "dielectric_ion_suppression_coefficient",
        "born_diameter_angstrom",
        "solvation_factor",
        "ionic_region_relative_permittivity",
        "electrolyte_model",
        "schreckenberg_dielectric_volume_m3_per_mol",
        "schreckenberg_dielectric_temperature_k",
        "zuber_ion_suppression_coefficient",
        "rueben_dipole_scaling_debye2",
        "rueben_polarizability_scaling_angstrom3",
        "rueben_correlation_integral_parameter",
        "permittivity_model",
    )
    native_parameters = dict(zip(names, _native_arguments(model._model), strict=True))
    for key in (
        "temperature_min_k",
        "temperature_max_k",
        "ion_mole_fraction_max",
        "pressure_min_pa",
        "pressure_max_pa",
    ):
        if math.isnan(float(configuration[key])):
            configuration[key] = None
    payload = {
        "schema": "epcsaft-resolved-native-model-v1",
        "native_parameters": native_parameters,
        "sdk": configuration,
    }
    target = Path(path)
    with target.open("x", encoding="utf-8") as stream:
        json.dump(payload, stream, allow_nan=False, separators=(",", ":"), sort_keys=True)
        stream.write("\n")


def _state_magnitude(value: object, unit: str, field: str) -> float:
    try:
        quantity = _quantity(value, unit, field=field)
    except ParameterError as exc:
        details = dict(exc.context)
        if exc.code == "invalid-value":
            rejected = details.get("value")
            if rejected is not None:
                try:
                    rejected_number = float(rejected)
                except ValueError:
                    pass
                else:
                    if not math.isfinite(rejected_number):
                        raise _eos_error(
                            "invalid-state",
                            field=field,
                            reason="nonfinite",
                            value=rejected,
                        ) from exc
            context = {"field": field, "reason": exc.code}
            for key in ("value", "value_type"):
                if key in details:
                    context[key] = details[key]
            raise _eos_error("invalid-state", **context) from exc
        raise _eos_error(
            "invalid-state",
            field=field,
            reason=exc.code,
            expected_unit=details.get("expected", unit),
        ) from exc
    magnitude = float(quantity.magnitude)
    if not math.isfinite(magnitude):
        raise _eos_error(
            "invalid-state",
            field=field,
            reason="nonfinite",
            value=magnitude,
        )
    if magnitude <= 0.0:
        raise _eos_error("invalid-state", field=field, reason="nonpositive", value=magnitude)
    return magnitude


@dataclass(frozen=True, slots=True)
class EosResult:
    """Immutable properties plus the resolved state and parameter identity."""

    component_ids: tuple[str, ...]
    temperature: pint.Quantity[float]
    mole_fractions: tuple[float, ...]
    parameter_fingerprint: str
    hard_chain: float
    dispersion: float
    association: float
    debye_huckel: float
    born: float
    bulk_relative_permittivity: float | None
    residual_helmholtz: float
    compressibility_factor: float
    pressure: pint.Quantity[float]
    residual_chemical_potential_over_rt: tuple[float, ...]
    log_fugacity_coefficient: tuple[float, ...] | None
    residual_entropy: pint.Quantity[float]
    residual_internal_energy: pint.Quantity[float]
    residual_enthalpy: pint.Quantity[float]
    residual_isochoric_heat_capacity: pint.Quantity[float]
    residual_isobaric_heat_capacity: pint.Quantity[float] | None
    molar_density: pint.Quantity[float]
    density_diagnostics: DensityDiagnostics | None


@dataclass(frozen=True, slots=True)
class DensityDiagnostics:
    bracket_lower: pint.Quantity[float]
    bracket_upper: pint.Quantity[float]
    pressure_residual: pint.Quantity[float]
    pressure_density_derivative: pint.Quantity[float]
    iterations: int
    stable: bool
    branch: str


__all__ = (
    "DensityDiagnostics",
    "EPCSAFT",
    "EosError",
    "EosResult",
    "Phase",
    "export_native_model",
    "native_sdk",
)
