import ctypes

from copy import deepcopy

import gc

import math

from pathlib import Path

import shutil

import subprocess

import time

import pytest

import epcsaft

import epcsaft._core as native_core
from epcsaft._core import native_sdk_v1 as _native_sdk_v1

from epcsaft import (
    Mixture,
    Parameters,
    export_native_model,
    native_sdk,
    unit_registry,
)

from tests._raw_cases import (
    ASSOCIATING_PURE,
    POLAR_BINARY,
    POLAR_PURE,
    electrolyte,
    neutral_binary,
    neutral_pure,
)


def three_site_association_case(
    values: tuple[float, float, float, float, float, float, float] = (
        1.4,
        3.2,
        220.0,
        2400.0,
        0.030,
        1800.0,
        0.025,
    ),
    *,
    site_order: tuple[str, str, str] = ("a", "b", "c"),
    multiplicities: tuple[int, int, int] = (1, 1, 2),
) -> dict[str, object]:
    segment_count, diameter, dispersion, energy_ab, volume_ab, energy_ac, volume_ac = values
    multiplicity = dict(zip(("a", "b", "c"), multiplicities, strict=True))
    return {
        "schema": "epcsaft.parameters",
        "schema_version": 1,
        "components": ("three-site-fluid",),
        "parameters": {
            "m": (segment_count,),
            "s": (diameter,),
            "e": (dispersion,),
            "k_ij": ((0.0,),),
            "sites": tuple(
                {
                    "component_id": "three-site-fluid",
                    "site_id": site_id,
                    "site_class": "generic",
                    "multiplicity": multiplicity[site_id],
                }
                for site_id in site_order
            ),
            "association": (
                {
                    "component_id_a": "three-site-fluid",
                    "site_id_a": "a",
                    "component_id_b": "three-site-fluid",
                    "site_id_b": "b",
                    "association_energy_over_k": energy_ab,
                    "association_volume": volume_ab,
                },
                {
                    "component_id_a": "three-site-fluid",
                    "site_id_a": "a",
                    "component_id_b": "three-site-fluid",
                    "site_id_b": "c",
                    "association_energy_over_k": energy_ac,
                    "association_volume": volume_ac,
                },
            ),
        },
        "options": {"permittivity_model": "none"},
        "validity": {"kind": "unknown"},
    }


CAPSULE_NAME = b"epcsaft.native_sdk.v1"

ABI_VERSION = 1

STATUS_OK = 0

STATUS_INVALID_ARGUMENT = 1

STATUS_UNSUPPORTED_MODEL = 2

STATUS_DOMAIN_ERROR = 3

STATUS_EVALUATION_ERROR = 4

STATUS_ABI_MISMATCH = 5

STATUS_CANCELLED = 6

STATUS_DEADLINE_EXCEEDED = 7


class PhaseBlockResult(ctypes.Structure):
    _fields_ = (
        ("struct_size", ctypes.c_uint32),
        ("status", ctypes.c_int32),
        ("helmholtz_over_rt_reference_amount", ctypes.c_double),
        ("gradient", ctypes.c_double * 2),
        ("hessian", ctypes.c_double * 4),
        ("third", ctypes.c_double * 8),
        ("pressure_pa", ctypes.c_double),
        ("chemical_potential_over_rt", ctypes.c_double),
        ("parameter_fingerprint", ctypes.c_char * 72),
        ("error", ctypes.c_char * 160),
    )


class ParameterizedPhaseBlockResult(ctypes.Structure):
    _fields_ = (
        ("struct_size", ctypes.c_uint32),
        ("status", ctypes.c_int32),
        ("helmholtz_over_rt_reference_amount", ctypes.c_double),
        ("gradient", ctypes.c_double * 5),
        ("hessian", ctypes.c_double * 25),
        ("pressure_pa", ctypes.c_double),
        ("chemical_potential_over_rt", ctypes.c_double),
        ("parameter_fingerprint", ctypes.c_char * 72),
        ("error", ctypes.c_char * 160),
    )


class MixturePhaseBlockResult(ctypes.Structure):
    _fields_ = (
        ("struct_size", ctypes.c_uint32),
        ("status", ctypes.c_int32),
        ("coordinate_count", ctypes.c_size_t),
        ("gradient_capacity", ctypes.c_size_t),
        ("hessian_capacity", ctypes.c_size_t),
        ("gradient", ctypes.POINTER(ctypes.c_double)),
        ("hessian", ctypes.POINTER(ctypes.c_double)),
        ("helmholtz_over_rt_reference_amount", ctypes.c_double),
        ("pressure_pa", ctypes.c_double),
        ("parameter_fingerprint", ctypes.c_char * 72),
        ("error", ctypes.c_char * 160),
    )


class AcceptedNativeSdkTable(ctypes.Structure):
    _fields_ = (
        ("abi_version", ctypes.c_uint32),
        ("table_size", ctypes.c_size_t),
        ("result_size", ctypes.c_size_t),
        ("model_context", ctypes.c_void_p),
        ("evaluate_pure_phase", ctypes.c_void_p),
        ("parameterized_result_size", ctypes.c_size_t),
        ("evaluate_pure_phase_parameters", ctypes.c_void_p),
    )


class FixedMixtureNativeSdkTable(ctypes.Structure):
    _fields_ = AcceptedNativeSdkTable._fields_ + (
        ("component_count", ctypes.c_size_t),
        ("mixture_result_size", ctypes.c_size_t),
        ("evaluate_mixture_phase", ctypes.c_void_p),
    )


class NativeSdkTable(ctypes.Structure):
    _fields_ = FixedMixtureNativeSdkTable._fields_ + (
        ("evaluate_mixture_phase_kij", ctypes.c_void_p),
    )


class ElectrolyteNativeSdkTable(ctypes.Structure):
    _fields_ = NativeSdkTable._fields_ + (
        ("component_ids", ctypes.POINTER(ctypes.c_char_p)),
        ("component_charges", ctypes.POINTER(ctypes.c_int32)),
        ("evaluate_electrolyte_phase", ctypes.c_void_p),
    )


class VolumeBoundsNativeSdkTable(ctypes.Structure):
    _fields_ = ElectrolyteNativeSdkTable._fields_ + (
        ("evaluate_molar_volume_bounds", ctypes.c_void_p),
    )


class PackingFractionNativeSdkTable(ctypes.Structure):
    _fields_ = VolumeBoundsNativeSdkTable._fields_ + (
        ("evaluate_packing_fraction", ctypes.c_void_p),
    )


class ApplicabilityNativeSdkTable(ctypes.Structure):
    _fields_ = PackingFractionNativeSdkTable._fields_ + (
        ("source_temperature_min_k", ctypes.c_double),
        ("source_temperature_max_k", ctypes.c_double),
        ("total_ion_mole_fraction_max", ctypes.c_double),
    )


class IonSolvationBornResult(ctypes.Structure):
    _fields_ = (
        ("struct_size", ctypes.c_uint32),
        ("status", ctypes.c_int32),
        ("solvation_gibbs_j_per_mol", ctypes.c_double),
        ("derivative_j_per_mol_per_angstrom", ctypes.c_double),
        ("reference_molality_mol_per_kg", ctypes.c_double),
        ("reference_convergence_error", ctypes.c_double),
        ("parameter_fingerprint", ctypes.c_char * 72),
        ("error", ctypes.c_char * 160),
    )


class ActiveBornNativeSdkTable(ctypes.Structure):
    _fields_ = ApplicabilityNativeSdkTable._fields_ + (
        ("ion_solvation_born_result_size", ctypes.c_size_t),
        ("evaluate_ion_solvation_born", ctypes.c_void_p),
    )


class AqueousMiacKijResult(ctypes.Structure):
    _fields_ = (
        ("struct_size", ctypes.c_uint32),
        ("status", ctypes.c_int32),
        ("log_mean_ionic_activity_coefficient_molality", ctypes.c_double),
        ("derivative", ctypes.c_double * 3),
        ("reference_molality_mol_per_kg", ctypes.c_double),
        ("reference_convergence_error", ctypes.c_double),
        ("reference_derivative_convergence_error", ctypes.c_double),
        ("parameter_fingerprint", ctypes.c_char * 72),
        ("error", ctypes.c_char * 160),
    )


class AqueousMiacKijNativeSdkTable(ctypes.Structure):
    _fields_ = ActiveBornNativeSdkTable._fields_ + (
        ("aqueous_miac_kij_result_size", ctypes.c_size_t),
        ("evaluate_aqueous_miac_kij", ctypes.c_void_p),
    )


class NeutralReferenceResult(ctypes.Structure):
    _fields_ = (
        ("struct_size", ctypes.c_uint32),
        ("status", ctypes.c_int32),
        ("component_count", ctypes.c_size_t),
        ("neutral_basis_row_count", ctypes.c_size_t),
        ("neutral_basis_capacity", ctypes.c_size_t),
        ("contraction_capacity", ctypes.c_size_t),
        ("reference_composition_capacity", ctypes.c_size_t),
        ("neutral_basis", ctypes.POINTER(ctypes.c_double)),
        ("log_fugacity_contractions", ctypes.POINTER(ctypes.c_double)),
        ("reference_composition", ctypes.POINTER(ctypes.c_double)),
        ("derivative_availability", ctypes.c_uint32),
        ("temperature_k", ctypes.c_double),
        ("pressure_pa", ctypes.c_double),
        ("solvent_molar_mass_kg_per_mol", ctypes.c_double),
        ("reference_amount_mol", ctypes.c_double),
        ("reference_number_density_mol_per_m3", ctypes.c_double),
        ("reference_molality_mol_per_kg", ctypes.c_double),
        ("reference_convergence_error", ctypes.c_double),
        ("parameter_fingerprint", ctypes.c_char * 72),
        ("helmholtz_basis_id", ctypes.c_char * 72),
        ("error", ctypes.c_char * 160),
    )


class NeutralReferenceNativeSdkTable(ctypes.Structure):
    _fields_ = AqueousMiacKijNativeSdkTable._fields_ + (
        ("neutral_reference_basis_row_count", ctypes.c_size_t),
        ("neutral_reference_result_size", ctypes.c_size_t),
        ("evaluate_neutral_reference", ctypes.c_void_p),
    )


class AqueousMiacSolvationFactorResult(ctypes.Structure):
    _fields_ = (
        ("struct_size", ctypes.c_uint32),
        ("status", ctypes.c_int32),
        ("log_mean_ionic_activity_coefficient_molality", ctypes.c_double),
        ("derivative", ctypes.c_double),
        ("reference_molality_mol_per_kg", ctypes.c_double),
        ("reference_convergence_error", ctypes.c_double),
        ("reference_derivative_convergence_error", ctypes.c_double),
        ("parameter_fingerprint", ctypes.c_char * 72),
        ("error", ctypes.c_char * 160),
    )


class AqueousMiacSolvationFactorNativeSdkTable(ctypes.Structure):
    _fields_ = NeutralReferenceNativeSdkTable._fields_ + (
        ("aqueous_miac_solvation_factor_result_size", ctypes.c_size_t),
        ("evaluate_aqueous_miac_solvation_factor", ctypes.c_void_p),
    )


class AqueousMiacBatchNativeSdkTable(ctypes.Structure):
    _fields_ = AqueousMiacSolvationFactorNativeSdkTable._fields_ + (
        ("evaluate_aqueous_miac_kij_batch", ctypes.c_void_p),
        ("evaluate_aqueous_miac_solvation_factor_batch", ctypes.c_void_p),
    )


class NativeEvaluationBudget(ctypes.Structure):
    _fields_ = (
        ("struct_size", ctypes.c_uint32),
        ("work_budget", ctypes.c_uint64),
        ("deadline_ns", ctypes.c_uint64),
    )


class AqueousMiacBoundedBatchNativeSdkTable(ctypes.Structure):
    _fields_ = AqueousMiacBatchNativeSdkTable._fields_ + (
        ("evaluation_budget_size", ctypes.c_size_t),
        ("evaluate_aqueous_miac_kij_batch_bounded", ctypes.c_void_p),
        ("evaluate_aqueous_miac_solvation_factor_batch_bounded", ctypes.c_void_p),
    )


class NativeCapabilityCoordinate(ctypes.Structure):
    _fields_ = (
        ("struct_size", ctypes.c_uint32),
        ("kind", ctypes.c_uint32),
        ("component_index", ctypes.c_int32),
        ("pair_component_index_a", ctypes.c_int32),
        ("pair_component_index_b", ctypes.c_int32),
        ("unit", ctypes.c_char_p),
    )


class NativeCapabilityDescriptor(ctypes.Structure):
    _fields_ = (
        ("struct_size", ctypes.c_uint32),
        ("schema_version", ctypes.c_uint32),
        ("capability", ctypes.c_uint32),
        ("parameter_family", ctypes.c_uint32),
        ("parameter_identity", ctypes.c_uint32),
        ("observation_contract", ctypes.c_uint32),
        ("model_domain", ctypes.c_uint32),
        ("tensor_layout", ctypes.c_uint32),
        ("derivative_order", ctypes.c_uint32),
        ("maturity", ctypes.c_uint32),
        ("authority_effect", ctypes.c_uint32),
        ("unsupported_status", ctypes.c_int32),
        ("domain_status", ctypes.c_int32),
        ("state_coordinate_count", ctypes.c_size_t),
        ("active_parameter_count", ctypes.c_size_t),
        ("coordinate_count", ctypes.c_size_t),
        ("coordinates", ctypes.POINTER(NativeCapabilityCoordinate)),
        ("component_count", ctypes.c_size_t),
        ("component_ids", ctypes.POINTER(ctypes.c_char_p)),
        ("temperature_min_k", ctypes.c_double),
        ("temperature_max_k", ctypes.c_double),
        ("parameter_fingerprint", ctypes.c_char * 72),
        ("topology_fingerprint", ctypes.c_char * 72),
        ("helmholtz_basis_id", ctypes.c_char * 72),
    )


class CapabilityNativeSdkTable(ctypes.Structure):
    _fields_ = AqueousMiacBoundedBatchNativeSdkTable._fields_ + (
        ("capability_count", ctypes.c_size_t),
        ("capabilities", ctypes.POINTER(NativeCapabilityDescriptor)),
    )


class LijNativeSdkTable(ctypes.Structure):
    _fields_ = CapabilityNativeSdkTable._fields_ + (
        ("evaluate_mixture_phase_lij", ctypes.c_void_p),
    )


class PureParameterNativeSdkTable(ctypes.Structure):
    _fields_ = LijNativeSdkTable._fields_ + (("evaluate_pure_phase_parameter", ctypes.c_void_p),)


class IonFractionSuppressionResult(ctypes.Structure):
    _fields_ = (
        ("struct_size", ctypes.c_uint32),
        ("status", ctypes.c_int32),
        ("bulk_relative_permittivity", ctypes.c_double),
        ("derivative", ctypes.c_double),
        ("salt_free_relative_permittivity", ctypes.c_double),
        ("parameter_fingerprint", ctypes.c_char * 72),
        ("error", ctypes.c_char * 160),
    )


class IonFractionSuppressionNativeSdkTable(ctypes.Structure):
    _fields_ = PureParameterNativeSdkTable._fields_ + (
        ("ion_fraction_suppression_result_size", ctypes.c_size_t),
        ("evaluate_ion_fraction_suppression", ctypes.c_void_p),
    )


class IonSolvationKijNativeSdkTable(ctypes.Structure):
    _fields_ = IonFractionSuppressionNativeSdkTable._fields_ + (
        ("ion_solvation_kij_result_size", ctypes.c_size_t),
        ("evaluate_ion_solvation_kij", ctypes.c_void_p),
    )


class IonSolvationIonicPermittivityNativeSdkTable(ctypes.Structure):
    _fields_ = IonSolvationKijNativeSdkTable._fields_ + (
        ("ion_solvation_ionic_permittivity_result_size", ctypes.c_size_t),
        ("evaluate_ion_solvation_ionic_permittivity", ctypes.c_void_p),
    )


class IonSolvationSolventPermittivityNativeSdkTable(ctypes.Structure):
    _fields_ = IonSolvationIonicPermittivityNativeSdkTable._fields_ + (
        ("ion_solvation_solvent_permittivity_result_size", ctypes.c_size_t),
        ("evaluate_ion_solvation_solvent_permittivity", ctypes.c_void_p),
    )


class ReservedAssociationNativeSdkTable(ctypes.Structure):
    _fields_ = IonSolvationSolventPermittivityNativeSdkTable._fields_ + (
        ("reserved_abi_size_0", ctypes.c_size_t),
        ("reserved_abi_slot_1", ctypes.c_void_p),
    )


class ElectrolytePhaseValueNativeSdkTable(ctypes.Structure):
    _fields_ = ReservedAssociationNativeSdkTable._fields_ + (
        ("electrolyte_phase_value_result_size", ctypes.c_size_t),
        ("evaluate_electrolyte_phase_value", ctypes.c_void_p),
    )


class ActiveParameterRequest(ctypes.Structure):
    _fields_ = (
        ("struct_size", ctypes.c_uint32),
        ("parameter_family", ctypes.c_uint32),
        ("parameter_identity", ctypes.c_uint32),
        ("component_index", ctypes.c_int32),
        ("pair_component_index_a", ctypes.c_int32),
        ("pair_component_index_b", ctypes.c_int32),
        ("value", ctypes.c_double),
        ("unit", ctypes.c_char_p),
    )


class NeutralReferenceDerivativeNativeSdkTable(ctypes.Structure):
    _fields_ = ElectrolytePhaseValueNativeSdkTable._fields_ + (
        ("neutral_reference_derivative_result_size", ctypes.c_size_t),
        ("evaluate_neutral_reference_derivatives", ctypes.c_void_p),
    )


class ReactingPhaseParameterResult(ctypes.Structure):
    _fields_ = (
        ("struct_size", ctypes.c_uint32),
        ("status", ctypes.c_int32),
        ("state_coordinate_count", ctypes.c_size_t),
        ("active_parameter_count", ctypes.c_size_t),
        ("state_gradient_capacity", ctypes.c_size_t),
        ("state_hessian_capacity", ctypes.c_size_t),
        ("parameter_gradient_capacity", ctypes.c_size_t),
        ("state_parameter_capacity", ctypes.c_size_t),
        ("packing_parameter_capacity", ctypes.c_size_t),
        ("molar_volume_bound_parameter_capacity", ctypes.c_size_t),
        ("pressure_parameter_capacity", ctypes.c_size_t),
        ("chemical_potential_parameter_capacity", ctypes.c_size_t),
        ("state_gradient", ctypes.POINTER(ctypes.c_double)),
        ("state_hessian", ctypes.POINTER(ctypes.c_double)),
        ("parameter_gradient", ctypes.POINTER(ctypes.c_double)),
        ("state_parameter_derivatives", ctypes.POINTER(ctypes.c_double)),
        ("packing_parameter_derivatives", ctypes.POINTER(ctypes.c_double)),
        (
            "molar_volume_bound_parameter_derivatives_m3_per_mol",
            ctypes.POINTER(ctypes.c_double),
        ),
        ("pressure_parameter_derivatives_pa", ctypes.POINTER(ctypes.c_double)),
        (
            "chemical_potential_parameter_derivatives_over_rt",
            ctypes.POINTER(ctypes.c_double),
        ),
        ("helmholtz_over_rt_reference_amount", ctypes.c_double),
        ("pressure_pa", ctypes.c_double),
        ("packing_fraction", ctypes.c_double),
        ("molar_volume_bounds_m3_per_mol", ctypes.c_double * 2),
        ("parameter_fingerprint", ctypes.c_char * 72),
        ("topology_fingerprint", ctypes.c_char * 72),
        ("helmholtz_basis_id", ctypes.c_char * 72),
        ("error", ctypes.c_char * 160),
        ("packing_state_gradient_capacity", ctypes.c_size_t),
        ("packing_state_hessian_capacity", ctypes.c_size_t),
        ("packing_state_gradient", ctypes.POINTER(ctypes.c_double)),
        ("packing_state_hessian", ctypes.POINTER(ctypes.c_double)),
    )


class ReactingSensitivityNativeSdkTable(ctypes.Structure):
    _fields_ = NeutralReferenceDerivativeNativeSdkTable._fields_ + (
        ("reacting_phase_parameter_result_size", ctypes.c_size_t),
        ("evaluate_reacting_phase_parameters", ctypes.c_void_p),
    )


class InversePackingGeometryResult(ctypes.Structure):
    _fields_ = (
        ("struct_size", ctypes.c_uint32),
        ("status", ctypes.c_int32),
        ("coordinate_count", ctypes.c_size_t),
        ("active_parameter_count", ctypes.c_size_t),
        ("gradient_capacity", ctypes.c_size_t),
        ("hessian_capacity", ctypes.c_size_t),
        ("gradient", ctypes.POINTER(ctypes.c_double)),
        ("hessian", ctypes.POINTER(ctypes.c_double)),
        ("volume_m3", ctypes.c_double),
        ("parameter_fingerprint", ctypes.c_char * 72),
        ("topology_fingerprint", ctypes.c_char * 72),
        ("error", ctypes.c_char * 160),
    )


class InversePackingNativeSdkTable(ctypes.Structure):
    _fields_ = ReactingSensitivityNativeSdkTable._fields_ + (
        ("inverse_packing_geometry_result_size", ctypes.c_size_t),
        ("evaluate_inverse_packing_geometry", ctypes.c_void_p),
    )


class FixedTopologyParameterSlot(ctypes.Structure):
    _fields_ = (
        ("struct_size", ctypes.c_uint32),
        ("parameter_family", ctypes.c_uint32),
        ("component_index", ctypes.c_int32),
        ("site_index_a", ctypes.c_int32),
        ("site_index_b", ctypes.c_int32),
        ("endpoint_component_index_a", ctypes.c_int32),
        ("endpoint_component_index_b", ctypes.c_int32),
        ("parameter_id", ctypes.c_char_p),
        ("unit", ctypes.c_char_p),
        ("lower_bound_exclusive", ctypes.c_double),
        ("upper_bound_exclusive", ctypes.c_double),
    )


class FixedTopologyAssociationDescriptor(ctypes.Structure):
    _fields_ = (
        ("struct_size", ctypes.c_uint32),
        ("schema_version", ctypes.c_uint32),
        ("observation_contract", ctypes.c_uint32),
        ("model_domain", ctypes.c_uint32),
        ("tensor_layout", ctypes.c_uint32),
        ("derivative_order", ctypes.c_uint32),
        ("maturity", ctypes.c_uint32),
        ("authority_effect", ctypes.c_uint32),
        ("unsupported_status", ctypes.c_int32),
        ("domain_status", ctypes.c_int32),
        ("state_coordinate_count", ctypes.c_size_t),
        ("slot_count", ctypes.c_size_t),
        ("slots", ctypes.POINTER(FixedTopologyParameterSlot)),
        ("component_count", ctypes.c_size_t),
        ("component_ids", ctypes.POINTER(ctypes.c_char_p)),
        ("site_count", ctypes.c_size_t),
        ("site_ids", ctypes.POINTER(ctypes.c_char_p)),
        ("site_component_indices", ctypes.POINTER(ctypes.c_int32)),
        ("site_multiplicities", ctypes.POINTER(ctypes.c_int32)),
        ("temperature_min_k", ctypes.c_double),
        ("temperature_max_k", ctypes.c_double),
        ("parameter_fingerprint", ctypes.c_char * 72),
        ("topology_fingerprint", ctypes.c_char * 72),
        ("artifact_fingerprint", ctypes.c_char * 72),
        ("helmholtz_basis_id", ctypes.c_char * 72),
        ("association_basis_id", ctypes.c_char * 72),
    )


class FixedTopologyParameterRequest(ctypes.Structure):
    _fields_ = (
        ("struct_size", ctypes.c_uint32),
        ("parameter_family", ctypes.c_uint32),
        ("component_index", ctypes.c_int32),
        ("site_index_a", ctypes.c_int32),
        ("site_index_b", ctypes.c_int32),
        ("unit", ctypes.c_char_p),
        ("value", ctypes.c_double),
    )


class FixedTopologyPhaseResult(ctypes.Structure):
    _fields_ = (
        ("struct_size", ctypes.c_uint32),
        ("status", ctypes.c_int32),
        ("coordinate_count", ctypes.c_size_t),
        ("active_parameter_count", ctypes.c_size_t),
        ("gradient_capacity", ctypes.c_size_t),
        ("hessian_capacity", ctypes.c_size_t),
        ("gradient", ctypes.POINTER(ctypes.c_double)),
        ("hessian", ctypes.POINTER(ctypes.c_double)),
        ("helmholtz_over_rt_reference_amount", ctypes.c_double),
        ("pressure_pa", ctypes.c_double),
        ("chemical_potential_over_rt", ctypes.c_double),
        ("parameter_fingerprint", ctypes.c_char * 72),
        ("topology_fingerprint", ctypes.c_char * 72),
        ("artifact_fingerprint", ctypes.c_char * 72),
        ("helmholtz_basis_id", ctypes.c_char * 72),
        ("error", ctypes.c_char * 160),
    )


class FixedTopologyNativeSdkTable(ctypes.Structure):
    _fields_ = InversePackingNativeSdkTable._fields_ + (
        ("reserved_abi_slot_2", ctypes.c_void_p),
        ("fixed_topology_association_result_size", ctypes.c_size_t),
        (
            "fixed_topology_association",
            ctypes.POINTER(FixedTopologyAssociationDescriptor),
        ),
        ("evaluate_fixed_topology_association", ctypes.c_void_p),
    )


class NativeSdkPrefix(ctypes.Structure):
    _fields_ = NativeSdkTable._fields_[:5]


FIXED_MIXTURE_PREFIX_SIZE = (
    FixedMixtureNativeSdkTable.evaluate_mixture_phase.offset + ctypes.sizeof(ctypes.c_void_p)
)

ACTIVE_KIJ_PREFIX_SIZE = NativeSdkTable.evaluate_mixture_phase_kij.offset + ctypes.sizeof(
    ctypes.c_void_p
)

ELECTROLYTE_PREFIX_SIZE = (
    ElectrolyteNativeSdkTable.evaluate_electrolyte_phase.offset + ctypes.sizeof(ctypes.c_void_p)
)

VOLUME_BOUNDS_PREFIX_SIZE = (
    VolumeBoundsNativeSdkTable.evaluate_molar_volume_bounds.offset + ctypes.sizeof(ctypes.c_void_p)
)

PACKING_FRACTION_PREFIX_SIZE = (
    PackingFractionNativeSdkTable.evaluate_packing_fraction.offset + ctypes.sizeof(ctypes.c_void_p)
)

APPLICABILITY_PREFIX_SIZE = (
    ApplicabilityNativeSdkTable.total_ion_mole_fraction_max.offset + ctypes.sizeof(ctypes.c_double)
)

ACTIVE_BORN_PREFIX_SIZE = (
    ActiveBornNativeSdkTable.evaluate_ion_solvation_born.offset + ctypes.sizeof(ctypes.c_void_p)
)

AQUEOUS_MIAC_KIJ_PREFIX_SIZE = (
    AqueousMiacKijNativeSdkTable.evaluate_aqueous_miac_kij.offset + ctypes.sizeof(ctypes.c_void_p)
)

NEUTRAL_REFERENCE_PREFIX_SIZE = (
    NeutralReferenceNativeSdkTable.evaluate_neutral_reference.offset
    + ctypes.sizeof(ctypes.c_void_p)
)

AQUEOUS_MIAC_BATCH_PREFIX_SIZE = (
    AqueousMiacBatchNativeSdkTable.evaluate_aqueous_miac_solvation_factor_batch.offset
    + ctypes.sizeof(ctypes.c_void_p)
)

AQUEOUS_MIAC_BOUNDED_BATCH_PREFIX_SIZE = (
    AqueousMiacBoundedBatchNativeSdkTable.evaluate_aqueous_miac_solvation_factor_batch_bounded.offset
    + ctypes.sizeof(ctypes.c_void_p)
)

CAPABILITY_PREFIX_SIZE = CapabilityNativeSdkTable.capabilities.offset + ctypes.sizeof(
    ctypes.c_void_p
)

ACTIVE_LIJ_PREFIX_SIZE = LijNativeSdkTable.evaluate_mixture_phase_lij.offset + ctypes.sizeof(
    ctypes.c_void_p
)

ION_FRACTION_SUPPRESSION_PREFIX_SIZE = (
    IonFractionSuppressionNativeSdkTable.evaluate_ion_fraction_suppression.offset
    + ctypes.sizeof(ctypes.c_void_p)
)

REACTING_SENSITIVITY_PREFIX_SIZE = (
    ReactingSensitivityNativeSdkTable.evaluate_reacting_phase_parameters.offset
    + ctypes.sizeof(ctypes.c_void_p)
)

INVERSE_PACKING_PREFIX_SIZE = (
    InversePackingNativeSdkTable.evaluate_inverse_packing_geometry.offset
    + ctypes.sizeof(ctypes.c_void_p)
)

FIXED_TOPOLOGY_ASSOCIATION_PREFIX_SIZE = (
    FixedTopologyNativeSdkTable.evaluate_fixed_topology_association.offset
    + ctypes.sizeof(ctypes.c_void_p)
)

EvaluatePurePhase = ctypes.CFUNCTYPE(
    ctypes.c_int,
    ctypes.c_void_p,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.POINTER(PhaseBlockResult),
)

EvaluateMixturePhase = ctypes.CFUNCTYPE(
    ctypes.c_int,
    ctypes.c_void_p,
    ctypes.c_double,
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_size_t,
    ctypes.c_double,
    ctypes.POINTER(MixturePhaseBlockResult),
)

EvaluateMixturePhaseKij = ctypes.CFUNCTYPE(
    ctypes.c_int,
    ctypes.c_void_p,
    ctypes.c_double,
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_size_t,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.POINTER(MixturePhaseBlockResult),
)

EvaluateMixturePhaseLij = EvaluateMixturePhaseKij

EvaluateReactingPhaseParameters = ctypes.CFUNCTYPE(
    ctypes.c_int,
    ctypes.c_void_p,
    ctypes.c_char_p,
    ctypes.c_char_p,
    ctypes.c_double,
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_size_t,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.POINTER(ActiveParameterRequest),
    ctypes.c_size_t,
    ctypes.POINTER(ReactingPhaseParameterResult),
)

EvaluateInversePackingGeometry = ctypes.CFUNCTYPE(
    ctypes.c_int,
    ctypes.c_void_p,
    ctypes.c_char_p,
    ctypes.c_char_p,
    ctypes.c_double,
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_size_t,
    ctypes.c_double,
    ctypes.POINTER(ActiveParameterRequest),
    ctypes.c_size_t,
    ctypes.POINTER(InversePackingGeometryResult),
)

EvaluateFixedTopologyAssociation = ctypes.CFUNCTYPE(
    ctypes.c_int,
    ctypes.c_void_p,
    ctypes.c_char_p,
    ctypes.c_char_p,
    ctypes.c_char_p,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.POINTER(FixedTopologyParameterRequest),
    ctypes.c_size_t,
    ctypes.POINTER(FixedTopologyPhaseResult),
)

EvaluatePackingFraction = ctypes.CFUNCTYPE(
    ctypes.c_int,
    ctypes.c_void_p,
    ctypes.c_char_p,
    ctypes.c_double,
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_size_t,
    ctypes.c_double,
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_size_t,
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_size_t,
)

EvaluateIonSolvationBorn = ctypes.CFUNCTYPE(
    ctypes.c_int,
    ctypes.c_void_p,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.POINTER(IonSolvationBornResult),
)

EvaluateAqueousMiacKij = ctypes.CFUNCTYPE(
    ctypes.c_int,
    ctypes.c_void_p,
    ctypes.c_char_p,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_size_t,
    ctypes.POINTER(AqueousMiacKijResult),
)

EvaluateNeutralReference = ctypes.CFUNCTYPE(
    ctypes.c_int,
    ctypes.c_void_p,
    ctypes.c_char_p,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.POINTER(NeutralReferenceResult),
)

EvaluateAqueousMiacSolvationFactor = ctypes.CFUNCTYPE(
    ctypes.c_int,
    ctypes.c_void_p,
    ctypes.c_char_p,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.POINTER(AqueousMiacSolvationFactorResult),
)

EvaluateIonFractionSuppression = ctypes.CFUNCTYPE(
    ctypes.c_int,
    ctypes.c_void_p,
    ctypes.c_char_p,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.POINTER(IonFractionSuppressionResult),
)

EvaluateAqueousMiacKijBatch = ctypes.CFUNCTYPE(
    ctypes.c_int,
    ctypes.c_void_p,
    ctypes.c_char_p,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_size_t,
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_size_t,
    ctypes.POINTER(AqueousMiacKijResult),
    ctypes.c_size_t,
)

EvaluateAqueousMiacSolvationFactorBatch = ctypes.CFUNCTYPE(
    ctypes.c_int,
    ctypes.c_void_p,
    ctypes.c_char_p,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_size_t,
    ctypes.c_double,
    ctypes.POINTER(AqueousMiacSolvationFactorResult),
    ctypes.c_size_t,
)

EvaluateAqueousMiacKijBatchBounded = ctypes.CFUNCTYPE(
    ctypes.c_int,
    ctypes.c_void_p,
    ctypes.c_char_p,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_size_t,
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_size_t,
    ctypes.POINTER(AqueousMiacKijResult),
    ctypes.c_size_t,
    ctypes.POINTER(NativeEvaluationBudget),
)

EvaluateAqueousMiacSolvationFactorBatchBounded = ctypes.CFUNCTYPE(
    ctypes.c_int,
    ctypes.c_void_p,
    ctypes.c_char_p,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_size_t,
    ctypes.c_double,
    ctypes.POINTER(AqueousMiacSolvationFactorResult),
    ctypes.c_size_t,
    ctypes.POINTER(NativeEvaluationBudget),
)


def _model(topology: str, component: str) -> Mixture:
    if topology == "associating":
        mapping = deepcopy(ASSOCIATING_PURE)
        mapping["components"] = (component,)
        for site in mapping["parameters"]["sites"]:
            site["component_id"] = component
        for pair in mapping["parameters"]["association"]:
            pair["component_id_a"] = component
            pair["component_id_b"] = component
        return Mixture(Parameters.from_dictionary(mapping))
    return Mixture(Parameters.from_dictionary(neutral_pure(component)))


def _mixture_model(component_ids: tuple[str, str] = ("methane", "ethane")) -> Mixture:
    return Mixture(Parameters.from_dictionary(neutral_binary(component_ids)))


def _electrolyte_model(
    component_ids: tuple[str, ...] = (
        "water",
        "sodium-cation",
        "chloride-anion",
    ),
) -> Mixture:
    return Mixture(Parameters.from_dictionary(electrolyte(component_ids)))


def _held_water_ionization_model(
    component_ids: tuple[str, str, str] = (
        "water",
        "hydronium-cation",
        "hydroxide-anion",
    ),
) -> Mixture:
    return Mixture(Parameters.from_dictionary(electrolyte(component_ids)))


def _table(model: Mixture) -> tuple[object, NativeSdkTable]:
    capsule = native_sdk(model)
    get_name = ctypes.pythonapi.PyCapsule_GetName
    get_name.argtypes = (ctypes.py_object,)
    get_name.restype = ctypes.c_char_p
    assert get_name(capsule) == CAPSULE_NAME

    get_pointer = ctypes.pythonapi.PyCapsule_GetPointer
    get_pointer.argtypes = (ctypes.py_object, ctypes.c_char_p)
    get_pointer.restype = ctypes.c_void_p
    address = get_pointer(capsule, CAPSULE_NAME)
    assert address
    fixed_table = ctypes.cast(address, ctypes.POINTER(FixedMixtureNativeSdkTable)).contents
    assert fixed_table.table_size >= ACTIVE_KIJ_PREFIX_SIZE
    return capsule, ctypes.cast(address, ctypes.POINTER(NativeSdkTable)).contents


def _capability_table(
    model: Mixture,
) -> tuple[object, CapabilityNativeSdkTable]:
    capsule = native_sdk(model)
    get_pointer = ctypes.pythonapi.PyCapsule_GetPointer
    get_pointer.argtypes = (ctypes.py_object, ctypes.c_char_p)
    get_pointer.restype = ctypes.c_void_p
    address = get_pointer(capsule, CAPSULE_NAME)
    assert address
    table = ctypes.cast(address, ctypes.POINTER(CapabilityNativeSdkTable)).contents
    assert table.table_size >= CAPABILITY_PREFIX_SIZE
    return capsule, table


def _ion_fraction_suppression_table(
    model: Mixture,
) -> tuple[object, IonFractionSuppressionNativeSdkTable]:
    capsule = native_sdk(model)
    get_pointer = ctypes.pythonapi.PyCapsule_GetPointer
    get_pointer.argtypes = (ctypes.py_object, ctypes.c_char_p)
    get_pointer.restype = ctypes.c_void_p
    address = get_pointer(capsule, CAPSULE_NAME)
    assert address
    table = ctypes.cast(address, ctypes.POINTER(IonFractionSuppressionNativeSdkTable)).contents
    assert table.table_size >= ctypes.sizeof(IonFractionSuppressionNativeSdkTable)
    return capsule, table


def _lij_table(model: Mixture) -> tuple[object, LijNativeSdkTable]:
    capsule = native_sdk(model)
    get_pointer = ctypes.pythonapi.PyCapsule_GetPointer
    get_pointer.argtypes = (ctypes.py_object, ctypes.c_char_p)
    get_pointer.restype = ctypes.c_void_p
    address = get_pointer(capsule, CAPSULE_NAME)
    assert address
    table = ctypes.cast(address, ctypes.POINTER(LijNativeSdkTable)).contents
    assert table.table_size >= ACTIVE_LIJ_PREFIX_SIZE
    return capsule, table


def _electrolyte_table(
    model: Mixture,
) -> tuple[object, ElectrolyteNativeSdkTable]:
    capsule = native_sdk(model)
    get_pointer = ctypes.pythonapi.PyCapsule_GetPointer
    get_pointer.argtypes = (ctypes.py_object, ctypes.c_char_p)
    get_pointer.restype = ctypes.c_void_p
    address = get_pointer(capsule, CAPSULE_NAME)
    assert address
    table = ctypes.cast(address, ctypes.POINTER(ElectrolyteNativeSdkTable)).contents
    assert table.table_size >= ELECTROLYTE_PREFIX_SIZE
    return capsule, table


def _reacting_sensitivity_table(
    model: Mixture,
) -> tuple[object, ReactingSensitivityNativeSdkTable]:
    capsule = native_sdk(model)
    get_pointer = ctypes.pythonapi.PyCapsule_GetPointer
    get_pointer.argtypes = (ctypes.py_object, ctypes.c_char_p)
    get_pointer.restype = ctypes.c_void_p
    address = get_pointer(capsule, CAPSULE_NAME)
    assert address
    table = ctypes.cast(address, ctypes.POINTER(ReactingSensitivityNativeSdkTable)).contents
    assert table.table_size >= REACTING_SENSITIVITY_PREFIX_SIZE
    return capsule, table


def _inverse_packing_table(
    model: Mixture,
) -> tuple[object, InversePackingNativeSdkTable]:
    capsule = native_sdk(model)
    get_pointer = ctypes.pythonapi.PyCapsule_GetPointer
    get_pointer.argtypes = (ctypes.py_object, ctypes.c_char_p)
    get_pointer.restype = ctypes.c_void_p
    address = get_pointer(capsule, CAPSULE_NAME)
    assert address
    prefix = ctypes.cast(address, ctypes.POINTER(ReactingSensitivityNativeSdkTable)).contents
    assert prefix.table_size >= INVERSE_PACKING_PREFIX_SIZE
    table = ctypes.cast(address, ctypes.POINTER(InversePackingNativeSdkTable)).contents
    return capsule, table


def _fixed_topology_table(
    model: Mixture,
) -> tuple[object, FixedTopologyNativeSdkTable]:
    capsule = native_sdk(model)
    get_pointer = ctypes.pythonapi.PyCapsule_GetPointer
    get_pointer.argtypes = (ctypes.py_object, ctypes.c_char_p)
    get_pointer.restype = ctypes.c_void_p
    address = get_pointer(capsule, CAPSULE_NAME)
    assert address
    table = ctypes.cast(address, ctypes.POINTER(FixedTopologyNativeSdkTable)).contents
    assert table.table_size >= FIXED_TOPOLOGY_ASSOCIATION_PREFIX_SIZE
    return capsule, table


def _packing_fraction_table(
    model: Mixture,
) -> tuple[object, PackingFractionNativeSdkTable]:
    capsule = native_sdk(model)
    get_pointer = ctypes.pythonapi.PyCapsule_GetPointer
    get_pointer.argtypes = (ctypes.py_object, ctypes.c_char_p)
    get_pointer.restype = ctypes.c_void_p
    address = get_pointer(capsule, CAPSULE_NAME)
    assert address
    table = ctypes.cast(address, ctypes.POINTER(PackingFractionNativeSdkTable)).contents
    assert table.table_size >= PACKING_FRACTION_PREFIX_SIZE
    return capsule, table


def _active_born_table(model: Mixture) -> tuple[object, ActiveBornNativeSdkTable]:
    capsule = native_sdk(model)
    get_pointer = ctypes.pythonapi.PyCapsule_GetPointer
    get_pointer.argtypes = (ctypes.py_object, ctypes.c_char_p)
    get_pointer.restype = ctypes.c_void_p
    address = get_pointer(capsule, CAPSULE_NAME)
    assert address
    table = ctypes.cast(address, ctypes.POINTER(ActiveBornNativeSdkTable)).contents
    assert table.table_size >= ACTIVE_BORN_PREFIX_SIZE
    return capsule, table


def _neutral_reference_table(
    model: Mixture,
) -> tuple[object, NeutralReferenceNativeSdkTable]:
    capsule = native_sdk(model)
    get_pointer = ctypes.pythonapi.PyCapsule_GetPointer
    get_pointer.argtypes = (ctypes.py_object, ctypes.c_char_p)
    get_pointer.restype = ctypes.c_void_p
    address = get_pointer(capsule, CAPSULE_NAME)
    assert address
    table = ctypes.cast(address, ctypes.POINTER(NeutralReferenceNativeSdkTable)).contents
    assert table.table_size >= NEUTRAL_REFERENCE_PREFIX_SIZE
    return capsule, table


def _aqueous_miac_batch_table(
    model: Mixture,
) -> tuple[object, AqueousMiacBatchNativeSdkTable]:
    capsule = native_sdk(model)
    get_pointer = ctypes.pythonapi.PyCapsule_GetPointer
    get_pointer.argtypes = (ctypes.py_object, ctypes.c_char_p)
    get_pointer.restype = ctypes.c_void_p
    address = get_pointer(capsule, CAPSULE_NAME)
    assert address
    table = ctypes.cast(address, ctypes.POINTER(AqueousMiacBatchNativeSdkTable)).contents
    assert table.table_size >= AQUEOUS_MIAC_BATCH_PREFIX_SIZE
    return capsule, table


def _aqueous_miac_bounded_batch_table(
    model: Mixture,
) -> tuple[object, AqueousMiacBoundedBatchNativeSdkTable]:
    capsule = native_sdk(model)
    get_pointer = ctypes.pythonapi.PyCapsule_GetPointer
    get_pointer.argtypes = (ctypes.py_object, ctypes.c_char_p)
    get_pointer.restype = ctypes.c_void_p
    address = get_pointer(capsule, CAPSULE_NAME)
    assert address
    table = ctypes.cast(address, ctypes.POINTER(AqueousMiacBoundedBatchNativeSdkTable)).contents
    assert table.table_size >= AQUEOUS_MIAC_BOUNDED_BATCH_PREFIX_SIZE
    return capsule, table


def _evaluate_aqueous_miac_solvation_factor(
    table: AqueousMiacSolvationFactorNativeSdkTable,
    fingerprint: str,
    molality: float,
    solvation_factor: float,
) -> AqueousMiacSolvationFactorResult:
    result = AqueousMiacSolvationFactorResult(
        struct_size=ctypes.sizeof(AqueousMiacSolvationFactorResult)
    )
    status = EvaluateAqueousMiacSolvationFactor(table.evaluate_aqueous_miac_solvation_factor)(
        table.model_context,
        fingerprint.encode(),
        298.15,
        100000.0,
        molality,
        solvation_factor,
        ctypes.byref(result),
    )
    assert status == result.status
    return result


def _evaluate_aqueous_miac_solvation_factor_batch(
    table: AqueousMiacBatchNativeSdkTable,
    fingerprint: str,
    molalities: tuple[float, ...],
    solvation_factor: float,
) -> tuple[AqueousMiacSolvationFactorResult, ...]:
    molality_values = (ctypes.c_double * len(molalities))(*molalities)
    results = (AqueousMiacSolvationFactorResult * len(molalities))(
        *(
            AqueousMiacSolvationFactorResult(
                struct_size=ctypes.sizeof(AqueousMiacSolvationFactorResult)
            )
            for _ in molalities
        )
    )
    status = EvaluateAqueousMiacSolvationFactorBatch(
        table.evaluate_aqueous_miac_solvation_factor_batch
    )(
        table.model_context,
        fingerprint.encode(),
        298.15,
        100000.0,
        molality_values,
        len(molalities),
        solvation_factor,
        results,
        len(molalities),
    )
    assert all(result.status == status for result in results)
    return tuple(results)


def _evaluate_neutral_reference(
    table: NeutralReferenceNativeSdkTable,
    fingerprint: str,
    *,
    temperature_k: float = 298.15,
    pressure_pa: float = 100000.0,
) -> tuple[NeutralReferenceResult, tuple[float, ...], tuple[float, ...], tuple[float, ...]]:
    component_count = table.component_count
    row_count = table.neutral_reference_basis_row_count
    basis = (ctypes.c_double * (component_count * row_count))()
    contractions = (ctypes.c_double * row_count)()
    reference_composition = (ctypes.c_double * component_count)()
    result = NeutralReferenceResult(
        struct_size=ctypes.sizeof(NeutralReferenceResult),
        component_count=component_count,
        neutral_basis_row_count=row_count,
        neutral_basis_capacity=component_count * row_count,
        contraction_capacity=row_count,
        reference_composition_capacity=component_count,
        neutral_basis=basis,
        log_fugacity_contractions=contractions,
        reference_composition=reference_composition,
    )
    status = EvaluateNeutralReference(table.evaluate_neutral_reference)(
        table.model_context,
        fingerprint.encode(),
        temperature_k,
        pressure_pa,
        ctypes.byref(result),
    )
    assert status == result.status
    return result, tuple(basis), tuple(contractions), tuple(reference_composition)


def _active_request(
    family: int,
    component: int,
    value: float,
    unit: bytes,
) -> ActiveParameterRequest:
    return ActiveParameterRequest(
        struct_size=ctypes.sizeof(ActiveParameterRequest),
        parameter_family=family,
        parameter_identity=2,
        component_index=component,
        pair_component_index_a=-1,
        pair_component_index_b=-1,
        value=value,
        unit=unit,
    )


def _evaluate_reacting_phase_parameters(
    table: ReactingSensitivityNativeSdkTable,
    fingerprint: str,
    topology_fingerprint: str,
    amounts: tuple[float, ...],
    volume_m3: float,
    requests: tuple[ActiveParameterRequest, ...],
    *,
    temperature_k: float = 298.15,
    packing_state_gradient_capacity: int | None = None,
    packing_state_hessian_capacity: int | None = None,
    provide_packing_state_gradient: bool = True,
    provide_packing_state_hessian: bool = True,
) -> tuple[
    ReactingPhaseParameterResult,
    tuple[float, ...],
    tuple[float, ...],
    tuple[float, ...],
    tuple[float, ...],
    tuple[float, ...],
    tuple[float, ...],
    tuple[float, ...],
    tuple[float, ...],
    tuple[float, ...],
    tuple[float, ...],
]:
    state_count = len(amounts) + 1
    parameter_count = len(requests)
    amounts_buffer = (ctypes.c_double * len(amounts))(*amounts)
    request_buffer = (ActiveParameterRequest * parameter_count)(*requests)
    state_gradient = (ctypes.c_double * state_count)(*(21.0 for _ in range(state_count)))
    state_hessian = (ctypes.c_double * (state_count * state_count))(
        *(22.0 for _ in range(state_count * state_count))
    )
    parameter_gradient = (ctypes.c_double * parameter_count)(
        *(23.0 for _ in range(parameter_count))
    )
    mixed = (ctypes.c_double * (state_count * parameter_count))(
        *(24.0 for _ in range(state_count * parameter_count))
    )
    packing = (ctypes.c_double * parameter_count)(*(25.0 for _ in range(parameter_count)))
    packing_state_gradient = (ctypes.c_double * state_count)(*(29.0 for _ in range(state_count)))
    packing_state_hessian = (ctypes.c_double * (state_count * state_count))(
        *(30.0 for _ in range(state_count * state_count))
    )
    volume_bounds = (ctypes.c_double * (2 * parameter_count))(
        *(26.0 for _ in range(2 * parameter_count))
    )
    pressure = (ctypes.c_double * parameter_count)(*(27.0 for _ in range(parameter_count)))
    chemical = (ctypes.c_double * (len(amounts) * parameter_count))(
        *(28.0 for _ in range(len(amounts) * parameter_count))
    )
    result = ReactingPhaseParameterResult(
        struct_size=ctypes.sizeof(ReactingPhaseParameterResult),
        state_coordinate_count=state_count,
        active_parameter_count=parameter_count,
        state_gradient_capacity=state_count,
        state_hessian_capacity=state_count * state_count,
        parameter_gradient_capacity=parameter_count,
        state_parameter_capacity=state_count * parameter_count,
        packing_parameter_capacity=parameter_count,
        molar_volume_bound_parameter_capacity=2 * parameter_count,
        pressure_parameter_capacity=parameter_count,
        chemical_potential_parameter_capacity=len(amounts) * parameter_count,
        state_gradient=state_gradient,
        state_hessian=state_hessian,
        parameter_gradient=parameter_gradient,
        state_parameter_derivatives=mixed,
        packing_parameter_derivatives=packing,
        molar_volume_bound_parameter_derivatives_m3_per_mol=volume_bounds,
        pressure_parameter_derivatives_pa=pressure,
        chemical_potential_parameter_derivatives_over_rt=chemical,
        packing_state_gradient_capacity=(
            state_count
            if packing_state_gradient_capacity is None
            else packing_state_gradient_capacity
        ),
        packing_state_hessian_capacity=(
            state_count * state_count
            if packing_state_hessian_capacity is None
            else packing_state_hessian_capacity
        ),
        packing_state_gradient=(packing_state_gradient if provide_packing_state_gradient else None),
        packing_state_hessian=(packing_state_hessian if provide_packing_state_hessian else None),
    )
    status = EvaluateReactingPhaseParameters(table.evaluate_reacting_phase_parameters)(
        table.model_context,
        fingerprint.encode(),
        topology_fingerprint.encode(),
        temperature_k,
        amounts_buffer,
        len(amounts),
        volume_m3,
        0.01,
        0.5,
        request_buffer,
        parameter_count,
        ctypes.byref(result),
    )
    assert status == result.status
    return (
        result,
        tuple(state_gradient),
        tuple(state_hessian),
        tuple(parameter_gradient),
        tuple(mixed),
        tuple(packing),
        tuple(packing_state_gradient),
        tuple(packing_state_hessian),
        tuple(volume_bounds),
        tuple(pressure),
        tuple(chemical),
    )


def _evaluate_inverse_packing_geometry(
    table: InversePackingNativeSdkTable,
    fingerprint: str,
    topology_fingerprint: str,
    amounts: tuple[float, ...],
    log_packing_fraction: float,
    requests: tuple[ActiveParameterRequest, ...],
    *,
    gradient_capacity: int | None = None,
    hessian_capacity: int | None = None,
    provide_gradient: bool = True,
    provide_hessian: bool = True,
) -> tuple[
    InversePackingGeometryResult,
    tuple[float, ...],
    tuple[float, ...],
]:
    coordinate_count = len(amounts) + 1 + len(requests)
    amounts_buffer = (ctypes.c_double * len(amounts))(*amounts)
    request_buffer = (ActiveParameterRequest * len(requests))(*requests) if requests else None
    gradient = (ctypes.c_double * coordinate_count)(*(31.0 for _ in range(coordinate_count)))
    hessian = (ctypes.c_double * (coordinate_count * coordinate_count))(
        *(32.0 for _ in range(coordinate_count * coordinate_count))
    )
    result = InversePackingGeometryResult(
        struct_size=ctypes.sizeof(InversePackingGeometryResult),
        coordinate_count=coordinate_count,
        active_parameter_count=len(requests),
        gradient_capacity=(coordinate_count if gradient_capacity is None else gradient_capacity),
        hessian_capacity=(
            coordinate_count * coordinate_count if hessian_capacity is None else hessian_capacity
        ),
        gradient=gradient if provide_gradient else None,
        hessian=hessian if provide_hessian else None,
    )
    status = EvaluateInversePackingGeometry(table.evaluate_inverse_packing_geometry)(
        table.model_context,
        fingerprint.encode(),
        topology_fingerprint.encode(),
        298.15,
        amounts_buffer,
        len(amounts),
        log_packing_fraction,
        request_buffer,
        len(requests),
        ctypes.byref(result),
    )
    assert status == result.status
    return result, tuple(gradient), tuple(hessian)


def _fixed_topology_requests(
    slots: tuple[FixedTopologyParameterSlot, ...],
    values: tuple[float, ...],
) -> tuple[FixedTopologyParameterRequest, ...]:
    return tuple(
        FixedTopologyParameterRequest(
            struct_size=ctypes.sizeof(FixedTopologyParameterRequest),
            parameter_family=slot.parameter_family,
            component_index=slot.component_index,
            site_index_a=slot.site_index_a,
            site_index_b=slot.site_index_b,
            unit=slot.unit,
            value=value,
        )
        for slot, value in zip(slots, values, strict=True)
    )


def _evaluate_fixed_topology(
    table: FixedTopologyNativeSdkTable,
    fingerprint: str,
    topology_fingerprint: str,
    temperature_k: float,
    amount_mol: float,
    volume_m3: float,
    requests: tuple[FixedTopologyParameterRequest, ...],
    *,
    expected_artifact_fingerprint: str | None = None,
    output_sentinel: float = 0.0,
) -> tuple[FixedTopologyPhaseResult, tuple[float, ...], tuple[float, ...]]:
    coordinate_count = 2 + len(requests)
    gradient = (ctypes.c_double * coordinate_count)(
        *(output_sentinel for _ in range(coordinate_count))
    )
    hessian = (ctypes.c_double * (coordinate_count * coordinate_count))(
        *(output_sentinel for _ in range(coordinate_count * coordinate_count))
    )
    request_array = (FixedTopologyParameterRequest * len(requests))(*requests)
    result = FixedTopologyPhaseResult(
        struct_size=ctypes.sizeof(FixedTopologyPhaseResult),
        coordinate_count=coordinate_count,
        active_parameter_count=len(requests),
        gradient_capacity=coordinate_count,
        hessian_capacity=coordinate_count * coordinate_count,
        gradient=gradient,
        hessian=hessian,
    )
    status = EvaluateFixedTopologyAssociation(table.evaluate_fixed_topology_association)(
        table.model_context,
        fingerprint.encode(),
        topology_fingerprint.encode(),
        (
            expected_artifact_fingerprint
            or table.fixed_topology_association.contents.artifact_fingerprint.decode()
        ).encode(),
        temperature_k,
        amount_mol,
        volume_m3,
        request_array,
        len(requests),
        ctypes.byref(result),
    )
    assert status == result.status
    return result, tuple(gradient), tuple(hessian)


def _evaluate_aqueous_miac_kij(
    table: AqueousMiacKijNativeSdkTable,
    fingerprint: str,
    molality: float,
    k_ij: tuple[float, float, float],
) -> AqueousMiacKijResult:
    coordinates = (ctypes.c_double * 3)(*k_ij)
    result = AqueousMiacKijResult(struct_size=ctypes.sizeof(AqueousMiacKijResult))
    status = EvaluateAqueousMiacKij(table.evaluate_aqueous_miac_kij)(
        table.model_context,
        fingerprint.encode(),
        298.15,
        100000.0,
        molality,
        coordinates,
        3,
        ctypes.byref(result),
    )
    assert status == result.status
    return result


def _evaluate_aqueous_miac_kij_batch(
    table: AqueousMiacBatchNativeSdkTable,
    fingerprint: str,
    molalities: tuple[float, ...],
    k_ij: tuple[float, float, float],
) -> tuple[AqueousMiacKijResult, ...]:
    molality_values = (ctypes.c_double * len(molalities))(*molalities)
    coordinates = (ctypes.c_double * 3)(*k_ij)
    results = (AqueousMiacKijResult * len(molalities))(
        *(AqueousMiacKijResult(struct_size=ctypes.sizeof(AqueousMiacKijResult)) for _ in molalities)
    )
    status = EvaluateAqueousMiacKijBatch(table.evaluate_aqueous_miac_kij_batch)(
        table.model_context,
        fingerprint.encode(),
        298.15,
        100000.0,
        molality_values,
        len(molalities),
        coordinates,
        3,
        results,
        len(molalities),
    )
    assert all(result.status == status for result in results)
    return tuple(results)


def _evaluate_active_born(
    table: ActiveBornNativeSdkTable,
    born_diameter_angstrom: float,
    *,
    pressure_pa: float = 100000.0,
) -> IonSolvationBornResult:
    result = IonSolvationBornResult(struct_size=ctypes.sizeof(IonSolvationBornResult))
    status = EvaluateIonSolvationBorn(table.evaluate_ion_solvation_born)(
        table.model_context,
        298.15,
        pressure_pa,
        born_diameter_angstrom,
        ctypes.byref(result),
    )
    assert status == result.status
    return result


def _evaluate(
    table: NativeSdkTable,
    temperature_k: float,
    amount_mol: float,
    volume_m3: float,
) -> PhaseBlockResult:
    result = PhaseBlockResult(struct_size=ctypes.sizeof(PhaseBlockResult))
    status = EvaluatePurePhase(table.evaluate_pure_phase)(
        table.model_context,
        temperature_k,
        amount_mol,
        volume_m3,
        ctypes.byref(result),
    )
    assert status == result.status
    return result


def _mixture_result(coordinate_count: int) -> tuple[MixturePhaseBlockResult, object, object]:
    gradient = (ctypes.c_double * coordinate_count)()
    hessian = (ctypes.c_double * (coordinate_count * coordinate_count))()
    result = MixturePhaseBlockResult(
        struct_size=ctypes.sizeof(MixturePhaseBlockResult),
        coordinate_count=coordinate_count,
        gradient_capacity=coordinate_count,
        hessian_capacity=coordinate_count * coordinate_count,
        gradient=gradient,
        hessian=hessian,
    )
    return result, gradient, hessian


def _evaluate_mixture(
    table: NativeSdkTable,
    temperature_k: float,
    amounts_mol: tuple[float, ...],
    volume_m3: float,
) -> tuple[MixturePhaseBlockResult, tuple[float, ...], tuple[float, ...]]:
    assert table.table_size >= FIXED_MIXTURE_PREFIX_SIZE
    assert table.evaluate_mixture_phase
    coordinate_count = len(amounts_mol) + 1
    result, gradient, hessian = _mixture_result(coordinate_count)
    amounts = (ctypes.c_double * len(amounts_mol))(*amounts_mol)
    status = EvaluateMixturePhase(table.evaluate_mixture_phase)(
        table.model_context,
        temperature_k,
        amounts,
        len(amounts_mol),
        volume_m3,
        ctypes.byref(result),
    )
    assert status == result.status
    return result, tuple(gradient), tuple(hessian)


def _evaluate_mixture_kij(
    table: NativeSdkTable,
    temperature_k: float,
    amounts_mol: tuple[float, ...],
    volume_m3: float,
    k_ij: float,
) -> tuple[MixturePhaseBlockResult, tuple[float, ...], tuple[float, ...]]:
    assert table.table_size >= ACTIVE_KIJ_PREFIX_SIZE
    assert table.evaluate_mixture_phase_kij
    coordinate_count = len(amounts_mol) + 2
    result, gradient, hessian = _mixture_result(coordinate_count)
    amounts = (ctypes.c_double * len(amounts_mol))(*amounts_mol)
    status = EvaluateMixturePhaseKij(table.evaluate_mixture_phase_kij)(
        table.model_context,
        temperature_k,
        amounts,
        len(amounts_mol),
        volume_m3,
        k_ij,
        ctypes.byref(result),
    )
    assert status == result.status
    return result, tuple(gradient), tuple(hessian)


def _evaluate_mixture_lij(
    table: LijNativeSdkTable,
    temperature_k: float,
    amounts_mol: tuple[float, ...],
    volume_m3: float,
    l_ij: float,
) -> tuple[MixturePhaseBlockResult, tuple[float, ...], tuple[float, ...]]:
    assert table.evaluate_mixture_phase_lij
    coordinate_count = len(amounts_mol) + 2
    result, gradient, hessian = _mixture_result(coordinate_count)
    amounts = (ctypes.c_double * len(amounts_mol))(*amounts_mol)
    status = EvaluateMixturePhaseLij(table.evaluate_mixture_phase_lij)(
        table.model_context,
        temperature_k,
        amounts,
        len(amounts_mol),
        volume_m3,
        l_ij,
        ctypes.byref(result),
    )
    assert status == result.status
    return result, tuple(gradient), tuple(hessian)


def _evaluate_electrolyte(
    table: ElectrolyteNativeSdkTable,
    temperature_k: float,
    amounts_mol: tuple[float, ...],
    volume_m3: float,
) -> tuple[MixturePhaseBlockResult, tuple[float, ...], tuple[float, ...]]:
    coordinate_count = len(amounts_mol) + 1
    result, gradient, hessian = _mixture_result(coordinate_count)
    amounts = (ctypes.c_double * len(amounts_mol))(*amounts_mol)
    status = EvaluateMixturePhase(table.evaluate_electrolyte_phase)(
        table.model_context,
        temperature_k,
        amounts,
        len(amounts_mol),
        volume_m3,
        ctypes.byref(result),
    )
    assert status == result.status
    return result, tuple(gradient), tuple(hessian)


__all__ = tuple(name for name in globals() if not name.startswith("__"))
