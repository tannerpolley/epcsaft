from __future__ import annotations

from copy import deepcopy

from epcsaft.eos import (
    _native_sdk_configuration,
    _python_native_artifact_fingerprint,
)
from tests._raw_cases import AQUEOUS_ELECTROLYTE

from ._support import *


def test_native_sdk_addition_preserves_the_original_table_prefix() -> None:
    capsule = native_sdk(_model("neutral", "methane"))
    get_pointer = ctypes.pythonapi.PyCapsule_GetPointer
    get_pointer.argtypes = (ctypes.py_object, ctypes.c_char_p)
    get_pointer.restype = ctypes.c_void_p
    address = get_pointer(capsule, CAPSULE_NAME)
    prefix = ctypes.cast(address, ctypes.POINTER(NativeSdkPrefix)).contents
    assert prefix.abi_version == ABI_VERSION
    assert prefix.table_size >= ctypes.sizeof(NativeSdkPrefix)
    assert prefix.result_size == ctypes.sizeof(PhaseBlockResult)
    assert prefix.model_context
    assert prefix.evaluate_pure_phase
    for field_name, _ in AcceptedNativeSdkTable._fields_:
        assert (
            getattr(NativeSdkTable, field_name).offset
            == getattr(AcceptedNativeSdkTable, field_name).offset
        )
    for field_name, _ in FixedMixtureNativeSdkTable._fields_:
        assert (
            getattr(NativeSdkTable, field_name).offset
            == getattr(FixedMixtureNativeSdkTable, field_name).offset
        )
    assert FIXED_MIXTURE_PREFIX_SIZE == ctypes.sizeof(FixedMixtureNativeSdkTable)
    assert NativeSdkTable.evaluate_mixture_phase_kij.offset == FIXED_MIXTURE_PREFIX_SIZE
    for field_name, _ in ElectrolyteNativeSdkTable._fields_:
        assert (
            getattr(VolumeBoundsNativeSdkTable, field_name).offset
            == getattr(ElectrolyteNativeSdkTable, field_name).offset
        )
    assert VolumeBoundsNativeSdkTable.evaluate_molar_volume_bounds.offset == ELECTROLYTE_PREFIX_SIZE
    for field_name, _ in VolumeBoundsNativeSdkTable._fields_:
        assert (
            getattr(PackingFractionNativeSdkTable, field_name).offset
            == getattr(VolumeBoundsNativeSdkTable, field_name).offset
        )
    assert VOLUME_BOUNDS_PREFIX_SIZE == ctypes.sizeof(VolumeBoundsNativeSdkTable)
    assert (
        PackingFractionNativeSdkTable.evaluate_packing_fraction.offset == VOLUME_BOUNDS_PREFIX_SIZE
    )
    for field_name, _ in PackingFractionNativeSdkTable._fields_:
        assert (
            getattr(ApplicabilityNativeSdkTable, field_name).offset
            == getattr(PackingFractionNativeSdkTable, field_name).offset
        )
    assert PACKING_FRACTION_PREFIX_SIZE == ctypes.sizeof(PackingFractionNativeSdkTable) == 128
    assert ApplicabilityNativeSdkTable.source_temperature_min_k.offset == 128
    assert ApplicabilityNativeSdkTable.source_temperature_max_k.offset == 136
    assert ApplicabilityNativeSdkTable.total_ion_mole_fraction_max.offset == 144
    assert APPLICABILITY_PREFIX_SIZE == ctypes.sizeof(ApplicabilityNativeSdkTable) == 152
    for field_name, _ in AqueousMiacBoundedBatchNativeSdkTable._fields_:
        assert (
            getattr(CapabilityNativeSdkTable, field_name).offset
            == getattr(AqueousMiacBoundedBatchNativeSdkTable, field_name).offset
        )
    assert AQUEOUS_MIAC_BOUNDED_BATCH_PREFIX_SIZE == 264
    assert CapabilityNativeSdkTable.capability_count.offset == 264
    assert CapabilityNativeSdkTable.capabilities.offset == 272
    assert ctypes.sizeof(CapabilityNativeSdkTable) == 280
    for field_name, _ in CapabilityNativeSdkTable._fields_:
        assert (
            getattr(LijNativeSdkTable, field_name).offset
            == getattr(CapabilityNativeSdkTable, field_name).offset
        )
    assert LijNativeSdkTable.evaluate_mixture_phase_lij.offset == 280
    assert ctypes.sizeof(LijNativeSdkTable) == 288
    for field_name, _ in LijNativeSdkTable._fields_:
        assert (
            getattr(PureParameterNativeSdkTable, field_name).offset
            == getattr(LijNativeSdkTable, field_name).offset
        )
    assert PureParameterNativeSdkTable.evaluate_pure_phase_parameter.offset == 288
    assert ctypes.sizeof(PureParameterNativeSdkTable) == 296
    for field_name, _ in PureParameterNativeSdkTable._fields_:
        assert (
            getattr(IonFractionSuppressionNativeSdkTable, field_name).offset
            == getattr(PureParameterNativeSdkTable, field_name).offset
        )
    assert IonFractionSuppressionNativeSdkTable.ion_fraction_suppression_result_size.offset == 296
    assert IonFractionSuppressionNativeSdkTable.evaluate_ion_fraction_suppression.offset == 304
    assert (
        ION_FRACTION_SUPPRESSION_PREFIX_SIZE
        == ctypes.sizeof(IonFractionSuppressionNativeSdkTable)
        == 312
    )


def test_native_sdk_header_preserves_prefix_and_appends_public_operations() -> None:
    compiler = shutil.which("cc")
    assert compiler is not None
    include_root = Path(__file__).parents[2] / "src" / "epcsaft" / "include"
    source = r"""
#include <stddef.h>
#include <epcsaft/native_sdk_v1.h>

_Static_assert(
    sizeof(epcsaft_fixed_topology_parameter_request_v1) == 40,
    "fixed-topology request ABI"
);
_Static_assert(
    offsetof(epcsaft_native_sdk_v1, fixed_topology_association_result_size) == 448,
    "fixed-topology result-size offset"
);
_Static_assert(
    offsetof(epcsaft_native_sdk_v1, fixed_topology_association) == 456,
    "fixed-topology descriptor offset"
);
_Static_assert(
    offsetof(epcsaft_native_sdk_v1, evaluate_fixed_topology_association) == 464,
    "fixed-topology callback offset"
);
_Static_assert(
    offsetof(epcsaft_native_sdk_v1, source_reference_transfer_result_size) == 472,
    "source-reference result-size offset"
);
_Static_assert(
    offsetof(epcsaft_native_sdk_v1, evaluate_source_reference_transfer) == 480,
    "source-reference callback offset"
);
_Static_assert(
    sizeof(epcsaft_fixed_topology_mixture_phase_result_v1) == 664,
    "fixed-topology mixture result ABI"
);
_Static_assert(
    offsetof(
        epcsaft_native_sdk_v1,
        fixed_topology_associating_mixture_result_size
    ) == 488,
    "fixed-topology mixture result-size offset"
);
_Static_assert(
    offsetof(
        epcsaft_native_sdk_v1,
        evaluate_fixed_topology_associating_mixture
    ) == 496,
    "fixed-topology mixture callback offset"
);
_Static_assert(sizeof(epcsaft_native_sdk_v1) == 504, "generic SDK table size");

int main(void) {
    epcsaft_evaluate_fixed_topology_association_v1 callback = 0;
    epcsaft_evaluate_fixed_topology_associating_mixture_v1 mixture_callback = 0;
    epcsaft_evaluate_source_reference_transfer_v1 source_callback = 0;
    return callback != 0 || source_callback != 0 || mixture_callback != 0;
}
"""
    completed = subprocess.run(
        (compiler, "-std=c11", f"-I{include_root}", "-x", "c", "-fsyntax-only", "-"),
        input=source,
        text=True,
        capture_output=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr


@pytest.mark.parametrize(
    ("temperature_min_k", "temperature_max_k", "ion_mole_fraction_max"),
    (
        (math.nan, math.nan, 0.38),
        (math.nan, 298.15, 0.38),
        (298.15, math.nan, 0.38),
        (299.0, 298.15, 0.38),
        (298.15, 298.15, 0.0),
        (298.15, 298.15, 1.01),
        (298.15, 298.15, math.inf),
    ),
)
def test_native_sdk_applicability_metadata_rejects_missing_or_malformed_limits(
    temperature_min_k: float,
    temperature_max_k: float,
    ion_mole_fraction_max: float,
) -> None:
    model = _electrolyte_model(("water", "sodium-cation", "chloride-anion"))
    with pytest.raises(ValueError):
        _native_sdk_v1(
            model._native,
            model.parameter_fingerprint,
            temperature_min_k,
            temperature_max_k,
            model.component_ids,
            ion_mole_fraction_max,
        )


def test_native_sdk_generic_phase__does_not_admit_fixed_callbacks_without_domain() -> None:
    mapping = deepcopy(AQUEOUS_ELECTROLYTE)
    mapping["validity"] = {"kind": "unknown"}

    configuration = _native_sdk_configuration(
        Mixture.from_dictionary(mapping),
        _python_native_artifact_fingerprint(),
    )

    assert configuration["phase_callbacks_admitted"] is True
    for name in (
        "active_born_admitted",
        "aqueous_kij_admitted",
        "ion_fraction_suppression_admitted",
        "ion_solvation_kij_admitted",
        "aqueous_solvation_factor_admitted",
        "neutral_reference_admitted",
    ):
        assert configuration[name] is False


def test_native_sdk_packing_fraction_fails_closed_for_abi_model_and_input_errors() -> None:
    model = _electrolyte_model()
    capsule, table = _packing_fraction_table(model)
    assert capsule is not None
    callback = EvaluatePackingFraction(table.evaluate_packing_fraction)
    amounts = (ctypes.c_double * 3)(0.98, 0.01, 0.01)
    value = ctypes.c_double(11.0)
    gradient = (ctypes.c_double * 4)(*(12.0 for _ in range(4)))
    hessian = (ctypes.c_double * 16)(*(13.0 for _ in range(16)))

    def call(
        *,
        fingerprint: bytes | None = model.parameter_fingerprint.encode(),
        temperature_k: float = 298.15,
        amount_values: object = amounts,
        amount_count: int = 3,
        volume_m3: float = 1.0e-3,
        gradient_capacity: int = 4,
        hessian_capacity: int = 16,
    ) -> int:
        return callback(
            table.model_context,
            fingerprint,
            temperature_k,
            amount_values,
            amount_count,
            volume_m3,
            ctypes.byref(value),
            gradient,
            gradient_capacity,
            hessian,
            hessian_capacity,
        )

    def assert_untouched() -> None:
        assert value.value == 11.0
        assert tuple(gradient) == (12.0,) * 4
        assert tuple(hessian) == (13.0,) * 16

    assert call(gradient_capacity=3) == STATUS_ABI_MISMATCH
    assert_untouched()
    assert call(hessian_capacity=15) == STATUS_ABI_MISMATCH
    assert_untouched()
    assert call(fingerprint=b"sha256:wrong") == STATUS_UNSUPPORTED_MODEL
    assert_untouched()
    assert call(temperature_k=math.nan) == STATUS_INVALID_ARGUMENT
    assert_untouched()
    assert call(temperature_k=294.0) == STATUS_DOMAIN_ERROR
    assert_untouched()
    assert call(amount_count=4) == STATUS_INVALID_ARGUMENT
    assert_untouched()
    assert call(amount_values=(ctypes.c_double * 3)(0.98, 0.01, 0.0)) == (STATUS_INVALID_ARGUMENT)
    assert_untouched()
    assert call(amount_values=(ctypes.c_double * 3)(0.96, 0.03, 0.01)) == (STATUS_DOMAIN_ERROR)
    assert_untouched()
    assert call(volume_m3=math.nan) == STATUS_INVALID_ARGUMENT
    assert_untouched()
    assert call(volume_m3=0.0) == STATUS_INVALID_ARGUMENT
    assert_untouched()
