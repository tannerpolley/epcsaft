import ctypes
import gc
import math

import pytest

import epcsaft
from epcsaft import EPCSAFT, EosError, ParameterBundle, native_sdk, unit_registry


CAPSULE_NAME = b"epcsaft.native_sdk.v1"
ABI_VERSION = 1
STATUS_OK = 0


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


class NativeSdkTable(ctypes.Structure):
    _fields_ = (
        ("abi_version", ctypes.c_uint32),
        ("table_size", ctypes.c_size_t),
        ("result_size", ctypes.c_size_t),
        ("model_context", ctypes.c_void_p),
        ("evaluate_pure_phase", ctypes.c_void_p),
        ("parameterized_result_size", ctypes.c_size_t),
        ("evaluate_pure_phase_parameters", ctypes.c_void_p),
    )


class NativeSdkPrefix(ctypes.Structure):
    _fields_ = NativeSdkTable._fields_[:5]


EvaluatePurePhase = ctypes.CFUNCTYPE(
    ctypes.c_int,
    ctypes.c_void_p,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.POINTER(PhaseBlockResult),
)

EvaluatePurePhaseParameters = ctypes.CFUNCTYPE(
    ctypes.c_int,
    ctypes.c_void_p,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.POINTER(ParameterizedPhaseBlockResult),
)


def _model(catalog: str, component: str) -> EPCSAFT:
    return EPCSAFT(ParameterBundle.from_catalog(catalog, version=1).select((component,)))


def _table(model: EPCSAFT) -> tuple[object, NativeSdkTable]:
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
    return capsule, ctypes.cast(address, ctypes.POINTER(NativeSdkTable)).contents


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


def _evaluate_parameters(
    table: NativeSdkTable,
    temperature_k: float,
    amount_mol: float,
    volume_m3: float,
    segment_count: float,
    segment_diameter_angstrom: float,
    dispersion_energy_over_k: float,
) -> ParameterizedPhaseBlockResult:
    assert table.table_size >= ctypes.sizeof(NativeSdkTable)
    assert table.evaluate_pure_phase_parameters
    result = ParameterizedPhaseBlockResult(struct_size=ctypes.sizeof(ParameterizedPhaseBlockResult))
    status = EvaluatePurePhaseParameters(table.evaluate_pure_phase_parameters)(
        table.model_context,
        temperature_k,
        amount_mol,
        volume_m3,
        segment_count,
        segment_diameter_angstrom,
        dispersion_energy_over_k,
        ctypes.byref(result),
    )
    assert status == result.status
    return result


@pytest.mark.parametrize(
    "catalog,component,temperature_k,density_mol_m3",
    (
        ("gross-2001-methane-ethane", "methane", 190.0, 1000.0),
        ("gross-2001-methane-ethane", "ethane", 250.0, 1000.0),
        ("gross-2001-propane", "propane", 300.0, 1000.0),
    ),
)
def test_native_sdk_phase_block_matches_public_thermodynamics(
    catalog: str,
    component: str,
    temperature_k: float,
    density_mol_m3: float,
) -> None:
    model = _model(catalog, component)
    capsule, table = _table(model)
    assert capsule is not None
    assert table.abi_version == ABI_VERSION
    assert table.table_size == ctypes.sizeof(NativeSdkTable)
    assert table.result_size == ctypes.sizeof(PhaseBlockResult)
    assert table.parameterized_result_size == ctypes.sizeof(ParameterizedPhaseBlockResult)
    assert table.model_context
    assert table.evaluate_pure_phase
    assert table.evaluate_pure_phase_parameters

    amount_mol = 1.0
    volume_m3 = amount_mol / density_mol_m3
    block = _evaluate(table, temperature_k, amount_mol, volume_m3)
    assert block.status == STATUS_OK
    assert block.error == b""
    assert block.parameter_fingerprint.decode() == model.parameter_fingerprint

    state = model.evaluate(
        temperature=temperature_k * unit_registry.kelvin,
        molar_density=density_mol_m3 * unit_registry.mole / unit_registry.meter**3,
        mole_fractions=(1.0,),
    )
    expected_helmholtz = amount_mol * (math.log(density_mol_m3) - 1.0 + state.residual_helmholtz)
    expected_chemical_potential = (
        math.log(density_mol_m3) + state.residual_chemical_potential_over_rt[0]
    )
    assert block.helmholtz_over_rt_reference_amount == pytest.approx(
        expected_helmholtz, rel=2.0e-12, abs=2.0e-12
    )
    assert block.pressure_pa == pytest.approx(
        state.pressure.to("pascal").magnitude, rel=2.0e-11, abs=1.0e-5
    )
    assert block.chemical_potential_over_rt == pytest.approx(
        expected_chemical_potential, rel=2.0e-11, abs=2.0e-11
    )


def test_native_sdk_derivative_tensors_match_directional_oracles() -> None:
    model = _model("gross-2001-methane-ethane", "methane")
    capsule, table = _table(model)
    del model
    gc.collect()
    assert capsule is not None
    temperature_k = 190.0
    point = (1.2, 1.2e-3)
    direction = (0.05, -3.0e-5)

    def shifted(scale: float) -> PhaseBlockResult:
        result = _evaluate(
            table,
            temperature_k,
            point[0] + scale * direction[0],
            point[1] + scale * direction[1],
        )
        assert result.status == STATUS_OK
        return result

    center = shifted(0.0)
    step = 2.0e-4
    lower = shifted(-step)
    upper = shifted(step)

    gradient_contraction = sum(center.gradient[index] * direction[index] for index in range(2))
    value_difference = (
        upper.helmholtz_over_rt_reference_amount - lower.helmholtz_over_rt_reference_amount
    ) / (2.0 * step)
    assert gradient_contraction == pytest.approx(value_difference, rel=2.0e-8)

    for row in range(2):
        hessian_vector = sum(
            center.hessian[row * 2 + column] * direction[column] for column in range(2)
        )
        gradient_difference = (upper.gradient[row] - lower.gradient[row]) / (2.0 * step)
        assert hessian_vector == pytest.approx(gradient_difference, rel=2.0e-7)

    def hessian_contraction(result: PhaseBlockResult) -> float:
        return sum(
            direction[row] * result.hessian[row * 2 + column] * direction[column]
            for row in range(2)
            for column in range(2)
        )

    third_contraction = sum(
        direction[i] * direction[j] * direction[k] * center.third[(i * 2 + j) * 2 + k]
        for i in range(2)
        for j in range(2)
        for k in range(2)
    )
    hessian_difference = (hessian_contraction(upper) - hessian_contraction(lower)) / (2.0 * step)
    assert third_contraction == pytest.approx(hessian_difference, rel=2.0e-6)


def test_native_sdk_addition_preserves_the_original_table_prefix() -> None:
    capsule = native_sdk(_model("gross-2001-methane-ethane", "methane"))
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


def test_native_sdk_parameterized_phase_block_has_exact_parameter_derivatives() -> None:
    model = _model("gross-2001-methane-ethane", "methane")
    _, table = _table(model)
    temperature_k = 150.0
    point = (1.0, 2.0e-3, 1.0, 3.7039, 150.03)
    direction = (0.01, -2.0e-5, 0.02, -0.03, 0.5)

    def shifted(scale: float) -> ParameterizedPhaseBlockResult:
        values = tuple(value + scale * delta for value, delta in zip(point, direction))
        result = _evaluate_parameters(table, temperature_k, *values)
        assert result.status == STATUS_OK
        return result

    center = shifted(0.0)
    fixed = _evaluate(table, temperature_k, point[0], point[1])
    assert center.helmholtz_over_rt_reference_amount == pytest.approx(
        fixed.helmholtz_over_rt_reference_amount, rel=2.0e-13, abs=2.0e-13
    )
    assert center.pressure_pa == pytest.approx(fixed.pressure_pa, rel=2.0e-13)
    assert center.chemical_potential_over_rt == pytest.approx(
        fixed.chemical_potential_over_rt, rel=2.0e-13
    )
    assert center.parameter_fingerprint.decode() == model.parameter_fingerprint

    step = 2.0e-5
    lower = shifted(-step)
    upper = shifted(step)
    value_direction = sum(center.gradient[index] * direction[index] for index in range(5))
    value_difference = (
        upper.helmholtz_over_rt_reference_amount - lower.helmholtz_over_rt_reference_amount
    ) / (2.0 * step)
    assert value_direction == pytest.approx(value_difference, rel=3.0e-8)

    for row in range(5):
        hessian_direction = sum(
            center.hessian[row * 5 + column] * direction[column] for column in range(5)
        )
        gradient_difference = (upper.gradient[row] - lower.gradient[row]) / (2.0 * step)
        assert hessian_direction == pytest.approx(gradient_difference, rel=3.0e-7)

    gas_constant = 8.31446261815324
    assert center.pressure_pa == pytest.approx(
        -gas_constant * temperature_k * center.gradient[1], rel=2.0e-13
    )
    assert center.chemical_potential_over_rt == pytest.approx(center.gradient[0], rel=2.0e-13)


def test_native_sdk_rejects_unsupported_models_and_invalid_states() -> None:
    binary = EPCSAFT(
        ParameterBundle.from_catalog("gross-2001-methane-ethane", version=1).select(
            ("methane", "ethane")
        )
    )
    with pytest.raises(EosError, match="unsupported-native-sdk"):
        native_sdk(binary)

    associating = _model("figiel-2025-reference-electrolytes", "ethanol")
    with pytest.raises(EosError, match="unsupported-native-sdk"):
        native_sdk(associating)

    model = _model("gross-2001-methane-ethane", "methane")
    gross_bundle = ParameterBundle.from_catalog("gross-2001-methane-ethane", version=1)
    gross = gross_bundle.select(("methane",))
    synthetic = ParameterBundle.from_records(
        bundle_id="synthetic-gross-copy",
        bundle_version=1,
        purpose="package-test-fixture",
        sources=gross_bundle.sources,
        domains=gross.domains,
        components=gross.components,
        singles=gross.records,
    )
    synthetic_model = EPCSAFT(synthetic.select(("methane",)))
    assert synthetic_model.parameter_fingerprint != model.parameter_fingerprint
    with pytest.raises(EosError, match="unadmitted-parameter-fingerprint"):
        native_sdk(synthetic_model)

    class SpoofedEPCSAFT(EPCSAFT):
        @property
        def parameter_fingerprint(self) -> str:
            return model.parameter_fingerprint

    spoofed_model = SpoofedEPCSAFT(synthetic.select(("methane",)))
    with pytest.raises(EosError, match="expected-exact-EPCSAFT-model"):
        native_sdk(spoofed_model)

    _, table = _table(model)
    for temperature, amount, volume in (
        (96.0, 1.0, 1.0e-3),
        (190.0, 0.0, 1.0e-3),
        (190.0, 1.0, 0.0),
        (math.nan, 1.0, 1.0e-3),
    ):
        result = _evaluate(table, temperature, amount, volume)
        assert result.status != STATUS_OK
        assert result.error

    for segment_count, segment_diameter, dispersion_energy in (
        (0.0, 3.7039, 150.03),
        (1.0, 0.0, 150.03),
        (1.0, 3.7039, 0.0),
        (1.0, 3.7039, math.nan),
    ):
        result = _evaluate_parameters(
            table,
            150.0,
            1.0,
            2.0e-3,
            segment_count,
            segment_diameter,
            dispersion_energy,
        )
        assert result.status != STATUS_OK
        assert result.error


def test_native_sdk_header_is_installed_package_data() -> None:
    header = (
        __import__("importlib.resources").resources.files(epcsaft)
        / "include"
        / "epcsaft"
        / "native_sdk_v1.h"
    )
    assert header.is_file()
