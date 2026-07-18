import ctypes
import gc
import math

import pytest

import epcsaft
from epcsaft import EPCSAFT, EosError, ParameterBundle, native_sdk, unit_registry


CAPSULE_NAME = b"epcsaft.native_sdk.v1"
ABI_VERSION = 1
STATUS_OK = 0
STATUS_INVALID_ARGUMENT = 1
STATUS_UNSUPPORTED_MODEL = 2
STATUS_DOMAIN_ERROR = 3
STATUS_ABI_MISMATCH = 5


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


class NativeSdkPrefix(ctypes.Structure):
    _fields_ = NativeSdkTable._fields_[:5]


FIXED_MIXTURE_PREFIX_SIZE = (
    FixedMixtureNativeSdkTable.evaluate_mixture_phase.offset + ctypes.sizeof(ctypes.c_void_p)
)
ACTIVE_KIJ_PREFIX_SIZE = NativeSdkTable.evaluate_mixture_phase_kij.offset + ctypes.sizeof(
    ctypes.c_void_p
)


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


def _model(catalog: str, component: str) -> EPCSAFT:
    return EPCSAFT(ParameterBundle.from_catalog(catalog, version=1).select((component,)))


def _mixture_model(component_ids: tuple[str, str] = ("methane", "ethane")) -> EPCSAFT:
    parameters = ParameterBundle.from_catalog("gross-2001-methane-ethane", version=1)
    return EPCSAFT(parameters.select(component_ids))


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
    fixed_table = ctypes.cast(address, ctypes.POINTER(FixedMixtureNativeSdkTable)).contents
    assert fixed_table.table_size >= ACTIVE_KIJ_PREFIX_SIZE
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


def test_native_sdk_mixture_phase_block_matches_public_thermodynamics_and_derivatives() -> None:
    model = _mixture_model()
    capsule, table = _table(model)
    assert capsule is not None
    assert table.component_count == 2
    assert table.mixture_result_size == ctypes.sizeof(MixturePhaseBlockResult)
    assert table.evaluate_mixture_phase

    temperature_k = 220.0
    point = (0.4, 0.6, 1.25e-4)
    direction = (0.02, -0.01, -2.0e-6)

    def shifted(
        scale: float,
    ) -> tuple[MixturePhaseBlockResult, tuple[float, ...], tuple[float, ...]]:
        return _evaluate_mixture(
            table,
            temperature_k,
            (point[0] + scale * direction[0], point[1] + scale * direction[1]),
            point[2] + scale * direction[2],
        )

    center, gradient, hessian = shifted(0.0)
    assert center.status == STATUS_OK
    assert center.error == b""
    assert center.coordinate_count == 3
    assert center.gradient_capacity == 3
    assert center.hessian_capacity == 9
    assert center.parameter_fingerprint.decode() == model.parameter_fingerprint

    total = point[0] + point[1]
    density = total / point[2]
    composition = (point[0] / total, point[1] / total)
    state = model.evaluate(
        temperature=temperature_k * unit_registry.kelvin,
        molar_density=density * unit_registry.mole / unit_registry.meter**3,
        mole_fractions=composition,
    )
    expected = sum(amount * (math.log(amount / point[2]) - 1.0) for amount in point[:2])
    expected += total * state.residual_helmholtz
    assert center.helmholtz_over_rt_reference_amount == pytest.approx(
        expected, rel=2.0e-12, abs=2.0e-12
    )
    assert center.pressure_pa == pytest.approx(
        state.pressure.to("pascal").magnitude, rel=2.0e-11, abs=1.0e-5
    )

    step = 2.0e-4
    lower, lower_gradient, _ = shifted(-step)
    upper, upper_gradient, _ = shifted(step)
    gradient_contraction = sum(value * delta for value, delta in zip(gradient, direction))
    value_difference = (
        upper.helmholtz_over_rt_reference_amount - lower.helmholtz_over_rt_reference_amount
    ) / (2.0 * step)
    assert gradient_contraction == pytest.approx(value_difference, rel=3.0e-8)
    for row in range(3):
        hessian_direction = sum(
            hessian[row * 3 + column] * direction[column] for column in range(3)
        )
        gradient_difference = (upper_gradient[row] - lower_gradient[row]) / (2.0 * step)
        assert hessian_direction == pytest.approx(gradient_difference, rel=3.0e-7)


def test_native_sdk_mixture_phase_block_is_component_permutation_invariant() -> None:
    forward_capsule, forward_table = _table(_mixture_model())
    reverse_capsule, reverse_table = _table(_mixture_model(("ethane", "methane")))
    assert forward_capsule is not None
    assert reverse_capsule is not None
    forward, forward_gradient, forward_hessian = _evaluate_mixture(
        forward_table, 220.0, (0.4, 0.6), 1.25e-4
    )
    reverse, reverse_gradient, reverse_hessian = _evaluate_mixture(
        reverse_table, 220.0, (0.6, 0.4), 1.25e-4
    )
    permutation = (1, 0, 2)
    assert reverse.helmholtz_over_rt_reference_amount == pytest.approx(
        forward.helmholtz_over_rt_reference_amount, rel=2.0e-13, abs=2.0e-13
    )
    assert reverse.pressure_pa == pytest.approx(forward.pressure_pa, rel=2.0e-13)
    for row in range(3):
        assert reverse_gradient[row] == pytest.approx(
            forward_gradient[permutation[row]], rel=2.0e-13, abs=2.0e-13
        )
        for column in range(3):
            assert reverse_hessian[row * 3 + column] == pytest.approx(
                forward_hessian[permutation[row] * 3 + permutation[column]],
                rel=2.0e-13,
                abs=2.0e-13,
            )


@pytest.mark.parametrize(
    ("temperature_k", "point", "direction"),
    (
        (220.0, (0.4, 0.6, 1.25e-4, 0.0), (0.02, -0.01, -2.0e-6, 0.01)),
        (280.0, (0.7, 0.3, 2.5e-4, -0.05), (-0.015, 0.02, 3.0e-6, -0.008)),
    ),
)
def test_native_sdk_active_kij_matches_fixed_and_directional_oracles(
    temperature_k: float,
    point: tuple[float, float, float, float],
    direction: tuple[float, float, float, float],
) -> None:
    model = _mixture_model()
    capsule, table = _table(model)
    assert capsule is not None
    assert table.evaluate_mixture_phase_kij

    def shifted(
        scale: float,
    ) -> tuple[MixturePhaseBlockResult, tuple[float, ...], tuple[float, ...]]:
        values = tuple(value + scale * delta for value, delta in zip(point, direction))
        return _evaluate_mixture_kij(
            table,
            temperature_k,
            (values[0], values[1]),
            values[2],
            values[3],
        )

    center, gradient, hessian = shifted(0.0)
    assert center.status == STATUS_OK
    assert center.error == b""
    assert center.coordinate_count == 4
    assert center.gradient_capacity == 4
    assert center.hessian_capacity == 16
    assert center.parameter_fingerprint.decode() == model.parameter_fingerprint
    assert center.pressure_pa == pytest.approx(
        -8.31446261815324 * temperature_k * gradient[2], rel=2.0e-13
    )

    if point[3] == 0.0:
        fixed, fixed_gradient, fixed_hessian = _evaluate_mixture(
            table, temperature_k, point[:2], point[2]
        )
        assert center.helmholtz_over_rt_reference_amount == pytest.approx(
            fixed.helmholtz_over_rt_reference_amount, rel=2.0e-13, abs=2.0e-13
        )
        assert center.pressure_pa == pytest.approx(fixed.pressure_pa, rel=2.0e-13)
        for row in range(3):
            assert gradient[row] == pytest.approx(fixed_gradient[row], rel=2.0e-13, abs=2.0e-13)
            for column in range(3):
                assert hessian[row * 4 + column] == pytest.approx(
                    fixed_hessian[row * 3 + column], rel=2.0e-13, abs=2.0e-13
                )

    step = 2.0e-4
    lower, lower_gradient, _ = shifted(-step)
    upper, upper_gradient, _ = shifted(step)
    gradient_direction = sum(value * delta for value, delta in zip(gradient, direction))
    value_difference = (
        upper.helmholtz_over_rt_reference_amount - lower.helmholtz_over_rt_reference_amount
    ) / (2.0 * step)
    assert gradient_direction == pytest.approx(value_difference, rel=3.0e-8)
    for row in range(4):
        hessian_direction = sum(
            hessian[row * 4 + column] * direction[column] for column in range(4)
        )
        gradient_difference = (upper_gradient[row] - lower_gradient[row]) / (2.0 * step)
        assert hessian_direction == pytest.approx(gradient_difference, rel=3.0e-7)


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
    associating = _model("figiel-2025-reference-electrolytes", "ethanol")
    with pytest.raises(EosError, match="unsupported-native-sdk"):
        native_sdk(associating)
    ionic_bundle = ParameterBundle.from_catalog("figiel-2025-reference-electrolytes", version=1)
    ionic = EPCSAFT(ionic_bundle.select(("water", "sodium-cation", "chloride-anion")))
    with pytest.raises(EosError, match="unsupported-native-sdk"):
        native_sdk(ionic)

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


def test_native_sdk_mixture_callback_rejects_malformed_buffers_and_scope() -> None:
    capsule, table = _table(_mixture_model())
    assert capsule is not None
    callback = EvaluateMixturePhase(table.evaluate_mixture_phase)
    amounts = (ctypes.c_double * 2)(0.4, 0.6)

    def call(
        result: MixturePhaseBlockResult,
        *,
        amount_pointer: object = amounts,
        amount_count: int = 2,
        temperature_k: float = 220.0,
    ) -> int:
        return callback(
            table.model_context,
            temperature_k,
            amount_pointer,
            amount_count,
            1.25e-4,
            ctypes.byref(result),
        )

    malformed: list[tuple[MixturePhaseBlockResult, dict[str, object], int]] = []
    result, _, _ = _mixture_result(3)
    malformed.append((result, {"amount_pointer": None}, STATUS_INVALID_ARGUMENT))
    result, _, _ = _mixture_result(3)
    malformed.append((result, {"amount_count": 1}, STATUS_INVALID_ARGUMENT))
    result, _, _ = _mixture_result(3)
    result.gradient = ctypes.POINTER(ctypes.c_double)()
    malformed.append((result, {}, STATUS_INVALID_ARGUMENT))
    result, _, _ = _mixture_result(3)
    result.hessian = ctypes.POINTER(ctypes.c_double)()
    malformed.append((result, {}, STATUS_INVALID_ARGUMENT))
    result, _, _ = _mixture_result(3)
    result.gradient_capacity = 2
    malformed.append((result, {}, STATUS_INVALID_ARGUMENT))
    result, _, _ = _mixture_result(3)
    result.hessian_capacity = 8
    malformed.append((result, {}, STATUS_INVALID_ARGUMENT))
    result, _, _ = _mixture_result(3)
    malformed.append((result, {"temperature_k": 96.0}, STATUS_DOMAIN_ERROR))

    for result, arguments, expected_status in malformed:
        gradient_pointer = ctypes.cast(result.gradient, ctypes.c_void_p).value
        hessian_pointer = ctypes.cast(result.hessian, ctypes.c_void_p).value
        status = call(result, **arguments)
        assert status == expected_status
        assert result.status == status
        assert result.error
        assert ctypes.cast(result.gradient, ctypes.c_void_p).value == gradient_pointer
        assert ctypes.cast(result.hessian, ctypes.c_void_p).value == hessian_pointer

    result, _, _ = _mixture_result(3)
    result.struct_size -= 1
    assert call(result) == STATUS_ABI_MISMATCH

    result, _, _ = _mixture_result(3)
    nonpositive_amounts = (ctypes.c_double * 2)(0.0, 1.0)
    assert call(result, amount_pointer=nonpositive_amounts) == STATUS_INVALID_ARGUMENT
    assert result.status == STATUS_INVALID_ARGUMENT
    assert result.error

    pure_capsule, pure_table = _table(_model("gross-2001-methane-ethane", "methane"))
    assert pure_capsule is not None
    pure_mixture_result, _, _ = _mixture_result(2)
    pure_amount = (ctypes.c_double * 1)(1.0)
    status = EvaluateMixturePhase(pure_table.evaluate_mixture_phase)(
        pure_table.model_context,
        190.0,
        pure_amount,
        1,
        1.0e-3,
        ctypes.byref(pure_mixture_result),
    )
    assert status == STATUS_UNSUPPORTED_MODEL
    assert pure_mixture_result.status == status
    assert pure_mixture_result.error

    mixture_pure_result = _evaluate(table, 220.0, 1.0, 1.25e-4)
    assert mixture_pure_result.status == STATUS_UNSUPPORTED_MODEL
    assert mixture_pure_result.error

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


def test_native_sdk_active_kij_rejects_malformed_buffers_models_and_nonfinite_inputs() -> None:
    capsule, table = _table(_mixture_model())
    assert capsule is not None
    callback = EvaluateMixturePhaseKij(table.evaluate_mixture_phase_kij)
    amounts = (ctypes.c_double * 2)(0.4, 0.6)

    def call(
        result: MixturePhaseBlockResult,
        *,
        amount_pointer: object = amounts,
        amount_count: int = 2,
        temperature_k: float = 220.0,
        volume_m3: float = 1.25e-4,
        k_ij: float = 0.0,
    ) -> int:
        return callback(
            table.model_context,
            temperature_k,
            amount_pointer,
            amount_count,
            volume_m3,
            k_ij,
            ctypes.byref(result),
        )

    malformed: list[tuple[MixturePhaseBlockResult, dict[str, object], int]] = []
    result, _, _ = _mixture_result(3)
    malformed.append((result, {}, STATUS_INVALID_ARGUMENT))
    result, _, _ = _mixture_result(4)
    result.gradient = ctypes.POINTER(ctypes.c_double)()
    malformed.append((result, {}, STATUS_INVALID_ARGUMENT))
    result, _, _ = _mixture_result(4)
    result.hessian_capacity = 15
    malformed.append((result, {}, STATUS_INVALID_ARGUMENT))
    result, _, _ = _mixture_result(4)
    malformed.append((result, {"amount_pointer": None}, STATUS_INVALID_ARGUMENT))
    result, _, _ = _mixture_result(4)
    malformed.append((result, {"amount_count": 1}, STATUS_INVALID_ARGUMENT))
    result, _, _ = _mixture_result(4)
    malformed.append((result, {"temperature_k": math.nan}, STATUS_INVALID_ARGUMENT))
    result, _, _ = _mixture_result(4)
    malformed.append((result, {"volume_m3": math.inf}, STATUS_INVALID_ARGUMENT))
    result, _, _ = _mixture_result(4)
    malformed.append((result, {"k_ij": math.nan}, STATUS_INVALID_ARGUMENT))
    result, _, _ = _mixture_result(4)
    malformed.append((result, {"k_ij": math.inf}, STATUS_INVALID_ARGUMENT))

    for result, arguments, expected_status in malformed:
        gradient_pointer = ctypes.cast(result.gradient, ctypes.c_void_p).value
        hessian_pointer = ctypes.cast(result.hessian, ctypes.c_void_p).value
        status = call(result, **arguments)
        assert status == expected_status
        assert result.status == status
        assert result.error
        assert ctypes.cast(result.gradient, ctypes.c_void_p).value == gradient_pointer
        assert ctypes.cast(result.hessian, ctypes.c_void_p).value == hessian_pointer

    result, _, _ = _mixture_result(4)
    result.struct_size -= 1
    assert call(result) == STATUS_ABI_MISMATCH

    finite_unbounded, _, _ = _mixture_result(4)
    assert call(finite_unbounded, k_ij=-0.2) == STATUS_OK
    assert finite_unbounded.status == STATUS_OK

    pure_capsule, pure_table = _table(_model("gross-2001-methane-ethane", "methane"))
    assert pure_capsule is not None
    pure_result, _, _ = _mixture_result(3)
    pure_amount = (ctypes.c_double * 1)(1.0)
    status = EvaluateMixturePhaseKij(pure_table.evaluate_mixture_phase_kij)(
        pure_table.model_context,
        190.0,
        pure_amount,
        1,
        1.0e-3,
        0.0,
        ctypes.byref(pure_result),
    )
    assert status == STATUS_UNSUPPORTED_MODEL
    assert pure_result.status == status
    assert pure_result.error


def test_native_sdk_header_is_installed_package_data() -> None:
    header = (
        __import__("importlib.resources").resources.files(epcsaft)
        / "include"
        / "epcsaft"
        / "native_sdk_v1.h"
    )
    assert header.is_file()
