from __future__ import annotations

import ctypes

import pytest
from parameter_cases import aqueous_parameters

from epcsaft import Mixture, Parameters, native_sdk
from epcsaft import regression
from epcsaft.regression import _native
from epcsaft.regression.workflow import _aqueous_kij_native_payload


def _models() -> tuple[Mixture, ...]:
    return tuple(
        Mixture(Parameters.from_dictionary(aqueous_parameters(("water", cation, anion))))
        for _, cation, anion, _ in regression.FIGIEL_AQUEOUS_KIJ_V1.salt_contracts
    )


class _NativeSdkV1(ctypes.Structure):
    _fields_ = [
        ("abi_version", ctypes.c_uint32),
        ("table_size", ctypes.c_size_t),
        ("result_size", ctypes.c_size_t),
        ("model_context", ctypes.c_void_p),
        ("evaluate_pure_phase", ctypes.c_void_p),
        ("parameterized_result_size", ctypes.c_size_t),
        ("evaluate_pure_phase_parameters", ctypes.c_void_p),
        ("component_count", ctypes.c_size_t),
        ("mixture_result_size", ctypes.c_size_t),
        ("evaluate_mixture_phase", ctypes.c_void_p),
        ("evaluate_mixture_phase_kij", ctypes.c_void_p),
        ("component_ids", ctypes.c_void_p),
        ("component_charges", ctypes.c_void_p),
        ("evaluate_electrolyte_phase", ctypes.c_void_p),
        ("evaluate_molar_volume_bounds", ctypes.c_void_p),
        ("evaluate_packing_fraction", ctypes.c_void_p),
        ("source_temperature_min_k", ctypes.c_double),
        ("source_temperature_max_k", ctypes.c_double),
        ("total_ion_mole_fraction_max", ctypes.c_double),
        ("ion_solvation_born_result_size", ctypes.c_size_t),
        ("evaluate_ion_solvation_born", ctypes.c_void_p),
        ("aqueous_miac_kij_result_size", ctypes.c_size_t),
        ("evaluate_aqueous_miac_kij", ctypes.c_void_p),
        ("neutral_reference_basis_row_count", ctypes.c_size_t),
        ("neutral_reference_result_size", ctypes.c_size_t),
        ("evaluate_neutral_reference", ctypes.c_void_p),
        ("aqueous_miac_solvation_factor_result_size", ctypes.c_size_t),
        ("evaluate_aqueous_miac_solvation_factor", ctypes.c_void_p),
        ("evaluate_aqueous_miac_kij_batch", ctypes.c_void_p),
        ("evaluate_aqueous_miac_solvation_factor_batch", ctypes.c_void_p),
        ("evaluation_budget_size", ctypes.c_size_t),
        ("evaluate_aqueous_miac_kij_batch_bounded", ctypes.c_void_p),
        (
            "evaluate_aqueous_miac_solvation_factor_batch_bounded",
            ctypes.c_void_p,
        ),
    ]


def _copied_provider_capsules(
    capsules: tuple[object, ...],
    *,
    clear_scalar: bool = False,
    truncate_before_bounded_batch: bool = False,
) -> tuple[tuple[object, ...], tuple[_NativeSdkV1, ...]]:
    get_pointer = ctypes.pythonapi.PyCapsule_GetPointer
    get_pointer.argtypes = (ctypes.py_object, ctypes.c_char_p)
    get_pointer.restype = ctypes.c_void_p
    new_capsule = ctypes.pythonapi.PyCapsule_New
    new_capsule.argtypes = (ctypes.c_void_p, ctypes.c_char_p, ctypes.c_void_p)
    new_capsule.restype = ctypes.py_object
    name = b"epcsaft.native_sdk.v1"
    tables: list[_NativeSdkV1] = []
    copied_capsules: list[object] = []
    for capsule in capsules:
        pointer = get_pointer(capsule, name)
        table = _NativeSdkV1()
        ctypes.memmove(ctypes.addressof(table), pointer, ctypes.sizeof(table))
        assert table.table_size >= ctypes.sizeof(table)
        table.table_size = ctypes.sizeof(table)
        if clear_scalar:
            table.evaluate_aqueous_miac_kij = None
        if truncate_before_bounded_batch:
            table.table_size = _NativeSdkV1.evaluate_aqueous_miac_kij_batch_bounded.offset
        tables.append(table)
        copied_capsules.append(new_capsule(ctypes.addressof(table), name, None))
    return tuple(copied_capsules), tuple(tables)


@pytest.mark.campaign
def test_aqueous_kij_residual_jacobian_maps_all_exact_provider_columns() -> None:
    specification = regression.FIGIEL_AQUEOUS_KIJ_V1
    models = _models()
    capsules = tuple(native_sdk(model) for model in models)
    payload = _aqueous_kij_native_payload(
        specification, tuple(model.parameter_fingerprint for model in models)
    )

    _, jacobian_native, rows_native = _native.evaluate_figiel_kij(
        capsules, payload, specification.published_parameters
    )
    jacobian = tuple(float(value) for value in jacobian_native)
    column_map = {salt: columns for salt, _, _, columns in specification.salt_contracts}
    for row_index, row_native in enumerate(rows_native):
        modeled = float(row_native[5])
        observed = float(row_native[3])
        local = tuple(float(value) for value in row_native[7])
        expected = [0.0] * 11
        for column, derivative in zip(column_map[str(row_native[1])], local, strict=True):
            expected[column] = -(modeled / observed) * derivative
        actual = jacobian[row_index * 11 : (row_index + 1) * 11]
        assert actual == pytest.approx(expected, rel=1.0e-14, abs=1.0e-14)


@pytest.mark.campaign
def test_aqueous_kij_exact_jacobian_matches_callback_value_direction() -> None:
    specification = regression.FIGIEL_AQUEOUS_KIJ_V1
    models = _models()
    capsules = tuple(native_sdk(model) for model in models)
    payload = _aqueous_kij_native_payload(
        specification, tuple(model.parameter_fingerprint for model in models)
    )
    parameters = specification.published_parameters
    direction = tuple((-1.0 if index % 2 else 1.0) * (index + 1.0) / 11.0 for index in range(11))
    _, jacobian_native, _ = _native.evaluate_figiel_kij(capsules, payload, parameters)
    jacobian = tuple(float(value) for value in jacobian_native)
    differences: list[tuple[float, ...]] = []
    for step in (1.0e-4, 5.0e-5):
        plus = tuple(
            value + step * component for value, component in zip(parameters, direction, strict=True)
        )
        minus = tuple(
            value - step * component for value, component in zip(parameters, direction, strict=True)
        )
        residuals_plus, _, _ = _native.evaluate_figiel_kij(capsules, payload, plus)
        residuals_minus, _, _ = _native.evaluate_figiel_kij(capsules, payload, minus)
        differences.append(
            tuple(
                (float(plus_value) - float(minus_value)) / (2.0 * step)
                for plus_value, minus_value in zip(residuals_plus, residuals_minus, strict=True)
            )
        )
    for row in range(164):
        exact = sum(jacobian[row * 11 + column] * direction[column] for column in range(11))
        coarse = differences[0][row]
        fine = differences[1][row]
        tolerance = max(
            1.0e-8,
            20.0 * abs(coarse - fine),
            2.0e-8 * abs(exact),
        )
        assert abs(exact - fine) <= tolerance


@pytest.mark.campaign
def test_public_aqueous_kij_fit_replays_retained_negative_result() -> None:
    result = regression.fit_figiel_aqueous_kij(models=_models())

    assert result.solver_converged
    assert result.numerically_converged
    assert result.physically_valid
    assert result.workflow_valid
    assert not result.scientifically_valid
    assert result.recovery_status == (
        "SOURCE_DESCRIBED_STAGED_RECOVERY_DID_NOT_REPRODUCE_PRINTED_TUPLE"
    )
    assert result.input_row_ids == result.evaluated_row_ids
    assert result.failed_row_ids == ()
    assert tuple(start.rank for start in result.starts) == (11, 11)
    assert result.start_parameter_max_abs_delta <= 1.0e-5
    assert result.published_parameter_max_abs_delta > 0.05
