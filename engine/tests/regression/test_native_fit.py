from __future__ import annotations

import ctypes
import gc
import sys
from functools import partial
from pathlib import Path

import epcsaft.regression._native as native
import pytest

from epcsaft import Mixture, Parameters, native_sdk
from epcsaft.regression import workflow
from epcsaft.regression.records import (
    ETHANE_SATURATION_FIT_V1,
    METHANE_SATURATION_FIT_V1,
    PROPANE_SATURATION_FIT_V1,
    load_pure_saturation_dataset,
)
from epcsaft.regression.workflow import _native_payload, fit_pure_saturation

SPECIFICATIONS = {
    "methane": METHANE_SATURATION_FIT_V1,
    "ethane": ETHANE_SATURATION_FIT_V1,
    "propane": PROPANE_SATURATION_FIT_V1,
}


class _ParameterizedResult(ctypes.Structure):
    _fields_ = (
        ("struct_size", ctypes.c_uint32),
        ("status", ctypes.c_int32),
        ("helmholtz", ctypes.c_double),
        ("gradient", ctypes.c_double * 5),
        ("hessian", ctypes.c_double * 25),
        ("pressure", ctypes.c_double),
        ("chemical_potential", ctypes.c_double),
        ("fingerprint", ctypes.c_char * 72),
        ("error", ctypes.c_char * 160),
    )


_ParameterizedCallback = ctypes.CFUNCTYPE(
    ctypes.c_int,
    ctypes.c_void_p,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.POINTER(_ParameterizedResult),
)


class _NativeSdkTable(ctypes.Structure):
    _fields_ = (
        ("abi_version", ctypes.c_uint32),
        ("table_size", ctypes.c_size_t),
        ("result_size", ctypes.c_size_t),
        ("model_context", ctypes.c_void_p),
        ("evaluate_pure_phase", ctypes.c_void_p),
        ("parameterized_result_size", ctypes.c_size_t),
        ("evaluate_pure_phase_parameters", ctypes.c_void_p),
    )


def _model(component_id: str) -> Mixture:
    packet = "gross-2001-propane" if component_id == "propane" else "gross-2001-methane-ethane"
    bundle = Path(__file__).parents[1] / "fixtures" / "data"
    return Mixture(
        Parameters.from_bundle(
            bundle / packet / "1" / "parameters",
            components=(component_id,),
        )
    )


def _capsule(component_id: str) -> object:
    return native_sdk(_model(component_id))


def _provider_table(capsule: object) -> _NativeSdkTable:
    capsule_pointer = ctypes.pythonapi.PyCapsule_GetPointer
    capsule_pointer.argtypes = (ctypes.py_object, ctypes.c_char_p)
    capsule_pointer.restype = ctypes.c_void_p
    pointer = capsule_pointer(capsule, b"epcsaft.native_sdk.v1")
    assert pointer
    return ctypes.cast(pointer, ctypes.POINTER(_NativeSdkTable)).contents


def _provider_phase(
    table: _NativeSdkTable,
    temperature_k: float,
    volume_m3: float,
    parameters: tuple[float, float, float],
) -> _ParameterizedResult:
    callback = _ParameterizedCallback(table.evaluate_pure_phase_parameters)
    result = _ParameterizedResult()
    result.struct_size = ctypes.sizeof(_ParameterizedResult)
    status = callback(
        table.model_context,
        temperature_k,
        1.0,
        volume_m3,
        *parameters,
        ctypes.byref(result),
    )
    assert status == result.status == 0
    return result


def _payload(component_id: str) -> tuple[object, ...]:
    dataset = load_pure_saturation_dataset(component_id)
    specification = SPECIFICATIONS[component_id]
    model = _model(component_id)
    return _native_payload(dataset, specification, model.parameter_fingerprint)


def test_methane_start_residuals_match_accepted_provider_anchor() -> None:
    residuals, jacobian, diagnostics, fingerprint = native.evaluate(
        _capsule("methane"), _payload("methane"), (0.0,) * 11
    )

    assert residuals[4:8] == pytest.approx(
        (-30.793189605316272, -0.043171211223460487, -0.64887647749645083, 0.0),
        rel=2.0e-13,
        abs=2.0e-13,
    )
    assert len(residuals) == 16
    assert len(jacobian) == 16 * 11
    assert len(diagnostics) == 4
    assert fingerprint.startswith("sha256:")


@pytest.mark.parametrize("component_id", ("methane", "ethane"))
def test_joint_pure_jacobian_is_the_exact_provider_hessian_chain_rule(
    component_id: str,
) -> None:
    specification = SPECIFICATIONS[component_id]
    dataset = load_pure_saturation_dataset(component_id)
    capsule = _capsule(component_id)
    residuals, jacobian, _, fingerprint = native.evaluate(
        capsule, _payload(component_id), (0.0,) * 11
    )
    table = _provider_table(capsule)

    expected = [0.0] * (16 * 11)
    for row_index, row in enumerate(dataset.training_rows):
        liquid_volume = specification.molar_mass_kg_per_mol / row.liquid_density_kg_m3
        vapor_volume = 8.31446261815324 * row.temperature_k / row.pressure_pa
        liquid = _provider_phase(table, row.temperature_k, liquid_volume, specification.start)
        vapor = _provider_phase(table, row.temperature_k, vapor_volume, specification.start)
        assert liquid.fingerprint.decode() == vapor.fingerprint.decode()
        assert liquid.fingerprint.decode() == fingerprint
        pressure_factor = 0.5 / row.pressure_pa
        residual_offset = 4 * row_index
        liquid_column = 3 + 2 * row_index
        vapor_column = liquid_column + 1
        for parameter, scale in enumerate(specification.parameter_scales):
            coordinate = 2 + parameter
            expected[(residual_offset + 0) * 11 + parameter] = (
                pressure_factor
                * -8.31446261815324
                * row.temperature_k
                * liquid.hessian[5 + coordinate]
                * scale
            )
            expected[(residual_offset + 1) * 11 + parameter] = (
                pressure_factor
                * -8.31446261815324
                * row.temperature_k
                * vapor.hessian[5 + coordinate]
                * scale
            )
            expected[(residual_offset + 2) * 11 + parameter] = (
                0.5 * (liquid.hessian[coordinate] - vapor.hessian[coordinate]) * scale
            )
        expected[(residual_offset + 0) * 11 + liquid_column] = (
            pressure_factor
            * -8.31446261815324
            * row.temperature_k
            * liquid.hessian[6]
            * liquid_volume
        )
        expected[(residual_offset + 1) * 11 + vapor_column] = (
            pressure_factor
            * -8.31446261815324
            * row.temperature_k
            * vapor.hessian[6]
            * vapor_volume
        )
        expected[(residual_offset + 2) * 11 + liquid_column] = (
            0.5 * liquid.hessian[1] * liquid_volume
        )
        expected[(residual_offset + 2) * 11 + vapor_column] = -0.5 * vapor.hessian[1] * vapor_volume
        expected[(residual_offset + 3) * 11 + liquid_column] = -0.5

    assert residuals
    assert tuple(jacobian) == pytest.approx(expected, rel=2.0e-14, abs=2.0e-14)


@pytest.mark.parametrize("component_id", ("methane", "ethane"))
def test_public_workflow_returns_strict_component_diagnostics(
    component_id: str,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    model = _model(component_id)
    dataset = load_pure_saturation_dataset(component_id)
    specification = SPECIFICATIONS[component_id]
    native_solve_general = native.solve_general
    general_solve_calls = 0

    def record_general_solve(*args: object) -> tuple[object, ...]:
        nonlocal general_solve_calls
        general_solve_calls += 1
        return native_solve_general(*args)

    monkeypatch.setattr(native, "solve_general", record_general_solve)
    assert not hasattr(native, "solve")
    result = fit_pure_saturation(
        model=model,
        dataset=dataset,
        specification=specification,
    )

    assert result.component_id == component_id
    assert result.solver_converged, result.failure_reasons
    assert result.numerically_converged, result.failure_reasons
    assert result.physically_valid, result.failure_reasons
    assert result.predictive_status == "NOT_ADJUDICATED_NO_APPROVED_HELD_OUT_CUTOFF"
    assert result.termination == "CONVERGENCE"
    assert result.solution_usable
    assert result.final_cost < result.initial_cost
    assert result.provider_fingerprint == model.parameter_fingerprint
    assert result.provider_fingerprint in result.compiled_problem_identity
    assert len(result.parameters) == 3
    assert any(abs(item.movement) > 1.0e-8 for item in result.parameters)
    assert not any(item.active_bound for item in result.parameters)
    assert result.jacobian.complete_columns
    assert result.jacobian.full_rank == 11
    assert result.jacobian.projected_parameter_rank == 3
    assert len(result.training_rows) == 4
    assert len(result.reporting_rows) == len(dataset.rows)
    assert tuple(row.temperature_k for row in result.reporting_rows) == tuple(
        row.temperature_k for row in dataset.rows
    )
    assert sum(row.training for row in result.reporting_rows) == 4
    assert all(row.physically_valid for row in result.reporting_rows if row.partition != "stress")
    assert tuple(row.partition for row in result.reporting_rows) == tuple(
        "training"
        if row.temperature_k in dataset.training_temperatures_k
        else "held_out"
        if row.temperature_k in dataset.held_out_temperatures_k
        else "stress"
        for row in dataset.rows
    )
    assert result.confirmation_parameter_scaled_max_delta <= 1.0e-5
    assert result.confirmation_cost_relative_delta <= 1.0e-8
    assert result.failure_reasons == ()
    assert general_solve_calls == 1


def test_propane_fit_preserves_distinct_statuses_at_the_frozen_pressure_gate() -> None:
    result = fit_pure_saturation(
        model=_model("propane"),
        dataset=load_pure_saturation_dataset("propane"),
        specification=PROPANE_SATURATION_FIT_V1,
    )

    assert result.solver_converged
    assert result.numerically_converged
    assert not result.physically_valid
    assert result.predictive_status == "NOT_ADJUDICATED_NO_APPROVED_HELD_OUT_CUTOFF"
    assert result.termination == "CONVERGENCE"
    assert result.confirmation_termination == "CONVERGENCE"
    assert result.iterations == 1090
    assert result.jacobian.full_rank == 11
    assert result.jacobian.projected_parameter_rank == 3
    assert not any(parameter.active_bound for parameter in result.parameters)
    assert result.confirmation_parameter_scaled_max_delta <= 1.0e-5
    assert result.confirmation_cost_relative_delta <= 1.0e-8
    failed_rows = tuple(row for row in result.reporting_rows if not row.physically_valid)
    assert "glos2004-propane-sat-120-k" in tuple(row.row_id for row in failed_rows)
    held_out_120 = next(row for row in failed_rows if row.row_id == "glos2004-propane-sat-120-k")
    assert held_out_120.partition == "held_out"
    assert (
        abs(held_out_120.raw_equilibrium_residuals[0]) / held_out_120.observed_pressure_pa
        > PROPANE_SATURATION_FIT_V1.reporting_pressure_scaled_residual_max
    )
    assert "training or reporting physical validity gate failed" in result.failure_reasons


def test_rank_deficient_parameter_jacobian_cannot_be_accepted(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    native_solve = native.solve_general

    def return_rank_deficient_result(*args: object) -> tuple[object, ...]:
        transport = list(native_solve(*args))
        transport[14] = 2
        return tuple(transport)

    monkeypatch.setattr(native, "solve_general", return_rank_deficient_result)
    result = fit_pure_saturation(
        model=_model("methane"),
        dataset=load_pure_saturation_dataset("methane"),
        specification=METHANE_SATURATION_FIT_V1,
    )

    assert result.jacobian.projected_parameter_rank == 2
    assert not result.solver_converged
    assert not result.numerically_converged
    assert not result.physically_valid
    assert (
        "training parameter Jacobian is rank deficient: 2 of 3 fitted parameter columns"
        in result.failure_reasons
    )


def test_rank_deficient_full_jacobian_cannot_be_accepted(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    native_solve = native.solve_general

    def return_rank_deficient_result(*args: object) -> tuple[object, ...]:
        transport = list(native_solve(*args))
        transport[11] = 10
        return tuple(transport)

    monkeypatch.setattr(native, "solve_general", return_rank_deficient_result)
    result = fit_pure_saturation(
        model=_model("methane"),
        dataset=load_pure_saturation_dataset("methane"),
        specification=METHANE_SATURATION_FIT_V1,
    )

    assert result.jacobian.full_rank == 10
    assert not result.solver_converged
    assert not result.numerically_converged
    assert not result.physically_valid


def test_low_cost_confirmation_agreement_uses_symmetric_relative_difference() -> None:
    result = fit_pure_saturation(
        model=_model("ethane"),
        dataset=load_pure_saturation_dataset("ethane"),
        specification=ETHANE_SATURATION_FIT_V1,
    )

    assert result.final_cost < 1.0
    assert 0.0 < result.confirmation_cost_relative_delta
    assert (
        result.confirmation_cost_relative_delta
        <= ETHANE_SATURATION_FIT_V1.confirmation_cost_relative_delta
    )


@pytest.mark.parametrize(
    "malformed_path",
    ("identity_text", "training_row", "payload_field", "reporting_row"),
)
def test_malformed_native_sequences_do_not_leak_references(
    malformed_path: str,
) -> None:
    capsule = _capsule("methane")
    payload = _payload("methane")
    dataset = load_pure_saturation_dataset("methane")
    reporting = tuple(workflow._row_payload(row) for row in dataset.rows)
    variables = (0.0,) * 11

    if malformed_path == "identity_text":
        identity = list(payload[0])
        identity[0] = object()
        tracked = tuple(identity)
        malformed_payload = (tracked, *payload[1:])
        call = partial(native.evaluate, capsule, malformed_payload, variables)
        expected_exception = ValueError
    elif malformed_path == "training_row":
        row = list(payload[1][0])
        row[0] = object()
        tracked = tuple(row)
        malformed_rows = (tracked, *payload[1][1:])
        malformed_payload = (payload[0], malformed_rows, *payload[2:])
        call = partial(native.evaluate, capsule, malformed_payload, variables)
        expected_exception = ValueError
    elif malformed_path == "payload_field":
        fields = list(payload)
        fields[2] = (object(), *payload[2][1:])
        tracked = tuple(fields)
        call = partial(native.evaluate, capsule, tracked, variables)
        expected_exception = ValueError
    else:
        row = list(reporting[0])
        row[0] = object()
        malformed_reporting = (tuple(row), *reporting[1:])
        tracked = tuple(malformed_reporting)
        call = partial(native.report_pure_saturation, capsule, payload, tracked, payload[2])
        expected_exception = ValueError

    gc.collect()
    reference_count = sys.getrefcount(tracked)
    for _ in range(32):
        with pytest.raises(expected_exception):
            call()
    gc.collect()

    assert sys.getrefcount(tracked) == reference_count


def test_generalized_workflow_preserves_accepted_methane_numerical_result() -> None:
    result = fit_pure_saturation(
        model=_model("methane"),
        dataset=load_pure_saturation_dataset("methane"),
        specification=METHANE_SATURATION_FIT_V1,
    )
    expected_reporting = (
        (100.0, 34626.07915160773, 436.84483289550474),
        (110.0, 88224.60866583801, 423.3969791365449),
        (120.0, 191083.41254773695, 409.3407614313209),
        (130.0, 366384.5067925305, 394.34823157350235),
        (140.0, 639981.9267634666, 377.95195367291694),
        (150.0, 1039603.4624909018, 359.43578939218077),
        (160.0, 1594405.452356535, 337.5806133970321),
        (170.0, 2334648.4404434026, 309.96652881184707),
        (180.0, 3290375.174877589, 270.4239126820564),
    )

    assert tuple(item.final for item in result.parameters) == pytest.approx(
        (0.9932081279826167, 3.717121437945618, 150.4888402511307),
        rel=2.0e-9,
        abs=2.0e-9,
    )
    assert result.initial_cost == pytest.approx(14340.021563034428, rel=2.0e-12)
    assert result.final_cost == pytest.approx(4.798586497669576e-6, rel=2.0e-9)
    assert result.jacobian.full_rank == 11
    assert result.jacobian.projected_parameter_rank == 3
    for observed, expected in zip(result.reporting_rows, expected_reporting, strict=True):
        assert (
            observed.temperature_k,
            observed.predicted_pressure_pa,
            observed.predicted_liquid_density_kg_m3,
        ) == pytest.approx(expected, rel=2.0e-9, abs=2.0e-9)


def test_generalized_workflow_preserves_accepted_ethane_numerical_result() -> None:
    result = fit_pure_saturation(
        model=_model("ethane"),
        dataset=load_pure_saturation_dataset("ethane"),
        specification=ETHANE_SATURATION_FIT_V1,
    )

    assert tuple(item.final for item in result.parameters) == pytest.approx(
        (1.6101710205193558, 3.524959232756593, 191.09459145171377),
        rel=1.0e-7,
        abs=1.0e-7,
    )
    assert result.initial_cost == pytest.approx(89326.40642623953, rel=2.0e-12)
    assert result.final_cost == pytest.approx(1.4622626154617253e-6, rel=2.0e-9)
    assert result.jacobian.full_rank == 11
    assert result.jacobian.projected_parameter_rank == 3
    assert result.predictive_status == ("NOT_ADJUDICATED_NO_APPROVED_HELD_OUT_CUTOFF")


def test_provider_callback_failure_is_returned_as_structured_fit_evidence(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    native_solve = native.solve_general

    def return_failure(*args: object) -> tuple[object, ...]:
        transport = list(native_solve(*args))
        transport[21] = "synthetic provider domain failure"
        return tuple(transport)

    monkeypatch.setattr(native, "solve_general", return_failure)

    result = fit_pure_saturation(
        model=_model("methane"),
        dataset=load_pure_saturation_dataset("methane"),
        specification=METHANE_SATURATION_FIT_V1,
    )

    assert not result.solver_converged
    assert any("synthetic provider domain failure" in reason for reason in result.failure_reasons)
