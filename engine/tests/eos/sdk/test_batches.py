from __future__ import annotations

from ._support import *


def test_native_sdk_aqueous_miac_bounded_batches_cancel_without_background_work() -> None:
    model = _electrolyte_model(("water", "sodium-cation", "bromide-anion"))
    _, table = _aqueous_miac_bounded_batch_table(model)
    assert table.evaluation_budget_size == ctypes.sizeof(NativeEvaluationBudget) == 24
    assert AqueousMiacBoundedBatchNativeSdkTable.evaluation_budget_size.offset == 240
    assert (
        AqueousMiacBoundedBatchNativeSdkTable.evaluate_aqueous_miac_kij_batch_bounded.offset == 248
    )
    assert (
        AqueousMiacBoundedBatchNativeSdkTable.evaluate_aqueous_miac_solvation_factor_batch_bounded.offset
        == 256
    )
    assert ctypes.sizeof(AqueousMiacBoundedBatchNativeSdkTable) == 264

    molalities = tuple(0.001 + (5.999 * index / 163.0) for index in range(164))
    molality_values = (ctypes.c_double * len(molalities))(*molalities)
    coordinates = (ctypes.c_double * 3)(-0.2, -0.3, 0.4)
    cancelled_results = (AqueousMiacKijResult * len(molalities))(
        *(AqueousMiacKijResult(struct_size=ctypes.sizeof(AqueousMiacKijResult)) for _ in molalities)
    )
    tiny_budget = NativeEvaluationBudget(
        struct_size=ctypes.sizeof(NativeEvaluationBudget),
        work_budget=0,
        deadline_ns=0,
    )
    started = time.perf_counter()
    cancelled_status = EvaluateAqueousMiacKijBatchBounded(
        table.evaluate_aqueous_miac_kij_batch_bounded
    )(
        table.model_context,
        model.parameter_fingerprint.encode(),
        298.15,
        100000.0,
        molality_values,
        len(molalities),
        coordinates,
        3,
        cancelled_results,
        len(molalities),
        ctypes.byref(tiny_budget),
    )
    cancellation_elapsed = time.perf_counter() - started
    assert cancelled_status == STATUS_CANCELLED
    assert all(result.status == STATUS_CANCELLED for result in cancelled_results)
    assert all(result.error for result in cancelled_results)
    assert cancellation_elapsed < 5.0

    deadline_results = (AqueousMiacKijResult * 1)(
        AqueousMiacKijResult(struct_size=ctypes.sizeof(AqueousMiacKijResult))
    )
    expired_budget = NativeEvaluationBudget(
        struct_size=ctypes.sizeof(NativeEvaluationBudget),
        work_budget=(1 << 64) - 1,
        deadline_ns=1,
    )
    deadline_status = EvaluateAqueousMiacKijBatchBounded(
        table.evaluate_aqueous_miac_kij_batch_bounded
    )(
        table.model_context,
        model.parameter_fingerprint.encode(),
        298.15,
        100000.0,
        (ctypes.c_double * 1)(1.0),
        1,
        coordinates,
        3,
        deadline_results,
        1,
        ctypes.byref(expired_budget),
    )
    assert deadline_status == STATUS_DEADLINE_EXCEEDED
    assert deadline_results[0].status == STATUS_DEADLINE_EXCEEDED
    assert deadline_results[0].error

    solvation_cancelled_results = (AqueousMiacSolvationFactorResult * len(molalities))(
        *(
            AqueousMiacSolvationFactorResult(
                struct_size=ctypes.sizeof(AqueousMiacSolvationFactorResult)
            )
            for _ in molalities
        )
    )
    solvation_budget = NativeEvaluationBudget(
        struct_size=ctypes.sizeof(NativeEvaluationBudget),
        work_budget=0,
        deadline_ns=0,
    )
    solvation_cancelled_status = EvaluateAqueousMiacSolvationFactorBatchBounded(
        table.evaluate_aqueous_miac_solvation_factor_batch_bounded
    )(
        table.model_context,
        model.parameter_fingerprint.encode(),
        298.15,
        100000.0,
        molality_values,
        len(molalities),
        1.35,
        solvation_cancelled_results,
        len(molalities),
        ctypes.byref(solvation_budget),
    )
    assert solvation_cancelled_status == STATUS_CANCELLED
    assert all(
        result.status == STATUS_CANCELLED and result.error for result in solvation_cancelled_results
    )

    bounded_results = (AqueousMiacKijResult * 1)(
        *(AqueousMiacKijResult(struct_size=ctypes.sizeof(AqueousMiacKijResult)) for _ in range(1))
    )
    sufficient_budget = NativeEvaluationBudget(
        struct_size=ctypes.sizeof(NativeEvaluationBudget),
        work_budget=(1 << 64) - 1,
        deadline_ns=0,
    )
    short_molalities = (ctypes.c_double * 1)(1.0)
    bounded_status = EvaluateAqueousMiacKijBatchBounded(
        table.evaluate_aqueous_miac_kij_batch_bounded
    )(
        table.model_context,
        model.parameter_fingerprint.encode(),
        298.15,
        100000.0,
        short_molalities,
        1,
        coordinates,
        3,
        bounded_results,
        1,
        ctypes.byref(sufficient_budget),
    )
    existing_results = _evaluate_aqueous_miac_kij_batch(
        ctypes.cast(ctypes.pointer(table), ctypes.POINTER(AqueousMiacBatchNativeSdkTable)).contents,
        model.parameter_fingerprint,
        (1.0,),
        (-0.2, -0.3, 0.4),
    )
    assert bounded_status == STATUS_OK
    for bounded, existing in zip(bounded_results, existing_results, strict=True):
        assert bytes(bounded) == bytes(existing)

    batch_table = ctypes.cast(
        ctypes.pointer(table), ctypes.POINTER(AqueousMiacBatchNativeSdkTable)
    ).contents
    scalar_kij = _evaluate_aqueous_miac_kij(
        batch_table,
        model.parameter_fingerprint,
        1.0,
        (-0.2, -0.3, 0.4),
    )
    assert bytes(existing_results[0]) == bytes(scalar_kij)

    batch_solvation = _evaluate_aqueous_miac_solvation_factor_batch(
        batch_table,
        model.parameter_fingerprint,
        (1.0,),
        1.35,
    )
    scalar_solvation = _evaluate_aqueous_miac_solvation_factor(
        batch_table,
        model.parameter_fingerprint,
        1.0,
        1.35,
    )
    assert bytes(batch_solvation[0]) == bytes(scalar_solvation)
