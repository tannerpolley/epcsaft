from __future__ import annotations

import math

import pytest

from epcsaft.equilibrium import (
    Constant,
    EvaluatedRow,
    ExpressionEvaluationError,
    LinearSum,
    Product,
    Ratio,
    RowReference,
    evaluate_expression,
)


def test_pressure_passthrough_preserves_owner_value_and_jacobian() -> None:
    result = evaluate_expression(
        RowReference("system_pressure"),
        [EvaluatedRow("system_pressure", 101_325.0, (2.5, -0.25))],
    )

    assert result.value == pytest.approx(101_325.0)
    assert result.jacobian == pytest.approx((2.5, -0.25))


def test_scaled_log_activity_expression_is_generic_ph_style_output() -> None:
    # pH = -log10(a_H+) = -ln(a_H+)/ln(10).  The logarithm is already an
    # owner-produced row; this tree only represents the declared scaling.
    result = evaluate_expression(
        Product(
            Constant(-1.0 / math.log(10.0)),
            RowReference("hydronium_log_activity"),
        ),
        [EvaluatedRow("hydronium_log_activity", math.log(1.0e-3), (0.5, -0.25))],
    )

    assert result.value == pytest.approx(3.0)
    assert result.jacobian == pytest.approx((-0.5 / math.log(10.0), 0.25 / math.log(10.0)))


def test_phase_recovery_is_a_species_neutral_amount_ratio() -> None:
    result = evaluate_expression(
        Ratio(RowReference("target_phase_amount"), RowReference("total_species_amount")),
        [
            EvaluatedRow("target_phase_amount", 0.35, (0.10, 0.0)),
            EvaluatedRow("total_species_amount", 1.40, (0.0, 0.20)),
        ],
    )

    assert result.value == pytest.approx(0.25)
    assert result.jacobian == pytest.approx((1.0 / 14.0, -1.0 / 28.0))


def test_linear_sum_and_product_use_exact_chain_rule() -> None:
    expression = LinearSum(
        terms=(
            (2.0, RowReference("x")),
            (-3.0, Product(RowReference("x"), RowReference("y"))),
        ),
        constant=4.0,
    )
    result = evaluate_expression(
        expression,
        [
            EvaluatedRow("x", 2.0, (1.0, 0.5)),
            EvaluatedRow("y", 3.0, (-2.0, 4.0)),
        ],
    )

    # 4 + 2x - 3xy = -10; d/dq = 2 dx - 3(y dx + x dy).
    assert result.value == pytest.approx(-10.0)
    assert result.jacobian == pytest.approx((5.0, -27.5))


@pytest.mark.parametrize(
    "rows, expression, message",
    [
        (
            [EvaluatedRow("x", 1.0, (1.0,)), EvaluatedRow("x", 2.0, (1.0,))],
            RowReference("x"),
            "duplicate row identity",
        ),
        (
            [EvaluatedRow("x", 1.0, (1.0,))],
            RowReference("missing"),
            "missing evaluated row",
        ),
        (
            [EvaluatedRow("x", None, None, status="failed")],
            RowReference("x"),
            "not evaluable",
        ),
    ],
)
def test_invalid_rows_and_references_fail_closed(
    rows: list[EvaluatedRow], expression: RowReference, message: str
) -> None:
    with pytest.raises(ExpressionEvaluationError, match=message):
        evaluate_expression(expression, rows)


def test_ratio_rejects_zero_or_near_zero_denominator() -> None:
    expression = Ratio(
        RowReference("numerator"), RowReference("denominator"), denominator_floor=1.0e-8
    )

    with pytest.raises(ExpressionEvaluationError, match="near zero"):
        evaluate_expression(
            expression,
            [
                EvaluatedRow("numerator", 1.0, (0.0,)),
                EvaluatedRow("denominator", 1.0e-9, (0.0,)),
            ],
        )


def test_nonfinite_rows_and_results_are_rejected() -> None:
    with pytest.raises(ExpressionEvaluationError, match="non-finite"):
        evaluate_expression(RowReference("x"), [EvaluatedRow("x", float("nan"), (1.0,))])

    with pytest.raises(ExpressionEvaluationError, match="intermediate"):
        evaluate_expression(
            Product(RowReference("x"), RowReference("y")),
            [
                EvaluatedRow("x", 1.0e308, (0.0,)),
                EvaluatedRow("y", 1.0e308, (0.0,)),
            ],
        )
