from __future__ import annotations

import math

import pytest

from epcsaft.equilibrium import _equilibrium

_PHASES = ("aqueous", "organic")
_SPECIES = ("A+", "B+", "X-", "BC", "AC", "solvent")
_CHARGES = (1.0, 1.0, -1.0, 0.0, 0.0, 0.0)
_REFERENCE = (
    (0.8, 0.2, 1.0, 0.0, 0.0, 0.5),
    (0.0, 0.0, 0.0, 0.3, 0.1, 1.0),
)
_BALANCES = (
    (1.0, 0.0, 0.0, 0.0, 1.0, 0.0),  # total A
    (0.0, 1.0, 0.0, 1.0, 0.0, 0.0),  # total B
    (0.0, 0.0, 0.0, 1.0, 1.0, 0.0),  # carrier
    (0.0, 0.0, 1.0, 0.0, 0.0, 0.0),  # chloride
    (0.0, 0.0, 0.0, 0.0, 0.0, 1.0),  # solvent
)


def _directions() -> tuple[tuple[float, ...], tuple[float, ...]]:
    dimension = len(_PHASES) * len(_SPECIES)
    exchange = [0.0] * dimension
    exchange[0] = -1.0  # A+ aqueous
    exchange[1] = 1.0  # B+ aqueous
    exchange[3 + len(_SPECIES)] = -1.0  # BC organic
    exchange[4 + len(_SPECIES)] = 1.0  # AC organic
    transfer = [0.0] * dimension
    transfer[5] = -1.0  # solvent aqueous
    transfer[5 + len(_SPECIES)] = 1.0  # solvent organic
    return tuple(exchange), tuple(transfer)


def _compile(
    phases: tuple[str, ...] = _PHASES,
    reference: tuple[tuple[float, ...], ...] = _REFERENCE,
    directions: tuple[tuple[float, ...], ...] = _directions(),
) -> dict[str, object]:
    return dict(
        _equilibrium._compile_phase_direction_system(
            phases,
            _SPECIES,
            _CHARGES,
            _BALANCES,
            reference,
            directions,
            (math.log(2.0), 0.0)[: len(directions)],
        )
    )


def _matmul(matrix: list[list[float]], vector: list[float]) -> list[float]:
    return [sum(value * vector[column] for column, value in enumerate(row)) for row in matrix]


def _transpose_product(matrix: list[list[float]]) -> list[list[float]]:
    if not matrix:
        return [[0.0] * 0 for _ in range(0)]
    rows = len(matrix)
    columns = len(matrix[0])
    return [
        [sum(matrix[row][i] * matrix[row][j] for row in range(rows)) for j in range(columns)]
        for i in range(columns)
    ]


def test_native_phase_direction_compiler_handles_exchange_transfer_and_permutation() -> None:
    compiled = _compile()
    directions = [list(row) for row in _directions()]
    invariant_rows = [list(row) for row in compiled["invariant_rows"]]  # type: ignore[index]
    reference = [value for row in _REFERENCE for value in row]
    invariant_totals = compiled["invariant_totals"]
    assert compiled["direction_rank"] == 2
    assert compiled["invariant_rank"] == 10
    assert compiled["invariant_orthonormal_residual_inf_norm"] < 2.0e-12
    assert compiled["invariant_direction_residual_inf_norm"] < 2.0e-12
    assert compiled["global_balance_residual_inf_norm"] < 2.0e-12
    assert compiled["phase_charge_residual_inf_norm"] < 2.0e-12
    assert compiled["reference_neutrality_residual_inf_norm"] < 2.0e-12
    assert _matmul(invariant_rows, reference) == pytest.approx(invariant_totals)

    potential = compiled["reference_potentials"]
    assert _matmul(directions, potential) == pytest.approx((-math.log(2.0), 0.0))
    assert compiled["potential_equation_residual_inf_norm"] < 2.0e-12

    # A row-space projector is basis-independent.  Under a phase permutation
    # it must transform by the same coordinate permutation, even though the
    # orthonormal basis A itself is only defined up to row rotations/signs.
    permutation = list(range(len(_SPECIES), 2 * len(_SPECIES))) + list(range(len(_SPECIES)))
    permuted_reference = tuple(_REFERENCE[index] for index in (1, 0))
    permuted_directions = tuple(tuple(row[index] for index in permutation) for row in _directions())
    permuted = _compile(
        phases=("organic", "aqueous"),
        reference=permuted_reference,
        directions=permuted_directions,
    )
    first_projector = _transpose_product(invariant_rows)
    second_projector = _transpose_product(
        [list(row) for row in permuted["invariant_rows"]]  # type: ignore[index]
    )
    expected = [
        [first_projector[permutation[i]][permutation[j]] for j in range(len(permutation))]
        for i in range(len(permutation))
    ]
    for actual_row, expected_row in zip(second_projector, expected, strict=True):
        assert actual_row == pytest.approx(expected_row, abs=2.0e-12)


@pytest.mark.parametrize(
    ("directions", "message"),
    [
        (((1.0,) * 12, (2.0,) * 12), "full row rank"),
        (((-1.0, 1.0, 0.0, 0.0, 0.0, 0.0) + (0.0,) * 6,), "global balance conservation"),
    ],
)
def test_native_phase_direction_compiler_rejects_rank_or_conservation_failures(
    directions: tuple[tuple[float, ...], ...], message: str
) -> None:
    with pytest.raises(ValueError, match=message):
        _compile(directions=directions)


@pytest.mark.parametrize(
    ("defect", "accepted"),
    ((5.0e-13, True), (5.0e-11, False)),
)
def test_native_phase_direction_conservation_boundary(
    defect: float, accepted: bool
) -> None:
    direction = list(_directions()[0])
    direction[0] += defect
    if accepted:
        _compile(directions=(tuple(direction),))
    else:
        with pytest.raises(ValueError, match="conservation"):
            _compile(directions=(tuple(direction),))
