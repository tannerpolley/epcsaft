from __future__ import annotations

import math

import pytest
from tests.eos.sdk._support import _evaluate_neutral_reference, _neutral_reference_table

from epcsaft import Mixture, unit_registry

_COMPONENTS = (
    "neutral-alpha",
    "neutral-beta",
    "cation-monovalent",
    "cation-secondary",
    "anion-monovalent",
    "anion-divalent",
)
_NINE_COMPONENTS = (
    "neutral-alpha",
    "neutral-beta",
    "neutral-gamma",
    "cation-monovalent",
    "cation-secondary",
    "anion-monovalent",
    "anion-secondary",
    "anion-tertiary",
    "anion-divalent",
)
_CHARGES = {
    "neutral-alpha": 0,
    "neutral-beta": 0,
    "neutral-gamma": 0,
    "cation-monovalent": 1,
    "cation-secondary": 1,
    "anion-monovalent": -1,
    "anion-secondary": -1,
    "anion-tertiary": -1,
    "anion-divalent": -2,
}


def _generic_electrolyte(
    component_ids: tuple[str, ...] = _COMPONENTS,
    *,
    variant: str = "associating",
    temperature_k: float = 298.15,
    pressure_interval_pa: tuple[float, float] = (50000.0, 200000.0),
) -> dict[str, object]:
    values = {
        "neutral-alpha": (0.0180153, 1.2047, None, 353.95, None, 1.5),
        "neutral-beta": (0.04401, 1.5, 3.3, 220.0, None, 1.2),
        "neutral-gamma": (0.06108, 2.0, 3.5, 260.0, None, 1.1),
        "cation-monovalent": (None, 1.0, 2.8232, 230.0, 3.445, 1.0),
        "cation-secondary": (None, 1.0, 2.9, 210.0, 3.6, 1.0),
        "anion-monovalent": (None, 1.0, 2.756, 170.0, 4.1, 1.0),
        "anion-secondary": (None, 1.0, 2.8, 180.0, 4.0, 1.0),
        "anion-tertiary": (None, 1.0, 2.9, 185.0, 4.1, 1.0),
        "anion-divalent": (None, 1.0, 3.0, 190.0, 4.2, 1.0),
    }
    charges = tuple(_CHARGES[component_id] for component_id in component_ids)
    parameters: dict[str, object] = {
        "mw": tuple(values[component_id][0] for component_id in component_ids),
        "m": tuple(values[component_id][1] for component_id in component_ids),
        "s": tuple(values[component_id][2] for component_id in component_ids),
        "e": tuple(values[component_id][3] for component_id in component_ids),
        "z": charges,
        "d_born": tuple(values[component_id][4] for component_id in component_ids),
        "f_solv": tuple(values[component_id][5] for component_id in component_ids),
        "k_ij": tuple(
            tuple(
                0.0
                if left == right
                or charges[left] * charges[right] > 0
                or charges[left] == charges[right] == 0
                else 0.8
                if charges[left] * charges[right] < 0
                else -0.3
                for right in range(len(component_ids))
            )
            for left in range(len(component_ids))
        ),
        "correlations": (
            {
                "component_id": "neutral-alpha",
                "family": "segment_diameter",
                "form": "constant-plus-sum-of-exponentials",
                "constant": 2.7927,
                "terms": (
                    {"amplitude": 10.11, "exponent_coefficient_per_k": -0.01775},
                    {"amplitude": -1.417, "exponent_coefficient_per_k": -0.01146},
                ),
            },
            {
                "component_id": "neutral-alpha",
                "family": "relative_permittivity",
                "form": "constant",
                "constant": 78.09,
            },
            {
                "component_id": "neutral-beta",
                "family": "relative_permittivity",
                "form": "constant",
                "constant": 45.0,
            },
        ),
    }
    if "neutral-gamma" in component_ids:
        correlations = parameters["correlations"]
        assert isinstance(correlations, tuple)
        parameters["correlations"] = correlations + (
            {
                "component_id": "neutral-gamma",
                "family": "relative_permittivity",
                "form": "constant",
                "constant": 30.0,
            },
        )
    if variant != "ordinary":
        parameters["sites"] = (
            {
                "component_id": "neutral-alpha",
                "site_id": "a",
                "site_class": "donor",
                "multiplicity": 1,
            },
            {
                "component_id": "neutral-alpha",
                "site_id": "b",
                "site_class": "acceptor",
                "multiplicity": 1,
            },
        )
        parameters["association"] = (
            {
                "component_id_a": "neutral-alpha",
                "site_id_a": "a",
                "component_id_b": "neutral-alpha",
                "site_id_b": "b",
                "association_energy_over_k": 2425.7,
                "association_volume": 0.04509,
            },
        )
    options: dict[str, object] = {
        "epsilon_r_ion": 8.0,
        "a_dh": 7.01,
        "permittivity_model": "ion-fraction-suppression",
    }
    if variant == "polar":
        parameters["dipole_moment"] = tuple(
            1.7 if component_id == "neutral-alpha" else None for component_id in component_ids
        )
        parameters["quadrupole_moment"] = tuple(
            3.8 if component_id == "neutral-beta" else None for component_id in component_ids
        )
        options["polar_formulation"] = "gross-vrabec-point-multipole"
    return {
        "schema": "epcsaft.parameters",
        "schema_version": 1,
        "components": component_ids,
        "parameters": parameters,
        "options": options,
        "validity": {
            "kind": "fit-range",
            "temperature_min_k": temperature_k,
            "temperature_max_k": temperature_k,
            "pressure_min_pa": pressure_interval_pa[0],
            "pressure_max_pa": pressure_interval_pa[1],
            "ion_mole_fraction_max": 0.38,
        },
    }


def _matrix_rank(rows: tuple[tuple[float, ...], ...]) -> int:
    matrix = [list(row) for row in rows]
    rank = 0
    for column in range(len(matrix[0])):
        pivot = next(
            (row for row in range(rank, len(matrix)) if abs(matrix[row][column]) > 1e-12),
            None,
        )
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        pivot_value = matrix[rank][column]
        matrix[rank] = [value / pivot_value for value in matrix[rank]]
        for row in range(len(matrix)):
            if row == rank:
                continue
            factor = matrix[row][column]
            matrix[row] = [
                value - factor * pivot_item
                for value, pivot_item in zip(matrix[row], matrix[rank], strict=True)
            ]
        rank += 1
    return rank


@pytest.mark.parametrize("variant", ("ordinary", "associating", "polar"))
def test_neutral_reference_admits_generic_multineutral_electrolytes(variant: str) -> None:
    model = Mixture.from_dictionary(_generic_electrolyte(variant=variant))
    capsule, table = _neutral_reference_table(model)

    assert capsule is not None
    assert table.neutral_reference_basis_row_count == len(_COMPONENTS) - 1
    assert table.neutral_reference_result_size > 0
    assert table.evaluate_neutral_reference

    result, flat_basis, contractions, composition = _evaluate_neutral_reference(
        table,
        model.parameter_fingerprint,
        pressure_pa=100000.0,
    )
    rows = tuple(
        flat_basis[offset : offset + len(_COMPONENTS)]
        for offset in range(0, len(flat_basis), len(_COMPONENTS))
    )

    assert result.status == 0
    assert composition == pytest.approx((0.5, 0.5, 0.0, 0.0, 0.0, 0.0), abs=0.0)
    reference_state = model.state(
        T=298.15 * unit_registry.kelvin,
        P=100000.0 * unit_registry.pascal,
        x=composition,
        phase="liquid",
    )
    assert reference_state._log_fugacity_coefficient is not None
    assert contractions[:2] == pytest.approx(
        reference_state._log_fugacity_coefficient[:2], abs=2.0e-12
    )
    assert _matrix_rank(rows) == len(rows)
    charges = tuple(_CHARGES[component_id] for component_id in _COMPONENTS)
    assert all(
        math.fsum(value * charge for value, charge in zip(row, charges, strict=True)) == 0.0
        for row in rows
    )
    assert all(math.isfinite(value) for value in composition + contractions)
    assert result.parameter_fingerprint.decode() == model.parameter_fingerprint


def test_neutral_reference_basis_is_component_permutation_safe() -> None:
    canonical = Mixture.from_dictionary(_generic_electrolyte())
    permutation = (
        "anion-divalent",
        "neutral-beta",
        "cation-monovalent",
        "cation-secondary",
        "neutral-alpha",
        "anion-monovalent",
    )
    permuted = Mixture.from_dictionary(_generic_electrolyte(permutation))

    _, canonical_table = _neutral_reference_table(canonical)
    _, canonical_basis, canonical_contractions, canonical_composition = _evaluate_neutral_reference(
        canonical_table, canonical.parameter_fingerprint
    )
    _, permuted_table = _neutral_reference_table(permuted)
    _, permuted_basis, permuted_contractions, permuted_composition = _evaluate_neutral_reference(
        permuted_table, permuted.parameter_fingerprint
    )

    inverse = tuple(permutation.index(component_id) for component_id in _COMPONENTS)
    component_count = len(_COMPONENTS)
    remapped_basis = tuple(
        tuple(permuted_basis[row * component_count + index] for index in inverse)
        for row in range(component_count - 1)
    )
    canonical_rows = tuple(
        canonical_basis[offset : offset + component_count]
        for offset in range(0, len(canonical_basis), component_count)
    )

    assert remapped_basis == canonical_rows
    assert permuted_contractions == pytest.approx(canonical_contractions, abs=3e-11)
    assert tuple(permuted_composition[index] for index in inverse) == canonical_composition
