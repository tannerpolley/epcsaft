from __future__ import annotations

import copy
import functools
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pytest
from chemical_equilibrium_cases import (
    base_system as _base_system,
)
from chemical_equilibrium_cases import (
    bind_records as _bind_record,
)
from chemical_equilibrium_cases import typed_problem as _typed_problem
from parameter_dictionaries import FIGIEL_REFERENCE_ELECTROLYTE_PARAMETERS

import epcsaft
from epcsaft.eos import _native_transport
from epcsaft import equilibrium
from epcsaft.equilibrium import _equilibrium

# gvbelov/Heterogeneous-Equilibrium commit c74b87545a3418262e8e38a7c7a2e31e1b12966a,
# test1.dat blob 5490caab3f7a3f8cdea0d1cb883cce4323902657, SHA-256 below.
# The element potentials and expected amounts are the independent 80-decimal
# calculation recorded in docs/research/2026-07-24-belov-aristova-chemical-equilibrium.md.
_BELOV_ELEMENT_POTENTIALS = (
    -10.0776517753450153959143507476192439107332024990741,
    44.8480568890243540815059574241188924984510845101041,
)
_BELOV_SOURCE_GIBBS = (
    -7.072468897041624,
    -28.765987200068952,
    -26.51300119486301,
    21.44955454258845,
    21.26761422175228,
    -34.35125913287487,
    -55.36542486137457,
    -17.10224518423043,
)
_HELD_WATER_IONIZATION_FINGERPRINT = (
    "sha256:8499a0cedeb7e8e34f7d70fbbe4c03180aea0018acae932e831736ce293e2aca"
)


def _packet_fingerprint(packet_root: Path) -> str:
    digest = hashlib.sha256()
    for path in sorted(item for item in packet_root.rglob("*") if item.is_file()):
        relative = path.relative_to(packet_root).as_posix().encode()
        digest.update(len(relative).to_bytes(8, "big"))
        digest.update(relative)
        digest.update(path.stat().st_size.to_bytes(8, "big"))
        digest.update(path.read_bytes())
    return digest.hexdigest()


@functools.cache
def _held_parameter_packet() -> Path:
    fixture = Path(__file__).parents[1] / "fixtures" / "data" / "held-cameretti-sadowski-2008"
    materialization = json.loads((fixture / "materialization.json").read_text())
    packet = fixture / "1"
    if _packet_fingerprint(packet) != materialization["packet_fingerprint"]:
        raise AssertionError("HELD parameter packet bytes changed")
    return packet


def _held_parameters(components: tuple[str, ...]) -> epcsaft.Parameters:
    return epcsaft.Parameters.from_bundle(
        _held_parameter_packet() / "parameters",
        components=components,
    )


def _water_ionization_source() -> dict[str, object]:
    path = Path(__file__).parents[2] / "data/reference/water_self_ionization_iapws_r11_07_2019.json"
    return json.loads(path.read_text())


def _iapws_p_kw(record: dict[str, object], temperature_k: float, density: float) -> float:
    correlation = record["correlation"]
    standard_state = record["standard_state"]
    assert isinstance(correlation, dict)
    assert isinstance(standard_state, dict)
    alpha = correlation["alpha"]
    beta = correlation["beta"]
    ideal = correlation["ideal_gas_p_kw_coefficients"]
    assert isinstance(alpha, list)
    assert isinstance(beta, list)
    assert isinstance(ideal, list)
    q_value = density * math.exp(
        alpha[0] + alpha[1] / temperature_k + alpha[2] * density ** (2.0 / 3.0) / temperature_k**2
    )
    p_kw_ideal = (
        ideal[0]
        + ideal[1] / temperature_k
        + ideal[2] / temperature_k**2
        + ideal[3] / temperature_k**3
    )
    return (
        -2.0
        * correlation["ion_coordination_number"]
        * (
            math.log10(1.0 + q_value)
            - q_value
            / (q_value + 1.0)
            * density
            * (beta[0] + beta[1] / temperature_k + beta[2] * density)
        )
        + p_kw_ideal
        + 2.0
        * math.log10(
            standard_state["standard_molality_mol_per_kg"]
            * standard_state["solvent_molar_mass_kg_per_mol"]
        )
    )


def _belov_aristova_gas_system() -> tuple[dict[str, object], tuple[float, ...]]:
    temperature_k = 2000.0
    pressure_pa = 100_000.0
    standard_pressure_pa = 101_325.0
    reactions = (
        (-2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        (-3.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0, -2.0, 1.0, 0.0, 0.0, 0.0),
        (-1.0, 0.0, 0.0, -1.0, 0.0, 1.0, 0.0, 0.0),
        (-2.0, 0.0, 0.0, -1.0, 0.0, 0.0, 1.0, 0.0),
        (-1.0, 0.0, 0.0, -2.0, 0.0, 0.0, 0.0, 1.0),
    )
    concentration_shift = math.log(8.31446261815324 * temperature_k / standard_pressure_pa)
    manufactured_gibbs = tuple(value + concentration_shift for value in _BELOV_SOURCE_GIBBS)
    ln_k = tuple(
        -sum(row[index] * manufactured_gibbs[index] for index in range(8)) for row in reactions
    )
    spec: dict[str, object] = {
        "species_ids": ("O", "O2", "O3", "C", "C2", "CO", "CO2", "C2O"),
        "charges": (0,) * 8,
        "provider_fingerprint": (
            "sha256:2f24904e5db64e96aa2ee96c53bd80863c9518b6c2ddd914e266181fe7459f84"
        ),
        "molar_masses_kg_per_mol": (
            0.016,
            0.032,
            0.048,
            0.012,
            0.024,
            0.028,
            0.044,
            0.040,
        ),
        "balance_matrix": (
            (0.0, 0.0, 0.0, 1.0, 2.0, 1.0, 1.0, 2.0),
            (1.0, 2.0, 3.0, 0.0, 0.0, 1.0, 2.0, 1.0),
        ),
        "reaction_matrix": reactions,
        "feed_amounts": (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0),
        "ln_k": ln_k,
        "temperature_k": temperature_k,
        "pressure_pa": pressure_pa,
    }
    _bind_record(spec)
    expected = (
        5.969938280702114e-17,
        5.249059580600665e-27,
        1.838054910648896e-47,
        1.7496236329525488e-5,
        4.996074090543112e-1,
        9.992323142113352e-1,
        4.453885591412746e-11,
        7.676856995870245e-4,
    )
    return spec, expected


def _amount_chart(
    charges: tuple[int, ...], coordinates: tuple[float, ...], trace_floor: float = 1.0e-12
) -> dict[str, object]:
    return _equilibrium._chemical_amount_chart(charges, coordinates, trace_floor)


def test_amount_chart_maps_neutral_logs_and_general_ionic_shares_exactly() -> None:
    neutral = _amount_chart((0, 0), (math.log(2.0), math.log(3.0)))
    assert neutral["amounts"] == pytest.approx((2.0, 3.0), abs=2.0e-15)
    assert neutral["charge_residual"] == 0.0

    charges = (0, 2, 1, -1, -2, 0)
    coordinates = (math.log(2.0), math.log(1.0 / 3.0), 0.0, math.log(3.0), math.log(4.0))
    ionic = _amount_chart(charges, coordinates)
    assert ionic["amounts"] == pytest.approx((3.0, 0.25, 1.5, 1.0, 0.5, 4.0), abs=2.0e-15)
    assert ionic["charge_residual"] == pytest.approx(0.0, abs=2.0e-15)


def test_amount_chart_has_exact_directional_first_and_second_derivatives() -> None:
    charges = (0, 2, 1, -1, -2, 0)
    center = (0.2, -0.7, 0.4, -0.3, 0.8)
    direction = (0.3, -0.2, 0.5, 0.1, -0.4)
    step = 2.0e-5

    lower = _amount_chart(
        charges,
        tuple(value - step * delta for value, delta in zip(center, direction, strict=True)),
    )
    result = _amount_chart(charges, center)
    upper = _amount_chart(
        charges,
        tuple(value + step * delta for value, delta in zip(center, direction, strict=True)),
    )

    dimension = len(center)
    jacobian = result["jacobian"]
    hessians = result["amount_hessians"]
    for species in range(len(charges)):
        first_directional = (upper["amounts"][species] - lower["amounts"][species]) / (2.0 * step)
        exact_first = sum(
            jacobian[species * dimension + column] * direction[column]
            for column in range(dimension)
        )
        assert first_directional == pytest.approx(exact_first, rel=3.0e-10, abs=3.0e-11)

        for row in range(dimension):
            jacobian_directional = (
                upper["jacobian"][species * dimension + row]
                - lower["jacobian"][species * dimension + row]
            ) / (2.0 * step)
            exact_second = sum(
                hessians[species * dimension * dimension + row * dimension + column]
                * direction[column]
                for column in range(dimension)
            )
            assert jacobian_directional == pytest.approx(exact_second, rel=2.0e-9, abs=5.0e-11)


def test_amount_chart_classifies_trace_floor_without_zeroing_species() -> None:
    trace_floor = 1.0e-10
    result = _amount_chart((0, 0), (math.log(0.5 * trace_floor), 0.0), trace_floor)

    assert result["minimum_amount"] == pytest.approx(0.5 * trace_floor, rel=2.0e-15)
    assert all(value > 0.0 for value in result["amounts"])


def test_amount_chart_rejects_one_sided_ionic_topology() -> None:
    with pytest.raises(ValueError):
        _amount_chart((1, 0), (0.0,))


def test_ionic_chart_pullback_spans_complete_physical_reaction_tangent() -> None:
    charges = np.asarray((0.0, 1.0, 1.0, -1.0, -1.0))
    chart = _amount_chart(tuple(charges.astype(int)), (0.2, -0.7, 0.4, -0.3))
    jacobian = np.asarray(chart["jacobian"], dtype=float).reshape(5, 4)
    balance = np.ones((1, 5))
    reactions = np.asarray(
        (
            (-2.0, 1.0, 0.0, 1.0, 0.0),
            (-2.0, 0.0, 1.0, 0.0, 1.0),
            (0.0, 1.0, -1.0, -1.0, 1.0),
        )
    )

    _, _, chart_null_rows = np.linalg.svd(balance @ jacobian)
    chart_feasible = jacobian @ chart_null_rows[1:].T
    chart_basis, _ = np.linalg.qr(chart_feasible)
    reaction_basis, _ = np.linalg.qr(reactions.T)
    principal_cosines = np.linalg.svd(
        chart_basis[:, :3].T @ reaction_basis[:, :3], compute_uv=False
    )

    assert np.linalg.matrix_rank(jacobian) == 4
    assert np.linalg.matrix_rank(chart_feasible) == 3
    assert np.linalg.matrix_rank(reactions) == 3
    assert np.max(np.abs(charges @ jacobian)) <= 2.0e-16
    assert np.max(np.abs(balance @ reactions.T)) <= 2.0e-16
    assert np.max(np.abs(charges @ reactions.T)) <= 2.0e-16
    assert principal_cosines == pytest.approx((1.0, 1.0, 1.0), abs=2.0e-15)

    chemical_potentials = 1.7 * balance[0] - 0.4 * charges
    assert reactions @ chemical_potentials == pytest.approx((0.0, 0.0, 0.0))
    chart_gradient = jacobian.T @ chemical_potentials
    chart_balance_gradient = jacobian.T @ balance[0]
    assert chart_gradient - 1.7 * chart_balance_gradient == pytest.approx(
        (0.0, 0.0, 0.0, 0.0), abs=2.0e-15
    )

    trace_chart = _amount_chart(
        tuple(charges.astype(int)),
        (0.0, -28.0, -27.0, math.log(2.0)),
        trace_floor=1.0e-18,
    )
    trace_jacobian = np.asarray(trace_chart["jacobian"], dtype=float).reshape(5, 4)
    _, _, trace_null_rows = np.linalg.svd(balance @ trace_jacobian)
    trace_pullback = trace_jacobian @ trace_null_rows[1:].T
    trace_singular_values = np.linalg.svd(trace_pullback, compute_uv=False)
    assert trace_singular_values[0] / trace_singular_values[-1] > 1.0e10


def test_incomplete_tangent_can_hide_a_failed_physical_affinity() -> None:
    balance = np.ones(5)
    charges = np.asarray((0.0, 1.0, 1.0, -1.0, -1.0))
    reactions = np.asarray(
        (
            (-2.0, 1.0, 0.0, 1.0, 0.0),
            (-2.0, 0.0, 1.0, 0.0, 1.0),
            (0.0, 1.0, -1.0, -1.0, 1.0),
        )
    )
    incomplete = np.vstack((balance, charges, reactions[:2]))
    _, _, right = np.linalg.svd(incomplete)
    chemical_potentials = right[-1]

    assert incomplete @ chemical_potentials == pytest.approx((0.0, 0.0, 0.0, 0.0), abs=8.0e-16)
    assert abs(reactions[2] @ chemical_potentials) > 0.5


def _manufactured_solve(
    spec: dict[str, object], options: dict[str, object] | None = None
) -> equilibrium.ChemicalEquilibriumResult:
    options = options or {}
    return equilibrium.chemical_equilibrium(
        equilibrium.IdealGasPhase(
            model_fingerprint=spec["provider_fingerprint"],
            reference_id="provider-helmholtz-coordinate-basis",
        ),
        spec["temperature_k"] * epcsaft.unit_registry.kelvin,
        spec["pressure_pa"] * epcsaft.unit_registry.pascal,
        _typed_problem(
            spec,
            minimum_amount_mol=options.get("trace_floor", 1.0e-12),
        ),
    )


@pytest.mark.parametrize(
    ("spec", "expected_amounts", "expected_volume"),
    (
        (_base_system(), (0.2, 0.8), None),
        (
            {
                **_base_system(),
                "balance_matrix": ((1.0, 2.0),),
                "reaction_matrix": ((-2.0, 1.0),),
                "molar_masses_kg_per_mol": (1.0, 2.0),
                "feed_amounts": (2.0, 0.0),
                "ln_k": (math.log(0.75),),
                "temperature_k": 300.0,
                "pressure_pa": 8.31446261815324 * 300.0,
            },
            (1.0, 0.5),
            1.5,
        ),
    ),
)
def test_manufactured_ideal_reactions_match_independent_analytic_states(
    spec: dict[str, object],
    expected_amounts: tuple[float, ...],
    expected_volume: float | None,
) -> None:
    _bind_record(spec)
    result = _manufactured_solve(spec)

    pressure_over_rt = spec["pressure_pa"] / (8.31446261815324 * spec["temperature_k"])
    volume = expected_volume or sum(expected_amounts) / pressure_over_rt
    diagnostics = result.diagnostics
    assert result.amounts_mol == pytest.approx(expected_amounts, rel=2.0e-8, abs=2.0e-10)
    assert result.volume_m3 == pytest.approx(volume, rel=2.0e-8)
    assert diagnostics.balance_inf_norm <= 1.0e-9
    assert diagnostics.charge_inf_norm <= 1.0e-12
    assert diagnostics.pressure_relative_residual <= 1.0e-8
    assert diagnostics.reaction_affinity_inf_norm <= 1.0e-7
    assert max(map(abs, diagnostics.reaction_affinity_residuals)) == pytest.approx(
        diagnostics.reaction_affinity_inf_norm
    )
    assert diagnostics.kkt_stationarity_inf_norm <= 1.0e-7
    assert max(map(abs, diagnostics.physical_stationarity_residuals)) == pytest.approx(
        diagnostics.kkt_stationarity_inf_norm
    )
    physical_rows, physical_columns = diagnostics.physical_constraint_shape
    derivative_dimension = len(diagnostics.derivative_coordinate_order)
    assert physical_columns == len(diagnostics.physical_stationarity_residuals)
    assert len(diagnostics.physical_equality_multipliers) == physical_rows
    assert len(diagnostics.physical_constraint_jacobian) == (physical_rows * physical_columns)
    physical_chart_rows, physical_chart_columns = diagnostics.physical_to_chart_shape
    assert physical_chart_rows == len(diagnostics.physical_lagrangian_gradient)
    assert physical_chart_columns == derivative_dimension
    assert len(diagnostics.physical_to_chart_jacobian) == (
        physical_chart_rows * physical_chart_columns
    )
    assert diagnostics.chart_physical_pullback_residual_inf_norm is not None
    assert diagnostics.chart_physical_pullback_residual_inf_norm <= 1.0e-10
    assert diagnostics.complementarity_inf_norm <= 1.0e-7
    assert diagnostics.objective_gradient_shape == (derivative_dimension,)
    assert diagnostics.constraint_jacobian_shape == (
        len(diagnostics.constraint_values),
        derivative_dimension,
    )
    assert diagnostics.lagrangian_hessian_shape == (
        derivative_dimension,
        derivative_dimension,
    )
    assert diagnostics.derivative_objective_basis == (
        "dimensionless_fixed_TP_A_plus_PV_plus_reference_over_RT"
    )
    assert diagnostics.derivative_constraint_basis == (
        "ordered_per_row_in_derivative_constraint_order"
    )
    assert len(diagnostics.derivative_constraint_order) == len(diagnostics.constraint_values)
    assert diagnostics.chart_stationarity_inf_norm <= 1.0e-7
    assert diagnostics.derivative_check_step is not None
    assert diagnostics.derivative_check_step > 0.0
    assert diagnostics.objective_gradient_check_relative_error is not None
    assert diagnostics.objective_gradient_check_relative_error <= 1.0e-5
    assert diagnostics.physical_objective_gradient_check_relative_error is not None
    assert diagnostics.physical_objective_gradient_check_relative_error <= 1.0e-5
    assert len(diagnostics.physical_objective_gradient) == (diagnostics.physical_to_chart_shape[0])
    assert diagnostics.constraint_jacobian_check_relative_error is not None
    assert diagnostics.constraint_jacobian_check_relative_error <= 1.0e-5
    assert diagnostics.lagrangian_hessian_check_relative_error is not None
    assert diagnostics.lagrangian_hessian_check_relative_error <= 2.0e-4
    assert diagnostics.reduced_hessian_check_relative_error is not None
    assert diagnostics.reduced_hessian_check_relative_error <= 2.0e-4
    assert diagnostics.derivative_check_worst_entry
    assert diagnostics.derivative_check_worst_relative_error is not None
    assert diagnostics.derivative_check_worst_analytic_value is not None
    assert diagnostics.derivative_check_worst_finite_difference_value is not None
    assert diagnostics.derivative_check_worst_step is not None
    assert diagnostics.kkt_root_shape == (
        diagnostics.kkt_dimension,
        diagnostics.kkt_dimension,
    )
    assert len(diagnostics.kkt_root_jacobian) == diagnostics.kkt_dimension**2
    assert diagnostics.kkt_root_jacobian_check_relative_error is not None
    assert diagnostics.kkt_root_jacobian_check_relative_error <= 2.0e-4
    assert sum(diagnostics.reduced_hessian_raw_inertia) == len(
        diagnostics.reduced_hessian_eigenvalues
    )
    reduced_dimension, full_dimension = diagnostics.reduced_hessian_nullspace_shape
    assert full_dimension == derivative_dimension
    assert reduced_dimension**2 == len(diagnostics.reduced_hessian)
    assert len(diagnostics.reduced_hessian_nullspace_basis) == (reduced_dimension * full_dimension)
    for row in range(reduced_dimension):
        for column in range(reduced_dimension):
            projected = math.fsum(
                diagnostics.reduced_hessian_nullspace_basis[row * full_dimension + left]
                * diagnostics.lagrangian_hessian[left * full_dimension + right]
                * diagnostics.reduced_hessian_nullspace_basis[column * full_dimension + right]
                for left in range(full_dimension)
                for right in range(full_dimension)
            )
            assert projected == pytest.approx(
                diagnostics.reduced_hessian[row * reduced_dimension + column],
                abs=2.0e-10,
            )
    assert diagnostics.minimum_amount_mol is not None
    assert diagnostics.trace_floor_mol == 1.0e-12
    assert diagnostics.minimum_amount_mol > diagnostics.trace_floor_mol


def test_manufactured_nonconvex_saddle_recovers_certified_lower_minimum() -> None:
    """A deterministic tangent displacement escapes a stationary saddle only."""

    spec = {**_base_system(), "ln_k": (0.0,)}
    _bind_record(spec)
    native = _equilibrium._chemical_solve_manufactured_nonconvex(
        spec,
        trace_floor=1.0e-12,
        max_iterations=200,
    )

    assert native["accepted"] is True
    assert native["amounts"][0] != pytest.approx(native["amounts"][1], abs=1.0e-5)
    pressure_over_rt = spec["pressure_pa"] / (8.31446261815324 * spec["temperature_k"])

    def objective(amounts: tuple[float, ...], volume: float) -> float:
        value = math.fsum(amount * (math.log(amount / volume) - 1.0) for amount in amounts)
        difference = amounts[0] - amounts[1]
        value += -2.3 * difference * difference + 2.0 * difference**4
        return value + pressure_over_rt * volume

    saddle_volume = 1.0 / pressure_over_rt
    saddle_coordinates = (
        math.log(0.5),
        math.log(0.5),
        math.log(saddle_volume),
    )
    lower = (
        math.log(0.1e-12),
        math.log(0.1e-12),
        math.log(saddle_volume) - 30.0,
    )
    upper = (0.0, 0.0, math.log(saddle_volume) + 30.0)
    tangent = (1.0 / math.sqrt(2.0), -1.0 / math.sqrt(2.0), 0.0)
    coarse_seed = _equilibrium._chemical_manufactured_recovery_displacement(
        saddle_coordinates, lower, upper, tangent, 1, 0
    )
    backtracked_seed = _equilibrium._chemical_manufactured_recovery_displacement(
        saddle_coordinates, lower, upper, tangent, 1, 3
    )
    assert objective(
        tuple(math.exp(value) for value in coarse_seed[:2]), saddle_volume
    ) > objective((0.5, 0.5), saddle_volume)
    assert objective(
        tuple(math.exp(value) for value in backtracked_seed[:2]), saddle_volume
    ) < objective((0.5, 0.5), saddle_volume)
    assert objective(tuple(native["amounts"]), native["volume_m3"]) < objective(
        (0.5, 0.5), saddle_volume
    )


def test_finite_but_inconsistent_callback_hessian_fails_spanning_audit() -> None:
    spec = {**_base_system(), "ln_k": (0.0,)}
    _bind_record(spec)

    native = _equilibrium._chemical_solve_manufactured_inconsistent_derivatives(
        spec,
        trace_floor=1.0e-12,
        max_iterations=200,
    )

    assert native["accepted"] is False
    assert native["derivative_check_worst_analytic_value"] is not None
    assert native["derivative_check_worst_finite_difference_value"] is not None
    assert native["derivative_check_worst_analytic_value"] != pytest.approx(
        native["derivative_check_worst_finite_difference_value"]
    )


def test_manufactured_nonconvex_search_retains_and_selects_distinct_basins() -> None:
    spec = {**_base_system(), "ln_k": (math.log(1.1),)}
    _bind_record(spec)

    native = _equilibrium._chemical_solve_manufactured_nonconvex(
        spec,
        trace_floor=1.0e-12,
        max_iterations=200,
    )

    search = native["search"]
    objectives = tuple(basin["objective"] for basin in search["basins"])
    assert len(objectives) == 2
    assert objectives[0] != pytest.approx(objectives[1], abs=1.0e-10)
    selected = search["selected_basin_ordinal"]
    assert search["selected_objective"] == pytest.approx(min(objectives))
    assert search["basins"][selected]["objective"] == pytest.approx(min(objectives))


@pytest.mark.parametrize(
    ("hessian", "inertia"),
    (
        ((2.0, 0.0, 0.0, 1.0), (2, 0, 0)),
        ((1.0, 0.0, 0.0, 0.0), (1, 1, 0)),
        ((1.0, 0.0, 0.0, -0.25), (1, 0, 1)),
    ),
)
def test_reduced_hessian_reports_spectrum_and_inertia(
    hessian: tuple[float, ...],
    inertia: tuple[int, int, int],
) -> None:
    evidence = _equilibrium._chemical_analyze_manufactured_reduced_hessian(hessian)

    assert tuple(evidence["inertia"]) == inertia
    assert tuple(evidence["raw_inertia"]) == inertia
    assert tuple(evidence["nullspace_basis"]) == pytest.approx((1.0, 0.0, 0.0, 1.0), abs=2.0e-13)
    assert tuple(evidence["eigenvalues"]) == pytest.approx(tuple(sorted(hessian[::3])), abs=2.0e-13)
    assert evidence["eigenvalue_tolerance"] == pytest.approx(
        1.0e-10 * max(1.0, *(abs(value) for value in hessian[::3]))
    )
    if inertia[2] > 0:
        assert evidence["curvature"] < 0.0


def test_reduced_hessian_separates_raw_spectrum_from_scaled_certification() -> None:
    evidence = _equilibrium._chemical_analyze_manufactured_reduced_hessian((1.0e20, 0.0, 0.0, 1.0))

    assert tuple(evidence["inertia"]) == (2, 0, 0)
    assert tuple(evidence["raw_inertia"]) == (1, 1, 0)


def test_reduced_hessian_nullspace_is_invariant_to_constraint_row_scaling() -> None:
    hessian = (2.0, 0.25, 0.0, 0.25, 1.5, 0.0, 0.0, 0.0, 0.75)
    baseline = _equilibrium._chemical_analyze_manufactured_reduced_hessian(
        hessian, (1.0, -2.0, 0.5), 1
    )
    scaled = _equilibrium._chemical_analyze_manufactured_reduced_hessian(
        hessian, (1.0e-18, -2.0e-18, 0.5e-18), 1
    )

    assert tuple(scaled["inertia"]) == tuple(baseline["inertia"])
    assert tuple(scaled["eigenvalues"]) == pytest.approx(
        tuple(baseline["eigenvalues"]), rel=2.0e-12, abs=2.0e-12
    )


def test_inactive_inequality_must_not_remove_negative_curvature_direction() -> None:
    hessian = (1.0, 0.0, 0.0, -1.0)
    jacobian = (1.0, 0.0, 0.0, 1.0)

    equality_tangent = _equilibrium._chemical_analyze_manufactured_reduced_hessian(
        hessian, jacobian[:2], 1
    )
    incorrectly_active_tangent = _equilibrium._chemical_analyze_manufactured_reduced_hessian(
        hessian, jacobian, 2
    )

    assert equality_tangent["curvature"] < 0.0
    assert tuple(equality_tangent["inertia"]) == (0, 0, 1)
    assert tuple(incorrectly_active_tangent["inertia"]) == (0, 0, 0)


def test_recovery_displacement_uses_mixed_sign_bound_room() -> None:
    variables = (0.0, 0.0)
    lower = (-1.0, -0.02)
    upper = (1.0, 1.0)
    direction = (1.0, -1.0)

    positive = _equilibrium._chemical_manufactured_recovery_displacement(
        variables, lower, upper, direction, 1
    )
    negative = _equilibrium._chemical_manufactured_recovery_displacement(
        variables, lower, upper, direction, -1
    )

    assert positive == pytest.approx((0.005, -0.005))
    assert negative == pytest.approx((-0.25, 0.25))
    assert all(lower[i] < positive[i] < upper[i] for i in range(2))
    assert all(lower[i] < negative[i] < upper[i] for i in range(2))


def test_balance_retraction_removes_raw_log_seed_balance_error() -> None:
    spec = _base_system()
    _bind_record(spec)
    pressure_over_rt = spec["pressure_pa"] / (8.31446261815324 * spec["temperature_k"])
    saddle_volume = 1.0 / pressure_over_rt
    variables = (math.log(0.5), math.log(0.5), math.log(saddle_volume))
    lower = (math.log(0.1e-12), math.log(0.1e-12), math.log(saddle_volume) - 30.0)
    upper = (0.0, 0.0, math.log(saddle_volume) + 30.0)
    raw = _equilibrium._chemical_manufactured_recovery_displacement(
        variables,
        lower,
        upper,
        (1.0 / math.sqrt(2.0), -1.0 / math.sqrt(2.0), 0.0),
        1,
        0,
    )
    raw_balance_error = abs(math.exp(raw[0]) + math.exp(raw[1]) - 1.0)
    assert raw_balance_error > 1.0e-4

    retracted = _equilibrium._chemical_retract_manufactured_balance(
        spec, raw, lower, upper, 1.0e-12
    )
    assert retracted
    assert abs(math.exp(retracted[0]) + math.exp(retracted[1]) - 1.0) <= 1.0e-9
    assert retracted[2] == variables[2]


def test_manufactured_reaction_reports_exact_conditioned_implicit_sensitivities() -> None:
    spec = _base_system()
    _bind_record(spec)
    native = _equilibrium._chemical_equilibrium(None, spec, None, None, 1.0e-12)

    sensitivities = native["sensitivities"]
    assert sensitivities["kkt_rank"] == sensitivities["kkt_dimension"] == 3
    assert sensitivities["condition_number_inf"] < 1.0e6

    amount_derivatives = sensitivities["amount_derivatives"]
    assert amount_derivatives[0] == pytest.approx((0.2, 0.8), abs=2.0e-9)
    assert amount_derivatives[1] == pytest.approx((-0.16, 0.16), abs=2.0e-9)
    assert amount_derivatives[2] == pytest.approx((0.0, 0.0), abs=2.0e-12)

    expected_volume = 8.31446261815324 * 350.0 / 200_000.0
    assert sensitivities["volume_derivatives"] == pytest.approx(
        (
            expected_volume,
            0.0,
            -expected_volume / 200_000.0,
        ),
        rel=2.0e-8,
        abs=2.0e-13,
    )

    steps = (1.0e-5, 1.0e-5, 5.0)
    for parameter, step in enumerate(steps):
        perturbed: list[dict[str, object]] = []
        for direction in (-1.0, 1.0):
            trial = copy.deepcopy(spec)
            if parameter == 0:
                trial["feed_amounts"] = (
                    1.0 + direction * step,
                    0.0,
                )
            elif parameter == 1:
                trial["ln_k"] = (
                    spec["ln_k"][0] + direction * step,  # type: ignore[index]
                )
            else:
                trial["pressure_pa"] = spec["pressure_pa"] + direction * step
            _bind_record(trial)
            perturbed.append(_equilibrium._chemical_equilibrium(None, trial, None, None, 1.0e-12))
        finite_difference_amounts = tuple(
            (perturbed[1]["amounts"][species] - perturbed[0]["amounts"][species]) / (2.0 * step)
            for species in range(2)
        )
        finite_difference_volume = (perturbed[1]["volume_m3"] - perturbed[0]["volume_m3"]) / (
            2.0 * step
        )
        assert amount_derivatives[parameter] == pytest.approx(
            finite_difference_amounts,
            rel=2.0e-6,
            abs=2.0e-10,
        )
        assert sensitivities["volume_derivatives"][parameter] == pytest.approx(
            finite_difference_volume,
            rel=2.0e-6,
            abs=2.0e-12,
        )


def test_manufactured_sensitivity_is_invariant_to_reaction_scaling() -> None:
    reaction_scale = 1.0e-11
    base = _base_system()
    _bind_record(base)
    scaled = copy.deepcopy(base)
    scaled["reaction_matrix"] = ((-reaction_scale, reaction_scale),)
    scaled["ln_k"] = (reaction_scale * math.log(4.0),)
    _bind_record(scaled)

    base_result = _equilibrium._chemical_equilibrium(None, base, None, None, 1.0e-12)[
        "sensitivities"
    ]
    scaled_result = _equilibrium._chemical_equilibrium(None, scaled, None, None, 1.0e-12)[
        "sensitivities"
    ]

    assert scaled_result["condition_number_inf"] == pytest.approx(
        base_result["condition_number_inf"], rel=2.0e-10
    )
    assert scaled_result["amount_derivatives"][0] == pytest.approx(
        base_result["amount_derivatives"][0], abs=2.0e-9
    )
    assert tuple(
        reaction_scale * value for value in scaled_result["amount_derivatives"][1]
    ) == pytest.approx(base_result["amount_derivatives"][1], abs=2.0e-9)
    assert scaled_result["amount_derivatives"][2] == pytest.approx(
        base_result["amount_derivatives"][2], abs=2.0e-12
    )


def test_manufactured_sensitivity_tracks_exact_species_permutation() -> None:
    base = _base_system()
    _bind_record(base)
    permutation = (1, 0)
    permuted = copy.deepcopy(base)
    for field in (
        "species_ids",
        "charges",
        "molar_masses_kg_per_mol",
        "feed_amounts",
    ):
        values = base[field]
        permuted[field] = tuple(values[index] for index in permutation)
    permuted["balance_matrix"] = tuple(
        tuple(row[index] for index in permutation)
        for row in base["balance_matrix"]  # type: ignore[union-attr]
    )
    permuted["reaction_matrix"] = tuple(
        tuple(row[index] for index in permutation)
        for row in base["reaction_matrix"]  # type: ignore[union-attr]
    )
    _bind_record(permuted)

    base_result = _equilibrium._chemical_equilibrium(None, base, None, None, 1.0e-12)[
        "sensitivities"
    ]
    permuted_result = _equilibrium._chemical_equilibrium(None, permuted, None, None, 1.0e-12)[
        "sensitivities"
    ]

    assert permuted_result["parameter_order"] == base_result["parameter_order"]
    assert permuted_result["condition_number_inf"] == pytest.approx(
        base_result["condition_number_inf"], rel=2.0e-10
    )
    for parameter in range(3):
        assert permuted_result["amount_derivatives"][parameter] == pytest.approx(
            tuple(base_result["amount_derivatives"][parameter][index] for index in permutation),
            abs=2.0e-9,
        )


def test_manufactured_reaction_with_ill_conditioned_kkt_has_no_sensitivity() -> None:
    balance_separation = 3.0e-6
    spec = _base_system()
    spec["species_ids"] = ("A", "B", "C")
    spec["charges"] = (0, 0, 0)
    spec["molar_masses_kg_per_mol"] = (1.0, 1.0, 1.0)
    spec["balance_matrix"] = ((1.0, 1.0, 1.0 + balance_separation),)
    spec["reaction_matrix"] = ((-1.0, 1.0, 0.0),)
    spec["feed_amounts"] = (0.5, 0.5, 1.0)
    spec["ln_k"] = (0.0,)
    spec["temperature_k"] = 300.0
    spec["pressure_pa"] = 8.31446261815324 * 300.0
    _bind_record(spec)

    native = _equilibrium._chemical_equilibrium(None, spec, None, None, 1.0e-12)

    sensitivities = native["sensitivities"]
    assert sensitivities["kkt_rank"] == sensitivities["kkt_dimension"] == 4
    assert sensitivities["condition_number_inf"] > 1.0e6
    assert native["accepted"] is False
    criteria = {record["name"]: record for record in native["numerical_criteria"]}
    assert (
        criteria["kkt_condition_number_inf"]["value"]
        > criteria["kkt_condition_number_inf"]["limit"]
    )


@pytest.mark.parametrize("trace_floor", (1.0e-50, 1.0e-55))
def test_belov_aristova_gas_restriction_resolves_extreme_positive_traces(
    trace_floor: float,
) -> None:
    spec, expected = _belov_aristova_gas_system()

    result = _manufactured_solve(spec, {"trace_floor": trace_floor})

    amounts = result.amounts_mol
    assert tuple(map(math.log, amounts)) == pytest.approx(
        tuple(map(math.log, expected)), abs=2.0e-6
    )
    balances = spec["balance_matrix"]
    assert tuple(
        sum(row[index] * amounts[index] for index in range(8))
        for row in balances  # type: ignore[union-attr]
    ) == pytest.approx((2.0, 1.0), abs=1.0e-9)
    assert (8.31446261815324 * 2000.0 * sum(amounts) / result.volume_m3) == pytest.approx(
        100_000.0, rel=1.0e-8
    )
    concentration_shift = math.log(8.31446261815324 * 2000.0 / 101_325.0)
    potentials = tuple(
        _BELOV_SOURCE_GIBBS[index]
        + concentration_shift
        + math.log(amounts[index] / result.volume_m3)
        for index in range(8)
    )
    reactions = spec["reaction_matrix"]
    assert tuple(
        sum(row[index] * potentials[index] for index in range(8))
        for row in reactions  # type: ignore[union-attr]
    ) == pytest.approx((0.0,) * 6, abs=1.0e-7)
    assert tuple(
        potentials[index]
        + balances[0][index] * _BELOV_ELEMENT_POTENTIALS[0]  # type: ignore[index]
        + balances[1][index] * _BELOV_ELEMENT_POTENTIALS[1]  # type: ignore[index]
        for index in range(8)
    ) == pytest.approx((0.0,) * 8, abs=1.0e-7)
    diagnostics = result.diagnostics
    assert diagnostics.balance_inf_norm <= 1.0e-9
    assert diagnostics.pressure_relative_residual <= 1.0e-8
    assert diagnostics.reaction_affinity_inf_norm <= 1.0e-7
    assert diagnostics.kkt_stationarity_inf_norm <= 1.0e-7


@pytest.mark.parametrize(
    "variant", ("reaction_basis", "species_order", "conservation_gauge", "feed_scale")
)
def test_belov_aristova_trace_solution_is_coordinate_invariant(variant: str) -> None:
    spec, expected = _belov_aristova_gas_system()
    options: dict[str, object] = {"trace_floor": 1.0e-50}

    if variant == "reaction_basis":
        reactions = spec["reaction_matrix"]
        ln_k = spec["ln_k"]
        spec["reaction_matrix"] = (
            tuple(
                reactions[0][index] + reactions[1][index]  # type: ignore[index]
                for index in range(8)
            ),
            *reactions[1:],  # type: ignore[index]
        )
        spec["ln_k"] = (ln_k[0] + ln_k[1], *ln_k[1:])  # type: ignore[index]
        _bind_record(spec)
    elif variant == "species_order":
        permutation = (7, 5, 3, 1, 6, 4, 2, 0)
        for field in (
            "species_ids",
            "charges",
            "molar_masses_kg_per_mol",
            "feed_amounts",
        ):
            values = spec[field]
            spec[field] = tuple(values[index] for index in permutation)  # type: ignore[index]
        spec["balance_matrix"] = tuple(
            tuple(row[index] for index in permutation)
            for row in spec["balance_matrix"]  # type: ignore[union-attr]
        )
        spec["reaction_matrix"] = tuple(
            tuple(row[index] for index in permutation)
            for row in spec["reaction_matrix"]  # type: ignore[union-attr]
        )
        expected = tuple(expected[index] for index in permutation)
    elif variant == "conservation_gauge":
        balances = spec["balance_matrix"]
        spec["balance_matrix"] = (
            tuple(
                balances[0][index] + balances[1][index]  # type: ignore[index]
                for index in range(8)
            ),
            tuple(
                2.0 * balances[0][index] - balances[1][index]  # type: ignore[index]
                for index in range(8)
            ),
        )
        _bind_record(spec)
    else:
        scale = 3.0
        spec["feed_amounts"] = tuple(
            scale * value
            for value in spec["feed_amounts"]  # type: ignore[union-attr]
        )
        expected = tuple(scale * value for value in expected)
        _bind_record(spec)

    result = _manufactured_solve(spec, options)

    assert tuple(map(math.log, result.amounts_mol)) == pytest.approx(
        tuple(map(math.log, expected)), abs=3.0e-6
    )
    assert result.diagnostics.reaction_affinity_inf_norm <= 1.0e-7


def test_manufactured_charged_reaction_uses_exact_electroneutral_chart() -> None:
    temperature_k = 300.0
    spec = {
        **_base_system(),
        "species_ids": ("A", "C+", "D-"),
        "charges": (0, 1, -1),
        "molar_masses_kg_per_mol": (2.0, 1.0, 1.0),
        "balance_matrix": ((2.0, 1.0, 1.0), (0.0, 1.0, -1.0)),
        "reaction_matrix": ((-1.0, 1.0, 1.0),),
        "feed_amounts": (1.0, 0.0, 0.0),
        "ln_k": (math.log(1.0 / 3.0),),
        "temperature_k": temperature_k,
        "pressure_pa": 8.31446261815324 * temperature_k,
    }
    _bind_record(spec)

    result = _manufactured_solve(spec)

    assert result.amounts_mol == pytest.approx((0.5, 0.5, 0.5), rel=3.0e-8)
    assert result.volume_m3 == pytest.approx(1.5, rel=3.0e-8)
    assert result.diagnostics.charge_inf_norm <= 2.0e-15


def test_manufactured_ultra_trace_charged_share_stays_differentiable() -> None:
    """Both charged softmax references and explicit shares stay interior."""

    temperature_k = 300.0
    trace_floor = 1.0e-18
    records = tuple(
        {
            "source_id": f"manufactured:reaction-{index}",
            "reference_id": "provider-helmholtz-coordinate-basis",
            "reaction_orientation": "products_positive",
            "conversion_id": "already-provider-basis",
            "dimensionless": True,
            "temperature_k": temperature_k,
            "pressure_pa": 8.31446261815324 * temperature_k,
        }
        for index in range(2)
    )

    def solve(log_ratio_constants: tuple[float, float], floor: float) -> dict[str, object]:
        spec = {
            **_base_system(),
            "species_ids": ("C1+", "C2+", "D1-", "D2-"),
            "charges": (1, 1, -1, -1),
            "molar_masses_kg_per_mol": (1.0, 1.0, 1.0, 1.0),
            "balance_matrix": ((1.0, 1.0, 1.0, 1.0), (1.0, 1.0, -1.0, -1.0)),
            "reaction_matrix": ((1.0, -1.0, 0.0, 0.0), (0.0, 0.0, 1.0, -1.0)),
            "feed_amounts": (1.0, 0.0, 1.0, 0.0),
            "ln_k": log_ratio_constants,
            "temperature_k": temperature_k,
            "pressure_pa": 8.31446261815324 * temperature_k,
            "equilibrium_constant_records": records,
        }
        spec["conserved_totals"] = (2.0, 0.0)
        return _equilibrium._chemical_equilibrium(None, spec, None, None, floor)

    explicit_trace = solve((-40.0, 0.0), trace_floor)
    assert explicit_trace["amounts"][0] > trace_floor

    reference_trace = solve((40.0, 40.0), trace_floor)
    assert reference_trace["amounts"][1] > trace_floor
    assert reference_trace["amounts"][3] > trace_floor

    rejected_trace = solve((40.0, 40.0), 1.0e-12)
    assert min(rejected_trace["amounts"]) <= 1.0e-12


def test_manufactured_equilibrium_is_gauge_scale_and_reaction_basis_invariant() -> None:
    base = _base_system()
    base["feed_amounts"] = (1.0, 0.0)
    _bind_record(base)
    plain = _manufactured_solve(base)
    gauge = copy.deepcopy(base)
    gauge["balance_matrix"] = ((3.25, 3.25),)
    _bind_record(gauge)
    gauged = _manufactured_solve(gauge)
    assert gauged.amounts_mol == pytest.approx(plain.amounts_mol, rel=2.0e-10)

    scaled = copy.deepcopy(base)
    scaled["feed_amounts"] = (7.0, 0.0)
    _bind_record(scaled)
    scaled_result = _manufactured_solve(scaled)
    assert scaled_result.amounts_mol == pytest.approx(
        tuple(7.0 * value for value in plain.amounts_mol), rel=3.0e-8
    )
    assert scaled_result.volume_m3 == pytest.approx(7.0 * plain.volume_m3, rel=3.0e-8)

    three = {
        **_base_system(),
        "species_ids": ("A", "B", "C"),
        "charges": (0, 0, 0),
        "molar_masses_kg_per_mol": (1.0, 1.0, 1.0),
        "balance_matrix": ((1.0, 1.0, 1.0),),
        "reaction_matrix": ((-1.0, 1.0, 0.0), (0.0, -1.0, 1.0)),
        "feed_amounts": (1.0, 0.0, 0.0),
        "ln_k": (math.log(2.0), math.log(3.0)),
    }
    _bind_record(three)
    first_basis = _manufactured_solve(three)
    changed_basis = copy.deepcopy(three)
    changed_basis["reaction_matrix"] = ((-1.0, 0.0, 1.0), (0.0, -1.0, 1.0))
    changed_basis["ln_k"] = (math.log(6.0), math.log(3.0))
    _bind_record(changed_basis)
    second_basis = _manufactured_solve(changed_basis)
    assert first_basis.amounts_mol == pytest.approx(second_basis.amounts_mol, rel=3.0e-9)


def test_manufactured_reaction_solve_has_no_feed_scaled_amount_cap() -> None:
    target = (0.5, 500.0, 0.5)
    volume = sum(target)
    spec = {
        **_base_system(),
        "species_ids": ("A", "B", "C"),
        "charges": (0, 0, 0),
        "molar_masses_kg_per_mol": (1001.0, 1.0, 1.0),
        "balance_matrix": ((1000.0, 1.0, 0.0), (1.0, 0.0, 1.0)),
        "reaction_matrix": ((-1.0, 1000.0, 1.0),),
        "feed_amounts": (1.0, 0.0, 0.0),
        "ln_k": (1000.0 * math.log(target[1] / volume),),
        "temperature_k": 300.0,
        "pressure_pa": 8.31446261815324 * 300.0,
    }
    _bind_record(spec)

    result = _manufactured_solve(spec, {"trace_floor": 0.1})

    assert result.amounts_mol == pytest.approx(target, rel=2.0e-7, abs=2.0e-8)


def test_manufactured_nlp_has_exact_directional_gradient_jacobian_and_hessian() -> None:
    spec = _base_system()
    _bind_record(spec)
    pressure_over_rt = spec["pressure_pa"] / (8.31446261815324 * spec["temperature_k"])
    center = (math.log(0.3), math.log(0.7), math.log(1.0 / pressure_over_rt))
    direction = (0.2, -0.4, 0.3)
    multipliers = (0.37,)
    step = 2.0e-5

    def evaluate(variables: tuple[float, ...]) -> dict[str, object]:
        return _equilibrium._chemical_evaluate_manufactured_nlp(spec, variables, multipliers)

    lower = evaluate(
        tuple(value - step * delta for value, delta in zip(center, direction, strict=True))
    )
    result = evaluate(center)
    upper = evaluate(
        tuple(value + step * delta for value, delta in zip(center, direction, strict=True))
    )

    objective_directional = (upper["objective"] - lower["objective"]) / (2.0 * step)
    assert objective_directional == pytest.approx(
        sum(
            result["objective_gradient"][index] * direction[index]
            for index in range(len(direction))
        ),
        rel=2.0e-9,
        abs=2.0e-10,
    )
    constraint_directional = (upper["constraints"][0] - lower["constraints"][0]) / (2.0 * step)
    assert constraint_directional == pytest.approx(
        sum(
            result["constraint_jacobian"][index] * direction[index]
            for index in range(len(direction))
        ),
        rel=2.0e-10,
        abs=2.0e-11,
    )
    dimension = len(direction)
    for row in range(dimension):
        gradient_directional = (
            upper["lagrangian_gradient"][row] - lower["lagrangian_gradient"][row]
        ) / (2.0 * step)
        hessian_directional = sum(
            result["lagrangian_hessian"][row * dimension + column] * direction[column]
            for column in range(dimension)
        )
        assert gradient_directional == pytest.approx(hessian_directional, rel=3.0e-9, abs=3.0e-10)


def test_manufactured_inverse_log_packing_has_exact_directional_pullback() -> None:
    spec = _base_system()
    _bind_record(spec)
    center = (math.log(0.3), math.log(0.7), 0.15)
    direction = (0.2, -0.4, 0.3)
    multipliers = (0.37,)
    step = 2.0e-5

    def evaluate(variables: tuple[float, ...]) -> dict[str, object]:
        return _equilibrium._chemical_evaluate_manufactured_inverse_log_packing_nlp(
            spec, variables, multipliers
        )

    lower = evaluate(
        tuple(value - step * delta for value, delta in zip(center, direction, strict=True))
    )
    result = evaluate(center)
    upper = evaluate(
        tuple(value + step * delta for value, delta in zip(center, direction, strict=True))
    )

    objective_directional = (upper["objective"] - lower["objective"]) / (2.0 * step)
    assert objective_directional == pytest.approx(
        sum(
            result["objective_gradient"][index] * direction[index]
            for index in range(len(direction))
        ),
        rel=3.0e-9,
        abs=3.0e-10,
    )
    constraint_directional = (upper["constraints"][0] - lower["constraints"][0]) / (2.0 * step)
    assert constraint_directional == pytest.approx(
        sum(
            result["constraint_jacobian"][index] * direction[index]
            for index in range(len(direction))
        ),
        rel=3.0e-10,
        abs=3.0e-11,
    )
    dimension = len(direction)
    for row in range(dimension):
        gradient_directional = (
            upper["lagrangian_gradient"][row] - lower["lagrangian_gradient"][row]
        ) / (2.0 * step)
        hessian_directional = sum(
            result["lagrangian_hessian"][row * dimension + column] * direction[column]
            for column in range(dimension)
        )
        assert gradient_directional == pytest.approx(hessian_directional, rel=5.0e-9, abs=5.0e-10)
    hessian = result["lagrangian_hessian"]
    for row in range(dimension):
        for column in range(dimension):
            scale = max(1.0, abs(hessian[row * dimension + column]))
            assert (
                abs(hessian[row * dimension + column] - hessian[column * dimension + row])
                <= 2.0e-13 * scale
            )
    assert result["volume_m3"] == pytest.approx(math.exp(-center[-1]), rel=2.0e-15)
    kkt_rhs = result["kkt_backtransform_rhs"]
    kkt_solution = result["kkt_backtransform_solution"]
    kkt_dimension = len(kkt_solution)
    assert kkt_dimension == dimension + 1
    for row in range(dimension):
        reconstructed = (
            sum(
                hessian[row * dimension + column] * kkt_solution[column]
                for column in range(dimension)
            )
            + result["constraint_jacobian"][row] * kkt_solution[-1]
        )
        assert reconstructed == pytest.approx(kkt_rhs[row], rel=2.0e-12, abs=2.0e-13)
    constraint_reconstructed = sum(
        result["constraint_jacobian"][column] * kkt_solution[column] for column in range(dimension)
    )
    assert constraint_reconstructed == pytest.approx(kkt_rhs[-1], rel=2.0e-12, abs=2.0e-13)
    zero_rhs = _equilibrium._chemical_evaluate_manufactured_inverse_log_packing_nlp(
        spec,
        center,
        multipliers,
        zero_kkt_rhs=True,
    )
    assert zero_rhs["kkt_backtransform_rhs"] == [0.0] * kkt_dimension
    assert zero_rhs["kkt_backtransform_solution"] == [0.0] * kkt_dimension


def test_manufactured_solver_rejects_trace_false_success() -> None:
    spec = _base_system()
    trace = copy.deepcopy(spec)
    trace["ln_k"] = (math.log(1.0e-10),)
    _bind_record(trace)
    with pytest.raises(equilibrium.ChemicalEquilibriumError) as failed:
        _manufactured_solve(trace, {"trace_floor": 1.0e-8})
    diagnostics = failed.value.diagnostics
    assert diagnostics.search.selected_basin_ordinal is None
    native = _equilibrium._chemical_equilibrium(None, trace, None, None, 1.0e-8)
    assert min(native["amounts"]) <= 1.0e-8


def test_manufactured_charged_solution_is_species_order_invariant() -> None:
    temperature_k = 300.0
    original = {
        **_base_system(),
        "species_ids": ("A", "C+", "D-"),
        "charges": (0, 1, -1),
        "molar_masses_kg_per_mol": (2.0, 1.0, 1.0),
        "balance_matrix": ((2.0, 1.0, 1.0), (0.0, 1.0, -1.0)),
        "reaction_matrix": ((-1.0, 1.0, 1.0),),
        "feed_amounts": (1.0, 0.0, 0.0),
        "ln_k": (math.log(1.0 / 3.0),),
        "temperature_k": temperature_k,
        "pressure_pa": 8.31446261815324 * temperature_k,
    }
    _bind_record(original)
    baseline = _manufactured_solve(original)
    permutation = (2, 0, 1)
    permuted = copy.deepcopy(original)
    for field in (
        "species_ids",
        "charges",
        "molar_masses_kg_per_mol",
        "feed_amounts",
    ):
        values = original[field]
        permuted[field] = tuple(values[index] for index in permutation)
    permuted["balance_matrix"] = tuple(
        tuple(row[index] for index in permutation) for row in original["balance_matrix"]
    )
    permuted["reaction_matrix"] = tuple(
        tuple(row[index] for index in permutation) for row in original["reaction_matrix"]
    )
    _bind_record(permuted)
    reordered = _manufactured_solve(permuted)
    inverse = tuple(permutation.index(index) for index in range(len(permutation)))
    assert tuple(reordered.amounts_mol[index] for index in inverse) == pytest.approx(
        baseline.amounts_mol, rel=4.0e-8
    )


def _figiel_provider_model(
    components: tuple[str, ...] = ("water", "sodium-cation", "chloride-anion"),
) -> epcsaft.Mixture:
    if components != tuple(FIGIEL_REFERENCE_ELECTROLYTE_PARAMETERS["components"]):
        raise ValueError("the compact Figiel test dictionary has one canonical component order")
    parameters = epcsaft.Parameters.from_dictionary(FIGIEL_REFERENCE_ELECTROLYTE_PARAMETERS)
    return epcsaft.Mixture(parameters)


@pytest.mark.parametrize("polar", (False, True), ids=("nonpolar", "full-polar"))
def test_provider_phase_block_matches_independent_finite_differences(
    polar: bool,
) -> None:
    mapping: dict[str, object] = {
        "schema": "epcsaft.parameters",
        "schema_version": 1,
        "components": ("solvent", "cation", "anion"),
        "parameters": {
            "mw": (0.018, None, None),
            "m": (1.2, 1.0, 1.0),
            "s": (2.8, 2.8, 2.8),
            "e": (350.0, 230.0, 170.0),
            "z": (0, 1, -1),
            "d_born": (None, 3.4, 4.1),
            "f_solv": (1.5, 1.0, 1.0),
            "epsilon_r": (78.0, None, None),
            "k_ij": (
                (0.0, -0.3, -0.3),
                (-0.3, 0.0, 0.8),
                (-0.3, 0.8, 0.0),
            ),
            "sites": (
                {
                    "component_id": "solvent",
                    "site_id": "a",
                    "site_class": "donor",
                    "multiplicity": 1,
                },
                {
                    "component_id": "solvent",
                    "site_id": "b",
                    "site_class": "acceptor",
                    "multiplicity": 1,
                },
            ),
            "association": (
                {
                    "component_id_a": "solvent",
                    "site_id_a": "a",
                    "component_id_b": "solvent",
                    "site_id_b": "b",
                    "association_energy_over_k": 2425.7,
                    "association_volume": 0.04509,
                },
            ),
        },
        "options": {
            "epsilon_r_ion": 8.0,
            "a_dh": 7.01,
            "permittivity_model": "ion-fraction-suppression",
        },
        "validity": {
            "kind": "reported-conditions",
            "temperature_min_k": 298.0,
            "temperature_max_k": 350.0,
            "pressure_min_pa": 10_000.0,
            "pressure_max_pa": 200_000.0,
            "ion_mole_fraction_max": 0.38,
        },
    }
    if polar:
        parameters = mapping["parameters"]
        options = mapping["options"]
        assert isinstance(parameters, dict)
        assert isinstance(options, dict)
        parameters["dipole_moment"] = (1.85, None, None)
        parameters["quadrupole_moment"] = (4.1, None, None)
        options["polar_formulation"] = "gross-vrabec-point-multipole"
    model = epcsaft.Mixture.from_dictionary(mapping)
    capsule = _native_transport(model)
    temperature_k = 313.15
    point = (0.98, 0.01, 0.01, 1.0e-3)

    def evaluate(values: tuple[float, float, float, float]) -> dict[str, object]:
        return _equilibrium._chemical_evaluate_provider_block(
            capsule,
            temperature_k,
            values[:3],
            values[3],
            model.parameter_fingerprint,
        )

    center = evaluate(point)
    gradient = tuple(center["gradient"])
    hessian = tuple(center["hessian"])
    directions = (
        ((1.0, 0.0, 0.0, 0.0), 2.0e-6),
        # The electrolyte phase model admits only electroneutral
        # states, so the ionic derivative is checked on its complete physical
        # tangent rather than through inadmissible single-ion perturbations.
        ((0.0, 1.0, 1.0, 0.0), 2.0e-7),
        ((0.0, 0.0, 0.0, 1.0), 2.0e-8),
    )
    for direction, step in directions:
        lower = tuple(value - step * delta for value, delta in zip(point, direction, strict=True))
        upper = tuple(value + step * delta for value, delta in zip(point, direction, strict=True))
        below = evaluate(lower)
        above = evaluate(upper)
        finite_difference_gradient = (float(above["value"]) - float(below["value"])) / (2.0 * step)
        analytic_directional_gradient = math.fsum(
            value * delta for value, delta in zip(gradient, direction, strict=True)
        )
        assert analytic_directional_gradient == pytest.approx(
            finite_difference_gradient, rel=3.0e-6, abs=3.0e-7
        )
        for row in range(4):
            finite_difference_hessian = (above["gradient"][row] - below["gradient"][row]) / (
                2.0 * step
            )
            analytic_directional_hessian = math.fsum(
                hessian[row * 4 + column] * direction[column] for column in range(4)
            )
            assert analytic_directional_hessian == pytest.approx(
                finite_difference_hessian, rel=3.0e-5, abs=3.0e-5
            )


def test_provider_certifies_multineutral_reactive_minimum() -> None:
    components = ("solvent", "cosolvent", "cation", "anion")
    model = epcsaft.Mixture.from_dictionary(
        {
            "schema": "epcsaft.parameters",
            "schema_version": 1,
            "components": components,
            "parameters": {
                "mw": (0.018, 0.030, None, None),
                "m": (1.2, 1.1, 1.0, 1.0),
                "s": (2.8, 3.0, 2.8, 2.8),
                "e": (350.0, 280.0, 230.0, 170.0),
                "z": (0, 0, 1, -1),
                "d_born": (None, None, 3.4, 4.1),
                "f_solv": (1.5, 1.2, 1.0, 1.0),
                "epsilon_r": (78.0, 40.0, None, None),
                "k_ij": tuple(tuple(0.0 for _ in components) for _ in components),
                "sites": (
                    {
                        "component_id": "solvent",
                        "site_id": "a",
                        "site_class": "donor",
                        "multiplicity": 1,
                    },
                    {
                        "component_id": "solvent",
                        "site_id": "b",
                        "site_class": "acceptor",
                        "multiplicity": 1,
                    },
                ),
                "association": (
                    {
                        "component_id_a": "solvent",
                        "site_id_a": "a",
                        "component_id_b": "solvent",
                        "site_id_b": "b",
                        "association_energy_over_k": 2425.7,
                        "association_volume": 0.04509,
                    },
                ),
            },
            "options": {
                "epsilon_r_ion": 8.0,
                "a_dh": 7.01,
                "permittivity_model": "ion-fraction-suppression",
            },
            "validity": {
                "kind": "reported-conditions",
                "temperature_min_k": 298.0,
                "temperature_max_k": 350.0,
                "pressure_min_pa": 1.0,
                "pressure_max_pa": 100_000_000.0,
            },
        }
    )
    phase = equilibrium.ProviderPhase(model, model.parameter_fingerprint, (1.0e-6, 0.74))
    target_amounts = (0.8, 0.2, 0.1, 0.1)
    target_volume = 0.025
    target = _equilibrium._chemical_evaluate_provider_block(
        _native_transport(model),
        313.15,
        target_amounts,
        target_volume,
        model.parameter_fingerprint,
    )
    target_gradient = tuple(float(value) for value in target["gradient"])
    target_ln_k = -target_gradient[0] + target_gradient[2] + target_gradient[3]
    problem = equilibrium.ChemicalEquilibriumProblem(
        species_ids=components,
        charges=(0, 0, 1, -1),
        molar_masses_kg_per_mol=(0.018, 0.030, 0.009, 0.009),
        balance_matrix=(
            (1.0, 0.0, 1.0, 0.0),
            (1.0, 0.0, 0.0, 1.0),
            (0.0, 1.0, 0.0, 0.0),
        ),
        conserved_totals=(0.9, 0.9, 0.2),
        reaction_matrix=((-1.0, 0.0, 1.0, 1.0),),
        feed_amounts_mol=(0.8, 0.2, 0.1, 0.1),
        equilibrium_constants=(
            equilibrium.ChemicalEquilibriumConstant(
                ln_value=target_ln_k,
                source_id="synthetic-multineutral-source",
                reference_id="provider-helmholtz-coordinate-basis",
                reaction_orientation="products_positive",
                conversion_id="already-provider-basis",
                dimensionless=True,
            ),
        ),
        strict_interior_amount_floor_mol=1.0e-30,
    )

    result = equilibrium.chemical_equilibrium(
        phase,
        313.15 * epcsaft.unit_registry.kelvin,
        float(target["pressure_pa"]) * epcsaft.unit_registry.pascal,
        problem,
    )

    assert min(result.amounts_mol) > problem.strict_interior_amount_floor_mol
    assert result.amounts_mol[1] == pytest.approx(0.2, abs=2.0e-10)
    assert result.amounts_mol[2] == pytest.approx(result.amounts_mol[3], abs=2.0e-10)


def _held_water_ionization_problem() -> tuple[
    epcsaft.Mixture, dict[str, object], dict[str, object]
]:
    source = _water_ionization_source()
    state = source["state"]
    standard_state = source["standard_state"]
    values = source["values"]
    assert isinstance(state, dict)
    assert isinstance(standard_state, dict)
    assert isinstance(values, dict)
    components = ("water", "hydronium-cation", "hydroxide-anion")
    parameters = _held_parameters(components)
    model = epcsaft.Mixture(parameters)
    assert model.parameter_fingerprint == _HELD_WATER_IONIZATION_FINGERPRINT
    # The installed EOS evaluates the declared pure-solvent molecular mass.
    conversion = math.log(standard_state["standard_molality_mol_per_kg"] * 0.01801528)
    spec: dict[str, object] = {
        "species_ids": components,
        "charges": (0, 1, -1),
        "provider_fingerprint": model.parameter_fingerprint,
        "molar_masses_kg_per_mol": (0.01801528, 0.01902322, 0.01700734),
        "balance_matrix": ((2.0, 3.0, 1.0), (1.0, 1.0, 1.0)),
        "reaction_matrix": ((-2.0, 1.0, 1.0),),
        "feed_amounts": (1.0, 0.0, 0.0),
        "ln_k": (values["ln_kw"],),
        "equilibrium_constant_records": (
            {
                "source_id": source["source"]["id"],
                "reference_id": standard_state["id"],
                "reaction_orientation": "products_positive",
                "conversion_id": "source-standard-state-to-provider-neutral-reference",
                "dimensionless": True,
                "temperature_k": state["temperature_k"],
                "pressure_pa": state["pressure_pa"],
            },
        ),
        "temperature_k": state["temperature_k"],
        "pressure_pa": state["pressure_pa"],
    }
    return (
        model,
        spec,
        {
            "id": standard_state["id"],
            "activity_scale_id": standard_state["activity_scale_id"],
            "log_activity_scale_factors": (0.0, conversion, conversion),
            "reference_pressure_pa": state["pressure_pa"],
            "source_reference_component_ids": components,
            "source_reference_solvent_composition": (1.0, 0.0, 0.0),
            "source_reference_activity_convention_id": ("molality-infinite-dilution-v1"),
            "source_reference_standard_molality_mol_per_kg": 1.0,
        },
    )


def _typed_standard_state(
    source_standard_state: dict[str, object],
) -> equilibrium.ChemicalStandardState:
    return equilibrium.ChemicalStandardState(
        id=source_standard_state["id"],
        activity_scale_id=source_standard_state["activity_scale_id"],
        log_activity_scale_factors=tuple(source_standard_state["log_activity_scale_factors"]),
        reference_pressure_pa=source_standard_state["reference_pressure_pa"],
        source_reference_component_ids=tuple(
            source_standard_state["source_reference_component_ids"]
        ),
        source_reference_solvent_composition=tuple(
            source_standard_state["source_reference_solvent_composition"]
        ),
        source_reference_activity_convention_id=source_standard_state[
            "source_reference_activity_convention_id"
        ],
        source_reference_standard_molality_mol_per_kg=(
            source_standard_state["source_reference_standard_molality_mol_per_kg"]
        ),
    )


def _provider_solve(
    model: epcsaft.Mixture,
    spec: dict[str, object],
    source_standard_state: dict[str, object] | None = None,
    sensitivity_request: equilibrium.ChemicalEquilibriumSensitivityRequest | None = None,
) -> equilibrium.ChemicalEquilibriumResult:
    standard_state = (
        None if source_standard_state is None else _typed_standard_state(source_standard_state)
    )
    return equilibrium.chemical_equilibrium(
        equilibrium.ProviderPhase(
            model=model,
            expected_parameter_fingerprint=spec["provider_fingerprint"],
            admissible_packing_fraction_interval=(1.0e-6, 0.74),
        ),
        spec["temperature_k"] * epcsaft.unit_registry.kelvin,
        spec["pressure_pa"] * epcsaft.unit_registry.pascal,
        _typed_problem(spec, source_standard_state=standard_state),
        sensitivity_request=sensitivity_request,
    )


def _provider_basis_spec(
    model: epcsaft.Mixture,
    source_spec: dict[str, object],
    source_standard_state: dict[str, object],
) -> tuple[dict[str, object], equilibrium.ChemicalEquilibriumResult]:
    transformed = _provider_solve(model, source_spec, source_standard_state)
    spec = copy.deepcopy(source_spec)
    spec["ln_k"] = transformed.ln_k_provider_basis
    spec["equilibrium_constant_records"] = tuple(
        {
            **record,
            "reference_id": "provider-helmholtz-coordinate-basis",
            "conversion_id": "already-provider-basis",
        }
        for record in spec["equilibrium_constant_records"]  # type: ignore[union-attr]
    )
    return spec, transformed


def _provider_native(
    model: epcsaft.Mixture,
    spec: dict[str, object],
    source_standard_state: dict[str, object] | None = None,
    packing_fraction_bounds: tuple[float, float] = (1.0e-6, 0.74),
) -> dict[str, object]:
    native_spec = copy.deepcopy(spec)
    if "conserved_totals" not in native_spec:
        balances = native_spec["balance_matrix"]
        feed = native_spec["feed_amounts"]
        native_spec["conserved_totals"] = tuple(
            math.fsum(
                row[index] * feed[index]  # type: ignore[index]
                for index in range(len(feed))  # type: ignore[arg-type]
            )
            for row in balances  # type: ignore[union-attr]
        )
    return _equilibrium._chemical_equilibrium(
        _native_transport(model),
        native_spec,
        source_standard_state,
        packing_fraction_bounds,
        1.0e-12,
    )


def test_iapws_2019_source_record_reproduces_release_equations() -> None:
    source = _water_ionization_source()
    state = source["state"]
    values = source["values"]
    assert isinstance(state, dict)
    assert isinstance(values, dict)
    assert _iapws_p_kw(source, 300.0, 1.0) == pytest.approx(13.906565, abs=5.0e-7)
    p_kw = _iapws_p_kw(
        source,
        state["temperature_k"],
        state["density_kg_per_m3"] / 1000.0,
    )
    assert p_kw == pytest.approx(values["p_kw"], abs=5.0e-13)
    assert -math.log(10.0) * p_kw == pytest.approx(values["ln_kw"], abs=1.0e-13)


def test_held_water_self_ionization_consumes_source_reference_and_provider() -> None:
    model, spec, source_standard_state = _held_water_ionization_problem()
    capsule = _native_transport(model)

    result = _provider_solve(
        model,
        spec,
        source_standard_state,
    )

    amounts = result.amounts_mol
    assert 2.0 * amounts[0] + 3.0 * amounts[1] + amounts[2] == pytest.approx(2.0, abs=2.0e-9)
    assert sum(amounts) == pytest.approx(1.0, abs=2.0e-9)
    assert amounts[1] == pytest.approx(amounts[2], abs=2.0e-12)
    final = _equilibrium._chemical_evaluate_provider_block(
        capsule,
        spec["temperature_k"],
        amounts,
        result.volume_m3,
        model.parameter_fingerprint,
    )
    provider_affinity = -2.0 * final["gradient"][0] + final["gradient"][1] + final["gradient"][2]
    assert provider_affinity - result.standard_offsets[0] == pytest.approx(
        spec["ln_k"][0], abs=2.0e-7
    )
    assert result.source_standard_state.id == ("iapws-molality-ions-mole-fraction-water")
    assert result.provider_reference_id == (
        "A_over_RT_reference_amount:n_ref=1mol:rho_ref=1mol_per_m3"
    )
    diagnostics = result.diagnostics
    assert diagnostics.reference_derivative_availability == 0
    assert diagnostics.reference_convergence_error <= 5.0e-5
    assert diagnostics.reference_representation_residual_inf_norm <= 1.0e-12
    assert (
        result.source_standard_state.activity_scale_id == source_standard_state["activity_scale_id"]
    )
    assert result.standard_offsets == pytest.approx((-21.20037599919826,), abs=2.0e-12)
    assert result.ln_k_provider_basis == pytest.approx((-53.42391841715451,), abs=2.0e-12)
    assert diagnostics.reaction_affinity_inf_norm <= 2.0e-12
    assert diagnostics.kkt_stationarity_inf_norm <= 2.0e-12
    assert result.sensitivity is None
    transfer = result.source_reference_transfer
    assert transfer is not None
    assert transfer.pressure_pa == spec["pressure_pa"]
    assert transfer.source_reference_pressure_pa == 100_000.0
    assert transfer.parameter_fingerprint == model.parameter_fingerprint
    for fingerprint in (
        transfer.topology_fingerprint,
        transfer.component_order_fingerprint,
        transfer.reference_state_fingerprint,
        transfer.domain_fingerprint,
        transfer.artifact_fingerprint,
    ):
        assert fingerprint.startswith("sha256:")
    assert result.source_standard_state.reference_pressure_pa == 100_000.0
    artifact = result.artifact_identity
    assert artifact.provider_record_sha256.startswith("sha256:")


def test_source_reference_pressure_is_distinct_from_trial_pressure() -> None:
    model, spec, source_standard_state = _held_water_ionization_problem()
    source_standard_state["reference_pressure_pa"] = 100_001.0

    value = _provider_solve(model, spec, source_standard_state)

    assert value.pressure_pa == 100_000.0
    assert value.source_reference_transfer is not None
    assert value.source_reference_transfer.pressure_pa == 100_000.0
    assert value.source_standard_state is not None
    assert value.source_standard_state.reference_pressure_pa == 100_001.0
    with pytest.raises(
        equilibrium.ChemicalEquilibriumError,
    ):
        _provider_solve(
            model,
            spec,
            source_standard_state,
            equilibrium.ChemicalEquilibriumSensitivityRequest(),
        )


def _held_active_parameter(
    family: str = "segment_diameter",
    value: float = 2.2740,
    unit: str = "angstrom",
) -> equilibrium.ChemicalEquilibriumActiveParameter:
    return equilibrium.ChemicalEquilibriumActiveParameter(
        family=family,
        identity="component",
        component_ids=("hydronium-cation",),
        value=value,
        unit=unit,
    )


def test_engine_model_basis_active_request_has_zero_lnk_cross_block() -> None:
    model, source_spec, source_standard_state = _held_water_ionization_problem()
    spec, _ = _provider_basis_spec(model, source_spec, source_standard_state)
    result = _provider_solve(
        model,
        spec,
        None,
        equilibrium.ChemicalEquilibriumSensitivityRequest(
            active_parameters=(_held_active_parameter(),)
        ),
    )
    assert result.sensitivity is not None


@pytest.mark.parametrize("variant", ("species_order", "reaction_orientation"))
def test_held_source_transform_is_coordinate_invariant(variant: str) -> None:
    model, spec, source_standard_state = _held_water_ionization_problem()
    if variant == "species_order":
        permutation = (2, 0, 1)
        components = tuple(spec["species_ids"][index] for index in permutation)
        parameters = _held_parameters(components)
        model = epcsaft.Mixture(parameters)
        for field in (
            "species_ids",
            "charges",
            "molar_masses_kg_per_mol",
            "feed_amounts",
        ):
            values = spec[field]
            spec[field] = tuple(values[index] for index in permutation)
        spec["balance_matrix"] = tuple(
            tuple(row[index] for index in permutation) for row in spec["balance_matrix"]
        )
        spec["reaction_matrix"] = tuple(
            tuple(row[index] for index in permutation) for row in spec["reaction_matrix"]
        )
        spec["provider_fingerprint"] = model.parameter_fingerprint
        scales = source_standard_state["log_activity_scale_factors"]
        source_standard_state["log_activity_scale_factors"] = tuple(
            scales[index] for index in permutation
        )
        reference_ids = source_standard_state["source_reference_component_ids"]
        source_standard_state["source_reference_component_ids"] = tuple(
            reference_ids[index] for index in permutation
        )
        solvent_composition = source_standard_state["source_reference_solvent_composition"]
        source_standard_state["source_reference_solvent_composition"] = tuple(
            solvent_composition[index] for index in permutation
        )
        water_index, cation_index, anion_index = (1, 2, 0)
    else:
        spec["reaction_matrix"] = tuple(
            tuple(-value for value in row) for row in spec["reaction_matrix"]
        )
        spec["ln_k"] = tuple(-value for value in spec["ln_k"])
        water_index, cation_index, anion_index = (0, 1, 2)

    result = _provider_solve(model, spec, source_standard_state)

    assert result.amounts_mol[water_index] == pytest.approx(0.9999999963727326, rel=2.0e-10)
    assert result.amounts_mol[cation_index] == pytest.approx(
        result.amounts_mol[anion_index], abs=2.0e-12
    )


def test_engine_model_manufactured_reaction_consumes_exact_phase_and_domain_blocks() -> None:
    model = _figiel_provider_model()
    capsule = _native_transport(model)
    temperature_k = 298.15
    target_amounts = (0.82, 0.08, 0.08)
    target_volume = 1.0e-3
    target = _equilibrium._chemical_evaluate_provider_block(
        capsule,
        temperature_k,
        target_amounts,
        target_volume,
        model.parameter_fingerprint,
    )
    ln_k = -target["gradient"][0] + target["gradient"][1] + target["gradient"][2]
    spec = {
        "species_ids": ("water", "sodium-cation", "chloride-anion"),
        "charges": (0, 1, -1),
        "molar_masses_kg_per_mol": (2.0, 1.0, 1.0),
        "provider_fingerprint": model.parameter_fingerprint,
        "balance_matrix": ((2.0, 1.0, 1.0), (0.0, 1.0, -1.0)),
        "reaction_matrix": ((-1.0, 1.0, 1.0),),
        "feed_amounts": (0.8, 0.1, 0.1),
        "ln_k": (ln_k,),
        "temperature_k": temperature_k,
        "pressure_pa": target["pressure_pa"],
    }
    _bind_record(spec)

    result = _provider_solve(model, spec)

    final = _equilibrium._chemical_evaluate_provider_block(
        capsule,
        temperature_k,
        result.amounts_mol,
        result.volume_m3,
        model.parameter_fingerprint,
    )
    assert 2.0 * result.amounts_mol[0] + result.amounts_mol[1] + result.amounts_mol[
        2
    ] == pytest.approx(1.8, abs=2.0e-9)
    assert -final["gradient"][0] + final["gradient"][1] + final["gradient"][2] == pytest.approx(
        ln_k, abs=2.0e-7
    )
    assert final["pressure_pa"] == pytest.approx(target["pressure_pa"], rel=1.0e-8)
    assert result.provider_parameter_fingerprint == model.parameter_fingerprint
    assert 1.0e-6 < result.diagnostics.packing_fraction < 0.74
    rejected = _provider_native(
        model,
        spec,
        packing_fraction_bounds=(0.73, 0.74),
    )
    assert rejected["accepted"] is False


def test_engine_model_pressure_sensitivity_matches_independent_resolves() -> None:
    model = _figiel_provider_model()
    capsule = _native_transport(model)
    temperature_k = 298.15
    target_amounts = (0.82, 0.08, 0.08)
    target_volume = 1.0e-3
    target = _equilibrium._chemical_evaluate_provider_block(
        capsule,
        temperature_k,
        target_amounts,
        target_volume,
        model.parameter_fingerprint,
    )
    spec = {
        "species_ids": ("water", "sodium-cation", "chloride-anion"),
        "charges": (0, 1, -1),
        "molar_masses_kg_per_mol": (2.0, 1.0, 1.0),
        "provider_fingerprint": model.parameter_fingerprint,
        "balance_matrix": ((2.0, 1.0, 1.0), (0.0, 1.0, -1.0)),
        "reaction_matrix": ((-1.0, 1.0, 1.0),),
        "feed_amounts": (0.8, 0.1, 0.1),
        "ln_k": (-target["gradient"][0] + target["gradient"][1] + target["gradient"][2],),
        "temperature_k": temperature_k,
        "pressure_pa": target["pressure_pa"],
    }
    _bind_record(spec)
    native = _provider_native(model, spec)
    sensitivities = native["sensitivities"]
    public = _provider_solve(
        model,
        spec,
        sensitivity_request=equilibrium.ChemicalEquilibriumSensitivityRequest(),
    )
    public_sensitivity = public.sensitivity
    assert public_sensitivity is not None
    assert sensitivities["parameter_fingerprint"] == model.parameter_fingerprint
    pressure_index = sensitivities["parameter_order"].index("pressure_pa")

    step_pa = 5.0
    perturbed: list[dict[str, object]] = []
    for direction in (-1.0, 1.0):
        trial = copy.deepcopy(spec)
        trial["pressure_pa"] = spec["pressure_pa"] + direction * step_pa
        _bind_record(trial)
        perturbed.append(_provider_native(model, trial))
    finite_difference_amounts = tuple(
        (perturbed[1]["amounts"][species] - perturbed[0]["amounts"][species]) / (2.0 * step_pa)
        for species in range(3)
    )
    finite_difference_volume = (perturbed[1]["volume_m3"] - perturbed[0]["volume_m3"]) / (
        2.0 * step_pa
    )
    assert sensitivities["amount_derivatives"][pressure_index] == pytest.approx(
        finite_difference_amounts,
        rel=2.0e-5,
        abs=2.0e-12,
    )
    assert sensitivities["volume_derivatives"][pressure_index] == pytest.approx(
        finite_difference_volume,
        rel=2.0e-5,
        abs=2.0e-13,
    )
