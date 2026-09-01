from __future__ import annotations

import math
from dataclasses import replace

import pytest

import epcsaft
from epcsaft import eos, equilibrium, regression
from tests.eos._cases import neutral_pure
from tests.equilibrium._cases import solved_pressure_case


def test_mixed_native_owner_fit_retains_order_and_both_jacobians() -> None:
    parameters, model, problem = solved_pressure_case()
    active = epcsaft.ActiveParameterSet(parameters, (parameters.parameter_specs[0].identity,))
    baseline = equilibrium.solve(model, problem, active_parameters=active)
    assert baseline.status == "evaluated", baseline.failure
    assert baseline.continuation_state is not None
    finite_state = baseline.continuation_state.phase(problem.phases[0].identity)
    finite_total = None
    for row, conserved in zip(
        problem.reaction_system.balance_matrix,
        problem.reaction_system.conserved_totals,
        strict=True,
    ):
        balance = math.fsum(
            coefficient * fraction
            for coefficient, fraction in zip(row, finite_state.mole_fractions, strict=True)
        )
        if abs(balance) > 1.0e-12:
            finite_total = conserved / balance
            break
    assert finite_total is not None
    anchored_phases = tuple(
        replace(
            phase,
            start=equilibrium.FinitePhaseStart(
                tuple(finite_total * fraction for fraction in state.mole_fractions),
                state.molar_volume_m3_per_mol,
            ),
        )
        if phase.amount_role == "finite"
        else phase
        for phase, state in zip(
            problem.phases,
            baseline.continuation_state.phases,
            strict=True,
        )
    )
    anchored = replace(
        problem,
        phases=anchored_phases,
        continuation_state=baseline.continuation_state,
    )
    liquid = baseline.phases[0]
    eos_owner = eos.compile_regression_evaluator(
        {
            "model": model,
            "family": "fixed-state",
            "temperature_k": 313.15,
            "molar_density_mol_m3": liquid.molar_density_mol_m3,
            "mole_fractions": liquid.mole_fractions,
            "row_identities": ("pressure_pa",),
        },
        active,
    )
    equilibrium_owner = equilibrium.compile_regression_evaluator(
        {"model": model, "problem": anchored}, active
    )
    eos_observation = model.state(
        T=313.15 * epcsaft.unit_registry.kelvin,
        rho=liquid.molar_density_mol_m3
        * epcsaft.unit_registry.mole
        / epcsaft.unit_registry.meter**3,
        x=liquid.mole_fractions,
    ).pressure.to("pascal").magnitude
    observations = (eos_observation,) + tuple(float(value) for value in baseline.values)
    result = regression.least_squares(
        (eos_owner, equilibrium_owner),
        observations,
        (active.values[0] + 0.05,),
        bounds=((active.values[0] - 5.0,), (active.values[0] + 5.0,)),
        affine_origin=(active.values[0] - 1.0,),
        affine_scale=(2.0,),
        residual_scales=(10.0,),
        maximum_iterations=20,
        maximum_elapsed_time_seconds=2.0,
    )
    assert result.usable, (result.message, result.failures)
    assert result.parameters == pytest.approx(active.values, abs=1.0e-8)
    assert result.predictions == pytest.approx(observations)
    assert result.residuals == pytest.approx((0.0,) * len(observations), abs=1.0e-8)
    assert tuple(
        (row["row_begin"], row["row_end"]) for row in result.evaluator_rows
    ) == ((0, 1), (1, len(observations)))
    assert all(math.isfinite(row[0]) for row in result.physical_jacobian)
    assert tuple(row[0] for row in result.optimizer_jacobian) == pytest.approx(
        tuple(row[0] * 2.0 / 10.0 for row in result.physical_jacobian)
    )
    assert result.residual_scale_metadata == "caller-defined discrepancy weights"
    assert result.covariance["status"] in {"available", "active_bound", "unavailable"}
    assert result.counters["uncached_jacobian_owner_attempts"] >= 4

    vector_scales = tuple(float(index + 1) for index in range(len(observations)))
    repeated = regression.curve_fit(
        (eos_owner, equilibrium_owner),
        observations,
        active.values,
        bounds=((active.values[0] - 5.0,), (active.values[0] + 5.0,)),
        affine_origin=(active.values[0] - 1.0,),
        affine_scale=(2.0,),
        residual_scales=vector_scales,
        maximum_iterations=4,
        maximum_elapsed_time_seconds=2.0,
    )
    assert repeated.usable
    assert repeated.parameters == pytest.approx(active.values)
    assert tuple(row[0] for row in repeated.optimizer_jacobian) == pytest.approx(
        tuple(
            row[0] * 2.0 / scale
            for row, scale in zip(repeated.physical_jacobian, vector_scales, strict=True)
        )
    )
    duplicate = regression.curve_fit(
        (eos_owner, equilibrium_owner),
        observations,
        active.values,
        bounds=((active.values[0] - 5.0,), (active.values[0] + 5.0,)),
        affine_origin=(active.values[0] - 1.0,),
        affine_scale=(2.0,),
        residual_scales=vector_scales,
        maximum_iterations=4,
        maximum_elapsed_time_seconds=2.0,
    )
    assert duplicate.parameters == repeated.parameters
    assert duplicate.predictions == repeated.predictions
    assert duplicate.physical_jacobian == repeated.physical_jacobian
    assert duplicate.optimizer_jacobian == repeated.optimizer_jacobian

    def fixed_prediction(parameter: float) -> float:
        margin = 1.0e-8 * step
        fixed = regression.least_squares(
            (eos_owner,),
            (eos_observation,),
            (parameter,),
            bounds=((parameter - margin,), (parameter + margin,)),
            maximum_iterations=2,
            maximum_elapsed_time_seconds=2.0,
        )
        assert fixed.usable, (fixed.message, fixed.failures)
        return fixed.predictions[0]

    step = 1.0e-5 * max(1.0, abs(active.values[0]))
    finite_difference = (
        fixed_prediction(active.values[0] + step)
        - fixed_prediction(active.values[0] - step)
    ) / (2.0 * step)
    assert result.physical_jacobian[0][0] == pytest.approx(
        finite_difference, rel=2.0e-5, abs=1.0e-6
    )

    active_bound = regression.least_squares(
        (eos_owner,),
        (eos_observation,),
        active.values,
        bounds=(active.values, (active.values[0] + 5.0,)),
        maximum_iterations=2,
        maximum_elapsed_time_seconds=2.0,
    )
    assert active_bound.usable
    assert active_bound.covariance["status"] == "active_bound"

    owner_limited = regression.least_squares(
        (eos_owner, equilibrium_owner),
        observations,
        active.values,
        maximum_uncached_owner_calls=1,
    )
    assert not owner_limited.usable
    assert owner_limited.parameters == owner_limited.predictions == ()
    assert owner_limited.failures[0]["code"] == "owner_call_limit"

    final_validation_limited = regression.least_squares(
        (eos_owner,),
        (eos_observation,),
        active.values,
        maximum_iterations=2,
        maximum_uncached_owner_calls=2,
    )
    assert not final_validation_limited.usable
    assert final_validation_limited.parameters == final_validation_limited.predictions == ()
    assert final_validation_limited.failures[-1]["evaluation_phase"] == "final_validation"
    assert final_validation_limited.failures[-1]["code"] == "owner_call_limit"

    invalid_initial = regression.least_squares(
        (eos_owner,),
        (eos_observation,),
        (active.values[0] + 2.0,),
        bounds=((active.values[0] - 1.0,), (active.values[0] + 1.0,)),
    )
    assert not invalid_initial.usable
    assert invalid_initial.parameters == invalid_initial.predictions == ()

    rejected = regression.least_squares(
        (eos_owner,),
        (eos_observation * 1000.0,),
        active.values,
        bounds=((-10000.0,), (10000.0,)),
        maximum_iterations=20,
        maximum_elapsed_time_seconds=2.0,
    )
    assert rejected.usable
    assert rejected.status == "limit_reached"
    assert rejected.counters["iteration_limit_reached"]
    assert rejected.counters["typed_trial_failures"] > 0
    assert any(
        failure["evaluation_phase"] == "ceres_trial"
        and failure["code"] == "eos_domain_rejection"
        for failure in rejected.failures
    )

    required_jacobian_failed = regression.least_squares(
        (eos_owner,),
        (eos_observation,),
        (active.values[0] + 0.05,),
        bounds=((active.values[0] - 5.0,), (active.values[0] + 5.0,)),
        maximum_iterations=20,
        maximum_elapsed_time_seconds=2.0,
        maximum_uncached_owner_calls=2,
    )
    assert not required_jacobian_failed.usable
    assert required_jacobian_failed.parameters == required_jacobian_failed.predictions == ()
    assert required_jacobian_failed.message == "a required owner Jacobian evaluation failed"
    assert required_jacobian_failed.failures[-1]["evaluation_phase"] == "ceres_trial"

    elapsed_limited = regression.least_squares(
        (eos_owner,),
        (eos_observation,),
        (active.values[0] + 0.05,),
        bounds=((active.values[0] - 5.0,), (active.values[0] + 5.0,)),
        maximum_iterations=1000,
        maximum_elapsed_time_seconds=1.0e-9,
    )
    assert elapsed_limited.usable
    assert elapsed_limited.status == "limit_reached"
    assert elapsed_limited.counters["elapsed_time_limit_reached"]

    starts = ((active.values[0] - 0.05,), (active.values[0] + 0.05,))

    def solve_start(start: tuple[float, ...]) -> regression.LeastSquaresResult:
        return regression.least_squares(
            (eos_owner,),
            (eos_observation,),
            start,
            bounds=((active.values[0] - 5.0,), (active.values[0] + 5.0,)),
            maximum_iterations=20,
            maximum_elapsed_time_seconds=2.0,
        )

    forward = {start: solve_start(start) for start in starts}
    reverse = {start: solve_start(start) for start in reversed(starts)}
    for start in starts:
        assert forward[start].parameters == reverse[start].parameters
        assert forward[start].predictions == reverse[start].predictions
        assert forward[start].physical_jacobian == reverse[start].physical_jacobian

    with pytest.raises(TypeError, match="requires an epcsaft EOS model"):
        eos.compile_regression_evaluator({}, active)


def test_pure_saturation_fit_preserves_owner_values_and_jacobian() -> None:
    parameters = epcsaft.Parameters.from_mapping(neutral_pure("propane"))
    model = epcsaft.Pure(parameters)
    active = epcsaft.ActiveParameterSet(
        parameters,
        (
            "component/propane/segment_count",
            "component/propane/segment_diameter",
        ),
    )
    block = model.observable_block(
        active,
        identity="eos.pure-saturation",
        T=300.0 * epcsaft.unit_registry.kelvin,
    )
    observations = tuple(float(value) for value in block.values if value is not None)[:3]
    owner = eos.compile_regression_evaluator(
        {
            "model": model,
            "family": "pure-saturation",
            "temperature_k": 300.0,
            "row_identities": (
                "saturation_pressure_pa",
                "vapor_molar_density_mol_per_m3",
                "liquid_molar_density_mol_per_m3",
            ),
        },
        active,
    )
    result = regression.least_squares(
        (owner,),
        observations,
        (active.values[0] + 0.01, active.values[1] - 0.01),
        bounds=((1.0, 1.0), (5.0, 5.0)),
        residual_scales=observations,
        maximum_iterations=20,
        maximum_elapsed_time_seconds=2.0,
    )

    assert result.usable
    assert result.parameters == pytest.approx(active.values, abs=5.0e-10)
    assert result.predictions == pytest.approx(observations)
    for actual, expected in zip(result.physical_jacobian, block.jacobian[:3], strict=True):
        assert actual == pytest.approx(expected)


def test_multiple_kij_coordinates_preserve_permutation_and_fit() -> None:
    parameters, model, problem = solved_pressure_case()
    baseline = equilibrium.solve(model, problem)
    assert baseline.status == "evaluated", baseline.failure
    liquid = baseline.phases[0]
    identities = (
        "pair/anion/cosolvent/k_ij",
        "pair/cation/solvent/k_ij",
        "pair/cosolvent/solvent/k_ij",
    )
    rows = ("pressure_pa",) + tuple(
        f"chemical_potential_over_rt/{component}" for component in model.component_ids
    )

    def solve(order: tuple[str, ...]) -> tuple[regression.LeastSquaresResult, tuple[tuple[float, ...], ...]]:
        active = epcsaft.ActiveParameterSet(parameters, order)
        block = model.observable_block(
            active,
            identity="eos.fixed-state",
            T=313.15 * epcsaft.unit_registry.kelvin,
            rho=liquid.molar_density_mol_m3
            * epcsaft.unit_registry.mole
            / epcsaft.unit_registry.meter**3,
            x=liquid.mole_fractions,
        )
        indexes = tuple(block.descriptor.row_identities.index(row) for row in rows)
        observations = tuple(
            float(block.values[index])
            for index in indexes
            if block.values[index] is not None
        )
        jacobian = tuple(
            tuple(block.jacobian[index] or ()) for index in indexes
        )
        owner = eos.compile_regression_evaluator(
            {
                "model": model,
                "family": "fixed-state",
                "temperature_k": 313.15,
                "molar_density_mol_m3": liquid.molar_density_mol_m3,
                "mole_fractions": liquid.mole_fractions,
                "row_identities": rows,
            },
            active,
        )
        result = regression.least_squares(
            (owner,),
            observations,
            tuple(value + 0.01 for value in active.values),
            bounds=((-0.5,) * 3, (0.5,) * 3),
            residual_scales=(abs(observations[0]), 1.0, 1.0, 1.0, 1.0),
            maximum_iterations=20,
            maximum_elapsed_time_seconds=2.0,
        )
        assert result.usable
        assert result.parameters == pytest.approx(active.values, abs=2.0e-11)
        for actual, expected in zip(result.physical_jacobian, jacobian, strict=True):
            assert actual == pytest.approx(expected)
        return result, jacobian

    ordered, _ = solve(identities)
    permuted, _ = solve(tuple(reversed(identities)))
    assert permuted.predictions == pytest.approx(ordered.predictions)
    for actual, expected in zip(
        permuted.physical_jacobian, ordered.physical_jacobian, strict=True
    ):
        assert actual == pytest.approx(tuple(reversed(expected)))
