from __future__ import annotations

import math
from dataclasses import replace

import pytest

import epcsaft
from epcsaft import eos, equilibrium, regression
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
    assert result.diagnostics["residual_scale_metadata"] == "caller-defined discrepancy weights"
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

    invalid_initial = regression.least_squares(
        (eos_owner,),
        (eos_observation,),
        (active.values[0] + 2.0,),
        bounds=((active.values[0] - 1.0,), (active.values[0] + 1.0,)),
    )
    assert not invalid_initial.usable
    assert invalid_initial.parameters == invalid_initial.predictions == ()
    with pytest.raises(TypeError, match="requires an epcsaft EOS model"):
        eos.compile_regression_evaluator({}, active)
