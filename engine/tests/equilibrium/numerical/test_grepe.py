from __future__ import annotations

import math
from dataclasses import replace

import pytest

import epcsaft
from epcsaft import equilibrium
from tests.equilibrium._cases import (
    finite_phase_case,
    ideal_uptake_case,
    provider_boundary_case,
    reverse_finite_phases,
    solved_pressure_case,
)


def test_finite_phases_retain_values_and_metrics_across_permutations() -> None:
    _, model, problem = finite_phase_case()
    forward = equilibrium.solve(model, problem)
    permuted_problem = reverse_finite_phases(problem)
    permuted = equilibrium.solve(model, permuted_problem)

    for result in (forward, permuted):
        assert result.status == "evaluated", result.failure
        evidence = dict(result.evidence)
        assert evidence["balance_inf_norm"] <= 1.0e-9
        assert evidence["pressure_relative_inf_norm"] <= 1.0e-8
        assert evidence["reaction_affinity_inf_norm"] <= 1.0e-6
        assert evidence["transfer_chemical_potential_inf_norm"] <= 1.0e-6
        assert evidence["kkt_stationarity_inf_norm"] <= 1.0e-7
        assert evidence["reduced_hessian_min_eigenvalue"] > 0.0
        assert evidence["local_evaluation_valid"] == 1
        assert evidence["equilibrium_value_valid"] == 1
        assert evidence["local_stability_status"] == "passed"
        assert evidence["phase_distinctness"] >= 1.0e-3
        assert evidence["declared_direction_rank"] == evidence["admissible_nullity"] == 3
        assert all(value is not None and value > 0.0 for value in result.values)
        for component in range(3):
            assert math.fsum(float(result.values[offset + component]) for offset in (0, 3)) == (
                pytest.approx(1.0, abs=1.0e-9)
            )
        for phase in result.phases:
            assert 1.0e-8 < phase.packing_fraction < 0.74
            assert phase.pressure_pa == pytest.approx(1.0e5, rel=1.0e-8)
            assert phase.mechanical_class == "strict_stable_local_minimum"

    forward_by_id = {
        phase.identity: state for phase, state in zip(problem.phases, forward.phases, strict=True)
    }
    permuted_by_id = {
        phase.identity: state
        for phase, state in zip(permuted_problem.phases, permuted.phases, strict=True)
    }
    for identity, expected in forward_by_id.items():
        assert permuted_by_id[identity].mole_fractions == pytest.approx(
            expected.mole_fractions, abs=2.0e-10
        )


def test_nonideal_finite_incipient_certifies_strict_vapor_stability() -> None:
    parameters, model, problem = solved_pressure_case()
    active = epcsaft.ActiveParameterSet(parameters, (parameters.parameter_specs[0].identity,))
    result = equilibrium.solve(model, problem, active_parameters=active)

    assert result.status == "evaluated", result.failure
    pressure, vapor_fraction, partial_pressure = result.values
    assert pressure is not None and 2.5e6 < pressure < 4.0e6
    assert vapor_fraction is not None and 0.0 < vapor_fraction < 1.0
    assert partial_pressure == pytest.approx(pressure * vapor_fraction, rel=2.0e-13)
    assert all(row.status == "evaluated" and row.jacobian is not None for row in result.rows)
    liquid, vapor = result.phases
    assert vapor.supported_component_ids == ("solvent", "cosolvent")
    assert liquid.chemical_potential_over_rt[:2] == pytest.approx(
        vapor.chemical_potential_over_rt, abs=1.0e-6
    )
    assert liquid.pressure_pa == pytest.approx(pressure, rel=1.0e-8)
    assert vapor.pressure_pa == pytest.approx(pressure, rel=1.0e-8)
    assert vapor.mechanical_class == "strict_stable_local_minimum"
    assert abs(liquid.mole_fractions[2] - liquid.mole_fractions[3]) <= 1.0e-9
    evidence = dict(result.evidence)
    assert evidence["equilibrium_value_valid"] == 1
    assert evidence["local_stability_status"] == "passed"
    assert evidence["balance_inf_norm"] <= 1.0e-9
    assert evidence["pressure_relative_inf_norm"] <= 1.0e-8
    assert evidence["reaction_affinity_inf_norm"] <= 1.0e-6
    assert math.isfinite(evidence["transfer_chemical_potential_inf_norm"])
    assert math.isfinite(evidence["kkt_stationarity_inf_norm"])
    assert evidence["bubble_log_closure_abs"] <= 2.0e-10


def test_vapor_stability_failure_retains_converged_value_payload() -> None:
    parameters, model, problem = solved_pressure_case()
    active = epcsaft.ActiveParameterSet(
        parameters,
        ("association/solvent/a/solvent/b/energy_over_k",),
        values=(10000.0,),
    )

    result = equilibrium.solve(model, problem, active_parameters=active)

    assert result.status == "non_evaluable"
    assert result.failure is not None and result.failure.code == "second_order_failure"
    assert result.physical_status == "failed"
    assert all(value is not None and math.isfinite(value) for value in result.values)
    assert all(value is None for value in result.jacobian)
    assert len(result.phases) == 2
    assert result.phases[1].mechanical_class == "stability_not_adjudicated"
    evidence = dict(result.evidence)
    assert evidence["equilibrium_value_valid"] == 1
    assert math.isfinite(evidence["bubble_log_closure_abs"])


@pytest.mark.parametrize(
    ("magnitude", "failure_code"),
    ((0.0, "singular_sensitivity"), (1.0e-10, "ill_conditioned_sensitivity")),
)
def test_sensitivity_failure_retains_complete_incipient_value_payload(
    magnitude: float,
    failure_code: str,
) -> None:
    parameters, model, problem = solved_pressure_case(ideal_vapor=True)
    active = epcsaft.ActiveParameterSet(parameters, (parameters.parameter_specs[0].identity,))
    baseline = equilibrium.solve(model, problem, active_parameters=active)
    assert baseline.status == "evaluated", baseline.failure
    pressure = float(baseline.values[0])
    feed = equilibrium.AffineFeed(
        problem.reaction_system.species_ids,
        problem.reaction_system.feed_amounts_mol,
        ((0.0,), (magnitude,), (0.0,), (0.0,)),
        (equilibrium.FeedCoordinate("sensitivity", (-1.0, 1.0)),),
    )
    boundary = equilibrium.IntensiveBoundary(
        "sensitivity-pressure", "sensitivity", None, "system.pressure", (), pressure, 0.0
    )
    candidate = replace(problem, feed=feed, intensive_boundaries=(boundary,))

    result = equilibrium.solve(model, candidate, active_parameters=active)

    assert result.failure is not None and result.failure.code == failure_code
    assert len(result.phases) == 2
    assert all(value is not None for value in result.values)
    assert all(value is None for value in result.jacobian)
    assert dict(result.evidence)["sensitivity_status"] == "unavailable"


def test_value_failure_is_typed_without_inventing_an_output() -> None:
    parameters, model, problem = solved_pressure_case(ideal_vapor=True)
    active = epcsaft.ActiveParameterSet(parameters, (parameters.parameter_specs[0].identity,))
    candidate = replace(
        problem,
        outputs=(replace(problem.outputs[2], coefficients=(1.0e308, 1.0e308)),),
    )

    result = equilibrium.solve(model, candidate, active_parameters=active)

    assert result.failure is not None and result.failure.code == "physical_rejection"
    assert result.values == (None,)
    assert result.jacobian == (None,)
    assert len(result.phases) == 2


def test_mixed_output_availability_is_retained_through_the_public_result() -> None:
    parameters, model, problem = solved_pressure_case(ideal_vapor=True)
    active = epcsaft.ActiveParameterSet(parameters, (parameters.parameter_specs[0].identity,))
    candidate = replace(
        problem,
        outputs=(
            replace(problem.outputs[2], coefficients=(1.0e308, 1.0e308)),
            problem.outputs[1],
        ),
    )

    result = equilibrium.solve(model, candidate, active_parameters=active)

    assert result.failure is not None and result.failure.code == "physical_rejection"
    assert result.values[0] is None
    assert result.values[1] is not None and 0.0 < result.values[1] < 1.0
    assert result.jacobian == (None, None)


def test_intensive_boundary_recovers_analytic_open_component_uptake() -> None:
    model, problem, expected = ideal_uptake_case()
    result = equilibrium.solve(model, problem)

    assert result.status == "evaluated", result.failure
    assert result.values == pytest.approx((expected,), rel=0.0, abs=2.0e-11)
    assert result.jacobian == ((),)
    assert result.phases[0].pressure_pa == pytest.approx(2.0e5, rel=1.0e-8)


@pytest.mark.slow
def test_intensive_boundary_total_jacobian_matches_public_parameter_resolves() -> None:
    parameters, model, problem, identity = provider_boundary_case()
    selected = epcsaft.ActiveParameterSet(parameters, (identity,))
    center = equilibrium.solve(model, problem, active_parameters=selected)
    step = 2.0e-5
    sides = tuple(
        equilibrium.solve(
            model,
            problem,
            active_parameters=epcsaft.ActiveParameterSet(
                parameters, (identity,), values=(selected.values[0] + sign * step,)
            ),
        )
        for sign in (-1.0, 1.0)
    )

    assert center.status == "evaluated", center.failure
    assert all(result.status == "evaluated" for result in sides)
    assert center.values[0] is not None and 0.0 < center.values[0] < 1.0
    assert center.jacobian[0] is not None
    exact = center.jacobian[0][0]
    centered = (float(sides[1].values[0]) - float(sides[0].values[0])) / (2.0 * step)
    assert abs(exact - centered) / max(abs(exact), abs(centered), 2.0e-9) <= 5.0e-5
