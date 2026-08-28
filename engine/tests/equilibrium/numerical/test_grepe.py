from __future__ import annotations

import math

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


def test_finite_phases_certify_balance_pressure_kkt_curvature_and_permutation() -> None:
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


@pytest.mark.slow
def test_finite_incipient_boundary_certifies_pressure_potentials_support_and_charge() -> None:
    _, model, problem = solved_pressure_case()
    result = equilibrium.solve(model, problem)

    assert result.status == "evaluated", result.failure
    pressure, vapor_fraction, partial_pressure = result.values
    assert pressure is not None and 2.5e6 < pressure < 4.0e6
    assert vapor_fraction is not None and 0.0 < vapor_fraction < 1.0
    assert partial_pressure == pytest.approx(pressure * vapor_fraction, rel=2.0e-13)
    liquid, vapor = result.phases
    assert vapor.supported_component_ids == ("solvent", "cosolvent")
    assert liquid.chemical_potential_over_rt[:2] == pytest.approx(
        vapor.chemical_potential_over_rt, abs=1.0e-6
    )
    assert liquid.pressure_pa == pytest.approx(pressure, rel=1.0e-8)
    assert vapor.pressure_pa == pytest.approx(pressure, rel=1.0e-8)
    assert abs(liquid.mole_fractions[2] - liquid.mole_fractions[3]) <= 1.0e-9
    evidence = dict(result.evidence)
    assert evidence["balance_inf_norm"] <= 1.0e-9
    assert evidence["pressure_relative_inf_norm"] <= 1.0e-8
    assert evidence["reaction_affinity_inf_norm"] <= 1.0e-6
    assert evidence["bubble_log_closure_abs"] <= 2.0e-10


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
