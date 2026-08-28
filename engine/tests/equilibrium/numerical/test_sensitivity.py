from __future__ import annotations

import pytest

import epcsaft
from epcsaft import equilibrium
from tests.eos._cases import AQUEOUS_ELECTROLYTE, NEUTRAL_TERNARY_IMMISCIBLE
from tests.equilibrium._cases import temperature_boundary_case


def _scaled_error(exact: float, finite_difference: float) -> float:
    return abs(exact - finite_difference) / max(
        abs(exact), abs(finite_difference), 2.0e-9
    )


def _solve_temperature_boundary(
    parameters: epcsaft.Parameters,
    model: epcsaft.Mixture,
    problem: equilibrium.GeneralReactiveEquilibriumProblem,
    identity: str,
    value: float,
) -> equilibrium.GeneralReactiveEquilibriumResult:
    return equilibrium.solve(
        model,
        problem,
        active_parameters=epcsaft.ActiveParameterSet(parameters, (identity,), values=(value,)),
    )


def test_coexistence_sensitivity_matches_centered_public_resolves() -> None:
    parameters = epcsaft.Parameters.from_mapping(AQUEOUS_ELECTROLYTE)
    model = epcsaft.Mixture(parameters)
    identity = "model/relative_permittivity/ion_fraction_suppression_coefficient"
    selected = epcsaft.ActiveParameterSet(parameters, (identity,))
    problem = equilibrium.PhaseEquilibriumProblem(
        equilibrium.FixedTCoexistence(
            298.15 * epcsaft.unit_registry.kelvin,
            equilibrium.DeclaredTopology(
                (
                    equilibrium.PhaseSpecification("liquid", equilibrium.AllComponents()),
                    equilibrium.PhaseSpecification("vapor", equilibrium.AllComponents()),
                )
            ),
            0,
            (0.96, 0.02, 0.02),
        )
    )

    def solve(value: float) -> equilibrium.PhaseEquilibriumResult:
        return equilibrium.solve(
            model,
            problem,
            active_parameters=epcsaft.ActiveParameterSet(parameters, (identity,), values=(value,)),
        )

    step = 1.0e-3
    center = solve(selected.values[0])
    lower = solve(selected.values[0] - step)
    upper = solve(selected.values[0] + step)
    assert isinstance(center.extent, equilibrium.IncipientBoundary)
    assert isinstance(lower.extent, equilibrium.IncipientBoundary)
    assert isinstance(upper.extent, equilibrium.IncipientBoundary)
    assert center.sensitivity is not None and center.sensitivity.status == "available"
    exact = center.sensitivity.output_parameter_jacobian[0][0]
    centered = (upper.extent.pressure_pa - lower.extent.pressure_pa) / (2.0 * step)
    assert _scaled_error(exact, centered) <= 5.0e-5


def test_declared_sensitivity_maps_identities_and_matches_centered_public_resolves() -> None:
    parameters = epcsaft.Parameters.from_mapping(NEUTRAL_TERNARY_IMMISCIBLE)
    model = epcsaft.Mixture(parameters)
    identity = "pair/alpha/beta/k_ij"
    selected = epcsaft.ActiveParameterSet(parameters, (identity,))
    problem = equilibrium.PhaseEquilibriumProblem(
        equilibrium.FixedTPFeed(
            300.0 * epcsaft.unit_registry.kelvin,
            1.0e5 * epcsaft.unit_registry.pascal,
            (1.0, 1.0, 1.0),
            equilibrium.DeclaredTopology(
                tuple(
                    equilibrium.PhaseSpecification("liquid", equilibrium.AllComponents())
                    for _ in range(3)
                )
            ),
        )
    )

    def solve(value: float) -> equilibrium.PhaseEquilibriumResult:
        return equilibrium.solve(
            model,
            problem,
            active_parameters=epcsaft.ActiveParameterSet(parameters, (identity,), values=(value,)),
        )

    def outputs(result: equilibrium.PhaseEquilibriumResult) -> tuple[float, ...]:
        assert isinstance(result.extent, equilibrium.FeedSplit)
        values = []
        for amount, phase in zip(result.extent.phase_amounts_mol, result.phases, strict=True):
            values.extend(amount * fraction for fraction in phase.mole_fractions)
            values.append(amount * phase.molar_volume_m3_per_mol)
        return tuple(values)

    step = 1.0e-5
    center = solve(selected.values[0])
    lower = outputs(solve(selected.values[0] - step))
    upper = outputs(solve(selected.values[0] + step))
    sensitivity = center.sensitivity
    assert sensitivity is not None and sensitivity.status == "available"
    assert sensitivity.parameters == selected.coordinates
    assert sensitivity.rank == sensitivity.dimension
    assert sensitivity.linear_solve_backward_error is not None
    assert sensitivity.linear_solve_backward_error <= 1.0e-8
    assert len(sensitivity.output_identities) == len(outputs(center))
    exact = tuple(row[0] for row in sensitivity.output_parameter_jacobian)
    centered = tuple(
        (right - left) / (2.0 * step) for left, right in zip(lower, upper, strict=True)
    )
    assert all(
        _scaled_error(actual, expected) <= 5.0e-5
        for actual, expected in zip(exact, centered, strict=True)
    )


@pytest.mark.slow
def test_solved_temperature_boundaries_certify_bubble_dew_and_total_derivatives() -> None:
    for reference_role in ("liquid", "vapor"):
        parameters, model, problem, identity = temperature_boundary_case(reference_role)
        selected = epcsaft.ActiveParameterSet(parameters, (identity,))
        step = 1.0e-5

        center = _solve_temperature_boundary(
            parameters, model, problem, identity, selected.values[0]
        )
        lower = _solve_temperature_boundary(
            parameters, model, problem, identity, selected.values[0] - step
        )
        upper = _solve_temperature_boundary(
            parameters, model, problem, identity, selected.values[0] + step
        )
        assert center.status == lower.status == upper.status == "evaluated"
        assert center.values[0] is not None and 295.0 < center.values[0] < 305.0
        for phase in center.phases:
            assert phase.pressure_pa == pytest.approx(1200.0, rel=1.0e-8)
        assert center.phases[1].chemical_potential_over_rt == pytest.approx(
            center.phases[0].chemical_potential_over_rt, abs=1.0e-6
        )
        for exact, low, high in zip(center.jacobian, lower.values, upper.values, strict=True):
            assert exact is not None and low is not None and high is not None
            centered = (high - low) / (2.0 * step)
            assert _scaled_error(exact[0], centered) <= 5.0e-5
