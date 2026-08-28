from __future__ import annotations

import pytest

import epcsaft
from epcsaft import equilibrium
from tests.eos._cases import STRONGLY_ASSOCIATING_BINARY


def _nonreactive_problem(
    component_ids: tuple[str, ...],
    amounts: tuple[float, ...],
) -> equilibrium.ChemicalEquilibriumProblem:
    return equilibrium.ChemicalEquilibriumProblem(
        component_ids,
        tuple(0 for _ in component_ids),
        tuple(0.05 for _ in component_ids),
        tuple(
            tuple(float(row == column) for column in range(len(component_ids)))
            for row in range(len(component_ids))
        ),
        amounts,
        (),
        amounts,
        (),
        1.0e-14,
    )


@pytest.mark.parametrize("reference_role", ["liquid", "vapor"])
def test_grepe_fixed_pressure_boundary_has_exact_bubble_and_dew_jacobians(
    reference_role: str,
) -> None:
    """Falsify role-specific dispatch or fixed-temperature regression partials."""

    parameters = epcsaft.Parameters.from_mapping(STRONGLY_ASSOCIATING_BINARY)
    model = epcsaft.Mixture(parameters)
    units = epcsaft.unit_registry
    parameter_identity = "pair/solute/solvent/k_ij"
    topology = equilibrium.DeclaredTopology(
        (
            equilibrium.PhaseSpecification("liquid", equilibrium.AllComponents()),
            equilibrium.PhaseSpecification("vapor", equilibrium.AllComponents()),
        )
    )
    bubble = equilibrium.solve(
        model,
        equilibrium.PhaseEquilibriumProblem(
            equilibrium.FixedTCoexistence(
                300.0 * units.kelvin, topology, 0, (0.5, 0.5)
            )
        ),
    )
    reference_index = 0 if reference_role == "liquid" else 1
    reference = (0.5, 0.5) if reference_index == 0 else bubble.phases[1].mole_fractions
    phases = tuple(
        equilibrium.ReactivePhase(
            phase.role,
            phase.role,
            "finite" if index == reference_index else "incipient",
            equilibrium.AllComponents(),
            equilibrium.ProviderModel((1.0e-9, 0.74), "installed-provider-eos"),
        )
        for index, phase in enumerate(topology.phases)
    )
    partner_index = 1 - reference_index
    problem = equilibrium.GeneralReactiveEquilibriumProblem(
        f"declared-{reference_role}-temperature",
        equilibrium.Solved(
            299.0 * units.kelvin,
            (295.0 * units.kelvin, 305.0 * units.kelvin),
            (296.0 * units.kelvin, 304.0 * units.kelvin),
        ),
        equilibrium.Fixed(bubble.extent.pressure_pa * units.pascal),
        phases,
        _nonreactive_problem(parameters.component_ids, reference),
        (),
        (
            equilibrium.EquilibriumOutput(
                "temperature", "system.temperature", "kelvin", "absolute-temperature"
            ),
            equilibrium.EquilibriumOutput(
                "partner-solute",
                "phase.mole_fraction",
                "dimensionless",
                "true-species-mole-fraction",
                phases[partner_index].identity,
                (1.0, 0.0),
            ),
        ),
        continuation_identity="declared-two-phase-temperature",
        feed=equilibrium.direct_amounts(parameters.component_ids, reference),
    )
    assert (
        equilibrium.general_reactive_equilibrium_problem_from_mapping(
            equilibrium.general_reactive_equilibrium_problem_to_mapping(problem)
        )
        == problem
    )

    def solve(value: float) -> equilibrium.GeneralReactiveEquilibriumResult:
        return equilibrium.solve(
            model,
            problem,
            active_parameters=epcsaft.ActiveParameterSet(
                parameters, (parameter_identity,), values=(value,)
            ),
        )

    center_value = epcsaft.ActiveParameterSet(
        parameters, (parameter_identity,)
    ).values[0]
    step = 1.0e-4
    center = solve(center_value)
    lower = solve(center_value - step)
    upper = solve(center_value + step)
    assert center.status == "evaluated", center.failure
    assert center.descriptor is not None
    assert center.descriptor.formulation == "temperature_boundary"
    assert center.values[0] == pytest.approx(300.0, abs=3.0e-8)
    for row, low, high in zip(center.jacobian, lower.values, upper.values, strict=True):
        assert row is not None and low is not None and high is not None
        assert row[0] == pytest.approx((high - low) / (2.0 * step), rel=2.0e-5)


def test_grepe_fixed_pressure_boundary_rejects_provider_domain() -> None:
    """Falsify accepting a certified root outside a declared Provider domain."""

    parameters = epcsaft.Parameters.from_mapping(STRONGLY_ASSOCIATING_BINARY)
    model = epcsaft.Mixture(parameters)
    units = epcsaft.unit_registry
    topology = equilibrium.DeclaredTopology(
        (
            equilibrium.PhaseSpecification("liquid", equilibrium.AllComponents()),
            equilibrium.PhaseSpecification("vapor", equilibrium.AllComponents()),
        )
    )
    bubble = equilibrium.solve(
        model,
        equilibrium.PhaseEquilibriumProblem(
            equilibrium.FixedTCoexistence(
                300.0 * units.kelvin, topology, 0, (0.5, 0.5)
            )
        ),
    )
    phases = tuple(
        equilibrium.ReactivePhase(
            phase.role,
            phase.role,
            "finite" if index == 0 else "incipient",
            equilibrium.AllComponents(),
            equilibrium.ProviderModel((1.0e-8, 1.0e-7), "installed-provider-eos"),
        )
        for index, phase in enumerate(topology.phases)
    )
    problem = equilibrium.GeneralReactiveEquilibriumProblem(
        "declared-provider-domain-rejection",
        equilibrium.Solved(
            299.0 * units.kelvin,
            (295.0 * units.kelvin, 305.0 * units.kelvin),
            (297.0 * units.kelvin,),
        ),
        equilibrium.Fixed(bubble.extent.pressure_pa * units.pascal),
        phases,
        _nonreactive_problem(parameters.component_ids, (0.5, 0.5)),
        (),
        (
            equilibrium.EquilibriumOutput(
                "temperature", "system.temperature", "kelvin", "absolute-temperature"
            ),
        ),
        continuation_identity="declared-provider-domain-rejection",
        feed=equilibrium.direct_amounts(parameters.component_ids, (0.5, 0.5)),
    )

    result = equilibrium.solve(model, problem)
    assert result.status == "non_evaluable"
    assert result.failure is not None
    assert result.failure.code == "provider_domain_rejection"
    assert result.provider_domain_status == "failed"
