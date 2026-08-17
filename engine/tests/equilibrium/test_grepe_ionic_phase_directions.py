from __future__ import annotations

import math
from dataclasses import replace
from pathlib import Path

import pytest

import epcsaft
from epcsaft import equilibrium


def test_electrolyte_finite_phases_preserve_exact_support_and_declared_directions() -> None:
    """Configured support must keep excluded ions exactly absent without losing charge."""

    component_ids = ("water", "sodium-cation", "chloride-anion")
    parameter_path = (
        Path(__file__).parents[1]
        / "fixtures/data/figiel-2025-reference-electrolytes/2/parameters/parameters.json"
    )
    parameters = epcsaft.Parameters.from_json(parameter_path, components=component_ids)
    model = epcsaft.Mixture(parameters)
    active = epcsaft.ActiveParameterSet(parameters, ("component/water/segment_count",))
    temperature_k = 298.15
    pressure_pa = 3_168.0
    starts = ((50.0, 1.0, 1.0), (5.5, 0.0, 0.0))
    phases = []
    start_states = []
    for identity, role, support, amounts in zip(
        ("liquid", "vapor"),
        ("liquid", "vapor"),
        (
            equilibrium.DeclaredComponents(tuple(reversed(component_ids))),
            equilibrium.DeclaredComponents(("water",)),
        ),
        starts,
        strict=True,
    ):
        total = math.fsum(amounts)
        state = model.state(
            T=temperature_k * epcsaft.unit_registry.kelvin,
            P=pressure_pa * epcsaft.unit_registry.pascal,
            x=tuple(amount / total for amount in amounts),
            phase=role,
        )
        start_states.append(state)
        molar_density = float(state.molar_density.to("mole / meter**3").magnitude)
        phases.append(
            equilibrium.ReactivePhase(
                identity,
                role,
                "finite",
                support,
                equilibrium.ProviderModel((1.0e-10, 0.74), "installed-provider"),
                start=equilibrium.FinitePhaseStart(amounts, 1.0 / molar_density),
            )
        )
    feed = tuple(math.fsum(row[index] for row in starts) for index in range(3))
    balances = ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))
    reaction_system = equilibrium.ChemicalEquilibriumProblem(
        component_ids,
        (0, 1, -1),
        (0.01801528, 0.02298977, 0.03545),
        balances,
        feed,
        (),
        feed,
        (),
        1.0e-14,
    )
    phase_reactions = equilibrium.PhaseReactionSystem(
        ("liquid", "vapor"),
        component_ids,
        (
            equilibrium.PhaseReactionRow(
                "water-transfer",
                (
                    equilibrium.PhaseSpeciesTerm("liquid", "water", -1.0),
                    equilibrium.PhaseSpeciesTerm("vapor", "water", 1.0),
                ),
                start_states[1].fugacity.ln_coefficient[0]
                - (
                    math.log(starts[0][0] / math.fsum(starts[0]))
                    + start_states[0].fugacity.ln_coefficient[0]
                ),
            ),
        ),
        (0.0, 1.0, -1.0),
        balances,
    )
    outputs = tuple(
        equilibrium.EquilibriumOutput(
            f"{phase_identity}-{component_id}-amount",
            "phase.amount",
            "mole",
            "manufactured-electrolyte-fixed-topology",
            phase_identity,
            tuple(float(index == component) for index in range(3)),
            support="nonnegative",
        )
        for phase_identity in ("liquid", "vapor")
        for component, component_id in enumerate(component_ids)
    )
    problem = equilibrium.GeneralReactiveEquilibriumProblem(
        "manufactured-electrolyte-finite-phases",
        equilibrium.Fixed(temperature_k * epcsaft.unit_registry.kelvin),
        equilibrium.Fixed(pressure_pa * epcsaft.unit_registry.pascal),
        tuple(phases),
        reaction_system,
        (),
        outputs,
        phase_reactions=phase_reactions,
    )
    assert equilibrium.general_reactive_equilibrium_problem_from_mapping(
        equilibrium.general_reactive_equilibrium_problem_to_mapping(problem)
    ) == problem

    result = equilibrium.phase_equilibrium(model, problem, active_parameters=active)

    assert result.status == "evaluated", result.failure
    assert result.owner_result is not None
    assert result.owner_result.declared_direction_rank == 1
    assert result.owner_result.admissible_nullity == 1
    assert result.owner_result.reaction_affinity_inf_norm < 1.0e-7
    assert result.owner_result.balance_inf_norm < 1.0e-9
    assert result.rows[4].value == 0.0
    assert result.rows[5].value == 0.0
    assert result.rows[4].jacobian == (0.0,)
    assert result.rows[5].jacobian == (0.0,)
    assert result.phases[0].mole_fractions[1] == pytest.approx(
        result.phases[0].mole_fractions[2], abs=2.0e-12
    )
    assert result.phases[1].mole_fractions[1:] == (0.0, 0.0)
    assert result.phases[0].supported_component_ids == component_ids
    assert result.phases[1].supported_component_ids == ("water",)
    assert result.phases[0].molar_density_mol_m3 > 10_000.0
    assert result.phases[1].molar_density_mol_m3 < 10.0

    rejected = equilibrium.phase_equilibrium(
        model,
        replace(
            problem,
            phase_reactions=replace(
                phase_reactions,
                reactions=(replace(phase_reactions.reactions[0], ln_equilibrium_constant=0.0),),
            ),
        ),
        active_parameters=active,
    )
    assert rejected.status == "non_evaluable"
    assert len(rejected.phases) == 2
    assert rejected.phases[1].mole_fractions[1:] == (0.0, 0.0)
    assert rejected.phases[1].mechanical_class == "rejected_configured_candidate"
