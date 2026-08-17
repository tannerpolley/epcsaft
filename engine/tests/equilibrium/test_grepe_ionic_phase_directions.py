from __future__ import annotations

import math
from pathlib import Path

import pytest

import epcsaft
from epcsaft import equilibrium


def test_electrolyte_finite_phases_preserve_phase_charge_and_declared_directions() -> None:
    """An ionic chart must remove charge invariants without losing fixed inventories."""

    component_ids = ("water", "sodium-cation", "chloride-anion")
    parameter_path = (
        Path(__file__).parents[1]
        / "fixtures/data/figiel-2025-reference-electrolytes/2/parameters/parameters.json"
    )
    model = epcsaft.Mixture(epcsaft.Parameters.from_json(parameter_path, components=component_ids))
    temperature_k = 298.15
    pressure_pa = 3_168.0
    starts = ((50.0, 0.999, 0.999), (5.5, 0.001, 0.001))
    phases = []
    for identity, role, amounts in zip(
        ("liquid", "vapor"), ("liquid", "vapor"), starts, strict=True
    ):
        total = math.fsum(amounts)
        state = model.state(
            T=temperature_k * epcsaft.unit_registry.kelvin,
            P=pressure_pa * epcsaft.unit_registry.pascal,
            x=tuple(amount / total for amount in amounts),
            phase=role,
        )
        molar_density = float(state.molar_density.to("mole / meter**3").magnitude)
        phases.append(
            equilibrium.ReactivePhase(
                identity,
                role,
                "finite",
                equilibrium.AllComponents(),
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
            support="positive",
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

    result = equilibrium.phase_equilibrium(model, problem)

    assert result.status == "evaluated", result.failure
    assert result.owner_result is not None
    assert result.owner_result.declared_direction_rank == 1
    assert result.owner_result.admissible_nullity == 1
    assert result.owner_result.reaction_affinity_inf_norm < 1.0e-7
    assert result.owner_result.balance_inf_norm < 1.0e-9
    assert result.rows[4].value == pytest.approx(starts[1][1], abs=2.0e-11)
    assert result.rows[5].value == pytest.approx(starts[1][2], abs=2.0e-11)
    for phase in result.phases:
        assert phase.mole_fractions[1] == pytest.approx(phase.mole_fractions[2], abs=2.0e-12)
    assert result.phases[0].molar_density_mol_m3 > 10_000.0
    assert result.phases[1].molar_density_mol_m3 < 10.0
