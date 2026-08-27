from __future__ import annotations

import math
from copy import deepcopy
from dataclasses import replace
from pathlib import Path

import pytest

import epcsaft
from epcsaft import equilibrium
from tests.eos._raw_cases import (
    HELD_STRATEGY_2_TERNARY,
    _component,
    _constant_correlation,
    _pair,
)


def _multineutral_electrolyte(*, organic_permittivity: float | None) -> dict[str, object]:
    mapping = deepcopy(HELD_STRATEGY_2_TERNARY)
    mapping["components"].append(
        _component(
            "neutral-organic",
            molar_mass=0.100,
            charge_number=0,
            segment_count=2.0,
            segment_diameter=3.5,
            dispersion_energy_over_k=250.0,
        )
    )
    for other in ("water", "lithium-cation", "chloride-anion"):
        mapping["pairs"].append(_pair(other, "neutral-organic", k_ij=0.0))
    if organic_permittivity is not None:
        mapping["correlations"].append(
            _constant_correlation("neutral-organic", "relative_permittivity", organic_permittivity)
        )
    return mapping


def test_dielectric_solvent_identity_is_scoped_to_charged_phase_support() -> None:
    """A neutral cosolute is not globally forced to be the electrolyte solvent."""

    temperature = 298.15 * epcsaft.unit_registry.kelvin
    parameters = epcsaft.Parameters.from_mapping(
        _multineutral_electrolyte(organic_permittivity=None)
    )
    model = epcsaft.Mixture(parameters)
    component_ids = ("water", "lithium-cation", "chloride-anion", "neutral-organic")
    parameter_identity = "component/water/segment_count"
    selected = epcsaft.ActiveParameterSet(parameters, (parameter_identity,))
    pressure = 1.0e5 * epcsaft.unit_registry.pascal
    starts = ((50.0, 1.0, 1.0, 0.5), (0.0, 0.0, 0.0, 5.0))
    supports = (
        equilibrium.DeclaredComponents(component_ids),
        equilibrium.DeclaredComponents(("neutral-organic",)),
    )
    phases = []
    states = []
    for identity, support, amounts in zip(
        ("aqueous", "organic"), supports, starts, strict=True
    ):
        total = math.fsum(amounts)
        state = model.state(
            T=temperature,
            P=pressure,
            x=tuple(amount / total for amount in amounts),
            phase="liquid",
        )
        states.append(state)
        molar_density = float(state.molar_density.to("mole / meter**3").magnitude)
        phases.append(
            equilibrium.ReactivePhase(
                identity,
                "liquid",
                "finite",
                support,
                equilibrium.ProviderModel((1.0e-10, 0.74), "installed-provider"),
                start=equilibrium.FinitePhaseStart(amounts, 1.0 / molar_density),
            )
        )
    assert states[0].bulk_relative_permittivity == pytest.approx(78.09)
    assert states[0].debye_huckel < 0.0
    assert states[1].debye_huckel == 0.0
    feed = tuple(math.fsum(row[index] for row in starts) for index in range(4))
    balances = tuple(tuple(float(row == column) for column in range(4)) for row in range(4))
    reaction_system = equilibrium.ChemicalEquilibriumProblem(
        component_ids,
        (0, 1, -1, 0),
        (0.01801528, 0.00694, 0.03545, 0.100),
        balances,
        feed,
        (),
        feed,
        (),
        1.0e-14,
    )
    transfer_ln_k = (
        math.log(starts[1][3] / math.fsum(starts[1]))
        + states[1].fugacity.ln_coefficient[3]
        - math.log(starts[0][3] / math.fsum(starts[0]))
        - states[0].fugacity.ln_coefficient[3]
    )
    phase_reactions = equilibrium.PhaseReactionSystem(
        ("aqueous", "organic"),
        component_ids,
        (
            equilibrium.PhaseReactionRow(
                "cosolute-transfer",
                (
                    equilibrium.PhaseSpeciesTerm(
                        "aqueous", "neutral-organic", -1.0
                    ),
                    equilibrium.PhaseSpeciesTerm(
                        "organic", "neutral-organic", 1.0
                    ),
                ),
                transfer_ln_k,
            ),
        ),
        (0.0, 1.0, -1.0, 0.0),
        balances,
    )
    outputs = tuple(
        equilibrium.EquilibriumOutput(
            f"organic-{component_id}-amount",
            "phase.amount",
            "mole",
            "manufactured-phase-scoped-solvent",
            "organic",
            selector,
            support="nonnegative",
        )
        for component_id, selector in (
            ("lithium", (0.0, 1.0, 0.0, 0.0)),
            ("cosolute", (0.0, 0.0, 0.0, 1.0)),
        )
    )
    problem = equilibrium.GeneralReactiveEquilibriumProblem(
        "manufactured-phase-scoped-solvent",
        equilibrium.Fixed(temperature),
        equilibrium.Fixed(pressure),
        tuple(phases),
        reaction_system,
        (),
        outputs,
        phase_reactions=phase_reactions,
    )

    def solve(value: float) -> equilibrium.GeneralReactiveEquilibriumResult:
        return equilibrium.solve(
            model,
            problem,
            active_parameters=epcsaft.ActiveParameterSet(
                parameters, (parameter_identity,), values=(value,)
            ),
        )

    step = 1.0e-5
    result = solve(selected.values[0])
    below = solve(selected.values[0] - step)
    above = solve(selected.values[0] + step)

    assert result.status == "evaluated", result.failure
    evidence = dict(result.evidence)
    assert evidence["balance_inf_norm"] < 1.0e-9
    assert evidence["reaction_affinity_inf_norm"] < 1.0e-7
    assert result.phases[0].supported_component_ids == component_ids
    assert result.phases[1].supported_component_ids == ("neutral-organic",)
    assert result.rows[0].value == 0.0
    assert result.rows[1].jacobian is not None
    numerical = (above.values[1] - below.values[1]) / (2.0 * step)
    assert numerical != 0.0
    assert result.rows[1].jacobian[0] == pytest.approx(numerical, rel=3.0e-6, abs=1.0e-8)

    ambiguous_model = epcsaft.Mixture(
        epcsaft.Parameters.from_mapping(
            _multineutral_electrolyte(organic_permittivity=4.0)
        )
    )
    no_solvent_problem = replace(
        problem,
        phases=(
            replace(
                problem.phases[0],
                support=equilibrium.DeclaredComponents(component_ids[1:]),
                start=replace(
                    problem.phases[0].start,
                    amounts_mol=(0.0, 1.0, 1.0, 0.5),
                ),
            ),
            problem.phases[1],
        ),
    )
    multisolvent = equilibrium.solve(ambiguous_model, problem)
    assert multisolvent.status == "evaluated", multisolvent.failure
    with pytest.raises(epcsaft.ParameterError, match="dielectric-solvent-required"):
        equilibrium.solve(model, no_solvent_problem)


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

    result = equilibrium.solve(model, problem, active_parameters=active)

    assert result.status == "evaluated", result.failure
    evidence = dict(result.evidence)
    assert evidence["declared_direction_rank"] == 1
    assert evidence["admissible_nullity"] == 1
    assert evidence["reaction_affinity_inf_norm"] < 1.0e-7
    assert evidence["balance_inf_norm"] < 1.0e-9
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

    rejected = equilibrium.solve(
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
