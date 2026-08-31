from __future__ import annotations

import math
from dataclasses import replace

import epcsaft
from epcsaft import equilibrium
from tests.eos._cases import (
    AQUEOUS_ELECTROLYTE,
    MULTINEUTRAL_REACTIVE,
    NEUTRAL_TERNARY_IMMISCIBLE,
    STRONGLY_ASSOCIATING_BINARY,
)


def nonreactive_system(
    component_ids: tuple[str, ...], amounts: tuple[float, ...]
) -> equilibrium.ChemicalEquilibriumProblem:
    size = len(component_ids)
    return equilibrium.ChemicalEquilibriumProblem(
        species_ids=component_ids,
        charges=(0,) * size,
        molar_masses_kg_per_mol=(0.05,) * size,
        balance_matrix=tuple(
            tuple(float(row == column) for column in range(size)) for row in range(size)
        ),
        conserved_totals=amounts,
        reaction_matrix=(),
        feed_amounts_mol=amounts,
        equilibrium_constants=(),
        strict_interior_amount_floor_mol=1.0e-14,
    )


def finite_phase_case() -> tuple[
    epcsaft.Parameters,
    epcsaft.Mixture,
    equilibrium.GeneralReactiveEquilibriumProblem,
]:
    parameters = epcsaft.Parameters.from_mapping(NEUTRAL_TERNARY_IMMISCIBLE)
    model = epcsaft.Mixture(parameters)
    feed = (1.0, 1.0, 1.0)
    starts = ((0.9, 0.55, 0.05), (0.1, 0.45, 0.95))
    phases = []
    for identity, amounts in zip(("phase-a", "phase-b"), starts, strict=True):
        total = math.fsum(amounts)
        state = model.state(
            T=300.0 * epcsaft.unit_registry.kelvin,
            P=1.0e5 * epcsaft.unit_registry.pascal,
            x=tuple(amount / total for amount in amounts),
            phase="liquid",
        )
        phases.append(
            equilibrium.ReactivePhase(
                identity,
                "liquid",
                "finite",
                equilibrium.AllComponents(),
                equilibrium.ProviderModel((1.0e-8, 0.74), "installed-provider-eos"),
                equilibrium.FinitePhaseStart(
                    amounts,
                    1.0 / float(state.molar_density.to("mole / meter**3").magnitude),
                ),
            )
        )
    reaction = nonreactive_system(model.component_ids, feed)
    directions = equilibrium.PhaseReactionSystem(
        tuple(phase.identity for phase in phases),
        model.component_ids,
        tuple(
            equilibrium.PhaseReactionRow(
                f"transfer-{component_id}",
                (
                    equilibrium.PhaseSpeciesTerm("phase-a", component_id, -1.0),
                    equilibrium.PhaseSpeciesTerm("phase-b", component_id, 1.0),
                ),
            )
            for component_id in model.component_ids
        ),
        (0.0, 0.0, 0.0),
        reaction.balance_matrix,
    )
    outputs = tuple(
        equilibrium.EquilibriumOutput(
            f"{phase.identity}-amount-{component}",
            "phase.amount",
            "mole",
            "true-species-amount",
            phase.identity,
            tuple(float(component == index) for index in range(3)),
            support="positive",
        )
        for phase in phases
        for component in range(3)
    )
    return (
        parameters,
        model,
        equilibrium.GeneralReactiveEquilibriumProblem(
            "finite-two-phase",
            equilibrium.Fixed(300.0 * epcsaft.unit_registry.kelvin),
            equilibrium.Fixed(1.0e5 * epcsaft.unit_registry.pascal),
            tuple(phases),
            reaction,
            (),
            outputs,
            phase_reactions=directions,
        ),
    )


def ionic_reaction_case(
    component_ids: tuple[str, ...] | None = None,
) -> tuple[epcsaft.Mixture, equilibrium.GeneralReactiveEquilibriumProblem]:
    parameters = epcsaft.Parameters.from_mapping(AQUEOUS_ELECTROLYTE)
    order = parameters.component_ids if component_ids is None else component_ids
    parameters = parameters.select(order)
    model = epcsaft.Mixture(parameters)
    charges = {"water": 0, "sodium-cation": 1, "chloride-anion": -1}
    feed = {"water": 0.0, "sodium-cation": 1.0, "chloride-anion": 1.0}
    stoichiometry = {"water": 2.0, "sodium-cation": -1.0, "chloride-anion": -1.0}
    amounts = tuple(feed[identity] for identity in order)
    charge_row = tuple(float(charges[identity]) for identity in order)
    reaction = equilibrium.ChemicalEquilibriumProblem(
        species_ids=order,
        charges=tuple(charges[identity] for identity in order),
        molar_masses_kg_per_mol=(1.0,) * 3,
        balance_matrix=((1.0,) * 3, charge_row),
        conserved_totals=(2.0, 0.0),
        reaction_matrix=(tuple(stoichiometry[identity] for identity in order),),
        feed_amounts_mol=amounts,
        equilibrium_constants=(
            equilibrium.ChemicalEquilibriumConstant(
                math.log(4.0),
                "manufactured-reaction",
                "provider-helmholtz-coordinate-basis",
                "products_positive",
                "already-provider-basis",
                True,
            ),
        ),
        strict_interior_amount_floor_mol=1.0e-14,
    )
    phase = equilibrium.ReactivePhase(
        "vapor",
        "vapor",
        "finite",
        equilibrium.AllComponents(),
        equilibrium.IdealGasModel("provider-helmholtz-coordinate-basis"),
        equilibrium.FinitePhaseStart(
            tuple(0.2 if identity == "water" else 0.9 for identity in order), 1.0
        ),
    )
    outputs = tuple(
        equilibrium.EquilibriumOutput(
            f"amount-{identity}",
            "phase.amount",
            "mole",
            "true-species-amount",
            phase.identity,
            tuple(float(identity == candidate) for candidate in order),
            support="positive",
        )
        for identity in order
    )
    return model, equilibrium.GeneralReactiveEquilibriumProblem(
        "analytic-ionic-reaction",
        equilibrium.Fixed(350.0 * epcsaft.unit_registry.kelvin),
        equilibrium.Fixed(2.0e5 * epcsaft.unit_registry.pascal),
        (phase,),
        reaction,
        (phase.identity,),
        outputs,
    )


def multineutral_reaction() -> equilibrium.ChemicalEquilibriumProblem:
    return equilibrium.ChemicalEquilibriumProblem(
        species_ids=("solvent", "cosolvent", "cation", "anion"),
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
                -269.6037501925047,
                "manufactured-reaction",
                "provider-helmholtz-coordinate-basis",
                "products_positive",
                "already-provider-basis",
                True,
            ),
        ),
        strict_interior_amount_floor_mol=1.0e-30,
    )


def solved_pressure_case(
    *, ideal_vapor: bool = False
) -> tuple[
    epcsaft.Parameters,
    epcsaft.Mixture,
    equilibrium.GeneralReactiveEquilibriumProblem,
]:
    parameters = epcsaft.Parameters.from_mapping(MULTINEUTRAL_REACTIVE)
    model = epcsaft.Mixture(parameters)
    reaction = multineutral_reaction()
    phases = (
        equilibrium.ReactivePhase(
            "finite-liquid",
            "liquid",
            "finite",
            equilibrium.AllComponents(),
            equilibrium.ProviderModel((1.0e-6, 0.70), "installed-provider-eos"),
        ),
        equilibrium.ReactivePhase(
            "incipient-vapor",
            "vapor",
            "incipient",
            equilibrium.DeclaredNeutralComponents(("solvent", "cosolvent")),
            equilibrium.IdealGasModel("provider-helmholtz-coordinate-basis")
            if ideal_vapor
            else equilibrium.ProviderModel((1.0e-6, 0.70), "installed-provider-eos"),
        ),
    )
    outputs = (
        equilibrium.EquilibriumOutput(
            "pressure", "system.pressure", "pascal", "true-species-pressure", support="positive"
        ),
        equilibrium.EquilibriumOutput(
            "vapor-solvent-fraction",
            "phase.mole_fraction",
            "dimensionless",
            "true-species-mole-fraction",
            phases[1].identity,
            (1.0, 0.0),
            support="positive",
        ),
        equilibrium.EquilibriumOutput(
            "vapor-solvent-partial-pressure",
            "phase.partial_pressure",
            "pascal",
            "true-species-partial-pressure",
            phases[1].identity,
            (1.0, 0.0),
            support="positive",
        ),
    )
    return (
        parameters,
        model,
        equilibrium.GeneralReactiveEquilibriumProblem(
            "ideal-finite-incipient-pressure"
            if ideal_vapor
            else "finite-incipient-pressure",
            equilibrium.Fixed(313.15 * epcsaft.unit_registry.kelvin),
            equilibrium.Solved(
                3.2e6 * epcsaft.unit_registry.pascal,
                (2.5e6 * epcsaft.unit_registry.pascal, 4.0e6 * epcsaft.unit_registry.pascal),
                (2.8e6 * epcsaft.unit_registry.pascal, 3.2e6 * epcsaft.unit_registry.pascal),
            ),
            phases,
            reaction,
            (phases[0].identity,),
            outputs,
            continuation_identity="ideal-finite-incipient-boundary"
            if ideal_vapor
            else "finite-incipient-boundary",
        ),
    )


def ideal_uptake_case() -> tuple[
    epcsaft.Mixture, equilibrium.GeneralReactiveEquilibriumProblem, float
]:
    model = epcsaft.Mixture(epcsaft.Parameters.from_mapping(AQUEOUS_ELECTROLYTE))
    component_ids = model.component_ids
    pressure_pa = 2.0e5
    closed_amount = 1.5
    expected = 0.375
    target = math.log(pressure_pa * expected / (closed_amount + expected))
    reaction = nonreactive_system(component_ids, (1.0, 0.25, 0.5))
    phase = equilibrium.ReactivePhase(
        "vapor",
        "vapor",
        "finite",
        equilibrium.AllComponents(),
        equilibrium.IdealGasModel("provider-helmholtz-coordinate-basis"),
    )
    problem = equilibrium.GeneralReactiveEquilibriumProblem(
        "analytic-open-component",
        equilibrium.Fixed(350.0 * epcsaft.unit_registry.kelvin),
        equilibrium.Fixed(pressure_pa * epcsaft.unit_registry.pascal),
        (phase,),
        reaction,
        (),
        (
            equilibrium.EquilibriumOutput(
                "uptake",
                "phase.amount",
                "mole",
                "true-species-amount",
                phase.identity,
                (0.0, 1.0, 0.0),
                support="positive",
            ),
        ),
        feed=equilibrium.AffineFeed(
            component_ids,
            (1.0, 0.0, 0.5),
            ((0.0,), (1.0,), (0.0,)),
            (equilibrium.FeedCoordinate("open-component", (0.05, 1.5)),),
        ),
        intensive_boundaries=(
            equilibrium.IntensiveBoundary(
                "reservoir-potential",
                "open-component",
                phase.identity,
                "phase.log_fugacity",
                (0.0, 1.0, 0.0),
                target,
                0.2,
            ),
        ),
    )
    return model, problem, expected


def provider_boundary_case() -> tuple[
    epcsaft.Parameters,
    epcsaft.Mixture,
    equilibrium.GeneralReactiveEquilibriumProblem,
    str,
]:
    parameters = epcsaft.Parameters.from_mapping(MULTINEUTRAL_REACTIVE)
    model = epcsaft.Mixture(parameters)
    reaction = multineutral_reaction()
    phase = equilibrium.ReactivePhase(
        "liquid",
        "liquid",
        "finite",
        equilibrium.AllComponents(),
        equilibrium.ProviderModel((1.0e-6, 0.74), "installed-provider-eos"),
    )
    problem = equilibrium.GeneralReactiveEquilibriumProblem(
        "open-component-response",
        equilibrium.Fixed(313.15 * epcsaft.unit_registry.kelvin),
        equilibrium.Fixed(119_199.78746707023 * epcsaft.unit_registry.pascal),
        (phase,),
        reaction,
        (phase.identity,),
        (
            equilibrium.EquilibriumOutput(
                "cation-fraction",
                "phase.mole_fraction",
                "dimensionless",
                "true-species-mole-fraction",
                phase.identity,
                (0.0, 0.0, 1.0, 0.0),
                support="positive",
            ),
        ),
        feed=equilibrium.AffineFeed(
            reaction.species_ids,
            (0.8, 0.1, 0.1, 0.1),
            ((0.0,), (1.0,), (0.0,), (0.0,)),
            (equilibrium.FeedCoordinate("open-neutral", (0.02, 0.5)),),
        ),
        intensive_boundaries=(
            equilibrium.IntensiveBoundary(
                "imposed-potential",
                "open-neutral",
                phase.identity,
                "phase.log_fugacity",
                (0.0, 1.0, 0.0, 0.0),
                15.0,
                0.2,
            ),
        ),
    )
    return parameters, model, problem, "pair/anion/solvent/k_ij"


def temperature_boundary_case(
    reference_role: str,
) -> tuple[
    epcsaft.Parameters,
    epcsaft.Mixture,
    equilibrium.GeneralReactiveEquilibriumProblem,
    str,
]:
    parameters = epcsaft.Parameters.from_mapping(STRONGLY_ASSOCIATING_BINARY)
    model = epcsaft.Mixture(parameters)
    reference_index = 0 if reference_role == "liquid" else 1
    reference = (0.5, 0.5) if reference_index == 0 else (0.01, 0.99)
    phases = tuple(
        equilibrium.ReactivePhase(
            role,
            role,
            "finite" if index == reference_index else "incipient",
            equilibrium.AllComponents(),
            equilibrium.ProviderModel((1.0e-9, 0.74), "installed-provider-eos"),
        )
        for index, role in enumerate(("liquid", "vapor"))
    )
    initial = 300.0 if reference_index == 0 else 297.0
    problem = equilibrium.GeneralReactiveEquilibriumProblem(
        f"{reference_role}-temperature-boundary",
        equilibrium.Solved(
            initial * epcsaft.unit_registry.kelvin,
            (295.0 * epcsaft.unit_registry.kelvin, 305.0 * epcsaft.unit_registry.kelvin),
            (297.0 * epcsaft.unit_registry.kelvin, 303.0 * epcsaft.unit_registry.kelvin),
        ),
        equilibrium.Fixed(1200.0 * epcsaft.unit_registry.pascal),
        phases,
        nonreactive_system(parameters.component_ids, reference),
        (),
        (
            equilibrium.EquilibriumOutput(
                "temperature", "system.temperature", "kelvin", "absolute-temperature"
            ),
            equilibrium.EquilibriumOutput(
                "partner-fraction",
                "phase.mole_fraction",
                "dimensionless",
                "true-species-mole-fraction",
                phases[1 - reference_index].identity,
                (1.0, 0.0),
            ),
        ),
        continuation_identity=f"{reference_role}-temperature-boundary",
        feed=equilibrium.direct_amounts(parameters.component_ids, reference),
    )
    return parameters, model, problem, "pair/solute/solvent/k_ij"


def reverse_finite_phases(
    problem: equilibrium.GeneralReactiveEquilibriumProblem,
) -> equilibrium.GeneralReactiveEquilibriumProblem:
    assert problem.phase_reactions is not None
    return replace(
        problem,
        phases=tuple(reversed(problem.phases)),
        phase_reactions=replace(
            problem.phase_reactions,
            phase_ids=tuple(reversed(problem.phase_reactions.phase_ids)),
        ),
    )
