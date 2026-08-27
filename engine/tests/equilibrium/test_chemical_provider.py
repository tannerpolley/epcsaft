from __future__ import annotations

import functools
import hashlib
import json
import math
from collections.abc import Sequence
from pathlib import Path
from typing import cast

import pytest

import epcsaft
from epcsaft import _core, equilibrium
from epcsaft.equilibrium import _grepe
from tests.eos._raw_cases import MULTINEUTRAL_REACTIVE, MULTISOLVENT_MULTIION

from .chemical_equilibrium_cases import typed_problem as _typed_problem


def _packet_fingerprint(packet_root: Path) -> str:
    digest = hashlib.sha256()
    for path in sorted(item for item in packet_root.rglob("*") if item.is_file()):
        relative = path.relative_to(packet_root).as_posix().encode()
        digest.update(len(relative).to_bytes(8, "big"))
        digest.update(relative)
        digest.update(path.stat().st_size.to_bytes(8, "big"))
        digest.update(path.read_bytes())
    return digest.hexdigest()


def _assert_provider_matches_direct_partials(
    parameters: epcsaft.Parameters,
    identities: tuple[str, ...],
    amounts: tuple[float, ...],
) -> None:
    model = epcsaft.Mixture(parameters)
    active = epcsaft.ActiveParameterSet(parameters, identities)
    volume = math.fsum(amounts) / 2000.0
    coordinates = model._active_coordinates(active)
    direct = model._native.evaluate_phase_derivatives(298.15, amounts, volume, coordinates, 2)
    provider = _core.evaluate_parameterized_phase(
        model._native_model(active),
        298.15,
        amounts,
        volume,
        1.0e-8,
        0.74,
    )

    state_count = len(amounts) + 1
    parameter_count = len(identities)
    joint_count = state_count + parameter_count
    direct_gradient = tuple(direct.gradient)
    direct_hessian = tuple(direct.hessian)
    state_hessian = tuple(
        direct_hessian[row * joint_count + column]
        for row in range(state_count)
        for column in range(state_count)
    )
    mixed = tuple(
        direct_hessian[row * joint_count + state_count + parameter]
        for row in range(state_count)
        for parameter in range(parameter_count)
    )

    def provider_values(key: str) -> tuple[float, ...]:
        return tuple(float(value) for value in cast(Sequence[object], provider[key]))

    assert tuple(cast(Sequence[object], provider["active_parameter_identities"])) == identities
    assert float(provider["value"]) == pytest.approx(
        direct.helmholtz_over_rt_reference_amount, rel=2.0e-12, abs=2.0e-12
    )
    assert provider_values("state_gradient") == pytest.approx(
        direct_gradient[:state_count], rel=2.0e-12, abs=2.0e-12
    )
    assert provider_values("state_hessian") == pytest.approx(
        state_hessian, rel=2.0e-12, abs=2.0e-12
    )
    assert provider_values("parameter_derivatives") == pytest.approx(
        direct_gradient[state_count:], rel=2.0e-12, abs=2.0e-12
    )
    assert provider_values("state_parameter_derivatives") == pytest.approx(
        mixed, rel=2.0e-12, abs=2.0e-12
    )
    assert provider_values("pressure_parameter_derivatives_pa") == pytest.approx(
        tuple(
            -8.31446261815324 * 298.15 * mixed[-parameter_count + parameter]
            for parameter in range(parameter_count)
        ),
        rel=2.0e-12,
        abs=2.0e-8,
    )
    assert provider_values("chemical_potential_parameter_derivatives_over_rt") == pytest.approx(
        mixed[: len(amounts) * parameter_count], rel=2.0e-12, abs=2.0e-12
    )


def test_provider_zero_coordinate_state_hessian_matches_value_owner() -> None:
    """Falsify a dummy-coordinate requirement in exact value/state derivatives."""

    model = epcsaft.Mixture(epcsaft.Parameters.from_mapping(MULTINEUTRAL_REACTIVE))
    amounts = (0.70, 0.20, 0.05, 0.05)
    volume = math.fsum(amounts) / 2000.0
    provider = _core.evaluate_parameterized_phase(
        model._native_model(), 313.15, amounts, volume, 1.0e-8, 0.74
    )
    state = model._state(
        T=313.15 * epcsaft.unit_registry.kelvin,
        rho=2000.0 * epcsaft.unit_registry.mole / epcsaft.unit_registry.meter**3,
        x=amounts,
    )
    gradient = tuple(float(value) for value in cast(Sequence[object], provider["state_gradient"]))
    hessian = tuple(float(value) for value in cast(Sequence[object], provider["state_hessian"]))
    expected_mu = tuple(
        math.log(amount / volume) + residual
        for amount, residual in zip(amounts, state.residual_chemical_potential_over_rt, strict=True)
    )
    assert tuple(cast(Sequence[object], provider["active_parameter_identities"])) == ()
    assert gradient[: len(amounts)] == pytest.approx(expected_mu, rel=2.0e-12, abs=2.0e-12)
    assert float(provider["pressure_pa"]) == pytest.approx(
        float(state.pressure.to("pascal").magnitude), rel=2.0e-12, abs=2.0e-8
    )
    assert len(hessian) == (len(amounts) + 1) ** 2
    assert all(math.isfinite(value) for value in hessian)


def test_provider_preserves_reordered_mixed_active_parameter_partials() -> None:
    """Falsify identity loss or family filtering at the EOS/Equilibrium seam."""

    _assert_provider_matches_direct_partials(
        epcsaft.Parameters.from_mapping(MULTISOLVENT_MULTIION),
        tuple(
            reversed(
                (
                    "component/solvent/segment_count",
                    "component/iodide-anion/dispersion_energy_over_k",
                    "component/solvent/segment_diameter",
                    "association/solvent/a/solvent/b/energy_over_k",
                    "pair/iodide-anion/solvent/k_ij",
                    "component/solvent/relative_permittivity",
                    "component/divalent-cation/born_diameter",
                    "component/cosolvent/solvation_factor",
                    "model/relative_permittivity/ionic_region_relative_permittivity",
                )
            )
        ),
        (0.60, 0.33, 0.02, 0.02, 0.01, 0.02),
    )


def test_provider_preserves_typed_correlation_coordinate_partials() -> None:
    """Falsify typed-correlation filtering at the EOS/Equilibrium seam."""

    _assert_provider_matches_direct_partials(
        epcsaft.Parameters.from_mapping(MULTINEUTRAL_REACTIVE),
        (
            "pair/anion/solvent/k_ij",
            "component/solvent/relative_permittivity/constant/constant",
        ),
        (0.70, 0.10, 0.10, 0.10),
    )


@pytest.mark.slow
def test_source_reference_uses_declared_multisolvent_multiion_basis() -> None:
    """Falsify equal-solvent and fixed-ion-position reference conventions."""

    component_ids = (
        "iodide-anion",
        "cosolvent",
        "divalent-cation",
        "solvent",
        "anion",
        "cation",
    )
    model = epcsaft.Mixture(
        epcsaft.Parameters.from_mapping(MULTISOLVENT_MULTIION, components=component_ids)
    )
    transfer = epcsaft.source_reference_transfer(
        model,
        epcsaft.SourceReferenceDeclaration(
            reference_state_id="multisolvent-multiion-fixture",
            component_ids=component_ids,
            solvent_mole_fractions=(0.0, 0.3, 0.0, 0.7, 0.0, 0.0),
            ion_pairs=(
                epcsaft.SourceReferenceIonPair(("cation", "anion"), (1, 1)),
                epcsaft.SourceReferenceIonPair(("cation", "iodide-anion"), (1, 1)),
                epcsaft.SourceReferenceIonPair(("divalent-cation", "anion"), (1, 2)),
            ),
            phase="liquid",
            reference_convention="molality-infinite-dilution",
            activity_convention_id="molality-infinite-dilution-v1",
            standard_molality_mol_per_kg=1.0,
            reference_pressure_pa=100_000.0,
            required_derivatives=(),
        ),
        T=298.15 * epcsaft.unit_registry.kelvin,
        P=100_000.0 * epcsaft.unit_registry.pascal,
    )

    assert transfer.source_composition == (0.0, 0.3, 0.0, 0.7, 0.0, 0.0)
    assert transfer.neutral_basis == (
        (0.0, 1.0, 0.0, 0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0, 1.0, 0.0, 0.0),
        (0.0, 0.0, 0.0, 0.0, 1.0, 1.0),
        (1.0, 0.0, 0.0, 0.0, 0.0, 1.0),
        (0.0, 0.0, 1.0, 0.0, 2.0, 0.0),
    )
    assert all(math.isfinite(value) for value in transfer.transfer_log_contractions)


@pytest.mark.slow
def test_source_reference_supports_mixed_pure_and_infinite_dilution_neutrals() -> None:
    """Falsify replacing a declared pure-neutral reference by the solvent blend."""

    component_ids = ("water", "hydronium-cation", "hydroxide-anion")
    model = epcsaft.Mixture(_held_parameters(component_ids))
    transfer = epcsaft.source_reference_transfer(
        model,
        epcsaft.SourceReferenceDeclaration(
            reference_state_id="mixed-pure-and-infinite-dilution-fixture",
            component_ids=component_ids,
            solvent_mole_fractions=(1.0, 0.0, 0.0),
            ion_pairs=(
                epcsaft.SourceReferenceIonPair(("hydronium-cation", "hydroxide-anion"), (1, 1)),
            ),
            phase="liquid",
            reference_convention="mole-fraction-mixed-reference",
            activity_convention_id="mole-fraction-mixed-reference-v1",
            standard_molality_mol_per_kg=1.0,
            reference_pressure_pa=100_000.0,
            required_derivatives=(),
            pure_components=(epcsaft.SourceReferencePureComponent("water", "liquid"),),
        ),
        T=298.15 * epcsaft.unit_registry.kelvin,
        P=100_000.0 * epcsaft.unit_registry.pascal,
    )
    pure_solvent = model.state(
        T=298.15 * epcsaft.unit_registry.kelvin,
        P=100_000.0 * epcsaft.unit_registry.pascal,
        x=(1.0, 0.0, 0.0),
        phase="liquid",
    )

    solvent_row = transfer.neutral_basis.index((1.0, 0.0, 0.0))
    assert transfer.log_fugacity_contractions[solvent_row] == pytest.approx(
        pure_solvent.fugacity.ln_coefficient[0], abs=2.0e-12
    )
    assert transfer.activity_scale_log_contractions == pytest.approx(
        (0.0,) * len(transfer.activity_scale_log_contractions), abs=0.0
    )


@functools.cache
def _held_parameter_packet() -> Path:
    fixture = Path(__file__).parents[1] / "fixtures" / "data" / "held-cameretti-sadowski-2008"
    materialization = json.loads((fixture / "materialization-v2.json").read_text())
    packet = fixture / "2"
    if _packet_fingerprint(packet) != materialization["packet_fingerprint"]:
        raise AssertionError("HELD parameter packet bytes changed")
    return packet


def _held_parameters(components: tuple[str, ...]) -> epcsaft.Parameters:
    return epcsaft.Parameters.from_json(
        _held_parameter_packet() / "parameters" / "parameters.json",
        components=components,
    )


def _water_ionization_source() -> dict[str, object]:
    path = Path(__file__).parents[2] / "data/reference/water_self_ionization_iapws_r11_07_2019.json"
    return json.loads(path.read_text())


def _iapws_p_kw(record: dict[str, object], temperature_k: float, density: float) -> float:
    correlation = record["correlation"]
    standard_state = record["standard_state"]
    assert isinstance(correlation, dict)
    assert isinstance(standard_state, dict)
    alpha = correlation["alpha"]
    beta = correlation["beta"]
    ideal = correlation["ideal_gas_p_kw_coefficients"]
    assert isinstance(alpha, list)
    assert isinstance(beta, list)
    assert isinstance(ideal, list)
    q_value = density * math.exp(
        alpha[0] + alpha[1] / temperature_k + alpha[2] * density ** (2.0 / 3.0) / temperature_k**2
    )
    p_kw_ideal = (
        ideal[0]
        + ideal[1] / temperature_k
        + ideal[2] / temperature_k**2
        + ideal[3] / temperature_k**3
    )
    return (
        -2.0
        * correlation["ion_coordination_number"]
        * (
            math.log10(1.0 + q_value)
            - q_value
            / (q_value + 1.0)
            * density
            * (beta[0] + beta[1] / temperature_k + beta[2] * density)
        )
        + p_kw_ideal
        + 2.0
        * math.log10(
            standard_state["standard_molality_mol_per_kg"]
            * standard_state["solvent_molar_mass_kg_per_mol"]
        )
    )


def _multineutral_problem() -> equilibrium.ChemicalEquilibriumProblem:
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
                ln_value=-269.6037501925047,
                source_id="synthetic-multineutral-source",
                reference_id="provider-helmholtz-coordinate-basis",
                reaction_orientation="products_positive",
                conversion_id="already-provider-basis",
                dimensionless=True,
            ),
        ),
        strict_interior_amount_floor_mol=1.0e-30,
    )


@pytest.mark.slow
def test_provider_certifies_multineutral_reactive_minimum() -> None:
    model = epcsaft.Mixture(epcsaft.Parameters.from_mapping(MULTINEUTRAL_REACTIVE))
    problem = _multineutral_problem()
    phase = equilibrium.ReactivePhase(
        "liquid",
        "liquid",
        "finite",
        equilibrium.AllComponents(),
        equilibrium.ProviderModel((1.0e-6, 0.74), "installed-provider-eos"),
    )
    result = equilibrium.solve(
        model,
        equilibrium.GeneralReactiveEquilibriumProblem(
            "multineutral-reactive-minimum",
            equilibrium.Fixed(313.15 * epcsaft.unit_registry.kelvin),
            equilibrium.Fixed(119_199.78746707023 * epcsaft.unit_registry.pascal),
            (phase,),
            problem,
            (phase.identity,),
            tuple(
                equilibrium.EquilibriumOutput(
                    f"amount-{index}",
                    "phase.amount",
                    "mole",
                    "true-species-amount",
                    phase.identity,
                    tuple(float(component == index) for component in range(4)),
                )
                for index in range(4)
            ),
        ),
    )

    assert result.status == "evaluated", result.failure
    amounts = tuple(float(value) for value in result.values)
    assert min(amounts) > problem.strict_interior_amount_floor_mol
    assert amounts[1] == pytest.approx(0.2, abs=2.0e-10)
    assert amounts[2] == pytest.approx(amounts[3], abs=2.0e-10)
    evidence = dict(result.evidence)
    assert evidence["reaction_affinity_inf_norm"] <= 2.0e-7
    assert evidence["kkt_stationarity_inf_norm"] <= 2.0e-7


def _held_water_ionization_problem() -> tuple[
    epcsaft.Mixture, dict[str, object], dict[str, object]
]:
    source = _water_ionization_source()
    state = source["state"]
    standard_state = source["standard_state"]
    values = source["values"]
    assert isinstance(state, dict)
    assert isinstance(standard_state, dict)
    assert isinstance(values, dict)
    components = ("water", "hydronium-cation", "hydroxide-anion")
    parameters = _held_parameters(components)
    model = epcsaft.Mixture(parameters)
    # The installed EOS evaluates the declared pure-solvent molecular mass.
    conversion = math.log(standard_state["standard_molality_mol_per_kg"] * 0.01801528)
    spec: dict[str, object] = {
        "species_ids": components,
        "charges": (0, 1, -1),
        "provider_fingerprint": model.parameter_fingerprint,
        "molar_masses_kg_per_mol": (0.01801528, 0.01902322, 0.01700734),
        "balance_matrix": ((2.0, 3.0, 1.0), (1.0, 1.0, 1.0)),
        "reaction_matrix": ((-2.0, 1.0, 1.0),),
        "feed_amounts": (1.0, 0.0, 0.0),
        "ln_k": (values["ln_kw"],),
        "equilibrium_constant_records": (
            {
                "source_id": source["source"]["id"],
                "reference_id": standard_state["id"],
                "reaction_orientation": "products_positive",
                "conversion_id": "source-standard-state-to-provider-neutral-reference",
                "dimensionless": True,
                "temperature_k": state["temperature_k"],
                "pressure_pa": state["pressure_pa"],
            },
        ),
        "temperature_k": state["temperature_k"],
        "pressure_pa": state["pressure_pa"],
    }
    return (
        model,
        spec,
        {
            "id": standard_state["id"],
            "activity_scale_id": standard_state["activity_scale_id"],
            "log_activity_scale_factors": (0.0, conversion, conversion),
            "reference_pressure_pa": state["pressure_pa"],
            "source_reference_component_ids": components,
            "source_reference_solvent_composition": (1.0, 0.0, 0.0),
            "source_reference_ion_pairs": (
                epcsaft.SourceReferenceIonPair(("hydronium-cation", "hydroxide-anion"), (1, 1)),
            ),
            "source_reference_phase": "liquid",
            "source_reference_convention": "molality-infinite-dilution",
            "source_reference_activity_convention_id": ("molality-infinite-dilution-v1"),
            "source_reference_standard_molality_mol_per_kg": 1.0,
        },
    )


def _typed_standard_state(
    source_standard_state: dict[str, object],
) -> equilibrium.ChemicalStandardState:
    return equilibrium.ChemicalStandardState(
        id=source_standard_state["id"],
        activity_scale_id=source_standard_state["activity_scale_id"],
        log_activity_scale_factors=tuple(source_standard_state["log_activity_scale_factors"]),
        reference_pressure_pa=source_standard_state["reference_pressure_pa"],
        source_reference_component_ids=tuple(
            source_standard_state["source_reference_component_ids"]
        ),
        source_reference_solvent_composition=tuple(
            source_standard_state["source_reference_solvent_composition"]
        ),
        source_reference_ion_pairs=tuple(source_standard_state["source_reference_ion_pairs"]),
        source_reference_phase=source_standard_state["source_reference_phase"],
        source_reference_convention=source_standard_state["source_reference_convention"],
        source_reference_activity_convention_id=source_standard_state[
            "source_reference_activity_convention_id"
        ],
        source_reference_standard_molality_mol_per_kg=(
            source_standard_state["source_reference_standard_molality_mol_per_kg"]
        ),
    )


def test_iapws_2019_source_record_reproduces_release_equations() -> None:
    source = _water_ionization_source()
    state = source["state"]
    values = source["values"]
    assert isinstance(state, dict)
    assert isinstance(values, dict)
    assert _iapws_p_kw(source, 300.0, 1.0) == pytest.approx(13.906565, abs=5.0e-7)
    p_kw = _iapws_p_kw(
        source,
        state["temperature_k"],
        state["density_kg_per_m3"] / 1000.0,
    )
    assert p_kw == pytest.approx(values["p_kw"], abs=5.0e-13)
    assert -math.log(10.0) * p_kw == pytest.approx(values["ln_kw"], abs=1.0e-13)


@pytest.mark.slow
def test_held_water_self_ionization_consumes_source_reference_and_provider() -> None:
    model, spec, source_standard_state = _held_water_ionization_problem()
    reaction = _typed_problem(
        spec,
        source_standard_state=_typed_standard_state(source_standard_state),
    )
    start_amounts = (1.0 - 2.0e-8, 1.0e-8, 1.0e-8)
    start_state = model.state(
        T=spec["temperature_k"] * epcsaft.unit_registry.kelvin,
        P=spec["pressure_pa"] * epcsaft.unit_registry.pascal,
        x=start_amounts,
        phase="liquid",
    )
    phase = equilibrium.ReactivePhase(
        "liquid",
        "liquid",
        "finite",
        equilibrium.AllComponents(),
        equilibrium.ProviderModel((1.0e-6, 0.74), "installed-provider-eos"),
        start=equilibrium.FinitePhaseStart(
            start_amounts,
            1.0 / float(start_state.molar_density.to("mole / meter**3").magnitude),
        ),
    )
    result = equilibrium.solve(
        model,
        equilibrium.GeneralReactiveEquilibriumProblem(
            "held-water-self-ionization",
            equilibrium.Fixed(spec["temperature_k"] * epcsaft.unit_registry.kelvin),
            equilibrium.Fixed(spec["pressure_pa"] * epcsaft.unit_registry.pascal),
            (phase,),
            reaction,
            (phase.identity,),
            tuple(
                equilibrium.EquilibriumOutput(
                    f"amount-{index}",
                    "phase.amount",
                    "mole",
                    "true-species-amount",
                    phase.identity,
                    tuple(float(component == index) for component in range(3)),
                )
                for index in range(3)
            ),
        ),
    )

    assert result.status == "evaluated", result.failure
    amounts = tuple(float(value) for value in result.values)
    assert 2.0 * amounts[0] + 3.0 * amounts[1] + amounts[2] == pytest.approx(2.0, abs=2.0e-9)
    assert sum(amounts) == pytest.approx(1.0, abs=2.0e-9)
    assert amounts[1] == pytest.approx(amounts[2], abs=2.0e-12)
    evidence = dict(result.evidence)
    assert evidence["source_reference_representation_residual_inf_norm"] <= 1.0e-12
    assert evidence["reaction_affinity_inf_norm"] <= 2.0e-12
    assert evidence["kkt_stationarity_inf_norm"] <= 2.0e-12


@pytest.mark.slow
def test_source_reference_pressure_column_matches_centered_resolves() -> None:
    """Keep the source-pressure transform inside the native total derivative."""

    temperature_k = 313.15
    pressure_pa = 119_199.78746707023
    model = epcsaft.Mixture(epcsaft.Parameters.from_mapping(MULTINEUTRAL_REACTIVE))
    base = _multineutral_problem()
    declaration = epcsaft.SourceReferenceDeclaration(
        reference_state_id="mixed-source-pressure-density-fixture",
        component_ids=base.species_ids,
        solvent_mole_fractions=(0.8, 0.2, 0.0, 0.0),
        ion_pairs=(epcsaft.SourceReferenceIonPair(("cation", "anion"), (1, 1)),),
        phase="liquid",
        reference_convention="mole-fraction-mixed-reference",
        activity_convention_id="mole-fraction-mixed-reference-v1",
        standard_molality_mol_per_kg=1.0,
        reference_pressure_pa=pressure_pa,
        required_derivatives=(),
    )
    transfer = epcsaft.source_reference_transfer(
        model,
        declaration,
        T=temperature_k * epcsaft.unit_registry.kelvin,
        P=pressure_pa * epcsaft.unit_registry.pascal,
    )
    solvent_row = transfer.neutral_basis.index((1.0, 0.0, 0.0, 0.0))
    ion_pair_row = transfer.neutral_basis.index((0.0, 0.0, 1.0, 1.0))
    offset = (
        -transfer.transfer_log_contractions[solvent_row]
        + transfer.transfer_log_contractions[ion_pair_row]
        + math.log(pressure_pa / (8.31446261815324 * temperature_k))
    )
    standard = equilibrium.ChemicalStandardState(
        id=declaration.reference_state_id,
        activity_scale_id=declaration.activity_convention_id,
        log_activity_scale_factors=(0.0, 0.0, 0.0, 0.0),
        reference_pressure_pa=pressure_pa,
        source_reference_component_ids=base.species_ids,
        source_reference_solvent_composition=(0.8, 0.2, 0.0, 0.0),
        source_reference_ion_pairs=declaration.ion_pairs,
        source_reference_phase="liquid",
        source_reference_convention="mole-fraction-mixed-reference",
        source_reference_activity_convention_id=declaration.activity_convention_id,
        source_reference_standard_molality_mol_per_kg=1.0,
    )
    reaction = equilibrium.ChemicalEquilibriumProblem(
        base.species_ids,
        base.charges,
        base.molar_masses_kg_per_mol,
        base.balance_matrix,
        base.conserved_totals,
        base.reaction_matrix,
        base.feed_amounts_mol,
        (
            equilibrium.ChemicalEquilibriumConstant(
                base.equilibrium_constants[0].ln_value - offset,
                "synthetic-multineutral-source",
                declaration.reference_state_id,
                "products_positive",
                "source-standard-state-to-provider-neutral-reference",
                True,
            ),
        ),
        base.strict_interior_amount_floor_mol,
        standard,
    )
    phase = equilibrium.ReactivePhase(
        "liquid",
        "liquid",
        "finite",
        equilibrium.AllComponents(),
        equilibrium.ProviderModel((1.0e-6, 0.74), "installed-provider-eos"),
    )
    problem = equilibrium.GeneralReactiveEquilibriumProblem(
        "multineutral-source-pressure-column",
        equilibrium.Fixed(temperature_k * epcsaft.unit_registry.kelvin),
        equilibrium.Fixed(pressure_pa * epcsaft.unit_registry.pascal),
        (phase,),
        reaction,
        (phase.identity,),
        (
            equilibrium.EquilibriumOutput(
                "solvent-amount",
                "phase.amount",
                "mole",
                "true-species-amount",
                phase.identity,
                (1.0, 0.0, 0.0, 0.0),
            ),
        ),
    )
    total = math.fsum(reaction.feed_amounts_mol)
    start_state = model.state(
        T=temperature_k * epcsaft.unit_registry.kelvin,
        P=pressure_pa * epcsaft.unit_registry.pascal,
        x=tuple(amount / total for amount in reaction.feed_amounts_mol),
        phase="liquid",
    )
    starts = (
        {
            "amounts_mol": reaction.feed_amounts_mol,
            "molar_volume_m3_per_mol": 1.0
            / float(start_state.molar_density.to("mole / meter**3").magnitude),
        },
    )
    outputs = tuple(
        {
            "selector": "phase.amount",
            "phase_index": 0,
            "coefficients": tuple(float(component == species) for component in range(4)),
            "solvent_mass_coefficients_kg_per_mol": (),
        }
        for species in range(4)
    )

    def solve_native(trial_pressure_pa: float) -> dict[str, object]:
        source = _grepe._reactive_finite_source_standard_state(
            model,
            problem,
            phase,
            temperature_k,
            trial_pressure_pa,
            None,
        )
        return _core._reactive_finite_phases(
            model._native_model(),
            _grepe._finite_phase_spec(problem, temperature_k, trial_pressure_pa),
            source,
            ((1.0e-6, 0.74),),
            starts,
            reaction.strict_interior_amount_floor_mol,
            False,
            False,
            None,
            ((0, 1, 2, 3),),
            outputs,
        )

    center = solve_native(pressure_pa)
    lower = solve_native(pressure_pa - 5.0)
    upper = solve_native(pressure_pa + 5.0)
    assert center["accepted"] and lower["accepted"] and upper["accepted"]
    pressure_column = center["input_order"].index("pressure_pa")
    width = len(center["input_order"])
    for species in range(4):
        analytic = center["output_input_derivatives"][species * width + pressure_column]
        numerical = (upper["output_values"][species] - lower["output_values"][species]) / 10.0
        assert abs(analytic - numerical) <= 5.0e-5 * max(abs(analytic), abs(numerical)) + 2.0e-9
