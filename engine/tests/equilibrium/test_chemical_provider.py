from __future__ import annotations

import copy
import functools
import hashlib
import json
import math
from collections.abc import Sequence
from pathlib import Path
from typing import cast

import pytest

import epcsaft
from epcsaft import equilibrium
from epcsaft.eos import _native_transport
from epcsaft.equilibrium import _equilibrium
from tests.eos._raw_cases import MULTINEUTRAL_REACTIVE, MULTISOLVENT_MULTIION

from .chemical_equilibrium_cases import bind_records as _bind_record
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
    provider = _equilibrium.evaluate_parameterized_phase(
        _native_transport(model, active),
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
    provider = _equilibrium.evaluate_parameterized_phase(
        _native_transport(model), 313.15, amounts, volume, 1.0e-8, 0.74
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


def _figiel_provider_model(
    components: tuple[str, ...] = ("water", "sodium-cation", "chloride-anion"),
) -> epcsaft.Mixture:
    fixture = (
        Path(__file__).parents[1]
        / "fixtures/data/figiel-2025-reference-electrolytes/2/parameters/parameters.json"
    )
    parameters = epcsaft.Parameters.from_json(fixture, components=components)
    return epcsaft.Mixture(parameters)


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
    phase = equilibrium.ProviderPhase(model, (1.0e-6, 0.74))
    problem = _multineutral_problem()

    result = equilibrium.chemical_equilibrium(
        phase,
        313.15 * epcsaft.unit_registry.kelvin,
        119_199.78746707023 * epcsaft.unit_registry.pascal,
        problem,
    )

    assert min(result.amounts_mol) > problem.strict_interior_amount_floor_mol
    assert result.amounts_mol[1] == pytest.approx(0.2, abs=2.0e-10)
    assert result.amounts_mol[2] == pytest.approx(result.amounts_mol[3], abs=2.0e-10)
    assert result.diagnostics.derivative_audits_performed
    assert result.diagnostics.reaction_affinity_inf_norm <= 2.0e-7
    assert result.diagnostics.kkt_stationarity_inf_norm <= 2.0e-7


def test_active_parameter_values_own_reaction_setup_geometry() -> None:
    """Falsify seeding an active solve from the inactive packing geometry."""

    parameters = epcsaft.Parameters.from_mapping(MULTINEUTRAL_REACTIVE)
    model = epcsaft.Mixture(parameters)
    identity = "component/solvent/segment_diameter"
    active = epcsaft.ActiveParameterSet(
        parameters,
        (identity,),
        values=(100.0 * epcsaft.unit_registry.angstrom,),
    )
    result = equilibrium.chemical_equilibrium(
        equilibrium.ProviderPhase(model, (1.0e-6, 0.74)),
        313.15 * epcsaft.unit_registry.kelvin,
        119_199.78746707023 * epcsaft.unit_registry.pascal,
        _multineutral_problem(),
        sensitivity_request=equilibrium.ChemicalEquilibriumSensitivityRequest(
            active_parameters=active
        ),
    )

    assert result.sensitivity is not None
    assert result.sensitivity.parameters[-1].identity == identity
    active_model = model._native.with_active_coordinates(model._active_coordinates(active))
    molar_bounds = active_model.molar_volume_bounds(
        313.15,
        result.mole_fractions,
        1.0e-6,
        0.74,
        math.nan,
    )
    total_amount = math.fsum(result.amounts_mol)
    molar_volume = result.volume_m3 / total_amount
    assert molar_bounds[0] < molar_volume < molar_bounds[1]


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


def _provider_solve(
    model: epcsaft.Mixture,
    spec: dict[str, object],
    source_standard_state: dict[str, object] | None = None,
    sensitivity_request: equilibrium.ChemicalEquilibriumSensitivityRequest | None = None,
) -> equilibrium.ChemicalEquilibriumResult:
    standard_state = (
        None if source_standard_state is None else _typed_standard_state(source_standard_state)
    )
    return equilibrium.chemical_equilibrium(
        equilibrium.ProviderPhase(
            model=model,
            admissible_packing_fraction_interval=(1.0e-6, 0.74),
        ),
        spec["temperature_k"] * epcsaft.unit_registry.kelvin,
        spec["pressure_pa"] * epcsaft.unit_registry.pascal,
        _typed_problem(spec, source_standard_state=standard_state),
        sensitivity_request=sensitivity_request,
    )


def _provider_native(
    model: epcsaft.Mixture,
    spec: dict[str, object],
    source_standard_state: dict[str, object] | None = None,
    packing_fraction_bounds: tuple[float, float] = (1.0e-6, 0.74),
) -> dict[str, object]:
    native_spec = copy.deepcopy(spec)
    if "conserved_totals" not in native_spec:
        balances = native_spec["balance_matrix"]
        feed = native_spec["feed_amounts"]
        native_spec["conserved_totals"] = tuple(
            math.fsum(
                row[index] * feed[index]  # type: ignore[index]
                for index in range(len(feed))  # type: ignore[arg-type]
            )
            for row in balances  # type: ignore[union-attr]
        )
    return _equilibrium._chemical_equilibrium(
        _native_transport(model),
        native_spec,
        source_standard_state,
        packing_fraction_bounds,
        1.0e-12,
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
    result = _provider_solve(
        model,
        spec,
        source_standard_state,
    )

    amounts = result.amounts_mol
    assert 2.0 * amounts[0] + 3.0 * amounts[1] + amounts[2] == pytest.approx(2.0, abs=2.0e-9)
    assert sum(amounts) == pytest.approx(1.0, abs=2.0e-9)
    assert amounts[1] == pytest.approx(amounts[2], abs=2.0e-12)
    assert result.source_standard_state.id == ("iapws-molality-ions-mole-fraction-water")
    assert result.provider_reference_id == (
        "A_over_RT_reference_amount:n_ref=1mol:rho_ref=1mol_per_m3"
    )
    diagnostics = result.diagnostics
    assert diagnostics.reference_derivative_availability == 0
    assert diagnostics.reference_convergence_error <= 5.0e-5
    assert diagnostics.reference_representation_residual_inf_norm <= 1.0e-12
    assert (
        result.source_standard_state.activity_scale_id == source_standard_state["activity_scale_id"]
    )
    assert result.standard_offsets == pytest.approx((-21.20037599919826,), abs=2.0e-12)
    assert result.ln_k_provider_basis == pytest.approx((-53.42391841715451,), abs=2.0e-12)
    assert diagnostics.reaction_affinity_inf_norm <= 2.0e-12
    assert diagnostics.kkt_stationarity_inf_norm <= 2.0e-12
    assert result.sensitivity is None
    transfer = result.source_reference_transfer
    assert transfer is not None
    assert transfer.pressure_pa == spec["pressure_pa"]
    assert transfer.source_reference_pressure_pa == 100_000.0
    assert transfer.parameter_fingerprint == model.parameter_fingerprint
    for fingerprint in (
        transfer.topology_fingerprint,
        transfer.component_order_fingerprint,
        transfer.reference_state_fingerprint,
        transfer.domain_fingerprint,
        transfer.artifact_fingerprint,
    ):
        assert fingerprint.startswith("sha256:")
    assert result.source_standard_state.reference_pressure_pa == 100_000.0
    artifact = result.artifact_identity
    assert artifact.provider_record_sha256.startswith("sha256:")


@pytest.mark.slow
def test_source_reference_pressure_is_distinct_from_trial_pressure() -> None:
    model, spec, source_standard_state = _held_water_ionization_problem()
    source_standard_state["reference_pressure_pa"] = 100_001.0

    value = _provider_solve(model, spec, source_standard_state)

    assert value.pressure_pa == 100_000.0
    assert value.source_reference_transfer is not None
    assert value.source_reference_transfer.pressure_pa == 100_000.0
    assert value.source_standard_state is not None
    assert value.source_standard_state.reference_pressure_pa == 100_001.0
    differentiated = _provider_solve(
        model,
        spec,
        source_standard_state,
        equilibrium.ChemicalEquilibriumSensitivityRequest(),
    )
    assert differentiated.sensitivity is not None
    assert differentiated.sensitivity.status == "available"
    assert differentiated.diagnostics.reference_derivative_availability == 2


@pytest.mark.slow
@pytest.mark.parametrize("variant", ("species_order", "reaction_orientation"))
def test_held_source_transform_is_coordinate_invariant(variant: str) -> None:
    model, spec, source_standard_state = _held_water_ionization_problem()
    if variant == "species_order":
        permutation = (2, 0, 1)
        components = tuple(spec["species_ids"][index] for index in permutation)
        parameters = _held_parameters(components)
        model = epcsaft.Mixture(parameters)
        for field in (
            "species_ids",
            "charges",
            "molar_masses_kg_per_mol",
            "feed_amounts",
        ):
            values = spec[field]
            spec[field] = tuple(values[index] for index in permutation)
        spec["balance_matrix"] = tuple(
            tuple(row[index] for index in permutation) for row in spec["balance_matrix"]
        )
        spec["reaction_matrix"] = tuple(
            tuple(row[index] for index in permutation) for row in spec["reaction_matrix"]
        )
        spec["provider_fingerprint"] = model.parameter_fingerprint
        scales = source_standard_state["log_activity_scale_factors"]
        source_standard_state["log_activity_scale_factors"] = tuple(
            scales[index] for index in permutation
        )
        reference_ids = source_standard_state["source_reference_component_ids"]
        source_standard_state["source_reference_component_ids"] = tuple(
            reference_ids[index] for index in permutation
        )
        solvent_composition = source_standard_state["source_reference_solvent_composition"]
        source_standard_state["source_reference_solvent_composition"] = tuple(
            solvent_composition[index] for index in permutation
        )
        water_index, cation_index, anion_index = (1, 2, 0)
    else:
        spec["reaction_matrix"] = tuple(
            tuple(-value for value in row) for row in spec["reaction_matrix"]
        )
        spec["ln_k"] = tuple(-value for value in spec["ln_k"])
        water_index, cation_index, anion_index = (0, 1, 2)

    result = _provider_solve(
        model,
        spec,
        source_standard_state,
        equilibrium.ChemicalEquilibriumSensitivityRequest(),
    )

    assert result.amounts_mol[water_index] == pytest.approx(0.9999999963727326, rel=2.0e-10)
    assert result.amounts_mol[cation_index] == pytest.approx(
        result.amounts_mol[anion_index], abs=2.0e-12
    )
    assert result.sensitivity is not None
    assert result.sensitivity.status == "available"
    assert result.sensitivity.pressure.status == "available"


def test_engine_model_provider_solve_enforces_balance_affinity_and_domain() -> None:
    model = _figiel_provider_model()
    temperature_k = 298.15
    pressure_pa = 1_037_015.835643739
    spec = {
        "species_ids": ("water", "sodium-cation", "chloride-anion"),
        "charges": (0, 1, -1),
        "molar_masses_kg_per_mol": (2.0, 1.0, 1.0),
        "provider_fingerprint": model.parameter_fingerprint,
        "balance_matrix": ((2.0, 1.0, 1.0), (0.0, 1.0, -1.0)),
        "reaction_matrix": ((-1.0, 1.0, 1.0),),
        "feed_amounts": (0.8, 0.1, 0.1),
        "ln_k": (-279.92590115413996,),
        "temperature_k": temperature_k,
        "pressure_pa": pressure_pa,
    }
    _bind_record(spec)

    result = _provider_solve(model, spec)

    assert 2.0 * result.amounts_mol[0] + result.amounts_mol[1] + result.amounts_mol[
        2
    ] == pytest.approx(1.8, abs=2.0e-9)
    assert result.amounts_mol[1] == pytest.approx(result.amounts_mol[2], abs=2.0e-12)
    assert result.diagnostics.reaction_affinity_inf_norm <= 2.0e-7
    assert result.diagnostics.pressure_relative_residual <= 1.0e-8
    assert result.provider_parameter_fingerprint == model.parameter_fingerprint
    assert 1.0e-6 < result.diagnostics.packing_fraction < 0.74
    rejected = _provider_native(
        model,
        spec,
        packing_fraction_bounds=(0.73, 0.74),
    )
    assert rejected["accepted"] is False


@pytest.mark.slow
def test_engine_model_pressure_sensitivity_matches_independent_resolves() -> None:
    model = _figiel_provider_model()
    temperature_k = 298.15
    pressure_pa = 1_037_015.835643739
    spec = {
        "species_ids": ("water", "sodium-cation", "chloride-anion"),
        "charges": (0, 1, -1),
        "molar_masses_kg_per_mol": (2.0, 1.0, 1.0),
        "provider_fingerprint": model.parameter_fingerprint,
        "balance_matrix": ((2.0, 1.0, 1.0), (0.0, 1.0, -1.0)),
        "reaction_matrix": ((-1.0, 1.0, 1.0),),
        "feed_amounts": (0.8, 0.1, 0.1),
        "ln_k": (-279.92590115413996,),
        "temperature_k": temperature_k,
        "pressure_pa": pressure_pa,
    }
    _bind_record(spec)
    native = _provider_native(model, spec)
    sensitivities = native["sensitivities"]
    public = _provider_solve(
        model,
        spec,
        sensitivity_request=equilibrium.ChemicalEquilibriumSensitivityRequest(),
    )
    public_sensitivity = public.sensitivity
    assert public_sensitivity is not None
    assert sensitivities["parameter_fingerprint"] == model.parameter_fingerprint
    pressure_index = sensitivities["parameter_order"].index("pressure_pa")

    step_pa = 5.0
    perturbed: list[dict[str, object]] = []
    for direction in (-1.0, 1.0):
        trial = copy.deepcopy(spec)
        trial["pressure_pa"] = spec["pressure_pa"] + direction * step_pa
        _bind_record(trial)
        perturbed.append(_provider_native(model, trial))
    finite_difference_amounts = tuple(
        (perturbed[1]["amounts"][species] - perturbed[0]["amounts"][species]) / (2.0 * step_pa)
        for species in range(3)
    )
    finite_difference_volume = (perturbed[1]["volume_m3"] - perturbed[0]["volume_m3"]) / (
        2.0 * step_pa
    )
    assert sensitivities["amount_derivatives"][pressure_index] == pytest.approx(
        finite_difference_amounts,
        rel=2.0e-5,
        abs=2.0e-12,
    )
    assert sensitivities["volume_derivatives"][pressure_index] == pytest.approx(
        finite_difference_volume,
        rel=2.0e-5,
        abs=2.0e-13,
    )
