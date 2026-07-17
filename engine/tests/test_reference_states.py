import math

import pytest

from epcsaft import EPCSAFT, ParameterBundle, unit_registry


def _model(*component_ids: str) -> EPCSAFT:
    bundle = ParameterBundle.from_catalog("figiel-2025-reference-electrolytes", version=1)
    return EPCSAFT(bundle.select(component_ids))


def test_water_nacl_reference_converges_on_an_electroneutral_path() -> None:
    model = _model("water", "sodium-cation", "chloride-anion")

    reference = model.reference_state(
        temperature=298.15 * unit_registry.kelvin,
        pressure=1.0 * unit_registry.bar,
        solvent_mole_fractions=(1.0, 0.0, 0.0),
        phase="liquid",
    )

    assert reference.component_ids == model.component_ids
    assert reference.ion_component_ids == ("sodium-cation", "chloride-anion")
    assert reference.formula_stoichiometry == (1, 1)
    assert reference.solvent_mole_fractions == (1.0, 0.0, 0.0)
    assert reference.reference_molality.to("mole / kilogram").magnitude <= 1.0e-10
    assert reference.convergence_error <= 5.0e-5
    assert all(math.isfinite(value) for value in reference.infinite_dilution_log_fugacity)


def test_activity_follows_figiel_reference_and_scale_equations() -> None:
    model = _model("water", "sodium-cation", "chloride-anion")
    reference = model.reference_state(
        temperature=298.15 * unit_registry.kelvin,
        pressure=1.0 * unit_registry.bar,
        solvent_mole_fractions=(1.0, 0.0, 0.0),
        phase="liquid",
    )

    activity = reference.at_molality(0.1 * unit_registry.mole / unit_registry.kilogram)
    state = model.evaluate(
        temperature=reference.temperature,
        pressure=reference.pressure,
        mole_fractions=activity.mole_fractions,
        phase="liquid",
    )
    assert state.log_fugacity_coefficient is not None
    expected_single = tuple(
        math.exp(state.log_fugacity_coefficient[index] - reference.infinite_dilution_log_fugacity[k])
        for k, index in enumerate((1, 2))
    )
    assert activity.single_ion_activity_coefficients == pytest.approx(expected_single)

    expected_mean_x = math.sqrt(expected_single[0] * expected_single[1])
    expected_mean_m = expected_mean_x / (1.0 + 0.0180153 * 0.1 * 2.0)
    assert activity.mean_ionic_activity_coefficient_mole_fraction == pytest.approx(expected_mean_x)
    assert activity.mean_ionic_activity_coefficient_molality == pytest.approx(expected_mean_m)
    assert math.fsum(activity.mole_fractions) == pytest.approx(1.0)
    assert activity.mole_fractions[1] == activity.mole_fractions[2]


def test_solvation_free_energy_follows_figiel_x_process() -> None:
    reference = _model("water", "sodium-cation", "chloride-anion").reference_state(
        temperature=298.15 * unit_registry.kelvin,
        pressure=1.0 * unit_registry.bar,
        solvent_mole_fractions=(1.0, 0.0, 0.0),
        phase="liquid",
    )

    gas_constant = 8.31446261815324
    expected = tuple(
        gas_constant * 298.15 * value for value in reference.infinite_dilution_log_fugacity
    )
    observed = tuple(
        value.to("joule / mole").magnitude for value in reference.solvation_free_energy
    )
    assert observed == pytest.approx(expected)


def test_transfer_free_energy_follows_figiel_ratio_and_reverses_sign() -> None:
    model = _model("water", "methanol", "sodium-cation", "bromide-anion")
    common = {
        "temperature": 298.15 * unit_registry.kelvin,
        "pressure": 1.0 * unit_registry.bar,
        "phase": "liquid",
    }
    water = model.reference_state(
        **common,
        solvent_mole_fractions=(1.0, 0.0, 0.0, 0.0),
    )
    methanol = model.reference_state(
        **common,
        solvent_mole_fractions=(0.0, 1.0, 0.0, 0.0),
    )

    forward = water.transfer_to(methanol)
    reverse = methanol.transfer_to(water)
    identity = water.transfer_to(water)
    scale = 8.31446261815324 * 298.15
    expected = tuple(
        scale * (target - source)
        for source, target in zip(
            water.infinite_dilution_log_fugacity,
            methanol.infinite_dilution_log_fugacity,
            strict=True,
        )
    )

    forward_values = tuple(value.to("joule / mole").magnitude for value in forward.values)
    reverse_values = tuple(value.to("joule / mole").magnitude for value in reverse.values)
    identity_values = tuple(value.to("joule / mole").magnitude for value in identity.values)
    assert forward.ion_component_ids == ("sodium-cation", "bromide-anion")
    assert forward_values == pytest.approx(expected)
    assert reverse_values == pytest.approx(tuple(-value for value in expected))
    assert identity_values == pytest.approx((0.0, 0.0), abs=1.0e-15)


def test_nacl_mean_ionic_activity_matches_hamer_wu_at_one_molal() -> None:
    reference = _model("water", "sodium-cation", "chloride-anion").reference_state(
        temperature=298.15 * unit_registry.kelvin,
        pressure=1.0 * unit_registry.bar,
        solvent_mole_fractions=(1.0, 0.0, 0.0),
        phase="liquid",
    )

    # Hamer and Wu (1972), DOI 10.1063/1.3253108, Table 16:
    # aqueous NaCl at 25 degC and m = 1.000 mol/kg has gamma_pm = 0.657.
    observed = reference.at_molality(
        1.0 * unit_registry.mole / unit_registry.kilogram
    ).mean_ionic_activity_coefficient_molality
    assert abs(observed - 0.657) <= 0.08


def test_nacl_activity_in_water_ethanol_matches_esteso() -> None:
    model = _model("water", "ethanol", "sodium-cation", "chloride-anion")
    water_moles = 0.4 / 0.0180153
    ethanol_moles = 0.6 / 0.046069
    solvent_moles = water_moles + ethanol_moles
    reference = model.reference_state(
        temperature=298.15 * unit_registry.kelvin,
        pressure=1.0 * unit_registry.bar,
        solvent_mole_fractions=(
            water_moles / solvent_moles,
            ethanol_moles / solvent_moles,
            0.0,
            0.0,
        ),
        phase="liquid",
    )

    # Esteso et al. (1989), J. Solution Chem. 18, Table I:
    # 60 mass % ethanol in the salt-free solvent, m = 0.500 mol/kg,
    # and gamma_pm,m = 0.415. The 0.08 gate is a model-accuracy tolerance,
    # not the experimental precision.
    observed = reference.at_molality(
        0.5 * unit_registry.mole / unit_registry.kilogram
    ).mean_ionic_activity_coefficient_molality
    assert abs(observed - 0.415) <= 0.08
