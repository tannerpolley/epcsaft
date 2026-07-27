import math
from pathlib import Path

import pytest

from epcsaft import (
    EPCSAFT,
    ElectrolyteReference,
    IonicActivity,
    ParameterBundle,
    TransferFreeEnergy,
    unit_registry,
)


def _model(*component_ids: str) -> EPCSAFT:
    bundle = ParameterBundle.from_catalog("figiel-2025-reference-electrolytes", version=1)
    return EPCSAFT(bundle.select(component_ids))


def _held_water_ionization_model() -> EPCSAFT:
    bundle = ParameterBundle.from_catalog("held-2008-water-self-ionization", version=1)
    return EPCSAFT(bundle.select(("water", "hydronium-cation", "hydroxide-anion")))


def test_held_water_ionization_reference_uses_the_native_sequence_owner() -> None:
    model = _held_water_ionization_model()

    reference = model.reference_state(
        temperature=298.15 * unit_registry.kelvin,
        pressure=100000.0 * unit_registry.pascal,
        solvent_mole_fractions=(1.0, 0.0, 0.0),
        phase="liquid",
    )
    pure_water = model.evaluate(
        temperature=reference.temperature,
        pressure=reference.pressure,
        mole_fractions=(1.0, 0.0, 0.0),
        phase="liquid",
    )

    assert reference.ion_component_ids == ("hydronium-cation", "hydroxide-anion")
    assert reference.formula_stoichiometry == (1, 1)
    assert reference.reference_molality.to("mole / kilogram").magnitude == 1.0e-12
    assert 0.0 <= reference.convergence_error <= 5.0e-5
    assert all(
        math.isfinite(value)
        for value in reference.infinite_dilution_log_fugacity_coefficient
    )
    assert pure_water.log_fugacity_coefficient is not None
    assert reference._pure_water_log_fugacity == pytest.approx(
        pure_water.log_fugacity_coefficient[0], abs=2.0e-12
    )


def test_water_nacl_reference_converges_on_an_electroneutral_path() -> None:
    model = _model("water", "sodium-cation", "chloride-anion")

    reference = model.reference_state(
        temperature=298.15 * unit_registry.kelvin,
        pressure=1.0 * unit_registry.bar,
        solvent_mole_fractions=(1.0, 0.0, 0.0),
        phase="liquid",
    )

    assert isinstance(reference, ElectrolyteReference)
    assert reference.component_ids == model.component_ids
    assert reference.ion_component_ids == ("sodium-cation", "chloride-anion")
    assert reference.formula_stoichiometry == (1, 1)
    assert reference.solvent_mole_fractions == (1.0, 0.0, 0.0)
    assert reference.reference_molality.to("mole / kilogram").magnitude <= 1.0e-10
    assert reference.convergence_error <= 5.0e-5
    assert all(
        math.isfinite(value) for value in reference.infinite_dilution_log_fugacity_coefficient
    )
    assert not hasattr(reference, "infinite_dilution_log_fugacity")
    assert reference.infinite_dilution_log_fugacity_coefficient == pytest.approx(
        (-154.0679654864907, -127.07688974125972),
        abs=2.0e-12,
    )


def test_activity_follows_figiel_reference_and_scale_equations() -> None:
    model = _model("water", "sodium-cation", "chloride-anion")
    reference = model.reference_state(
        temperature=298.15 * unit_registry.kelvin,
        pressure=1.0 * unit_registry.bar,
        solvent_mole_fractions=(1.0, 0.0, 0.0),
        phase="liquid",
    )

    activity = reference.at_molality(0.1 * unit_registry.mole / unit_registry.kilogram)
    assert isinstance(activity, IonicActivity)
    assert activity.component_ids == model.component_ids
    state = model.evaluate(
        temperature=reference.temperature,
        pressure=reference.pressure,
        mole_fractions=activity.mole_fractions,
        phase="liquid",
    )
    assert state.log_fugacity_coefficient is not None
    expected_single = tuple(
        math.exp(
            state.log_fugacity_coefficient[index]
            - reference.infinite_dilution_log_fugacity_coefficient[k]
        )
        for k, index in enumerate((1, 2))
    )
    assert activity.single_ion_activity_coefficients == pytest.approx(expected_single)

    expected_mean_x = math.sqrt(expected_single[0] * expected_single[1])
    expected_mean_m = expected_mean_x / (1.0 + 0.0180153 * 0.1 * 2.0)
    pure_solvent = model.evaluate(
        temperature=reference.temperature,
        pressure=reference.pressure,
        mole_fractions=reference.solvent_mole_fractions,
        phase="liquid",
    )
    assert pure_solvent.log_fugacity_coefficient is not None
    log_solvent_activity = (
        math.log(activity.mole_fractions[0])
        + state.log_fugacity_coefficient[0]
        - pure_solvent.log_fugacity_coefficient[0]
    )
    expected_osmotic = -log_solvent_activity / (2.0 * 0.1 * 0.0180153)
    assert activity.mean_ionic_activity_coefficient_mole_fraction == pytest.approx(expected_mean_x)
    assert activity.mean_ionic_activity_coefficient_molality == pytest.approx(expected_mean_m)
    assert activity.practical_osmotic_coefficient_molality == pytest.approx(expected_osmotic)
    assert math.fsum(activity.mole_fractions) == pytest.approx(1.0)
    assert activity.mole_fractions[1] == activity.mole_fractions[2]


def test_aqueous_osmotic_coefficient_approaches_ideal_dilution_limit() -> None:
    reference = _model("water", "sodium-cation", "chloride-anion").reference_state(
        temperature=298.15 * unit_registry.kelvin,
        pressure=0.101 * unit_registry.megapascal,
        solvent_mole_fractions=(1.0, 0.0, 0.0),
        phase="liquid",
    )

    observed = reference.at_molality(
        1.0e-6 * unit_registry.mole / unit_registry.kilogram
    ).practical_osmotic_coefficient_molality

    assert observed is not None
    assert abs(observed - 1.0) <= 2.0e-3


@pytest.mark.parametrize("cation", ("sodium-cation", "potassium-cation"))
def test_source_complete_iodides_expose_aqueous_osmotic_coefficient(cation: str) -> None:
    bundle = ParameterBundle.from_catalog("figiel-2025-aqueous-iodides", version=1)
    reference = EPCSAFT(bundle.select(("water", cation, "iodide-anion"))).reference_state(
        temperature=298.15 * unit_registry.kelvin,
        pressure=1.0 * unit_registry.bar,
        solvent_mole_fractions=(1.0, 0.0, 0.0),
        phase="liquid",
    )

    observed = reference.at_molality(
        0.1 * unit_registry.mole / unit_registry.kilogram
    ).practical_osmotic_coefficient_molality

    assert observed is not None
    assert math.isfinite(observed) and observed > 0.0


def test_solvation_free_energy_follows_figiel_x_process() -> None:
    reference = _model("water", "sodium-cation", "chloride-anion").reference_state(
        temperature=298.15 * unit_registry.kelvin,
        pressure=1.0 * unit_registry.bar,
        solvent_mole_fractions=(1.0, 0.0, 0.0),
        phase="liquid",
    )

    gas_constant = 8.31446261815324
    expected = tuple(
        gas_constant * 298.15 * value
        for value in reference.infinite_dilution_log_fugacity_coefficient
    )
    observed = tuple(
        value.to("joule / mole").magnitude for value in reference.solvation_free_energy
    )
    assert observed == pytest.approx(expected)


def test_reference_sequence_algorithm_has_one_native_owner() -> None:
    package_root = Path(__file__).parents[1] / "src" / "epcsaft"
    python_source = (package_root / "reference.py").read_text()
    native_source = (package_root / "_native" / "electrolyte_reference.hpp").read_text()
    eos_source = (package_root / "_native" / "eos.cpp").read_text()

    assert "_REFERENCE_MOLALITIES" not in python_source
    assert "_REFERENCE_CONVERGENCE_TOLERANCE" not in python_source
    assert "def _formula_composition" not in python_source
    assert "for molality in" not in python_source
    assert "REFERENCE_MOLALITIES" in native_source
    assert "REFERENCE_CONVERGENCE_TOLERANCE" in native_source
    assert "formula_composition" in native_source
    assert "electrolyte_reference_sequence" in native_source
    assert eos_source.count("detail::electrolyte_reference_sequence(") == 4
    for delegate in (
        python_source,
        eos_source,
        (package_root / "_native" / "module.cpp").read_text(),
        (package_root / "_native" / "native_sdk.cpp").read_text(),
    ):
        assert "1.0e-6, 1.0e-8, 1.0e-10, 1.0e-12" not in delegate


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
    assert isinstance(forward, TransferFreeEnergy)
    scale = 8.31446261815324 * 298.15
    expected = tuple(
        scale * (target - source)
        for source, target in zip(
            water.infinite_dilution_log_fugacity_coefficient,
            methanol.infinite_dilution_log_fugacity_coefficient,
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
    activity = reference.at_molality(0.5 * unit_registry.mole / unit_registry.kilogram)
    assert abs(activity.mean_ionic_activity_coefficient_molality - 0.415) <= 0.08
    assert activity.practical_osmotic_coefficient_molality is None
