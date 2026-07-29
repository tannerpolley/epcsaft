from dataclasses import FrozenInstanceError
import os
from pathlib import Path
from statistics import median
from time import perf_counter

import pytest

import epcsaft
from epcsaft import Mixture, ParameterError, Parameters, State, unit_registry as u


def test_neutral_density_state_exposes_the_frontend_and_residual_bases() -> None:
    parameters = Parameters.from_catalog(
        "gross-2001-methane-ethane",
        version=1,
        components=("methane", "ethane"),
    )
    mixture = Mixture(parameters)

    state = mixture.state(
        T=250.0 * u.kelvin,
        rho=8_000.0 * u.mole / u.meter**3,
        x=(0.4, 0.6),
    )

    assert isinstance(state, State)
    assert parameters.component_ids == mixture.component_ids == state.component_ids
    assert state.hard_chain == pytest.approx(0.8611345413201740, rel=5.0e-12)
    assert state.dispersion == pytest.approx(-1.9282055648091894, rel=5.0e-12)
    assert state.ares(basis="reduced") == pytest.approx(-1.0670710234890153, rel=5.0e-12)
    assert state.ares().to("joule / mole").magnitude == pytest.approx(
        -1.0670710234890153 * 8.31446261815324 * 250.0,
        rel=5.0e-12,
    )
    assert state.compressibility_factor == pytest.approx(0.19067336346364548, rel=5.0e-12)
    assert state.pressure.to("pascal").magnitude == pytest.approx(
        3_170_693.1055920523,
        rel=5.0e-12,
    )
    assert not hasattr(epcsaft, "ParameterBundle")
    assert not hasattr(epcsaft, "ParameterSet")
    assert not hasattr(epcsaft, "EPCSAFT")
    assert not hasattr(epcsaft, "EosResult")
    with pytest.raises(FrozenInstanceError):
        state.temperature = 300.0 * u.kelvin


def test_pressure_state_exposes_caloric_and_fugacity_projections() -> None:
    mixture = Mixture(
        Parameters.from_catalog(
            "gross-2001-methane-ethane",
            components=("methane", "ethane"),
        )
    )

    state = mixture.state(
        T=250.0 * u.kelvin,
        P=1.0 * u.bar,
        x=(0.4, 0.6),
    )

    assert state.molar_density.to("mole / meter ** 3").magnitude == pytest.approx(
        48.506191136225354,
        rel=5.0e-12,
    )
    assert state.ures(basis="reduced") == pytest.approx(-0.014998504370230486)
    assert state.hres(basis="reduced") == pytest.approx(-0.023188162428916304)
    assert state.gres(basis="reduced") == pytest.approx(-0.016391161478456207)
    assert state.sres(basis="reduced") == pytest.approx(-0.006797000950460048)
    assert state.cvres(basis="reduced") == pytest.approx(0.0068262704173949335)
    assert state.cpres(basis="reduced") == pytest.approx(0.03730187183349076)
    assert state.fugacity is not None
    assert state.fugacity.ln_coefficient == pytest.approx(
        (-0.002435251323893545, -0.011989472356340493)
    )
    assert state.fugacity.coefficient == pytest.approx(
        (0.9975677114950541, 0.9880821149835884)
    )
    assert tuple(value.to("pascal").magnitude for value in state.fugacity.value) == pytest.approx(
        (39902.7085, 59284.9269),
        rel=2.0e-9,
    )
    assert state.activity is None


def test_associating_state_retains_the_catalog_property_anchor() -> None:
    mixture = Mixture(
        Parameters.from_catalog(
            "figiel-2025-reference-electrolytes",
            version=1,
            components=("ethanol",),
        )
    )

    state = mixture.state(
        T=298.15 * u.kelvin,
        rho=17_000.0 * u.mole / u.meter**3,
        x=(1.0,),
    )

    assert state.hard_chain == pytest.approx(5.313169940020586, rel=2.0e-11)
    assert state.dispersion == pytest.approx(-8.162193153374211, rel=2.0e-11)
    assert state.association == pytest.approx(-4.73937962762838, rel=2.0e-11)
    assert state.debye_huckel == 0.0
    assert state.born == 0.0


def test_explicit_electrolyte_reference_attaches_activity_to_the_state() -> None:
    mixture = Mixture(
        Parameters.from_catalog(
            "figiel-2025-reference-electrolytes",
            version=1,
            components=("water", "sodium-cation", "chloride-anion"),
        )
    )
    reference = mixture.reference(
        T=298.15 * u.kelvin,
        P=1.0 * u.bar,
        solvent_x=(1.0, 0.0, 0.0),
        phase="liquid",
    )

    state = reference.state(molality=0.1 * u.mole / u.kilogram)

    assert reference.infinite_dilution_log_fugacity_coefficient == pytest.approx(
        (-154.0679654864907, -127.07688974125972),
        abs=2.0e-12,
    )
    assert state.activity is not None
    assert state.activity.single_ion == pytest.approx(
        (0.7886051541872077, 0.7806970382125156)
    )
    assert state.activity.mean_mole_fraction == pytest.approx(0.7846411333807813)
    assert state.activity.mean_molality == pytest.approx(0.7818241739725079)
    assert state.activity.osmotic_coefficient == pytest.approx(0.936012244631617)
    assert state.debye_huckel == pytest.approx(-0.0007777581052593353)
    assert state.born == pytest.approx(-0.5103672648190744)
    assert state.association == pytest.approx(-5.3492574554560335)
    assert state.bulk_relative_permittivity == pytest.approx(76.17297208704825)
    assert tuple(
        value.to("joule / mole").magnitude for value in reference.solvation_free_energy
    ) == pytest.approx(
        (-381927.8660792742, -315018.14912410337),
        rel=5.0e-12,
    )

    # Hamer and Wu (1972), Table 16, aqueous NaCl at 25 degC and 1 mol/kg.
    one_molal = reference.state(molality=1.0 * u.mole / u.kilogram)
    assert one_molal.activity is not None
    assert abs(one_molal.activity.mean_molality - 0.657) <= 0.08


def test_transfer_free_energy_is_antisymmetric_and_zero_for_identity() -> None:
    mixture = Mixture(
        Parameters.from_catalog(
            "figiel-2025-reference-electrolytes",
            version=1,
            components=("water", "methanol", "sodium-cation", "bromide-anion"),
        )
    )
    common = {
        "T": 298.15 * u.kelvin,
        "P": 1.0 * u.bar,
        "phase": "liquid",
    }
    water = mixture.reference(
        **common,
        solvent_x=(1.0, 0.0, 0.0, 0.0),
    )
    methanol = mixture.reference(
        **common,
        solvent_x=(0.0, 1.0, 0.0, 0.0),
    )

    forward = water.transfer_to(methanol)
    reverse = methanol.transfer_to(water)
    identity = water.transfer_to(water)
    forward_values = tuple(value.to("joule / mole").magnitude for value in forward.values)
    reverse_values = tuple(value.to("joule / mole").magnitude for value in reverse.values)

    assert forward_values == pytest.approx(
        (8777.23988451225, 12597.383341393788),
        rel=5.0e-12,
    )
    assert reverse_values == pytest.approx(tuple(-value for value in forward_values))
    assert tuple(value.to("joule / mole").magnitude for value in identity.values) == pytest.approx(
        (0.0, 0.0),
        abs=1.0e-15,
    )


def test_parameter_bundle_template_uses_one_profile_driven_schema(tmp_path: Path) -> None:
    target = tmp_path / "electrolyte-parameters"

    Parameters.create_bundle_template(
        target,
        components=("solvent", "cation", "anion"),
        profile="associating-electrolyte",
        correlations=True,
    )

    assert tuple(sorted(path.name for path in target.iterdir())) == (
        "association.csv",
        "bundle.toml",
        "components.csv",
        "correlations.toml",
        "model.csv",
        "pair.csv",
        "single.csv",
        "sites.csv",
    )
    manifest = (target / "bundle.toml").read_text(encoding="utf-8")
    assert 'purpose = "user-provided"' in manifest
    assert "relative_permittivity_formulation" not in manifest
    with pytest.raises(ParameterError):
        Parameters.create_bundle_template(
            target,
            components=("solvent",),
            profile="neutral",
        )


def test_parameters_round_trip_through_the_strict_bundle_schema(tmp_path: Path) -> None:
    original = Parameters.from_catalog(
        "gross-2001-methane-ethane",
        version=1,
        components=("ethane", "methane"),
    )
    path = tmp_path / "parameters"

    original.to_bundle(path)
    restored = Parameters.from_bundle(
        path,
        components=("ethane", "methane"),
    )

    assert restored.component_ids == original.component_ids
    assert restored.bundle_fingerprint == original.bundle_fingerprint
    assert restored.fingerprint == original.fingerprint


def _timings(call, *, repeats: int = 7) -> tuple[float, float, float]:
    call()
    samples = []
    for _ in range(repeats):
        started = perf_counter()
        call()
        samples.append(perf_counter() - started)
    return median(samples), min(samples), max(samples)


@pytest.mark.skipif(
    os.environ.get("EPCSAFT_TIMING") != "1",
    reason="set EPCSAFT_TIMING=1 to report repeated frontend timings",
)
def test_frontend_unique_computational_path_timings() -> None:
    neutral = Mixture(
        Parameters.from_catalog(
            "gross-2001-methane-ethane",
            components=("methane", "ethane"),
        )
    )
    electrolyte = Mixture(
        Parameters.from_catalog(
            "figiel-2025-reference-electrolytes",
            version=1,
            components=("water", "methanol", "sodium-cation", "bromide-anion"),
        )
    )
    common_reference = {
        "T": 298.15 * u.kelvin,
        "P": 1.0 * u.bar,
        "phase": "liquid",
    }
    water = electrolyte.reference(
        **common_reference,
        solvent_x=(1.0, 0.0, 0.0, 0.0),
    )
    methanol = electrolyte.reference(
        **common_reference,
        solvent_x=(0.0, 1.0, 0.0, 0.0),
    )
    calls = {
        "density_state": lambda: neutral.state(
            T=250.0 * u.kelvin,
            rho=8_000.0 * u.mole / u.meter**3,
            x=(0.4, 0.6),
        ),
        "pressure_state": lambda: neutral.state(
            T=250.0 * u.kelvin,
            P=1.0 * u.bar,
            x=(0.4, 0.6),
        ),
        "reference": lambda: electrolyte.reference(
            **common_reference,
            solvent_x=(1.0, 0.0, 0.0, 0.0),
        ),
        "activity_state": lambda: water.state(molality=0.1 * u.mole / u.kilogram),
        "solvation_projection": lambda: water.solvation_free_energy,
        "transfer_projection": lambda: water.transfer_to(methanol),
    }

    report = {name: _timings(call) for name, call in calls.items()}
    for name, (middle, fastest, slowest) in report.items():
        print(
            f"{name}: median={middle:.6g}s min={fastest:.6g}s max={slowest:.6g}s"
        )
        assert 0.0 <= fastest <= middle <= slowest
