import math
from decimal import Decimal
from pathlib import Path
from typing import Literal

import pytest

from epcsaft import EPCSAFT, ParameterBundle, ParameterError, unit_registry as u
from epcsaft.records import (
    AssociationParameterRecord,
    ComponentRecord,
    ConstantCorrelation,
    ConstantPlusSumOfExponentialsCorrelation,
    ExponentialTerm,
    ModelParameterRecord,
    PairParameterRecord,
    SingleParameterRecord,
    SiteRecord,
    SourceRecord,
    ValidityDomain,
)


def _metadata(*, purpose: Literal["package-test-fixture", "user-provided"]) -> ParameterBundle:
    return ParameterBundle.from_records(
        bundle_id="user-propane",
        bundle_version=1,
        purpose=purpose,
        sources=(
            SourceRecord(
                "user-source",
                "User-supplied parameter source",
                "user-provided model input",
            ),
        ),
        domains=(ValidityDomain("user-domain", "unknown"),),
        components=(ComponentRecord("propane", name="propane", formula="C3H8"),),
    )


def test_user_provided_bundle_purpose_is_fingerprinted() -> None:
    user_bundle = _metadata(purpose="user-provided")
    fixture_bundle = _metadata(purpose="package-test-fixture")

    assert user_bundle.purpose == "user-provided"
    assert user_bundle.fingerprint != fixture_bundle.fingerprint


def _all_record_families_bundle() -> ParameterBundle:
    provenance = {
        "source_id": "user-source",
        "locator": "User parameter worksheet",
        "domain_id": "user-domain",
    }
    return ParameterBundle.from_records(
        bundle_id="user-all-records",
        bundle_version=1,
        purpose="user-provided",
        sources=(
            SourceRecord(
                "user-source",
                "User-supplied parameters",
                "user-provided model input",
                doi="10.0000/example",
            ),
        ),
        domains=(
            ValidityDomain(
                "user-domain",
                "reported-conditions",
                temperature_min=Decimal("250.0000000000000000001") * u.kelvin,
                temperature_max=350 * u.kelvin,
                pressure_min=1 * u.bar,
                pressure_max=100 * u.bar,
                ion_mole_fraction_max=Decimal("0.2000000000000000001"),
            ),
        ),
        components=(
            ComponentRecord("methane", name="methane", formula="CH4", aliases=("CH4",)),
            ComponentRecord("water", name="water", formula="H2O"),
        ),
        singles=(
            SingleParameterRecord(
                "methane-segment-count",
                "methane",
                "segment_count",
                1.0,
                **provenance,
            ),
            SingleParameterRecord(
                "methane-dispersion-energy",
                "methane",
                "dispersion_energy_over_k",
                150 * u.kelvin,
                **provenance,
            ),
            SingleParameterRecord(
                "methane-charge",
                "methane",
                "charge_number",
                0,
                **provenance,
            ),
        ),
        pairs=(
            PairParameterRecord(
                "methane-water-kij",
                "methane",
                "water",
                "k_ij",
                0.01,
                **provenance,
            ),
        ),
        sites=(
            SiteRecord("water-site-a", "water", "a", "a", 1, **provenance),
            SiteRecord("water-site-b", "water", "b", "b", 1, **provenance),
        ),
        associations=(
            AssociationParameterRecord(
                "water-association-energy",
                "water",
                "a",
                "water",
                "b",
                "association_energy_over_k",
                1200 * u.kelvin,
                **provenance,
            ),
            AssociationParameterRecord(
                "water-association-volume",
                "water",
                "a",
                "water",
                "b",
                "association_volume",
                0.03,
                **provenance,
            ),
        ),
        correlations=(
            ConstantCorrelation(
                "water-relative-permittivity",
                "water",
                "relative_permittivity",
                Decimal("78.3000000000000000001"),
                **provenance,
            ),
            ConstantPlusSumOfExponentialsCorrelation(
                "water-segment-diameter",
                "water",
                "segment_diameter",
                3.0 * u.angstrom,
                (
                    ExponentialTerm(
                        Decimal("0.1000000000000000001") * u.angstrom,
                        Decimal("-0.0100000000000000001") / u.kelvin,
                    ),
                ),
                **provenance,
            ),
        ),
        models=(
            ModelParameterRecord(
                "user-ionic-permittivity",
                "ionic_region_relative_permittivity",
                8.0,
                **provenance,
            ),
            ModelParameterRecord(
                "relative-permittivity-formulation",
                "relative_permittivity_formulation",
                "none",
                **provenance,
            ),
        ),
    )


def _propane_bundle() -> ParameterBundle:
    provenance = {
        "source_id": "user-source",
        "locator": "User parameter worksheet, propane row",
        "domain_id": "propane-domain",
    }
    return ParameterBundle.from_records(
        bundle_id="user-propane",
        bundle_version=1,
        purpose="user-provided",
        sources=(
            SourceRecord(
                "user-source",
                "User-supplied propane parameters",
                "user-provided model input",
            ),
        ),
        domains=(
            ValidityDomain(
                "propane-domain",
                "fit-range",
                temperature_min=85 * u.kelvin,
                temperature_max=523 * u.kelvin,
            ),
        ),
        components=(ComponentRecord("propane", name="propane", formula="C3H8"),),
        singles=(
            SingleParameterRecord(
                "propane-molar-mass",
                "propane",
                "molar_mass",
                0.044096 * u.kilogram / u.mole,
                **provenance,
            ),
            SingleParameterRecord(
                "propane-segment-count",
                "propane",
                "segment_count",
                2.0020,
                **provenance,
            ),
            SingleParameterRecord(
                "propane-segment-diameter",
                "propane",
                "segment_diameter",
                3.6184 * u.angstrom,
                **provenance,
            ),
            SingleParameterRecord(
                "propane-dispersion-energy",
                "propane",
                "dispersion_energy_over_k",
                208.11 * u.kelvin,
                **provenance,
            ),
        ),
        models=(
            ModelParameterRecord(
                "relative-permittivity-formulation",
                "relative_permittivity_formulation",
                "none",
                **provenance,
            ),
        ),
    )


def test_user_bundle_round_trip_preserves_every_record_form(tmp_path: Path) -> None:
    bundle = _all_record_families_bundle()
    target = tmp_path / "user-bundle"

    bundle.to_path(target)
    loaded = ParameterBundle.from_path(target)

    assert loaded.fingerprint == bundle.fingerprint
    assert {type(record) for record in loaded.records} == {
        AssociationParameterRecord,
        ConstantCorrelation,
        ConstantPlusSumOfExponentialsCorrelation,
        ModelParameterRecord,
        PairParameterRecord,
        SingleParameterRecord,
        SiteRecord,
    }
    assert {path.name for path in target.iterdir()} == {
        "association.csv",
        "bundle.toml",
        "components.csv",
        "correlations.toml",
        "model.csv",
        "pair.csv",
        "single.csv",
        "sites.csv",
    }


def test_reloaded_user_bundle_evaluates_a_state(tmp_path: Path) -> None:
    target = tmp_path / "user-propane"
    _propane_bundle().to_path(target)

    loaded = ParameterBundle.from_path(target)
    model = EPCSAFT(loaded.select(("propane",)))
    state = model.evaluate(
        temperature=300 * u.kelvin,
        pressure=1 * u.bar,
        mole_fractions=(1.0,),
        phase="vapor",
    )

    assert model.component_ids == ("propane",)
    assert math.isclose(state.pressure.to("bar").magnitude, 1.0, rel_tol=1e-12)
    assert math.isfinite(state.residual_helmholtz)


def test_to_path_refuses_to_overwrite_an_existing_path(tmp_path: Path) -> None:
    target = tmp_path / "existing"
    target.mkdir()

    with pytest.raises(ParameterError) as exc_info:
        _metadata(purpose="user-provided").to_path(target)

    assert exc_info.value.code == "invalid-value"
    assert not tuple(target.iterdir())


def test_to_path_leaves_no_target_after_a_failed_write(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    target = tmp_path / "user-bundle"

    def fail_rename(self: Path, destination: Path) -> Path:
        raise OSError("injected write failure")

    monkeypatch.setattr(Path, "rename", fail_rename)

    with pytest.raises(ParameterError) as exc_info:
        _metadata(purpose="user-provided").to_path(target)

    assert exc_info.value.code == "invalid-value"
    assert not target.exists()
