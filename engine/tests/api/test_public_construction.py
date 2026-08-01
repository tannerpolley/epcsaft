from __future__ import annotations

from pathlib import Path

import pytest

from epcsaft import Mixture, ParameterError, Parameters

from tests._raw_cases import NEUTRAL_BINARY


def test_parameters__exposes_only_dictionary_and_bundle_path_loaders(tmp_path: Path) -> None:
    assert not hasattr(Parameters, "from_catalog")
    assert not hasattr(Parameters, "from_records")

    source = Parameters.from_dictionary(NEUTRAL_BINARY)
    bundle = tmp_path / "external-bundle"
    source.to_bundle(bundle)
    restored = Parameters.from_bundle(bundle)

    assert restored.component_ids == source.component_ids
    assert restored.fingerprint == source.fingerprint


def test_mixture_from_dictionary__delegates_to_strict_parameter_loader() -> None:
    model = Mixture.from_dictionary(NEUTRAL_BINARY)

    assert model.component_ids == ("methane", "ethane")
    assert model.parameter_fingerprint == Parameters.from_dictionary(NEUTRAL_BINARY).fingerprint


def test_bundle_template__is_created_without_catalog_lookup(tmp_path: Path) -> None:
    target = tmp_path / "authoring-template"

    Parameters.create_bundle_template(
        target,
        components=("component-a", "component-b"),
        profile="neutral",
    )

    assert (target / "bundle.toml").is_file()
    assert sorted(path.name for path in target.iterdir()) == [
        "bundle.toml",
        "components.csv",
        "model.csv",
        "single.csv",
    ]


def test_installed_catalog_purpose__normalizes_metadata_without_lookup(
    tmp_path: Path,
) -> None:
    parameters = Parameters.from_dictionary(NEUTRAL_BINARY)
    external_bundle = tmp_path / "external-bundle"
    legacy_bundle = tmp_path / "legacy-data-v1-bundle"
    parameters.to_bundle(external_bundle)
    parameters.to_bundle(legacy_bundle)

    for bundle, purpose in (
        (external_bundle, "external-bundle"),
        (legacy_bundle, "installed-catalog"),
    ):
        manifest = bundle / "bundle.toml"
        manifest.write_text(
            manifest.read_text().replace(
                'purpose = "user-provided"', f'purpose = "{purpose}"'
            )
        )

    external = Parameters.from_bundle(external_bundle)
    legacy = Parameters.from_bundle(legacy_bundle)
    normalized_bundle = tmp_path / "normalized-bundle"
    legacy.to_bundle(normalized_bundle)

    assert legacy.fingerprint == external.fingerprint
    assert legacy.component_ids == external.component_ids
    assert legacy.records == external.records
    assert 'purpose = "external-bundle"' in (
        normalized_bundle / "bundle.toml"
    ).read_text()


def test_external_bundle_purpose__rejects_unknown(tmp_path: Path) -> None:
    parameters = Parameters.from_dictionary(NEUTRAL_BINARY)
    bundle = tmp_path / "external-bundle"
    parameters.to_bundle(bundle)
    manifest = bundle / "bundle.toml"

    manifest.write_text(
        manifest.read_text().replace('purpose = "user-provided"', 'purpose = "typo"')
    )
    with pytest.raises(ParameterError, match="invalid-value"):
        Parameters.from_bundle(bundle)
