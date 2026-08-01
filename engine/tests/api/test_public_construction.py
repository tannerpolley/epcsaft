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


@pytest.mark.parametrize("purpose", ("installed-catalog", "typo"))
def test_external_bundle_purpose__rejects_legacy_catalog_and_unknown(
    tmp_path: Path,
    purpose: str,
) -> None:
    parameters = Parameters.from_dictionary(NEUTRAL_BINARY)
    bundle = tmp_path / "external-bundle"
    parameters.to_bundle(bundle)
    manifest = bundle / "bundle.toml"
    original = manifest.read_text()

    manifest.write_text(
        original.replace('purpose = "user-provided"', f'purpose = "{purpose}"')
    )
    with pytest.raises(ParameterError, match="invalid-value"):
        Parameters.from_bundle(bundle)
