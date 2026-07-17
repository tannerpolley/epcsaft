from dataclasses import FrozenInstanceError

import pytest

from epcsaft import ParameterBundle, ParameterError, unit_registry
from epcsaft._resolve import _resolve_model
from epcsaft.records import (
    AssociationParameterRecord,
    ComponentRecord,
    ConstantCorrelation,
    ModelParameterRecord,
    PairParameterRecord,
    SingleParameterRecord,
    SiteRecord,
    SourceRecord,
    ValidityDomain,
)


PROVENANCE = {
    "source_id": "source-a",
    "locator": "source fixture",
    "domain_id": "domain-a",
}


def _single(record_id: str, component_id: str, family: str, value: object):
    return SingleParameterRecord(
        record_id, component_id, family, value, **PROVENANCE
    )


def _base_singles(component_id: str, *, molar_mass: float = 0.018):
    return (
        _single(
            f"{component_id}-molar-mass",
            component_id,
            "molar_mass",
            molar_mass * unit_registry.kilogram / unit_registry.mole,
        ),
        _single(f"{component_id}-m", component_id, "segment_count", 1),
        _single(
            f"{component_id}-sigma",
            component_id,
            "segment_diameter",
            3 * unit_registry.angstrom,
        ),
        _single(
            f"{component_id}-epsilon",
            component_id,
            "dispersion_energy_over_k",
            200 * unit_registry.kelvin,
        ),
    )


def _bundle(
    components: tuple[ComponentRecord, ...],
    *,
    singles: tuple[SingleParameterRecord, ...],
    pairs: tuple[PairParameterRecord, ...] = (),
    sites: tuple[SiteRecord, ...] = (),
    associations: tuple[AssociationParameterRecord, ...] = (),
    correlations: tuple[ConstantCorrelation, ...] = (),
    models: tuple[ModelParameterRecord, ...] = (),
):
    return ParameterBundle.from_records(
        bundle_id="resolution-fixture",
        bundle_version=1,
        purpose="package-test-fixture",
        sources=(SourceRecord("source-a", "Source A", "test fixture"),),
        domains=(ValidityDomain("domain-a", "unknown"),),
        components=components,
        singles=singles,
        pairs=pairs,
        sites=sites,
        associations=associations,
        correlations=correlations,
        models=models,
    ).select(tuple(component.component_id for component in components))


def test_neutral_associating_and_ionic_topologies_resolve_immutably() -> None:
    neutral = ParameterBundle.from_catalog(
        "gross-2001-methane-ethane", version=1
    ).select(("methane", "ethane"))
    neutral_model = _resolve_model(neutral)
    assert neutral_model.topology == "neutral"
    assert neutral_model.component_ids == ("methane", "ethane")
    assert neutral_model.site_ids == ()
    assert neutral_model.charge_number == ()

    water = ComponentRecord("water")
    sites = (
        SiteRecord("water-site-a", "water", "a", "a", 1, **PROVENANCE),
        SiteRecord("water-site-b", "water", "b", "b", 1, **PROVENANCE),
    )
    association = (
        AssociationParameterRecord(
            "water-ab-energy",
            "water",
            "a",
            "water",
            "b",
            "association_energy_over_k",
            2400 * unit_registry.kelvin,
            **PROVENANCE,
        ),
        AssociationParameterRecord(
            "water-ab-volume",
            "water",
            "a",
            "water",
            "b",
            "association_volume",
            0.03,
            **PROVENANCE,
        ),
    )
    associating_model = _resolve_model(
        _bundle((water,), singles=_base_singles("water"), sites=sites, associations=association)
    )
    assert associating_model.topology == "associating"
    assert associating_model.site_ids == (("water", "a"), ("water", "b"))
    assert len(associating_model.association_pairs) == 1

    components = (water, ComponentRecord("sodium"), ComponentRecord("chloride"))
    singles = (
        *_base_singles("water"),
        *_base_singles("sodium", molar_mass=0.023),
        *_base_singles("chloride", molar_mass=0.0355),
        _single("water-charge", "water", "charge_number", 0),
        _single("sodium-charge", "sodium", "charge_number", 1),
        _single("chloride-charge", "chloride", "charge_number", -1),
        _single("sodium-born", "sodium", "born_diameter", 3 * unit_registry.angstrom),
        _single("chloride-born", "chloride", "born_diameter", 4 * unit_registry.angstrom),
        _single("water-solvation", "water", "solvation_factor", 1.4),
        _single("sodium-solvation", "sodium", "solvation_factor", 1),
        _single("chloride-solvation", "chloride", "solvation_factor", 1),
    )
    pairs = tuple(
        PairParameterRecord(
            f"{a}-{b}-kij", a, b, "k_ij", 0, **PROVENANCE
        )
        for a, b in (("water", "sodium"), ("water", "chloride"), ("sodium", "chloride"))
    )
    correlations = (
        ConstantCorrelation(
            "water-relative-permittivity",
            "water",
            "relative_permittivity",
            78.3,
            **PROVENANCE,
        ),
    )
    models = (
        ModelParameterRecord(
            "ionic-region-permittivity",
            "ionic_region_relative_permittivity",
            8,
            **PROVENANCE,
        ),
        ModelParameterRecord(
            "dielectric-suppression",
            "dielectric_ion_suppression_coefficient",
            7.01,
            **PROVENANCE,
        ),
    )
    ionic_model = _resolve_model(
        _bundle(
            components,
            singles=singles,
            pairs=pairs,
            sites=sites,
            associations=association,
            correlations=correlations,
            models=models,
        )
    )
    assert ionic_model.topology == "ionic-ssm-ds"
    assert ionic_model.charge_number == (0, 1, -1)
    assert ionic_model.ionic_region_relative_permittivity == 8.0
    assert ionic_model.dielectric_ion_suppression_coefficient == 7.01

    without_unused_ion_masses = tuple(
        record
        for record in singles
        if not (
            record.family == "molar_mass"
            and record.component_id in {"sodium", "chloride"}
        )
    )
    lean_model = _resolve_model(
        _bundle(
            components,
            singles=without_unused_ion_masses,
            pairs=pairs,
            sites=sites,
            associations=association,
            correlations=correlations,
            models=models,
        )
    )
    assert lean_model.molar_mass_kg_per_mol == (0.018, None, None)
    with pytest.raises(FrozenInstanceError):
        ionic_model.topology = "neutral"


@pytest.mark.parametrize(
    "removed_family",
    (
        "molar_mass",
        "segment_count",
        "segment_diameter",
        "dispersion_energy_over_k",
        "charge_number",
        "born_diameter",
        "solvation_factor",
        "relative_permittivity",
        "ionic_region_relative_permittivity",
        "dielectric_ion_suppression_coefficient",
        "association_energy_over_k",
        "association_volume",
        "k_ij",
    ),
)
def test_partial_ionic_topologies_are_rejected(removed_family: str) -> None:
    # Task 4 supplies the physical fixture. This compact synthetic topology
    # exercises completeness without claiming model-output evidence.
    components = (ComponentRecord("water"), ComponentRecord("sodium"), ComponentRecord("chloride"))
    singles = [
        *(_base_singles("water")),
        *(_base_singles("sodium", molar_mass=0.023)),
        *(_base_singles("chloride", molar_mass=0.0355)),
        _single("water-charge", "water", "charge_number", 0),
        _single("sodium-charge", "sodium", "charge_number", 1),
        _single("chloride-charge", "chloride", "charge_number", -1),
        _single("sodium-born", "sodium", "born_diameter", 3 * unit_registry.angstrom),
        _single("chloride-born", "chloride", "born_diameter", 4 * unit_registry.angstrom),
        _single("water-solvation", "water", "solvation_factor", 1.4),
        _single("sodium-solvation", "sodium", "solvation_factor", 1),
        _single("chloride-solvation", "chloride", "solvation_factor", 1),
    ]
    sites = [
        SiteRecord("water-site-a", "water", "a", "a", 1, **PROVENANCE),
        SiteRecord("water-site-b", "water", "b", "b", 1, **PROVENANCE),
    ]
    associations = [
        AssociationParameterRecord(
            "water-ab-energy", "water", "a", "water", "b",
            "association_energy_over_k", 2400 * unit_registry.kelvin, **PROVENANCE
        ),
        AssociationParameterRecord(
            "water-ab-volume", "water", "a", "water", "b",
            "association_volume", 0.03, **PROVENANCE
        ),
    ]
    pairs = [
        PairParameterRecord(f"{a}-{b}-kij", a, b, "k_ij", 0, **PROVENANCE)
        for a, b in (("water", "sodium"), ("water", "chloride"), ("sodium", "chloride"))
    ]
    correlations = [
        ConstantCorrelation(
            "water-relative-permittivity", "water", "relative_permittivity", 78.3, **PROVENANCE
        )
    ]
    models = [
        ModelParameterRecord(
            "ionic-region-permittivity", "ionic_region_relative_permittivity", 8, **PROVENANCE
        ),
        ModelParameterRecord(
            "dielectric-suppression", "dielectric_ion_suppression_coefficient", 7.01, **PROVENANCE
        ),
    ]

    def without_family(records):
        return [record for record in records if record.family != removed_family]

    singles = without_family(singles)
    associations = without_family(associations)
    pairs = without_family(pairs)
    correlations = without_family(correlations)
    models = without_family(models)

    with pytest.raises(ParameterError, match="missing-record|invalid-value"):
        _resolve_model(
            _bundle(
                components,
                singles=tuple(singles),
                pairs=tuple(pairs),
                sites=tuple(sites),
                associations=tuple(associations),
                correlations=tuple(correlations),
                models=tuple(models),
            )
        )


def test_source_backed_figiel_catalog_resolves_without_unpublished_defaults() -> None:
    parameters = ParameterBundle.from_catalog(
        "figiel-2025-water-sodium-chloride", version=1
    ).select(("water", "sodium-cation", "chloride-anion"))
    model = _resolve_model(parameters)

    assert model.topology == "ionic-ssm-ds"
    assert model.charge_number == (0, 1, -1)
    assert model.molar_mass_kg_per_mol == (0.0180153, None, None)
    assert model.k_ij == (0.0, -0.3, -0.3, -0.3, 0.0, 0.8, -0.3, 0.8, 0.0)
    assert model.l_ij == (None,) * 9
    assert model.k_hb_ij == (None,) * 9
    assert model.ionic_region_relative_permittivity == 8.0
    assert model.dielectric_ion_suppression_coefficient == 7.01
