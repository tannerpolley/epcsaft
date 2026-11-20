import json

import pytest
from epcsaft import Mixture, ParameterError, Parameters, _core, equilibrium as eq
from tests.packets import parameters_json

TEMPERATURE_K, PRESSURE_PA = 298.15, 100000.0
WATER, SODIUM, MAGNESIUM, CHLORIDE, METHANOL = (
    "water", "sodium-cation", "magnesium-cation", "chloride-anion", "methanol",
)


def brine_model(*components: str) -> Mixture:
    path = parameters_json("figiel-2025-chloride-activity-refit", 2)
    return Mixture(Parameters.from_json(path, components=components))


def test_salt_recipe_molalities() -> None:
    water, sodium, magnesium, chloride = components = (WATER, SODIUM, MAGNESIUM, CHLORIDE)
    model = brine_model(*components)
    salts = {"NaCl": ({sodium: 1, chloride: 1}, 0.5), "MgCl2": ({magnesium: 1, chloride: 2}, 0.25)}
    solvent_mass = 2.0
    feed = eq.Molalities.from_salts(model, salts, solvent=water, solvent_mass_kg=solvent_mass)
    direct = eq.Amounts(dict(zip(components, (solvent_mass / model.parameters.molar_masses[0], 1, 0.5, 2), strict=True)))

    def solve(feed: eq.Feed) -> eq.EquilibriumResult:
        result = eq.solve_equilibrium(model, eq.Problem(
            phases=[eq.Phase("aqueous", kind="liquid")], T=TEMPERATURE_K, P=PRESSURE_PA, feed=feed,
        ))
        assert result.success, result.message
        return result

    mixed, explicit = solve(feed), solve(direct)
    assert mixed.amounts == explicit.amounts
    assert mixed.molar_densities == explicit.molar_densities
    assert mixed.mole_fractions == explicit.mole_fractions
    amounts = [mixed.amounts[0] * x for x in mixed.mole_fractions]
    mass = amounts[0] * model.parameters.molar_masses[0]
    charge = amounts[1] + 2 * amounts[2] - amounts[3]
    # Roundoff from reconstructing inventories from the common solved phase.
    assert amounts == pytest.approx(list(direct.values.values()), rel=1e-10)
    assert mass == pytest.approx(solvent_mass, rel=1e-10)
    assert charge == pytest.approx(0.0, abs=1e-10)
    with pytest.raises(ParameterError, match="MgCl formula unit must be electroneutral"):
        eq.Molalities.from_salts(model, {"MgCl": ({magnesium: 1, chloride: 1}, 0.25)},
                                 solvent=water, solvent_mass_kg=solvent_mass)


@pytest.mark.parametrize("solute_molalities", [None, (0.0, 0.5, 0.25, 1.0)])
def test_activity_rejects_mismatched_reference_ions(solute_molalities: tuple[float, ...] | None) -> None:
    model = brine_model(WATER, SODIUM, MAGNESIUM, CHLORIDE)
    solvent = (1.0, 0.0, 0.0, 0.0)
    sodium = model.reference(TEMPERATURE_K, PRESSURE_PA, solvent_fractions=solvent,
                             ions=(SODIUM, CHLORIDE), stoichiometry=(1, 1))
    # The public wrapper normally pairs its own ion indices and reference.
    with pytest.raises(ParameterError, match="ordered ion indices"):
        _core.electrolyte_activity(model._model, TEMPERATURE_K, PRESSURE_PA, solvent, (2, 3), (1, 2),
                                   0.25, sodium._native, solute_molalities)


def test_full_solute_activity_requires_exactly_pure_solvent() -> None:
    path = parameters_json("figiel-2025-chloride-activity-refit", 2)
    document = json.loads(path.read_text())
    document["topology"].setdefault("no_cross_association", []).append({
        "component_id_a": WATER, "component_id_b": METHANOL,
        "provenance": {"locator": "Manufactured solvent admission test"},
    })
    model = Mixture(Parameters.from_mapping(document, components=(WATER, SODIUM, CHLORIDE, METHANOL)))
    reference = model.reference(TEMPERATURE_K, PRESSURE_PA, solvent_fractions=(1.0, 0.0, 0.0, 5e-13),
                                ions=(SODIUM, CHLORIDE), stoichiometry=(1, 1))
    with pytest.raises(ParameterError, match="one pure solvent"):
        reference.activity(0.5, solute_molalities=(0.0, 0.5, 0.5, 0.0))


def test_full_solute_activity_uses_cation_molality_in_either_ion_order() -> None:
    model = brine_model(WATER, SODIUM, MAGNESIUM, CHLORIDE)
    solvent = (1.0, 0.0, 0.0, 0.0)
    molalities = (0.0, 0.5, 0.25, 1.0)
    forward = model.reference(TEMPERATURE_K, PRESSURE_PA, solvent_fractions=solvent,
                              ions=(SODIUM, CHLORIDE), stoichiometry=(1, 1))
    reverse = model.reference(TEMPERATURE_K, PRESSURE_PA, solvent_fractions=solvent,
                              ions=(CHLORIDE, SODIUM), stoichiometry=(1, 1))
    expected = forward.activity(0.5, solute_molalities=molalities)
    actual = reverse.activity(0.5, solute_molalities=molalities)
    assert actual.mean_ionic_molality_scale == pytest.approx(expected.mean_ionic_molality_scale, rel=1e-10)
    assert actual.reference_solvent_fractions == expected.reference_solvent_fractions == list(solvent)
    with pytest.raises(ParameterError, match="match the selected formula cation molality"):
        reverse.activity(1.0, solute_molalities=molalities)
