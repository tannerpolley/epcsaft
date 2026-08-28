from __future__ import annotations

import math

import pytest

from epcsaft import Activity, EosError, Mixture, Parameters, TransferFreeEnergy
from epcsaft import unit_registry as u
from tests.eos._cases import AQUEOUS_ELECTROLYTE, MULTISOLVENT_MULTIION


def test_reference_projection_activity_and_identity_transfer() -> None:
    model = Mixture(Parameters.from_mapping(AQUEOUS_ELECTROLYTE))
    reference = model.reference(
        T=298.15 * u.kelvin,
        P=1.0 * u.bar,
        solvent_x=(1.0, 0.0, 0.0),
        phase="liquid",
        ion_component_ids=("sodium-cation", "chloride-anion"),
        formula_stoichiometry=(1, 1),
        reference_convention="molality-infinite-dilution",
    )
    state = reference.state(molality=0.1 * u.mole / u.kilogram)
    finite = model.state(
        T=reference.temperature,
        P=reference.pressure,
        x=state.mole_fractions,
        phase="liquid",
    )

    assert isinstance(state.activity, Activity)
    assert finite.fugacity is not None
    projected = tuple(
        finite.fugacity.ln_coefficient[index] - infinite
        for index, infinite in zip(
            (1, 2),
            reference.infinite_dilution_log_fugacity_coefficient,
            strict=True,
        )
    )
    assert state.activity.single_ion == pytest.approx(
        tuple(math.exp(value) for value in projected),
        rel=2.0e-12,
        abs=2.0e-12,
    )
    assert state.activity.mean_mole_fraction == pytest.approx(
        math.exp(math.fsum(projected) / 2.0), rel=2.0e-12
    )
    assert state.activity.osmotic_coefficient is not None
    assert math.isfinite(state.activity.osmotic_coefficient)

    identity = reference.transfer_to(reference)
    assert isinstance(identity, TransferFreeEnergy)
    assert tuple(value.to("joule / mole").magnitude for value in identity.values) == (
        0.0,
        0.0,
    )


def test_reference_rejects_missing_ambiguous_and_incompatible_declarations() -> None:
    model = Mixture(Parameters.from_mapping(MULTISOLVENT_MULTIION))
    arguments = {
        "T": 298.15 * u.kelvin,
        "P": 1.0 * u.bar,
        "solvent_x": (0.6, 0.4, 0.0, 0.0, 0.0, 0.0),
        "phase": "liquid",
        "formula_stoichiometry": (1, 1),
        "reference_convention": "molality-infinite-dilution",
    }

    with pytest.raises(EosError) as missing:
        model.reference(
            **arguments,
            ion_component_ids=("missing-cation", "anion"),
        )
    assert missing.value.code == "invalid-reference-state"

    with pytest.raises(EosError) as ambiguous_roles:
        model.reference(
            **arguments,
            ion_component_ids=("anion", "cation"),
        )
    assert ambiguous_roles.value.code == "invalid-reference-state"

    monovalent = model.reference(
        **arguments,
        ion_component_ids=("cation", "anion"),
    )
    divalent = model.reference(
        **(arguments | {"formula_stoichiometry": (1, 2)}),
        ion_component_ids=("divalent-cation", "iodide-anion"),
    )
    with pytest.raises(EosError) as incompatible:
        monovalent.transfer_to(divalent)
    assert incompatible.value.code == "invalid-reference-state"
