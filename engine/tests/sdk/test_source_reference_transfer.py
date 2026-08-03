from __future__ import annotations

from dataclasses import replace
import math

import pytest

from epcsaft import (
    EosError,
    Mixture,
    SourceReferenceDeclaration,
    source_reference_transfer,
    unit_registry,
)
from tests.sdk.test_neutral_reference import (
    _CHARGES,
    _COMPONENTS,
    _generic_electrolyte,
    _matrix_rank,
)


def _declaration(
    component_ids: tuple[str, ...] = _COMPONENTS,
    *,
    required_derivatives: tuple[str, ...] = (),
) -> SourceReferenceDeclaration:
    return SourceReferenceDeclaration(
        reference_state_id="declared-solvent-molality-infinite-dilution-v1",
        component_ids=component_ids,
        solvent_mole_fractions=tuple(
            1.0 if component_id == "neutral-alpha" else 0.0
            for component_id in component_ids
        ),
        activity_convention_id="molality-infinite-dilution-v1",
        standard_molality_mol_per_kg=1.0,
        reference_pressure_pa=100000.0,
        required_derivatives=required_derivatives,
    )


@pytest.mark.parametrize("variant", ("ordinary", "associating", "polar"))
def test_source_reference_transfer_uses_declared_solvent_for_every_model(
    variant: str,
) -> None:
    model = Mixture.from_dictionary(_generic_electrolyte(variant=variant))

    result = source_reference_transfer(
        model,
        _declaration(),
        T=298.15 * unit_registry.kelvin,
        P=75000.0 * unit_registry.pascal,
    )
    at_source_reference_pressure = source_reference_transfer(
        model,
        _declaration(),
        T=298.15 * unit_registry.kelvin,
        P=100000.0 * unit_registry.pascal,
    )

    assert result.component_ids == _COMPONENTS
    assert result.status == "ok"
    assert result.domain_status == "admitted"
    assert result.convergence_status == "converged"
    assert result.reference_state_id == _declaration().reference_state_id
    assert result.activity_convention_id == "molality-infinite-dilution-v1"
    assert result.source_composition == (1.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    assert len(result.neutral_basis) == len(_COMPONENTS) - 1
    assert _matrix_rank(result.neutral_basis) == len(result.neutral_basis)
    charges = tuple(_CHARGES[component_id] for component_id in _COMPONENTS)
    assert all(
        math.fsum(value * charge for value, charge in zip(row, charges, strict=True))
        == 0.0
        for row in result.neutral_basis
    )
    assert all(math.isfinite(value) for value in result.transfer_log_contractions)
    assert result.activity_scale_log_contractions[0] == 0.0
    assert result.activity_scale_log_contractions[1] == pytest.approx(
        math.log(0.0180153), abs=2.0e-12
    )
    assert result.pressure_pa == 75000.0
    assert result.source_reference_pressure_pa == 100000.0
    assert any(
        abs(trial - source) > 1.0e-8
        for trial, source in zip(
            result.transfer_log_contractions,
            at_source_reference_pressure.transfer_log_contractions,
            strict=True,
        )
    )
    assert result.source_pressure_interval_pa == (50000.0, 200000.0)
    assert result.parameter_fingerprint == model.parameter_fingerprint
    assert result.topology_fingerprint.startswith("sha256:")
    assert result.component_order_fingerprint.startswith("sha256:")
    assert result.reference_state_fingerprint.startswith("sha256:")
    assert result.domain_fingerprint.startswith("sha256:")
    assert result.artifact_fingerprint.startswith("sha256:")


def test_source_reference_transfer_is_component_permutation_safe() -> None:
    permutation = (
        "anion-divalent",
        "neutral-beta",
        "cation-monovalent",
        "cation-secondary",
        "neutral-alpha",
        "anion-monovalent",
    )
    canonical = source_reference_transfer(
        Mixture.from_dictionary(_generic_electrolyte()),
        _declaration(),
        T=298.15 * unit_registry.kelvin,
        P=100000.0 * unit_registry.pascal,
    )
    permuted = source_reference_transfer(
        Mixture.from_dictionary(_generic_electrolyte(permutation)),
        _declaration(permutation),
        T=298.15 * unit_registry.kelvin,
        P=100000.0 * unit_registry.pascal,
    )
    inverse = tuple(permutation.index(component_id) for component_id in _COMPONENTS)

    assert tuple(
        tuple(row[index] for index in inverse) for row in permuted.neutral_basis
    ) == canonical.neutral_basis
    assert permuted.transfer_log_contractions == pytest.approx(
        canonical.transfer_log_contractions, abs=6.0e-11
    )


def test_source_reference_transfer_fails_closed() -> None:
    model = Mixture.from_dictionary(_generic_electrolyte())

    for pressure_pa in (50000.0, 200000.0):
        assert source_reference_transfer(
            model,
            _declaration(),
            T=298.15 * unit_registry.kelvin,
            P=pressure_pa * unit_registry.pascal,
        ).pressure_pa == pressure_pa
    for pressure_pa in (49999.99, 200000.01):
        with pytest.raises(EosError, match="provider-domain-rejection"):
            source_reference_transfer(
                model,
                _declaration(),
                T=298.15 * unit_registry.kelvin,
                P=pressure_pa * unit_registry.pascal,
            )

    charged_source = replace(
        _declaration(),
        solvent_mole_fractions=(0.9, 0.0, 0.1, 0.0, 0.0, 0.0),
    )
    with pytest.raises(EosError, match="invalid-source-composition"):
        source_reference_transfer(
            model,
            charged_source,
            T=298.15 * unit_registry.kelvin,
            P=100000.0 * unit_registry.pascal,
        )
    with pytest.raises(EosError, match="derivative-unavailable"):
        source_reference_transfer(
            model,
            _declaration(required_derivatives=("pressure",)),
            T=298.15 * unit_registry.kelvin,
            P=100000.0 * unit_registry.pascal,
        )
    with pytest.raises(EosError, match="nonfinite-source-reference-state"):
        source_reference_transfer(
            model,
            replace(_declaration(), reference_pressure_pa=math.nan),
            T=298.15 * unit_registry.kelvin,
            P=100000.0 * unit_registry.pascal,
        )
    with pytest.raises(EosError, match="derivative-unavailable"):
        source_reference_transfer(
            model,
            _declaration(required_derivatives=("unknown",)),
            T=298.15 * unit_registry.kelvin,
            P=100000.0 * unit_registry.pascal,
        )
    with pytest.raises(EosError, match="component-order-mismatch"):
        source_reference_transfer(
            model,
            replace(_declaration(), component_ids=tuple(reversed(_COMPONENTS))),
            T=298.15 * unit_registry.kelvin,
            P=100000.0 * unit_registry.pascal,
        )
    with pytest.raises(EosError, match="unsupported-source-reference-transfer"):
        source_reference_transfer(
            model,
            replace(_declaration(), activity_convention_id="unsupported"),
            T=298.15 * unit_registry.kelvin,
            P=100000.0 * unit_registry.pascal,
        )


def test_source_reference_declaration_detaches_mutable_sequences() -> None:
    component_ids = list(_COMPONENTS)
    solvent = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    declaration = SourceReferenceDeclaration(
        reference_state_id="declared-solvent-molality-infinite-dilution-v1",
        component_ids=component_ids,  # type: ignore[arg-type]
        solvent_mole_fractions=solvent,  # type: ignore[arg-type]
        activity_convention_id="molality-infinite-dilution-v1",
        standard_molality_mol_per_kg=1.0,
        reference_pressure_pa=100000.0,
        required_derivatives=[],  # type: ignore[arg-type]
    )
    component_ids[0] = "mutated"
    solvent[0] = 0.0

    result = source_reference_transfer(
        Mixture.from_dictionary(_generic_electrolyte()),
        declaration,
        T=298.15 * unit_registry.kelvin,
        P=75000.0 * unit_registry.pascal,
    )

    assert result.component_ids == _COMPONENTS
    assert result.source_composition[0] == 1.0
