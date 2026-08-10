from __future__ import annotations

import pytest

from epcsaft import ActiveParameterSet, Mixture, Parameters
from epcsaft import unit_registry as u
from epcsaft.equilibrium import (
    PhaseEquilibriumError,
    PhaseEquilibriumProblem,
    PhaseHint,
    phase_equilibrium,
    tp_flash,
)
from tests.eos._raw_cases import AQUEOUS_ELECTROLYTE, NEUTRAL_BINARY

_TEMPERATURE = 300.0 * u.kelvin
_PRESSURE = 1.0e5 * u.pascal
_FEED = (0.5, 0.5)
_TOTAL_AMOUNT_MOL = 2.0
_AMOUNT_FEED = tuple(_TOTAL_AMOUNT_MOL * value for value in _FEED)


def test_phase_equilibrium_discovery_and_role_constraint_preserve_held() -> None:
    """Falsify the claim that the renamed route preserves incumbent HELD.

    Oracle: the retained public ``tp_flash`` result at the same fixed T, P,
    normalized composition, and Parameters. The amount-feed result and its
    exact sensitivity matrices obey the declared extensive-coordinate scale.
    A matching role constraint must retain the same state; a different topology
    must fail closed. Failure means the adapter changed the incumbent HELD
    calculation, broke its sensitivity basis, or treated a declared topology as
    an initialization hint.
    """
    parameters = Parameters.from_mapping(NEUTRAL_BINARY)
    active = ActiveParameterSet(
        parameters,
        (
            "pair/ethane/methane/k_ij",
            "component/methane/dispersion_energy_over_k",
        ),
    )
    model = Mixture(parameters)
    retained = tp_flash(
        model,
        _TEMPERATURE,
        _PRESSURE,
        _FEED,
        active_parameters=active,
    )
    problem = PhaseEquilibriumProblem(_TEMPERATURE, _PRESSURE, _AMOUNT_FEED)
    discovered = phase_equilibrium(model, problem, active_parameters=active)

    assert discovered.diagnostics.route == "held"
    assert discovered.diagnostics.globality_claim == "not_established"
    assert discovered.diagnostics.hints_applied is False
    assert discovered.total_free_energy_over_rt == pytest.approx(
        _TOTAL_AMOUNT_MOL * retained.total_free_energy_over_rt,
        rel=2.0e-12,
        abs=2.0e-12,
    )
    assert len(discovered.phases) == len(retained.phases)
    for phase, expected in zip(discovered.phases, retained.phases, strict=True):
        assert phase.amount_mol == pytest.approx(
            _TOTAL_AMOUNT_MOL * expected.amount_mol,
            abs=2.0e-12,
        )
        assert phase.mole_fractions == pytest.approx(
            expected.mole_fractions, rel=2.0e-12, abs=2.0e-12
        )
        assert phase.volume_m3 == pytest.approx(
            _TOTAL_AMOUNT_MOL * expected.volume_m3,
            rel=2.0e-12,
            abs=2.0e-12,
        )
        assert phase.pressure_pa == pytest.approx(expected.pressure_pa, rel=2.0e-12, abs=2.0e-12)

    selected_sensitivity = discovered.sensitivity
    incumbent_sensitivity = retained.sensitivity
    assert selected_sensitivity is not None
    assert incumbent_sensitivity is not None
    assert selected_sensitivity.parameters == incumbent_sensitivity.parameters
    assert (
        selected_sensitivity.residual_parameter_jacobian
        == incumbent_sensitivity.residual_parameter_jacobian
    )
    for selected_row, incumbent_row in zip(
        selected_sensitivity.residual_state_jacobian,
        incumbent_sensitivity.residual_state_jacobian,
        strict=True,
    ):
        assert selected_row == pytest.approx(
            tuple(value / _TOTAL_AMOUNT_MOL for value in incumbent_row),
            rel=2.0e-12,
            abs=2.0e-12,
        )
    for selected_rows, incumbent_rows in (
        (
            selected_sensitivity.state_parameter_jacobian,
            incumbent_sensitivity.state_parameter_jacobian,
        ),
        (
            selected_sensitivity.output_parameter_jacobian,
            incumbent_sensitivity.output_parameter_jacobian,
        ),
    ):
        for selected_row, incumbent_row in zip(
            selected_rows,
            incumbent_rows,
            strict=True,
        ):
            assert selected_row == pytest.approx(
                tuple(_TOTAL_AMOUNT_MOL * value for value in incumbent_row),
                rel=2.0e-12,
                abs=2.0e-12,
            )

    roles = tuple(phase.role for phase in discovered.phases)
    constrained = phase_equilibrium(
        model,
        PhaseEquilibriumProblem(
            _TEMPERATURE,
            _PRESSURE,
            _AMOUNT_FEED,
            mode="constrained",
            phase_roles=roles,
        ),
        active_parameters=active,
    )
    assert constrained.phases == discovered.phases
    assert constrained.diagnostics.requested_roles == roles
    assert constrained.diagnostics.returned_roles == roles

    hints = (
        PhaseHint(
            roles[0],
            tuple(reversed(discovered.phases[0].mole_fractions)),
            1.5 * discovered.phases[0].molar_density_mol_m3,
        ),
        PhaseHint(
            "liquid" if roles[0] == "vapor" else "vapor",
            discovered.phases[0].mole_fractions,
            0.5 * discovered.phases[0].molar_density_mol_m3,
        ),
    )
    hinted_results = tuple(
        phase_equilibrium(
            model,
            PhaseEquilibriumProblem(
                _TEMPERATURE,
                _PRESSURE,
                _AMOUNT_FEED,
                hints=ordered_hints,
            ),
            active_parameters=active,
        )
        for ordered_hints in (hints, tuple(reversed(hints)))
    )
    for hinted in hinted_results:
        assert hinted.phases == discovered.phases
        assert hinted.total_free_energy_over_rt == discovered.total_free_energy_over_rt
        assert hinted.sensitivity == discovered.sensitivity
        assert hinted.diagnostics.hint_count == 2
        assert hinted.diagnostics.hints_applied is False
    assert hinted_results[0].diagnostics == hinted_results[1].diagnostics

    rejected_roles = ("liquid" if roles[0] == "vapor" else "vapor",)
    with pytest.raises(PhaseEquilibriumError) as rejected:
        phase_equilibrium(
            model,
            PhaseEquilibriumProblem(
                _TEMPERATURE,
                _PRESSURE,
                _AMOUNT_FEED,
                mode="constrained",
                phase_roles=rejected_roles,
                hints=tuple(reversed(hints)),
            ),
            active_parameters=active,
        )
    assert rejected.value.diagnostics.failure_stage == "phase_selection"
    assert rejected.value.best_candidate == discovered.phases


def test_phase_role_classification_gap_fails_closed() -> None:
    """Falsify admission of a phase without a bounded physical role.

    Oracle: retained HELD at this explicit selected ``k_ij`` coordinate returns
    one mechanically admitted phase whose independently recomputed
    compressibility factor lies strictly inside the declared 0.25--0.45 role
    gap. The candidate must return no phase topology and mark physical status
    not adjudicated. Failure means a provisional liquid/vapor label escaped its
    stated numerical evidence boundary.
    """
    parameters = Parameters.from_mapping(NEUTRAL_BINARY)
    active = ActiveParameterSet(
        parameters,
        ("pair/ethane/methane/k_ij",),
        values=(-0.8,),
    )
    model = Mixture(parameters)
    retained = tp_flash(
        model,
        320.0 * u.kelvin,
        20.0e6 * u.pascal,
        _FEED,
        active_parameters=active,
    )
    assert len(retained.phases) == 1
    phase = retained.phases[0]
    compressibility = (
        20.0e6
        * phase.volume_m3
        / (phase.amount_mol * 8.31446261815324 * 320.0)
    )
    assert 0.25 < compressibility < 0.45

    with pytest.raises(PhaseEquilibriumError) as rejected:
        phase_equilibrium(
            model,
            PhaseEquilibriumProblem(
                320.0 * u.kelvin,
                20.0e6 * u.pascal,
                _FEED,
            ),
            active_parameters=active,
        )
    assert rejected.value.diagnostics.failure_reason == (
        "phase_role_classification_inconclusive"
    )
    assert rejected.value.diagnostics.physical_status == "not_adjudicated"
    assert rejected.value.diagnostics.returned_roles == ()
    assert rejected.value.best_candidate == ()


def test_charged_native_rejection_retains_held2_route_identity() -> None:
    """Falsify loss of the selected charged formulation on native failure.

    Oracle: a non-electroneutral feed is rejected by the retained electrolyte
    controller. The candidate failure must remain attributed to HELD2 even
    though the protected ``tp_flash`` failure payload predates typed route
    preservation for pre-result exceptions.
    """
    model = Mixture(Parameters.from_mapping(AQUEOUS_ELECTROLYTE))
    with pytest.raises(PhaseEquilibriumError) as rejected:
        phase_equilibrium(
            model,
            PhaseEquilibriumProblem(
                298.15 * u.kelvin,
                1.0e5 * u.pascal,
                (0.97, 0.02, 0.01),
            ),
        )
    assert rejected.value.diagnostics.route == "held2"
    assert rejected.value.diagnostics.failure_stage == "held2"
