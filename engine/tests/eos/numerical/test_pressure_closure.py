from __future__ import annotations

import pytest

from epcsaft import EosError, Mixture, Parameters, Pure
from epcsaft import unit_registry as u
from tests.eos._cases import NEUTRAL_BINARY, neutral_pure


@pytest.mark.parametrize("phase", ("vapor", "liquid"))
def test_fixed_pressure_reconstructs_selected_stable_branch(phase: str) -> None:
    model = Mixture(Parameters.from_mapping(NEUTRAL_BINARY))
    temperature = 200.0
    pressure = 1.0e5
    composition = (0.4, 0.6)
    solved = model.state(
        T=temperature * u.kelvin,
        P=pressure * u.pascal,
        x=composition,
        phase=phase,
    )
    diagnostics = solved.density_diagnostics
    assert diagnostics is not None
    assert diagnostics.stable is True
    assert diagnostics.branch == phase
    assert abs(solved.pressure.to("pascal").magnitude - pressure) / pressure <= 1.0e-8

    explicit = model.state(
        T=temperature * u.kelvin,
        rho=solved.molar_density,
        x=composition,
    )
    assert explicit.ares(basis="reduced") == pytest.approx(
        solved.ares(basis="reduced"), rel=2.0e-13, abs=2.0e-13
    )
    assert explicit.residual_chemical_potential_over_rt == pytest.approx(
        solved.residual_chemical_potential_over_rt,
        rel=2.0e-12,
        abs=2.0e-12,
    )


def test_fixed_pressure_unbracketed_and_domain_requests_fail_typed() -> None:
    model = Pure(Parameters.from_mapping(neutral_pure("propane")))

    with pytest.raises(EosError) as unbracketed:
        model.state(T=300.0 * u.kelvin, P=1.0e12 * u.pascal, phase="vapor")
    assert unbracketed.value.code == "evaluation-failed"

    with pytest.raises(EosError) as invalid_domain:
        model.state(T=300.0 * u.kelvin, P=1.0e5 * u.pascal, phase="solid")
    assert invalid_domain.value.code == "invalid-state"
