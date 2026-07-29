import pytest

from epcsaft import Mixture, EosError, _ParameterBundle, unit_registry


def _methane() -> Mixture:
    return Mixture(
        _ParameterBundle.from_catalog("gross-2001-methane-ethane", version=1).select(("methane",))
    )


def test_pressure_state_returns_one_certified_root_and_diagnostics() -> None:
    model = _methane()
    result = model.state(
        T=250.0 * unit_registry.kelvin,
        P=10.0 * unit_registry.bar,
        x=(1.0,),
    )

    assert result.pressure.to("bar").magnitude == pytest.approx(10.0, abs=1.0e-8)
    assert result.molar_density.to("mole / meter ** 3").magnitude > 0.0
    diagnostics = result.density_diagnostics
    assert diagnostics is not None
    assert diagnostics.stable is True
    assert diagnostics.bracket_lower < result.molar_density < diagnostics.bracket_upper
    assert abs(diagnostics.pressure_residual.to("pascal").magnitude) < 1.0e-3
    assert diagnostics.pressure_density_derivative.to("pascal * meter ** 3 / mole").magnitude > 0.0
    assert diagnostics.iterations > 0
    assert diagnostics.branch == "single"

    density = result.molar_density.to("mole / meter ** 3").magnitude
    step = density * 1.0e-4

    def pressure_at(value: float) -> float:
        return (
            model.state(
                T=250.0 * unit_registry.kelvin,
                rho=value * unit_registry.mole / unit_registry.meter**3,
                x=(1.0,),
            )
            .pressure.to("pascal")
            .magnitude
        )

    finite_difference = (pressure_at(density + step) - pressure_at(density - step)) / (2.0 * step)
    assert diagnostics.pressure_density_derivative.to(
        "pascal * meter ** 3 / mole"
    ).magnitude == pytest.approx(finite_difference, rel=2.0e-7)


def test_three_root_state_requires_an_explicit_branch() -> None:
    model = _methane()
    state = {
        "T": 150.0 * unit_registry.kelvin,
        "P": 10.0 * unit_registry.bar,
        "x": (1.0,),
    }
    with pytest.raises(EosError):
        model.state(**state)

    vapor = model.state(**state, phase="vapor")
    liquid = model.state(**state, phase="liquid")

    assert vapor.molar_density < liquid.molar_density
    assert vapor.density_diagnostics is not None
    assert liquid.density_diagnostics is not None
    assert vapor.density_diagnostics.branch == "vapor"
    assert liquid.density_diagnostics.branch == "liquid"
    assert vapor.density_diagnostics.stable
    assert liquid.density_diagnostics.stable


def test_phase_request_rejects_a_single_root() -> None:
    with pytest.raises(EosError, match="requested density branch is unavailable"):
        _methane().state(
            T=250.0 * unit_registry.kelvin,
            P=10.0 * unit_registry.bar,
            x=(1.0,),
            phase="liquid",
        )


def test_root_topology_loss_is_reported_without_branch_relabelling() -> None:
    model = _methane()
    two_root_state = {
        "T": 190.0 * unit_registry.kelvin,
        "P": 45.0 * unit_registry.bar,
        "x": (1.0,),
    }
    vapor = model.state(**two_root_state, phase="vapor")
    liquid = model.state(**two_root_state, phase="liquid")

    assert vapor.density_diagnostics is not None
    assert liquid.density_diagnostics is not None
    assert vapor.density_diagnostics.branch == "vapor"
    assert liquid.density_diagnostics.branch == "liquid"
    assert vapor.molar_density < liquid.molar_density

    one_root_state = {
        "T": 191.0 * unit_registry.kelvin,
        "P": 46.0 * unit_registry.bar,
        "x": (1.0,),
    }
    result = model.state(**one_root_state)
    assert result.density_diagnostics is not None
    assert result.density_diagnostics.branch == "single"
    assert result.density_diagnostics.stable
    assert abs(
        result.density_diagnostics.pressure_residual.to("pascal").magnitude
    ) < 1.0e-3
    for requested_phase in ("vapor", "liquid"):
        with pytest.raises(EosError, match="requested density branch is unavailable"):
            model.state(**one_root_state, phase=requested_phase)


def test_density_and_pressure_inputs_are_mutually_exclusive() -> None:
    model = _methane()
    common = {
        "T": 250.0 * unit_registry.kelvin,
        "x": (1.0,),
    }
    density = 1000.0 * unit_registry.mole / unit_registry.meter**3
    pressure = 10.0 * unit_registry.bar

    with pytest.raises(EosError):
        model.state(**common)
    with pytest.raises(EosError):
        model.state(**common, rho=density, P=pressure)
    with pytest.raises(EosError):
        model.state(**common, rho=density, phase="vapor")
    with pytest.raises(EosError):
        model.state(**common, P=pressure, phase="supercritical")
    with pytest.raises(EosError):
        model.state(**common, P=1.0 * unit_registry.meter)
    with pytest.raises(EosError):
        model.state(**common, P=-1.0 * unit_registry.bar)

    explicit = model.state(**common, rho=density)
    assert explicit.molar_density == density
    assert explicit.density_diagnostics is None


def test_pressure_closure_runs_the_complete_electrolyte_eos() -> None:
    model = Mixture(
        _ParameterBundle.from_catalog("figiel-2025-reference-electrolytes", version=1).select(
            ("water", "sodium-cation", "chloride-anion")
        )
    )
    result = model.state(
        T=298.15 * unit_registry.kelvin,
        P=1.0 * unit_registry.bar,
        x=(0.9, 0.05, 0.05),
        phase="liquid",
    )

    assert result.pressure.to("bar").magnitude == pytest.approx(1.0, abs=1.0e-9)
    assert result.association != 0.0
    assert result.debye_huckel != 0.0
    assert result.born != 0.0
    assert result.density_diagnostics is not None
    assert result.density_diagnostics.branch == "liquid"
    assert result.density_diagnostics.stable


def test_ethanol_nabr_density_matches_held_measurement() -> None:
    model = Mixture(
        _ParameterBundle.from_catalog("figiel-2025-reference-electrolytes", version=1).select(
            ("ethanol", "sodium-cation", "bromide-anion")
        )
    )
    ethanol_moles = 1.0 / 0.046069
    salt_moles = 0.1422
    total_moles = ethanol_moles + 2.0 * salt_moles
    mole_fractions = (
        ethanol_moles / total_moles,
        salt_moles / total_moles,
        salt_moles / total_moles,
    )
    result = model.state(
        T=298.15 * unit_registry.kelvin,
        P=1.0 * unit_registry.bar,
        x=mole_fractions,
        phase="liquid",
    )
    molar_masses = (0.046069, 0.02298, 0.0799)
    mass_density = result.molar_density.to("mole / meter ** 3").magnitude * sum(
        fraction * molar_mass
        for fraction, molar_mass in zip(mole_fractions, molar_masses, strict=True)
    )

    # Held et al. (2012), DOI 10.1016/j.ces.2011.09.040, Table 1:
    # ethanol + NaBr at 298.15 K, ambient pressure, m = 0.1422 mol/kg,
    # rho = 795.96 kg/m3. One bar represents ambient pressure; the 20 kg/m3
    # gate is a model-accuracy tolerance, not the reported instrument uncertainty.
    assert abs(mass_density - 795.96) <= 20.0
