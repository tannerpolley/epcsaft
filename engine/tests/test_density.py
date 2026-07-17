import pytest

from epcsaft import EPCSAFT, EosError, ParameterBundle, unit_registry


def _methane() -> EPCSAFT:
    return EPCSAFT(
        ParameterBundle.from_catalog(
            "gross-2001-methane-ethane", version=1
        ).select(("methane",))
    )


def test_pressure_state_returns_one_certified_root_and_diagnostics() -> None:
    model = _methane()
    result = model.evaluate(
        temperature=250.0 * unit_registry.kelvin,
        pressure=10.0 * unit_registry.bar,
        mole_fractions=(1.0,),
    )

    assert result.pressure.to("bar").magnitude == pytest.approx(10.0, abs=1.0e-8)
    assert result.molar_density.to("mole / meter ** 3").magnitude > 0.0
    diagnostics = result.density_diagnostics
    assert diagnostics is not None
    assert diagnostics.stable is True
    assert diagnostics.bracket_lower < result.molar_density < diagnostics.bracket_upper
    assert abs(diagnostics.pressure_residual.to("pascal").magnitude) < 1.0e-3
    assert diagnostics.pressure_density_derivative.to(
        "pascal * meter ** 3 / mole"
    ).magnitude > 0.0
    assert diagnostics.iterations > 0
    assert diagnostics.branch == "single"
    assert diagnostics.rejection_reason is None

    density = result.molar_density.to("mole / meter ** 3").magnitude
    step = density * 1.0e-4

    def pressure_at(value: float) -> float:
        return model.evaluate(
            temperature=250.0 * unit_registry.kelvin,
            molar_density=value * unit_registry.mole / unit_registry.meter**3,
            mole_fractions=(1.0,),
        ).pressure.to("pascal").magnitude

    finite_difference = (pressure_at(density + step) - pressure_at(density - step)) / (
        2.0 * step
    )
    assert diagnostics.pressure_density_derivative.to(
        "pascal * meter ** 3 / mole"
    ).magnitude == pytest.approx(finite_difference, rel=2.0e-7)


def test_three_root_state_requires_an_explicit_branch() -> None:
    model = _methane()
    state = {
        "temperature": 150.0 * unit_registry.kelvin,
        "pressure": 10.0 * unit_registry.bar,
        "mole_fractions": (1.0,),
    }
    with pytest.raises(EosError):
        model.evaluate(**state)

    vapor = model.evaluate(**state, phase="vapor")
    liquid = model.evaluate(**state, phase="liquid")

    assert vapor.molar_density < liquid.molar_density
    assert vapor.density_diagnostics is not None
    assert liquid.density_diagnostics is not None
    assert vapor.density_diagnostics.branch == "vapor"
    assert liquid.density_diagnostics.branch == "liquid"
    assert vapor.density_diagnostics.stable
    assert liquid.density_diagnostics.stable


def test_phase_hint_does_not_relabel_a_single_root() -> None:
    result = _methane().evaluate(
        temperature=250.0 * unit_registry.kelvin,
        pressure=10.0 * unit_registry.bar,
        mole_fractions=(1.0,),
        phase="liquid",
    )

    assert result.density_diagnostics is not None
    assert result.density_diagnostics.branch == "single"


def test_density_and_pressure_inputs_are_mutually_exclusive() -> None:
    model = _methane()
    common = {
        "temperature": 250.0 * unit_registry.kelvin,
        "mole_fractions": (1.0,),
    }
    density = 1000.0 * unit_registry.mole / unit_registry.meter**3
    pressure = 10.0 * unit_registry.bar

    with pytest.raises(EosError):
        model.evaluate(**common)
    with pytest.raises(EosError):
        model.evaluate(**common, molar_density=density, pressure=pressure)
    with pytest.raises(EosError):
        model.evaluate(**common, molar_density=density, phase="vapor")

    explicit = model.evaluate(**common, molar_density=density)
    assert explicit.molar_density == density
    assert explicit.density_diagnostics is None


def test_pressure_closure_runs_the_complete_electrolyte_eos() -> None:
    model = EPCSAFT(
        ParameterBundle.from_catalog(
            "figiel-2025-water-sodium-chloride", version=1
        ).select(("water", "sodium-cation", "chloride-anion"))
    )
    result = model.evaluate(
        temperature=298.15 * unit_registry.kelvin,
        pressure=1.0 * unit_registry.bar,
        mole_fractions=(0.9, 0.05, 0.05),
        phase="liquid",
    )

    assert result.pressure.to("bar").magnitude == pytest.approx(1.0, abs=1.0e-9)
    assert result.association != 0.0
    assert result.debye_huckel != 0.0
    assert result.born_ssm_ds != 0.0
    assert result.density_diagnostics is not None
    assert result.density_diagnostics.branch == "liquid"
    assert result.density_diagnostics.stable
