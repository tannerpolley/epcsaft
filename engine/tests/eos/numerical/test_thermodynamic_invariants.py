from __future__ import annotations

import math

import pytest
from tests.eos._raw_cases import (
    AQUEOUS_ELECTROLYTE,
    ASSOCIATING_PURE,
    DIPOLAR_PERMITTIVITY_ELECTROLYTE,
    NEUTRAL_BINARY,
    neutral_binary,
    neutral_pure,
)

from epcsaft import Activity, Mixture, Parameters, TransferFreeEnergy
from epcsaft import unit_registry as u

AVOGADRO = 6.02214076e23
BOLTZMANN = 1.380649e-23
ELEMENTARY_CHARGE = 1.602176634e-19
GAS_CONSTANT = 8.31446261815324
VACUUM_PERMITTIVITY = 8.8541878128e-12


def _first_derivative(function, value: float, step: float) -> float:
    return (
        -function(value + 2.0 * step)
        + 8.0 * function(value + step)
        - 8.0 * function(value - step)
        + function(value - 2.0 * step)
    ) / (12.0 * step)


def _second_derivative(function, value: float, step: float) -> float:
    return (
        -function(value + 2.0 * step)
        + 16.0 * function(value + step)
        - 30.0 * function(value)
        + 16.0 * function(value - step)
        - function(value - 2.0 * step)
    ) / (12.0 * step * step)


def _state(mapping, *, x, density=1000.0, temperature=298.15):
    return Mixture(Parameters.from_dictionary(mapping)).state(
        T=temperature * u.kelvin,
        rho=density * u.mole / u.meter**3,
        x=x,
    )


@pytest.mark.parametrize(
    ("mapping", "composition"),
    (
        (NEUTRAL_BINARY, (0.4, 0.6)),
        (AQUEOUS_ELECTROLYTE, (0.98, 0.01, 0.01)),
    ),
    ids=("neutral-binary", "electrolyte"),
)
def test_residual_helmholtz__equals_sum_of_active_contributions(mapping, composition) -> None:
    state = _state(mapping, x=composition)

    assert state.ares(basis="reduced") == pytest.approx(
        state.hard_chain
        + state.dispersion
        + state.association
        + state.ares_polar
        + state.debye_huckel
        + state.born,
        rel=2.0e-14,
        abs=2.0e-14,
    )
    assert all(
        math.isfinite(value)
        for value in (
            state.compressibility_factor,
            state.pressure.to("pascal").magnitude,
            *state.residual_chemical_potential_over_rt,
        )
    )


def test_neutral_model__inactive_electrolyte_terms_are_exactly_zero() -> None:
    state = _state(NEUTRAL_BINARY, x=(0.4, 0.6))

    assert state.association == 0.0
    assert state.ares_polar == 0.0
    assert state.debye_huckel == 0.0
    assert state.born == 0.0
    assert state.bulk_relative_permittivity is None


def test_component_permutation__preserves_scalar_thermodynamics() -> None:
    forward = _state(NEUTRAL_BINARY, x=(0.4, 0.6))
    reverse = _state(neutral_binary(("ethane", "methane")), x=(0.6, 0.4))

    assert reverse.ares(basis="reduced") == pytest.approx(forward.ares(basis="reduced"))
    assert reverse.pressure.to("pascal").magnitude == pytest.approx(
        forward.pressure.to("pascal").magnitude
    )


def test_pressure_identity__matches_compressibility_definition() -> None:
    model = Mixture(Parameters.from_dictionary(NEUTRAL_BINARY))
    state = model.state(
        T=298.15 * u.kelvin,
        rho=40.0 * u.mole / u.meter**3,
        x=(0.4, 0.6),
    )
    expected = state.compressibility_factor * 40.0 * 8.31446261815324 * 298.15

    assert state.pressure.to("pascal").magnitude == pytest.approx(expected, rel=2.0e-14)


@pytest.mark.parametrize(
    ("temperature", "density", "composition", "hard_chain", "dispersion"),
    (
        (200.0, 1000.0, (1.0, 0.0), 0.06290111429923406, -0.16906790398822233),
        (250.0, 8000.0, (0.4, 0.6), 0.8611345413201740, -1.9282055648091894),
    ),
    ids=("pure-limit", "binary"),
)
def test_neutral_contributions__match_independent_gross_anchors(
    temperature, density, composition, hard_chain, dispersion
) -> None:
    state = Mixture.from_dictionary(NEUTRAL_BINARY).state(
        T=temperature * u.kelvin,
        rho=density * u.mole / u.meter**3,
        x=composition,
    )

    assert state.hard_chain == pytest.approx(hard_chain, rel=5.0e-12, abs=5.0e-13)
    assert state.dispersion == pytest.approx(dispersion, rel=5.0e-12, abs=5.0e-13)


def test_association__matches_independent_pure_2b_closed_form() -> None:
    temperature = 300.0
    density = 10000.0
    sigma = 3.0
    diameter = sigma * (1.0 - 0.12 * math.exp(-3.0 * 200.0 / temperature))
    number_density = density * AVOGADRO * 1.0e-30
    eta = math.pi * number_density * diameter**3 / 6.0
    contact = (
        1.0 / (1.0 - eta)
        + 1.5 * diameter * number_density * math.pi * diameter**2 / 6.0 / (1.0 - eta) ** 2
        + 0.5 * diameter**2 * (number_density * math.pi * diameter**2 / 6.0) ** 2 / (1.0 - eta) ** 3
    )
    delta = contact * math.expm1(2400.0 / temperature) * sigma**3 * 0.03
    site_fraction = 2.0 / (1.0 + math.sqrt(1.0 + 4.0 * number_density * delta))
    expected = 2.0 * (math.log(site_fraction) - 0.5 * site_fraction + 0.5)

    state = _state(ASSOCIATING_PURE, x=(1.0,), density=density, temperature=temperature)

    assert state.association == pytest.approx(expected, rel=2.0e-11, abs=2.0e-12)


def test_debye_huckel__matches_independent_si_oracle() -> None:
    temperature = 298.15
    density = 40.0
    fractions = (0.9, 0.05, 0.05)
    charges = (0, 1, -1)
    diameters = (
        2.7927
        + 10.11 * math.exp(-0.01775 * temperature)
        - 1.417 * math.exp(-0.01146 * temperature),
        2.8232,
        2.7560,
    )
    permittivity = 78.09 / (1.0 + 7.01 * 0.1)
    charge_square_sum = math.fsum(
        fraction * charge * charge for fraction, charge in zip(fractions, charges, strict=True)
    )
    kappa = math.sqrt(
        density
        * AVOGADRO
        * ELEMENTARY_CHARGE**2
        * charge_square_sum
        / (VACUUM_PERMITTIVITY * permittivity * BOLTZMANN * temperature)
    )
    weighted_finite_size = 0.0
    for fraction, charge, diameter in zip(fractions, charges, diameters, strict=True):
        if charge == 0:
            continue
        reduced = kappa * 0.88 * diameter * 1.0e-10
        chi = (
            3.0
            / reduced**3
            * (1.5 + math.log1p(reduced) - 2.0 * (1.0 + reduced) + 0.5 * (1.0 + reduced) ** 2)
        )
        weighted_finite_size += fraction * charge * charge * chi
    expected = (
        -kappa
        * ELEMENTARY_CHARGE**2
        * weighted_finite_size
        / (12.0 * math.pi * VACUUM_PERMITTIVITY * permittivity * BOLTZMANN * temperature)
    )

    state = _state(AQUEOUS_ELECTROLYTE, x=fractions, density=density)

    assert state.debye_huckel == pytest.approx(expected, rel=2.0e-11, abs=2.0e-13)


def test_corrected_born__matches_independent_two_region_oracle() -> None:
    temperature = 298.15
    fractions = (0.999998, 0.000001, 0.000001)
    charges = (0, 1, -1)
    diameters = (0.0, 3.445, 4.100)
    solvation = (1.5, 1.0, 1.0)
    mixed_solvation = math.fsum(
        fraction * factor for fraction, factor in zip(fractions, solvation, strict=True)
    )
    bulk_permittivity = 78.09 / (1.0 + 7.01 * math.fsum(fractions[1:]))
    radial_sum = 0.0
    for fraction, charge, diameter in zip(fractions, charges, diameters, strict=True):
        if charge == 0:
            continue
        born_diameter = diameter * 1.0e-10
        shell_diameter = born_diameter * (1.0 + (mixed_solvation - 1.0) / abs(charge))
        radial_sum += (
            fraction
            * charge
            * charge
            * (
                (1.0 - 1.0 / bulk_permittivity) / shell_diameter
                + (1.0 - 1.0 / 8.0) * (1.0 / born_diameter - 1.0 / shell_diameter)
            )
        )
    expected = (
        -(ELEMENTARY_CHARGE**2)
        * radial_sum
        / (4.0 * math.pi * VACUUM_PERMITTIVITY * BOLTZMANN * temperature)
    )

    state = _state(AQUEOUS_ELECTROLYTE, x=fractions, density=40000.0)

    assert state.born == pytest.approx(expected, rel=3.0e-12, abs=2.0e-12)


def test_permittivity__matches_independent_ion_fraction_and_dipolar_equations() -> None:
    fractions = (0.9, 0.05, 0.05)
    density = 55000.0
    ion_fraction_state = _state(AQUEOUS_ELECTROLYTE, x=fractions, density=density)
    assert ion_fraction_state.bulk_relative_permittivity == pytest.approx(
        78.09 / (1.0 + 7.01 * 0.1), rel=2.0e-12
    )

    number_density = density * AVOGADRO * 1.0e-30
    beta = 1.0 / (BOLTZMANN * 298.15)
    effective = 4.0 * math.pi / 9.0 * number_density * fractions[0] * (beta * 1.0e-19 + 3.0)
    integral = 1.0 + math.expm1(-effective)
    expected = 1.0 + 3.0 * effective * (
        1.0 + effective + ((17.0 / 16.0) * integral - 1.0) * effective**2
    )
    dipolar = _state(DIPOLAR_PERMITTIVITY_ELECTROLYTE, x=fractions, density=density)
    assert dipolar.bulk_relative_permittivity == pytest.approx(expected, rel=2.0e-12)


def test_density_closure__returns_certified_single_root_and_slope() -> None:
    model = Mixture.from_dictionary(neutral_pure("methane"))
    state = model.state(T=250.0 * u.kelvin, P=10.0 * u.bar, x=(1.0,))
    diagnostics = state.density_diagnostics

    assert diagnostics is not None
    assert diagnostics.stable is True
    assert diagnostics.branch == "single"
    assert abs(diagnostics.pressure_residual.to("pascal").magnitude) < 1.0e-3
    density = state.molar_density.to("mole / meter ** 3").magnitude
    step = density * 1.0e-4
    pressures = tuple(
        model.state(
            T=250.0 * u.kelvin,
            rho=value * u.mole / u.meter**3,
            x=(1.0,),
        )
        .pressure.to("pascal")
        .magnitude
        for value in (density - step, density + step)
    )
    slope = (pressures[1] - pressures[0]) / (2.0 * step)
    assert diagnostics.pressure_density_derivative.to(
        "pascal * meter ** 3 / mole"
    ).magnitude == pytest.approx(slope, rel=2.0e-7)


def test_residual_chemical_potentials__match_total_amount_derivatives() -> None:
    model = Mixture.from_dictionary(NEUTRAL_BINARY)
    temperature = 250.0
    density = 8000.0
    fractions = (0.4, 0.6)
    volume = 1.0 / density
    state = _state(NEUTRAL_BINARY, x=fractions, density=density, temperature=temperature)
    expected = []
    step = 2.0e-5
    for component in range(2):

        def total_helmholtz(delta: float) -> float:
            amounts = [*fractions]
            amounts[component] += delta
            total = math.fsum(amounts)
            perturbed = model.state(
                T=temperature * u.kelvin,
                rho=(total / volume) * u.mole / u.meter**3,
                x=tuple(value / total for value in amounts),
            )
            return total * perturbed.ares(basis="reduced")

        expected.append(_first_derivative(total_helmholtz, 0.0, step))

    assert state.residual_chemical_potential_over_rt == pytest.approx(
        expected, rel=2.0e-7, abs=2.0e-8
    )
    assert state.fugacity is not None
    assert state.fugacity.ln_coefficient == pytest.approx(
        tuple(value - math.log(state.compressibility_factor) for value in expected),
        rel=2.0e-7,
        abs=2.0e-8,
    )


def test_residual_caloric_properties__follow_helmholtz_identities() -> None:
    temperature = 250.0
    density = 1000.0
    fractions = (0.4, 0.6)
    state = _state(NEUTRAL_BINARY, x=fractions, density=density, temperature=temperature)

    def helmholtz_at(value: float) -> float:
        return _state(NEUTRAL_BINARY, x=fractions, density=density, temperature=value).ares(
            basis="reduced"
        )

    temperature_step = temperature * 2.0e-4
    a_t = _first_derivative(helmholtz_at, temperature, temperature_step)
    a_tt = _second_derivative(helmholtz_at, temperature, temperature_step)

    def z_at_temperature(value: float) -> float:
        return _state(
            NEUTRAL_BINARY, x=fractions, density=density, temperature=value
        ).compressibility_factor

    def z_at_density(value: float) -> float:
        return _state(
            NEUTRAL_BINARY, x=fractions, density=value, temperature=temperature
        ).compressibility_factor

    z = state.compressibility_factor
    z_t = _first_derivative(z_at_temperature, temperature, temperature_step)
    z_rho = _first_derivative(z_at_density, density, density * 2.0e-4)
    cv = -2.0 * temperature * a_t - temperature**2 * a_tt
    expected = {
        "ures": -temperature * a_t,
        "hres": -temperature * a_t + z - 1.0,
        "sres": -state.ares(basis="reduced") - temperature * a_t,
        "cvres": cv,
        "cpres": cv + (z + temperature * z_t) ** 2 / (z + density * z_rho) - 1.0,
    }

    for method, value in expected.items():
        assert getattr(state, method)(basis="reduced") == pytest.approx(
            value, rel=3.0e-5, abs=3.0e-5
        )


def test_electrolyte_reference__matches_public_projection_identities() -> None:
    model = Mixture.from_dictionary(AQUEOUS_ELECTROLYTE)
    reference = model.reference(
        T=298.15 * u.kelvin,
        P=1.0 * u.bar,
        solvent_x=(1.0, 0.0, 0.0),
        phase="liquid",
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
    log_activity = tuple(
        finite.fugacity.ln_coefficient[index] - infinite
        for index, infinite in zip(
            (1, 2),
            reference.infinite_dilution_log_fugacity_coefficient,
            strict=True,
        )
    )
    assert state.activity.single_ion == pytest.approx(
        tuple(math.exp(value) for value in log_activity)
    )
    assert state.activity.mean_mole_fraction == pytest.approx(
        math.exp(math.fsum(log_activity) / 2.0)
    )
    assert state.activity.osmotic_coefficient is not None
    assert math.isfinite(state.activity.osmotic_coefficient)

    scale = GAS_CONSTANT * 298.15
    assert tuple(
        value.to("joule / mole").magnitude for value in reference.solvation_free_energy
    ) == pytest.approx(
        tuple(scale * value for value in reference.infinite_dilution_log_fugacity_coefficient)
    )
    identity = reference.transfer_to(reference)
    assert isinstance(identity, TransferFreeEnergy)
    assert tuple(value.to("joule / mole").magnitude for value in identity.values) == (0.0, 0.0)
