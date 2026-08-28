from __future__ import annotations

import math
from copy import deepcopy

import pytest

from epcsaft import Mixture, Parameters, Pure
from epcsaft import unit_registry as u
from tests.eos._cases import (
    AQUEOUS_ELECTROLYTE,
    ASSOCIATING_PURE,
    DIPOLAR_PERMITTIVITY_ELECTROLYTE,
    NEUTRAL_BINARY,
)

AVOGADRO = 6.02214076e23
BOLTZMANN = 1.380649e-23
ELEMENTARY_CHARGE = 1.602176634e-19
VACUUM_PERMITTIVITY = 8.8541878128e-12

DISPERSION_A_COEFFICIENTS = (
    (0.9105631445, -0.3084016918, -0.0906148351),
    (0.6361281449, 0.1860531159, 0.4527842806),
    (2.6861347891, -2.5030047259, 0.5962700728),
    (-26.547362491, 21.419793629, -1.7241829131),
    (97.759208784, -65.255885330, -4.1302112531),
    (-159.59154087, 83.318680481, 13.776631870),
    (91.297774084, -33.746922930, -8.6728470368),
)
DISPERSION_B_COEFFICIENTS = (
    (0.7240946941, -0.5755498075, 0.0976883116),
    (2.2382791861, 0.6995095521, -0.2557574982),
    (-4.0025849485, 3.8925673390, -9.1558561530),
    (-21.003576815, -17.215471648, 20.642075974),
    (26.855641363, 192.67226447, -38.804430052),
    (206.55133841, -161.82646165, 93.626774077),
    (-355.60235612, -165.20769346, -29.666905585),
)


def _state(mapping, *, composition, density=1000.0, temperature=298.15):
    parameters = Parameters.from_mapping(mapping)
    if len(parameters.component_ids) == 1:
        return Pure(parameters).state(
            T=temperature * u.kelvin,
            rho=density * u.mole / u.meter**3,
        )
    return Mixture(parameters).state(
        T=temperature * u.kelvin,
        rho=density * u.mole / u.meter**3,
        x=composition,
    )


def _hard_chain_and_dispersion(mapping, temperature, density, composition):
    components = mapping["components"]
    component_ids = [component["component_id"] for component in components]

    def coefficient(component, family):
        return next(
            item["value"]["magnitude"]
            for item in component["coefficients"]
            if item["family"] == family
        )

    def pair_coefficient(left, right, family):
        if left == right:
            return 0.0
        left, right = sorted((left, right))
        pair = next(
            item
            for item in mapping["pairs"]
            if item["component_id_a"] == left and item["component_id_b"] == right
        )
        return next(
            (
                item["value"]["magnitude"]
                for item in pair["coefficients"]
                if item["family"] == family
            ),
            0.0,
        )

    segment_count = [coefficient(component, "segment_count") for component in components]
    sigma = [coefficient(component, "segment_diameter") for component in components]
    epsilon = [
        coefficient(component, "dispersion_energy_over_k") for component in components
    ]
    diameter = [
        size * (1.0 - 0.12 * math.exp(-3.0 * energy / temperature))
        for size, energy in zip(sigma, epsilon, strict=True)
    ]
    number_density = density * AVOGADRO * 1.0e-30
    zeta = tuple(
        math.pi
        / 6.0
        * number_density
        * math.fsum(
            fraction * segments * size**order
            for fraction, segments, size in zip(
                composition, segment_count, diameter, strict=True
            )
        )
        for order in range(4)
    )
    eta = zeta[3]
    one_minus_eta = 1.0 - eta
    hard_sphere = (
        3.0 * zeta[1] * zeta[2] / one_minus_eta
        + zeta[2] ** 3 / (eta * one_minus_eta**2)
        + (zeta[2] ** 3 / eta**2 - zeta[0]) * math.log(one_minus_eta)
    ) / zeta[0]
    mean_segments = math.fsum(
        fraction * segments
        for fraction, segments in zip(composition, segment_count, strict=True)
    )
    hard_chain = mean_segments * hard_sphere - math.fsum(
        fraction
        * (segments - 1.0)
        * math.log(
            1.0 / one_minus_eta
            + 1.5 * size * zeta[2] / one_minus_eta**2
            + 0.5 * size**2 * zeta[2] ** 2 / one_minus_eta**3
        )
        for fraction, segments, size in zip(
            composition, segment_count, diameter, strict=True
        )
    )

    chain_one = (mean_segments - 1.0) / mean_segments
    chain_two = (mean_segments - 1.0) * (mean_segments - 2.0) / mean_segments**2
    integral_one = math.fsum(
        (base + chain_one * first + chain_two * second) * eta**order
        for order, (base, first, second) in enumerate(DISPERSION_A_COEFFICIENTS)
    )
    integral_two = math.fsum(
        (base + chain_one * first + chain_two * second) * eta**order
        for order, (base, first, second) in enumerate(DISPERSION_B_COEFFICIENTS)
    )
    correction = 1.0 / (
        1.0
        + mean_segments * (8.0 * eta - 2.0 * eta**2) / one_minus_eta**4
        + (1.0 - mean_segments)
        * (20.0 * eta - 27.0 * eta**2 + 12.0 * eta**3 - 2.0 * eta**4)
        / ((one_minus_eta * (2.0 - eta)) ** 2)
    )
    first_moment = 0.0
    second_moment = 0.0
    for i, left in enumerate(component_ids):
        for j, right in enumerate(component_ids):
            cross_sigma = 0.5 * (sigma[i] + sigma[j]) * (
                1.0 - pair_coefficient(left, right, "l_ij")
            )
            cross_epsilon = math.sqrt(epsilon[i] * epsilon[j]) * (
                1.0 - pair_coefficient(left, right, "k_ij")
            )
            weight = (
                composition[i]
                * composition[j]
                * segment_count[i]
                * segment_count[j]
                * cross_sigma**3
            )
            first_moment += weight * cross_epsilon / temperature
            second_moment += weight * (cross_epsilon / temperature) ** 2
    dispersion = (
        -2.0 * math.pi * number_density * integral_one * first_moment
        - math.pi
        * number_density
        * mean_segments
        * correction
        * integral_two
        * second_moment
    )
    return hard_chain, dispersion


@pytest.mark.parametrize(
    ("mapping", "composition"),
    (
        (NEUTRAL_BINARY, (0.4, 0.6)),
        (AQUEOUS_ELECTROLYTE, (0.98, 0.01, 0.01)),
    ),
    ids=("neutral", "electrolyte"),
)
def test_residual_helmholtz_closes_over_retained_contributions(mapping, composition) -> None:
    state = _state(mapping, composition=composition)

    assert state.ares(basis="reduced") == pytest.approx(
        state.hard_chain
        + state.dispersion
        + state.association
        + state.debye_huckel
        + state.born,
        rel=2.0e-14,
        abs=2.0e-14,
    )


@pytest.mark.parametrize(
    ("temperature", "density", "composition"),
    (
        (200.0, 1000.0, (1.0, 0.0)),
        (250.0, 8000.0, (0.4, 0.6)),
    ),
    ids=("pure-limit", "binary"),
)
def test_hard_chain_and_dispersion_match_source_equations(
    temperature, density, composition
) -> None:
    # Gross and Sadowski, Ind. Eng. Chem. Res. 40 (2001) 1244-1260,
    # Appendix A, Eqs. (A.6), (A.8)-(A.17), with universal coefficients from
    # Table 1. Basis: a = Ares/(N k_B T) = Ares,m/(R T); rho is converted from
    # mol/m3 to particles/angstrom3. The calculation reads only the public case
    # mapping and does not call an Engine contribution routine.
    hard_chain, dispersion = _hard_chain_and_dispersion(
        NEUTRAL_BINARY, temperature, density, composition
    )
    state = _state(
        NEUTRAL_BINARY,
        composition=composition,
        density=density,
        temperature=temperature,
    )

    assert state.hard_chain == pytest.approx(hard_chain, rel=5.0e-12, abs=5.0e-13)
    assert state.dispersion == pytest.approx(dispersion, rel=5.0e-12, abs=5.0e-13)


def test_association_matches_two_site_closed_form() -> None:
    temperature = 300.0
    density = 10000.0
    sigma = 3.0
    diameter = sigma * (1.0 - 0.12 * math.exp(-3.0 * 200.0 / temperature))
    number_density = density * AVOGADRO * 1.0e-30
    eta = math.pi * number_density * diameter**3 / 6.0
    contact = (
        1.0 / (1.0 - eta)
        + 1.5
        * diameter
        * number_density
        * math.pi
        * diameter**2
        / 6.0
        / (1.0 - eta) ** 2
        + 0.5
        * diameter**2
        * (number_density * math.pi * diameter**2 / 6.0) ** 2
        / (1.0 - eta) ** 3
    )
    delta = contact * math.expm1(2400.0 / temperature) * sigma**3 * 0.03
    site_fraction = 2.0 / (1.0 + math.sqrt(1.0 + 4.0 * number_density * delta))
    expected = 2.0 * (math.log(site_fraction) - 0.5 * site_fraction + 0.5)

    state = _state(
        ASSOCIATING_PURE,
        composition=(1.0,),
        density=density,
        temperature=temperature,
    )

    assert state.association == pytest.approx(expected, rel=2.0e-11, abs=2.0e-12)


def test_debye_huckel_matches_independent_si_equation() -> None:
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
        fraction * charge * charge
        for fraction, charge in zip(fractions, charges, strict=True)
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
        chi = 3.0 / reduced**3 * (
            1.5
            + math.log1p(reduced)
            - 2.0 * (1.0 + reduced)
            + 0.5 * (1.0 + reduced) ** 2
        )
        weighted_finite_size += fraction * charge * charge * chi
    expected = (
        -kappa
        * ELEMENTARY_CHARGE**2
        * weighted_finite_size
        / (12.0 * math.pi * VACUUM_PERMITTIVITY * permittivity * BOLTZMANN * temperature)
    )

    state = _state(AQUEOUS_ELECTROLYTE, composition=fractions, density=density)

    assert state.debye_huckel == pytest.approx(expected, rel=2.0e-11, abs=2.0e-13)


@pytest.mark.parametrize(
    ("c_shell", "c_dielectric"),
    ((0.0, 0.0), (1.0, 0.0), (1.0, 1.0)),
    ids=("base", "shell", "shell-dielectric"),
)
def test_born_matches_independent_nested_equations(c_shell, c_dielectric) -> None:
    temperature = 298.15
    fractions = (0.90, 0.05, 0.05)
    charges = (0, 1, -1)
    diameters = (0.0, 3.445, 4.100)
    solvation = (1.5, 1.0, 1.0)
    mapping = deepcopy(AQUEOUS_ELECTROLYTE)
    electrolyte = next(
        family for family in mapping["model_families"] if family["kind"] == "electrolyte"
    )
    electrolyte["c_shell"] = c_shell
    electrolyte["c_dielectric"] = c_dielectric
    state = _state(mapping, composition=fractions, density=40000.0)

    mixed_solvation = math.fsum(
        fraction * factor for fraction, factor in zip(fractions, solvation, strict=True)
    )
    bulk_permittivity = 78.09 / (1.0 + 7.01 * math.fsum(fractions[1:]))
    radial_sum = 0.0
    for fraction, charge, diameter in zip(fractions, charges, diameters, strict=True):
        if charge == 0:
            continue
        born_diameter = diameter * 1.0e-10
        shell_diameter = born_diameter * (
            1.0 + c_shell * (mixed_solvation - 1.0) / abs(charge)
        )
        radial_sum += fraction * charge * charge * (
            (1.0 - 1.0 / bulk_permittivity) / shell_diameter
            + c_dielectric
            * (1.0 - 1.0 / 8.0)
            * (1.0 / born_diameter - 1.0 / shell_diameter)
        )
    expected = (
        -(ELEMENTARY_CHARGE**2)
        * radial_sum
        / (4.0 * math.pi * VACUUM_PERMITTIVITY * BOLTZMANN * temperature)
    )

    assert state.born == pytest.approx(expected, rel=3.0e-12, abs=2.0e-12)


def test_relative_permittivity_matches_both_retained_equations() -> None:
    fractions = (0.9, 0.05, 0.05)
    density = 55000.0
    ion_fraction = _state(AQUEOUS_ELECTROLYTE, composition=fractions, density=density)
    assert ion_fraction.bulk_relative_permittivity == pytest.approx(
        78.09 / (1.0 + 7.01 * 0.1), rel=2.0e-12
    )

    number_density = density * AVOGADRO * 1.0e-30
    beta = 1.0 / (BOLTZMANN * 298.15)
    effective = 4.0 * math.pi / 9.0 * number_density * fractions[0] * (
        beta * 1.0e-19 + 3.0
    )
    integral = 1.0 + math.expm1(-effective)
    expected = 1.0 + 3.0 * effective * (
        1.0 + effective + ((17.0 / 16.0) * integral - 1.0) * effective**2
    )

    dipolar_density = _state(
        DIPOLAR_PERMITTIVITY_ELECTROLYTE,
        composition=fractions,
        density=density,
    )

    assert dipolar_density.bulk_relative_permittivity == pytest.approx(
        expected, rel=2.0e-12
    )


def test_dipolar_density_permittivity_preserves_density_and_temperature_derivatives() -> None:
    fractions = (0.9, 0.05, 0.05)
    density = 55000.0
    temperature = 298.15
    mapping = deepcopy(DIPOLAR_PERMITTIVITY_ELECTROLYTE)
    mapping["domains"][0]["temperature_min"]["magnitude"] = 297.0
    mapping["domains"][0]["temperature_max"]["magnitude"] = 299.0

    def helmholtz(*, rho=density, T=temperature) -> float:
        return _state(
            mapping,
            composition=fractions,
            density=rho,
            temperature=T,
        ).ares(basis="reduced")

    density_step = density * 1.0e-5
    density_derivative = (
        helmholtz(rho=density + density_step)
        - helmholtz(rho=density - density_step)
    ) / (2.0 * density_step)
    temperature_step = temperature * 1.0e-4
    temperature_derivative = (
        helmholtz(T=temperature + temperature_step)
        - helmholtz(T=temperature - temperature_step)
    ) / (2.0 * temperature_step)
    state = _state(
        mapping,
        composition=fractions,
        density=density,
        temperature=temperature,
    )

    assert state.compressibility_factor == pytest.approx(
        1.0 + density * density_derivative, rel=2.0e-6, abs=2.0e-8
    )
    assert state.ures(basis="reduced") == pytest.approx(
        -temperature * temperature_derivative, rel=2.0e-6, abs=2.0e-8
    )
