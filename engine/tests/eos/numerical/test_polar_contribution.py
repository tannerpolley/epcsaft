from __future__ import annotations

import math
from copy import deepcopy

import pytest

from epcsaft import (
    Mixture,
    Parameters,
    Pure,
)
from epcsaft import (
    unit_registry as u,
)
from tests.eos._raw_cases import (
    DIPOLAR_PERMITTIVITY_ELECTROLYTE,
    POLAR_BINARY,
    POLAR_PURE,
    component_values,
    pair_matrix,
    set_binary_pair_value,
    set_component_values,
)

COEFFICIENTS = {
    "dd_a0": (0.3043504, -0.1358588, 1.4493329, 0.3556977, -2.0653308),
    "dd_a1": (0.9534641, -1.8396383, 2.0131180, -7.3724958, 8.2374135),
    "dd_a2": (-1.1610080, 4.5258607, 0.9751222, -12.281038, 5.9397575),
    "dd_b0": (0.2187939, -1.1896431, 1.1626889, 0.0, 0.0),
    "dd_b1": (-0.5873164, 1.2489132, -0.5085280, 0.0, 0.0),
    "dd_b2": (3.4869576, -14.915974, 15.372022, 0.0, 0.0),
    "dd_c0": (-0.0646774, 0.1975882, -0.8087562, 0.6902849, 0.0),
    "dd_c1": (-0.9520876, 2.9924258, -2.3802636, -0.2701261, 0.0),
    "dd_c2": (-0.6260979, 1.2924686, 1.6542783, -3.4396744, 0.0),
    "qq_a0": (1.2378308, 2.4355031, 1.6330905, -1.6118152, 6.9771185),
    "qq_a1": (1.2854109, -11.465615, 22.086893, 7.4691383, -17.197772),
    "qq_a2": (1.7942954, 0.7695103, 7.2647923, 94.486699, -77.148458),
    "qq_b0": (0.4542718, -4.5016264, 3.5858868, 0.0, 0.0),
    "qq_b1": (-0.8137340, 10.064030, -10.876631, 0.0, 0.0),
    "qq_b2": (6.8682675, -5.1732238, -17.240207, 0.0, 0.0),
    "qq_c0": (-0.5000437, 6.5318692, -16.014780, 14.425970, 0.0),
    "qq_c1": (2.0002094, -6.7838658, 20.383246, -10.895984, 0.0),
    "qq_c2": (3.1358271, 7.2475888, 3.0759478, 0.0, 0.0),
    "dq_a0": (0.6970950, -0.6335541, 2.9455090, -1.4670273, 0.0),
    "dq_a1": (-0.6734593, -1.4258991, 4.1944139, 1.0266216, 0.0),
    "dq_a2": (0.6703408, -4.3384718, 7.2341684, 0.0, 0.0),
    "dq_b0": (-0.4840383, 1.9704055, -2.1185727, 0.0, 0.0),
    "dq_b1": (0.6765101, -3.0138675, 0.4674266, 0.0, 0.0),
    "dq_b2": (-1.1675601, 2.1348843, 0.0, 0.0, 0.0),
    "dq_c0": (7.846431, 33.42700, 4.689111, 0.0, 0.0),
    "dq_c1": (-20.72202, -58.63904, -1.764887, 0.0, 0.0),
    "dq_alpha": 1.19374,
}

AVOGADRO = 6.02214076e23
BOLTZMANN = 1.380649e-23
VACUUM_PERMITTIVITY = 8.8541878128e-12
DEBYE_C_M = 1.0e-21 / 299792458.0
ANGSTROM_M = 1.0e-10


def _state(mapping=POLAR_BINARY):
    model = Mixture(Parameters.from_mapping(mapping))
    return model, model.state(
        T=315.0 * u.kelvin,
        rho=4200.0 * u.mole / u.meter**3,
        x=(0.35, 0.65),
    )


def _chain(values, chain):
    zero, one, two = values
    first = (chain - 1.0) / chain
    return zero + first * one + first * (chain - 2.0) / chain * two


def _poly(values, eta):
    return sum(value * eta**order for order, value in enumerate(values))


def _integral(family, kind, chain, eta, epsilon_over_temperature=0.0):
    if kind == "j2":
        values = tuple(
            _chain(
                tuple(COEFFICIENTS[f"{family}_a{column}"][order] for column in range(3)),
                chain,
            )
            + epsilon_over_temperature
            * _chain(
                tuple(COEFFICIENTS[f"{family}_b{column}"][order] for column in range(3)),
                chain,
            )
            for order in range(5)
        )
    elif family == "dq":
        values = tuple(
            COEFFICIENTS["dq_c0"][order] + (chain - 1.0) / chain * COEFFICIENTS["dq_c1"][order]
            for order in range(5)
        )
    else:
        values = tuple(
            _chain(
                tuple(COEFFICIENTS[f"{family}_c{column}"][order] for column in range(3)),
                chain,
            )
            for order in range(5)
        )
    return _poly(values, eta)


def _reduced_moments(mapping):
    segments = component_values(mapping, "segment_count")
    sigma = component_values(mapping, "segment_diameter")
    epsilon = component_values(mapping, "dispersion_energy_over_k")
    dipoles = component_values(mapping, "dipole_moment")
    quadrupoles = component_values(mapping, "quadrupole_moment")
    mu2 = []
    q2 = []
    for m, sig, eps, dipole, quadrupole in zip(
        segments, sigma, epsilon, dipoles, quadrupoles, strict=True
    ):
        denominator = 4.0 * math.pi * VACUUM_PERMITTIVITY * m * BOLTZMANN * eps
        mu = 0.0 if dipole is None else dipole * DEBYE_C_M
        q = 0.0 if quadrupole is None else quadrupole * DEBYE_C_M * ANGSTROM_M
        mu2.append(mu * mu / (denominator * (sig * ANGSTROM_M) ** 3))
        q2.append(q * q / (denominator * (sig * ANGSTROM_M) ** 5))
    return tuple(mu2), tuple(q2)


def _paper_terms(
    mapping,
    temperature,
    density,
    fractions,
    *,
    dd_cap=True,
    dq_cap=False,
    dq_diameter_powers=(3, 5),
):
    segments = component_values(mapping, "segment_count")
    sigma = component_values(mapping, "segment_diameter")
    epsilon = component_values(mapping, "dispersion_energy_over_k")
    kij = pair_matrix(mapping, "k_ij")
    mu2, q2 = _reduced_moments(mapping)
    number_density = density * AVOGADRO / 1.0e30
    diameter = tuple(
        sig * (1.0 - 0.12 * math.exp(-3.0 * eps / temperature))
        for sig, eps in zip(sigma, epsilon, strict=True)
    )
    eta = (
        math.pi
        / 6.0
        * number_density
        * math.fsum(x * m * d**3 for x, m, d in zip(fractions, segments, diameter, strict=True))
    )

    dd_pair = dd_triplet = qq_pair = qq_triplet = 0.0
    dq_pair = dq_triplet = 0.0
    for i in range(len(segments)):
        for j in range(len(segments)):
            sigma_ij = 0.5 * (sigma[i] + sigma[j])
            epsilon_ij = (1.0 - kij[i][j]) * math.sqrt(epsilon[i] * epsilon[j])
            if mu2[i] and mu2[j]:
                chain = math.sqrt(segments[i] * segments[j])
                if dd_cap:
                    chain = min(2.0, chain)
                dd_pair += (
                    fractions[i]
                    * fractions[j]
                    * epsilon[i]
                    * epsilon[j]
                    / temperature**2
                    * sigma[i] ** 3
                    * sigma[j] ** 3
                    / sigma_ij**3
                    * mu2[i]
                    * mu2[j]
                    * _integral("dd", "j2", chain, eta, epsilon_ij / temperature)
                )
            if q2[i] and q2[j]:
                chain = math.sqrt(segments[i] * segments[j])
                qq_pair += (
                    fractions[i]
                    * fractions[j]
                    * epsilon[i]
                    * epsilon[j]
                    / temperature**2
                    * sigma[i] ** 5
                    * sigma[j] ** 5
                    / sigma_ij**7
                    * q2[i]
                    * q2[j]
                    * _integral("qq", "j2", chain, eta, epsilon_ij / temperature)
                )
            if mu2[i] and q2[j]:
                chain = math.sqrt(segments[i] * segments[j])
                dq_pair += (
                    fractions[i]
                    * fractions[j]
                    * epsilon[i]
                    * epsilon[j]
                    / temperature**2
                    * sigma[i] ** dq_diameter_powers[0]
                    * sigma[j] ** dq_diameter_powers[1]
                    / sigma_ij**5
                    * mu2[i]
                    * q2[j]
                    * _integral("dq", "j2", chain, eta, epsilon_ij / temperature)
                )
            for k in range(len(segments)):
                sigma_ik = 0.5 * (sigma[i] + sigma[k])
                sigma_jk = 0.5 * (sigma[j] + sigma[k])
                common = (
                    fractions[i]
                    * fractions[j]
                    * fractions[k]
                    * epsilon[i]
                    * epsilon[j]
                    * epsilon[k]
                    / temperature**3
                )
                if mu2[i] and mu2[j] and mu2[k]:
                    chain = (segments[i] * segments[j] * segments[k]) ** (1.0 / 3.0)
                    if dd_cap:
                        chain = min(2.0, chain)
                    dd_triplet += (
                        common
                        * sigma[i] ** 3
                        * sigma[j] ** 3
                        * sigma[k] ** 3
                        / (sigma_ij * sigma_ik * sigma_jk)
                        * mu2[i]
                        * mu2[j]
                        * mu2[k]
                        * _integral("dd", "j3", chain, eta)
                    )
                if q2[i] and q2[j] and q2[k]:
                    chain = (segments[i] * segments[j] * segments[k]) ** (1.0 / 3.0)
                    qq_triplet += (
                        common
                        * sigma[i] ** 5
                        * sigma[j] ** 5
                        * sigma[k] ** 5
                        / (sigma_ij**3 * sigma_ik**3 * sigma_jk**3)
                        * q2[i]
                        * q2[j]
                        * q2[k]
                        * _integral("qq", "j3", chain, eta)
                    )
                moment_product = (
                    mu2[i] * mu2[j] * q2[k] + COEFFICIENTS["dq_alpha"] * mu2[i] * q2[j] * q2[k]
                )
                if moment_product:
                    chain = (segments[i] * segments[j] * segments[k]) ** (1.0 / 3.0)
                    if dq_cap:
                        chain = min(2.0, chain)
                    dq_triplet += (
                        common
                        * sigma[i] ** 4
                        * sigma[j] ** 4
                        * sigma[k] ** 4
                        / (sigma_ij**2 * sigma_ik**2 * sigma_jk**2)
                        * moment_product
                        * _integral("dq", "j3", chain, eta)
                    )
    return (
        -math.pi * number_density * dd_pair,
        -4.0 * math.pi**2 / 3.0 * number_density**2 * dd_triplet,
        -9.0 * math.pi / 16.0 * number_density * qq_pair,
        9.0 * math.pi**2 / 16.0 * number_density**2 * qq_triplet,
        -9.0 * math.pi / 4.0 * number_density * dq_pair,
        -(number_density**2) * dq_triplet,
    )


def _pade_sum(terms):
    return math.fsum(
        0.0 if second == 0.0 and third == 0.0 else second * second / (second - third)
        for second, third in zip(terms[::2], terms[1::2], strict=True)
    )


def test_polar_terms__match_primary_paper_oracle_and_tables() -> None:
    _, state = _state()
    expected = _paper_terms(POLAR_BINARY, 315.0, 4200.0, (0.35, 0.65))
    assert state.ares_polar == pytest.approx(_pade_sum(expected), rel=2.0e-13)

    pure_model = Pure(Parameters.from_mapping(POLAR_PURE))
    pure_expected = _paper_terms(POLAR_PURE, 315.0, 4200.0, (1.0,))
    pure_state = pure_model.state(
        T=315.0 * u.kelvin,
        rho=4200.0 * u.mole / u.meter**3,
    )
    assert pure_state.ares_polar == pytest.approx(_pade_sum(pure_expected), rel=2.0e-13)


@pytest.mark.parametrize(
    "moment_family",
    ("dipole_moment", "quadrupole_moment"),
)
def test_pure_family_terms__match_primary_paper_oracles(moment_family) -> None:
    mapping = deepcopy(POLAR_PURE)
    other = "quadrupole_moment" if moment_family == "dipole_moment" else "dipole_moment"
    set_component_values(mapping, other, (None,))
    expected = _paper_terms(mapping, 315.0, 4200.0, (1.0,))
    state = Pure(Parameters.from_mapping(mapping)).state(
        T=315.0 * u.kelvin,
        rho=4200.0 * u.mole / u.meter**3,
    )
    assert state.ares_polar == pytest.approx(_pade_sum(expected), rel=2.0e-13)


def test_polar_chain_rules__cap_dd_and_leave_dq_uncapped_above_two() -> None:
    mapping = deepcopy(POLAR_PURE)
    set_component_values(mapping, "segment_count", (3.25,))
    observed = Pure(Parameters.from_mapping(mapping)).state(
        T=315.0 * u.kelvin,
        rho=4200.0 * u.mole / u.meter**3,
    ).ares_polar
    expected = _paper_terms(mapping, 315.0, 4200.0, (1.0,))
    uncapped_dd = _paper_terms(mapping, 315.0, 4200.0, (1.0,), dd_cap=False)
    capped_dq = _paper_terms(mapping, 315.0, 4200.0, (1.0,), dq_cap=True)

    assert observed == pytest.approx(_pade_sum(expected), rel=2.0e-13)
    assert observed != pytest.approx(_pade_sum(uncapped_dd), rel=2.0e-13)
    assert observed != pytest.approx(_pade_sum(capped_dq), rel=2.0e-13)


def test_polar_formulation__ignores_diameter_interaction() -> None:
    without_lij = deepcopy(POLAR_BINARY)
    set_binary_pair_value(without_lij, "l_ij", 0.0)

    _, reference = _state()
    _, changed = _state(without_lij)

    assert changed.ares_polar == pytest.approx(reference.ares_polar, rel=0.0, abs=0.0)


@pytest.mark.parametrize(
    ("moment_family", "moments"),
    (
        ("dipole_moment", (1.7, 2.1)),
        ("quadrupole_moment", (3.6, 4.1)),
    ),
)
def test_dd_and_qq_unlike_mixing__matches_public_oracle(moment_family, moments) -> None:
    mapping = deepcopy(POLAR_BINARY)
    set_component_values(mapping, "dipole_moment", (None, None))
    set_component_values(mapping, "quadrupole_moment", (None, None))
    set_component_values(mapping, moment_family, moments)
    no_kij = deepcopy(mapping)
    set_binary_pair_value(no_kij, "k_ij", 0.0)
    no_lij = deepcopy(mapping)
    set_binary_pair_value(no_lij, "l_ij", 0.0)

    def public_ares(values):
        return Mixture(Parameters.from_mapping(values)).state(
            T=315.0 * u.kelvin,
            rho=4200.0 * u.mole / u.meter**3,
            x=(0.35, 0.65),
        ).ares_polar

    reference = public_ares(mapping)
    changed_kij = public_ares(no_kij)
    changed_lij = public_ares(no_lij)
    assert reference == pytest.approx(
        _pade_sum(_paper_terms(mapping, 315.0, 4200.0, (0.35, 0.65))), rel=2.0e-13
    )
    assert changed_kij == pytest.approx(
        _pade_sum(_paper_terms(no_kij, 315.0, 4200.0, (0.35, 0.65))), rel=2.0e-13
    )
    assert changed_kij != pytest.approx(reference)
    assert changed_lij == pytest.approx(reference, rel=0.0, abs=0.0)


def test_polar_activation_and_dilute_limit_are_publicly_observable() -> None:
    inactive = deepcopy(POLAR_BINARY)
    set_component_values(inactive, "dipole_moment", (None, None))
    set_component_values(inactive, "quadrupole_moment", (None, None))
    _, inactive_state = _state(inactive)
    assert inactive_state.ares_polar == 0.0

    dd_only = deepcopy(POLAR_BINARY)
    set_component_values(dd_only, "quadrupole_moment", (None, None))
    _, dd_state = _state(dd_only)
    assert dd_state.ares_polar == pytest.approx(
        _pade_sum(_paper_terms(dd_only, 315.0, 4200.0, (0.35, 0.65))), rel=2.0e-13
    )

    qq_only = deepcopy(POLAR_BINARY)
    set_component_values(qq_only, "dipole_moment", (None, None))
    _, qq_state = _state(qq_only)
    assert qq_state.ares_polar == pytest.approx(
        _pade_sum(_paper_terms(qq_only, 315.0, 4200.0, (0.35, 0.65))), rel=2.0e-13
    )

    dilute = Mixture(Parameters.from_mapping(POLAR_BINARY)).state(
        T=315.0 * u.kelvin,
        rho=1.0 * u.mole / u.meter**3,
        x=(0.35, 0.65),
    )
    assert math.isfinite(dilute.ares_polar)
    assert 0.0 < abs(dilute.ares_polar) < abs(_state()[1].ares_polar)

    dielectric_dipole = Mixture(
        Parameters.from_mapping(DIPOLAR_PERMITTIVITY_ELECTROLYTE)
    ).state(
        T=298.15 * u.kelvin,
        rho=1000.0 * u.mole / u.meter**3,
        x=(0.98, 0.01, 0.01),
    )
    assert dielectric_dipole.ares_polar == 0.0


def test_polar_identity_permutation_and_cppad_properties_are_consistent() -> None:
    changed_moment = deepcopy(POLAR_BINARY)
    set_component_values(changed_moment, "dipole_moment", (1.8, None))
    assert (
        Parameters.from_mapping(changed_moment).fingerprint
        != Parameters.from_mapping(POLAR_BINARY).fingerprint
    )

    reverse = deepcopy(POLAR_BINARY)
    reverse["components"] = list(reversed(reverse["components"]))
    forward_model, forward = _state()
    reverse_model = Mixture(Parameters.from_mapping(reverse))
    reversed_state = reverse_model.state(
        T=315.0 * u.kelvin,
        rho=4200.0 * u.mole / u.meter**3,
        x=(0.65, 0.35),
    )
    assert reversed_state.ares_polar == pytest.approx(forward.ares_polar, rel=2.0e-14)
    assert reverse_model.parameter_fingerprint != forward_model.parameter_fingerprint

    closure_reference = forward_model.state(
        T=315.0 * u.kelvin,
        rho=100.0 * u.mole / u.meter**3,
        x=(0.35, 0.65),
    )
    closure = forward_model.state(
        T=315.0 * u.kelvin,
        P=closure_reference.pressure,
        x=(0.35, 0.65),
        phase="vapor",
    )
    assert closure.molar_density.to("mole / meter ** 3").magnitude == pytest.approx(
        100.0, rel=2.0e-10
    )
    assert closure.ares_polar == pytest.approx(closure_reference.ares_polar, rel=2.0e-10)

    temperature_step = 315.0e-4
    density_step = 4200.0e-4

    def ares(temperature, density, fractions=(0.35, 0.65)):
        return forward_model.state(
            T=temperature * u.kelvin,
            rho=density * u.mole / u.meter**3,
            x=fractions,
        ).ares(basis="reduced")

    a_temperature = (
        ares(315.0 + temperature_step, 4200.0) - ares(315.0 - temperature_step, 4200.0)
    ) / (2.0 * temperature_step)
    a_density = (ares(315.0, 4200.0 + density_step) - ares(315.0, 4200.0 - density_step)) / (
        2.0 * density_step
    )
    assert forward.sres(basis="reduced") == pytest.approx(
        -forward.ares(basis="reduced") - 315.0 * a_temperature,
        rel=2.0e-7,
        abs=2.0e-8,
    )
    assert forward.compressibility_factor == pytest.approx(
        1.0 + 4200.0 * a_density, rel=2.0e-7, abs=2.0e-8
    )

    volume = 1.0 / 4200.0
    expected_chemical_potential = []
    for component in range(2):
        amount_step = 2.0e-5

        def total_ares(delta, component=component):
            amounts = [0.35, 0.65]
            amounts[component] += delta
            total = math.fsum(amounts)
            return total * ares(
                315.0,
                total / volume,
                tuple(amount / total for amount in amounts),
            )

        expected_chemical_potential.append(
            (total_ares(amount_step) - total_ares(-amount_step)) / (2.0 * amount_step)
        )
    assert forward.residual_chemical_potential_over_rt == pytest.approx(
        expected_chemical_potential, rel=3.0e-7, abs=3.0e-8
    )
