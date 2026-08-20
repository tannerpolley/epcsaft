from __future__ import annotations

import math
from copy import deepcopy

import pytest

from epcsaft import (
    ActiveParameterSet,
    Activity,
    EosError,
    Mixture,
    ParameterError,
    Parameters,
    Pure,
    TransferFreeEnergy,
)
from epcsaft import unit_registry as u
from tests.eos._raw_cases import (
    AQUEOUS_ELECTROLYTE,
    ASSOCIATING_PURE,
    DIPOLAR_PERMITTIVITY_ELECTROLYTE,
    HELD_STRATEGY_2_TERNARY,
    MULTINEUTRAL_REACTIVE,
    MULTISOLVENT_MULTIION,
    NEUTRAL_BINARY,
    _component,
    _document,
    _model_coefficient,
    _model_families,
    _pair,
    neutral_binary,
    neutral_pure,
)

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
    parameters = Parameters.from_mapping(mapping)
    if len(parameters.component_ids) == 1:
        return Pure(parameters).state(
            T=temperature * u.kelvin,
            rho=density * u.mole / u.meter**3,
        )
    return Mixture(parameters).state(
        T=temperature * u.kelvin,
        rho=density * u.mole / u.meter**3,
        x=x,
    )


def _revised_debye_huckel(mapping):
    candidate = deepcopy(mapping)
    electrolyte = next(
        family for family in candidate["model_families"] if family["kind"] == "electrolyte"
    )
    electrolyte["choice"] = "revised-debye-huckel"
    electrolyte.pop("c_shell")
    electrolyte.pop("c_dielectric")
    permittivity = next(
        family for family in candidate["model_families"] if family["kind"] == "permittivity"
    )
    permittivity["choice"] = "solvent-only"
    candidate["model_coefficients"] = []
    for component in candidate["components"]:
        charged = component["fixed"]["charge_number"]["value"]["magnitude"] != 0
        component["coefficients"] = [
            coefficient
            for coefficient in component["coefficients"]
            if coefficient["family"] not in {"born_diameter", "solvation_factor"}
            and not (charged and coefficient["family"] == "relative_permittivity")
        ]
    return candidate


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


def test_fully_dissociated_debye_huckel__retains_declared_opposite_charge_dispersion() -> None:
    state = Mixture(Parameters.from_mapping(HELD_STRATEGY_2_TERNARY)).state(
        T=298.15 * u.kelvin,
        P=1.0 * u.bar,
        x=(0.96, 0.02, 0.02),
        phase="liquid",
    )
    suppressed_mapping = deepcopy(HELD_STRATEGY_2_TERNARY)
    pair = next(
        record
        for record in suppressed_mapping["pairs"]
        if {record["component_id_a"], record["component_id_b"]}
        == {"lithium-cation", "chloride-anion"}
    )
    pair["coefficients"][0]["value"]["magnitude"] = 1.0
    suppressed = Mixture(Parameters.from_mapping(suppressed_mapping)).state(
        T=298.15 * u.kelvin,
        P=1.0 * u.bar,
        x=(0.96, 0.02, 0.02),
        phase="liquid",
    )

    assert state.dispersion == pytest.approx(-8.936297809884787, rel=2.0e-13)
    assert suppressed.dispersion == pytest.approx(-8.93489498474385, rel=2.0e-13)
    assert state.fugacity is not None
    assert suppressed.fugacity is not None
    assert state.fugacity.ln_coefficient[1:] == pytest.approx(
        (-3.0464389064816944, 2.9426852991208703), rel=2.0e-12
    )
    assert state.fugacity.ln_coefficient[1:] != pytest.approx(
        suppressed.fugacity.ln_coefficient[1:]
    )


def test_component_permutation__preserves_scalar_thermodynamics() -> None:
    forward = _state(NEUTRAL_BINARY, x=(0.4, 0.6))
    reverse = _state(neutral_binary(("ethane", "methane")), x=(0.6, 0.4))

    assert reverse.ares(basis="reduced") == pytest.approx(forward.ares(basis="reduced"))
    assert reverse.pressure.to("pascal").magnitude == pytest.approx(
        forward.pressure.to("pascal").magnitude
    )


def test_pressure_identity__matches_compressibility_definition() -> None:
    model = Mixture(Parameters.from_mapping(NEUTRAL_BINARY))
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
    state = Mixture(Parameters.from_mapping(NEUTRAL_BINARY)).state(
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


@pytest.mark.parametrize(
    ("c_shell", "c_dielectric"),
    ((0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)),
    ids=("uncorrected", "shell-only", "shell-and-dielectric", "redundant"),
)
def test_born_corrections_match_independent_nested_equations(
    c_shell: float, c_dielectric: float
) -> None:
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

    try:
        state = _state(mapping, x=fractions, density=40000.0)
    except ParameterError:
        state = None
    if (c_shell, c_dielectric) == (0.0, 1.0):
        assert state is None
        return
    assert state is not None

    mixed_solvation = math.fsum(
        fraction * factor for fraction, factor in zip(fractions, solvation, strict=True)
    )
    bulk_permittivity = 78.09 / (1.0 + 7.01 * math.fsum(fractions[1:]))
    radial_sum = 0.0
    for fraction, charge, diameter in zip(fractions, charges, diameters, strict=True):
        if charge == 0:
            continue
        born_diameter = diameter * 1.0e-10
        shell_diameter = born_diameter * (1.0 + c_shell * (mixed_solvation - 1.0) / abs(charge))
        radial_sum += (
            fraction
            * charge
            * charge
            * (
                (1.0 - 1.0 / bulk_permittivity) / shell_diameter
                + c_dielectric * (1.0 - 1.0 / 8.0) * (1.0 / born_diameter - 1.0 / shell_diameter)
            )
        )
    expected = (
        -(ELEMENTARY_CHARGE**2)
        * radial_sum
        / (4.0 * math.pi * VACUUM_PERMITTIVITY * BOLTZMANN * temperature)
    )

    assert state.born == pytest.approx(expected, rel=3.0e-12, abs=2.0e-12)


def test_original_born_accepts_solvent_only_permittivity() -> None:
    """Bülow's original Born form uses only neutral-solvent permittivities."""
    mapping = deepcopy(AQUEOUS_ELECTROLYTE)
    next(
        family
        for family in mapping["model_families"]
        if family["kind"] == "permittivity"
    )["choice"] = "solvent-only"
    mapping["model_coefficients"] = [
        coefficient
        for coefficient in mapping["model_coefficients"]
        if coefficient["family"] != "ion_fraction_suppression_coefficient"
    ]
    electrolyte = next(
        family
        for family in mapping["model_families"]
        if family["kind"] == "electrolyte"
    )
    electrolyte["c_shell"] = 0.0
    electrolyte["c_dielectric"] = 0.0

    state = _state(mapping, x=(0.90, 0.05, 0.05), density=40000.0)

    assert state.bulk_relative_permittivity == pytest.approx(78.09, rel=2.0e-13)
    assert math.isfinite(state.debye_huckel)
    assert math.isfinite(state.born)


@pytest.mark.parametrize(
    "declaration",
    (
        "neutral-packing-diameter",
        "neutral-debye-huckel-diameter",
        "neutral-born-diameter",
        "nonpositive-packing-diameter",
        "nonpositive-debye-huckel-diameter",
        "nonpositive-born-diameter",
        "nonpositive-solvation-factor",
        "nonpositive-ionic-permittivity",
        "unused-charged-relative-permittivity",
    ),
)
def test_invalid_born_physical_declarations_do_not_resolve(declaration: str) -> None:
    mapping = deepcopy(AQUEOUS_ELECTROLYTE)
    if declaration == "unused-charged-relative-permittivity":
        ion = next(
            component
            for component in mapping["components"]
            if component["component_id"] == "sodium-cation"
        )
        ion["coefficients"].extend(
            _component("sodium-cation", relative_permittivity=10.0)["coefficients"]
        )
        ion["coefficients"].sort(key=lambda item: str(item["family"]))
        family = "relative_permittivity"
    else:
        family = declaration.removeprefix("neutral-").removeprefix("nonpositive-").replace("-", "_")
    if declaration.startswith("neutral-"):
        ion = next(
            component
            for component in mapping["components"]
            if component["component_id"] == "sodium-cation"
        )
        source = deepcopy(
            next(
                coefficient
                for coefficient in ion["coefficients"]
                if coefficient["family"] == family
            )
        )
        source["identity"] = f"component/water/{family}"
        solvent = next(
            component for component in mapping["components"] if component["component_id"] == "water"
        )
        solvent["coefficients"].append(source)
        solvent["coefficients"].sort(key=lambda item: str(item["family"]))
    elif declaration == "unused-charged-relative-permittivity":
        pass
    elif family == "ionic_permittivity":
        coefficient = next(
            item
            for item in mapping["model_coefficients"]
            if item["family"] == "ionic_region_relative_permittivity"
        )
        coefficient["value"]["magnitude"] = 0.0
    else:
        component_id = "water" if family == "solvation_factor" else "sodium-cation"
        component = next(
            item for item in mapping["components"] if item["component_id"] == component_id
        )
        coefficient = next(item for item in component["coefficients"] if item["family"] == family)
        coefficient["value"]["magnitude"] = 0.0

    try:
        model = Mixture(Parameters.from_mapping(mapping))
    except ParameterError:
        model = None
    assert model is None


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

    component_mixing = deepcopy(AQUEOUS_ELECTROLYTE)
    permittivity_family = next(
        family for family in component_mixing["model_families"] if family["kind"] == "permittivity"
    )
    permittivity_family["choice"] = "component-permittivity-mixing"
    component_mixing["model_coefficients"] = [
        coefficient
        for coefficient in component_mixing["model_coefficients"]
        if coefficient["family"] != "ion_fraction_suppression_coefficient"
    ]
    for component_id, permittivity in (
        ("sodium-cation", 10.0),
        ("chloride-anion", 5.0),
    ):
        component = next(
            item for item in component_mixing["components"] if item["component_id"] == component_id
        )
        component["coefficients"].extend(
            _component(component_id, relative_permittivity=permittivity)["coefficients"]
        )
        component["coefficients"].sort(key=lambda item: str(item["family"]))
    mixed = _state(component_mixing, x=fractions, density=density)
    assert mixed.bulk_relative_permittivity == pytest.approx(
        fractions[0] * 78.09 + fractions[1] * 10.0 + fractions[2] * 5.0,
        rel=2.0e-12,
    )


def test_component_permittivity__admits_exact_pure_ionic_support() -> None:
    components = [
        _component(
            component_id,
            molar_mass=molar_mass,
            charge_number=charge,
            segment_count=segment_count,
            segment_diameter=segment_diameter,
            dispersion_energy_over_k=dispersion_energy,
            relative_permittivity=relative_permittivity,
            packing_diameter=segment_diameter,
            debye_huckel_diameter=segment_diameter,
            born_diameter=segment_diameter,
            solvation_factor=1.0,
        )
        for (
            component_id,
            molar_mass,
            charge,
            segment_count,
            segment_diameter,
            dispersion_energy,
            relative_permittivity,
        ) in (
            ("cation", 0.111, 1, 1.5, 3.6, 210.0, 2.0),
            ("anion", 0.115, -1, 2.0, 3.8, 240.0, 4.0),
        )
    ]
    mapping = _document(
        "pure-ionic-fluid-fixture",
        components,
        pairs=[_pair("cation", "anion", k_ij=0.0)],
        model_families=_model_families(
            electrolyte="born",
            permittivity="component-permittivity-mixing",
            born_corrections=(0.0, 0.0),
        ),
        model_coefficients=[
            _model_coefficient("ionic_region_relative_permittivity", 1.0)
        ],
    )

    state = _state(mapping, x=(0.5, 0.5), density=5000.0)
    all_pair_mapping = deepcopy(mapping)
    electrolyte = next(
        family
        for family in all_pair_mapping["model_families"]
        if family["kind"] == "electrolyte"
    )
    electrolyte["ionic_dispersion"] = "all-ion-pairs"
    all_pair = _state(all_pair_mapping, x=(0.5, 0.5), density=5000.0)

    assert state.bulk_relative_permittivity == pytest.approx(3.0, rel=2.0e-14)
    assert math.isfinite(state.debye_huckel)
    assert math.isfinite(state.born)
    assert math.isfinite(state.pressure.to("pascal").magnitude)
    assert state.dispersion == pytest.approx(-1.1002931545314838, rel=2.0e-13)
    assert all_pair.dispersion == pytest.approx(-2.3098829509931313, rel=2.0e-13)
    assert all_pair.debye_huckel == pytest.approx(state.debye_huckel, rel=2.0e-14)


def test_excluded_same_sign_ion_pairs_resolve_to_unit_k_ij() -> None:
    mapping = _revised_debye_huckel(MULTISOLVENT_MULTIION)
    explicit = deepcopy(mapping)
    explicit["pairs"].extend(
        (
            _pair("cation", "divalent-cation", k_ij=0.25),
            _pair("anion", "iodide-anion", k_ij=0.35),
        )
    )
    composition = (0.445, 0.38, 0.05, 0.05, 0.025, 0.05)
    identity = "pair/anion/divalent-cation/k_ij"

    def evaluate(candidate):
        parameters = Parameters.from_mapping(candidate)
        model = Mixture(parameters)
        state = model.state(
            T=298.15 * u.kelvin,
            rho=2000.0 * u.mole / u.meter**3,
            x=composition,
        )
        block = model.observable_block(
            ActiveParameterSet(parameters, (identity,)),
            identity="eos.fixed-state",
            T=298.15 * u.kelvin,
            rho=2000.0 * u.mole / u.meter**3,
            x=composition,
        )
        return state, block

    resolved_state, resolved_block = evaluate(mapping)
    explicit_state, explicit_block = evaluate(explicit)
    all_pair = deepcopy(explicit)
    electrolyte = next(
        family for family in all_pair["model_families"] if family["kind"] == "electrolyte"
    )
    electrolyte["ionic_dispersion"] = "all-ion-pairs"
    all_pair_state = _state(all_pair, x=composition, density=2000.0)
    all_pair_parameters = Parameters.from_mapping(all_pair)
    all_pair_block = Mixture(all_pair_parameters).observable_block(
        ActiveParameterSet(all_pair_parameters, ("pair/cation/divalent-cation/k_ij",)),
        identity="eos.fixed-state",
        T=298.15 * u.kelvin,
        rho=2000.0 * u.mole / u.meter**3,
        x=composition,
    )

    assert resolved_state.dispersion == explicit_state.dispersion
    assert resolved_block.values == explicit_block.values
    assert resolved_block.jacobian == explicit_block.jacobian
    assert all_pair_state.dispersion < resolved_state.dispersion
    assert any(row[0] != 0.0 for row in all_pair_block.jacobian)

    explicit_parameters = Parameters.from_mapping(explicit)
    explicit_model = Mixture(explicit_parameters)
    with pytest.raises(EosError, match="inactive-by-selected-model"):
        explicit_model.observable_block(
            ActiveParameterSet(
                explicit_parameters,
                ("pair/cation/divalent-cation/k_ij",),
            ),
            identity="eos.fixed-state",
            T=298.15 * u.kelvin,
            rho=2000.0 * u.mole / u.meter**3,
            x=composition,
        )


def test_density_closure__returns_certified_single_root_and_slope() -> None:
    model = Pure(Parameters.from_mapping(neutral_pure("methane")))
    state = model.state(T=250.0 * u.kelvin, P=10.0 * u.bar)
    diagnostics = state.density_diagnostics

    assert diagnostics is not None
    assert diagnostics.stable is True
    assert abs(diagnostics.pressure_residual.to("pascal").magnitude) < 1.0e-3
    density = state.molar_density.to("mole / meter ** 3").magnitude
    step = density * 1.0e-4
    pressures = tuple(
        model.state(
            T=250.0 * u.kelvin,
            rho=value * u.mole / u.meter**3,
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
    model = Mixture(Parameters.from_mapping(NEUTRAL_BINARY))
    temperature = 250.0
    density = 8000.0
    fractions = (0.4, 0.6)
    volume = 1.0 / density
    state = _state(NEUTRAL_BINARY, x=fractions, density=density, temperature=temperature)
    expected = []
    step = 2.0e-5
    for component in range(2):

        def total_helmholtz(delta: float, component_index: int = component) -> float:
            amounts = [*fractions]
            amounts[component_index] += delta
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


def test_electrolyte_reference_matches_projection_identities() -> None:
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
    assert state.activity.reference_convention == "molality-infinite-dilution"
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
    assert identity.reference_convention == "molality-infinite-dilution"
    assert tuple(value.to("joule / mole").magnitude for value in identity.values) == (0.0, 0.0)

    reference_arguments = {
        "T": 298.15 * u.kelvin,
        "P": 1.0 * u.bar,
        "solvent_x": (1.0, 0.0, 0.0),
        "phase": "liquid",
        "ion_component_ids": ("sodium-cation", "chloride-anion"),
        "formula_stoichiometry": (1, 1),
        "reference_convention": "molality-infinite-dilution",
    }
    for invalid in (
        {"formula_stoichiometry": (1, 2)},
        {"reference_convention": "mole-fraction-pure-solvent"},
    ):
        try:
            rejected = model.reference(**(reference_arguments | invalid))
        except EosError:
            rejected = None
        assert rejected is None


def test_revised_debye_huckel_retains_opposite_charge_dispersion() -> None:
    """The revised Held form admits solvent-ion and cation-anion k_ij values."""

    mapping = _revised_debye_huckel(AQUEOUS_ELECTROLYTE)

    model = Mixture(Parameters.from_mapping(mapping))
    state = model.state(
        T=298.15 * u.kelvin,
        P=1.0 * u.bar,
        x=(0.98, 0.01, 0.01),
        phase="liquid",
    )

    assert state.fugacity is not None
    assert all(math.isfinite(value) for value in state.fugacity.ln_coefficient)

    suppressed = deepcopy(mapping)
    ion_pair = next(
        pair
        for pair in suppressed["pairs"]
        if {pair["component_id_a"], pair["component_id_b"]} == {"sodium-cation", "chloride-anion"}
    )
    ion_pair["coefficients"][0]["value"]["magnitude"] = 1.0
    suppressed_state = Mixture(Parameters.from_mapping(suppressed)).state(
        T=298.15 * u.kelvin,
        P=1.0 * u.bar,
        x=(0.98, 0.01, 0.01),
        phase="liquid",
    )

    assert suppressed_state.fugacity is not None
    assert state.fugacity.ln_coefficient != pytest.approx(
        suppressed_state.fugacity.ln_coefficient,
        rel=1.0e-8,
        abs=1.0e-10,
    )


def test_revised_debye_huckel_uses_declared_multisolvent_permittivity() -> None:
    """The revised no-Born model admits a specified neutral solvent mixture."""

    mapping = _revised_debye_huckel(MULTINEUTRAL_REACTIVE)

    composition = (0.50, 0.40, 0.05, 0.05)
    parameters = Parameters.from_mapping(mapping)
    model = Mixture(parameters)
    state = model.state(
        T=313.15 * u.kelvin,
        P=1.0 * u.bar,
        x=composition,
        phase="liquid",
    )
    expected = (
        composition[0] * 0.018 * 78.0 + composition[1] * 0.030 * 40.0
    ) / (composition[0] * 0.018 + composition[1] * 0.030)

    assert state.bulk_relative_permittivity == pytest.approx(expected, rel=2.0e-13)

    assert state.fugacity is not None
    assert all(math.isfinite(value) for value in state.fugacity.ln_coefficient)


def test_solvent_only_permittivity_tracks_temperature_and_excludes_neutral_solute() -> None:
    """Temperature correlations and dielectric participation are both explicit."""

    mapping = _revised_debye_huckel(MULTINEUTRAL_REACTIVE)
    mapping["components"].append(
        _component(
            "neutral-solute",
            molar_mass=0.044,
            charge_number=0,
            segment_count=1.4,
            segment_diameter=3.2,
            dispersion_energy_over_k=210.0,
        )
    )
    for component_id in ("solvent", "cosolvent", "cation", "anion"):
        mapping["pairs"].append(_pair(component_id, "neutral-solute", k_ij=0.0))

    def correlation(component_id: str, coefficients: tuple[float, ...]) -> dict[str, object]:
        correlation_id = (
            f"component/{component_id}/relative_permittivity/"
            "constant-plus-polynomial-in-temperature"
        )
        return {
            "correlation_id": correlation_id,
            "component_id": component_id,
            "family": "relative_permittivity",
            "form": "constant-plus-polynomial-in-temperature",
            "independent_variables": ["temperature"],
            "constant": {
                "identity": f"{correlation_id}/constant",
                "value": {"magnitude": coefficients[0], "unit": "dimensionless"},
            },
            "reference_temperature": {"magnitude": 300.0, "unit": "kelvin"},
            "terms": [
                {
                    "power": power,
                    "coefficient": {
                        "identity": f"{correlation_id}/term-{power - 1}/coefficient",
                        "value": {
                            "magnitude": coefficient,
                            "unit": f"1 / kelvin ** {power}",
                        },
                    },
                }
                for power, coefficient in enumerate(coefficients[1:], start=1)
            ],
            "provenance": deepcopy(mapping["correlations"][0]["provenance"]),
        }

    mapping["correlations"] = [
        correlation("solvent", (78.0, -0.3)),
        correlation("cosolvent", (40.0, -0.1, 0.001)),
    ]
    composition = (0.45, 0.35, 0.05, 0.05, 0.10)
    temperature = 320.0
    parameters = Parameters.from_mapping(mapping)
    model = Mixture(parameters)
    state = model.state(
        T=temperature * u.kelvin,
        rho=2000.0 * u.mole / u.meter**3,
        x=composition,
    )
    solvent_epsilon = 78.0 - 0.3 * 20.0
    cosolvent_epsilon = 40.0 - 0.1 * 20.0 + 0.001 * 20.0**2
    expected = (
        composition[0] * 0.018 * solvent_epsilon
        + composition[1] * 0.030 * cosolvent_epsilon
    ) / (composition[0] * 0.018 + composition[1] * 0.030)

    assert state.bulk_relative_permittivity == pytest.approx(expected, rel=2.0e-13)

    identity = (
        "component/solvent/relative_permittivity/"
        "constant-plus-polynomial-in-temperature/term-0/coefficient"
    )
    exact = model.observable_block(
        ActiveParameterSet(parameters, (identity,)),
        identity="eos.fixed-state",
        T=temperature * u.kelvin,
        rho=2000.0 * u.mole / u.meter**3,
        x=composition,
    )

    def values(slope: float) -> tuple[float, ...]:
        candidate = deepcopy(mapping)
        correlation_mapping = next(
            item
            for item in candidate["correlations"]
            if item["component_id"] == "solvent"
        )
        correlation_mapping["terms"][0]["coefficient"]["value"]["magnitude"] = slope
        candidate_state = Mixture(Parameters.from_mapping(candidate)).state(
            T=temperature * u.kelvin,
            rho=2000.0 * u.mole / u.meter**3,
            x=composition,
        )
        return (
            float(candidate_state.ares(basis="reduced")),
            float(candidate_state.pressure.to("pascal").magnitude),
            *candidate_state.residual_chemical_potential_over_rt,
        )

    step = 1.0e-6
    lower = values(-0.3 - step)
    upper = values(-0.3 + step)
    numerical = tuple(
        (high - low) / (2.0 * step) for low, high in zip(lower, upper, strict=True)
    )
    for row, expected_derivative in enumerate(numerical):
        assert exact.jacobian[row][0] == pytest.approx(
            expected_derivative, rel=6.0e-5, abs=4.0e-7
        )


def test_multisolvent_multiion_reference_is_component_order_invariant() -> None:
    solvent_by_id = {"solvent": 0.6, "cosolvent": 0.4}
    ion_pair = ("divalent-cation", "iodide-anion")
    stoichiometry = (1, 2)

    def evaluate(component_ids: tuple[str, ...]) -> tuple[tuple[float, ...], tuple[float, ...]]:
        parameters = Parameters.from_mapping(
            MULTISOLVENT_MULTIION,
            components=component_ids,
        )
        model = Mixture(parameters)
        reference = model.reference(
            T=298.15 * u.kelvin,
            P=1.0 * u.bar,
            solvent_x=tuple(solvent_by_id.get(component_id, 0.0) for component_id in component_ids),
            phase="liquid",
            ion_component_ids=ion_pair,
            formula_stoichiometry=stoichiometry,
            reference_convention="molality-infinite-dilution",
        )
        state = reference.state(molality=0.05 * u.mole / u.kilogram)
        assert state.activity is not None
        assert state.activity.reference_convention == "molality-infinite-dilution"
        return reference.infinite_dilution_log_fugacity_coefficient, state.activity.single_ion

    canonical_ids = tuple(
        component["component_id"] for component in MULTISOLVENT_MULTIION["components"]
    )
    permuted_ids = (
        "iodide-anion",
        "cosolvent",
        "divalent-cation",
        "solvent",
        "anion",
        "cation",
    )
    canonical = evaluate(canonical_ids)
    permuted = evaluate(permuted_ids)

    assert permuted[0] == pytest.approx(canonical[0], rel=2.0e-10, abs=2.0e-10)
    assert permuted[1] == pytest.approx(canonical[1], rel=2.0e-10, abs=2.0e-10)
