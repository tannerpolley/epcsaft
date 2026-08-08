from __future__ import annotations

from collections.abc import Callable
from copy import deepcopy

import pytest

from epcsaft import ActiveParameterSet, EosError, Mixture, Parameters, Pure
from epcsaft import unit_registry as u
from tests.eos._raw_cases import (
    _PROVENANCE,
    AQUEOUS_ELECTROLYTE,
    ASSOCIATING_PURE,
    MULTISOLVENT_MULTIION,
    POLAR_BINARY,
    POLAR_PURE,
    _component,
    _document,
    _model_families,
    _pair,
    neutral_pure,
)
from tests.equilibrium.parameter_dictionaries import PURE_SATURATION_PARAMETERS

Values = tuple[float, ...]


def _set_parameter(mapping: dict[str, object], identity: str, value: float) -> None:
    pending: list[object] = [mapping]
    while pending:
        item = pending.pop()
        if isinstance(item, dict):
            if item.get("identity") == identity:
                quantity = item["value"]
                assert isinstance(quantity, dict)
                quantity["magnitude"] = value
                return
            pending.extend(item.values())
        elif isinstance(item, list):
            pending.extend(item)
    raise AssertionError(f"parameter identity not found: {identity}")


def _parameter_value(mapping: dict[str, object], identity: str) -> float:
    parameters = Parameters.from_mapping(mapping)
    value = next(spec.value for spec in parameters.parameter_specs if spec.identity == identity)
    return float(getattr(value, "magnitude", value))


def _perturbed(mapping: dict[str, object], identity: str, value: float) -> dict[str, object]:
    changed = deepcopy(mapping)
    _set_parameter(changed, identity, value)
    return changed


def _central_difference(
    mapping: dict[str, object],
    identity: str,
    step: float,
    evaluate: Callable[[dict[str, object]], Values],
) -> Values:
    value = _parameter_value(mapping, identity)
    lower = evaluate(_perturbed(mapping, identity, value - step))
    upper = evaluate(_perturbed(mapping, identity, value + step))
    return tuple((high - low) / (2.0 * step) for low, high in zip(lower, upper, strict=True))


def _fixed_state_values(
    mapping: dict[str, object],
    temperature_k: float,
    density: float,
    composition: tuple[float, ...],
) -> Values:
    parameters = Parameters.from_mapping(mapping)
    model = Pure(parameters) if len(composition) == 1 else Mixture(parameters)
    state = (
        model.state(T=temperature_k * u.kelvin, rho=density * u.mole / u.meter**3)
        if isinstance(model, Pure)
        else model.state(
            T=temperature_k * u.kelvin,
            rho=density * u.mole / u.meter**3,
            x=composition,
        )
    )
    return (
        float(state.ares(basis="reduced")),
        float(state.pressure.to("pascal").magnitude),
        *state.residual_chemical_potential_over_rt,
    )


def _neutral_associating_mixture() -> dict[str, object]:
    mapping = deepcopy(ASSOCIATING_PURE)
    mapping["document_id"] = "neutral-associating-mixture-fixture"
    components = mapping["components"]
    pairs = mapping["pairs"]
    assert isinstance(components, list)
    assert isinstance(pairs, list)
    components.append(
        _component(
            "inert-fluid",
            segment_count=1.4,
            segment_diameter=3.4,
            dispersion_energy_over_k=180.0,
        )
    )
    pairs.append(_pair("associating-fluid", "inert-fluid", k_ij=0.04))
    return mapping


def _neutral_ternary_pairs() -> dict[str, object]:
    return _document(
        "neutral-ternary-pair-fixture",
        [
            _component(
                "alpha",
                segment_count=1.1,
                segment_diameter=3.1,
                dispersion_energy_over_k=180.0,
            ),
            _component(
                "beta",
                segment_count=1.6,
                segment_diameter=3.5,
                dispersion_energy_over_k=230.0,
            ),
            _component(
                "gamma",
                segment_count=2.0,
                segment_diameter=3.9,
                dispersion_energy_over_k=290.0,
            ),
        ],
        pairs=[
            _pair("alpha", "beta", k_ij=0.03, l_ij=-0.02),
            _pair("alpha", "gamma", k_ij=0.07, l_ij=0.04),
            _pair("beta", "gamma", k_ij=-0.05, l_ij=0.01),
        ],
    )


def _association_edge(
    component: str,
    site_a: str,
    site_b: str,
    energy: float,
    volume: float,
) -> dict[str, object]:
    endpoint_a, endpoint_b = sorted(((component, site_a), (component, site_b)))
    prefix = f"association/{endpoint_a[0]}/{endpoint_a[1]}/{endpoint_b[0]}/{endpoint_b[1]}"
    return {
        "endpoint_a": {"component_id": endpoint_a[0], "site_id": endpoint_a[1]},
        "endpoint_b": {"component_id": endpoint_b[0], "site_id": endpoint_b[1]},
        "energy_over_k": {
            "identity": f"{prefix}/energy_over_k",
            "value": {"magnitude": energy, "unit": "kelvin"},
        },
        "volume": {
            "identity": f"{prefix}/volume",
            "value": {"magnitude": volume, "unit": "dimensionless"},
        },
        "source": {"kind": "explicit", "provenance": [deepcopy(_PROVENANCE)]},
    }


def _combined_edge(
    endpoint_a: tuple[str, str],
    endpoint_b: tuple[str, str],
    input_a: tuple[str, str, str],
    input_b: tuple[str, str, str],
) -> dict[str, object]:
    self_edges = tuple(
        f"association/{component}/{site_a}/{component}/{site_b}"
        for component, site_a, site_b in (input_a, input_b)
    )
    return {
        "endpoint_a": {"component_id": endpoint_a[0], "site_id": endpoint_a[1]},
        "endpoint_b": {"component_id": endpoint_b[0], "site_id": endpoint_b[1]},
        "source": {
            "kind": "combining-rule",
            "rule_id": "arithmetic-energy-geometric-volume",
            "inputs": {
                "energy_over_k": [f"{prefix}/energy_over_k" for prefix in self_edges],
                "volume": [f"{prefix}/volume" for prefix in self_edges],
            },
            "provenance": [deepcopy(_PROVENANCE)],
        },
    }


def _multi_role_associating_binary() -> dict[str, object]:
    pair = _pair("alpha", "beta", k_ij=0.02, l_ij=0.01)
    coefficients = pair["coefficients"]
    assert isinstance(coefficients, list)
    for role_a, role_b, value in (
        ("acceptor", "donor", 0.08),
        ("effective", "effective", -0.04),
    ):
        identity = f"pair/alpha/beta/association/{role_a}/{role_b}/k_hb_ij"
        coefficients.append(
            {
                "identity": identity,
                "family": "k_hb_ij",
                "value": {"magnitude": value, "unit": "dimensionless"},
                "provenance": deepcopy(_PROVENANCE),
            }
        )
    sites = [
        {
            "component_id": component,
            "site_id": site,
            "site_role": role,
            "multiplicity": 1,
            "provenance": deepcopy(_PROVENANCE),
        }
        for component in ("alpha", "beta")
        for site, role in (
            ("a", "acceptor"),
            ("d", "donor"),
            ("e1", "effective"),
            ("e2", "effective"),
        )
    ]
    edges = [
        _association_edge("alpha", "a", "d", 2200.0, 0.030),
        _association_edge("beta", "a", "d", 1800.0, 0.025),
        _association_edge("alpha", "e1", "e2", 1500.0, 0.020),
        _association_edge("beta", "e1", "e2", 1300.0, 0.018),
        _combined_edge(
            ("alpha", "a"),
            ("beta", "d"),
            ("alpha", "a", "d"),
            ("beta", "a", "d"),
        ),
        _combined_edge(
            ("alpha", "d"),
            ("beta", "a"),
            ("alpha", "a", "d"),
            ("beta", "a", "d"),
        ),
        _combined_edge(
            ("alpha", "e1"),
            ("beta", "e2"),
            ("alpha", "e1", "e2"),
            ("beta", "e1", "e2"),
        ),
        _combined_edge(
            ("alpha", "e2"),
            ("beta", "e1"),
            ("alpha", "e1", "e2"),
            ("beta", "e1", "e2"),
        ),
    ]
    return _document(
        "multi-role-associating-binary-fixture",
        [
            _component(
                "alpha",
                segment_count=1.3,
                segment_diameter=3.2,
                dispersion_energy_over_k=210.0,
            ),
            _component(
                "beta",
                segment_count=1.7,
                segment_diameter=3.6,
                dispersion_energy_over_k=250.0,
            ),
        ],
        pairs=[pair],
        model_families=_model_families(association="general-site"),
        topology={"presets": [], "sites": sites, "edges": edges},
    )


@pytest.mark.parametrize(
    ("mapping", "temperature_k", "density", "composition", "identities", "steps"),
    (
        (
            neutral_pure("propane"),
            300.0,
            4000.0,
            (1.0,),
            (
                "component/propane/segment_count",
                "component/propane/segment_diameter",
                "component/propane/dispersion_energy_over_k",
            ),
            (1.0e-5, 1.0e-5, 1.0e-3),
        ),
        (
            POLAR_PURE,
            315.0,
            4200.0,
            (1.0,),
            (
                "component/polar-fluid/segment_count",
                "component/polar-fluid/segment_diameter",
                "component/polar-fluid/dispersion_energy_over_k",
            ),
            (1.0e-5, 1.0e-5, 1.0e-3),
        ),
        (
            _neutral_associating_mixture(),
            300.0,
            3500.0,
            (0.6, 0.4),
            (
                "association/associating-fluid/a/associating-fluid/b/energy_over_k",
                "association/associating-fluid/a/associating-fluid/b/volume",
            ),
            (1.0e-3, 1.0e-6),
        ),
        (
            AQUEOUS_ELECTROLYTE,
            298.15,
            2000.0,
            (0.98, 0.01, 0.01),
            (
                "association/water/a/water/b/energy_over_k",
                "association/water/a/water/b/volume",
            ),
            (1.0e-3, 1.0e-6),
        ),
    ),
    ids=(
        "base",
        "polar-fixed-moment",
        "neutral-mixture-association",
        "electrolyte-temperature-dependent-association",
    ),
)
def test_fixed_state_selected_parameter_derivatives_match_independent_differences(
    mapping: dict[str, object],
    temperature_k: float,
    density: float,
    composition: tuple[float, ...],
    identities: tuple[str, ...],
    steps: tuple[float, ...],
) -> None:
    parameters = Parameters.from_mapping(mapping)
    active = ActiveParameterSet(parameters, identities)
    model = Pure(parameters) if len(composition) == 1 else Mixture(parameters)
    exact = (
        model.observable_block(
            active,
            identity="eos.fixed-state",
            T=temperature_k * u.kelvin,
            rho=density * u.mole / u.meter**3,
        )
        if isinstance(model, Pure)
        else model.observable_block(
            active,
            identity="eos.fixed-state",
            T=temperature_k * u.kelvin,
            rho=density * u.mole / u.meter**3,
            x=composition,
        )
    )
    for column, (identity, step) in enumerate(zip(identities, steps, strict=True)):
        numerical = _central_difference(
            mapping,
            identity,
            step,
            lambda candidate: _fixed_state_values(candidate, temperature_k, density, composition),
        )
        for row, expected in enumerate(numerical):
            assert exact.jacobian[row][column] == pytest.approx(expected, rel=3.0e-5, abs=2.0e-7)


def test_canonical_kij_coordinate_is_exact_and_component_order_independent() -> None:
    """Falsify the selected k_ij Jacobian and unordered-pair identity together."""

    identity = "pair/dipolar/quadrupolar/k_ij"
    trial = 0.12
    step = 1.0e-5
    parameters = Parameters.from_mapping(POLAR_BINARY)
    active = ActiveParameterSet(parameters, (identity,), values=(trial,))
    exact = Mixture(parameters).observable_block(
        active,
        identity="eos.fixed-state",
        T=300.0 * u.kelvin,
        rho=5000.0 * u.mole / u.meter**3,
        x=(0.4, 0.6),
    )
    trial_mapping = _perturbed(POLAR_BINARY, identity, trial)
    numerical = _central_difference(
        trial_mapping,
        identity,
        step,
        lambda candidate: _fixed_state_values(candidate, 300.0, 5000.0, (0.4, 0.6)),
    )
    # h=1e-5 balances centered O(h^2) truncation with double-precision
    # cancellation; 3e-5 relative tolerance spans reduced and pascal rows.
    for row, expected in enumerate(numerical):
        assert exact.jacobian[row][0] == pytest.approx(expected, rel=3.0e-5, abs=2.0e-7)

    reordered_parameters = Parameters.from_mapping(
        POLAR_BINARY, components=("quadrupolar", "dipolar")
    )
    reordered = Mixture(reordered_parameters).observable_block(
        ActiveParameterSet(reordered_parameters, (identity,), values=(trial,)),
        identity="eos.fixed-state",
        T=300.0 * u.kelvin,
        rho=5000.0 * u.mole / u.meter**3,
        x=(0.6, 0.4),
    )
    # Reordering changes only component-labelled rows; the 2e-10 tolerance is
    # tight enough to expose an incorrect directed write while allowing AD roundoff.
    for canonical_row, reordered_row in enumerate((0, 1, 3, 2)):
        assert reordered.values[reordered_row] == pytest.approx(
            exact.values[canonical_row], rel=2.0e-12, abs=2.0e-12
        )
        assert reordered.jacobian[reordered_row][0] == pytest.approx(
            exact.jacobian[canonical_row][0], rel=2.0e-10, abs=2.0e-10
        )


def test_mixed_ternary_pair_coordinates_are_independent_and_ordered() -> None:
    """Falsify reciprocal mixed-pair columns and their declared ordering."""

    mapping = _neutral_ternary_pairs()
    identities = (
        "pair/alpha/beta/k_ij",
        "pair/beta/gamma/k_ij",
        "pair/alpha/gamma/l_ij",
        "pair/beta/gamma/l_ij",
    )
    composition = (0.2, 0.3, 0.5)
    parameters = Parameters.from_mapping(mapping)
    exact = Mixture(parameters).observable_block(
        ActiveParameterSet(parameters, identities),
        identity="eos.fixed-state",
        T=315.0 * u.kelvin,
        rho=4300.0 * u.mole / u.meter**3,
        x=composition,
    )

    # h=1e-5 balances centered O(h^2) error and cancellation across the
    # dimensionless Helmholtz/chemical rows and the pressure row in pascals.
    # A mismatch means a selected coordinate wrote the wrong pair or column.
    for column, identity in enumerate(identities):
        numerical = _central_difference(
            mapping,
            identity,
            1.0e-5,
            lambda candidate: _fixed_state_values(candidate, 315.0, 4300.0, composition),
        )
        for row, expected in enumerate(numerical):
            assert exact.jacobian[row][column] == pytest.approx(expected, rel=4.0e-5, abs=3.0e-7)

    reordered_parameters = Parameters.from_mapping(mapping, components=("gamma", "alpha", "beta"))
    reordered = Mixture(reordered_parameters).observable_block(
        ActiveParameterSet(reordered_parameters, tuple(reversed(identities))),
        identity="eos.fixed-state",
        T=315.0 * u.kelvin,
        rho=4300.0 * u.mole / u.meter**3,
        x=(0.5, 0.2, 0.3),
    )
    for canonical_row, reordered_row in enumerate((0, 1, 3, 4, 2)):
        assert reordered.values[reordered_row] == pytest.approx(
            exact.values[canonical_row], rel=2.0e-12, abs=2.0e-12
        )
        for canonical_column in range(len(identities)):
            reordered_column = len(identities) - canonical_column - 1
            assert reordered.jacobian[reordered_row][reordered_column] == pytest.approx(
                exact.jacobian[canonical_row][canonical_column],
                rel=2.0e-10,
                abs=2.0e-10,
            )


def test_role_aware_cross_association_corrections_apply_exact_chain() -> None:
    """Falsify role-aware energy correction and its implicit EOS chain rule."""

    mapping = _multi_role_associating_binary()
    identities = (
        "pair/alpha/beta/k_ij",
        "pair/alpha/beta/l_ij",
        "pair/alpha/beta/association/acceptor/donor/k_hb_ij",
        "pair/alpha/beta/association/effective/effective/k_hb_ij",
    )
    composition = (0.45, 0.55)
    parameters = Parameters.from_mapping(mapping)
    # Independent combining-rule oracle: energies are kelvin and volumes are
    # dimensionless; 1e-13 relative tolerance covers Decimal-to-double
    # conversion only. A failure means k_HB altered volume or missed an
    # eligible multiplicity-equivalent edge.
    cross_edges = {
        (edge.endpoint_a, edge.endpoint_b): (
            float(edge.energy_over_k.to("kelvin").magnitude),
            float(edge.volume),
        )
        for edge in parameters.topology.edges
        if edge.endpoint_a[0] != edge.endpoint_b[0]
    }
    for endpoints in (
        (("alpha", "a"), ("beta", "d")),
        (("alpha", "d"), ("beta", "a")),
    ):
        energy, volume = cross_edges[endpoints]
        assert energy == pytest.approx(2000.0 * (1.0 - 0.08), rel=1.0e-13, abs=1.0e-12)
        assert volume == pytest.approx((0.030 * 0.025) ** 0.5, rel=1.0e-13, abs=1.0e-15)
    for endpoints in (
        (("alpha", "e1"), ("beta", "e2")),
        (("alpha", "e2"), ("beta", "e1")),
    ):
        energy, volume = cross_edges[endpoints]
        assert energy == pytest.approx(1400.0 * (1.0 - (-0.04)), rel=1.0e-13, abs=1.0e-12)
        assert volume == pytest.approx((0.020 * 0.018) ** 0.5, rel=1.0e-13, abs=1.0e-15)

    exact = Mixture(parameters).observable_block(
        ActiveParameterSet(parameters, identities),
        identity="eos.fixed-state",
        T=310.0 * u.kelvin,
        rho=2800.0 * u.mole / u.meter**3,
        x=composition,
    )

    # The same h=1e-5/tolerance pair spans reduced Helmholtz and chemical rows
    # plus pressure in pascals. A mismatch means the taped correction omitted
    # an edge or failed the site-fraction implicit chain rule.
    for column, identity in enumerate(identities):
        numerical = _central_difference(
            mapping,
            identity,
            1.0e-5,
            lambda candidate: _fixed_state_values(candidate, 310.0, 2800.0, composition),
        )
        for row, expected in enumerate(numerical):
            assert exact.jacobian[row][column] == pytest.approx(expected, rel=4.0e-5, abs=3.0e-7)

    reordered_parameters = Parameters.from_mapping(mapping, components=("beta", "alpha"))
    reordered = Mixture(reordered_parameters).observable_block(
        ActiveParameterSet(reordered_parameters, tuple(reversed(identities))),
        identity="eos.fixed-state",
        T=310.0 * u.kelvin,
        rho=2800.0 * u.mole / u.meter**3,
        x=(0.55, 0.45),
    )
    for canonical_row, reordered_row in enumerate((0, 1, 3, 2)):
        assert reordered.values[reordered_row] == pytest.approx(
            exact.values[canonical_row], rel=2.0e-12, abs=2.0e-12
        )
        for canonical_column in range(len(identities)):
            reordered_column = len(identities) - canonical_column - 1
            assert reordered.jacobian[reordered_row][reordered_column] == pytest.approx(
                exact.jacobian[canonical_row][canonical_column],
                rel=2.0e-10,
                abs=2.0e-10,
            )


@pytest.mark.parametrize(
    ("mapping", "temperature_k", "identities", "steps"),
    (
        (
            PURE_SATURATION_PARAMETERS["propane"],
            300.0,
            (
                "component/propane/segment_count",
                "component/propane/segment_diameter",
                "component/propane/dispersion_energy_over_k",
            ),
            (1.0e-5, 1.0e-5, 1.0e-3),
        ),
        (
            ASSOCIATING_PURE,
            300.0,
            (
                "association/associating-fluid/a/associating-fluid/b/energy_over_k",
                "association/associating-fluid/a/associating-fluid/b/volume",
            ),
            (1.0e-3, 1.0e-6),
        ),
    ),
    ids=("base", "association"),
)
def test_saturation_implicit_parameter_derivatives_match_resolved_boundaries(
    mapping: dict[str, object],
    temperature_k: float,
    identities: tuple[str, ...],
    steps: tuple[float, ...],
) -> None:
    parameters = Parameters.from_mapping(mapping)
    exact = Pure(parameters).observable_block(
        ActiveParameterSet(parameters, identities),
        identity="eos.pure-saturation",
        T=temperature_k * u.kelvin,
    )
    branches = {branch.identity: branch for branch in exact.saturation_branches}
    for row, identity in ((1, "vapor"), (2, "liquid")):
        branch = branches[identity]
        density = exact.values[row]
        assert density is not None
        assert branch.density_bracket_mol_per_m3[0] < density < branch.density_bracket_mol_per_m3[1]
        assert branch.pressure_relative_residual is not None
        assert branch.pressure_relative_residual <= 1.0e-8
        assert branch.residual_density_derivative_m3_per_mol is not None
        assert branch.residual_density_derivative_m3_per_mol > 0.0

    def values(candidate: dict[str, object]) -> Values:
        result = Pure(Parameters.from_mapping(candidate)).saturation(T=temperature_k * u.kelvin)
        return (
            result.saturation_pressure_pa,
            result.vapor.molar_density_mol_m3,
            result.liquid.molar_density_mol_m3,
        )

    for column, (identity, step) in enumerate(zip(identities, steps, strict=True)):
        numerical = _central_difference(mapping, identity, step, values)
        for row, expected in enumerate(numerical):
            assert exact.jacobian[row][column] == pytest.approx(expected, rel=3.0e-5, abs=2.0e-6)


def test_saturation_trial_geometry_uses_its_active_density_domain() -> None:
    mapping = PURE_SATURATION_PARAMETERS["propane"]
    identity = "component/propane/segment_diameter"
    trial = 2.5
    parameters = Parameters.from_mapping(mapping)
    exact = Pure(parameters).observable_block(
        ActiveParameterSet(parameters, (identity,), values=(trial,)),
        identity="eos.pure-saturation",
        T=300.0 * u.kelvin,
    )
    trial_mapping = _perturbed(mapping, identity, trial)

    def values(candidate: dict[str, object]) -> Values:
        result = Pure(Parameters.from_mapping(candidate)).saturation(T=300.0 * u.kelvin)
        return (
            result.saturation_pressure_pa,
            result.vapor.molar_density_mol_m3,
            result.liquid.molar_density_mol_m3,
        )

    numerical = _central_difference(trial_mapping, identity, 1.0e-5, values)
    assert exact.values[:3] == pytest.approx(values(trial_mapping), rel=2.0e-12, abs=2.0e-12)
    for row, expected in enumerate(numerical):
        derivative = exact.jacobian[row]
        assert derivative is not None
        assert derivative[0] == pytest.approx(expected, rel=3.0e-5, abs=2.0e-6)


def test_fixed_pressure_fugacity_derivatives_include_density_closure() -> None:
    mapping = neutral_pure("propane")
    identities = (
        "component/propane/segment_diameter",
        "component/propane/dispersion_energy_over_k",
    )
    trial_values = (5.0, 220.0)
    trial_mapping = _perturbed(
        _perturbed(mapping, identities[0], trial_values[0]),
        identities[1],
        trial_values[1],
    )
    parameters = Parameters.from_mapping(mapping)
    exact = Pure(parameters).observable_block(
        ActiveParameterSet(parameters, identities, values=trial_values),
        identity="eos.fixed-pressure-fugacity",
        T=300.0 * u.kelvin,
        P=1.0e5 * u.pascal,
        phase="vapor",
    )

    def value(candidate: dict[str, object]) -> Values:
        state = Pure(Parameters.from_mapping(candidate)).state(
            T=300.0 * u.kelvin, P=1.0e5 * u.pascal, phase="vapor"
        )
        assert state.fugacity is not None
        return state.fugacity.ln_coefficient

    for column, (identity, step) in enumerate(zip(identities, (1.0e-5, 1.0e-3), strict=True)):
        numerical = _central_difference(trial_mapping, identity, step, value)
        assert exact.jacobian[0][column] == pytest.approx(numerical[0], rel=3.0e-5, abs=2.0e-8)
    assert exact.values == pytest.approx(value(trial_mapping), rel=2.0e-12, abs=2.0e-12)


def test_active_geometry_rejects_close_packing_before_taping() -> None:
    mapping = neutral_pure("propane")
    identity = "component/propane/segment_diameter"
    parameters = Parameters.from_mapping(mapping)
    with pytest.raises(EosError, match="active packing fraction is outside") as captured:
        Pure(parameters).observable_block(
            ActiveParameterSet(parameters, (identity,), values=(7.4,)),
            identity="eos.fixed-state",
            T=300.0 * u.kelvin,
            rho=4000.0 * u.mole / u.meter**3,
        )
    assert captured.value.code == "unsupported-parameter-derivative"


def test_mixed_aqueous_active_set_matches_differences_and_preserves_axes() -> None:
    """Falsify the generic tape across scopes and Born input families."""

    identities = (
        "component/solvent/segment_count",
        "component/iodide-anion/dispersion_energy_over_k",
        "component/solvent/segment_diameter",
        "association/solvent/a/solvent/b/energy_over_k",
        "pair/iodide-anion/solvent/k_ij",
        "component/solvent/relative_permittivity",
        "component/divalent-cation/born_diameter",
        "component/cosolvent/solvation_factor",
        "model/relative_permittivity/ionic_region_relative_permittivity",
    )
    steps = (
        1.0e-5,
        1.0e-3,
        1.0e-5,
        1.0e-3,
        1.0e-5,
        1.0e-4,
        1.0e-5,
        1.0e-5,
        1.0e-5,
    )
    composition = (0.60, 0.33, 0.02, 0.02, 0.01, 0.02)
    parameters = Parameters.from_mapping(MULTISOLVENT_MULTIION)
    exact = Mixture(parameters).observable_block(
        ActiveParameterSet(parameters, identities),
        identity="eos.fixed-state",
        T=298.15 * u.kelvin,
        rho=2000.0 * u.mole / u.meter**3,
        x=composition,
    )

    # Each column is independently reconstructed through Parameters. The
    # scale-specific steps balance centered O(h^2) truncation and cancellation;
    # the tolerance spans reduced rows and pressure in pascals.
    for column, (identity, step) in enumerate(zip(identities, steps, strict=True)):
        numerical = _central_difference(
            MULTISOLVENT_MULTIION,
            identity,
            step,
            lambda candidate: _fixed_state_values(candidate, 298.15, 2000.0, composition),
        )
        for row, expected in enumerate(numerical):
            assert exact.jacobian[row][column] == pytest.approx(expected, rel=5.0e-5, abs=3.0e-7)

    reversed_identities = tuple(reversed(identities))
    reordered_parameters = Parameters.from_mapping(
        MULTISOLVENT_MULTIION,
        components=(
            "iodide-anion",
            "cosolvent",
            "divalent-cation",
            "solvent",
            "anion",
            "cation",
        ),
    )
    reordered = Mixture(reordered_parameters).observable_block(
        ActiveParameterSet(reordered_parameters, reversed_identities),
        identity="eos.fixed-state",
        T=298.15 * u.kelvin,
        rho=2000.0 * u.mole / u.meter**3,
        x=(0.02, 0.33, 0.01, 0.60, 0.02, 0.02),
    )
    for canonical_row, reordered_row in enumerate((0, 1, 5, 3, 7, 6, 4, 2)):
        assert reordered.values[reordered_row] == pytest.approx(
            exact.values[canonical_row], rel=2.0e-12, abs=2.0e-12
        )
        for canonical_column in range(len(identities)):
            reordered_column = len(identities) - canonical_column - 1
            assert reordered.jacobian[reordered_row][reordered_column] == pytest.approx(
                exact.jacobian[canonical_row][canonical_column],
                rel=2.0e-10,
                abs=2.0e-10,
            )
