from __future__ import annotations

import copy
import math

import pytest

import epcsaft
from epcsaft import equilibrium
from epcsaft.equilibrium import _equilibrium
from epcsaft.equilibrium._api import _chemical_equilibrium_impl
from tests.eos._raw_cases import NEUTRAL_TERNARY_IMMISCIBLE, STRONGLY_ASSOCIATING_BINARY

from .chemical_equilibrium_cases import base_system as _base_system
from .chemical_equilibrium_cases import bind_records as _bind_record
from .chemical_equilibrium_cases import typed_problem as _typed_problem

_BELOV_ELEMENT_POTENTIALS = (
    -10.0776517753450153959143507476192439107332024990741,
    44.8480568890243540815059574241188924984510845101041,
)
_BELOV_SOURCE_GIBBS = (
    -7.072468897041624,
    -28.765987200068952,
    -26.51300119486301,
    21.44955454258845,
    21.26761422175228,
    -34.35125913287487,
    -55.36542486137457,
    -17.10224518423043,
)


def _belov_aristova_gas_system() -> tuple[dict[str, object], tuple[float, ...]]:
    temperature_k = 2000.0
    pressure_pa = 100_000.0
    standard_pressure_pa = 101_325.0
    reactions = (
        (-2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        (-3.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0, -2.0, 1.0, 0.0, 0.0, 0.0),
        (-1.0, 0.0, 0.0, -1.0, 0.0, 1.0, 0.0, 0.0),
        (-2.0, 0.0, 0.0, -1.0, 0.0, 0.0, 1.0, 0.0),
        (-1.0, 0.0, 0.0, -2.0, 0.0, 0.0, 0.0, 1.0),
    )
    concentration_shift = math.log(8.31446261815324 * temperature_k / standard_pressure_pa)
    manufactured_gibbs = tuple(value + concentration_shift for value in _BELOV_SOURCE_GIBBS)
    ln_k = tuple(
        -sum(row[index] * manufactured_gibbs[index] for index in range(8)) for row in reactions
    )
    spec: dict[str, object] = {
        "species_ids": ("O", "O2", "O3", "C", "C2", "CO", "CO2", "C2O"),
        "charges": (0,) * 8,
        "provider_fingerprint": (
            "sha256:2f24904e5db64e96aa2ee96c53bd80863c9518b6c2ddd914e266181fe7459f84"
        ),
        "molar_masses_kg_per_mol": (
            0.016,
            0.032,
            0.048,
            0.012,
            0.024,
            0.028,
            0.044,
            0.040,
        ),
        "balance_matrix": (
            (0.0, 0.0, 0.0, 1.0, 2.0, 1.0, 1.0, 2.0),
            (1.0, 2.0, 3.0, 0.0, 0.0, 1.0, 2.0, 1.0),
        ),
        "reaction_matrix": reactions,
        "feed_amounts": (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0),
        "ln_k": ln_k,
        "temperature_k": temperature_k,
        "pressure_pa": pressure_pa,
    }
    _bind_record(spec)
    expected = (
        5.969938280702114e-17,
        5.249059580600665e-27,
        1.838054910648896e-47,
        1.7496236329525488e-5,
        4.996074090543112e-1,
        9.992323142113352e-1,
        4.453885591412746e-11,
        7.676856995870245e-4,
    )
    return spec, expected


def _manufactured_solve(
    spec: dict[str, object], options: dict[str, object] | None = None
) -> equilibrium.ChemicalEquilibriumResult:
    options = options or {}
    return _chemical_equilibrium_impl(
        equilibrium.IdealGasPhase(
            model_fingerprint=spec["provider_fingerprint"],
            reference_id="provider-helmholtz-coordinate-basis",
        ),
        spec["temperature_k"] * epcsaft.unit_registry.kelvin,
        spec["pressure_pa"] * epcsaft.unit_registry.pascal,
        _typed_problem(
            spec,
            minimum_amount_mol=options.get("trace_floor", 1.0e-12),
        ),
    )


@pytest.mark.parametrize(
    ("spec", "expected_amounts", "expected_volume"),
    (
        (_base_system(), (0.2, 0.8), None),
        (
            {
                **_base_system(),
                "balance_matrix": ((1.0, 2.0),),
                "reaction_matrix": ((-2.0, 1.0),),
                "molar_masses_kg_per_mol": (1.0, 2.0),
                "feed_amounts": (2.0, 0.0),
                "ln_k": (math.log(0.75),),
                "temperature_k": 300.0,
                "pressure_pa": 8.31446261815324 * 300.0,
            },
            (1.0, 0.5),
            1.5,
        ),
    ),
)
def test_manufactured_ideal_reactions_match_independent_analytic_states(
    spec: dict[str, object],
    expected_amounts: tuple[float, ...],
    expected_volume: float | None,
) -> None:
    _bind_record(spec)
    result = _manufactured_solve(spec)

    pressure_over_rt = spec["pressure_pa"] / (8.31446261815324 * spec["temperature_k"])
    volume = expected_volume or sum(expected_amounts) / pressure_over_rt
    diagnostics = result.diagnostics
    assert result.amounts_mol == pytest.approx(expected_amounts, rel=2.0e-8, abs=2.0e-10)
    assert result.volume_m3 == pytest.approx(volume, rel=2.0e-8)
    assert diagnostics.balance_inf_norm <= 1.0e-9
    assert diagnostics.charge_inf_norm <= 1.0e-12
    assert diagnostics.pressure_relative_residual <= 1.0e-8
    assert diagnostics.reaction_affinity_inf_norm <= 1.0e-7
    assert max(map(abs, diagnostics.reaction_affinity_residuals)) == pytest.approx(
        diagnostics.reaction_affinity_inf_norm
    )
    assert diagnostics.kkt_stationarity_inf_norm <= 1.0e-7
    assert max(map(abs, diagnostics.physical_stationarity_residuals)) == pytest.approx(
        diagnostics.kkt_stationarity_inf_norm
    )
    physical_rows, physical_columns = diagnostics.physical_constraint_shape
    derivative_dimension = len(diagnostics.derivative_coordinate_order)
    assert physical_columns == len(diagnostics.physical_stationarity_residuals)
    assert len(diagnostics.physical_equality_multipliers) == physical_rows
    assert len(diagnostics.physical_constraint_jacobian) == (physical_rows * physical_columns)
    physical_chart_rows, physical_chart_columns = diagnostics.physical_to_chart_shape
    assert physical_chart_rows == len(diagnostics.physical_lagrangian_gradient)
    assert physical_chart_columns == derivative_dimension
    assert len(diagnostics.physical_to_chart_jacobian) == (
        physical_chart_rows * physical_chart_columns
    )
    assert diagnostics.chart_physical_pullback_residual_inf_norm is not None
    assert diagnostics.chart_physical_pullback_residual_inf_norm <= 1.0e-10
    assert diagnostics.complementarity_inf_norm <= 1.0e-7
    assert diagnostics.objective_gradient_shape == (derivative_dimension,)
    assert diagnostics.constraint_jacobian_shape == (
        len(diagnostics.constraint_values),
        derivative_dimension,
    )
    assert diagnostics.lagrangian_hessian_shape == (
        derivative_dimension,
        derivative_dimension,
    )
    assert diagnostics.derivative_objective_basis == (
        "dimensionless_fixed_TP_A_plus_PV_plus_reference_over_RT"
    )
    assert diagnostics.derivative_constraint_basis == (
        "ordered_per_row_in_derivative_constraint_order"
    )
    assert len(diagnostics.derivative_constraint_order) == len(diagnostics.constraint_values)
    assert diagnostics.chart_stationarity_inf_norm <= 1.0e-7
    assert diagnostics.derivative_check_step is not None
    assert diagnostics.derivative_check_step > 0.0
    assert diagnostics.objective_gradient_check_relative_error is not None
    assert diagnostics.objective_gradient_check_relative_error <= 1.0e-5
    assert diagnostics.physical_objective_gradient_check_relative_error is not None
    assert diagnostics.physical_objective_gradient_check_relative_error <= 1.0e-5
    assert len(diagnostics.physical_objective_gradient) == (diagnostics.physical_to_chart_shape[0])
    assert diagnostics.constraint_jacobian_check_relative_error is not None
    assert diagnostics.constraint_jacobian_check_relative_error <= 1.0e-5
    assert diagnostics.lagrangian_hessian_check_relative_error is not None
    assert diagnostics.lagrangian_hessian_check_relative_error <= 2.0e-4
    assert diagnostics.reduced_hessian_check_relative_error is not None
    assert diagnostics.reduced_hessian_check_relative_error <= 2.0e-4
    assert diagnostics.derivative_check_worst_entry
    assert diagnostics.derivative_check_worst_relative_error is not None
    assert diagnostics.derivative_check_worst_analytic_value is not None
    assert diagnostics.derivative_check_worst_finite_difference_value is not None
    assert diagnostics.derivative_check_worst_step is not None
    assert diagnostics.kkt_root_shape == (
        diagnostics.kkt_dimension,
        diagnostics.kkt_dimension,
    )
    assert len(diagnostics.kkt_root_jacobian) == diagnostics.kkt_dimension**2
    assert diagnostics.kkt_root_jacobian_check_relative_error is not None
    assert diagnostics.kkt_root_jacobian_check_relative_error <= 2.0e-4
    assert sum(diagnostics.reduced_hessian_raw_inertia) == len(
        diagnostics.reduced_hessian_eigenvalues
    )
    reduced_dimension, full_dimension = diagnostics.reduced_hessian_nullspace_shape
    assert full_dimension == derivative_dimension
    assert reduced_dimension**2 == len(diagnostics.reduced_hessian)
    assert len(diagnostics.reduced_hessian_nullspace_basis) == (reduced_dimension * full_dimension)
    for row in range(reduced_dimension):
        for column in range(reduced_dimension):
            projected = math.fsum(
                diagnostics.reduced_hessian_nullspace_basis[row * full_dimension + left]
                * diagnostics.lagrangian_hessian[left * full_dimension + right]
                * diagnostics.reduced_hessian_nullspace_basis[column * full_dimension + right]
                for left in range(full_dimension)
                for right in range(full_dimension)
            )
            assert projected == pytest.approx(
                diagnostics.reduced_hessian[row * reduced_dimension + column],
                abs=2.0e-10,
            )
    assert diagnostics.minimum_amount_mol is not None
    assert diagnostics.trace_floor_mol == 1.0e-12
    assert diagnostics.minimum_amount_mol > diagnostics.trace_floor_mol

    generic_spec = copy.deepcopy(spec)
    generic_spec["species_ids"] = ("solute", "solvent")
    generic = equilibrium.solve(
        epcsaft.Mixture(epcsaft.Parameters.from_mapping(STRONGLY_ASSOCIATING_BINARY)),
        equilibrium.GeneralReactiveEquilibriumProblem(
            "manufactured-ideal-reaction",
            equilibrium.Fixed(spec["temperature_k"] * epcsaft.unit_registry.kelvin),
            equilibrium.Fixed(spec["pressure_pa"] * epcsaft.unit_registry.pascal),
            (
                equilibrium.ReactivePhase(
                    "vapor",
                    "vapor",
                    "finite",
                    equilibrium.AllComponents(),
                    equilibrium.IdealGasModel("provider-helmholtz-coordinate-basis"),
                ),
            ),
            _typed_problem(generic_spec),
            ("vapor",),
            tuple(
                equilibrium.EquilibriumOutput(
                    f"amount-{index}",
                    "phase.amount",
                    "mole",
                    "species amount",
                    "vapor",
                    tuple(float(component == index) for component in range(2)),
                )
                for index in range(2)
            ),
        ),
    )
    assert generic.status == "evaluated"
    assert generic.values == pytest.approx(expected_amounts, rel=2.0e-8, abs=2.0e-10)
    evidence = dict(generic.evidence)
    assert evidence["balance_inf_norm"] <= 1.0e-9
    assert evidence["electroneutrality_inf"] <= 1.0e-12
    assert evidence["pressure_relative_inf_norm"] <= 1.0e-8
    assert evidence["reaction_affinity_inf_norm"] <= 1.0e-7
    assert evidence["kkt_stationarity_inf_norm"] <= 1.0e-7


def test_manufactured_reaction_reports_exact_conditioned_implicit_sensitivities() -> None:
    spec = _base_system()
    _bind_record(spec)
    native = _equilibrium._chemical_equilibrium(None, spec, None, None, 1.0e-12)

    sensitivities = native["sensitivities"]
    assert sensitivities["kkt_rank"] == sensitivities["kkt_dimension"] == 3
    assert sensitivities["condition_number_inf"] < 1.0e6

    amount_derivatives = sensitivities["amount_derivatives"]
    assert amount_derivatives[0] == pytest.approx((0.2, 0.8), abs=2.0e-9)
    assert amount_derivatives[1] == pytest.approx((-0.16, 0.16), abs=2.0e-9)
    assert amount_derivatives[2] == pytest.approx((0.0, 0.0), abs=2.0e-12)

    expected_volume = 8.31446261815324 * 350.0 / 200_000.0
    assert sensitivities["volume_derivatives"] == pytest.approx(
        (
            expected_volume,
            0.0,
            -expected_volume / 200_000.0,
        ),
        rel=2.0e-8,
        abs=2.0e-13,
    )

    steps = (1.0e-5, 1.0e-5, 5.0)
    for parameter, step in enumerate(steps):
        perturbed: list[dict[str, object]] = []
        for direction in (-1.0, 1.0):
            trial = copy.deepcopy(spec)
            if parameter == 0:
                trial["feed_amounts"] = (
                    1.0 + direction * step,
                    0.0,
                )
            elif parameter == 1:
                trial["ln_k"] = (
                    spec["ln_k"][0] + direction * step,  # type: ignore[index]
                )
            else:
                trial["pressure_pa"] = spec["pressure_pa"] + direction * step
            _bind_record(trial)
            perturbed.append(_equilibrium._chemical_equilibrium(None, trial, None, None, 1.0e-12))
        finite_difference_amounts = tuple(
            (perturbed[1]["amounts"][species] - perturbed[0]["amounts"][species]) / (2.0 * step)
            for species in range(2)
        )
        finite_difference_volume = (perturbed[1]["volume_m3"] - perturbed[0]["volume_m3"]) / (
            2.0 * step
        )
        assert amount_derivatives[parameter] == pytest.approx(
            finite_difference_amounts,
            rel=2.0e-6,
            abs=2.0e-10,
        )
        assert sensitivities["volume_derivatives"][parameter] == pytest.approx(
            finite_difference_volume,
            rel=2.0e-6,
            abs=2.0e-12,
        )


def test_manufactured_sensitivity_is_invariant_to_reaction_scaling() -> None:
    reaction_scale = 1.0e-11
    base = _base_system()
    _bind_record(base)
    scaled = copy.deepcopy(base)
    scaled["reaction_matrix"] = ((-reaction_scale, reaction_scale),)
    scaled["ln_k"] = (reaction_scale * math.log(4.0),)
    _bind_record(scaled)

    base_result = _equilibrium._chemical_equilibrium(None, base, None, None, 1.0e-12)[
        "sensitivities"
    ]
    scaled_result = _equilibrium._chemical_equilibrium(None, scaled, None, None, 1.0e-12)[
        "sensitivities"
    ]

    assert scaled_result["condition_number_inf"] == pytest.approx(
        base_result["condition_number_inf"], rel=2.0e-10
    )
    assert scaled_result["amount_derivatives"][0] == pytest.approx(
        base_result["amount_derivatives"][0], abs=2.0e-9
    )
    assert tuple(
        reaction_scale * value for value in scaled_result["amount_derivatives"][1]
    ) == pytest.approx(base_result["amount_derivatives"][1], abs=2.0e-9)
    assert scaled_result["amount_derivatives"][2] == pytest.approx(
        base_result["amount_derivatives"][2], abs=2.0e-12
    )


def test_manufactured_sensitivity_tracks_exact_species_permutation() -> None:
    base = _base_system()
    _bind_record(base)
    permutation = (1, 0)
    permuted = copy.deepcopy(base)
    for field in (
        "species_ids",
        "charges",
        "molar_masses_kg_per_mol",
        "feed_amounts",
    ):
        values = base[field]
        permuted[field] = tuple(values[index] for index in permutation)
    permuted["balance_matrix"] = tuple(
        tuple(row[index] for index in permutation)
        for row in base["balance_matrix"]  # type: ignore[union-attr]
    )
    permuted["reaction_matrix"] = tuple(
        tuple(row[index] for index in permutation)
        for row in base["reaction_matrix"]  # type: ignore[union-attr]
    )
    _bind_record(permuted)

    base_result = _equilibrium._chemical_equilibrium(None, base, None, None, 1.0e-12)[
        "sensitivities"
    ]
    permuted_result = _equilibrium._chemical_equilibrium(None, permuted, None, None, 1.0e-12)[
        "sensitivities"
    ]

    assert permuted_result["parameter_order"] == base_result["parameter_order"]
    assert permuted_result["condition_number_inf"] == pytest.approx(
        base_result["condition_number_inf"], rel=2.0e-10
    )
    for parameter in range(3):
        assert permuted_result["amount_derivatives"][parameter] == pytest.approx(
            tuple(base_result["amount_derivatives"][parameter][index] for index in permutation),
            abs=2.0e-9,
        )


def test_manufactured_reaction_with_ill_conditioned_kkt_has_no_sensitivity() -> None:
    balance_separation = 3.0e-6
    spec = _base_system()
    spec["species_ids"] = ("A", "B", "C")
    spec["charges"] = (0, 0, 0)
    spec["molar_masses_kg_per_mol"] = (1.0, 1.0, 1.0)
    spec["balance_matrix"] = ((1.0, 1.0, 1.0 + balance_separation),)
    spec["reaction_matrix"] = ((-1.0, 1.0, 0.0),)
    spec["feed_amounts"] = (0.5, 0.5, 1.0)
    spec["ln_k"] = (0.0,)
    spec["temperature_k"] = 300.0
    spec["pressure_pa"] = 8.31446261815324 * 300.0
    _bind_record(spec)

    native = _equilibrium._chemical_equilibrium(None, spec, None, None, 1.0e-12)

    sensitivities = native["sensitivities"]
    assert sensitivities["kkt_rank"] == sensitivities["kkt_dimension"] == 4
    assert sensitivities["condition_number_inf"] > 1.0e6
    assert native["accepted"] is False
    criteria = {record["name"]: record for record in native["numerical_criteria"]}
    assert (
        criteria["kkt_condition_number_inf"]["value"]
        > criteria["kkt_condition_number_inf"]["limit"]
    )

    parameters = epcsaft.Parameters.from_mapping(NEUTRAL_TERNARY_IMMISCIBLE)
    spec["species_ids"] = parameters.component_ids
    phase = equilibrium.ReactivePhase(
        "vapor",
        "vapor",
        "finite",
        equilibrium.AllComponents(),
        equilibrium.IdealGasModel("provider-helmholtz-coordinate-basis"),
    )
    result = equilibrium.solve(
        epcsaft.Mixture(parameters),
        equilibrium.GeneralReactiveEquilibriumProblem(
            "ill-conditioned-ideal-reaction",
            equilibrium.Fixed(spec["temperature_k"] * epcsaft.unit_registry.kelvin),
            equilibrium.Fixed(spec["pressure_pa"] * epcsaft.unit_registry.pascal),
            (phase,),
            _typed_problem(spec),
            (phase.identity,),
            (
                equilibrium.EquilibriumOutput(
                    "amount",
                    "phase.amount",
                    "mole",
                    "true-species-amount",
                    phase.identity,
                    (1.0, 0.0, 0.0),
                ),
            ),
        ),
    )
    assert result.status == "non_evaluable"
    assert result.failure is not None
    assert result.failure.code == "singular_sensitivity"


@pytest.mark.parametrize("trace_floor", (1.0e-50, 1.0e-55))
def test_belov_aristova_gas_restriction_resolves_extreme_positive_traces(
    trace_floor: float,
) -> None:
    spec, expected = _belov_aristova_gas_system()

    result = _manufactured_solve(spec, {"trace_floor": trace_floor})

    amounts = result.amounts_mol
    assert tuple(map(math.log, amounts)) == pytest.approx(
        tuple(map(math.log, expected)), abs=2.0e-6
    )
    balances = spec["balance_matrix"]
    assert tuple(
        sum(row[index] * amounts[index] for index in range(8))
        for row in balances  # type: ignore[union-attr]
    ) == pytest.approx((2.0, 1.0), abs=1.0e-9)
    assert (8.31446261815324 * 2000.0 * sum(amounts) / result.volume_m3) == pytest.approx(
        100_000.0, rel=1.0e-8
    )
    concentration_shift = math.log(8.31446261815324 * 2000.0 / 101_325.0)
    potentials = tuple(
        _BELOV_SOURCE_GIBBS[index]
        + concentration_shift
        + math.log(amounts[index] / result.volume_m3)
        for index in range(8)
    )
    reactions = spec["reaction_matrix"]
    assert tuple(
        sum(row[index] * potentials[index] for index in range(8))
        for row in reactions  # type: ignore[union-attr]
    ) == pytest.approx((0.0,) * 6, abs=1.0e-7)
    assert tuple(
        potentials[index]
        + balances[0][index] * _BELOV_ELEMENT_POTENTIALS[0]  # type: ignore[index]
        + balances[1][index] * _BELOV_ELEMENT_POTENTIALS[1]  # type: ignore[index]
        for index in range(8)
    ) == pytest.approx((0.0,) * 8, abs=1.0e-7)
    diagnostics = result.diagnostics
    assert diagnostics.balance_inf_norm <= 1.0e-9
    assert diagnostics.pressure_relative_residual <= 1.0e-8
    assert diagnostics.reaction_affinity_inf_norm <= 1.0e-7
    assert diagnostics.kkt_stationarity_inf_norm <= 1.0e-7


@pytest.mark.parametrize(
    "variant", ("reaction_basis", "species_order", "conservation_gauge", "feed_scale")
)
def test_belov_aristova_trace_solution_is_coordinate_invariant(variant: str) -> None:
    spec, expected = _belov_aristova_gas_system()
    options: dict[str, object] = {"trace_floor": 1.0e-50}

    if variant == "reaction_basis":
        reactions = spec["reaction_matrix"]
        ln_k = spec["ln_k"]
        spec["reaction_matrix"] = (
            tuple(
                reactions[0][index] + reactions[1][index]  # type: ignore[index]
                for index in range(8)
            ),
            *reactions[1:],  # type: ignore[index]
        )
        spec["ln_k"] = (ln_k[0] + ln_k[1], *ln_k[1:])  # type: ignore[index]
        _bind_record(spec)
    elif variant == "species_order":
        permutation = (7, 5, 3, 1, 6, 4, 2, 0)
        for field in (
            "species_ids",
            "charges",
            "molar_masses_kg_per_mol",
            "feed_amounts",
        ):
            values = spec[field]
            spec[field] = tuple(values[index] for index in permutation)  # type: ignore[index]
        spec["balance_matrix"] = tuple(
            tuple(row[index] for index in permutation)
            for row in spec["balance_matrix"]  # type: ignore[union-attr]
        )
        spec["reaction_matrix"] = tuple(
            tuple(row[index] for index in permutation)
            for row in spec["reaction_matrix"]  # type: ignore[union-attr]
        )
        expected = tuple(expected[index] for index in permutation)
    elif variant == "conservation_gauge":
        balances = spec["balance_matrix"]
        spec["balance_matrix"] = (
            tuple(
                balances[0][index] + balances[1][index]  # type: ignore[index]
                for index in range(8)
            ),
            tuple(
                2.0 * balances[0][index] - balances[1][index]  # type: ignore[index]
                for index in range(8)
            ),
        )
        _bind_record(spec)
    else:
        scale = 3.0
        spec["feed_amounts"] = tuple(
            scale * value
            for value in spec["feed_amounts"]  # type: ignore[union-attr]
        )
        expected = tuple(scale * value for value in expected)
        _bind_record(spec)

    result = _manufactured_solve(spec, options)

    assert tuple(map(math.log, result.amounts_mol)) == pytest.approx(
        tuple(map(math.log, expected)), abs=3.0e-6
    )
    assert result.diagnostics.reaction_affinity_inf_norm <= 1.0e-7


def test_manufactured_charged_reaction_uses_exact_electroneutral_chart() -> None:
    temperature_k = 300.0
    spec = {
        **_base_system(),
        "species_ids": ("A", "C+", "D-"),
        "charges": (0, 1, -1),
        "molar_masses_kg_per_mol": (2.0, 1.0, 1.0),
        "balance_matrix": ((2.0, 1.0, 1.0), (0.0, 1.0, -1.0)),
        "reaction_matrix": ((-1.0, 1.0, 1.0),),
        "feed_amounts": (1.0, 0.0, 0.0),
        "ln_k": (math.log(1.0 / 3.0),),
        "temperature_k": temperature_k,
        "pressure_pa": 8.31446261815324 * temperature_k,
    }
    _bind_record(spec)

    result = _manufactured_solve(spec)

    assert result.amounts_mol == pytest.approx((0.5, 0.5, 0.5), rel=3.0e-8)
    assert result.volume_m3 == pytest.approx(1.5, rel=3.0e-8)
    assert result.diagnostics.charge_inf_norm <= 2.0e-15


def test_manufactured_ultra_trace_charged_share_stays_differentiable() -> None:
    """Both charged softmax references and explicit shares stay interior."""

    temperature_k = 300.0
    trace_floor = 1.0e-18
    records = tuple(
        {
            "source_id": f"manufactured:reaction-{index}",
            "reference_id": "provider-helmholtz-coordinate-basis",
            "reaction_orientation": "products_positive",
            "conversion_id": "already-provider-basis",
            "dimensionless": True,
            "temperature_k": temperature_k,
            "pressure_pa": 8.31446261815324 * temperature_k,
        }
        for index in range(2)
    )

    def solve(log_ratio_constants: tuple[float, float], floor: float) -> dict[str, object]:
        spec = {
            **_base_system(),
            "species_ids": ("C1+", "C2+", "D1-", "D2-"),
            "charges": (1, 1, -1, -1),
            "molar_masses_kg_per_mol": (1.0, 1.0, 1.0, 1.0),
            "balance_matrix": ((1.0, 1.0, 1.0, 1.0), (1.0, 1.0, -1.0, -1.0)),
            "reaction_matrix": ((1.0, -1.0, 0.0, 0.0), (0.0, 0.0, 1.0, -1.0)),
            "feed_amounts": (1.0, 0.0, 1.0, 0.0),
            "ln_k": log_ratio_constants,
            "temperature_k": temperature_k,
            "pressure_pa": 8.31446261815324 * temperature_k,
            "equilibrium_constant_records": records,
        }
        spec["conserved_totals"] = (2.0, 0.0)
        return _equilibrium._chemical_equilibrium(None, spec, None, None, floor)

    explicit_trace = solve((-40.0, 0.0), trace_floor)
    assert explicit_trace["amounts"][0] > trace_floor

    reference_trace = solve((40.0, 40.0), trace_floor)
    assert reference_trace["amounts"][1] > trace_floor
    assert reference_trace["amounts"][3] > trace_floor

    rejected_trace = solve((40.0, 40.0), 1.0e-12)
    assert min(rejected_trace["amounts"]) <= 1.0e-12


def test_manufactured_equilibrium_is_gauge_scale_and_reaction_basis_invariant() -> None:
    base = _base_system()
    base["feed_amounts"] = (1.0, 0.0)
    _bind_record(base)
    plain = _manufactured_solve(base)
    gauge = copy.deepcopy(base)
    gauge["balance_matrix"] = ((3.25, 3.25),)
    _bind_record(gauge)
    gauged = _manufactured_solve(gauge)
    assert gauged.amounts_mol == pytest.approx(plain.amounts_mol, rel=2.0e-10)

    scaled = copy.deepcopy(base)
    scaled["feed_amounts"] = (7.0, 0.0)
    _bind_record(scaled)
    scaled_result = _manufactured_solve(scaled)
    assert scaled_result.amounts_mol == pytest.approx(
        tuple(7.0 * value for value in plain.amounts_mol), rel=3.0e-8
    )
    assert scaled_result.volume_m3 == pytest.approx(7.0 * plain.volume_m3, rel=3.0e-8)

    three = {
        **_base_system(),
        "species_ids": ("A", "B", "C"),
        "charges": (0, 0, 0),
        "molar_masses_kg_per_mol": (1.0, 1.0, 1.0),
        "balance_matrix": ((1.0, 1.0, 1.0),),
        "reaction_matrix": ((-1.0, 1.0, 0.0), (0.0, -1.0, 1.0)),
        "feed_amounts": (1.0, 0.0, 0.0),
        "ln_k": (math.log(2.0), math.log(3.0)),
    }
    _bind_record(three)
    first_basis = _manufactured_solve(three)
    changed_basis = copy.deepcopy(three)
    changed_basis["reaction_matrix"] = ((-1.0, 0.0, 1.0), (0.0, -1.0, 1.0))
    changed_basis["ln_k"] = (math.log(6.0), math.log(3.0))
    _bind_record(changed_basis)
    second_basis = _manufactured_solve(changed_basis)
    assert first_basis.amounts_mol == pytest.approx(second_basis.amounts_mol, rel=3.0e-9)


def test_manufactured_reaction_solve_has_no_feed_scaled_amount_cap() -> None:
    target = (0.5, 500.0, 0.5)
    volume = sum(target)
    spec = {
        **_base_system(),
        "species_ids": ("A", "B", "C"),
        "charges": (0, 0, 0),
        "molar_masses_kg_per_mol": (1001.0, 1.0, 1.0),
        "balance_matrix": ((1000.0, 1.0, 0.0), (1.0, 0.0, 1.0)),
        "reaction_matrix": ((-1.0, 1000.0, 1.0),),
        "feed_amounts": (1.0, 0.0, 0.0),
        "ln_k": (1000.0 * math.log(target[1] / volume),),
        "temperature_k": 300.0,
        "pressure_pa": 8.31446261815324 * 300.0,
    }
    _bind_record(spec)

    result = _manufactured_solve(spec, {"trace_floor": 0.1})

    assert result.amounts_mol == pytest.approx(target, rel=2.0e-7, abs=2.0e-8)


def test_manufactured_solver_rejects_trace_false_success() -> None:
    spec = _base_system()
    trace = copy.deepcopy(spec)
    trace["ln_k"] = (math.log(1.0e-10),)
    _bind_record(trace)
    with pytest.raises(equilibrium.ChemicalEquilibriumError) as failed:
        _manufactured_solve(trace, {"trace_floor": 1.0e-8})
    diagnostics = failed.value.diagnostics
    assert diagnostics.search.selected_basin_ordinal is None
    native = _equilibrium._chemical_equilibrium(None, trace, None, None, 1.0e-8)
    assert min(native["amounts"]) <= 1.0e-8


def test_manufactured_charged_solution_is_species_order_invariant() -> None:
    temperature_k = 300.0
    original = {
        **_base_system(),
        "species_ids": ("A", "C+", "D-"),
        "charges": (0, 1, -1),
        "molar_masses_kg_per_mol": (2.0, 1.0, 1.0),
        "balance_matrix": ((2.0, 1.0, 1.0), (0.0, 1.0, -1.0)),
        "reaction_matrix": ((-1.0, 1.0, 1.0),),
        "feed_amounts": (1.0, 0.0, 0.0),
        "ln_k": (math.log(1.0 / 3.0),),
        "temperature_k": temperature_k,
        "pressure_pa": 8.31446261815324 * temperature_k,
    }
    _bind_record(original)
    baseline = _manufactured_solve(original)
    permutation = (2, 0, 1)
    permuted = copy.deepcopy(original)
    for field in (
        "species_ids",
        "charges",
        "molar_masses_kg_per_mol",
        "feed_amounts",
    ):
        values = original[field]
        permuted[field] = tuple(values[index] for index in permutation)
    permuted["balance_matrix"] = tuple(
        tuple(row[index] for index in permutation) for row in original["balance_matrix"]
    )
    permuted["reaction_matrix"] = tuple(
        tuple(row[index] for index in permutation) for row in original["reaction_matrix"]
    )
    _bind_record(permuted)
    reordered = _manufactured_solve(permuted)
    inverse = tuple(permutation.index(index) for index in range(len(permutation)))
    assert tuple(reordered.amounts_mol[index] for index in inverse) == pytest.approx(
        baseline.amounts_mol, rel=4.0e-8
    )
