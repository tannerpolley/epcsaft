from __future__ import annotations

from collections.abc import Sequence
from typing import cast

import numpy as np
import pytest
from scipy.optimize import least_squares

from epcsaft import ActiveParameterSet, Mixture, Parameters, _core
from epcsaft import unit_registry as u
from epcsaft.equilibrium import (
    AllComponents,
    DeclaredTopology,
    FixedTPFeed,
    PhaseEquilibriumProblem,
    PhaseSpecification,
    solve,
)
from tests.eos._raw_cases import AQUEOUS_ELECTROLYTE, NEUTRAL_BINARY

_TEMPERATURE_K = 243.58
_PRESSURE_PA = 3.949e6
_FEED = (0.48815, 0.51185)
_PARAMETER_IDENTITIES = (
    "pair/ethane/methane/k_ij",
    "component/methane/dispersion_energy_over_k",
)
_PARAMETER_STEPS = (1.0e-4, 1.0e-2)
_INITIAL_STATE = (
    0.14985423556818814,
    0.34549510864416604,
    3.276689695199133e-5,
    1.808100839084689e-4,
)


def _provider_block(
    native_model: object,
    amounts: tuple[float, float],
    volume: float,
) -> dict[str, object]:
    return _core.evaluate_parameterized_phase(
        native_model,
        _TEMPERATURE_K,
        amounts,
        volume,
        1.0e-6,
        0.74,
    )


def _solve_two_phase_state(
    parameters: Parameters,
    identities: tuple[str, ...],
    values: tuple[float, ...],
    initial: tuple[float, ...] = _INITIAL_STATE,
) -> tuple[tuple[tuple[float, float], tuple[float, float]], tuple[float, float]]:
    active = ActiveParameterSet(parameters, identities, values=values)
    native_model = Mixture(parameters)._native_model(active)

    def residual(state: np.ndarray[tuple[int], np.dtype[np.float64]]) -> np.ndarray:
        first_amounts = (float(state[0]), float(state[1]))
        second_amounts = (
            _FEED[0] - first_amounts[0],
            _FEED[1] - first_amounts[1],
        )
        first = _provider_block(native_model, first_amounts, float(state[2]))
        second = _provider_block(native_model, second_amounts, float(state[3]))
        first_gradient = cast(Sequence[float], first["state_gradient"])
        second_gradient = cast(Sequence[float], second["state_gradient"])
        return np.asarray(
            (
                float(first_gradient[0]) - float(second_gradient[0]),
                float(first_gradient[1]) - float(second_gradient[1]),
                float(first["pressure_pa"]) / _PRESSURE_PA - 1.0,
                float(second["pressure_pa"]) / _PRESSURE_PA - 1.0,
            )
        )

    solved = least_squares(
        residual,
        np.asarray(initial),
        bounds=(
            (1.0e-8, 1.0e-8, 1.0e-8, 1.0e-8),
            (_FEED[0] - 1.0e-8, _FEED[1] - 1.0e-8, 1.0e-2, 1.0e-2),
        ),
        jac="2-point",
        diff_step=1.0e-5,
        x_scale="jac",
        ftol=1.0e-14,
        xtol=1.0e-14,
        gtol=1.0e-14,
        max_nfev=100,
    )
    assert solved.success
    assert max(abs(value) for value in residual(solved.x)) <= 2.0e-11
    first_amounts = (float(solved.x[0]), float(solved.x[1]))
    return (
        (
            first_amounts,
            (_FEED[0] - first_amounts[0], _FEED[1] - first_amounts[1]),
        ),
        (float(solved.x[2]), float(solved.x[3])),
    )


def _sensitivity(
    parameters: Parameters,
    identities: tuple[str, ...],
    values: tuple[float, ...],
    amounts: tuple[tuple[float, float], tuple[float, float]],
    volumes: tuple[float, float],
    feed: tuple[float, float] = _FEED,
) -> dict[str, object]:
    active = ActiveParameterSet(parameters, identities, values=values)
    return _core.evaluate_solved_state_sensitivity(
        Mixture(parameters)._native_model(active),
        _TEMPERATURE_K,
        _PRESSURE_PA,
        feed,
        amounts,
        volumes,
    )


def _matrix(
    sensitivity: dict[str, object],
    name: str,
    rows: int,
    columns: int,
) -> tuple[tuple[float, ...], ...]:
    values = tuple(float(value) for value in cast(Sequence[object], sensitivity[name]))
    assert len(values) == rows * columns
    return tuple(values[row * columns : (row + 1) * columns] for row in range(rows))


def _outputs(
    amounts: tuple[tuple[float, float], tuple[float, float]],
    volumes: tuple[float, float],
) -> tuple[float, ...]:
    return tuple(value for phase in range(2) for value in (*amounts[phase], volumes[phase]))


def test_exact_two_phase_sensitivity_matches_independent_centered_resolves() -> None:
    """Falsify exact neutral two-phase totals, ordering, and refusal boundaries.

    Oracle: independently converged nonlinear defining-equation solves and
    centered parameter resolves. Compared quantities are dimensionless mu/RT
    and P/P_target residuals, amount derivatives in mol per parameter unit,
    and volume derivatives in m^3 per parameter unit. Solver residuals use
    2e-11 absolute tolerance; total derivatives use 2e-6 relative and 2e-11
    absolute tolerance; declared permutations use 2e-12 relative/absolute;
    the linear backward error is bounded by 1e-13. These tolerances cover
    double-precision roundoff, nonlinear solve error, and O(h^2) centered
    truncation at the explicit steps. Failure means the implicit equations,
    parameter identity/order, output permutation, or nondifferentiable-state
    refusal is scientifically incorrect.
    """
    parameters = Parameters.from_mapping(NEUTRAL_BINARY)
    selected = ActiveParameterSet(parameters, _PARAMETER_IDENTITIES)
    values = selected.values
    amounts, volumes = _solve_two_phase_state(parameters, _PARAMETER_IDENTITIES, values)
    sensitivity = _sensitivity(
        parameters,
        _PARAMETER_IDENTITIES,
        values,
        amounts,
        volumes,
    )

    assert sensitivity["status"] == "available"
    assert sensitivity["failure_reason"] == ""
    assert sensitivity["derivative_order"] == 1
    assert sensitivity["higher_order_status"] == "unavailable"
    assert sensitivity["higher_order_reason"] == ("higher_order_total_derivatives_not_implemented")
    assert sensitivity["tensor_layout"] == "row-major"
    assert tuple(cast(Sequence[object], sensitivity["parameter_identities"])) == (
        _PARAMETER_IDENTITIES
    )
    assert tuple(cast(Sequence[object], sensitivity["residual_identities"])) == (
        "phase[0].mu[methane]-phase[1].mu[methane]",
        "phase[0].mu[ethane]-phase[1].mu[ethane]",
        "phase[0].pressure_over_target_minus_one",
        "phase[1].pressure_over_target_minus_one",
    )
    assert int(sensitivity["rank"]) == int(sensitivity["dimension"]) == 4
    assert float(sensitivity["scaled_condition_number_2"]) < 1.0e8
    assert float(sensitivity["linear_solve_backward_error"]) <= 1.0e-13
    exact = _matrix(sensitivity, "output_parameter_jacobian", 6, 2)

    for parameter, step in enumerate(_PARAMETER_STEPS):
        lower = list(values)
        upper = list(values)
        lower[parameter] -= step
        upper[parameter] += step
        lower_state = _solve_two_phase_state(
            parameters,
            _PARAMETER_IDENTITIES,
            tuple(lower),
            (*amounts[0], *volumes),
        )
        upper_state = _solve_two_phase_state(
            parameters,
            _PARAMETER_IDENTITIES,
            tuple(upper),
            (*amounts[0], *volumes),
        )
        centered = tuple(
            (upper_value - lower_value) / (2.0 * step)
            for lower_value, upper_value in zip(
                _outputs(*lower_state), _outputs(*upper_state), strict=True
            )
        )
        assert tuple(row[parameter] for row in exact) == pytest.approx(
            centered, rel=2.0e-6, abs=2.0e-11
        )

    phase_permuted = _sensitivity(
        parameters,
        _PARAMETER_IDENTITIES,
        values,
        tuple(reversed(amounts)),
        tuple(reversed(volumes)),
    )
    phase_permuted_jacobian = _matrix(phase_permuted, "output_parameter_jacobian", 6, 2)
    expected_phase_permutation = exact[3:] + exact[:3]
    for observed, expected in zip(phase_permuted_jacobian, expected_phase_permutation, strict=True):
        assert observed == pytest.approx(expected, rel=2.0e-12, abs=2.0e-12)

    reversed_parameters = parameters.select(("ethane", "methane"))
    reversed_identities = tuple(reversed(_PARAMETER_IDENTITIES))
    reordered = _sensitivity(
        reversed_parameters,
        reversed_identities,
        tuple(reversed(values)),
        tuple((phase[1], phase[0]) for phase in amounts),
        volumes,
        tuple(reversed(_FEED)),
    )
    reordered_jacobian = _matrix(reordered, "output_parameter_jacobian", 6, 2)
    expected_reordering = tuple(
        tuple(reversed(row)) for row in (exact[1], exact[0], exact[2], exact[4], exact[3], exact[5])
    )
    for observed, expected in zip(reordered_jacobian, expected_reordering, strict=True):
        assert observed == pytest.approx(expected, rel=2.0e-12, abs=2.0e-12)

    unavailable = _sensitivity(
        parameters,
        _PARAMETER_IDENTITIES,
        values,
        (
            (_FEED[0] / 2.0, _FEED[1] / 2.0),
            (_FEED[0] / 2.0, _FEED[1] / 2.0),
        ),
        (1.0e-4, 1.0e-4),
    )
    assert unavailable["status"] == "unavailable"
    assert unavailable["failure_reason"] == "nonunique_coalesced_phases"

    unsolved = _sensitivity(
        parameters,
        _PARAMETER_IDENTITIES,
        values,
        amounts,
        (volumes[0] * 1.01, volumes[1]),
    )
    assert unsolved["status"] == "unavailable"
    assert unsolved["failure_reason"] == "defining_equations_not_satisfied"

    boundary = _sensitivity(
        parameters,
        _PARAMETER_IDENTITIES,
        values,
        (
            (0.0, amounts[0][1]),
            (_FEED[0], amounts[1][1]),
        ),
        volumes,
    )
    assert boundary["status"] == "unavailable"
    assert boundary["failure_reason"] == "phase_boundary_not_differentiable"

    public = solve(
        Mixture(parameters),
        PhaseEquilibriumProblem(
            FixedTPFeed(
                300.0 * u.kelvin,
                1.0e5 * u.pascal,
                (0.5, 0.5),
                DeclaredTopology(
                    (PhaseSpecification("vapor", AllComponents()),)
                ),
            )
        ),
        active_parameters=selected,
    )
    assert public.sensitivity is not None
    assert public.sensitivity.status == "available"
    assert public.sensitivity.parameters == selected.coordinates
    assert public.sensitivity.derivative_order == 1
    assert public.sensitivity.higher_order_status == "unavailable"
    assert public.sensitivity.tensor_layout == "row-major"


def test_one_phase_association_and_model_sensitivity_match_centered_resolves() -> None:
    """Falsify ionic one-phase association/model totals and topology refusal.

    Oracle: independent centered pressure re-solves in log-volume. The
    compared volume derivatives are m^3/K for association energy_over_k and
    m^3 per dimensionless suppression coefficient; pressure is scaled by its
    target. Root residuals use 1e-12 absolute tolerance and total derivatives
    use 2e-6 relative plus 1e-12 m^3-per-parameter absolute tolerance. The
    bounds cover double-precision solve error and O(h^2) truncation at steps
    1e-3 K and 1e-4, respectively. Failure means association/model EOS
    partials were not carried through the solved state, or charged multiphase
    constraints were overclaimed before the unified constrained owner exists.
    """
    temperature_k = 298.15
    feed = (0.98, 0.01, 0.01)
    volume_m3 = 5.0e-4
    identities = (
        "association/water/a/water/b/energy_over_k",
        "model/relative_permittivity/ion_fraction_suppression_coefficient",
    )
    steps = (1.0e-3, 1.0e-4)
    parameters = Parameters.from_mapping(AQUEOUS_ELECTROLYTE)
    selected = ActiveParameterSet(parameters, identities)

    def native_model(values: tuple[float, ...]) -> object:
        active = ActiveParameterSet(parameters, identities, values=values)
        return Mixture(parameters)._native_model(active)

    def phase(provider: object, volume: float) -> dict[str, object]:
        return _core.evaluate_parameterized_phase(
            provider,
            temperature_k,
            feed,
            volume,
            1.0e-6,
            0.74,
        )

    base = native_model(selected.values)
    pressure_pa = float(phase(base, volume_m3)["pressure_pa"])
    sensitivity = _core.evaluate_solved_state_sensitivity(
        base,
        temperature_k,
        pressure_pa,
        feed,
        (feed,),
        (volume_m3,),
    )
    assert sensitivity["status"] == "available"
    exact = _matrix(sensitivity, "output_parameter_jacobian", 4, 2)

    for parameter, step in enumerate(steps):
        volumes = []
        for sign in (-1.0, 1.0):
            values = list(selected.values)
            values[parameter] += sign * step
            provider = native_model(tuple(values))
            solved = least_squares(
                lambda log_volume, provider=provider: np.asarray(
                    (
                        float(phase(provider, float(np.exp(log_volume[0])))["pressure_pa"])
                        / pressure_pa
                        - 1.0,
                    )
                ),
                np.asarray((np.log(volume_m3),)),
                ftol=1.0e-14,
                xtol=1.0e-14,
                gtol=1.0e-14,
                max_nfev=50,
            )
            assert solved.success
            assert abs(float(solved.fun[0])) <= 1.0e-12
            volumes.append(float(np.exp(solved.x[0])))
        centered = (volumes[1] - volumes[0]) / (2.0 * step)
        assert exact[3][parameter] == pytest.approx(centered, rel=2.0e-6, abs=1.0e-12)

    charged_multiphase = _core.evaluate_solved_state_sensitivity(
        base,
        temperature_k,
        pressure_pa,
        feed,
        (
            tuple(amount / 2.0 for amount in feed),
            tuple(amount / 2.0 for amount in feed),
        ),
        (volume_m3 / 4.0, 3.0 * volume_m3 / 4.0),
    )
    assert charged_multiphase["status"] == "unavailable"
    assert charged_multiphase["failure_reason"] == ("unsupported_charged_multiphase_constraints")
