from __future__ import annotations

import json
import math
from copy import deepcopy
from pathlib import Path

import pytest

import epcsaft
from epcsaft import _core
from epcsaft import equilibrium
from epcsaft.regression import fit
from tests.regression._cases import (
    fixed_state_fit_case,
    fixed_state_observation,
    full_space_fit_case,
    pressure,
    scaled_error,
    solved_pressure_fit_case,
    with_parameter,
)


class _NativePayloadCaptured(Exception):
    pass


def _native_full_space_payload(
    mapping: dict[str, object], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> tuple[object, ...]:
    """Capture the existing resolved native request without a test-only API."""

    payload: dict[str, tuple[object, ...]] = {}
    native_fit = _core._fit_resolved

    def capture(*arguments: object) -> object:
        payload["arguments"] = arguments
        raise _NativePayloadCaptured

    monkeypatch.setattr(_core, "_fit_resolved", capture)
    with pytest.raises(_NativePayloadCaptured):
        fit(mapping, base_path=tmp_path)
    monkeypatch.setattr(_core, "_fit_resolved", native_fit)
    return payload["arguments"]


def test_maximum_iterations_is_nonusable_and_retains_final_diagnostics(
    tmp_path: Path,
) -> None:
    """Falsify max-iteration acceptance or loss of an inspectable final state."""

    reduced_mapping, *_ = fixed_state_fit_case(tmp_path)
    reduced_mapping["solver"]["controls"]["maximum_iterations"] = 1
    reduced = fit(reduced_mapping, base_path=tmp_path)
    assert all(start.solver_status == "maximum_iterations_exceeded" for start in reduced.starts)
    assert all(not start.usable for start in reduced.starts)
    assert reduced.best_usable_start is None
    for start in reduced.starts:
        kkt = start.kkt_diagnostics
        assert kkt["retained"]
        assert len(kkt["coordinate_lower_bound_multipliers"]) == len(
            start.final_affine
        )
        assert len(kkt["coordinate_upper_bound_multipliers"]) == len(
            start.final_affine
        )
        assert kkt["parameter_constraint_multipliers"] == []
        assert all(
            math.isfinite(kkt[name])
            for name in (
                "objective",
                "primal_feasibility_inf_norm",
                "stationarity_inf_norm",
                "complementarity_inf_norm",
                "dual_sign_violation_inf_norm",
            )
        )

    full_mapping, *_ = full_space_fit_case(tmp_path)
    full_mapping["solver"]["controls"]["maximum_iterations"] = 1
    stalled = fit(full_mapping, base_path=tmp_path).starts[0]
    assert stalled.solver_status == "maximum_iterations_exceeded"
    assert not stalled.usable
    assert all(math.isfinite(value) for value in stalled.final_physical)
    diagnostic = stalled.owner_diagnostics[0]
    assert diagnostic["solver_status"]
    assert diagnostic["local_variables"]
    certificate = diagnostic["certificate"]
    assert set(certificate) == {
        "balance_inf_norm",
        "pressure_relative_inf_norm",
        "reaction_affinity_inf_norm",
        "transfer_chemical_potential_inf_norm",
        "kkt_stationarity_inf_norm",
        "phase_distinctness",
        "minimum_amount_mol",
        "minimum_bound_margin",
        "reduced_hessian_min_eigenvalue",
        "reduced_hessian_scale",
        "reduced_hessian_threshold",
    }
    assert all(value is not None for value in certificate.values())
    assert stalled.owner_failures[0].reason == diagnostic["failure_reason"]
    kkt = stalled.kkt_diagnostics
    assert kkt["retained"]
    assert len(kkt["coordinate_lower_bound_multipliers"]) == len(
        stalled.final_affine
    )
    assert len(kkt["coordinate_upper_bound_multipliers"]) == len(
        stalled.final_affine
    )
    assert kkt["parameter_constraint_multipliers"] == []
    assert all(
        math.isfinite(kkt[name])
        for name in (
            "objective",
            "primal_feasibility_inf_norm",
            "stationarity_inf_norm",
            "complementarity_inf_norm",
            "dual_sign_violation_inf_norm",
        )
    )


def test_finite_owner_full_space_recovers_truth_with_complete_certificate(
    tmp_path: Path,
) -> None:
    """Falsify finite-owner recovery, derivatives, KKT evidence, or mixed-owner order."""

    mapping, parameters, model, problem, identity, truth_value, truth = full_space_fit_case(
        tmp_path
    )
    result = fit(mapping, base_path=tmp_path)
    accepted = result.best_usable_start
    assert accepted is not None

    coordinate_scale = 0.2
    recovery_limit = 1.0e-6
    observation_relative_tolerance = 2.0e-6
    assert (
        max(abs(value) for value in accepted.affine_jacobian) * recovery_limit
        < observation_relative_tolerance
    )
    assert abs(accepted.final_physical[0] - truth_value) / coordinate_scale <= recovery_limit
    assert accepted.predictions == pytest.approx(
        truth.values, rel=observation_relative_tolerance, abs=2.0e-8
    )
    assert accepted.owner_failures == ()
    assert accepted.batch_outcomes[0].status == "evaluated"
    assert result.accounting.dropped_rows == 0
    assert result.accounting.rejected_row_trials == 0
    assert result.accounting.failed_row_trials == 0

    step = 1.0e-5
    centered_values = []
    for value in (accepted.final_physical[0] - step, accepted.final_physical[0] + step):
        solved = equilibrium.solve(
            model,
            problem,
            active_parameters=epcsaft.ActiveParameterSet(parameters, (identity,), values=(value,)),
        )
        assert solved.status == "evaluated" and solved.values[0] is not None
        centered_values.append(float(solved.values[0]))
    centered = (centered_values[1] - centered_values[0]) / (2.0 * step)
    row = result.specification.rows[0]
    owner_derivative = accepted.physical_jacobian[0] * row.scale / row.multiplier
    assert scaled_error(owner_derivative, centered) <= 5.0e-5

    assert len(accepted.owner_receipts) == 1
    receipt = json.loads(accepted.owner_receipts[0].canonical_bytes_utf8)
    evidence = dict(receipt["evidence"])
    assert evidence["balance_inf_norm"] <= 1.0e-9
    assert evidence["pressure_relative_inf_norm"] <= 1.0e-8
    assert evidence["reaction_affinity_inf_norm"] <= 1.0e-6
    assert evidence["kkt_stationarity_inf_norm"] <= 1.0e-7
    expected_branch = json.loads(json.dumps(truth.branch_signature))
    assert receipt["branch_signature"] == expected_branch
    assert receipt["certificate_scope"] == truth.certificate_scope
    assert receipt["globality_claim"] == truth.globality_claim
    assert {name: outcome["status"] for name, outcome in receipt["outcomes"].items()} == {
        "solver_convergence": "passed",
        "numerical_convergence": "passed",
        "physical_certification": "passed",
        "scientific_fit_acceptance": "not_applicable",
    }
    assert result.evaluated_rows[0].identity == "manufactured-equilibrium-row"
    assert result.evaluated_rows[0].unit == problem.outputs[0].unit
    assert result.evaluated_rows[0].basis == problem.outputs[0].basis

    assert isinstance(problem.temperature, equilibrium.Fixed)
    phase = truth.phases[0]
    temperature = float(problem.temperature.value.to("kelvin").magnitude)
    direct_pressure = pressure(
        with_parameter(parameters, identity, truth_value, "mixed-owner-truth"),
        temperature,
        phase.molar_density_mol_m3,
        phase.mole_fractions,
    )
    direct = fixed_state_observation(
        (
            {
                "temperature_k": temperature,
                "molar_density_mol_per_m3": phase.molar_density_mol_m3,
                "pressure": direct_pressure,
            },
        ),
        phase.mole_fractions,
        identity="mixed-direct-owner",
    )
    observations = mapping["observations"]
    assert isinstance(observations, list)
    mapping["observations"] = [direct, *observations]

    mixed_result = fit(mapping, base_path=tmp_path)
    mixed = mixed_result.best_usable_start
    assert mixed is not None
    assert mixed.predictions == pytest.approx(
        (direct_pressure, float(truth.values[0])), rel=2.0e-6, abs=2.0e-8
    )
    assert tuple(row.selector for row in mixed_result.specification.rows) == (
        "pressure",
        "cation-fraction",
    )
    assert tuple(outcome.status for outcome in mixed.batch_outcomes) == (
        "evaluated",
        "evaluated",
    )


@pytest.mark.slow
@pytest.mark.parametrize("ideal_vapor", [True, False], ids=("ideal", "provider"))
def test_solved_pressure_full_space_derivatives_match_independent_resolves(
    tmp_path: Path,
    ideal_vapor: bool,
) -> None:
    """Falsify solved-pressure liquid response, vapor outputs, or residual scaling."""

    mapping, parameters, model, problem, identity, truth_value, _ = solved_pressure_fit_case(
        tmp_path, ideal_vapor=ideal_vapor
    )
    result = fit(mapping, base_path=tmp_path)
    accepted = result.best_usable_start
    assert accepted is not None
    assert result.evaluation_strategy == "simultaneous_sparse_full_space"
    assert accepted.final_physical == pytest.approx((truth_value,), abs=3.0e-6)
    assert accepted.owner_failures == ()
    assert result.accounting.failed_row_trials == 0
    assert result.accounting.rejected_row_trials == 0
    assert tuple((row.selector, row.residual, row.scale) for row in result.specification.rows) == (
        ("pressure", "log_ratio", 0.05),
        ("vapor-solvent-fraction", "scaled_difference", 0.01),
        ("vapor-solvent-partial-pressure", "log_ratio", 0.05),
    )

    center_value = accepted.final_physical[0]
    step = 1.0e-5
    sides = tuple(
        equilibrium.solve(
            model,
            problem,
            active_parameters=epcsaft.ActiveParameterSet(
                parameters, (identity,), values=(center_value + sign * step,)
            ),
        )
        for sign in (-1.0, 1.0)
    )
    assert all(side.status == "evaluated" for side in sides)
    centered = tuple(
        (float(high) - float(low)) / (2.0 * step)
        for low, high in zip(sides[0].values, sides[1].values, strict=True)
    )
    owner_derivatives = tuple(
        residual_derivative
        * row.scale
        * (prediction if row.residual == "log_ratio" else 1.0)
        / row.multiplier
        for residual_derivative, row, prediction in zip(
            accepted.physical_jacobian,
            result.specification.rows,
            accepted.predictions,
            strict=True,
        )
    )
    for exact, reference in zip(owner_derivatives, centered, strict=True):
        assert scaled_error(exact, reference) <= 5.0e-5

    if ideal_vapor:
        center = equilibrium.solve(
            model,
            problem,
            active_parameters=epcsaft.ActiveParameterSet(
                parameters, (identity,), values=(center_value,)
            ),
        )
        assert center.status == "evaluated"
        temperature_k = float(problem.temperature.value.to("kelvin").magnitude)
        liquid = center.phases[0]
        fugacities = tuple(
            math.exp(liquid.chemical_potential_over_rt[index])
            * 8.31446261815324
            * temperature_k
            for index in (0, 1)
        )
        pressure_pa = math.fsum(fugacities)
        reference = (pressure_pa, fugacities[0] / pressure_pa, fugacities[0])
        assert center.values == pytest.approx(reference, rel=2.0e-10, abs=2.0e-8)


@pytest.mark.slow
@pytest.mark.parametrize("temperature_k", (313.15, 353.15, 393.15))
def test_solved_pressure_reaction_coordinate_certifies_full_space_path(
    tmp_path: Path,
    temperature_k: float,
) -> None:
    """Falsify reaction-coordinate derivatives, branch identity, or full-space evidence."""

    identity = "reaction:R1:correlation:b_k"
    mapping, parameters, model, problem, identity, _, _ = solved_pressure_fit_case(
        tmp_path,
        ideal_vapor=False,
        temperature_k=temperature_k,
        reaction_correlation=True,
        parameter_identity=identity,
        truth_value=-100.0,
        bounds=(-200.0, 0.0),
        origin=-100.0,
        scale=100.0,
    )
    result = fit(mapping, base_path=tmp_path)
    accepted = result.best_usable_start
    assert accepted is not None
    assert result.evaluation_strategy == "simultaneous_sparse_full_space"
    assert -200.0 <= accepted.final_physical[0] <= 0.0
    assert accepted.owner_failures == ()
    assert accepted.trial_failure_code is None
    assert accepted.batch_outcomes[0].status == "evaluated"
    assert all(row.evaluation_status == "evaluated" for row in result.evaluated_rows)
    assert tuple(row.selector for row in result.specification.rows) == (
        "pressure",
        "vapor-solvent-fraction",
        "vapor-solvent-partial-pressure",
        "liquid-cation-fraction",
    )

    center = equilibrium.solve(
        model,
        problem,
        active_parameters=epcsaft.ActiveParameterSet(
            parameters, (identity,), values=accepted.final_physical
        ),
    )
    assert center.status == "evaluated" and center.continuation_state is not None
    assert problem.continuation_identity == "finite-incipient-boundary"
    assert problem.continuation_state is not None
    for fingerprint in (
        problem.continuation_state.problem_fingerprint,
        problem.continuation_state.reaction_fingerprint,
        problem.continuation_state.topology_fingerprint,
        problem.continuation_state.parameter_fingerprint,
        center.continuation_state.problem_fingerprint,
        center.continuation_state.reaction_fingerprint,
        center.continuation_state.topology_fingerprint,
        center.continuation_state.parameter_fingerprint,
    ):
        assert fingerprint.startswith("sha256:")
    step = 10.0 if temperature_k == 393.15 else 2.0
    sides = tuple(
        equilibrium.solve(
            model,
            problem,
            active_parameters=epcsaft.ActiveParameterSet(
                parameters,
                (identity,),
                values=(accepted.final_physical[0] + sign * step,),
            ),
        )
        for sign in (-1.0, 1.0)
    )
    assert all(side.status == "evaluated" for side in sides)
    centered = tuple(
        (float(high) - float(low)) / (2.0 * step)
        for low, high in zip(sides[0].values, sides[1].values, strict=True)
    )
    owner_derivatives = tuple(
        residual_derivative
        * row.scale
        * (prediction if row.residual == "log_ratio" else 1.0)
        / row.multiplier
        for residual_derivative, row, prediction in zip(
            accepted.physical_jacobian,
            result.specification.rows,
            accepted.predictions,
            strict=True,
        )
    )
    assert all(
        scaled_error(actual, reference) <= 5.0e-5
        for actual, reference in zip(owner_derivatives, centered, strict=True)
    )

    receipt = json.loads(accepted.owner_receipts[0].canonical_bytes_utf8)
    branch = receipt["branch_signature"]
    expected_branch = json.loads(json.dumps(center.branch_signature))
    assert branch[0] == expected_branch[0]
    assert branch[1].startswith("sha256:")
    assert branch[2] == pytest.approx(expected_branch[2], rel=2.0e-10)
    assert branch[3] == pytest.approx(expected_branch[3], abs=2.0e-12)
    assert branch[4] == pytest.approx(expected_branch[4], abs=2.0e-12)
    assert tuple(coordinate.identity for coordinate in result.specification.coordinates) == (
        identity,
    )
    evidence = dict(receipt["evidence"])
    assert evidence["pressure_relative_inf_norm"] <= 1.0e-8
    assert evidence["reaction_affinity_inf_norm"] <= 1.0e-6
    assert accepted.kkt_diagnostics["retained"]
    assert accepted.kkt_diagnostics["primal_feasibility_inf_norm"] <= 1.0e-7
    assert accepted.kkt_diagnostics["stationarity_inf_norm"] <= 1.0e-7
    assert accepted.kkt_diagnostics["complementarity_inf_norm"] <= 1.0e-7
    assert accepted.kkt_diagnostics["dual_sign_violation_inf_norm"] <= 1.0e-7


@pytest.mark.slow
def test_solved_pressure_full_space_certifies_a_closed_declared_parameter_bound(
    tmp_path: Path,
) -> None:
    """Falsify affine-bound rounding or final certification at a closed bound."""

    identity = "pair/cosolvent/solvent/k_ij"
    mapping, _, _, _, _, truth_value, _ = solved_pressure_fit_case(
        tmp_path,
        ideal_vapor=False,
        parameter_identity=identity,
        truth_value=0.01,
        bounds=(-0.01, 0.01),
        origin=0.0,
        scale=0.01,
    )
    result = fit(mapping, base_path=tmp_path)
    accepted = result.best_usable_start
    assert accepted is not None
    assert result.evaluation_strategy == "simultaneous_sparse_full_space"
    assert accepted.final_physical == pytest.approx((truth_value,), abs=2.0e-8)
    assert accepted.owner_failures == ()


@pytest.mark.slow
def test_solved_pressure_full_space_rejects_only_physical_domain_trials(
    tmp_path: Path,
) -> None:
    """Falsify Ipopt restoration crashes when a trial leaves the EOS domain."""

    identity = "pair/cosolvent/solvent/k_ij"
    mapping, *_ = solved_pressure_fit_case(
        tmp_path,
        ideal_vapor=False,
        parameter_identity=identity,
        truth_value=0.01,
        bounds=(-0.01, 0.01),
        origin=0.0,
        scale=0.01,
    )
    coordinate = mapping["active"]["coordinates"][identity]
    coordinate.update({"bounds": [-1.0, 1.0], "scale": 1.0})
    observation = mapping["observations"][0]
    observation["targets"] = observation["targets"][:1]
    observation["targets"][0]["observed"] *= 100.0
    mapping["solver"]["controls"]["maximum_iterations"] = 20

    result = fit(mapping, base_path=tmp_path)
    accepted = result.best_usable_start
    assert result.evaluation_strategy == "simultaneous_sparse_full_space"
    assert accepted is not None
    assert accepted.solver_status == "solve_succeeded"
    assert accepted.owner_failures == ()
    assert accepted.final_physical[0] == pytest.approx(0.0157182582, abs=2.0e-8)


@pytest.mark.slow
def test_full_space_retains_a_typed_non_evaluable_reaction_trial(tmp_path: Path) -> None:
    """Falsify loss of a typed failed start from the full-space result."""

    mapping, *_ = solved_pressure_fit_case(
        tmp_path,
        ideal_vapor=False,
        temperature_k=393.15,
        reaction_correlation=True,
        reaction_b_upper=30_000.0,
        parameter_identity="reaction:R1:correlation:b_k",
        truth_value=-100.0,
        bounds=(-200.0, 30_000.0),
        origin=-100.0,
        scale=100.0,
    )
    mapping["solver"]["multistart"]["starts"] = [
        {"identity": "non-evaluable", "values": [30_000.0]}
    ]

    result = fit(mapping, base_path=tmp_path)
    refused = result.starts[0]
    assert result.evaluation_strategy == "simultaneous_sparse_full_space"
    assert not refused.usable
    assert refused.trial_failure_code == "branch_loss"
    assert refused.non_evaluable_count >= 1
    assert refused.owner_failures[0].code == refused.trial_failure_code
    assert result.accounting.non_evaluable_trials >= 1


def test_full_space_transports_constructor_and_selection_failures(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Falsify escaped construction failures or unbounded selected output rows."""

    mapping, *_ = full_space_fit_case(tmp_path)
    arguments = _native_full_space_payload(mapping, tmp_path, monkeypatch)
    native_fit = _core._fit_resolved
    coordinates = deepcopy(arguments[2])
    assert isinstance(coordinates, list)
    coordinates[0].update(
        origin=2.0, scale=3.0, lower_affine=-1.0, upper_affine=1.0
    )
    controls = deepcopy(arguments[4])
    assert isinstance(controls, dict)
    controls["starts"] = [
        {"identity": "owner-contract", "affine": [-1.0]},
        {"identity": "provider-domain", "affine": [1.0]},
    ]
    constructor = native_fit(
        arguments[0],
        arguments[1],
        coordinates,
        arguments[3],
        controls,
        arguments[5],
        arguments[6],
        "simultaneous_sparse_full_space",
    )
    assert constructor["preflight_failure"] is None
    assert [start["trial_failure_code"] for start in constructor["starts"]] == [
        "owner_contract_failure",
        "provider_domain_rejection",
    ]
    assert all(
        start["batch_outcomes"]
        == [
            {
                "index": 0,
                "row_count": 1,
                "status": "rejected",
                "failure_code": start["trial_failure_code"],
            }
        ]
        for start in constructor["starts"]
    )

    experiments = [dict(experiment) for experiment in arguments[6]]
    assert isinstance(experiments, list)
    experiments[0]["selected_rows"] = [1]
    for strategy in ("automatic", "simultaneous_sparse_full_space"):
        selection = native_fit(
            arguments[0],
            arguments[1],
            arguments[2],
            arguments[3],
            arguments[4],
            arguments[5],
            experiments,
            strategy,
        )
        assert selection["evaluation_strategy"] == "simultaneous_sparse_full_space"
        assert selection["starts"][0]["trial_failure_code"] == "owner_contract_failure"
        assert selection["starts"][0]["batch_outcomes"][0]["status"] == "rejected"


def test_full_space_without_a_final_iterate_keeps_batches_not_evaluated(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Falsify fabricated full-space penalties after an invalid owner objective."""

    mapping, *_ = full_space_fit_case(tmp_path)
    arguments = _native_full_space_payload(mapping, tmp_path, monkeypatch)
    observations = deepcopy(arguments[3])
    assert isinstance(observations, list)
    observations[0]["targets"][0].update(observed=1.0, scale=1.0e-300)
    native = _core._fit_resolved(
        arguments[0],
        arguments[1],
        arguments[2],
        observations,
        arguments[4],
        arguments[5],
        arguments[6],
        "simultaneous_sparse_full_space",
    )
    start = native["starts"][0]
    assert start["solver_status"] == "invalid_number_detected"
    assert start["trial_failure_code"] == "native_residual_failure"
    assert start["failure_reason"] == "ipopt_final_iterate_unavailable"
    assert start["batch_outcomes"] == [
        {
            "index": 0,
            "row_count": 1,
            "status": "not_evaluated",
            "failure_code": "native_residual_failure",
        }
    ]
    assert all(math.isnan(value) for value in start["predictions"])
    assert all(math.isnan(value) for value in start["residuals"])
    assert start["affine_jacobian"] == [0.0]


@pytest.mark.slow
def test_full_space_finalizes_later_rejected_batches(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Falsify first-rejection truncation of final full-space batch outcomes."""

    first_mapping, *_ = solved_pressure_fit_case(
        tmp_path,
        ideal_vapor=False,
        temperature_k=333.15,
        reaction_correlation=True,
        reaction_b_upper=30_000.0,
        parameter_identity="reaction:R1:correlation:b_k",
        truth_value=-100.0,
        bounds=(-200.0, 30_000.0),
        origin=-100.0,
        scale=100.0,
    )
    second_mapping, *_ = solved_pressure_fit_case(
        tmp_path,
        ideal_vapor=False,
        temperature_k=393.15,
        reaction_correlation=True,
        reaction_b_upper=30_000.0,
        parameter_identity="reaction:R1:correlation:b_k",
        truth_value=-100.0,
        bounds=(-200.0, 30_000.0),
        origin=-100.0,
        scale=100.0,
    )
    first = _native_full_space_payload(first_mapping, tmp_path, monkeypatch)
    second = _native_full_space_payload(second_mapping, tmp_path, monkeypatch)
    coordinates = deepcopy(first[2])
    controls = deepcopy(first[4])
    assert isinstance(controls, dict)
    controls["starts"] = [{"identity": "rejected-finalization", "affine": [301.0]}]
    controls["maximum_iterations"] = 1
    native = _core._fit_resolved(
        first[0],
        first[1],
        coordinates,
        [*first[3], *second[3]],
        controls,
        first[5],
        [*first[6], *second[6]],
        "simultaneous_sparse_full_space",
    )
    outcomes = native["starts"][0]["batch_outcomes"]
    assert [(item["index"], item["status"], item["failure_code"]) for item in outcomes] == [
        (0, "rejected", "physical_rejection"),
        (1, "rejected", "physical_rejection"),
    ]
    assert native["starts"][0]["owner_failures"] == [
        {
            "code": "physical_rejection",
            "reason": "reactive_finite_phase_certificate_failed",
            "count": 2,
        }
    ]
