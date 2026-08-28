from __future__ import annotations

import copy
import hashlib
import json
import math
from collections.abc import Mapping, Sequence
from dataclasses import replace
from typing import Any, cast
from unittest.mock import patch

import pytest

import epcsaft
from epcsaft import _core, equilibrium, regression
from epcsaft.equilibrium._grepe import _native_general_reactive_equilibrium_inputs
from epcsaft.regression._equilibrium_owner import _equilibrium_owner_batches
from epcsaft.regression._model import ParameterFit
from tests.eos._cases import MULTINEUTRAL_REACTIVE
from tests.equilibrium.test_chemical_provider import _multineutral_problem
from tests.regression.numerical._monitoring import call_without_python_frames


def _source_standard_problem(
    model: epcsaft.Mixture,
    problem: equilibrium.ChemicalEquilibriumProblem,
) -> equilibrium.ChemicalEquilibriumProblem:
    temperature_k = 313.15
    pressure_pa = 119_199.78746707023
    declaration = epcsaft.SourceReferenceDeclaration(
        reference_state_id="multineutral-full-space-source",
        component_ids=problem.species_ids,
        solvent_mole_fractions=(0.8, 0.2, 0.0, 0.0),
        ion_pairs=(epcsaft.SourceReferenceIonPair(("cation", "anion"), (1, 1)),),
        phase="liquid",
        reference_convention="mole-fraction-mixed-reference",
        activity_convention_id="mole-fraction-mixed-reference-v1",
        standard_molality_mol_per_kg=1.0,
        reference_pressure_pa=pressure_pa,
        required_derivatives=(),
    )
    transfer = epcsaft.source_reference_transfer(
        model,
        declaration,
        T=temperature_k * epcsaft.unit_registry.kelvin,
        P=pressure_pa * epcsaft.unit_registry.pascal,
    )
    solvent = transfer.neutral_basis.index((1.0, 0.0, 0.0, 0.0))
    ion_pair = transfer.neutral_basis.index((0.0, 0.0, 1.0, 1.0))
    offset = (
        -transfer.transfer_log_contractions[solvent]
        + transfer.transfer_log_contractions[ion_pair]
        + math.log(pressure_pa / (8.31446261815324 * temperature_k))
    )
    constant = problem.equilibrium_constants[0]
    return replace(
        problem,
        equilibrium_constants=(
            replace(
                constant,
                ln_value=constant.ln_value - offset,
                reference_id=declaration.reference_state_id,
                conversion_id="source-standard-state-to-provider-neutral-reference",
            ),
        ),
        source_standard_state=equilibrium.ChemicalStandardState(
            id=declaration.reference_state_id,
            activity_scale_id=declaration.activity_convention_id,
            log_activity_scale_factors=(0.0, 0.0, 0.0, 0.0),
            reference_pressure_pa=pressure_pa,
            source_reference_component_ids=problem.species_ids,
            source_reference_solvent_composition=(0.8, 0.2, 0.0, 0.0),
            source_reference_ion_pairs=declaration.ion_pairs,
            source_reference_phase="liquid",
            source_reference_convention="mole-fraction-mixed-reference",
            source_reference_activity_convention_id=declaration.activity_convention_id,
            source_reference_standard_molality_mol_per_kg=1.0,
        ),
    )


@pytest.mark.slow
def test_full_space_regression_recovers_shared_parameter_and_reactive_state(tmp_path) -> None:
    """Run the retained GREPE derivative, N=1, and generic N=8 acceptance case.

    Benchmark from ``engine/`` with::

        /usr/bin/time -f 'wall_time_seconds=%e' env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 ../.venv/bin/python -m pytest -q -m slow tests/regression/numerical/test_full_space_regression.py
    """
    parameters = epcsaft.Parameters.from_mapping(MULTINEUTRAL_REACTIVE)
    model = epcsaft.Mixture(parameters)
    reaction = _source_standard_problem(model, _multineutral_problem())
    temperature = 313.15 * epcsaft.unit_registry.kelvin
    pressure = 119_199.78746707023 * epcsaft.unit_registry.pascal
    phase_id = "multineutral-liquid"
    output = equilibrium.EquilibriumOutput(
        "cation-mole-fraction",
        "phase.mole_fraction",
        "dimensionless",
        "true-species-liquid-mole-fraction",
        phase_id,
        (0.0, 0.0, 1.0, 0.0),
        support="positive",
    )
    total = math.fsum(reaction.feed_amounts_mol)
    state = model.state(
        T=temperature,
        P=pressure,
        x=tuple(amount / total for amount in reaction.feed_amounts_mol),
        phase="liquid",
    )
    phase = equilibrium.ReactivePhase(
        phase_id,
        "liquid",
        "finite",
        equilibrium.AllComponents(),
        equilibrium.ProviderModel((1.0e-6, 0.74), "installed-provider-eos"),
        start=equilibrium.FinitePhaseStart(
            reaction.feed_amounts_mol,
            1.0 / float(state.molar_density.to("mole / meter**3").magnitude),
        ),
    )
    unanchored = equilibrium.GeneralReactiveEquilibriumProblem(
        "manufactured-full-space",
        equilibrium.Fixed(temperature),
        equilibrium.Fixed(pressure),
        (phase,),
        reaction,
        (phase_id,),
        (output,),
    )
    anchor = equilibrium.solve(model, unanchored)
    assert anchor.status == "evaluated", anchor.failure
    assert anchor.continuation_state is not None
    problem = replace(
        unanchored,
        continuation_identity="multineutral-liquid-branch",
        continuation_state=anchor.continuation_state,
        phases=(
            replace(
                phase,
                start=None,
            ),
        ),
    )
    parameter_identity = "component/solvent/segment_diameter"
    truth_parameter = 2.9
    truth = equilibrium.solve(
        model,
        problem,
        active_parameters=epcsaft.ActiveParameterSet(
            parameters, (parameter_identity,), values=(truth_parameter,)
        ),
    )
    assert truth.status == "evaluated", truth.failure
    observed = float(truth.values[0])

    payload = json.dumps(parameters.to_mapping(), separators=(",", ":"), sort_keys=True).encode()
    parameter_path = tmp_path / "parameters.json"
    parameter_path.write_bytes(payload)
    mapping = {
        "schema_version": 2,
        "parameters": {
            "path": str(parameter_path),
            "sha256": f"sha256:{hashlib.sha256(payload).hexdigest()}",
        },
        "model": {"kind": "mixture"},
        "active": {
            "identities": [parameter_identity],
            "coordinates": {
                parameter_identity: {"bounds": [2.5, 3.1], "origin": 2.8, "scale": 0.2}
            },
            "ties": [],
        },
        "observations": [
            {
                "family": "equilibrium",
                "identity": problem.identity,
                "request": equilibrium.general_reactive_equilibrium_problem_to_mapping(problem),
                "targets": [
                    {
                        "identity": "manufactured-cation",
                        "prediction_identity": output.identity,
                        "observed": observed,
                        "unit": output.unit,
                        "residual": "scaled_difference",
                        "scale": observed,
                        "multiplier": 1.0,
                        "role": "diagnostic",
                        "measurement_classification": "direct",
                        "basis": output.basis,
                        "source_identity": "manufactured-full-space",
                        "source_hash": f"sha256:{'2' * 64}",
                        "included": True,
                        "drop_reason": None,
                        "aggregate_identity": None,
                        "covariance_identity": None,
                    }
                ],
                "multiplier": 1.0,
            }
        ],
        "solver": {
            "multistart": {
                "starts": [{"identity": "origin", "values": [2.8]}],
                "affine_bound_margin": 1.0e-7,
            },
            "controls": {"maximum_iterations": 50, "optimizer_tolerance": 1.0e-8},
        },
    }

    resolved = ParameterFit.from_mapping(mapping).resolve(parameters, base_path=tmp_path)
    _, batches, *_ = _equilibrium_owner_batches(resolved)
    batch = batches[0]
    active = epcsaft.ActiveParameterSet(parameters, (parameter_identity,), values=(2.8,))
    _, native_problem, fingerprint = _native_general_reactive_equilibrium_inputs(
        batch.model, batch.problem, active_parameters=active
    )
    initialized = _core._initialize_grepe_local(
        batch.model._native, native_problem, active, fingerprint
    ).equilibrium
    local = tuple(initialized.local_variables)
    lower = tuple(initialized.local_lower_bounds)
    upper = tuple(initialized.local_upper_bounds)

    def evaluate(parameter: float, values: tuple[float, ...]) -> object:
        trial = epcsaft.ActiveParameterSet(parameters, (parameter_identity,), values=(parameter,))
        return _core._evaluate_grepe_local(
            batch.model._native,
            native_problem,
            trial,
            fingerprint,
            values,
            lower,
            upper,
        ).equilibrium

    center = evaluate(2.8, local)
    assert center.accepted
    assert len(local) == len(center.kkt_root_values)
    assert reaction.source_standard_state is not None

    constraint_columns: list[tuple[float, ...]] = []
    output_columns: list[tuple[float, ...]] = []
    parameter_step = 1.0e-5
    parameter_below = evaluate(2.8 - parameter_step, local)
    parameter_above = evaluate(2.8 + parameter_step, local)
    constraint_columns.append(
        tuple(
            (high - low) / (2.0 * parameter_step)
            for low, high in zip(
                parameter_below.kkt_root_values,
                parameter_above.kkt_root_values,
                strict=True,
            )
        )
    )
    output_columns.append(
        tuple(
            (high - low) / (2.0 * parameter_step)
            for low, high in zip(
                parameter_below.output_values,
                parameter_above.output_values,
                strict=True,
            )
        )
    )
    input_column = tuple(center.input_order).index(parameter_identity)
    input_width = len(center.input_order)
    parameter_column = tuple(center.parameter_order).index(parameter_identity)
    parameter_width = len(center.parameter_order)
    exact_constraint_columns = [
        tuple(
            center.kkt_input_jacobian[row * input_width + input_column] for row in range(len(local))
        )
    ]
    exact_output_columns = [
        tuple(
            center.output_parameter_partials[row * parameter_width + parameter_column]
            for row in range(len(center.output_values))
        )
    ]
    for column, value in enumerate(local):
        step = 1.0e-5 * max(1.0, abs(value))
        below_values = tuple(
            item - step if index == column else item for index, item in enumerate(local)
        )
        above_values = tuple(
            item + step if index == column else item for index, item in enumerate(local)
        )
        below = evaluate(2.8, below_values)
        above = evaluate(2.8, above_values)
        constraint_columns.append(
            tuple(
                (high - low) / (2.0 * step)
                for low, high in zip(below.kkt_root_values, above.kkt_root_values, strict=True)
            )
        )
        output_columns.append(
            tuple(
                (high - low) / (2.0 * step)
                for low, high in zip(below.output_values, above.output_values, strict=True)
            )
        )
        exact_constraint_columns.append(
            tuple(center.kkt_root_jacobian[row * len(local) + column] for row in range(len(local)))
        )
        exact_output_columns.append(
            tuple(
                center.output_local_derivatives[row * len(local) + column]
                for row in range(len(center.output_values))
            )
        )

    assert all(
        math.isfinite(value)
        for block in (
            center.kkt_root_values,
            center.kkt_root_jacobian,
            center.kkt_input_jacobian,
            center.output_values,
            center.output_local_derivatives,
            center.output_parameter_partials,
        )
        for value in block
    )
    for exact_columns, centered_columns in (
        (exact_constraint_columns, constraint_columns),
        (exact_output_columns, output_columns),
    ):
        for exact, centered in zip(exact_columns, centered_columns, strict=True):
            for expected, actual in zip(exact, centered, strict=True):
                scale = max(2.0e-9, abs(expected), abs(actual))
                assert abs(expected - actual) / scale <= 5.0e-5

    native_fit = _core._fit_resolved
    captured: dict[str, object] = {}
    frames: tuple[object, ...] = ()
    native_arguments: tuple[object, ...] = ()

    def monitored_fit(*args: object) -> object:
        nonlocal frames, native_arguments
        native_arguments = args
        value, frames = call_without_python_frames(native_fit, *args)
        captured.update(value)
        return value

    with patch.object(_core, "_fit_resolved", side_effect=monitored_fit):
        result = regression.fit(mapping, base_path=tmp_path)
    assert frames == ()
    assert captured["evaluation_strategy"] == "simultaneous_sparse_full_space"
    malformed_experiments = [
        dict(cast(Mapping[str, object], experiment))
        for experiment in cast(Sequence[object], native_arguments[6])
    ]
    malformed_experiments[0]["selected_rows"] = ()
    with pytest.raises(ValueError, match="one row per observation target"):
        native_fit(*native_arguments[:6], malformed_experiments)
    accepted = result.best_usable_start
    assert accepted is not None
    assert accepted.final_physical == pytest.approx((truth_parameter,), abs=2.0e-7)
    assert accepted.predictions == pytest.approx((observed,), rel=2.0e-6, abs=2.0e-8)
    assert accepted.owner_failures == ()

    retained_provenance: dict[str, object] | None = None

    def assert_owner_receipts(start: object, count: int) -> None:
        nonlocal retained_provenance
        receipts = cast(Sequence[object], cast(Any, start).owner_receipts)
        assert len(receipts) == count
        for opaque in receipts:
            receipt = json.loads(cast(Any, opaque).canonical_bytes_utf8)
            assert receipt["certificate_scope"] == truth.certificate_scope
            assert receipt["branch_signature"] == json.loads(json.dumps(truth.branch_signature))
            provenance = {
                key: receipt[key]
                for key in (
                    "descriptor_fingerprint",
                    "problem_fingerprint",
                    "reaction_fingerprint",
                    "topology_fingerprint",
                    "reference_fingerprint",
                    "domain_fingerprint",
                    "formulation",
                    "certificate_scope",
                    "globality_claim",
                    "resolved_policy",
                    "branch_signature",
                )
            }
            if retained_provenance is None:
                retained_provenance = provenance
            else:
                assert provenance == retained_provenance
            assert {key: value["status"] for key, value in receipt["outcomes"].items()} == {
                "solver_convergence": "passed",
                "numerical_convergence": "passed",
                "physical_certification": "passed",
                "scientific_fit_acceptance": "not_applicable",
            }
            evidence = dict(receipt["evidence"])
            assert evidence["balance_inf_norm"] < 1.0e-9
            assert evidence["pressure_relative_inf_norm"] < 1.0e-8
            assert evidence["reaction_affinity_inf_norm"] < 1.0e-6
            assert evidence["kkt_stationarity_inf_norm"] < 1.0e-7
            mole_fractions = receipt["phases"][0]["mole_fractions"]
            assert abs(mole_fractions[2] - mole_fractions[3]) < 1.0e-9

    assert_owner_receipts(accepted, 1)
    assert all(outcome.status == "evaluated" for outcome in accepted.batch_outcomes)
    assert result.accounting.dropped_rows == 0
    assert result.accounting.rejected_row_trials == 0
    assert result.accounting.failed_row_trials == 0

    eight_mapping = copy.deepcopy(mapping)
    eight_mapping["observations"] = []
    for index in range(8):
        observation = copy.deepcopy(mapping["observations"][0])
        observation["identity"] = f"manufactured-full-space-{index}"
        observation["targets"][0]["identity"] = f"manufactured-cation-{index}"
        observation["targets"][0]["source_identity"] = f"manufactured-full-space-{index}"
        eight_mapping["observations"].append(observation)
    captured_eight: dict[str, object] = {}
    eight_frames: tuple[object, ...] = ()

    def capture_eight(*args: object) -> object:
        nonlocal eight_frames
        value, eight_frames = call_without_python_frames(native_fit, *args)
        captured_eight.update(value)
        return value

    with patch.object(_core, "_fit_resolved", side_effect=capture_eight):
        eight_result = regression.fit(eight_mapping, base_path=tmp_path)
    assert eight_frames == ()
    assert captured_eight["evaluation_strategy"] == "simultaneous_sparse_full_space"
    eight_accepted = eight_result.best_usable_start
    assert eight_accepted is not None
    assert eight_accepted.final_physical == pytest.approx((truth_parameter,), abs=2.0e-7)
    assert eight_accepted.predictions == pytest.approx((observed,) * 8, rel=2.0e-6, abs=2.0e-8)
    assert eight_accepted.owner_failures == ()
    assert_owner_receipts(eight_accepted, 8)
    assert all(outcome.status == "evaluated" for outcome in eight_accepted.batch_outcomes)
    assert eight_result.accounting.dropped_rows == 0
    assert eight_result.accounting.rejected_row_trials == 0
    assert eight_result.accounting.failed_row_trials == 0
