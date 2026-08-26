from __future__ import annotations

import json
import os
from copy import deepcopy

import pytest

import epcsaft
from epcsaft import equilibrium, regression
from epcsaft.regression._equilibrium_owner import _equilibrium_owner_batches
from epcsaft.regression._fit import _fit_owner_batches
from epcsaft.regression._model import ParameterFit
from tests.eos._raw_cases import MULTINEUTRAL_REACTIVE
from tests.equilibrium.test_chemical_provider import _multineutral_problem
from tests.equilibrium.test_homogeneous_reactive_observations import (
    _observation_problem,
    _reactive_fit_mapping,
    _target_mapping,
)


@pytest.mark.slow
def test_full_space_regression_recovers_shared_parameter_and_reactive_states(
    tmp_path,
) -> None:
    """Falsify a Regression route that does not solve shared GREPE states simultaneously."""

    parameters = epcsaft.Parameters.from_mapping(MULTINEUTRAL_REACTIVE)
    model = epcsaft.Mixture(parameters)
    problem = _multineutral_problem()
    phase = equilibrium.ProviderPhase(model, (1.0e-6, 0.74))
    owner_problem = _observation_problem(phase, problem)
    parameter_identity = "pair/anion/solvent/k_ij"
    truth_parameter = 0.08
    row = equilibrium.HomogeneousReactiveObservationRow(
        "cation-mole-fraction",
        "mole_fraction",
        (0.0, 0.0, 1.0, 0.0),
        support="positive",
    )
    truth_active = epcsaft.ActiveParameterSet(
        parameters, (parameter_identity,), values=(truth_parameter,)
    )
    truth = equilibrium.evaluate_homogeneous_reactive_observations(
        owner_problem,
        313.15 * epcsaft.unit_registry.kelvin,
        119_199.78746707023 * epcsaft.unit_registry.pascal,
        equilibrium.homogeneous_reactive_observation_descriptor(
            owner_problem, (row,), active_parameters=truth_active
        ),
        active_parameters=truth_active,
    )
    assert truth.status == "evaluated", truth.failure
    observed = float(truth.values[0])
    mapping = _reactive_fit_mapping(
        tmp_path,
        parameters,
        owner_problem,
        (row,),
        (
            _target_mapping(
                "manufactured-cation-0",
                row.identity,
                observed,
                row.unit,
                "scaled_difference",
                observed,
            ),
        ),
        parameter_identity,
        lower=-0.05,
        upper=0.2,
        origin=0.0,
        scale=0.1,
        starts=(("origin", 0.0),),
    )
    experiment_count = int(os.environ.get("EPCSAFT_ISSUE51_EXPERIMENTS", "1"))
    template = mapping["observations"][0]
    observations = []
    for index in range(experiment_count):
        observation = deepcopy(template)
        observation["identity"] = f"manufactured-batch-{index}"
        observation["request"]["identity"] = f"manufactured-problem-{index}"
        observation["targets"][0]["identity"] = f"manufactured-cation-{index}"
        observations.append(observation)
    mapping["observations"] = observations
    mapping["solver"]["controls"] = {
        "maximum_iterations": 50,
        "optimizer_tolerance": 1.0e-8,
    }

    result = regression.fit(mapping, base_path=tmp_path)

    accepted = result.best_usable_start
    assert accepted is not None
    assert accepted.final_physical == pytest.approx((truth_parameter,), abs=2.0e-7)
    assert accepted.predictions == pytest.approx((observed,) * experiment_count, abs=2.0e-10)
    assert len(accepted.owner_receipts) == experiment_count
    assert accepted.cost <= 1.0e-14
    if experiment_count == 1:
        resolved = ParameterFit.from_mapping(mapping).resolve(parameters, base_path=tmp_path)
        _, _, owner, residuals, receipt = _equilibrium_owner_batches(resolved)
        reduced = _fit_owner_batches(resolved, owner, residuals, receipt).best_usable_start
        assert reduced is not None
        assert accepted.final_physical == pytest.approx(reduced.final_physical, abs=2.0e-8)
        assert accepted.cost == pytest.approx(reduced.cost, abs=2.0e-14)
        assert accepted.predictions == pytest.approx(reduced.predictions, abs=2.0e-10)
        assert accepted.physical_jacobian == pytest.approx(
            reduced.physical_jacobian, rel=2.0e-8, abs=2.0e-10
        )
        for full_receipt, reduced_receipt in zip(
            accepted.owner_receipts, reduced.owner_receipts, strict=True
        ):
            full_payload = json.loads(full_receipt.canonical_bytes_utf8)
            reduced_payload = json.loads(reduced_receipt.canonical_bytes_utf8)
            for key in (
                "certificate_scope",
                "globality_claim",
                "outcomes",
                "provider_domain_status",
                "resolved_policy",
                "branch_signature",
            ):
                assert full_payload[key] == reduced_payload[key]
    if "EPCSAFT_ISSUE51_EXPERIMENTS" in os.environ:
        print(f"N={experiment_count} iterations={accepted.iterations} cost={accepted.cost:.12g}")
