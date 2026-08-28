from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

import epcsaft
from epcsaft import equilibrium, regression
from tests.eos._cases import STRONGLY_ASSOCIATING_BINARY
from tests.equilibrium.test_grepe_temperature_boundary import _nonreactive_problem


def test_regression_recovers_binary_parameter_from_declared_bubble_temperature(
    tmp_path: Path,
) -> None:
    """Falsify a GREPE temperature row that bypasses Regression or loses its exact J."""

    parameters = epcsaft.Parameters.from_mapping(STRONGLY_ASSOCIATING_BINARY)
    model = epcsaft.Mixture(parameters)
    units = epcsaft.unit_registry
    parameter_identity = "pair/solute/solvent/k_ij"
    truth_value = -0.02
    truth_active = epcsaft.ActiveParameterSet(
        parameters, (parameter_identity,), values=(truth_value,)
    )
    topology = equilibrium.DeclaredTopology(
        (
            equilibrium.PhaseSpecification("liquid", equilibrium.AllComponents()),
            equilibrium.PhaseSpecification("vapor", equilibrium.AllComponents()),
        )
    )
    reference = (0.5, 0.5)
    truth = equilibrium.solve(
        model,
        equilibrium.PhaseEquilibriumProblem(
            equilibrium.FixedTCoexistence(
                300.0 * units.kelvin, topology, 0, reference
            )
        ),
        active_parameters=truth_active,
    )
    phases = (
        equilibrium.ReactivePhase(
            "liquid",
            "liquid",
            "finite",
            equilibrium.AllComponents(),
            equilibrium.ProviderModel((1.0e-9, 0.74), "installed-provider-eos"),
        ),
        equilibrium.ReactivePhase(
            "vapor",
            "vapor",
            "incipient",
            equilibrium.AllComponents(),
            equilibrium.ProviderModel((1.0e-9, 0.74), "installed-provider-eos"),
        ),
    )
    request = equilibrium.GeneralReactiveEquilibriumProblem(
        "regression-declared-bubble-temperature",
        equilibrium.Solved(
            299.0 * units.kelvin,
            (295.0 * units.kelvin, 305.0 * units.kelvin),
            (297.0 * units.kelvin, 303.0 * units.kelvin),
        ),
        equilibrium.Fixed(truth.extent.pressure_pa * units.pascal),
        phases,
        _nonreactive_problem(parameters.component_ids, reference),
        (),
        (
            equilibrium.EquilibriumOutput(
                "bubble-temperature",
                "system.temperature",
                "kelvin",
                "absolute-temperature",
            ),
        ),
        continuation_identity="declared-bubble-temperature",
        feed=equilibrium.direct_amounts(parameters.component_ids, reference),
    )
    parameter_bytes = json.dumps(
        parameters.to_mapping(), separators=(",", ":"), sort_keys=True
    ).encode()
    parameter_path = tmp_path / "parameters.json"
    parameter_path.write_bytes(parameter_bytes)
    mapping = {
        "schema_version": 2,
        "parameters": {
            "path": str(parameter_path),
            "sha256": f"sha256:{hashlib.sha256(parameter_bytes).hexdigest()}",
        },
        "model": {"kind": "mixture"},
        "active": {
            "identities": [parameter_identity],
            "coordinates": {
                parameter_identity: {
                    "bounds": [-0.05, 0.02],
                    "origin": -0.04,
                    "scale": 0.1,
                }
            },
            "ties": [],
        },
        "observations": [
            {
                "family": "equilibrium",
                "identity": "declared-bubble-temperature-fit",
                "request": equilibrium.general_reactive_equilibrium_problem_to_mapping(request),
                "targets": [
                    {
                        "identity": "manufactured-bubble-temperature",
                        "prediction_identity": "bubble-temperature",
                        "observed": 300.0,
                        "unit": "kelvin",
                        "residual": "scaled_difference",
                        "scale": 1.0,
                        "multiplier": 1.0,
                        "role": "validation",
                        "measurement_classification": "direct",
                        "basis": "absolute-temperature",
                        "source_identity": "manufactured-strongly-associating-binary",
                        "source_hash": f"sha256:{'4' * 64}",
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
                "starts": [{"identity": "off-truth", "values": [-0.04]}],
                "affine_bound_margin": 1.0e-7,
            },
            "controls": {
                "maximum_iterations": 10,
                "optimizer_tolerance": 1.0e-8,
            },
        },
    }

    result = regression.fit(mapping, base_path=tmp_path)
    accepted = result.best_usable_start
    assert accepted is not None, result.to_mapping()
    assert accepted.final_physical == pytest.approx((truth_value,), abs=2.0e-7)
    assert accepted.predictions == pytest.approx((300.0,), abs=2.0e-8)
    assert accepted.owner_receipts
