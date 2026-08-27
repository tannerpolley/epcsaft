from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

import epcsaft
from epcsaft import equilibrium, regression
from tests.eos._raw_cases import MULTINEUTRAL_REACTIVE


@pytest.mark.slow
def test_generic_regression_recovers_parameter_from_imposed_potential_uptake(
    tmp_path: Path,
) -> None:
    """Falsify a Regression path that cannot consume implicit GREPE outputs."""

    parameters = epcsaft.Parameters.from_mapping(MULTINEUTRAL_REACTIVE)
    model = epcsaft.Mixture(parameters)
    component_ids = parameters.component_ids
    reaction = equilibrium.ChemicalEquilibriumProblem(
        species_ids=component_ids,
        charges=(0, 0, 1, -1),
        molar_masses_kg_per_mol=(0.018, 0.030, 0.009, 0.009),
        balance_matrix=(
            (1.0, 0.0, 0.0, 0.0),
            (0.0, 1.0, 0.0, 0.0),
            (0.0, 0.0, 1.0, 0.0),
            (0.0, 0.0, 0.0, 1.0),
        ),
        conserved_totals=(0.8, 0.2, 0.1, 0.1),
        reaction_matrix=(),
        feed_amounts_mol=(0.8, 0.2, 0.1, 0.1),
        equilibrium_constants=(),
        strict_interior_amount_floor_mol=1.0e-14,
    )
    parameter_identity = "pair/anion/solvent/k_ij"
    truth_parameter = 0.08
    truth_active = epcsaft.ActiveParameterSet(
        parameters, (parameter_identity,), values=(truth_parameter,)
    )
    phase = equilibrium.ReactivePhase(
        "configured-liquid",
        "liquid",
        "finite",
        equilibrium.AllComponents(),
        equilibrium.ProviderModel((1.0e-6, 0.74), "installed-provider-eos"),
    )
    temperature = equilibrium.Fixed(313.15 * epcsaft.unit_registry.kelvin)
    pressure = equilibrium.Fixed(119_199.78746707023 * epcsaft.unit_registry.pascal)
    truth_parameters = parameters.with_values(
        {parameter_identity: truth_parameter},
        source_id="mixed-regression-truth",
        locator="project-authored numerical fixture",
    )
    direct_composition = (0.6, 0.2, 0.1, 0.1)
    direct_density = 500.0
    direct_pressure = float(
        epcsaft.Mixture(truth_parameters)
        .state(
            T=313.15 * epcsaft.unit_registry.kelvin,
            rho=direct_density * epcsaft.unit_registry.mole / epcsaft.unit_registry.meter**3,
            x=direct_composition,
        )
        .pressure.to("pascal")
        .magnitude
    )
    potential = equilibrium.EquilibriumOutput(
        "cosolvent-log-fugacity",
        "phase.log_fugacity",
        "dimensionless",
        "natural-log-pascal-true-species",
        phase.identity,
        (0.0, 1.0, 0.0, 0.0),
    )
    fixed_truth = equilibrium.GeneralReactiveEquilibriumProblem(
        "nonreactive-open-component-truth",
        temperature,
        pressure,
        (phase,),
        reaction,
        (),
        (potential,),
        feed=equilibrium.direct_amounts(component_ids, (0.8, 0.2, 0.1, 0.1)),
    )
    truth = equilibrium.solve(
        model,
        fixed_truth,
        active_parameters=truth_active,
    )
    assert truth.status == "evaluated", truth.failure
    problem = equilibrium.GeneralReactiveEquilibriumProblem(
        "nonreactive-open-component-fit",
        temperature,
        pressure,
        (phase,),
        reaction,
        (),
        (
            equilibrium.EquilibriumOutput(
                "cosolvent-uptake",
                "phase.amount",
                "mole",
                "true-species-amount",
                phase.identity,
                (0.0, 1.0, 0.0, 0.0),
                support="positive",
            ),
        ),
        feed=equilibrium.AffineFeed(
            component_ids,
            (0.8, 0.05, 0.1, 0.1),
            ((0.0,), (1.0,), (0.0,), (0.0,)),
            (equilibrium.FeedCoordinate("open-neutral", (0.02, 0.5)),),
        ),
        intensive_boundaries=(
            equilibrium.IntensiveBoundary(
                "imposed-neutral-potential",
                "open-neutral",
                phase.identity,
                "phase.log_fugacity",
                potential.coefficients,
                float(truth.values[0]),
                0.15,
            ),
        ),
    )
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
                parameter_identity: {
                    "bounds": [-0.05, 0.15],
                    "origin": 0.0,
                    "scale": 0.1,
                }
            },
            "ties": [],
        },
        "observations": [
            {
                "family": "fixed_state_pressure",
                "composition": list(direct_composition),
                "source": {
                    "kind": "table",
                    "path": "inline-mixed-direct-pressure",
                    "sha256": f"sha256:{hashlib.sha256(json.dumps([{'temperature_k': 313.15, 'molar_density_mol_per_m3': direct_density, 'pressure': direct_pressure}], separators=(',', ':'), sort_keys=True).encode()).hexdigest()}",
                    "citation": "Project-authored mixed-owner numerical fixture.",
                    "domain": {
                        "temperature_k": [312.0, 314.0],
                        "molar_density_mol_per_m3": [499.0, 501.0],
                    },
                    "rows": [
                        {
                            "temperature_k": 313.15,
                            "molar_density_mol_per_m3": direct_density,
                            "pressure": direct_pressure,
                        }
                    ],
                },
                "targets": ["pressure"],
                "residuals": {"pressure": "log_ratio"},
                "units": {"pressure": "pascal"},
                "basis": {"pressure": "absolute"},
                "reference_convention": {"pressure": "absolute"},
                "transformation": {"pressure": {"input": "identity", "output": "identity"}},
                "multiplier": 1.0,
            },
            {
                "family": "equilibrium",
                "identity": "open-component-fit-batch",
                "request": equilibrium.general_reactive_equilibrium_problem_to_mapping(problem),
                "targets": [
                    {
                        "identity": "manufactured-open-component-uptake",
                        "prediction_identity": "cosolvent-uptake",
                        "observed": 0.2,
                        "unit": "mole",
                        "residual": "scaled_difference",
                        "scale": 0.1,
                        "multiplier": 1.0,
                        "role": "diagnostic",
                        "measurement_classification": "direct",
                        "basis": "true-species-amount",
                        "source_identity": "manufactured-open-component",
                        "source_hash": f"sha256:{'4' * 64}",
                        "included": True,
                        "drop_reason": None,
                        "aggregate_identity": None,
                        "covariance_identity": None,
                    }
                ],
                "multiplier": 1.0,
            },
        ],
        "solver": {
            "multistart": {
                "starts": [{"identity": "origin", "values": [0.0]}],
                "affine_bound_margin": 1.0e-7,
            },
            "controls": {
                "maximum_iterations": 12,
                "optimizer_tolerance": 1.0e-8,
            },
        },
    }
    default_mapping = json.loads(json.dumps(mapping))
    del default_mapping["solver"]["controls"]
    default_controls = regression.ParameterFit.from_mapping(default_mapping).solver.controls
    assert default_controls.maximum_iterations == 25
    assert default_controls.optimizer_tolerance == 1.0e-6

    result = regression.fit(mapping, base_path=tmp_path)

    accepted = result.best_usable_start
    assert accepted is not None
    assert accepted.final_physical == pytest.approx((truth_parameter,), abs=3.0e-7)
    assert accepted.predictions == pytest.approx((direct_pressure, 0.2), rel=2.0e-12, abs=3.0e-8)
    assert accepted.owner_receipts
    assert tuple(row.selector for row in result.evaluated_rows) == (
        "pressure",
        "cosolvent-uptake",
    )
    fit_receipt = result.specification.to_mapping()["fit"]
    assert fit_receipt["fit"]["solver"]["controls"] == mapping["solver"]["controls"]
