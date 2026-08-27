from __future__ import annotations

import hashlib
import json
import math
from dataclasses import replace
from unittest.mock import patch

import pytest

import epcsaft
from epcsaft import _core, equilibrium, regression
from epcsaft.regression._equilibrium_owner import _equilibrium_owner_batches, _full_space_batches
from epcsaft.regression._model import ParameterFit
from tests.eos._raw_cases import MULTINEUTRAL_REACTIVE
from tests.equilibrium.test_chemical_provider import _multineutral_problem


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
    full_space = _full_space_batches(resolved, batches)
    assert full_space is not None
    evaluator, setups = full_space
    setup = setups[0]
    point = [0.0, *setup["initial"]]
    step = 1.0e-4
    for checked_column in (0, 1):
        below = evaluator(
            0,
            [
                value - step if index == checked_column else value
                for index, value in enumerate(point)
            ],
        )
        above = evaluator(
            0,
            [
                value + step if index == checked_column else value
                for index, value in enumerate(point)
            ],
        )
        difference = [
            (high - low) / (2.0 * step)
            for low, high in zip(below["constraints"], above["constraints"], strict=True)
        ]
        exact = [0.0] * len(difference)
        for row, column, value in zip(
            setup["jacobian_rows"],
            setup["jacobian_columns"],
            evaluator(0, point)["constraint_jacobian"],
            strict=True,
        ):
            if column == checked_column:
                exact[row] = value
        assert exact == pytest.approx(difference, rel=2.0e-5, abs=2.0e-7)

    with patch.object(_core, "_fit_resolved", wraps=_core._fit_resolved) as native_fit:
        result = regression.fit(mapping, base_path=tmp_path)
    assert native_fit.call_args.args[7] is not None
    accepted = result.best_usable_start
    assert accepted is not None
    assert accepted.final_physical == pytest.approx((truth_parameter,), abs=2.0e-7)
    assert accepted.predictions == pytest.approx((observed,), abs=2.0e-10)
