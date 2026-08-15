from __future__ import annotations

import hashlib
import json
from dataclasses import replace
from pathlib import Path

import pytest

import epcsaft
from epcsaft import equilibrium, regression
from epcsaft.eos import _native_transport
from epcsaft.equilibrium import _equilibrium as _native_equilibrium
from epcsaft.equilibrium import _grepe
from tests.eos._raw_cases import NEUTRAL_TERNARY_IMMISCIBLE

_BASIS = "manufactured-reactive-lle-basis"
_PARAMETER_IDENTITY = "pair/alpha/beta/k_ij"


def _fixture() -> tuple[
    epcsaft.Parameters,
    epcsaft.Mixture,
    equilibrium.GeneralReactiveEquilibriumProblem,
]:
    parameters = epcsaft.Parameters.from_mapping(NEUTRAL_TERNARY_IMMISCIBLE)
    model = epcsaft.Mixture(parameters)
    seed = equilibrium.phase_equilibrium(
        model,
        equilibrium.PhaseEquilibriumProblem(
            equilibrium.FixedTPFeed(
                300.0 * epcsaft.unit_registry.kelvin,
                1.0e5 * epcsaft.unit_registry.pascal,
                (1.0, 1.0, 1.0),
                equilibrium.DeclaredTopology(
                    tuple(
                        equilibrium.PhaseSpecification("liquid", equilibrium.AllComponents())
                        for _ in range(3)
                    )
                ),
            )
        ),
    )
    aqueous, organic = seed.phases[:2]
    starts = (aqueous.mole_fractions, organic.mole_fractions)
    feed = tuple(starts[0][index] + starts[1][index] for index in range(3))
    ln_k = aqueous.chemical_potential_over_rt[1] - aqueous.chemical_potential_over_rt[0]
    reaction = equilibrium.ChemicalEquilibriumProblem(
        species_ids=model.component_ids,
        charges=(0, 0, 0),
        molar_masses_kg_per_mol=(0.03, 0.03, 0.03),
        balance_matrix=((1.0, 1.0, 0.0), (0.0, 0.0, 1.0)),
        conserved_totals=(feed[0] + feed[1], feed[2]),
        reaction_matrix=((-1.0, 1.0, 0.0),),
        feed_amounts_mol=feed,
        equilibrium_constants=(
            equilibrium.ChemicalEquilibriumConstant(
                ln_k,
                "manufactured-reactive-lle-oracle",
                "provider-helmholtz-coordinate-basis",
                "products_positive",
                "already-provider-basis",
                True,
            ),
        ),
        strict_interior_amount_floor_mol=1.0e-12,
    )
    phases = tuple(
        equilibrium.ReactivePhase(
            identity,
            "liquid",
            "finite",
            equilibrium.AllComponents(),
            equilibrium.ProviderModel((1.0e-8, 0.74), "installed-provider-eos"),
            start=equilibrium.FinitePhaseStart(
                amounts,
                state.molar_volume_m3_per_mol,
            ),
        )
        for identity, amounts, state in zip(
            ("aqueous", "organic"),
            starts,
            (aqueous, organic),
            strict=True,
        )
    )
    outputs = tuple(
        equilibrium.EquilibriumOutput(
            f"amount-{phase.identity}-{component_id}",
            "phase.amount",
            "mole",
            _BASIS,
            phase.identity,
            tuple(float(index == component) for index in range(3)),
            support="positive",
        )
        for phase in phases
        for component, component_id in enumerate(model.component_ids)
    ) + (
        equilibrium.EquilibriumOutput(
            "aqueous-alpha-mole-fraction",
            "phase.mole_fraction",
            "dimensionless",
            _BASIS,
            "aqueous",
            (1.0, 0.0, 0.0),
            support="positive",
        ),
        equilibrium.EquilibriumOutput(
            "organic-molar-density",
            "phase.molar_density",
            "mole / meter**3",
            _BASIS,
            "organic",
            support="positive",
        ),
        equilibrium.EquilibriumOutput(
            "aqueous-pressure",
            "phase.pressure",
            "pascal",
            _BASIS,
            "aqueous",
            support="positive",
        ),
        equilibrium.EquilibriumOutput(
            "aqueous-beta-chemical-potential",
            "phase.chemical_potential_over_rt",
            "dimensionless",
            _BASIS,
            "aqueous",
            (0.0, 1.0, 0.0),
        ),
    )
    return (
        parameters,
        model,
        equilibrium.GeneralReactiveEquilibriumProblem(
            "manufactured-reactive-lle",
            equilibrium.Fixed(300.0 * epcsaft.unit_registry.kelvin),
            equilibrium.Fixed(1.0e5 * epcsaft.unit_registry.pascal),
            phases,
            reaction,
            ("aqueous", "organic"),
            outputs,
        ),
    )


def _reaction_mapping(
    problem: equilibrium.ChemicalEquilibriumProblem,
) -> dict[str, object]:
    return {
        "species_ids": list(problem.species_ids),
        "charges": list(problem.charges),
        "molar_masses_kg_per_mol": list(problem.molar_masses_kg_per_mol),
        "balance_matrix": [list(row) for row in problem.balance_matrix],
        "conserved_totals": list(problem.conserved_totals),
        "reaction_matrix": [list(row) for row in problem.reaction_matrix],
        "feed_amounts_mol": list(problem.feed_amounts_mol),
        "equilibrium_constants": [
            [
                item.ln_value,
                item.source_id,
                item.reference_id,
                item.reaction_orientation,
                item.conversion_id,
                item.dimensionless,
            ]
            for item in problem.equilibrium_constants
        ],
        "strict_interior_amount_floor_mol": problem.strict_interior_amount_floor_mol,
        "source_standard_state": None,
    }


def _request_mapping(
    problem: equilibrium.GeneralReactiveEquilibriumProblem,
) -> dict[str, object]:
    return {
        "identity": problem.identity,
        "temperature": {"role": "fixed", "value": 300.0, "unit": "kelvin"},
        "pressure": {"role": "fixed", "value": 1.0e5, "unit": "pascal"},
        "phases": [
            {
                "identity": phase.identity,
                "fluid_role": phase.fluid_role,
                "amount_role": phase.amount_role,
                "support": {"kind": "all_components", "component_ids": []},
                "model": {
                    "kind": "provider",
                    "reference_id": "installed-provider-eos",
                    "admissible_packing_fraction_interval": list(
                        phase.model.admissible_packing_fraction_interval
                    ),
                },
                "continuation": None,
                "start": {
                    "amounts_mol": list(phase.start.amounts_mol),
                    "molar_volume_m3_per_mol": phase.start.molar_volume_m3_per_mol,
                },
            }
            for phase in problem.phases
            if isinstance(phase.model, equilibrium.ProviderModel) and phase.start is not None
        ],
        "reaction_system": _reaction_mapping(problem.reaction_system),
        "reaction_phase_ids": list(problem.reaction_phase_ids),
        "outputs": [
            {
                "identity": output.identity,
                "selector": output.selector,
                "unit": output.unit,
                "basis": output.basis,
                "phase_identity": output.phase_identity,
                "coefficients": list(output.coefficients),
                "solvent_mass_coefficients_kg_per_mol": [],
                "support": output.support,
                "censor_limit": output.censor_limit,
                "aggregate_identity": output.aggregate_identity,
                "covariance_identity": output.covariance_identity,
            }
            for output in problem.outputs
        ],
        "continuation": None,
    }


def test_reactive_lle_is_one_certified_coupled_minimum_and_label_covariant() -> None:
    """Falsify staged speciation, missing transfer equilibrium, or label ordering."""

    parameters, model, problem = _fixture()
    active = epcsaft.ActiveParameterSet(parameters, (_PARAMETER_IDENTITY,))
    result = equilibrium.phase_equilibrium(model, problem, active_parameters=active)

    assert result.status == "evaluated", result.failure
    assert result.numerical_status == "passed"
    assert result.physical_status == "passed"
    assert result.provider_domain_status == "passed"
    assert result.descriptor is not None
    assert result.descriptor._owner_kind == "reactive_lle"
    assert result.certificate_scope == "fixed_TP_two_finite_reacting_liquids_local_minimum"
    assert result.globality_claim == "not_established"
    assert len(result.phases) == 2
    phase_amounts = {
        phase.identity: tuple(
            result.rows[phase_index * 3 + component].value for component in range(3)
        )
        for phase_index, phase in enumerate(problem.phases)
    }
    assert all(value is not None for values in phase_amounts.values() for value in values)
    reconstructed_species = tuple(
        sum(float(values[component]) for values in phase_amounts.values()) for component in range(3)
    )
    reconstructed_conserved = (
        reconstructed_species[0] + reconstructed_species[1],
        reconstructed_species[2],
    )
    assert reconstructed_conserved == pytest.approx(
        problem.reaction_system.conserved_totals, abs=2.0e-10
    )
    aqueous, organic = result.phases
    for left, right in zip(
        aqueous.chemical_potential_over_rt,
        organic.chemical_potential_over_rt,
        strict=True,
    ):
        assert left == pytest.approx(right, abs=3.0e-8)
    for phase in result.phases:
        assert (
            -phase.chemical_potential_over_rt[0] + phase.chemical_potential_over_rt[1]
        ) == pytest.approx(0.0, abs=3.0e-8)
        assert phase.pressure_pa == pytest.approx(1.0e5, rel=7.0e-9)
    assert (
        max(
            abs(left - right)
            for left, right in zip(aqueous.mole_fractions, organic.mole_fractions, strict=True)
        )
        > 0.5
    )

    observations = tuple(
        equilibrium.GeneralEquilibriumOwnerObservation(output.identity, output.identity)
        for output in problem.outputs
    )
    receipt = equilibrium.serialize_general_reactive_equilibrium_receipt(
        result,
        batch_identity="manufactured/reactive-lle",
        observations=observations,
    )
    assert equilibrium.validate_general_reactive_equilibrium_receipt(
        receipt,
        descriptor=result.descriptor,
        problem=problem,
        batch_identity="manufactured/reactive-lle",
        active_parameters=active,
        observations=observations,
    ).startswith("sha256:")
    certificate = json.loads(receipt)["owner_receipt"]["certificate"]
    assert certificate["declared_direction_rank"] == certificate["admissible_nullity"] == 4
    assert certificate["reduced_hessian_min_eigenvalue"] > 0.0

    permuted = replace(
        problem,
        phases=tuple(reversed(problem.phases)),
        reaction_phase_ids=tuple(reversed(problem.reaction_phase_ids)),
    )
    swapped = equilibrium.phase_equilibrium(model, permuted, active_parameters=active)
    assert swapped.status == "evaluated", swapped.failure
    states = {
        phase.identity: state for phase, state in zip(problem.phases, result.phases, strict=True)
    }
    swapped_states = {
        phase.identity: state for phase, state in zip(permuted.phases, swapped.phases, strict=True)
    }
    for identity in states:
        assert swapped_states[identity].mole_fractions == pytest.approx(
            states[identity].mole_fractions, abs=2.0e-10
        )


def test_reactive_lle_total_parameter_jacobian_matches_independent_resolves() -> None:
    """Falsify a partial-only LLE Jacobian or fixed-pressure derivative leakage."""

    parameters, model, problem = _fixture()
    center_value = 0.2
    step = 1.0e-5

    def solve(value: float) -> equilibrium.GeneralReactiveEquilibriumResult:
        result = equilibrium.phase_equilibrium(
            model,
            problem,
            active_parameters=epcsaft.ActiveParameterSet(
                parameters,
                (_PARAMETER_IDENTITY,),
                values=(value,),
            ),
        )
        assert result.status == "evaluated", result.failure
        return result

    center = solve(center_value)
    below = solve(center_value - step)
    above = solve(center_value + step)
    for row, low, high in zip(center.rows, below.values, above.values, strict=True):
        assert row.jacobian is not None and low is not None and high is not None
        numerical = (high - low) / (2.0 * step)
        scale = max(1.0, abs(row.jacobian[0]), abs(numerical))
        assert abs(row.jacobian[0] - numerical) / scale < 3.0e-6


@pytest.mark.slow
def test_reactive_lle_exact_nlp_derivatives_match_centered_oracles() -> None:
    """Falsify wrong objective, constraint, or Lagrangian derivatives supplied to IPOPT."""

    parameters, model, problem = _fixture()
    active = epcsaft.ActiveParameterSet(parameters, (_PARAMETER_IDENTITY,))
    native = _native_equilibrium._reactive_finite_phases(
        _native_transport(model, active),
        _grepe._reactive_lle_spec(problem, 300.0, 1.0e5),
        tuple(
            phase.model.admissible_packing_fraction_interval
            for phase in problem.phases
            if isinstance(phase.model, equilibrium.ProviderModel)
        ),
        tuple(
            {
                "amounts_mol": phase.start.amounts_mol,
                "molar_volume_m3_per_mol": phase.start.molar_volume_m3_per_mol,
            }
            for phase in problem.phases
            if isinstance(phase.start, equilibrium.FinitePhaseStart)
        ),
        problem.reaction_system.strict_interior_amount_floor_mol,
        True,
        True,
    )

    assert native["accepted"], native
    assert native["derivative_audits_performed"]
    assert native["derivative_audit_status"] == "passed"
    assert native["objective_gradient_audit_max_relative_error"] < 2.0e-6
    assert native["constraint_jacobian_audit_max_relative_error"] < 2.0e-8
    assert native["lagrangian_hessian_audit_max_relative_error"] < 2.0e-5


def test_reactive_lle_provider_domain_failure_preserves_failure_axis() -> None:
    """Falsify collapsing a Provider source-domain rejection into nonconvergence."""

    _, model, problem = _fixture()
    first_phase = problem.phases[0]
    assert isinstance(first_phase.model, equilibrium.ProviderModel)
    outside_source_domain = replace(
        problem,
        phases=(
            replace(
                first_phase,
                model=replace(
                    first_phase.model,
                    admissible_packing_fraction_interval=(1.0e-8, 1.0e-7),
                ),
            ),
            problem.phases[1],
        ),
    )
    result = equilibrium.phase_equilibrium(model, outside_source_domain)

    assert result.status == "non_evaluable"
    assert result.failure is not None
    assert result.failure.code == "provider_domain_rejection"
    assert result.solver_status == "not_started"
    assert result.numerical_status == "not_adjudicated"
    assert result.physical_status == "not_adjudicated"
    assert result.provider_domain_status == "failed"


def test_reactive_lle_strict_mapping_and_unsupported_charge_fail_closed() -> None:
    parameters, model, problem = _fixture()
    with pytest.raises(ValueError, match="both a start and continuation"):
        replace(
            problem.phases[0],
            continuation_identity="ambiguous-local-anchor",
            branch_policy="local_certified_role_selected_state",
        )
    mapping = _request_mapping(problem)
    loaded = equilibrium.general_reactive_equilibrium_problem_from_mapping(mapping)
    assert loaded == problem

    mapping["phases"][0]["start"]["amounts_mol"][0] = 0.0
    with pytest.raises(ValueError, match="positive amounts"):
        equilibrium.general_reactive_equilibrium_problem_from_mapping(mapping)

    charged = replace(
        problem,
        reaction_system=replace(problem.reaction_system, charges=(1, 1, 0)),
    )
    rejected = equilibrium.phase_equilibrium(model, charged)
    assert rejected.status == "non_evaluable"
    assert rejected.failure is not None
    assert rejected.failure.code == "unsupported_topology"
    assert rejected.numerical_status == "unsupported"
    assert rejected.physical_status == "unsupported"
    assert rejected.provider_domain_status == "unsupported"
    assert "charged reactive LLE" in rejected.failure.diagnostic
    with pytest.raises(ValueError, match="charged reactive LLE"):
        equilibrium.general_reactive_equilibrium_descriptor(model, charged)


def test_generic_regression_recovers_reactive_lle_parameter(tmp_path: Path) -> None:
    """Falsify any hidden LLE-specific Regression objective or owner receipt."""

    parameters, model, problem = _fixture()
    true_value = 0.202
    truth = equilibrium.phase_equilibrium(
        model,
        problem,
        active_parameters=epcsaft.ActiveParameterSet(
            parameters,
            (_PARAMETER_IDENTITY,),
            values=(true_value,),
        ),
    )
    assert truth.status == "evaluated", truth.failure
    observed = truth.values[0]
    assert observed is not None
    parameter_payload = json.dumps(
        parameters.to_mapping(), separators=(",", ":"), sort_keys=True
    ).encode()
    parameter_path = tmp_path / "reactive-lle-parameters.json"
    parameter_path.write_bytes(parameter_payload)
    fit_mapping = {
        "schema_version": 2,
        "parameters": {
            "path": str(parameter_path),
            "sha256": f"sha256:{hashlib.sha256(parameter_payload).hexdigest()}",
        },
        "model": {"kind": "mixture"},
        "active": {
            "identities": [_PARAMETER_IDENTITY],
            "coordinates": {
                _PARAMETER_IDENTITY: {
                    "bounds": [0.19, 0.21],
                    "origin": 0.2,
                    "scale": 0.005,
                }
            },
            "ties": [],
        },
        "observations": [
            {
                "family": "equilibrium",
                "identity": "manufactured/reactive-lle-regression",
                "request": _request_mapping(problem),
                "targets": [
                    {
                        "identity": "manufactured/reactive-lle-alpha-amount",
                        "prediction_identity": problem.outputs[0].identity,
                        "observed": observed,
                        "unit": "mole",
                        "residual": "scaled_difference",
                        "scale": observed,
                        "multiplier": 1.0,
                        "role": "diagnostic",
                        "measurement_classification": "manufactured",
                        "basis": _BASIS,
                        "source_identity": "manufactured-reactive-lle-oracle",
                        "source_hash": f"sha256:{'3' * 64}",
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
                "starts": [{"identity": "offset-start", "values": [0.198]}],
                "affine_bound_margin": 1.0e-7,
            }
        },
    }
    fitted = regression.fit(fit_mapping, base_path=tmp_path)
    accepted = fitted.best_usable_start
    assert accepted is not None, fitted.to_mapping()
    assert accepted.final_physical == pytest.approx((true_value,), abs=2.0e-8)
    assert accepted.predictions == pytest.approx((observed,), rel=2.0e-9)
    assert len(accepted.owner_receipts) == 1
    assert accepted.owner_receipts[0].digest.startswith("sha256:")
