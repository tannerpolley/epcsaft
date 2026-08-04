from __future__ import annotations

import math
from dataclasses import asdict, replace

import pytest

from epcsaft import Mixture, Parameters, unit_registry
from epcsaft.regression import (
    AcquisitionClass,
    AffineParameterTransform,
    AssociationParameterIdentity,
    ComponentParameterIdentity,
    ConfirmationControls,
    FixedCompositionMixtureDensityObservation,
    FixedTopologyAssociationCapability,
    ObjectiveContract,
    ObservationDataset,
    ObservationPartition,
    PairParameterIdentity,
    ParameterFamily,
    ParameterRequest,
    RankControls,
    RowProvenance,
    SolverControls,
    SourceInput,
    parameter_capabilities,
    prepare_fit,
)

_COMPONENTS = ("component-a", "component-b", "component-c")
_MOLAR_MASSES = (0.200, 0.400, 0.100)


def _pair_matrix(value: float) -> tuple[tuple[float, ...], ...]:
    return ((0.0, value), (value, 0.0))


def _neutral_model(*, k_ij: float = 0.0, l_ij: float = 0.0) -> Mixture:
    return Mixture(
        Parameters.from_dictionary(
            {
                "schema": "epcsaft.parameters",
                "schema_version": 1,
                "components": _COMPONENTS[:2],
                "parameters": {
                    "mw": (0.020, 0.040),
                    "m": (1.1, 1.4),
                    "s": (3.1, 3.4),
                    "e": (180.0, 220.0),
                    "k_ij": _pair_matrix(k_ij),
                    "l_ij": _pair_matrix(l_ij),
                },
                "options": {"permittivity_model": "none"},
                "validity": {
                    "kind": "reported-conditions",
                    "temperature_min_k": 200.0,
                    "temperature_max_k": 400.0,
                    "pressure_min_pa": 1.0,
                    "pressure_max_pa": 100.0e6,
                },
            }
        )
    )


def _associating_model(
    *,
    segment_count: float = 2.1,
    segment_diameter: float = 3.7,
    dispersion_energy_over_k: float = 310.0,
    association_energy_over_k: float = 1900.0,
) -> Mixture:
    return Mixture(
        Parameters.from_dictionary(
            {
                "schema": "epcsaft.parameters",
                "schema_version": 1,
                "components": _COMPONENTS,
                "parameters": {
                    "mw": _MOLAR_MASSES,
                    "m": (1.4, segment_count, 1.2),
                    "s": (3.2, segment_diameter, 3.0),
                    "e": (220.0, dispersion_energy_over_k, 180.0),
                    "k_ij": ((0.0, 0.0, 0.0),) * 3,
                    "sites": (
                        {
                            "component_id": "component-b",
                            "site_id": "donor",
                            "site_class": "generic",
                            "multiplicity": 1,
                        },
                        {
                            "component_id": "component-b",
                            "site_id": "acceptor",
                            "site_class": "generic",
                            "multiplicity": 1,
                        },
                    ),
                    "association": (
                        {
                            "component_id_a": "component-b",
                            "site_id_a": "donor",
                            "component_id_b": "component-b",
                            "site_id_b": "acceptor",
                            "association_energy_over_k": association_energy_over_k,
                            "association_volume": 0.025,
                        },
                    ),
                },
                "options": {"permittivity_model": "none"},
                "validity": {
                    "kind": "reported-conditions",
                    "temperature_min_k": 250.0,
                    "temperature_max_k": 360.0,
                    "pressure_min_pa": 1.0,
                    "pressure_max_pa": 100.0e6,
                },
            }
        )
    )


def _rows(
    target: Mixture,
    *,
    component_ids: tuple[str, ...],
    molar_masses: tuple[float, ...],
    states: tuple[tuple[float, float, tuple[float, ...]], ...],
    source_id: str,
) -> tuple[FixedCompositionMixtureDensityObservation, ...]:
    rows = []
    for index, (temperature, volume, fractions) in enumerate(states):
        pressure = float(
            target.state(
                T=temperature * unit_registry.kelvin,
                rho=(1.0 / volume) * unit_registry.mole / unit_registry.meter**3,
                x=fractions,
            )
            .pressure.to("pascal")
            .magnitude
        )
        molar_mass = math.fsum(
            fraction * mass for fraction, mass in zip(fractions, molar_masses, strict=True)
        )
        density = molar_mass / volume
        rows.append(
            FixedCompositionMixtureDensityObservation(
                row_id=f"{source_id}-{index}",
                source_id=source_id,
                source_locator=f"manufactured-row-{index}",
                component_ids=component_ids,
                mole_fractions=fractions,
                temperature_k=temperature,
                pressure_pa=pressure,
                density_kg_per_m3=density,
                molar_mass_kg_per_mol=molar_mass,
                pressure_scale_pa=pressure,
                density_scale_kg_per_m3=density,
                volume_origin_m3_per_mol=volume,
                volume_start_m3_per_mol=volume,
                volume_bounds_m3_per_mol=(0.8 * volume, 1.2 * volume),
                partition=ObservationPartition.TRAINING,
            )
        )
    return tuple(rows)


def _prepare(
    model: Mixture,
    rows: tuple[FixedCompositionMixtureDensityObservation, ...],
    requests: tuple[ParameterRequest, ...],
    starts: tuple[tuple[float, ...], ...],
):
    records = []
    provenance = {}
    for row in rows:
        record = asdict(row)
        record["partition"] = row.partition.value
        records.append(record)
        provenance[row.row_id] = RowProvenance(
            AcquisitionClass.DIRECT_MEASUREMENT,
            "unique manufactured row",
            "included",
            "not a critical-region exclusion",
            "not censored",
            "no outlier rule applied",
        )
    dataset = ObservationDataset.from_records(
        FixedCompositionMixtureDensityObservation,
        tuple(records),
        source=SourceInput(
            rows[0].source_id,
            "Anonymous manufactured fixed-composition mixture-density fixture",
            "tests/test_mixture_density_regression.py",
            "0" * 64,
            "public EOS states; no row selection",
            "SI; ordered mole fractions; mixture mass-density basis",
            "generic regression mechanics only",
            "observed pressure and density magnitudes",
        ),
        objective=ObjectiveContract(
            "fixed_composition_mixture_density",
            "native_scaled_least_squares",
            "observation_residual_scales",
            "independent_no_covariance",
            "squared",
            (),
            "fail_fit",
        ),
        row_provenance=provenance,
    )
    return prepare_fit(
        model,
        datasets=(dataset,),
        parameters=requests,
        parameter_slot_indices=tuple(range(len(requests))),
        start_vectors=starts,
        solver=SolverControls(150, 8.0, 1.0e-12, 1.0e-12, 1.0e-12),
        rank=RankControls(1.0e12),
        confirmation=ConfirmationControls(1.0e-6, 1.0e-8),
    )


def _assert_public_jacobian(prepared, physical: tuple[float, ...]) -> None:
    lifted = (0.0,) * len(prepared.problem.training_observations)
    evaluation = prepared.evaluate(physical, lifted_solver_point=lifted)
    parameter_count = len(physical)
    variable_count = parameter_count + len(lifted)
    assert evaluation == prepared.evaluate(physical, lifted_solver_point=lifted)
    assert evaluation.physical_parameter_point == physical
    assert evaluation.solver_parameter_point == tuple(
        coordinate.transform.to_solver(value)
        for coordinate, value in zip(prepared.problem.parameters, physical, strict=True)
    )
    assert evaluation.lifted_physical_point == tuple(
        row.volume_origin_m3_per_mol for row in prepared.problem.training_observations
    )
    assert evaluation.lifted_solver_point == lifted
    assert evaluation.fitted_variable_identities == prepared.problem.parameters
    assert evaluation.jacobian_layout == "row_major"
    assert len(evaluation.residual_vector) == 2 * len(lifted)
    assert len(evaluation.jacobian) == len(evaluation.residual_vector) * variable_count
    assert evaluation.lifted_variable_ids == tuple(
        f"{row.row_id}:volume_m3_per_mol" for row in prepared.problem.training_observations
    )
    assert tuple(item.component for item in evaluation.residual_identities) == (
        "pressure",
        "density",
    ) * len(lifted)
    assert evaluation.jacobian_diagnostics.full_rank == variable_count
    assert evaluation.jacobian_diagnostics.projected_parameter_rank == parameter_count
    assert math.isfinite(evaluation.jacobian_diagnostics.full_condition_number)
    assert math.isfinite(evaluation.jacobian_diagnostics.projected_parameter_condition_number)
    assert evaluation.parameter_fingerprints == tuple(
        dict.fromkeys(
            coordinate.provider_parameter_fingerprint for coordinate in prepared.problem.parameters
        )
    )
    assert evaluation.topology_fingerprints == tuple(
        dict.fromkeys(
            coordinate.provider_topology_fingerprint for coordinate in prepared.problem.parameters
        )
    )
    assert evaluation.capability_artifact_fingerprints == tuple(
        dict.fromkeys(
            coordinate.provider_artifact_fingerprint
            for coordinate in prepared.problem.parameters
            if coordinate.provider_artifact_fingerprint is not None
        )
    )
    assert evaluation.preparation_fingerprint == prepared.preparation_fingerprint
    assert evaluation.dataset_fingerprints == tuple(
        dataset.provenance_sha256 for dataset in prepared.datasets
    )
    assert evaluation.source_fingerprints == tuple(
        dataset.source.source_artifact_sha256 for dataset in prepared.datasets
    )
    assert evaluation.observation_order_fingerprint.startswith("sha256:")

    step = 1.0e-6
    for column in range(variable_count):
        lower_physical = list(physical)
        upper_physical = list(physical)
        lower_lifted = list(lifted)
        upper_lifted = list(lifted)
        if column < parameter_count:
            coordinate = prepared.problem.parameters[column]
            solver = evaluation.solver_parameter_point[column]
            lower_physical[column] = coordinate.transform.to_physical(solver - step)
            upper_physical[column] = coordinate.transform.to_physical(solver + step)
        else:
            lower_lifted[column - parameter_count] -= step
            upper_lifted[column - parameter_count] += step
        lower = prepared.evaluate(
            tuple(lower_physical), lifted_solver_point=tuple(lower_lifted)
        ).residual_vector
        upper = prepared.evaluate(
            tuple(upper_physical), lifted_solver_point=tuple(upper_lifted)
        ).residual_vector
        for row, (lower_value, upper_value) in enumerate(zip(lower, upper, strict=True)):
            assert evaluation.jacobian[row * variable_count + column] == pytest.approx(
                (upper_value - lower_value) / (2.0 * step),
                rel=4.0e-5,
                abs=4.0e-6,
            )


@pytest.mark.parametrize(
    ("family", "target_value"),
    ((ParameterFamily.K_IJ, 0.03), (ParameterFamily.L_IJ, 0.02)),
)
def test_neutral_pair_mixture_density_recovery_and_public_jacobian(
    family: ParameterFamily, target_value: float
) -> None:
    target = _neutral_model(
        k_ij=target_value if family is ParameterFamily.K_IJ else 0.0,
        l_ij=target_value if family is ParameterFamily.L_IJ else 0.0,
    )
    rows = _rows(
        target,
        component_ids=_COMPONENTS[:2],
        molar_masses=(0.020, 0.040),
        states=tuple(
            (temperature, volume, (0.35, 0.65))
            for temperature, volume in zip(
                (220.0, 240.0, 260.0, 280.0),
                (4.0e-4, 5.0e-4, 6.0e-4, 7.0e-4),
                strict=True,
            )
        ),
        source_id=f"manufactured-{family.value}",
    )
    prepared = _prepare(
        _neutral_model(),
        rows,
        (
            ParameterRequest(
                family,
                PairParameterIdentity(*_COMPONENTS[:2]),
                AffineParameterTransform(0.0, 0.05),
                -0.1,
                0.1,
            ),
        ),
        ((0.0,), (0.01,)),
    )
    assert prepared.preflight().ready
    _assert_public_jacobian(prepared, (0.0,))
    result = prepared.fit()
    assert result.workflow_valid
    assert result.parameters[0].final == pytest.approx(target_value, abs=2.0e-8)


def test_associating_mixture_recovers_three_generic_component_parameters() -> None:
    targets = (2.25, 3.82, 345.0)
    target = _associating_model(
        segment_count=targets[0],
        segment_diameter=targets[1],
        dispersion_energy_over_k=targets[2],
    )
    compositions = (
        (0.20, 0.60, 0.20),
        (0.25, 0.55, 0.20),
        (0.30, 0.50, 0.20),
        (0.20, 0.65, 0.15),
        (0.25, 0.60, 0.15),
        (0.30, 0.55, 0.15),
    )
    rows = _rows(
        target,
        component_ids=_COMPONENTS,
        molar_masses=_MOLAR_MASSES,
        states=tuple(
            (temperature, volume, composition)
            for temperature, volume, composition in zip(
                (275.0, 290.0, 305.0, 320.0, 335.0, 350.0),
                (1.8e-3, 2.0e-3, 2.2e-3, 2.4e-3, 2.6e-3, 2.8e-3),
                compositions,
                strict=True,
            )
        ),
        source_id="manufactured-associating-principal-three",
    )
    model = _associating_model(
        segment_count=2.1,
        segment_diameter=3.70,
        dispersion_energy_over_k=325.0,
    )
    requests = tuple(
        ParameterRequest(
            family,
            ComponentParameterIdentity("component-b"),
            transform,
            *bounds,
        )
        for family, transform, bounds in (
            (ParameterFamily.SEGMENT_COUNT, AffineParameterTransform(2.0, 0.2), (1.5, 3.0)),
            (ParameterFamily.SEGMENT_DIAMETER, AffineParameterTransform(3.6, 0.2), (3.0, 4.5)),
            (
                ParameterFamily.DISPERSION_ENERGY_OVER_K,
                AffineParameterTransform(300.0, 40.0),
                (200.0, 500.0),
            ),
        )
    )
    prepared = _prepare(
        model,
        rows,
        requests,
        ((2.1, 3.70, 325.0), (2.15, 3.75, 335.0)),
    )
    preflight = prepared.preflight()
    assert preflight.ready
    assert preflight.starts[0].projected_parameter_rank == 3
    _assert_public_jacobian(prepared, (2.1, 3.70, 325.0))
    result = prepared.fit()
    assert result.workflow_valid
    assert result.failed_row_count == 0
    assert result.jacobian.projected_parameter_rank == 3
    assert all(item.active_bound is None for item in result.parameters)
    assert tuple(item.final for item in result.parameters) == pytest.approx(
        targets, rel=2.0e-6, abs=2.0e-6
    )
    assert result.to_json_bytes(prepared=prepared) == result.to_json_bytes(prepared=prepared)


def test_associating_mixture_accepts_a_descriptor_owned_association_slot() -> None:
    model = _associating_model(association_energy_over_k=1900.0)
    target = _associating_model(association_energy_over_k=2300.0)
    rows = _rows(
        target,
        component_ids=_COMPONENTS,
        molar_masses=_MOLAR_MASSES,
        states=tuple(
            (temperature, 2.0e-3, (0.25, 0.55, 0.20))
            for temperature in (280.0, 300.0, 320.0, 340.0)
        ),
        source_id="manufactured-association-slot",
    )
    capability = next(
        item
        for item in parameter_capabilities(model)
        if isinstance(item, FixedTopologyAssociationCapability)
    )
    slot = next(
        item
        for item in capability.slots
        if item.family is ParameterFamily.ASSOCIATION_ENERGY_OVER_K
        and isinstance(item.identity, AssociationParameterIdentity)
    )
    prepared = _prepare(
        model,
        rows,
        (
            ParameterRequest(
                ParameterFamily.ASSOCIATION_ENERGY_OVER_K,
                slot.identity,
                AffineParameterTransform(1900.0, 200.0),
                1500.0,
                2600.0,
            ),
        ),
        ((1800.0,), (1900.0,)),
    )
    assert prepared.preflight().ready
    result = prepared.fit()
    assert result.workflow_valid
    assert result.parameters[0].final == pytest.approx(2300.0, abs=3.0e-4)

    with pytest.raises(ValueError, match="mole_fractions must sum to one"):
        replace(rows[0], mole_fractions=(0.25, 0.50, 0.20))
    with pytest.raises(RuntimeError, match="^eos_domain_failure:"):
        prepared.evaluate(
            (1900.0,),
            lifted_solver_point=(math.log(2.0),) * len(rows),
        )
    reversed_rows = tuple(
        replace(
            row,
            component_ids=tuple(reversed(row.component_ids)),
            mole_fractions=tuple(reversed(row.mole_fractions)),
        )
        for row in rows
    )
    incompatible = _prepare(
        model,
        reversed_rows,
        (
            ParameterRequest(
                ParameterFamily.ASSOCIATION_ENERGY_OVER_K,
                slot.identity,
                AffineParameterTransform(1900.0, 200.0),
                1500.0,
                2600.0,
            ),
        ),
        ((1800.0,), (1900.0,)),
    ).preflight()
    assert not incompatible.ready
    assert all(
        "does not match the installed fixed-topology association descriptor" in start.failure_reason
        for start in incompatible.starts
    )
