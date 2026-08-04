from __future__ import annotations

import pytest
from parameter_cases import neutral_parameters

from epcsaft import Mixture, Parameters
from epcsaft.regression import (
    AcquisitionClass,
    AffineParameterTransform,
    ConfirmationControls,
    FixedCompositionVleObservation,
    ObjectiveContract,
    ObservationDataset,
    PairParameterIdentity,
    ParameterFamily,
    ParameterRequest,
    RankControls,
    RowProvenance,
    SolverControls,
    SourceInput,
    prepare_fit,
    usability,
)

SHA = "a" * 64


def _source() -> SourceInput:
    return SourceInput(
        source_id="doi:example",
        citation="Example et al. (2026)",
        durable_locator="https://doi.org/10.0000/example",
        source_artifact_sha256=SHA,
        transformation_record="bar to Pa; no exclusions",
        units_and_bases="T/K, P/Pa, mole fractions",
        use_basis="source-backed test data",
        residual_scale_rationale="pressure by observed pressure; mu/RT",
    )


def _objective() -> ObjectiveContract:
    return ObjectiveContract(
        residual_family="fixed_composition_vle",
        interpretation="native_scaled_least_squares",
        row_weighting="observation_residual_scales",
        covariance_interpretation="independent_no_covariance",
        loss="squared",
        loss_parameters=(),
        failed_row_treatment="fail_fit",
    )


def _records() -> tuple[dict[str, object], ...]:
    return (
        {
            "row_id": "row-1",
            "source_id": "doi:example",
            "source_locator": "table-1:row-1",
            "component_ids": ["methane", "ethane"],
            "temperature_k": 203.22,
            "pressure_pa": 2_124_000.0,
            "liquid_mole_fraction_first": 0.3653,
            "vapor_mole_fraction_first": 0.8667,
            "pressure_scale_pa": 2_124_000.0,
            "chemical_potential_scales": [1.0, 1.0],
            "liquid_volume_origin_m3_per_mol": 6.0e-5,
            "liquid_volume_start_m3_per_mol": 6.5e-5,
            "liquid_volume_bounds_m3_per_mol": [2.0e-5, 1.0e-4],
            "vapor_volume_origin_m3_per_mol": 9.0e-4,
            "vapor_volume_start_m3_per_mol": 1.0e-3,
            "vapor_volume_bounds_m3_per_mol": [1.0e-4, 1.0e-2],
            "partition": "training",
        },
    )


def _dataset() -> ObservationDataset:
    return ObservationDataset.from_records(
        FixedCompositionVleObservation,
        _records(),
        source=_source(),
        objective=_objective(),
        row_provenance={
            "row-1": RowProvenance(
                acquisition=AcquisitionClass.DIRECT_MEASUREMENT,
                duplicate_decision="unique source row",
                exclusion_decision="included",
                critical_region_decision="outside critical region",
                censoring_decision="not censored",
                outlier_decision="retained; no outlier rule applied",
            )
        },
    )


def _model() -> Mixture:
    return Mixture(Parameters.from_dictionary(neutral_parameters(("methane", "ethane"))))


def _prepared(dataset: ObservationDataset | None = None):
    model = _model()
    request = ParameterRequest(
        family=ParameterFamily.K_IJ,
        identity=PairParameterIdentity("methane", "ethane"),
        transform=AffineParameterTransform(0.0, 0.01),
        lower_bound=-0.15,
        upper_bound=0.10,
    )
    return prepare_fit(
        model,
        datasets=(dataset or _dataset(),),
        parameters=(request,),
        parameter_slot_indices=(0,),
        start_vectors=((0.0,), (-0.05,)),
        solver=SolverControls(
            maximum_iterations=200,
            maximum_solver_time_seconds=30.0,
            function_tolerance=1.0e-12,
            gradient_tolerance=1.0e-12,
            parameter_tolerance=1.0e-12,
        ),
        rank=RankControls(maximum_condition_number=1.0e10),
        confirmation=ConfirmationControls(
            parameter_scaled_max_delta=1.0e-5,
            cost_relative_delta=1.0e-8,
        ),
    )


def test_prepare_and_preflight_resolve_capability_and_exact_rank(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    prepared = _prepared()
    problem = prepared.problem
    support = usability.parameter_capabilities(prepared.model)
    capability = next(
        item for item in support if getattr(item, "family", None) is ParameterFamily.K_IJ
    )

    assert problem.parameters[0].capability_id == capability.capability_id
    assert problem.parameters[0].provider_parameter_fingerprint == (
        capability.parameter_fingerprint
    )
    assert problem.parameters[0].unit == "1"
    assert problem.maximum_iterations == 200
    assert problem.maximum_condition_number == 1.0e10
    assert capability.installed_ready
    assert capability.unsupported_reason == ""

    calls: list[tuple[float, ...]] = []

    def evaluate(problem: object, model: object, variables: tuple[float, ...]):
        calls.append(variables)
        # R=4, N=1, Q=2; all three columns are independent.
        return (
            (0.0, 0.0, 0.0, 0.0),
            (1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0),
        )

    monkeypatch.setattr(usability, "_evaluate_parameters", evaluate)
    report = prepared.preflight()

    assert len(calls) == 2
    assert (report.fitted_parameter_count, report.lifted_variable_count) == (1, 2)
    assert report.residual_count == 4
    assert report.starts[0].full_rank == 3
    assert report.starts[0].projected_parameter_rank == 1
    assert report.starts[0].fitted_bounds_status == ("interior",)
    assert report.starts[0].lifted_bounds_status == ("interior", "interior")
    assert report.ready
    assert report.reasons == ()


def test_public_scalar_evaluation_matches_centered_difference_oracle() -> None:
    prepared = _prepared()
    evaluation = prepared.evaluate((0.0,))
    solver_point = (
        *evaluation.solver_parameter_point,
        *evaluation.lifted_solver_point,
    )
    step = 1.0e-6

    assert evaluation.jacobian_layout == "row_major"
    assert evaluation.jacobian_diagnostics.full_rank == len(solver_point)
    for column in range(len(solver_point)):
        lower = list(solver_point)
        upper = list(solver_point)
        lower[column] -= step
        upper[column] += step
        lower_evaluation = prepared.evaluate(
            (prepared.problem.parameters[0].transform.to_physical(lower[0]),),
            lifted_solver_point=tuple(lower[1:]),
        )
        upper_evaluation = prepared.evaluate(
            (prepared.problem.parameters[0].transform.to_physical(upper[0]),),
            lifted_solver_point=tuple(upper[1:]),
        )
        for row, (lower_residual, upper_residual) in enumerate(
            zip(
                lower_evaluation.residual_vector,
                upper_evaluation.residual_vector,
                strict=True,
            )
        ):
            centered_difference = (upper_residual - lower_residual) / (2.0 * step)
            assert evaluation.jacobian[row * len(solver_point) + column] == pytest.approx(
                centered_difference, rel=2.0e-5, abs=2.0e-6
            )


@pytest.mark.parametrize(
    ("error", "reason"),
    (
        (
            ValueError("installed capability derivative changed"),
            "unavailable_derivatives:",
        ),
        (
            RuntimeError("EOS does not expose the associating-mixture interface"),
            "unavailable_derivatives:",
        ),
        (RuntimeError("density iteration left the EOS domain"), "eos_domain_failure:"),
        (
            RuntimeError("observation payload has wrong field count"),
            "malformed_rows:",
        ),
        (TypeError("row payload shape changed"), "malformed_rows:"),
        (ValueError("component order changed"), "incompatible_contract:"),
        (ArithmeticError("nonfinite exact Jacobian"), "nonfinite_evaluation:"),
    ),
)
def test_preflight_preserves_failure_classes(
    monkeypatch: pytest.MonkeyPatch,
    error: Exception,
    reason: str,
) -> None:
    def fail(*_args: object) -> None:
        raise error

    monkeypatch.setattr(usability, "_evaluate_parameters", fail)
    report = _prepared().preflight()

    assert not report.ready
    assert all(start.failure_reason.startswith(reason) for start in report.starts)
    with pytest.raises(RuntimeError, match=f"^{reason}"):
        _prepared().evaluate((0.0,))


@pytest.mark.parametrize(
    ("residuals", "jacobian", "reason"),
    (
        ((0.0,) * 4, (0.0,) * 11, "unavailable_derivatives:"),
        ((0.0, 0.0, 0.0, float("nan")), (0.0,) * 12, "nonfinite_evaluation:"),
        ((0.0,) * 4, (0.0,) * 11 + (float("nan"),), "nonfinite_evaluation:"),
    ),
)
def test_public_evaluation_rejects_incomplete_or_nonfinite_payload(
    monkeypatch: pytest.MonkeyPatch,
    residuals: tuple[float, ...],
    jacobian: tuple[float, ...],
    reason: str,
) -> None:
    monkeypatch.setattr(
        usability,
        "_evaluate_parameters",
        lambda *_: (residuals, jacobian),
    )

    with pytest.raises(RuntimeError, match=f"^{reason}"):
        _prepared().evaluate((0.0,))


def test_public_evaluation_classifies_malformed_and_diagnostic_failures(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(
        usability,
        "_evaluate_parameters",
        lambda *_: (None, (0.0,) * 12),
    )
    with pytest.raises(RuntimeError, match="^malformed_evaluation:"):
        _prepared().evaluate((0.0,))

    monkeypatch.setattr(
        usability,
        "_evaluate_parameters",
        lambda *_: ((0.0,) * 4, (0.0,) * 12),
    )

    def fail(*_args: object) -> None:
        raise ValueError("diagnostic payload mismatch")

    monkeypatch.setattr(usability._native, "diagnose_jacobian", fail)
    with pytest.raises(RuntimeError, match="^diagnostic_failure:"):
        _prepared().evaluate((0.0,))


@pytest.fixture(scope="module")
def installed_fit():
    prepared = _prepared()
    return prepared, prepared.preflight(), prepared.fit()
