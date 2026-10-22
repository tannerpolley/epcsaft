"""Measurement declarations, with or without a joint correlation, and the native worker transport preserve the same numerical problem."""
from typing import Any

from epcsaft import regression as reg
from epcsaft.equilibrium import Amounts, EOS_STANDARD_STATE_ID, Phase, Problem, Reaction, ReactionLogPolynomial
from tests.native_inputs import ETHANE, parameters


def test_shared_measurement_transport_and_worker_equivalence(evidence: Any) -> None:
    p = parameters((ETHANE,) * 2, pairs=[("k_ij", (0, 1), 0.0)])
    rows = []
    uncertainty = reg.MeasurementUncertainty()
    groups = []
    for i, (temperature, value) in enumerate(((300.0, 1.52), (320.0, 1.59), (340.0, 1.72))):
        problem = Problem(T=temperature, P=5e6, phases=[Phase("liquid", kind="liquid")],
            feed=Amounts({"s0": 0.5, "s1": 0.5}), reactions=[Reaction({"s0": -1.0, "s1": 1.0}, name="isomerisation",
                correlation=ReactionLogPolynomial(0.0, 0.0, 0.0, 0.0045, 300.0, 250.0, 450.0, EOS_STANDARD_STATE_ID, e=900.0, f=1e-7))])
        row = reg.observation(p, "pressure", problem, observed=value, form="difference", scale=1.0, mass_basis=False)
        row.kind, row.reaction = reg.ObservationKind.LogEquilibriumConstant, "isomerisation"
        rows.append(row)
        groups.append(reg.measurement_group(p, [i], [reg.measured_input(p, "temperature", sigma=2.0,
            bounds=(250.0, 450.0), unit="K")]))
    uncertainty.groups, uncertainty.output_sigma = groups, [0.02] * 3
    coordinate = reg.reaction_coordinate("isomerisation", "d", origin=0.0, scale=0.001, bounds=(0.001, 0.01))
    fitted = [reg.fit(p, [coordinate], rows, uncertainty=uncertainty,
        controls=reg.FitControls(maximum_iterations=100, maximum_elapsed_time_seconds=600.0, workers=workers))
        for workers in (1, 2)]
    a, b = fitted
    assert a.usable and b.usable and a.status == b.status == "converged"
    assert a.reconciliation is not None and b.reconciliation is not None
    assert a.physical == b.physical and a.predictions == b.predictions and a.final_cost == b.final_cost
    assert a.optimizer_jacobian == b.optimizer_jacobian
    assert a.reconciliation.values == b.reconciliation.values and a.reconciliation.jacobian == b.reconciliation.jacobian
    assert a.reconciliation.column_groups == b.reconciliation.column_groups == [0, 1, 2]
    assert a.reconciliation.column_inputs == b.reconciliation.column_inputs == [0, 0, 0]
    assert a.covariance.optimizer == b.covariance.optimizer and a.covariance.variance_factor is None
    evidence.observe("measurement-workers", physical=a.physical, reconciled=a.reconciliation.values,
        output_cost=a.reconciliation.output_cost, input_cost=a.reconciliation.input_cost, workers=(1, 2))
    named = reg.measured_input(p, "mole_fraction", species="s1", sigma=0.001, bounds=(0.1, 0.9), unit="mol/mol")
    chart = reg.measurement_group(p, [0], [named], dependent_species=["s0"], feed_total_mol=1.0)
    assert chart.inputs[0].species == 1 and chart.dependent_species == [0]
    # A declared correlation over (y, T) is admitted in the parent and whitened there; workers return identical rows.
    uncertainty.correlation = [0.5 ** abs(i - j) for i in range(6) for j in range(6)]
    a, b = (reg.fit(p, [coordinate], rows, uncertainty=uncertainty,
        controls=reg.FitControls(maximum_iterations=100, maximum_elapsed_time_seconds=600.0, workers=workers))
        for workers in (1, 2))
    assert a.status == b.status == "converged" and a.physical == b.physical and a.final_cost == b.final_cost
    assert a.reconciliation is not None and b.reconciliation is not None
    assert a.reconciliation.jacobian == b.reconciliation.jacobian and "joint" in a.covariance.assumption
    uncertainty.covariance = []
    refused = reg.fit(p, [coordinate], rows, uncertainty=uncertainty,
        controls=reg.FitControls(maximum_iterations=100, maximum_elapsed_time_seconds=600.0, workers=2))
    assert refused.status == "invalid_setup" and "both" in refused.message
    uncertainty.correlation = None
    refused = reg.fit(p, [coordinate], rows, uncertainty=uncertainty,
        controls=reg.FitControls(maximum_iterations=100, maximum_elapsed_time_seconds=600.0, workers=2))
    assert refused.status == "invalid_setup" and "dimension" in refused.message
