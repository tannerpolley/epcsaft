"""Named-phase permittivity at the public regression and fit-worker boundaries."""

import pytest

from epcsaft import Mixture, ParameterError, Parameters, Pure
from epcsaft import regression as reg
from epcsaft.equilibrium import Amounts, Phase, Problem

from tests.native_inputs import ETHANE, parameters


def salt_parameters(alpha: float) -> Parameters:
    return parameters(
        ((18.0153, 1.2047, 2.7927, 353.95), (22.98977, 1.0, 2.8232, 230.0), (35.4527, 1.0, 2.756, 170.0)),
        pairs=(("k_ij", (0, 1), -0.3), ("k_ij", (0, 2), 0.0), ("k_ij", (1, 2), 0.8),
               ("association_energy_over_k", (0, 1), 2425.7), ("association_volume", (0, 1), 0.04509)),
        sites=((0, 0), (1, 1)),
        charges=(0, 1, -1),
        extra=(
            ("packing_diameter", (1,), 2.484416, "angstrom"),
            ("packing_diameter", (2,), 2.42528, "angstrom"),
            ("debye_huckel_diameter", (1,), 2.484416, "angstrom"),
            ("debye_huckel_diameter", (2,), 2.42528, "angstrom"),
            ("dielectric_volume", (0,), 0.3777, "liter / mole"),
            ("dielectric_temperature", (0,), 1403.0, "kelvin"),
            ("ion_suppression_coefficient", (1,), alpha, "dimensionless"),
            ("ion_suppression_coefficient", (2,), alpha, "dimensionless"),
        ),
        choices=("fully-dissociated-debye-huckel", "exclude-same-sign-ion-pairs",
                 "solvent-density-temperature-ion-specific-suppression"),
    )


def salt_problem(pressure: float) -> Problem:
    return Problem(phases=[Phase("L", kind="liquid")], T=298.15, P=pressure,
                   feed=Amounts({"s0": 0.96, "s1": 0.02, "s2": 0.02}))


@pytest.mark.parametrize("kind", ["phase_density", "bulk_relative_permittivity"])
def test_suppression_denominator_refuses_public_prediction_and_fit(kind: str) -> None:
    p = salt_parameters(-30)
    model = Mixture(p)
    row = reg.observation(p, kind, salt_problem(101325), observed=1, form="difference", phase="L",
                          scale=1, mass_basis=False if kind == "phase_density" else None)
    prediction = reg.predict(model, row)
    gradient = reg.differentiate(model, row)
    cause = "invalid permittivity suppression denominator"
    assert prediction.status == gradient.status == reg.ObservationStatus.CentralSolveUnavailable
    assert prediction.value is None and not gradient.value
    assert prediction.message == gradient.message and cause in prediction.message

    coordinate = reg.coordinate(p, [("ion_suppression_coefficient", "s1", 1)], origin=0, scale=1, bounds=(-100, 10))
    result = reg.fit(p, [coordinate], [row], controls=reg.FitControls(
        maximum_iterations=2, maximum_elapsed_time_seconds=60, workers=1))
    assert not result.usable and not result.predictions
    if result.status == "invalid_setup":
        assert result.message == f"observation 0: central_solve_unavailable: {cause}"
        assert result.residual_evaluations == 0
    else:
        assert result.status == "initial_evaluation_failed"
        failure = result.failures[0]
        assert failure.observation == 0 and failure.status == reg.ObservationStatus.CentralSolveUnavailable
        assert cause in failure.message


def test_permittivity_form_is_required_at_admission() -> None:
    p = parameters((ETHANE,))
    problem = Problem(phases=[Phase("L", kind="liquid")], T=180, P=1e6, feed=Amounts({"s0": 1}))
    row = reg.observation(p, "bulk_relative_permittivity", problem, observed=1, form="difference", phase="L", scale=1)
    model = Pure(p)
    message = "bulk relative permittivity requires a declared permittivity form"
    with pytest.raises(ParameterError, match=message):
        reg.predict(model, row)
    with pytest.raises(ParameterError, match=message):
        reg.differentiate(model, row)
    coordinate = reg.coordinate(p, [("dispersion_energy_over_k", "s0", 1)], origin=0, scale=1, bounds=(100, 300))
    result = reg.fit(p, [coordinate], [row], controls=reg.FitControls(
        maximum_iterations=2, maximum_elapsed_time_seconds=60, workers=1))
    assert not result.usable and result.status == "invalid_setup" and message in result.message
    assert result.residual_evaluations == 0


def test_solved_permittivity_on_fit_worker() -> None:
    p = salt_parameters(2)
    first, second = salt_problem(101325), salt_problem(215000)
    rows = [
        reg.observation(p, "phase_density", first, observed=1, form="difference", phase="L", scale=1000, mass_basis=False),
        reg.observation(p, "bulk_relative_permittivity", second, observed=1, form="difference", phase="L", scale=3.2),
        reg.observation(p, "phase_density", second, observed=1, form="difference", phase="L", scale=1000, mass_basis=False),
    ]
    model = Mixture(p)
    for row in rows:
        value = reg.predict(model, row).value
        assert value is not None
        row.observed = value
    coordinate = reg.coordinate(p, [("ion_suppression_coefficient", "s1", 1)], origin=0, scale=2, bounds=(0, 10))
    # The second eligible group and its density row belong to process one.
    serial, parallel = (
        reg.fit(p, [coordinate], rows, controls=reg.FitControls(
            maximum_iterations=10, maximum_elapsed_time_seconds=60, workers=workers))
        for workers in (1, 2)
    )
    assert serial.usable and parallel.usable
    assert len(serial.predictions) == len(serial.optimizer_jacobian) == len(rows)
    assert serial.optimizer_jacobian[1] != 0
    assert serial.predictions == parallel.predictions
    assert serial.optimizer_jacobian == parallel.optimizer_jacobian
