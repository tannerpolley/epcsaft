import math
from collections.abc import Callable, Sequence
from typing import Any, cast

import pytest

from epcsaft import (
    DensityActionDirection,
    Mixture,
    ParameterError,
    Parameters,
    PressureActionDirection,
    PropertyActionStatus,
    PropertyObservable,
    Pure,
)
from tests.native_inputs import mixed_difference, packet_metadata
from tests.packets import parameters_json


def test_public_residual_actions_match_independent_stencils(evidence: Any) -> None:
    path = parameters_json("gross-2001-methane-ethane", 3)
    components = ("methane", "ethane")
    model = Mixture(Parameters.from_json(path, components=components))
    temperature = 300.0
    density = 5000.0
    composition = (0.4, 0.6)
    central = model.state(temperature, rho=density, x=composition)
    observable = PropertyObservable.ResidualEnthalpy
    raw_density_directions = (
        (1.0, 0.0, (0.0, 0.0)),
        (0.0, 20.0, (0.0, 0.0)),
        (0.0, 0.0, (0.02, -0.02)),
    )
    raw_pressure_directions = (
        (1.0, 0.0, (0.0, 0.0)),
        (0.0, 1000.0, (0.0, 0.0)),
        (0.0, 0.0, (0.02, -0.02)),
    )
    density_directions = tuple(DensityActionDirection(*item) for item in raw_density_directions)
    pressure_directions = tuple(PressureActionDirection(*item) for item in raw_pressure_directions)
    metadata = packet_metadata(path, components)
    criterion = {
        "type": "scalar",
        "rule": "symmetric_close",
        "atol": 2e-8,
        "rtol": 2e-5,
        "scale": "max(abs(model), abs(reference))",
        "combining_rule": "max(atol, rtol*scale)",
    }

    for chart, coordinate, raw, directions in (
        ("density", density, raw_density_directions, density_directions),
        ("pressure", central.pressure, raw_pressure_directions, pressure_directions),
    ):
        for order in range(1, 4):
            evaluate: Callable[[float, float, tuple[float, ...]], float]
            if chart == "density":
                result = model.density_actions(
                    temperature, density, x=composition,
                    observables=(observable,),
                    directions=cast(
                        Sequence[DensityActionDirection], directions[:order]
                    ),
                ).results[0]

                def evaluate_density(
                    t: float, rho: float, x: tuple[float, ...]
                ) -> float:
                    return model.state(t, rho=rho, x=x).residual_enthalpy

                evaluate = evaluate_density
            else:
                result = model.pressure_actions(
                    temperature, central.pressure, x=composition, phase="vapor",
                    anchor=central, observables=(observable,),
                    directions=cast(
                        Sequence[PressureActionDirection], directions[:order]
                    ),
                ).results[0]

                def evaluate_pressure(
                    t: float, pressure: float, x: tuple[float, ...]
                ) -> float:
                    return model.state(
                        t, P=pressure, x=x, phase="vapor", anchor=central,
                    ).residual_enthalpy

                evaluate = evaluate_pressure

            assert result.status == PropertyActionStatus.Available
            assert result.value is not None and result.action is not None
            assert math.isclose(result.value, central.residual_enthalpy, rel_tol=1e-12)
            steps = (0.2, 0.1) if chart == "density" else (0.1, 0.05)
            for step in steps:
                reference = mixed_difference(
                    evaluate, temperature, coordinate, composition, raw[:order], step,
                )
                evidence.compare(
                    f"public_actions/{chart}/order_{order}/step_{step}",
                    quantity="residual enthalpy directional action", model=result.action,
                    reference=reference, unit="J/mol", criterion=criterion,
                    inputs={
                        "T_K": temperature, "chart_coordinate": coordinate,
                        "composition": list(composition), "directions": raw[:order],
                        "dimensionless_step": step, "_parameters": metadata,
                    },
                    basis="Relationship", method="independent centered mixed stencil",
                    claim="The public action equals the native-owned directional contraction.",
                    coordinates={"chart": chart, "order": order},
                    quantity_latex=r"D_1\cdots D_k h^{res}",
                    inputs_latex=r"T,(\rho\;\mathrm{or}\;P),\mathbf{x}",
                )
                assert math.isclose(result.action, reference, rel_tol=2e-5, abs_tol=2e-8)


def test_public_action_availability_and_pure_identity() -> None:
    path = parameters_json("gross-2001-methane-ethane", 3)
    observable = PropertyObservable.ResidualEnthalpy
    mixture = Mixture(Parameters.from_json(path, components=("methane", "ethane")))
    def inserted(eps: float) -> Any:
        return mixture.density_actions(
            300.0, 5000.0, x=(1.0 - eps, eps), observables=(observable,),
            directions=(DensityActionDirection(0.0, 0.0, (-0.01, 0.01)),),
        ).results[0]

    # A neutral insertion is its exact infinite-dilution limit (#73): f(0) = 2 f(eps) - f(2 eps) + O(eps^2).
    support = inserted(0.0)
    assert support.status == PropertyActionStatus.Available
    assert support.value is not None and math.isfinite(support.value)
    assert math.isclose(support.action, 2 * inserted(1e-7).action - inserted(2e-7).action, rel_tol=1e-9)

    partial = mixture.density_actions(
        300.0, 5000.0, x=(0.4, 0.6),
        observables=(PropertyObservable.IdealEnthalpy, observable),
        directions=(DensityActionDirection(1.0, 0.0, (0.0, 0.0)),),
    ).results
    assert partial[0].value is None and partial[0].action is None
    assert partial[0].status == PropertyActionStatus.MissingThermochemistry
    assert partial[1].status == PropertyActionStatus.Available
    assert partial[1].value is not None and partial[1].action is not None
    assert math.isfinite(partial[1].value) and math.isfinite(partial[1].action)

    pure = Pure(Parameters.from_json(path, components=("ethane",))).density_actions(
        300.0, 5000.0, observables=(observable,),
        directions=(DensityActionDirection(1.0),),
    ).results[0]
    endpoint = mixture.density_actions(
        300.0, 5000.0, x=(0.0, 1.0), observables=(observable,),
        directions=(DensityActionDirection(1.0, 0.0, (0.0, 0.0)),),
    ).results[0]
    assert pure.value is not None and pure.action is not None
    assert endpoint.value is not None and endpoint.action is not None
    assert math.isclose(pure.value, endpoint.value, rel_tol=1e-12)
    assert math.isclose(pure.action, endpoint.action, rel_tol=1e-12)


def test_public_action_rejects_invalid_directions() -> None:
    path = parameters_json("gross-2001-methane-ethane", 3)
    model = Mixture(Parameters.from_json(path, components=("methane", "ethane")))
    observable = PropertyObservable.ResidualEnthalpy
    valid = DensityActionDirection(1.0, 0.0, (0.0, 0.0))
    invalid = (
        (),
        (valid,) * 4,
        (DensityActionDirection(0.0, 0.0, (0.0,)),),
        (DensityActionDirection(math.inf, 0.0, (0.0, 0.0)),),
        (DensityActionDirection(0.0, 0.0, (0.1, 0.0)),),
    )
    for directions in invalid:
        with pytest.raises(ParameterError):
            model.density_actions(
                300.0, 5000.0, x=(0.4, 0.6),
                observables=(observable,), directions=directions,
            )
    for observables in ((), (observable, observable)):
        with pytest.raises(ParameterError):
            model.density_actions(
                300.0, 5000.0, x=(0.4, 0.6),
                observables=observables, directions=(valid,),
            )

    ionic = Mixture(Parameters.from_json(
        parameters_json("figiel-2025-reference-electrolytes", 3),
        components=("water", "lithium-cation", "chloride-anion"),
    ))
    with pytest.raises(ParameterError, match="charge balance"):
        ionic.density_actions(
            298.15, 55000.0, x=(0.98, 0.01, 0.01), observables=(observable,),
            directions=(DensityActionDirection(0.0, 0.0, (0.0, 0.001, -0.001)),),
        )
