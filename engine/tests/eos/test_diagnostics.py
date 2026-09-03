from functools import partial
from typing import Callable

import pytest

from epcsaft import Convergence, Mixture, Parameters, Stop, trace

from tests.packets import parameters_json


def test_trace_and_convergence_round_trip() -> None:
    parameters = Parameters.from_json(
        parameters_json("gross-2001-methane-ethane", 2),
        components=("methane", "ethane"),
    )
    mixture = Mixture(parameters)
    with trace() as captured:
        state = mixture.state(
            300.0, P=1.0e5, x=(0.5, 0.5), phase="liquid"
        )
    convergence = state.density_convergence
    assert convergence is not None
    assert isinstance(convergence, Convergence)
    root = next(span for span in captured.spans if span.parent is None)
    assert root.convergence is not None
    assert root.convergence.residual == convergence.residual
    assert root.convergence.iterations == convergence.iterations
    assert root.convergence.stop == convergence.stop
    if convergence.stop == Stop.converged:
        criterion = convergence.criterion
        assert criterion is not None
        assert convergence.residual <= criterion.value
    assert "DENSITY_ROOT_RESIDUAL" in repr(convergence)


def test_reference_requires_stoichiometry() -> None:
    parameters = Parameters.from_json(
        parameters_json("figiel-2025-reference-electrolytes", 3),
        components=("water", "lithium-cation", "chloride-anion"),
    )
    with pytest.raises(TypeError):
        reference: Callable[..., object] = getattr(Mixture(parameters), "reference")
        call = partial(
            reference,
            298.15,
            1.0e5,
            solvent_fractions=(1.0, 0.0, 0.0),
            ions=("lithium-cation", "chloride-anion"),
        )
        call()
