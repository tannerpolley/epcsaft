import gc
import sys

import pytest

from epcsaft import Mixture, ParameterError, Parameters, _core
from epcsaft.equilibrium import (
    Amounts,
    Molalities,
    Phase,
    Problem,
    Reaction,
    problem_from_mapping,
)

from tests.packets import parameters_json


def mixture() -> Mixture:
    return Mixture(Parameters.from_json(
        parameters_json("gross-2001-methane-ethane", 2),
        components=("methane", "ethane"),
    ))


def problem(feed: Amounts) -> Problem:
    return Problem(phases=[Phase("L")], T=300.0, P=1.0e5, feed=feed)


def test_compiled_system_keeps_model_alive() -> None:
    model = mixture()
    native_model = model._model
    references = sys.getrefcount(native_model)
    compiled = _core.compile_problem(native_model, problem(Amounts({"methane": 1.0}))._native(model))
    assert sys.getrefcount(native_model) == references + 1
    del compiled
    gc.collect()
    assert sys.getrefcount(native_model) == references


def test_unknown_species_and_amount_status_are_rejected_by_name() -> None:
    model = mixture()
    with pytest.raises(ParameterError, match="unknown species methnae"):
        problem(Amounts({"methnae": 1.0}))._native(model)
    with pytest.raises(ParameterError, match="unknown amount status fixed"):
        problem_from_mapping({
            "phases": [{"name": "L", "amount": {"status": "fixed"}}],
            "T": 300.0,
            "P": 1.0e5,
            "feed": {"basis": "amounts", "values": {"methane": 1.0}},
        })


def test_physical_arguments_are_required() -> None:
    with pytest.raises(TypeError, match="equilibrium_constant"):
        Reaction({"methane": -1.0})  # type: ignore[call-arg]
    with pytest.raises(TypeError, match="solvent_mass_kg"):
        Molalities({"methane": 1.0}, solvent=("methane",))  # type: ignore[call-arg]
