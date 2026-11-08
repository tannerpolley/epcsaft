"""#1103: a vector-of-class field returns copies; a list read earlier once dangled after the setter replaced the vector."""
from collections.abc import Callable
from typing import Any

import pytest

from epcsaft import _core, ThermochemistryPacket
from tests.packets import PACKETS

RECORDS = PACKETS / "nasa-glenn-2002-hda-species/1/parameters/nasa9-records.json"


def declared(owner: type, field: str, element: type, attribute: str, value: Any) -> Callable[[], tuple[Any, str, Callable[[Any], Any]]]:
    def make() -> tuple[Any, str, Callable[[Any], Any]]:
        item, holder = element(), owner()
        setattr(item, attribute, value)
        setattr(holder, field, [item])
        return holder, field, lambda items: getattr(items[0], attribute)
    return make


def thermochemistry() -> tuple[Any, str, Callable[[Any], Any]]:
    return ThermochemistryPacket(RECORDS).select(["toluene"]), "component_intervals", lambda items: items[0][0].lower


CASES = {
    "EquilibriumProblem.phases": declared(_core.EquilibriumProblem, "phases", _core.PhaseSpec, "name", "L"),
    "EquilibriumProblem.reactions": declared(_core.EquilibriumProblem, "reactions", _core.ReactionSpec, "name", "r"),
    "EquilibriumProblem.solids": declared(_core.EquilibriumProblem, "solids", _core.SolidSpec, "name", "s"),
    "StabilityOptions.trials": declared(_core.StabilityOptions, "trials", _core.PhaseSpec, "name", "t"),
    "Coordinate.targets": declared(_core.Coordinate, "targets", _core.CoordinateTarget, "coefficient", 2.5),
    "CoordinateSet.coordinates": declared(_core.CoordinateSet, "coordinates", _core.Coordinate, "scale", 2.5),
    "MeasurementGroup.inputs": declared(_core.MeasurementGroup, "inputs", _core.MeasuredInput, "sigma", 2.5),
    "MeasurementUncertainty.groups": declared(_core.MeasurementUncertainty, "groups", _core.MeasurementGroup, "feed_total_mol", 2.5),
    "ThermochemistryRecord.component_intervals": thermochemistry,
}


@pytest.mark.parametrize("make", CASES.values(), ids=CASES)
def test_elements_read_before_reassignment_stay_valid(make: Callable[[], tuple[Any, str, Callable[[Any], Any]]]) -> None:
    holder, field, read = make()
    earlier = getattr(holder, field)
    expected = read(earlier)
    setattr(holder, field, getattr(holder, field))
    setattr(holder, field, earlier)
    assert read(earlier) == expected == read(getattr(holder, field))
