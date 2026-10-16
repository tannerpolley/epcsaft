"""#569: one NASA9 thermochemistry packet (NASA/TP-2002-211556 records, nasa-glenn-2002-hda-species/1), loaded once.
Its native preparation must equal, bit for bit, the records #73's HDA script built (each coefficient times
R_source / R_Engine, intervals [lower, upper) with the last closed), for every species order and repeated selection."""
import copy
import json
from pathlib import Path
from typing import Any

import pytest

from epcsaft import (
    GAS_CONSTANT_J_PER_MOL_K, IdealCorrelation, IdealInterval, IdealNASA9, ParameterError, Parameters, Pure,
    ThermochemistryPacket, ThermochemistryRecord, ThermochemistryStatus,
)
from tests.packets import PACKETS

RECORDS = PACKETS / "nasa-glenn-2002-hda-species/1/parameters/nasa9-records.json"


def former(components: tuple[str, ...]) -> ThermochemistryRecord:
    """The former analyses/2026-hda-two-flash/scripts/hda_inputs.py construction at d0e992c2 (#73)."""
    records = json.loads(RECORDS.read_text())
    scale, p0 = records["gas_constant_J_per_mol_K"] / GAS_CONSTANT_J_PER_MOL_K, records["standard_pressure_Pa"]
    by_id = {item["species_id"]: sorted(item["intervals"], key=lambda i: i["lower_K"]) for item in records["species"]}
    return ThermochemistryRecord(p0, [[IdealInterval(
        item["lower_K"], item["upper_K"], True, index == len(by_id[name]) - 1,
        IdealCorrelation(IdealNASA9([a * scale for a in item["coefficients"]]), p0))
        for index, item in enumerate(by_id[name])] for name in components])


def fields(record: ThermochemistryRecord) -> list[Any]:
    return [record.standard_pressure] + [(i.lower, i.upper, i.lower_inclusive, i.upper_inclusive, i.correlation.source_pressure,
                                          list(i.correlation.form.coefficients)) for s in record.component_intervals for i in s]


def test_prepared_packet_equals_the_former_records_for_every_selection() -> None:
    packet = ThermochemistryPacket(RECORDS)
    species = packet.species_ids
    assert species == ("hydrogen", "methane", "benzene", "toluene")
    assert sum(len(s) for s in packet.select(species).component_intervals) == 9
    for selection in (species, species[::-1], ("toluene", "benzene"), ("methane",), species):
        assert fields(packet.select(selection)) == fields(former(selection))  # == on floats: bitwise here


def test_interval_ends_are_half_open_with_the_last_closed() -> None:
    packet = ThermochemistryPacket(RECORDS)
    parameters = Parameters.from_json(PACKETS / "hda-pc-saft/1/parameters/parameters.json", components=["toluene"])
    low, high = packet.select(["toluene"]).component_intervals[0]

    def cp(temperature: float, intervals: list[IdealInterval]) -> Any:
        state = Pure(parameters, thermochemistry=ThermochemistryRecord(1e5, [intervals])).state(temperature, rho=10.0)
        return state.ideal_isobaric_heat_capacity, state.thermochemistry_status

    closed_low = IdealInterval(low.lower, low.upper, True, True, low.correlation)
    assert cp(1000.0, [low, high]) == cp(1000.0, [high]) != cp(1000.0, [closed_low])
    assert cp(6000.0, [low, high])[1] == ThermochemistryStatus.Available
    assert cp(6000.000001, [low, high])[1] == ThermochemistryStatus.OutsideInterval


@pytest.mark.parametrize("change, cause", [
    (lambda d: d.update(schema_version=2), "form or version"),
    (lambda d: d.update(coefficient_order=d["coefficient_order"][::-1]), "form or version"),
    (lambda d: d["species"].append(copy.deepcopy(d["species"][0])), "duplicate"),
    (lambda d: d.update(gas_constant_J_per_mol_K=0.0), "gas constant"),
    (lambda d: d["species"][0]["intervals"][0]["coefficients"].__setitem__(2, float("nan")), "nonfinite"),
    (lambda d: d["species"][0]["intervals"][1].update(lower_K=1001.0), "overlapping or gapped"),
    (lambda d: d["species"][0].update(intervals=[]), "missing ideal intervals"),
])
def test_malformed_packets_refuse(tmp_path: Path, change: Any, cause: str) -> None:
    document = json.loads(RECORDS.read_text())
    change(document)
    (tmp_path / "records.json").write_text(json.dumps(document))
    with pytest.raises(ParameterError, match=cause):
        ThermochemistryPacket(tmp_path / "records.json")
    with pytest.raises(ParameterError, match="no species ethanol"):
        ThermochemistryPacket(RECORDS).select(["ethanol"])
