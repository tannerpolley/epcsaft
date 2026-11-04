"""Retained association failure and diagnostic transport cases."""
import copy
from typing import Any, cast

import pytest

from epcsaft import EosError, Parameters, PressureActionDirection, PropertyObservable, Pure
from epcsaft.parameters import document_mapping
from tests.packets import PACKETS


def water_parameters() -> Parameters:
    return Parameters.from_packets(PACKETS, [("held-cameretti-sadowski-2008/4", ["water"])])


def test_association_conditioning_is_readable() -> None:
    model = Pure(water_parameters())
    batch = model.pressure_actions(298.15, 1e5, phase="liquid", observables=[PropertyObservable.ResidualEnthalpy],
                                   directions=[PressureActionDirection(temperature=1.0)])
    rank = batch.diagnostics.site_conditioning
    assert rank is not None and rank.rank == rank.dimension == 2
    assert rank.basis == "implicit_site_fractions"
    assert "rank=2" in repr(rank) and "dimension=2" in repr(rank)


def test_singular_association_names_declared_sites() -> None:
    document = cast(dict[str, Any], copy.deepcopy(document_mapping(water_parameters())))
    document["topology"]["edges"][0]["energy_over_k"]["value"]["magnitude"] = 100000.0
    model = Pure(Parameters.from_mapping(document, components=["water"]))
    with pytest.raises(EosError, match="association Jacobian is singular") as caught:
        model.state(298.15, P=1e5, phase="liquid")
    assert "water/a" in str(caught.value) and "water/b" in str(caught.value)


def test_nonassociating_conditioning_is_absent() -> None:
    params = Parameters.from_packets(PACKETS, [("gross-2001-propane/3", ["propane"])])
    batch = Pure(params).pressure_actions(300.0, 1e5, phase="vapor",
        observables=[PropertyObservable.ResidualEnthalpy], directions=[PressureActionDirection(temperature=1.0)])
    assert batch.diagnostics.site_conditioning is None


def test_asl_thermochemistry_interval_message(capfd: pytest.CaptureFixture[str]) -> None:
    from pathlib import Path
    from pyomo.environ import ExternalFunction  # type: ignore[import-untyped]
    from epcsaft import ThermochemistryPacket, _core
    from epcsaft.eos import property_asl_record

    params = Parameters.from_packets(PACKETS, [("gross-2001-methane-ethane/3", ["methane"])])
    thermal = ThermochemistryPacket(PACKETS / "nasa-glenn-2002-hda-species/1/parameters/nasa9-records.json").select(["methane"])
    record = property_asl_record(Pure(params, thermochemistry=thermal))
    library = str(Path(_core.__file__).with_name("libepcsaft_asl.so"))
    fn = ExternalFunction(library=library, function="epcsaft_pressure_property")
    fn.evaluate_fgh((20000.0, 1e5, 1.0, record, "1", "1"), fgh=2)
    message = capfd.readouterr().err
    assert "OutsideInterval" in message and "thermochemistry interval" in message
    assert "20000" in message and ":1" in message and "site rank 0" not in message


def test_density_structure_message_retains_context() -> None:
    from epcsaft import DensityStructureError

    document = cast(dict[str, Any], copy.deepcopy(document_mapping(water_parameters())))
    coefficients = document["components"][0]["coefficients"]
    energy = next(item for item in coefficients if item["identity"] == "component/water/dispersion_energy_over_k")
    energy["value"]["magnitude"] = 100000.0
    model = Pure(Parameters.from_mapping(document, components=["water"]))
    with pytest.raises(DensityStructureError) as caught:
        model.state(298.15, P=1e5, phase="liquid")
    message = str(caught.value)
    assert "could not resolve the density structure" in message
    assert "enumeration interval=[" in message and "unresolved:" in message
    assert "T=" in message and "P=" in message and "curvature order" in message


def test_diagnostic_reprs_show_values() -> None:
    from epcsaft import _core
    from epcsaft.equilibrium import SolvedStateObservable, SolvedStateObservableKind

    estimate = _core.NumericalEstimate()
    estimate.absolute = 0.125
    assert "absolute=0.125" in repr(estimate) and "available=False" in repr(estimate)
    observable = SolvedStateObservable(SolvedStateObservableKind.Pressure, phase=2, component=1)
    assert "phase=2" in repr(observable) and "component=1" in repr(observable)
