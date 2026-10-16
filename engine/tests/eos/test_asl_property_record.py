"""The ASL property record: fitted parameters write back into their start document, and every thermochemistry form
exports. The callback equals the native pressure actions bitwise (value, T/P gradient and Hessian) at 300 K, 1e6 Pa."""
from __future__ import annotations

import csv
from pathlib import Path
from types import SimpleNamespace
from typing import Any, cast

import pytest
from pyomo.environ import ExternalFunction  # type: ignore[import-untyped]

from epcsaft import IdealDIPPR107, IdealNASA7, IdealNASA9, IdealPolynomial, IdealShomate, Mixture, Parameters, Pure
from epcsaft import _core, regression
from epcsaft.eos import PressureActionDirection, PropertyObservable, property_asl_record
from epcsaft.equilibrium import Amounts, Free, Phase, Pinned, Problem
from tests.native_inputs import p3_forms, record
from tests.packets import PACKETS, parameters_json
from tests.regression.test_temperature_laws import _document

SELECTIONS = Path(__file__).resolve().parents[3] / "analyses/2026-dielectric-alternatives/inputs/selections"
BRINE = ("water", "sodium-cation", "chloride-anion")


def entries(document: _core.ParameterDocument) -> list[tuple[Any, ...]]:
    return [(e.identity, e.field, list(e.index), e.magnitude.hex(), e.unit) for e in document.entries]


def changed(document: _core.ParameterDocument, values: dict[str, float]) -> _core.ParameterDocument:
    items = document.entries
    for item in items:
        item.magnitude = values.get(item.identity, item.magnitude)
    document.entries = items
    return document


def assert_callback_equals_native(model: Pure) -> None:
    library = str(Path(_core.__file__).parent / "libepcsaft_asl.so")
    callback, text = ExternalFunction(library=library, function="epcsaft_pressure_property"), property_asl_record(model)
    observables = [PropertyObservable.TotalEnthalpy, PropertyObservable.LogFugacityCoefficient]
    t, p = PressureActionDirection(1.0, 0.0, [0.0]), PressureActionDirection(0.0, 1.0, [0.0])
    for phase, name in enumerate(("liquid", "vapor")):
        first, second = ([model.pressure_actions(300.0, 1.0e6, phase=name, observables=observables, directions=d)
                          for d in group] for group in (([t], [p]), ([t, t], [t, p], [p, p])))
        native = [[first[0].molar_density, *(b.density_action[0] for b in first), *(b.density_mixed_action for b in second)]]
        native += [[first[0].results[k].value, *(b.results[k].action for b in first + second)] for k in (0, 1)]
        for row in range(3):
            value, gradient, hessian = callback.evaluate_fgh([300.0, 1.0e6, 1.0, text, str(row), str(phase)], fgh=2)
            assert [float(v).hex() for v in (value, *gradient[:2], *hessian[:3])] == [None if v is None else v.hex() for v in native[row]]


def test_fitted_propane_exports() -> None:
    """F1: the P10 fit (Gross 2001 start, NIST SRD 69 propane at 150-360 K) writes back and exports bitwise."""
    start = Parameters.from_json(parameters_json("gross-2001-propane", 3), components=["propane"])
    rows = [r for r in csv.DictReader((PACKETS / "nist-srd69-pure-saturation-properties/1/observations"
                                       / "propane-reference-saturation.csv").open()) if 150.0 <= float(r["T_K"]) <= 360.0]
    observations = []
    for r in rows:
        saturation = Problem(phases=[Phase("L", kind="liquid"), Phase("V", kind="vapor", amount=Pinned(0.0))],
                             T=float(r["T_K"]), P=Free(float(r["p_sat_Pa"])), feed=Amounts({"propane": 1.0}))
        observations += [regression.observation(start, "pressure", saturation, observed=float(r["p_sat_Pa"]), form="log_ratio"),
                         regression.observation(start, "phase_density", saturation, observed=float(r["rho_sat_liq_kg_m3"]),
                                                form="log_ratio", mass_basis=True, phase="L")]
    coordinates = [regression.coordinate(start, [(name, "propane", 1.0)], origin=0.0, scale=scale, bounds=bounds)
                   for name, scale, bounds in (("segment_count", 1.0, (1.0, 5.0)), ("segment_diameter", 1.0, (2.5, 5.0)),
                                               ("dispersion_energy_over_k", 100.0, (100.0, 400.0)))]
    result = regression.fit(start, coordinates, observations,
                            controls=regression.FitControls(maximum_iterations=100, maximum_elapsed_time_seconds=600.0))
    assert result.usable
    fitted = regression.fitted_parameters(start, result)
    assert isinstance(fitted._document, _core.ParameterDocument) and isinstance(start._document, _core.ParameterDocument)
    fields = ("segment_count", "sigma_constant", "dispersion_energy_over_k")
    written = {e.field: e.magnitude.hex() for e in fitted._document.entries if e.field in fields}
    assert written == {field: value.hex() for field, value in zip(fields, result.physical, strict=True)}
    assert [e for e in entries(fitted._document) if e[1] not in fields] == [e for e in entries(start._document) if e[1] not in fields]
    assert_callback_equals_native(Pure(fitted, thermochemistry=record(p3_forms()["nasa9"], upper=450.0)))


def test_write_back_across_families() -> None:
    """F1b: water-first Born selections and a temperature-law document round-trip, changed values rebuild bitwise, and
    a slot with no entry is refused."""
    def load(name: str) -> _core.ParameterDocument:
        parameters = Parameters.from_mapping(_document()) if name == "law" else Parameters.from_json(SELECTIONS / f"{name}.json", components=BRINE)
        return cast(_core.ParameterDocument, parameters._document)

    suppression = {"pair/sodium-cation/water/k_ij": -0.25, "association/water/a/water/b/energy_over_k": 2400.3,
                   "model/relative_permittivity/ionic_region_relative_permittivity": 9.5,
                   "model/relative_permittivity/ion_fraction_suppression_coefficient": 6.5}
    # 0.4884 L/mol does not survive a unit round trip, so an unfitted entry must not be rewritten; a fitted 0.302 L/mol
    # rebuilds only through value / 1e-3, not value * 1e3.
    volume = "component/water/dielectric_volume"
    # A law term's index is its law and power, not a slot; its fitted a0 and a1 must reach their own entries.
    law = {"manufactured/k01/a0": -0.061, "manufactured/k01/a1": 1.7e-4, "manufactured/k01/s": 12.5}
    for name, base, values in (("ion-fraction-suppression", {}, suppression), ("density-temperature", {volume: 0.4884}, {volume: 0.302}),
                               ("law", {}, law)):
        document, target = changed(load(name), base), changed(load(name), values)
        assert entries(_core.fitted_document(document, _core.build_parameters(document))) == entries(document)
        written = _core.fitted_document(document, _core.build_parameters(target))
        assert entries(_core.fitted_document(target, _core.build_parameters(written))) == entries(target)
        assert [e for e in entries(written) if e[1] != "dielectric_volume"] == [e for e in entries(target) if e[1] != "dielectric_volume"]

    # A coefficient may carry its law's identity; only the coefficient entry is written back, never the law declaration.
    def collided(a1: float) -> _core.ParameterDocument:
        document = load("law")
        items = document.entries
        for item in items:
            if item.identity == "manufactured/k01/a1":
                item.identity, item.magnitude = "manufactured/k01", a1
        document.entries = items
        return document
    target = collided(1.7e-4)
    assert entries(_core.fitted_document(collided(1e-4), _core.build_parameters(target))) == entries(target)

    target = changed(load("ion-fraction-suppression"), suppression)
    lij = _core.ParameterEntry()
    lij.identity, lij.field, lij.index, lij.magnitude, lij.unit = "pair/sodium-cation/water/l_ij", "l_ij", [1, 0], 0.01, "dimensionless"
    target.entries = [*target.entries, lij]
    refusal = "fitted l_ij slot 1 has no document entry"
    with pytest.raises(ValueError, match=refusal):
        _core.fitted_document(load("ion-fraction-suppression"), _core.build_parameters(target))
    start = Parameters.from_json(SELECTIONS / "ion-fraction-suppression.json", components=BRINE)
    fit = cast(regression.FitResult, SimpleNamespace(parameters=_core.build_parameters(target)))
    model = Mixture(regression.fitted_parameters(start, fit))
    assert model.state(298.15, P=1.0e5, x={"water": 0.98, "sodium-cation": 0.01, "chloride-anion": 0.01}, phase="liquid").pressure > 0.0
    with pytest.raises(ValueError, match="ASL record: the start document cannot represent the fit: " + refusal):
        property_asl_record(model)


@pytest.mark.parametrize("form", [
    IdealPolynomial([70.1, 0.113, -2.31e-4, 1.72e-7], 298.0, 1234.5, 271.3),
    IdealShomate([-4.22, 306.3, -158.6, 32.15, 0.0875, -114.2, 169.3, -104.7], -104680.0),
    IdealNASA7([0.9335, 0.02642, 6.106e-6, -2.198e-8, 9.515e-12, -13958.5, 19.2017]),
    IdealNASA9([-243314.4, 4656.27, -29.3947, 0.118895, -1.37631e-4, 8.81482e-8, -2.34299e-11, -35403.35, 184.175]),
    IdealDIPPR107([59.192, 176.56, 1626.5, 95.23, 726.2], 300.5, 2345.6, 270.1),
], ids=["polynomial", "shomate", "nasa7", "nasa9", "dippr107"])
def test_every_thermochemistry_form_exports(form: Any) -> None:
    """F2: every coefficient and anchor nonzero and distinct, so a dropped or permuted slot changes the callback."""
    model = Pure(Parameters.from_json(parameters_json("gross-2001-propane", 3), components=["propane"]),
                 thermochemistry=record(form, upper=450.0))
    assert all((model.state(T, P=1.0e3, phase="vapor").ideal_isobaric_heat_capacity or 0.0) > 0.0 for T in (250.0, 350.0, 450.0))
    assert_callback_equals_native(model)
