"""The ASL property record: fitted parameters write back into their start document, and every thermochemistry form
exports. The callback equals the native pressure actions bitwise (value, T/P gradient and Hessian) at 300 K, 1e6 Pa.
Charged records (Held 2008 Debye-Hueckel, Figiel 2025 Born) at 298.15 K and 1e5 Pa: values equal the direct Engine,
derivatives along charge-neutral directions match central differences (verification only), non-electroneutral weights
are refused, and an IDAES mixer closes its energy balance through Ipopt."""
from __future__ import annotations

import csv
import os
import shutil
import subprocess
from pathlib import Path
from types import SimpleNamespace
from typing import Any, cast

import numpy as np
import pytest
from pyomo.environ import ConcreteModel, Constraint, ExternalFunction, Objective, Var, value  # type: ignore[import-untyped]

from epcsaft import IdealDIPPR107, IdealNASA7, IdealNASA9, IdealPolynomial, IdealShomate, Mixture, Parameters, Pure
from epcsaft import _core, regression
from epcsaft.eos import ParameterError, PressureActionDirection, PropertyObservable, property_asl_record
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


def test_external_ampl_caller_exits_normally(tmp_path: Path) -> None:
    """The external AMPL solver unloads its library before its thread-exit notifications."""
    executable = shutil.which("ipopt") or str(Path.home() / ".idaes/bin/ipopt")
    if not Path(executable).is_file():
        pytest.skip("external AMPL Ipopt executable unavailable")
    model = Pure(Parameters.from_json(parameters_json("gross-2001-propane", 3), components=["propane"]))
    problem = ConcreteModel()
    problem.pressure = Var(initialize=1.0e6, bounds=(0.9e6, 1.1e6))
    problem.fixed_pressure = Constraint(expr=problem.pressure == 1.0e6)
    library = str(Path(_core.__file__).with_name("libepcsaft_asl.so"))
    problem.density = ExternalFunction(library=library,
                                      function="epcsaft_pressure_property")
    problem.objective = Objective(expr=problem.density(300.0, problem.pressure, 1.0,
                                                         property_asl_record(model), "0", "1"))
    problem.write(str(tmp_path / "pressure.nl"), format="nl")
    (tmp_path / "ipopt.opt").write_text("linear_solver mumps\n")
    result = subprocess.run([executable, "pressure.nl", "-AMPL"], cwd=tmp_path, env={**os.environ, "AMPLFUNC": library},
                            capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stdout + result.stderr
    solution = (tmp_path / "pressure.sol").read_text()
    assert "Optimal Solution Found" in solution


CHARGED = {"NaCl": (("held-cameretti-sadowski-2008", 3), BRINE, [0.96, 0.02, 0.02], [0.0, 1.0, 1.0]),
           "MgCl2": (("held-cameretti-sadowski-2008", 3), ("water", "magnesium-cation", "chloride-anion"), [0.97, 0.01, 0.02], [0.0, 1.0, 2.0]),
           "Born": (("figiel-2025-chloride-activity-refit", 2), BRINE, [0.96, 0.02, 0.02], [0.0, 1.0, 1.0])}


def charged(name: str) -> Mixture:
    """A numerical identity thermochemistry record (cp 30 J/(mol K)); no physical ion enthalpy is claimed."""
    form = IdealPolynomial([30.0, 0.0], 300.0, 1000.0, 50.0)
    return Mixture(Parameters.from_json(parameters_json(*CHARGED[name][0]), components=CHARGED[name][1]), thermochemistry=record(form, form, form))


def aqueous_package(model: Mixture) -> Any:
    from idaes.core import FlowsheetBlock  # type: ignore[import-untyped]
    from idaes.models.properties.modular_properties.base.generic_property import GenericParameterBlock  # type: ignore[import-untyped]
    from epcsaft.idaes import property_configuration
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.props = GenericParameterBlock(**property_configuration(model, phases=("liquid",), temperature_ref=298.15, pressure_ref=1.0e5,
        state_bounds={"flow_mol": (0, 1, 10), "temperature": (273.15, 298.15, 373.15), "pressure": (5e4, 1e5, 1e6)}))
    return m


def property_fgh(text: str, q: Any, row: int, fgh: int) -> tuple[float, Any, Any]:
    from epcsaft.idaes import PropertyExternalFunction, asl_library
    callback = PropertyExternalFunction(library=str(asl_library()), function="epcsaft_pressure_property")
    v, g, h = callback.evaluate_fgh([*q, text, str(row), "0"], fgh=fgh)
    n = len(q)
    hessian = np.array([[h[min(i, j) + max(i, j) * (max(i, j) + 1) // 2] for j in range(n)] for i in range(n)]) if fgh == 2 else None
    return float(v), np.array(g[:n]) if fgh else None, hessian


@pytest.mark.parametrize("name", CHARGED)
def test_charged_values_equal_engine(name: str) -> None:
    """F1: every row at every order, and the IDAES aqueous block (record order, not component_list order), equal the
    direct Engine bitwise."""
    from epcsaft.idaes import PressureProperties
    model, w = charged(name), CHARGED[name][2]
    batch = model.pressure_actions(298.15, 1.0e5, x=w, phase="liquid", directions=[PressureActionDirection(0.0, 0.0, [0.0] * 3)],
                                   observables=[PropertyObservable.TotalEnthalpy, PropertyObservable.LogFugacityCoefficient])
    direct = [cast(float, batch.molar_density), *(cast(float, r.value) for r in batch.results)]
    text = property_asl_record(model)
    assert [[property_fgh(text, [298.15, 1.0e5, *w], row, fgh)[0].hex() for fgh in (0, 1, 2)] for row in range(5)] == [[v.hex()] * 3 for v in direct]
    m = aqueous_package(model)
    m.fs.sb = m.fs.props.build_state_block([0], defined_state=True)
    sb = m.fs.sb[0]
    sb.temperature.fix(298.15)
    sb.pressure.fix(1.0e5)
    for j, v in zip(model.component_ids, w, strict=True):
        sb.mole_frac_phase_comp["liquid", j].set_value(v)
    x = {j: sb.mole_frac_phase_comp["liquid", j] for j in model.component_ids}
    got = [value(sb.dens_mol_phase["liquid"]), value(sb.enth_mol_phase["liquid"]),
           *(value(PressureProperties._at(sb, "liquid", 2 + k, sb.temperature, sb.pressure, x, True)) for k in range(3))]
    assert [v.hex() for v in got] == [v.hex() for v in direct]


@pytest.mark.parametrize("name", CHARGED)
def test_charged_derivatives_match_central_differences(name: str) -> None:
    """F2, verification only: |D - A| <= tau |A| + 10 eps (|f+| + |f-|) / 2h along T, P, water and salt; tau 1e-6 for
    the gradient, 1e-5 for the Hessian from differences of the gradient. A mixed (composition, P) entry is differenced
    along its composition direction: a 10 Pa difference of a composition gradient is roundoff-limited."""
    model, w, salt = charged(name), CHARGED[name][2], CHARGED[name][3]
    text, q = property_asl_record(model), np.array([298.15, 1.0e5, *w])
    directions = {"T": (np.eye(5)[0], 1e-3), "P": (np.eye(5)[1], 10.0), "water": (np.eye(5)[2], 1e-5), "salt": (np.array([0, 0, *salt]), 1e-5)}
    eps = float(np.finfo(float).eps)
    for row in range(5):
        _, gradient, hessian = property_fgh(text, q, row, 2)
        for b, (d, h) in directions.items():
            plus, minus = property_fgh(text, q + h * d, row, 1), property_fgh(text, q - h * d, row, 1)
            pairs = [(gradient @ d, plus[0], minus[0], 1e-6)]
            pairs += [(a @ hessian @ d, a @ plus[1], a @ minus[1], 1e-5) for k, (a, _) in directions.items() if b != "P" or k == "P"]
            for exact, fp, fm, tau in pairs:
                assert abs((fp - fm) / (2 * h) - exact) <= tau * abs(exact) + 10 * eps * (abs(fp) + abs(fm)) / (2 * h), (row, b)


def test_non_electroneutral_weights_are_refused() -> None:
    """F3: every order and the endpoint audit refuse a charge-unbalanced composition; charged VLE is refused."""
    from epcsaft.idaes import check_state, property_configuration
    model, bad = charged("NaCl"), [0.96, 0.03, 0.01]
    message = "non-electroneutral composition: charge sum 0.0199"
    for fgh in (0, 1, 2):
        with pytest.raises(RuntimeError, match="epcsaft_property:invalid_input:" + message):
            property_fgh(property_asl_record(model), [298.15, 1.0e5, *bad], 2, fgh)
    m = aqueous_package(model)
    m.fs.sb = m.fs.props.build_state_block([0], defined_state=True)
    sb = m.fs.sb[0]
    sb.temperature.fix(298.15)
    sb.pressure.fix(1.0e5)
    for j, v in zip(model.component_ids, bad, strict=True):
        sb.mole_frac_phase_comp["liquid", j].set_value(v)
    with pytest.raises(ParameterError, match=message):
        check_state(sb, "liquid")
    with pytest.raises(ValueError, match="a charged record admits only"):
        property_configuration(model, phases=("liquid", "vapor"), temperature_ref=298.15, pressure_ref=1.0e5, state_bounds={})


def test_charged_mixer_ipopt_round_trip(tmp_path: Path) -> None:
    """F4: an IDAES Mixer of two MgCl2 brines solves with no adapter refusal at any iterate, an exact outlet
    composition and closed energy balance. The ion fractions start above Ipopt's 0.01 bound push."""
    from idaes.core.solvers import get_solver  # type: ignore[import-untyped]
    from idaes.models.unit_models import Mixer, MomentumMixingType  # type: ignore[import-untyped]
    from epcsaft.idaes import asl_library, solve_with_diagnostics
    solver = get_solver()
    if not solver.available(exception_flag=False):
        pytest.skip("IDAES Ipopt executable unavailable")
    os.environ["AMPLFUNC"] = "\n".join(filter(None, (os.environ.get("AMPLFUNC"), str(asl_library()))))
    model = charged("MgCl2")
    m = aqueous_package(model)
    m.fs.mix = Mixer(property_package=m.fs.props, inlet_list=["in1", "in2"], momentum_mixing_type=MomentumMixingType.equality)
    inlets = (((0.94, 0.02, 0.04), 298.15), ((0.955, 0.015, 0.03), 303.15))
    for port, (x, T) in zip((m.fs.mix.in1, m.fs.mix.in2), inlets, strict=True):
        port.flow_mol.fix(1.0)
        port.temperature.fix(T)
        port.pressure.fix(1.0e5)
        for j, v in zip(model.component_ids, x, strict=True):
            port.mole_frac_comp[0, j].fix(v)
    m.fs.mix.in2.pressure.unfix()  # pressure equality sets it; fixing both inlets over-specifies the mixer
    states = (m.fs.mix.in1_state[0], m.fs.mix.in2_state[0], m.fs.mix.mixed_state[0])
    for state, (x, _) in zip(states, (*inlets, inlets[0]), strict=True):
        for j, v in zip(model.component_ids, x, strict=True):
            state.mole_frac_comp[j].set_value(v)
            state.mole_frac_phase_comp["liquid", j].set_value(v)
    states[2].flow_mol.set_value(2.0)
    log = tmp_path / "ipopt.log"
    solve_with_diagnostics(solver, m, logfile=log, endpoints=((states[2], "liquid"),))
    assert "epcsaft_property failure" not in log.read_text()
    assert [value(states[2].mole_frac_comp[j]) for j in model.component_ids] == pytest.approx([0.9475, 0.0175, 0.035], abs=1e-10)
    flows = []
    for state in states:
        batch = model.pressure_actions(value(state.temperature), value(state.pressure), phase="liquid",
                                       x=[value(state.mole_frac_comp[j]) for j in model.component_ids],
                                       directions=[PressureActionDirection(0.0, 0.0, [0.0] * 3)], observables=[PropertyObservable.TotalEnthalpy])
        flows.append(value(state.flow_mol) * cast(float, batch.results[0].value))
    assert abs(flows[0] + flows[1] - flows[2]) <= 1e-8 * sum(map(abs, flows))
