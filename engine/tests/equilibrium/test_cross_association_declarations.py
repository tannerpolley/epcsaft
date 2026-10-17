"""#569: explicit cross-association admission (full EOS §5). P11 is ethanol/water at 351 K, 2B/2B, k_ij = -0.03
(Gross-Sadowski 2002 Table 1 pure values), whose bubble pressures the Engine at 517a3450 retained in the 2026-10-15
audit probe (build/audit-2026-10-15/opus/probes/p11): without cross edges it found seven bubble points and refused
x_ethanol 0.3, 0.5 and 0.7 as locally unstable; with Wolbach-Sandler cross edges it found all ten."""
import ctypes
import math
from pathlib import Path
from typing import Any

import pytest

from epcsaft import (
    IdealCorrelation, IdealInterval, IdealNASA9, Mixture, ParameterError, Parameters, ThermochemistryRecord, _core,
)
from epcsaft.eos import property_asl_record
from epcsaft.equilibrium import MoleFractions, Phase, bubble_point, solve_equilibrium

# x_ethanol: (P / Pa, y_ethanol) without cross edges (None: refused) and with Wolbach-Sandler edges.
RETAINED = {
    0.3: (None, (83271.34990375616, 0.5413045652524671)),
    0.5: (None, (89525.09777215665, 0.634311475445005)),
    0.7: (None, (94561.68351836465, 0.7569685835733344)),
    0.8: ((163110.48974928935, 0.5814939139688228), (96434.64276910479, 0.8298302556612382)),
    0.85: ((161904.58289442613, 0.5888270983822638), (97192.18546472151, 0.8692382577119703)),
    0.88: ((160226.73267754345, 0.5979597675935466), (97588.0978176885, 0.8938623260594497)),
    0.9: ((158425.19709870123, 0.6072633895056476), (97827.4417765058, 0.9106939508769847)),
    0.92: ((155763.58897099853, 0.6206263289099221), (98047.04813682604, 0.9278627356823218)),
    0.95: ((148882.86971319484, 0.6549514169377036), (98339.39168574553, 0.9542576236325229)),
    0.98: ((133098.2773987926, 0.7401918871661524), (98587.24270050864, 0.9814366438314484)),
}
PURE = {"ethanol": (0.04607, 2.3827, 3.1771, 198.24, 2653.4, 0.032384),  # kg/mol, m, sigma/A, eps/K, eps_AB/K, kappa
        "water": (0.018015, 1.0656, 3.0007, 366.51, 2500.7, 0.034868)}
NONE = {"component_id_a": "ethanol", "component_id_b": "water", "provenance": {"locator": "P11: no cross edge"}}


def value(magnitude: float, unit: str) -> dict[str, Any]:
    return {"magnitude": magnitude, "unit": unit}


def document(none: list[dict[str, Any]] | None = None, rule: bool = False) -> dict[str, Any]:
    components, sites, edges = [], [], []
    for c, (mass, m, sigma, eps, energy, volume) in PURE.items():
        components.append({"component_id": c, "fixed": {"molar_mass": {"value": value(mass, "kilogram / mole")},
                                                         "charge_number": {"value": value(0, "elementary-charge")}},
                           "coefficients": [{"identity": f"{c}/{f}", "family": f, "value": value(v, u)} for f, v, u in (
                               ("segment_count", m, "dimensionless"), ("segment_diameter", sigma, "angstrom"),
                               ("dispersion_energy_over_k", eps, "kelvin"))]})
        sites += [{"component_id": c, "site_id": s, "multiplicity": 1} for s in "ab"]
        edges.append({"endpoint_a": {"component_id": c, "site_id": "a"}, "endpoint_b": {"component_id": c, "site_id": "b"},
                      "energy_over_k": {"identity": f"{c}/energy", "value": value(energy, "kelvin")},
                      "volume": {"identity": f"{c}/volume", "value": value(volume, "dimensionless")},
                      "source": {"kind": "explicit"}})
    for a, b in (("ethanol", "water"), ("water", "ethanol")) if rule else ():
        edges.append({"endpoint_a": {"component_id": a, "site_id": "a"}, "endpoint_b": {"component_id": b, "site_id": "b"},
                      "source": {"kind": "combining-rule", "rule_id": "wolbach-sandler", "inputs": {
                          "energy_over_k": [f"{a}/energy", f"{b}/energy"], "volume": [f"{a}/volume", f"{b}/volume"]}}})
    families = [{"kind": "base", "choice": "pc-saft"}, {"kind": "association", "choice": "general-site"},
                {"kind": "electrolyte", "choice": "none"}, {"kind": "permittivity", "choice": "none"}]
    pairs = [{"component_id_a": "ethanol", "component_id_b": "water",
              "coefficients": [{"identity": "k_ij", "family": "k_ij", "value": value(-0.03, "dimensionless")}]}]
    return {"schema": "epcsaft.parameters", "schema_version": 1, "components": components, "pairs": pairs,
            "model_families": families, "topology": {"sites": sites, "edges": edges, "no_cross_association": none or []}}


def close(actual: float, expected: float) -> bool:
    return abs(actual - expected) <= 1e-10 * max(abs(expected), 1.0)


@pytest.mark.parametrize("rule", [False, True])
def test_p11_declared_cross_association_preserves_retained_bubble_points(rule: bool) -> None:
    mixture = Mixture(Parameters.from_mapping(document(None if rule else [NONE], rule)))
    for x, retained in RETAINED.items():
        expected = retained[rule]
        result = solve_equilibrium(mixture, bubble_point(liquid=Phase("L", kind="liquid"), T=351.0,
                                                         feed=MoleFractions({"ethanol": x, "water": 1.0 - x})))
        if expected is None:
            assert not result.success and result.message.startswith("local phase instability"), (x, result.message)
        else:
            assert result.success and result.pressure is not None and close(result.pressure, expected[0]), (x, result.pressure)
            assert close(result.mole_fractions[2], expected[1]), (x, result.mole_fractions)


def test_declarations_refuse_before_evaluation_and_ignore_other_selections() -> None:
    with pytest.raises(ParameterError, match="missing cross association or no_cross_association: ethanol/water"):
        Parameters.from_mapping(document())
    with pytest.raises(ParameterError, match="no_cross_association contradicts a cross edge: ethanol/water"):
        Parameters.from_mapping(document([NONE], rule=True))
    reversed_none = {**NONE, "component_id_a": "water", "component_id_b": "ethanol"}
    for bad in ([NONE, reversed_none], [{**NONE, "component_id_b": "ethanol"}], [{**NONE, "component_id_b": "methanol"}],
                [{"component_id_a": "ethanol", "component_id_b": "water"}]):
        with pytest.raises(ParameterError, match=r"no_cross_association\[\d\]"):
            Parameters.from_mapping(document(bad), components=["water"])
    Parameters.from_mapping(document([NONE]), components=["water"])  # declaration outside the selection
    forward, backward = (Mixture(Parameters.from_mapping(document([n]), components=order)).state(
        351.0, rho=15000.0, x=(0.8, 0.2) if order[0] == "ethanol" else (0.2, 0.8))
        for n, order in ((NONE, ["ethanol", "water"]), (reversed_none, ["water", "ethanol"])))
    assert close(forward.pressure, backward.pressure)
    # The native owner refuses a direct document's duplicate, self-pair and nonzero declarations.
    native = Parameters.from_mapping(document([NONE]))._document
    assert isinstance(native, _core.ParameterDocument)
    for index, magnitude, cause in (((1, 0), 0.0, "duplicate"), ((0, 0), 0.0, "invalid component pair"),
                                    ((0, 1), 1.0, "nonzero declaration")):
        entry = next(e for e in native.entries if e.field == "no_cross_association")
        changed = _core.ParameterEntry()
        changed.identity, changed.field, changed.index, changed.magnitude, changed.unit = (
            "probe", entry.field, list(index), magnitude, entry.unit)
        trial = _core.ParameterDocument()
        for name in ("component_ids", "site_component", "site_multiplicity", "base_choice", "association_choice",
                     "electrolyte_choice", "ionic_dispersion_choice", "permittivity_choice", "polar_choice"):
            setattr(trial, name, getattr(native, name))
        trial.entries = [e for e in native.entries if e.field != "no_cross_association" or cause == "duplicate"] + [changed]
        with pytest.raises(ParameterError, match=cause):
            _core.build_parameters(trial)


def test_asl_record_transports_the_declaration_to_the_same_native_owner() -> None:
    """The epcsaft_pressure_property record carries the declaration entry; its liquid density equals the Python
    state's, and the same record without the entry is refused by native admission inside the ASL library."""
    from pyomo.core.base.external import _ARGLIST  # type: ignore[import-untyped]
    from pyomo.environ import ExternalFunction  # type: ignore[import-untyped]

    interval = IdealInterval(200.0, 1000.0, True, True, IdealCorrelation(IdealNASA9([0.0, 0.0, 4.0] + [0.0] * 6), 1e5))
    mixture = Mixture(Parameters.from_mapping(document([NONE])), thermochemistry=ThermochemistryRecord(1e5, [[interval]] * 2))
    external = ExternalFunction(library=str(Path(_core.__file__).parent / "libepcsaft_asl.so"), function="epcsaft_pressure_property")
    external.load_library()
    function = external._known_functions["epcsaft_pressure_property"][0]

    def density(record: str) -> tuple[float, str]:
        packed = _ARGLIST([351.0, 1e5, 0.8, 0.2, record, "0", "0"], fgh=0)
        return function(ctypes.byref(packed)), packed.Errmsg.decode() if packed.Errmsg else ""

    transported, message = density(property_asl_record(mixture))
    assert message == "" and close(transported, mixture.state(351.0, P=1e5, x=(0.8, 0.2), phase="liquid").molar_density)
    native = mixture.parameters._document
    assert isinstance(native, _core.ParameterDocument)
    native.entries = [e for e in native.entries if e.field != "no_cross_association"]
    refused, message = density(property_asl_record(mixture))
    assert not math.isfinite(refused) and "missing cross association or no_cross_association: ethanol/water" in message
