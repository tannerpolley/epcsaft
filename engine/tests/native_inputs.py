"""Explicit inputs of the retained native fixtures, restated through the public package (#59).

`ROWS` holds the C++ values recorded by oracle/record_native_rows.py for these same inputs.
"""
import math
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from epcsaft import (
    GAS_CONSTANT_J_PER_MOL_K as R,
    IdealCorrelation,
    IdealDIPPR107,
    IdealInterval,
    IdealNASA7,
    IdealNASA9,
    IdealPolynomial,
    IdealShomate,
    Parameters,
    ThermochemistryRecord,
    _core,
)

ROWS = {
    line.split("\t")[1]: float(line.split("\t")[3])
    for line in (Path(__file__).parent / "oracle" / "native-rows.tsv").read_text().splitlines()[1:]
}
# The retained same-owner transport rule: max(1e-12, 1e-12 max(|model|, |reference|)).
SAME = {"type": "scalar", "rule": "symmetric_close", "atol": 1e-12, "rtol": 1e-12}


def same(evidence: Any, name: str, model: float | None, reference: float | None = None) -> None:
    reference = ROWS[name] if reference is None else reference
    evidence.compare(name, model=model, reference=reference, criterion=SAME)
    assert model is not None and math.isclose(model, reference, rel_tol=1e-12, abs_tol=1e-12), name


def parameters(
    species: Sequence[tuple[float, float, float, float]],  # g/mol, segments, sigma / angstrom, epsilon / K
    pairs: Sequence[tuple[str, tuple[int, int], float]] = (),  # k_ij and association site-pair entries
    sites: tuple[Sequence[int], Sequence[int]] = ((), ()),
    ions: dict[str, Any] | None = None,
) -> Parameters:
    document = _core.ParameterDocument()
    document.component_ids = [f"s{i}" for i in range(len(species))]
    document.site_component, document.site_multiplicity = list(sites[0]), list(sites[1])
    document.electrolyte_choice = document.permittivity_choice = "none"
    rows: list[tuple[str, tuple[int, ...], float, str]] = list(ions["entries"]) if ions else []
    if ions:
        document.electrolyte_choice, document.ionic_dispersion_choice, document.permittivity_choice = ions["choices"]
    for i, (mass, segments, sigma, epsilon) in enumerate(species):
        rows += [("molar_mass", (i,), mass, "gram / mole"), ("segment_count", (i,), segments, "dimensionless"),
                 ("sigma_constant", (i,), sigma, "angstrom"), ("dispersion_energy_over_k", (i,), epsilon, "kelvin")]
        if not ions:
            rows.append(("charge_number", (i,), 0.0, "elementary-charge"))
    rows += [(field, index, value, "kelvin" if field == "association_energy_over_k" else "dimensionless")
             for field, index, value in pairs]
    entries = []
    for field, index, magnitude, unit in rows:
        item = _core.ParameterEntry()
        item.identity, item.field, item.index = f"{field}/{index}", field, list(index)
        item.magnitude, item.unit = magnitude, unit
        entries.append(item)
    document.entries = entries
    site_ids = tuple((document.component_ids[c], f"site{k}") for k, c in enumerate(sites[0]))
    return Parameters(_core.build_parameters(document), tuple(document.component_ids), site_ids)


ETHANE = (30.07, 1.6069, 3.5206, 191.42)


def associating_binary() -> Parameters:
    """properties_support::associating_binary(), the surface example's amine-water-like pair."""
    edges = (((0, 1), 2383.4744, 0.118488), ((0, 3), 2093.8222, 0.1540467874409226),
             ((1, 2), 2093.8222, 0.1540467874409226), ((2, 3), 1804.17, 0.2039))
    pairs: list[tuple[str, tuple[int, int], float]] = [("k_ij", (0, 1), -0.052)]
    for index, energy, volume in edges:
        pairs += [("association_energy_over_k", index, energy), ("association_volume", index, volume)]
    return parameters(((61.08, 4.5354, 2.6019, 204.0438), (18.01528, 2.1945, 2.229, 141.66)), pairs,
                      ((0, 0, 1, 1), (2, 1, 2, 2)))


def p3_forms() -> dict[str, Any]:
    """properties_support::p3_forms(): cp = 30 J/(mol K), h(300 K) = 1000 J/mol, s(300 K) = 50 J/(mol K)."""
    nasa = (1000.0 / R - 30.0 * 300.0 / R, 50.0 / R - 30.0 * math.log(300.0) / R)
    return {
        "polynomial": IdealPolynomial([30.0, 0.0], 300.0, 1000.0, 50.0),
        "shomate": IdealShomate([30.0, 0.0, 0.0, 0.0, 0.0, -9.0, 50.0 - 30.0 * math.log(0.3), 0.0], 1000.0),
        "nasa7": IdealNASA7([30.0 / R, 0.0, 0.0, 0.0, 0.0, *nasa]),
        "nasa9": IdealNASA9([0.0, 0.0, 30.0 / R, 0.0, 0.0, 0.0, 0.0, *nasa]),
        "dippr107": IdealDIPPR107([0.0, 30.0, 0.0, 0.0, 0.0], 300.0, 1000.0, 50.0),
    }


def record(*forms: Any, lower: float = 250.0, upper: float = 350.0) -> ThermochemistryRecord:
    """One closed interval per component at 1 bar, in component order."""
    return ThermochemistryRecord(1.0e5, [[IdealInterval(lower, upper, True, True, IdealCorrelation(form, 1.0e5))]
                                         for form in forms])
