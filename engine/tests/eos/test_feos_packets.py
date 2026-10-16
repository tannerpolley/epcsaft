"""#569: the four FeOs packets (feos-org/feos af84e3f9, data/tools/import_feos.py) against their sources and the
Engine's association conventions. Input verification only; no physical-validation claim."""
import csv
import importlib.util
import json
import math
import sys
from decimal import Decimal
from typing import Any

import pytest

from epcsaft import Mixture, Parameters
from tests.packets import PACKETS

ROOT = PACKETS.parents[1]
SOURCE_RECORDS = {"feos-gross-sadowski-2001": 78, "feos-gross-sadowski-2002": 18, "feos-esper-2023": 1842, "feos-rehner-2023": 7848}


def document(packet: str) -> dict[str, Any]:
    result: dict[str, Any] = json.loads((PACKETS / packet / "1/parameters/parameters.json").read_text())
    return result


def coefficients(component: dict[str, Any]) -> dict[str, float]:
    found = {item["family"]: item["value"]["magnitude"] for item in component["coefficients"]}
    return found | {"molar_mass": component["fixed"]["molar_mass"]["value"]["magnitude"]}


def test_every_source_record_is_admitted_or_refused_once() -> None:
    for packet, count in SOURCE_RECORDS.items():
        rows = list(csv.DictReader((PACKETS / packet / "1/sources/record-accounting.csv").open()))
        assert [int(r["source_index"]) for r in rows] == list(range(count))
        assert all(r["disposition"] == "admitted" or r["cause"] for r in rows)
        if packet != "feos-rehner-2023":
            assert len(document(packet)["components"]) == sum(r["disposition"] == "admitted" for r in rows)
    refused = [r for r in csv.DictReader((PACKETS / "feos-rehner-2023/1/sources/record-accounting.csv").open())
               if r["disposition"] == "refused"]
    assert [r["source_index"] for r in refused] == ["863", "872", "885", "898"]  # negative fitted eps^AiBj


def test_gross_2001_values_equal_the_existing_packet_bitwise() -> None:
    spelling = {"n-propylbenzene": "n-propyl-benzene", "propyl ethanoate": "n-propyl-ethanoate",
                "propyl methanoate": "n-propyl-methanoate", "propyl propanoate": "n-propyl-propanoate"}
    old_packet = json.loads((PACKETS / "gross-sadowski-2001/3/parameters/parameters.json").read_text())
    existing = {c["component_id"]: coefficients(c) for c in old_packet["components"]}
    compared = 0
    for component in document("feos-gross-sadowski-2001")["components"]:
        name = component["component_id"]
        old = existing[spelling.get(name, name.replace(" ", "-").replace(",", "-"))]
        for family, magnitude in coefficients(component).items():
            assert magnitude == old[family], (name, family)  # == on floats: bitwise
            compared += 1
    assert compared == 312


def test_gross_2002_equals_table_1_and_declares_its_rules() -> None:
    table = ROOT / "analyses/reproduction/2002-gross/results/tables/table-001/source/table-001.csv"
    rows = [r for r in csv.reader(table.open()) if len(r) > 6 and r[1]][1:]
    packet = document("feos-gross-sadowski-2002")
    components = {c["component_id"]: c for c in packet["components"]}
    self_edges = {e["endpoint_a"]["component_id"]: e for e in packet["topology"]["edges"] if e["source"]["kind"] == "explicit"}
    compared = 0
    for name, mass, m, sigma, epsilon, kappa, energy, *_ in rows:
        found = coefficients(components[name])
        edge = self_edges[name]
        expected = (float(Decimal(mass) / 1000), float(m), float(sigma), float(epsilon), float(kappa), float(energy))
        actual = (found["molar_mass"], found["segment_count"], found["segment_diameter"], found["dispersion_energy_over_k"],
                  edge["volume"]["value"]["magnitude"], edge["energy_over_k"]["value"]["magnitude"])
        assert actual == expected, name
        compared += len(expected)
    assert (len(rows), compared) == (18, 108)
    rules = [e for e in packet["topology"]["edges"] if e["source"]["kind"] == "combining-rule"]
    assert len(rules) == 18 * 17 and {e["source"]["rule_id"] for e in rules} == {"wolbach-sandler"}


def test_esper_maps_nonpolar_dipolar_self_and_acceptor_only_records() -> None:
    packet = document("feos-esper-2023")
    sites: dict[str, list[tuple[str, int]]] = {}
    for site in packet["topology"]["sites"]:
        sites.setdefault(site["component_id"], []).append((site["site_id"], site["multiplicity"]))
    edges = {e["endpoint_a"]["component_id"] for e in packet["topology"]["edges"]}
    components = {c["component_id"]: coefficients(c) for c in packet["components"]}
    assert (components["methane"]["dipole_moment"], components["methane"]["dipole_count"], sites.get("methane")) == (0, 0, None)
    assert components["acetone"]["dipole_count"] == 1 and components["acetone"]["dipole_moment"] > 0
    assert sites["ethanol"] == [("a", 1), ("b", 1)] and "ethanol" in edges
    assert sites["acetaldehyde"] == [("b", 1)] and "acetaldehyde" not in edges
    assert len(edges) == 387 and sum(set(s) == {("b", n) for _, n in s} for s in sites.values()) == 703


def test_rehner_fitted_cross_volume_keeps_the_source_strength() -> None:
    """Rehner 2023 Eq. 17 builds the cross volume on (sigma_i sigma_j)^(3/2); the Engine's direct edge on
    ((sigma_i + sigma_j)/2)^3. First pinned induced row: acetaldehyde (acceptor only) / 2-propen-1-ol."""
    packet = document("feos-rehner-2023")
    sigma = {c["component_id"]: coefficients(c)["segment_diameter"] for c in packet["components"]}
    edge = next(e for e in packet["topology"]["edges"]
                if {e["endpoint_a"]["component_id"], e["endpoint_b"]["component_id"]} == {"acetaldehyde", "2-propen-1-ol"})
    assert (edge["endpoint_a"]["component_id"], edge["energy_over_k"]["value"]["magnitude"]) == ("2-propen-1-ol", 2319.833028322245)
    s_a, s_b = sigma["acetaldehyde"], sigma["2-propen-1-ol"]
    engine = edge["volume"]["value"]["magnitude"] * ((s_a + s_b) / 2) ** 3
    assert math.isclose(engine, 0.04087 * (s_a * s_b) ** 1.5, rel_tol=1e-10)


def test_converted_direct_edge_reproduces_the_wolbach_sandler_strength() -> None:
    """The same Eq. 17 strength through a direct edge (converted volume, Eq. 18 energy) and through the rule."""
    packet = document("feos-gross-sadowski-2002")
    ids = ("methanol", "water")
    self_edges = {e["endpoint_a"]["component_id"]: e for e in packet["topology"]["edges"] if e["source"]["kind"] == "explicit"}
    sigma = {c["component_id"]: coefficients(c)["segment_diameter"] for c in packet["components"]}
    k = {c: self_edges[c]["volume"]["value"]["magnitude"] for c in ids}
    e = {c: self_edges[c]["energy_over_k"]["value"]["magnitude"] for c in ids}
    volume = math.sqrt(k[ids[0]] * k[ids[1]]) * (sigma[ids[0]] * sigma[ids[1]]) ** 1.5 / ((sigma[ids[0]] + sigma[ids[1]]) / 2) ** 3
    states = []
    for direct in (False, True):
        mapping = dict(packet, pairs=[{"component_id_a": ids[0], "component_id_b": ids[1], "coefficients": [
            {"identity": "k_ij", "family": "k_ij", "value": {"magnitude": 0.0, "unit": "dimensionless"}}]}])
        edges = [x for x in packet["topology"]["edges"] if x["source"]["kind"] == "explicit"
                 or {x["endpoint_a"]["component_id"], x["endpoint_b"]["component_id"]} == set(ids)]
        if direct:
            edges = [x for x in edges if x["source"]["kind"] == "explicit"] + [{
                "endpoint_a": {"component_id": a, "site_id": "a"}, "endpoint_b": {"component_id": b, "site_id": "b"},
                "energy_over_k": {"identity": f"{a}/{b}/e", "value": {"magnitude": (e[a] + e[b]) / 2, "unit": "kelvin"}},
                "volume": {"identity": f"{a}/{b}/k", "value": {"magnitude": volume, "unit": "dimensionless"}},
                "source": {"kind": "explicit"}} for a, b in (ids, ids[::-1])]
        mapping["topology"] = dict(packet["topology"], edges=edges)
        state = Mixture(Parameters.from_mapping(mapping, components=list(ids))).state(330.0, rho=3.0e4, x=(0.4, 0.6))
        states.append([state.pressure, *state.residual_chemical_potential_over_rt])
    assert all(math.isclose(a, b, rel_tol=1e-10, abs_tol=1e-10) for a, b in zip(*states, strict=True))


def test_unsupported_source_terms_refuse() -> None:
    spec = importlib.util.spec_from_file_location("import_feos", ROOT / "data/tools/import_feos.py")
    assert spec is not None and spec.loader is not None
    recipe = importlib.util.module_from_spec(spec)
    sys.path.insert(0, str(ROOT / "data/tools"))
    spec.loader.exec_module(recipe)
    record = {"identifier": {"name": "x", "cas": "0"}, "molarweight": Decimal("44.0"), "m": Decimal("2"),
              "sigma": Decimal("3"), "epsilon_k": Decimal("200"), "mu": Decimal("0")}
    with pytest.raises(recipe.Refused, match="unconsumed field q"):  # a planted quadrupole is never dropped
        recipe.pure(record | {"q": Decimal("4.4")}, "planted", "planted", True)
    with pytest.raises(recipe.Refused, match="association parameters without a donor/acceptor pair"):
        recipe.pure(record | {"association_sites": [{"nb": Decimal(1), "kappa_ab": Decimal("0.01")}]}, "planted", "planted", True)
