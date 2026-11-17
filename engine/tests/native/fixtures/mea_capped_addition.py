"""Regenerate mea_capped_addition.inc: python mea_capped_addition.py <_core extension> <output>.

The committed .inc is checked against this output by tools/check-affected.sh --python; no thermodynamics here."""
import copy
import importlib.util
import sys
import types
import gzip
import hashlib
import json
from pathlib import Path
from typing import Any

root = Path(__file__).resolve().parents[4]
package = types.ModuleType("epcsaft")
package.__dict__["__path__"] = [str(root / "engine/src/epcsaft")]
sys.modules["epcsaft"] = package
spec = importlib.util.spec_from_file_location("epcsaft._core", sys.argv[1])
assert spec is not None and spec.loader is not None
core = importlib.util.module_from_spec(spec)
sys.modules["epcsaft._core"] = core
spec.loader.exec_module(core)
from epcsaft.parameters import Parameters  # noqa: E402
raw = gzip.decompress((root / "tools/fixtures/mea-f1.json.gz").read_bytes())
identity = json.loads((root / "tools/fixtures/mea-f1-identity.json").read_text())
assert hashlib.sha256(raw).hexdigest() == identity["packet_sha256"]
packet = json.loads(raw)
packet["parameters"]["topology"]["no_cross_association"] = [{
    "component_id_a": "carbon-dioxide", "component_id_b": "monoethanolamine",
    "provenance": {"source_id": "mea-f1", "locator": "tools/fixtures/mea-f1.json.gz: no CO2/MEA edge"}}]
d = Parameters.from_mapping(packet["parameters"])._document
assert d is not None and not isinstance(d, str)
out = ["using namespace epcsaft; using namespace epcsaft::equilibrium; using namespace epcsaft::regression;", "ParameterDocument document;"]


def literal(value: Any) -> str:
    if isinstance(value, str):
        return json.dumps(value)
    if isinstance(value, bool):
        return str(value).lower()
    if isinstance(value, list):
        return "{" + ",".join(literal(v) for v in value) + "}"
    return repr(value)


for field in ("component_ids", "site_component", "site_multiplicity", "base_choice", "association_choice",
              "electrolyte_choice", "ionic_dispersion_choice", "permittivity_choice", "polar_choice"):
    out.append(f"document.{field}={literal(getattr(d, field))};")
for entry in d.entries:
    values = [entry.identity, entry.field, list(entry.index), entry.magnitude, entry.unit]
    out.append(f"document.entries.push_back({literal(values)});")
counter = 0


def declaration(value: Any) -> str:
    global counter
    if value is None:
        return "std::nullopt"
    if isinstance(value, list):
        return "{" + ",".join(declaration(v) for v in value) + "}"
    if not isinstance(value, dict):
        return literal(value)
    if "enum" in value:
        return str(value["enum"]) + "::" + str(value["name"])
    counter += 1
    name = f"input_{counter}"
    out.append(f"{value['type']} {name}{{}};")
    for key, operand in value["fields"].items():
        if operand is None:
            continue
        key = {"neutral_reference": "source_reference", "row_tolerance": "solver_options.values[epcsaft::policy::row]"}.get(key, key)
        expression = declaration(operand)
        out.append(f"{name}.{key}={expression};")
    return name


observation = copy.deepcopy(packet["observations"][141])
fields = observation["fields"]
bubble = copy.deepcopy(packet["observations"][0]["fields"]["problem"])
for key in ("temperature", "reactions", "neutral_reference", "reference"):
    if key in fields["problem"]["fields"]:
        bubble["fields"][key] = copy.deepcopy(fields["problem"]["fields"][key])
bubble["fields"]["pressure"]["fields"]["value"] = fields["pressure"]
bubble["fields"]["feed"] = copy.deepcopy(fields["final_problem"]["fields"]["feed"])
fields.update(kind={"enum": "ObservationKind", "name": "SaturationCappedAdditionEnthalpy"},
              final_problem=bubble, basis_amount=1., loading_cap=1.2, loading_bounds=[.05, 1.2],
              loading_grid=[.05, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1., 1.1, 1.2])
row = declaration(observation)
out.append(f"auto row={row};")
thermal = declaration(packet["thermochemistry"])
out.append(f"auto thermal={thermal};")
coords = declaration([packet["coordinates"][i] for i in (0, 6)])
out.append(f"CoordinateSet coordinates; coordinates.coordinates={coords};")
out.append("auto parameters=build_parameters(document);")
Path(sys.argv[2]).write_text("\n".join(out) + "\n")
