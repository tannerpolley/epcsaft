"""Mechanically transport the retained MEA declarations to native C++ inputs."""
import hashlib
import importlib.util
import json
from pathlib import Path

import epcsaft
from epcsaft import _core, regression as reg

ROOT = Path(__file__).resolve().parents[3]
LEAF = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "analyses/2026-fitting-state-reuse/inputs/mea-initial.json"
OUT = ROOT / "build/issue-510"
spec = importlib.util.spec_from_file_location("retained_declarations", SOURCE.parents[1] / "scripts/replay.py")
source = importlib.util.module_from_spec(spec)
spec.loader.exec_module(source)
data = json.loads(SOURCE.read_text())
documents = []
original = _core.build_parameters


def capture(document):
    documents.append(document)
    return original(document)


_core.build_parameters = capture
try:
    parameters = epcsaft.Parameters.from_mapping(source.declared(data))
finally:
    _core.build_parameters = original
lines = ["// Generated declarations; edit the retained source, not this file.", "FitRequest mea_request() {"]
serial = 0


def expression(value):
    global serial
    if isinstance(value, bool):
        return str(value).lower()
    if isinstance(value, (int, float)):
        return repr(value)
    if isinstance(value, str):
        return json.dumps(value)
    if isinstance(value, (list, tuple)):
        return "{" + ",".join(expression(item) for item in value) + "}"
    if hasattr(value, "name") and hasattr(value, "value"):
        return type(value).__name__ + "::" + value.name
    name = "declaration_" + str(serial)
    serial += 1
    lines.append(type(value).__name__ + " " + name + ";")
    for field in dir(value):
        # solver_options is the read-only view of row_tolerance, which the native solver options own (#486).
        if field.startswith("_") or field == "solver_options" or callable(getattr(value, field)) or getattr(value, field) is None:
            continue
        target = {"neutral_reference": "source_reference",
                  "row_tolerance": "solver_options.values[epcsaft::policy::row]"}.get(field, field)
        assigned = expression(getattr(value, field))
        lines.append(name + "." + target + "=" + assigned + ";")
    return name


document = expression(documents[0])
coordinate = reg.coordinate(parameters, [("k_ij", ("water", "carbamate-anion"), 1.0)],
    **{key: data["coordinates"][0][key] for key in ("origin", "scale")}, bounds=(-0.5, 0.5))
coordinate_name = expression(coordinate)
observations = []
for state in data["states"][:16]:
    problem = source.declaration(state["problem"])
    observations.append(expression(reg.observation(parameters, "partial_pressure", problem, observed=1.0,
        form="difference", scale=1.0, mass_basis=False, phase=problem.phases[1].name, species=["carbon-dioxide"])))
lines += ["FitRequest request; request.parameters=build_parameters(" + document + ");",
          "request.coordinates.coordinates={" + coordinate_name + "};", "request.observations={" + ",".join(observations) + "};", "return request;", "}"]
if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    header = OUT / "mea-declarations.hpp"
    header.write_text("\n".join(lines) + "\n")
    (LEAF / "inputs/provenance.json").write_text(json.dumps({"source": str(SOURCE.relative_to(ROOT)),
        "source_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(), "generated": str(header.relative_to(ROOT)),
        "generated_sha256": hashlib.sha256(header.read_bytes()).hexdigest(), "states": [s["identity"] for s in data["states"][:16]],
        "component_ids": parameters.component_ids, "base": json.loads((LEAF / "inputs/campaign-freeze.json").read_text())["comparison_base"]}, separators=(",", ":")) + "\n")


def declarations():
    return parameters, coordinate, [source.declaration(state["problem"]) for state in data["states"][:16]]
