"""Compile the retained PC-SAFT source data to bit-exact C++ arrays; no evaluation logic."""
import json
import sys
from pathlib import Path

doc = json.loads(Path(sys.argv[1]).read_text())
assert len(doc["intervals"]) == 10 and len(doc["w_edges"]) == 11
assert len(doc["Ttilde_crit"]) == 11
assert doc["w_edges"][0] == 1 / 64 and doc["w_edges"][-1] == 1.0
assert all(a < b for a, b in zip(doc["w_edges"], doc["w_edges"][1:]))
segments = []
offsets = []
nodes = []

def collection(values, lower, upper):
    assert values[0]["xmin"] == lower and values[-1]["xmax"] == upper
    assert all(a["xmax"] == b["xmin"] for a, b in zip(values, values[1:]))
    assert all(len(v["coef"]) == 17 and v["xmin"] < v["xmax"] for v in values)
    offsets.append(len(segments))
    segments.extend([[v["xmin"], v["xmax"], *v["coef"]] for v in values])

assert doc["Ttilde_crit"][0]["xmin"] <= 1 / 64
collection(doc["Ttilde_crit"], doc["Ttilde_crit"][0]["xmin"], 1.0)
for i, interval in enumerate(doc["intervals"]):
    assert interval["wnodes"][0] == doc["w_edges"][i + 1]
    assert interval["wnodes"][-1] == doc["w_edges"][i]
    assert all(a > b for a, b in zip(interval["wnodes"], interval["wnodes"][1:]))
    assert len(interval["wnodes"]) == 17
    nodes.append(interval["wnodes"])
    for phase in ("liquid", "vapor"):
        assert len(interval[phase]) == 17
        for values in interval[phase]:
            collection(values, 0.0, 1.0)
offsets.append(len(segments))

def array(name, rows, width):
    values = ["{" + ",".join(float(v).hex() for v in row) + "}" for row in rows]
    return "inline constexpr double " + name + "[][" + str(width) + "] = {\n" + ",\n".join(values) + "};\n"

text = "// Generated from engine/data/superancillary/pcsaft.json; NIST notice accompanies that file.\n#pragma once\nnamespace epcsaft::superancillary {\n"
text += array("segments", segments, 19)
text += array("nodes", nodes, 17)
text += array("w_edges", [doc["w_edges"]], 11)
text += array("lower_temperature", [doc["Ttilde_min_coefficients"]], 2)
text += "inline constexpr unsigned offsets[] = {" + ",".join(map(str, offsets)) + "};\n}\n"
Path(sys.argv[2]).write_text(text)
