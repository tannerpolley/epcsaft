"""Diagnose the native-unit signal floor without changing the specified derivative checks."""
import copy
import csv
import json
from pathlib import Path

from epcsaft import Parameters, Pure, regression
from epcsaft.eos import IdealCorrelation, IdealInterval, IdealPolynomial, ThermochemistryRecord, GAS_CONSTANT_J_PER_MOL_K

SET = Path(__file__).resolve().parents[1]
source = json.loads((SET.parents[1] / "data/packets/held-cameretti-sadowski-2008/3/parameters/parameters.json").read_text())
base = Parameters.from_mapping(source, components=["water"])
record = ThermochemistryRecord(1e5, [[IdealInterval(150, 500, True, True,
    IdealCorrelation(IdealPolynomial([4*GAS_CONSTANT_J_PER_MOL_K], 300, 0, 0), 1e5))]])
chart = [("segment_count", 5e-4), ("segment_diameter", 5e-4), ("dispersion_energy_over_k", 0.05),
         ("association_energy_over_k", 1.0), ("association_volume", 1e-5)]
coordinates = [regression.coordinate(base, [(family, (("water", "a"), ("water", "b")) if family.startswith("association") else "water", 1)],
    origin=0, scale=1) for family, _ in chart]
set_ = regression.CoordinateSet()
set_.coordinates = coordinates
model = Pure(base, thermochemistry=record, active_parameters=set_.active())


def shifted(family, shift):
    mapping = copy.deepcopy(source)
    if family.startswith("association"):
        edge = next(e for e in mapping["topology"]["edges"] if e["endpoint_a"]["component_id"] == "water")
        node = edge["energy_over_k" if family == "association_energy_over_k" else "volume"]
    elif family == "segment_diameter":
        node = next(c for c in mapping["correlations"] if c.get("component_id") == "water"
                    and c["family"] == "segment_diameter")["constant"]
    else:
        component = next(c for c in mapping["components"] if c["component_id"] == "water")
        node = next(c for c in component["coefficients"] if c["family"] == family)
    node["value"]["magnitude"] += shift
    return Pure(Parameters.from_mapping(mapping, components=["water"]), thermochemistry=record)


rows = []
for T, P in [(298.15, 1e5), (400, 1e6), (420, 1e6), (440, 1e6), (450, 1e6)]:
    item = regression.isobaric_heat_capacity(base, T, P, {"water": 1}, phase="liquid", observed=1,
                                            form="log_ratio")
    exact = regression.differentiate(model, item).value
    for k, (family, h) in enumerate(chart):
        coarse = (regression.predict(shifted(family, h), item).value
                  - regression.predict(shifted(family, -h), item).value)/(2*h)
        fine = (regression.predict(shifted(family, h/2), item).value
                - regression.predict(shifted(family, -h/2), item).value)/h
        richardson = (4*fine-coarse)/3
        rows.append(dict(T_K=T, P_Pa=P, family=family, step_native=h, exact=exact[k], coarse=coarse,
            fine=fine, richardson=richardson, signal=abs(fine), signal_floor=1e-2,
            fine_error=abs(exact[k]-fine), richardson_error=abs(exact[k]-richardson)))
        print(T, family, exact[k], fine, "signal_floor_met", abs(fine)>=1e-2, flush=True)
with (SET / "results/signal-probe.csv").open("w", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
