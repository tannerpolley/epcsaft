"""Reactive MEA liquid cold-start cases (#154, #157): the problems exactly as the applications declare them.

ABSORBER: one declared liquid, Amounts feed of apparent CO2/MEA/water, no composition guess
(MEA-Absorption-Column 4f2d12e9 `reactive_bundle.engine_liquid` / `casadi_reactive.EnginePhase._problem`).
MEA: packet speciation states through MEA's own request mapping (calibration-misfit e522dc4).
"""
from __future__ import annotations

import hashlib
import json
import math
import os
import sys
from pathlib import Path

import epcsaft
from epcsaft import equilibrium as q

ROOT = Path(__file__).resolve().parents[3]
CLONES = Path(os.environ.get("DOWNSTREAM_CLONES", ROOT / "downstream/.worktrees"))  # read-only application inputs
ABSORBER = CLONES / "MEA-Absorption-Column-157-inputs"  # absorber 4f2d12e9
MEA = CLONES / "MEA-Thermodynamics-calibration-misfit"  # MEA e522dc4
DATA = ABSORBER / "src/mea_absorption_column/data/epcsaft_datasets/MEA_greenfield_exploratory"
FIELDS = ("a", "b", "c", "d", "reference_temperature", "temperature_min", "temperature_max", "standard_state_id")

# #157 reproduction rows: T [K], P [Pa], apparent CO2/MEA/water feed [mol].
ABSORBER_STATES = [
    (329.39, 110900.0, (2.43232, 9.72928, 76.0902)),
    (329.39, 110901.0, (2.43232, 9.72928, 76.0902)),
    (318.15, 111100.0, (1.991416, 9.7292821, 76.9751925)),
    (318.15, 111110.0, (1.991416, 9.7292821, 76.9751925)),
    (318.15, 111200.0, (1.991416, 9.7292821, 76.9751925)),
    (318.15, 110900.0, (2.550636, 9.7292821, 76.9751925)),
]
MEA_STATES = ("Bottinger2008_state_042", "Bottinger2008_state_046", "Matin2012_state_016")


def absorber_thermochemistry(ids):
    """The absorber's physical ideal-gas records through its own reader (`reactive_bundle.ideal_gas_record`)."""
    import importlib.util  # noqa: PLC0415
    spec = importlib.util.spec_from_file_location(
        "reactive_bundle", ABSORBER / "src/mea_absorption_column/Thermodynamics/reactive_bundle.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.ideal_gas_record(DATA / "ideal-gas-thermochemistry.json", ids)


def absorber_liquid():
    """(model, problem factory) for the absorber's nine-species reactive liquid."""
    record = json.loads((DATA / "engine-reactions.json").read_text())
    mapping = json.loads((DATA / "parameters.json").read_text())
    assert hashlib.sha256((DATA / "parameters.json").read_bytes()).hexdigest() == record["parameter_sha256"]
    for family in mapping["model_families"]:
        if family.get("kind") == "electrolyte" and family.get("choice") == "born":
            for key, value in record["born_runtime_defaults"].items():
                family.setdefault(key, value)
    reactions = [q.Reaction(row["stoichiometry"], name=row["name"],
                            correlation=q.ReactionLogPolynomial(*(row["correlation"][k] for k in FIELDS)),
                            reference=q.ReactionReference(getattr(q.ReferenceSourceBasis, row["reference"]["source_basis"]),
                                                          reference_pressure_pa=row["reference"]["reference_pressure_pa"]))
                 for row in record["reactions"]]
    ids = [c["component_id"] for c in mapping["components"]]
    neutral = q.NeutralReference(**record["neutral_reference"])
    model = epcsaft.Mixture(epcsaft.Parameters.from_mapping(mapping), thermochemistry=absorber_thermochemistry(ids))

    def problem(temperature, pressure, feed):
        return q.Problem(phases=[q.Phase("liquid", kind="liquid")], T=temperature, P=pressure,
                         feed=q.Amounts(dict(zip(ids[:3], feed))), reactions=reactions, neutral_reference=neutral)

    return model, problem, ids


def mea_evaluation():
    """MEA's shared_evaluation module, its packet and the pre-refit / refit-C parameter mappings."""
    sys.argv = [sys.argv[0], str(MEA), "timing"]
    sys.path.insert(0, str(ROOT / "analyses/2026-mea-reactive-fit-timing/scripts"))
    import run as timing  # noqa: PLC0415 - #152's loader: record hash, packet, reactions, refit C
    base = epcsaft.Mixture(epcsaft.Parameters.from_mapping(timing.mapping()))
    refit = epcsaft.Mixture(epcsaft.Parameters.from_mapping(
        timing.mapping(dict(zip(timing.REFIT_C["parameters"], timing.REFIT_C["start"])))))
    return timing, {"pre_refit": base, "refit_c": refit}


def relative_composition_difference(left, right) -> float:
    return max(abs(a - b) / max(abs(a), abs(b), 1e-300) for a, b in zip(left, right))


def finite(value) -> float | None:
    return value if value is not None and math.isfinite(value) else None
