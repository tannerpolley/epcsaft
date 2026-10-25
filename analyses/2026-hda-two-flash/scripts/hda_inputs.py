"""HDA packet inputs (#73): the hda-pc-saft parameter document and the NASA9 thermochemistry."""

import json
from pathlib import Path

import epcsaft
from epcsaft.equilibrium import Free, MoleFractions, Phase, Pinned, Problem

ROOT = Path(__file__).resolve().parents[3]
ANALYSIS = Path(__file__).resolve().parents[1]
PARAMETERS = ROOT / "data/packets/hda-pc-saft/1/parameters/parameters.json"
NASA9 = ROOT / "data/packets/nasa-glenn-2002-hda-species/1/parameters/nasa9-records.json"
PAIR, P, X = ("benzene", "toluene"), 105000.0, [0.5, 0.5]  # the #73 equimolar binary path


def parameters(components, pure_override=None):
    """The packet document for `components`; `pure_override[name][family]` replaces a pure value."""
    document = json.loads(PARAMETERS.read_text())
    for record in document["components"]:
        for item in record["coefficients"]:
            item["value"]["magnitude"] = (pure_override or {}).get(record["component_id"], {}).get(
                item["family"], item["value"]["magnitude"])
    return epcsaft.Parameters.from_mapping(document, components=list(components))


# NASA9 records at 1 bar, loaded once; native preparation scales the coefficients by R_source/R_engine so cp, h and
# s keep their source values. Intervals are [lower, upper) except the last, which is closed.
thermochemistry = epcsaft.ThermochemistryPacket(NASA9).select


def boundary_problem(kind, guess):
    """Bubble ("bubble": vapor amount pinned at 0) or dew temperature of the equimolar pair at P."""
    liquid = Phase("L", kind="liquid", composition_guess=X, amount=Pinned(0.0) if kind == "dew" else None)
    vapor = Phase("V", kind="vapor", composition_guess=X, amount=Pinned(0.0) if kind == "bubble" else None)
    return Problem(phases=[liquid, vapor], T=Free(guess), P=P, feed=MoleFractions(dict(zip(PAIR, X))))
