"""Pure-fluid input for the retained reduced-row comparison; no calculation."""
import json
from pathlib import Path

from epcsaft import Pure
from tests.native_inputs import parameters


def declarations():
    owner = Path(__file__).resolve().parents[2] / '2026-competitor-cross-check/inputs/cases.json'
    component = json.loads(owner.read_text())['models']['C1']['components'][0]
    return Pure(parameters([(1000 * component['molar_mass'], component['m'], component['sigma'], component['epsilon'])]))
