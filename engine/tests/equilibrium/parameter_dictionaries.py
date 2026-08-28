from __future__ import annotations

from tests.eos._cases import (
    neutral_pure,
    set_component_values,
)

METHANE_PARAMETERS = neutral_pure("methane")
ETHANE_PARAMETERS = neutral_pure("ethane")
PROPANE_PARAMETERS = neutral_pure("propane")
set_component_values(PROPANE_PARAMETERS, "molar_mass", (0.044096,))

PURE_SATURATION_PARAMETERS = {
    "methane": METHANE_PARAMETERS,
    "ethane": ETHANE_PARAMETERS,
    "propane": PROPANE_PARAMETERS,
}
