"""Panel contracts shared by Figiel data generation and rendering."""

from __future__ import annotations


FIGURE_8_PANELS = (
    ("a)", "LiBr", 5.0, 3.0),
    ("b)", "NaI", 1.5, 1.125),
    ("c)", "NaBr", 1.5, 1.125),
)

FIGURE_9_PANELS = (
    ("a)", "NaBr", "water-methanol", "methanol", "x_Methanol", 3.0),
    ("b)", "NaBr", "water-ethanol", "ethanol", "x_Ethanol", 4.0),
    ("c)", "NaCl", "water-methanol", "methanol", "x_Methanol", 2.0),
    ("d)", "NaCl", "water-ethanol", "ethanol", "x_Ethanol", 2.0),
)
