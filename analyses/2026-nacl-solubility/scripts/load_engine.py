"""Use this checkout's installed Engine and retain its loaded extension identity."""
from pathlib import Path
from epcsaft import _core

ROOT = Path(__file__).resolve().parents[3]
LIB = Path(_core.__file__).resolve()
assert ROOT in LIB.parents, "Prepare this checkout's .venv with README Checks"
