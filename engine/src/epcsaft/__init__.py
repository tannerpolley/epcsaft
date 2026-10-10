import hashlib
import json
import subprocess
from importlib import metadata
from pathlib import Path
from urllib.parse import unquote, urlparse

from . import _core
from .eos import (
    ActiveParameter,
    DensityActionDirection,
    DensityBranch,
    DensityDiagnostics,
    EosError,
    GAS_CONSTANT_J_PER_MOL_K,
    IdealCorrelation,
    IdealDIPPR107,
    IdealInterval,
    IdealNASA7,
    IdealNASA9,
    IdealPolynomial,
    IdealShomate,
    Mixture,
    ParameterActionDirection,
    ParameterFamily,
    PressureActionDirection,
    PropertyActionBatch,
    PropertyActionDiagnostics,
    PropertyActionResult,
    PropertyActionStatus,
    PropertyObservable,
    Pure,
    State,
    ThermochemistryRecord,
    ThermochemistryStatus,
)
from .parameters import ParameterError, Parameters
from .reference import Activity, Reference

trace = _core.Capture
Stop = _core.Stop
Convergence = _core.Convergence
ConvergenceError = _core.ConvergenceError
tolerances = _core.tolerances
DeclarationError = _core.DeclarationError


def provenance() -> tuple[Path, str] | tuple[str, bool]:
    """Identify the running Engine: a wheel install's file path and SHA-256 (its direct_url.json), or
    for the editable install the commit and whether engine/ is clean, with no uncommitted or untracked
    files (the import rebuilt the extension from that tree)."""
    url = json.loads(metadata.distribution("epcsaft").read_text("direct_url.json") or "{}")
    if url.get("dir_info", {}).get("editable"):
        git = ["git", "-C", str(Path(__file__).resolve().parents[2])]
        commit, changes = (subprocess.run([*git, *args], capture_output=True, text=True, check=True).stdout
                           for args in (["rev-parse", "HEAD"], ["status", "--porcelain", "--", "."]))
        return commit.strip(), not changes
    path = Path(unquote(urlparse(url.get("url", "")).path)).resolve()
    if not url.get("url", "").startswith("file:") or not path.is_file():
        raise FileNotFoundError(f"epcsaft was not installed from a wheel file that still exists: {url!r}")
    return path, hashlib.sha256(path.read_bytes()).hexdigest()

__all__ = (
    "ActiveParameter",
    "Activity",
    "Convergence",
    "ConvergenceError",
    "DeclarationError",
    "DensityActionDirection",
    "DensityBranch",
    "DensityDiagnostics",
    "EosError",
    "GAS_CONSTANT_J_PER_MOL_K",
    "IdealCorrelation",
    "IdealDIPPR107",
    "IdealInterval",
    "IdealNASA7",
    "IdealNASA9",
    "IdealPolynomial",
    "IdealShomate",
    "Mixture",
    "ParameterActionDirection",
    "ParameterError",
    "ParameterFamily",
    "Parameters",
    "PressureActionDirection",
    "PropertyActionBatch",
    "PropertyActionDiagnostics",
    "PropertyActionResult",
    "PropertyActionStatus",
    "PropertyObservable",
    "Pure",
    "Reference",
    "State",
    "Stop",
    "ThermochemistryRecord",
    "ThermochemistryStatus",
    "provenance",
    "tolerances",
    "trace",
)
