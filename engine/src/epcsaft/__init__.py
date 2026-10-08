import hashlib
import json
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


def installed_wheel() -> tuple[Path, str]:
    """Return the path and SHA-256 of the wheel file this package was installed from (its direct_url.json)."""
    url = json.loads(metadata.distribution("epcsaft").read_text("direct_url.json") or "{}").get("url", "")
    path = Path(unquote(urlparse(url).path)).resolve()
    if not url.startswith("file:") or not path.is_file():
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
    "installed_wheel",
    "tolerances",
    "trace",
)
