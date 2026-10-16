from pathlib import Path

from . import _core
from .eos import (
    ActiveParameter,
    DensityActionDirection,
    DensityBranch,
    DensityDiagnostics,
    EosError,
    DensityStructureError,
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
    ThermochemistryPacket,
    ThermochemistryRecord,
    ThermochemistryStatus,
)
from .parameters import ParameterError, Parameters
from .reference import Activity, Reference, SolventActivity

trace = _core.Capture
Stop = _core.Stop
Convergence = _core.Convergence
ConvergenceError = _core.ConvergenceError
tolerances = _core.tolerances
DeclarationError = _core.DeclarationError


def _repr(record: type, *fields: str) -> None:
    record.__repr__ = lambda self: f"{record.__name__}({', '.join(f'{n}={getattr(self, n)!r}' for n in fields)})"  # type: ignore


def _phase_mole_fractions(self: _core.EquilibriumResult, phase: int) -> list[float]:
    index = range(self.phase_count)[phase]  # IndexError outside the declared phases
    n = len(self.mole_fractions) // self.phase_count
    return self.mole_fractions[index * n : (index + 1) * n]


_repr(_core.State, "pressure", "molar_density", "compressibility_factor", "density_branch", "residual_helmholtz")
_repr(_core.EquilibriumResult, "success", "temperature", "pressure", "amounts", "mole_fractions")
_repr(_core.ContinuationResult, "success", "stop", "reached_path", "accepted_steps", "message")
_repr(_core.FitResult, "status", "usable", "final_cost", "iterations", "physical")
_repr(_core.ElectrolyteActivity, "mean_ionic_molality_scale", "osmotic_coefficient", "single_ion")
setattr(_core.EquilibriumResult, "phase_mole_fractions", _phase_mole_fractions)


def provenance() -> tuple[Path, str] | tuple[str, bool]:
    """Identify the running Engine: a wheel install's file path and SHA-256 (its direct_url.json), or
    for the editable install the commit and whether engine/ is clean, with no uncommitted or untracked
    files (the import rebuilt the extension from that tree)."""
    import hashlib
    import json
    import subprocess
    from importlib import metadata
    from urllib.parse import unquote, urlparse
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
    "DensityStructureError",
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
    "SolventActivity",
    "State",
    "Stop",
    "ThermochemistryPacket",
    "ThermochemistryRecord",
    "ThermochemistryStatus",
    "provenance",
    "tolerances",
    "trace",
)
