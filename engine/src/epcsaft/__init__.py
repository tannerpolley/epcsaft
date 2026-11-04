"""Declare fluid models and read native thermodynamic results.

Temperatures are in K, pressures in Pa and molar densities in mol/m³.
Caloric state energies use J/mol; residual Helmholtz energies divided by RT
are dimensionless. Equilibrium amounts use mol.
Dimensionless logarithms use the natural logarithm unless named ``log10``.
``Mixture`` and ``Pure`` own state evaluation; ``equilibrium`` declares
phase and reaction problems, and ``regression`` declares parameter estimation.
Native result types and exceptions are re-exported from ``_core``.
"""

import logging
import os
import sys
from importlib.metadata import version
from pathlib import Path

__version__ = version("epcsaft")

_dll_directories: list[object] = []
if sys.platform == "win32":
    for _directory in (Path(__file__).resolve().parent, Path(sys.prefix) / "Library" / "bin"):
        if _directory.is_dir():
            _dll_directories.append(os.add_dll_directory(str(_directory)))

from . import _core  # noqa: E402 -- Windows DLL directories must be registered before native imports.
from .eos import (  # noqa: E402
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
from .parameters import ParameterError, Parameters  # noqa: E402
from .reference import Activity, Reference, SolventActivity  # noqa: E402
from ._log import log as _log  # noqa: E402

ReactionThermochemistry = _core.ReactionThermochemistry
reaction_thermochemistry_density = _core.reaction_thermochemistry_density

Stop = _core.Stop
Convergence = _core.Convergence
ConvergenceError = _core.ConvergenceError
tolerances = _core.tolerances
DeclarationError = _core.DeclarationError


def trace(diagnostics: bool | None = None, ipopt_log: str | None = None, ipopt_print_level: int | None = None) -> _core.Capture:
    """Capture native spans and events, and optionally an Ipopt log, inside a ``with`` block.

    Parameters
    ----------
    diagnostics : bool, optional
        Record detailed observations as well as work spans (default: no).
    ipopt_log : str, optional
        Absolute path of a file receiving every Ipopt run in the block (default: none).
    ipopt_print_level : int, optional
        Ipopt journal level of that file (J_MAIN is capped at 5). Default: 12 when the
        ``epcsaft`` logger is enabled for DEBUG, else 5 (one summary line per iteration).

    Returns
    -------
    capture : Capture
        Context manager; read ``spans``, ``events`` and ``trace`` after it closes.

    Notes
    -----
    Capture is diagnostic and changes derivative work (ALG-TERMINATION); the
    ``epcsaft`` logger itself never opens one.
    """
    if ipopt_print_level is None:
        ipopt_print_level = 12 if _log.isEnabledFor(logging.DEBUG) else 5
    return _core.Capture(bool(diagnostics), ipopt_log or "", ipopt_print_level)


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
    """Identify the installed Engine and the source used by its editable build.

    Returns
    -------
    identity : tuple
        For an editable install, the Git commit string and a Boolean stating
        whether ``engine/`` has no uncommitted or untracked files. For a wheel,
        its local ``Path`` and the SHA-256 digest of the wheel's current bytes.

    Raises
    ------
    FileNotFoundError
        A wheel was not installed from a local file, or that file no longer exists.
    importlib.metadata.PackageNotFoundError
        The ``epcsaft`` distribution metadata is unavailable.
    subprocess.CalledProcessError
        Git cannot read the editable checkout's commit or working tree.

    Notes
    -----
    The import rebuilds the editable extension when needed. This identity does
    not establish numerical agreement or physical validity. Wheel identity is
    read from ``direct_url.json``; malformed JSON and file read errors propagate.
    """
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
    "ReactionThermochemistry",
    "reaction_thermochemistry_density",
    "provenance",
    "tolerances",
    "trace",
)
