from __future__ import annotations

import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[5]
PAPER_ROOT = Path(__file__).resolve().parents[2]
PARAMETERS_PATH = ROOT / "data/packets/figiel-2025-reference-electrolytes/2/parameters/parameters.json"
GATES_PATH = PAPER_ROOT / "inputs/docs" / "legacy-derived-gates.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(ROOT).as_posix()
