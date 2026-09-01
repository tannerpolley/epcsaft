"""Fail closed unless the active interpreter uses the declared unified Engine wheel."""

from __future__ import annotations

import argparse
import json
import re
from importlib import import_module, metadata
from pathlib import Path
from urllib.parse import unquote, urlparse
from zipfile import ZipFile

PUBLIC_MODULES = ("epcsaft", "epcsaft.equilibrium", "epcsaft.regression")
NATIVE_ENTRYPOINTS = (
    "_compile_eos_regression_evaluator",
    "_compile_equilibrium_regression_evaluator",
    "_solve_grepe",
    "_solve_least_squares",
)
NATIVE_MEMBER_PREFIX = "epcsaft/_core."
TRANSIENT_WHEEL_RUNPATH = b"/wheel/platlib/lib"


def normalized_distribution_name(name: str) -> str:
    return re.sub(r"[-_.]+", "-", name).lower()


def verify(wheel: Path) -> None:
    wheel = wheel.resolve()
    if wheel.suffix != ".whl" or not wheel.is_file():
        raise SystemExit(f"Engine wheel does not exist: {wheel}")

    with ZipFile(wheel) as archive:
        members = archive.namelist()
        matches = [
            member
            for member in members
            if member.startswith("epcsaft/") and member.endswith(".so")
        ]
        if len(matches) != 1 or not matches[0].startswith(NATIVE_MEMBER_PREFIX):
            raise SystemExit(f"expected only the epcsaft._core extension, found {matches}")
        if TRANSIENT_WHEEL_RUNPATH in archive.read(matches[0]):
            raise SystemExit(f"transient build runpath embedded in {matches[0]}")

    epcsaft_distributions = [
        distribution
        for distribution in metadata.distributions()
        if normalized_distribution_name(distribution.metadata.get("Name", "")).startswith(
            "epcsaft"
        )
    ]
    names = [distribution.metadata.get("Name", "") for distribution in epcsaft_distributions]
    if names != ["epcsaft"]:
        raise SystemExit(f"expected only epcsaft, found {names}")

    distribution = epcsaft_distributions[0]
    direct_url = json.loads(distribution.read_text("direct_url.json") or "{}")
    if direct_url.get("dir_info", {}).get("editable", False):
        raise SystemExit("editable Engine install is forbidden")
    installed_from = Path(unquote(urlparse(direct_url.get("url", "")).path)).resolve()
    if installed_from != wheel:
        raise SystemExit(f"unexpected Engine install source: {installed_from}")

    installed_root = Path(distribution.locate_file("")).resolve()
    origins = []
    for name in PUBLIC_MODULES:
        module = import_module(name)
        origin = Path(module.__file__).resolve()
        if installed_root not in origin.parents:
            raise SystemExit(f"{name} is outside the installed wheel: {origin}")
        origins.append(f"{name}={origin.relative_to(installed_root)}")
    core = import_module("epcsaft._core")
    if any(not callable(getattr(core, name, None)) for name in NATIVE_ENTRYPOINTS):
        raise SystemExit(f"unified Engine extension lacks {NATIVE_ENTRYPOINTS}")

    print(f"epcsaft {distribution.version} installed from {wheel}")
    print("module origins: " + ", ".join(origins))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--wheel", type=Path, required=True)
    arguments = parser.parse_args()
    verify(arguments.wheel)


if __name__ == "__main__":
    main()
