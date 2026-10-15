"""Fail closed unless the active interpreter uses the declared unified Engine wheel.

Keep it working on older wheels: time-key-solves and refresh-main-wheel verify wheels from
earlier commits, so use only the wheel's files and metadata, not newer epcsaft functions."""

from __future__ import annotations

import argparse
import base64
import csv
import hashlib
import json
import re
from importlib import import_module, metadata
from pathlib import Path
from urllib.parse import unquote, urlparse
from zipfile import ZipFile

PUBLIC_MODULES = ("epcsaft", "epcsaft.equilibrium")
NATIVE_ENTRYPOINTS = (
    "build_parameters",
    "HelmholtzModel",
    "infinite_dilution_reference",
    "electrolyte_activity",
    "compile_problem",
    "sensitivity",
    "solve_equilibrium",
)
NATIVE_MEMBER_PREFIX = "epcsaft/_core."
ASL_MEMBER = "epcsaft/libepcsaft_asl.so"
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
        extensions = [member for member in matches if member.startswith(NATIVE_MEMBER_PREFIX)]
        if (len(extensions) != 1 or ASL_MEMBER not in matches
                or len(matches) != 2 + ("epcsaft/libepcsaft_native.so" in matches)):
            raise SystemExit(
                f"expected one epcsaft._core extension and {ASL_MEMBER}, found {matches}"
            )
        for member in matches:
            if TRANSIENT_WHEEL_RUNPATH in archive.read(member):
                raise SystemExit(f"transient build runpath embedded in {member}")
        (record_name,) = [member for member in members if member.endswith(".dist-info/RECORD")]
        record = {
            path: digest
            for path, digest, _ in csv.reader(archive.read(record_name).decode().splitlines())
            if digest
        }

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
    # The installed bytes must be the named wheel's bytes, not just its path.
    for path, digest in record.items():
        algorithm, expected = digest.split("=", 1)
        installed = installed_root / path
        if not installed.is_file():
            raise SystemExit(f"{installed} is missing; wheel RECORD lists {path}")
        actual = base64.urlsafe_b64encode(
            hashlib.new(algorithm, installed.read_bytes()).digest()
        ).rstrip(b"=").decode()
        if actual != expected:
            raise SystemExit(
                f"{installed} differs from {wheel.name} RECORD: {algorithm}={actual}, "
                f"expected {digest}"
            )
    extra = sorted(
        str(path.relative_to(installed_root))
        for path in (installed_root / "epcsaft").rglob("*")
        if path.is_file()
        and "__pycache__" not in path.parts
        and path.relative_to(installed_root).as_posix() not in record
    )
    if extra:
        raise SystemExit(f"installed epcsaft files absent from {wheel.name} RECORD: {extra}")
    installed_asl = Path(distribution.locate_file(ASL_MEMBER)).resolve()
    if not installed_asl.is_file() or installed_root not in installed_asl.parents:
        raise SystemExit(f"ASL library is outside the installed wheel: {installed_asl}")
    origins = []
    for name in PUBLIC_MODULES:
        module = import_module(name)
        origin = Path(module.__file__).resolve()
        if installed_root not in origin.parents:
            raise SystemExit(f"{name} is outside the installed wheel: {origin}")
        origins.append(f"{name}={origin.relative_to(installed_root)}")
    core = import_module("epcsaft._core")
    missing = [name for name in NATIVE_ENTRYPOINTS if not callable(getattr(core, name, None))]
    if missing:
        raise SystemExit(f"current Engine extension lacks {missing}")

    print(f"epcsaft {distribution.version} installed from {wheel}; {len(record)} RECORD hashes match")
    print("module origins: " + ", ".join(origins))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--wheel", type=Path, required=True)
    arguments = parser.parse_args()
    verify(arguments.wheel)


if __name__ == "__main__":
    main()
