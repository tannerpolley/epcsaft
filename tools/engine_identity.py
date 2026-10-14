"""Identify the loaded Engine images in an installed or prepared environment."""

import hashlib
import json
import os
import subprocess
import sys
from importlib import metadata
from pathlib import Path
from urllib.parse import unquote, urlparse


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def loaded_identity():
    from epcsaft import _core

    core = Path(_core.__file__).resolve(strict=True)
    direct = json.loads(metadata.distribution("epcsaft").read_text("direct_url.json") or "{}")
    source = Path(unquote(urlparse(direct.get("url", "")).path))
    wheel = source if source.suffix == ".whl" and source.is_file() else None
    mappings = Path("/proc/self/maps").read_text().splitlines()
    for line in mappings:
        fields = line.split(maxsplit=5)
        if (len(fields) == 6 and (core.name in fields[5] or "libepcsaft_native.so" in fields[5])
                and (fields[5].endswith(" (deleted)") or Path(fields[5]).stat().st_ino != int(fields[4]))):
            raise RuntimeError("loaded Engine image was replaced; identify it in a fresh interpreter")
    native = {Path(line.split(maxsplit=5)[5]).resolve() for line in mappings
              if len(line.split(maxsplit=5)) == 6 and line.split(maxsplit=5)[5].endswith("/libepcsaft_native.so")}
    if len(native) > 1:
        raise RuntimeError(f"multiple loaded native libraries: {native}")
    numerical = next(iter(native), core)
    options = core.parent / "meson-info/intro-buildoptions.json"
    editable = direct.get("dir_info", {}).get("editable", False)
    revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=source.parent, text=True).strip() if editable else None
    difference = subprocess.check_output(["git", "diff", "HEAD", "--", "engine"], cwd=source.parent) if editable else b""
    return {"interpreter": sys.executable, "python": sys.version, "direct_url": direct,
                "wheel": str(wheel) if wheel else None, "wheel_sha256": digest(wheel) if wheel else None,
                "source_revision": revision, "source_diff_sha256": hashlib.sha256(difference).hexdigest() if editable else None,
                "build_configuration": str(options) if options.is_file() else None,
                "build_configuration_sha256": digest(options) if options.is_file() else None,
                "extension": str(core), "extension_sha256": digest(core),
                "native_library": str(numerical), "native_library_sha256": digest(numerical)}


def asl_library(identity=None):
    identity = identity or loaded_identity()
    path = Path(os.environ.get("EPCSAFT_ASL_LIBRARY", str(Path(identity["extension"]).parent / "libepcsaft_asl.so")))
    return path.resolve(strict=True)


if __name__ == "__main__":
    saved_stdout = os.dup(1)
    os.dup2(2, 1)
    try:
        identity = loaded_identity()
        sys.stdout.flush()
    finally:
        os.dup2(saved_stdout, 1)
        os.close(saved_stdout)
    print(json.dumps(identity, sort_keys=True))
