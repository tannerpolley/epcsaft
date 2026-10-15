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
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def prepared_python(value, environment=None):
    """Resolve a Python path, venv, or prepared root without creating an environment."""
    env = dict(os.environ if environment is None else environment)
    path = Path(value).expanduser().absolute()
    if path.is_dir():
        root = path
        marker = root / "python"
        if marker.is_file() and not marker.is_symlink():
            target = Path(marker.read_text().strip())
            path = target if target.is_absolute() else root / target
        else:
            path = next((p for p in (root / "venv/bin/python", root / "bin/python") if p.is_file()), path)
    path = Path(os.path.abspath(path))
    if not path.is_file() or not os.access(path, os.X_OK):
        raise ValueError(f"prepared input must name a Python interpreter or root: {value}")
    return path, env


def _source_root(direct):
    selected = os.environ.get("EPCSAFT_SOURCE_ROOT")
    if selected:
        return Path(selected).expanduser().resolve(strict=True)
    url = urlparse(direct.get("url", ""))
    if direct.get("dir_info", {}).get("editable") and url.scheme == "file":
        path = Path(unquote(url.path)).resolve(strict=True)
        return path.parent if path.name == "engine" else path
    return None


def _source_identity(root):
    if root is None:
        return None, None
    revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip()
    diff = subprocess.check_output(["git", "diff", "--binary", "HEAD"], cwd=root)
    untracked = subprocess.check_output(
        ["git", "ls-files", "--others", "--exclude-standard", "-z"], cwd=root
    ).split(b"\0")
    digestor = hashlib.sha256(diff)
    for name in filter(None, untracked):
        path = root / os.fsdecode(name)
        if path.is_file():
            digestor.update(name + b"\0" + path.read_bytes())
    return revision, digestor.hexdigest()


def _build_dir(core):
    value = os.environ.get("EPCSAFT_BUILD_DIR")
    if value:
        path = Path(value).expanduser().resolve(strict=True)
        if not (path / "meson-info/intro-buildoptions.json").is_file():
            raise RuntimeError(f"prepared Meson identity is missing: {path}")
        return path
    return next((p for p in core.parents if (p / "meson-info/intro-buildoptions.json").is_file()), None)


def _build_source(build):
    files = json.loads((build / "meson-info/intro-buildsystem_files.json").read_text())
    roots = {Path(path).resolve().parents[1] for path in files
             if Path(path).name == "meson.build" and Path(path).parent.name == "engine"}
    if len(roots) != 1:
        raise RuntimeError(f"Meson build does not identify one Engine source root: {build}")
    return roots.pop()


def _build_identity(build, source):
    if build is None:
        return None, None, None
    names = ("meson-info/intro-buildoptions.json", "meson-info/intro-compilers.json",
             "meson-info/intro-buildsystem_files.json",
             "meson-info/intro-targets.json", "meson-private/cmd_line.txt", "build.ninja", "toolchain-key")
    files = {name: digest(build / name) for name in names if (build / name).is_file()}
    configuration = hashlib.sha256(json.dumps(files, sort_keys=True).encode()).hexdigest()
    objects = [entry["output"] for entry in json.loads((build / "compile_commands.json").read_text())
               if "engine/native/" in entry["file"]]
    dependencies = subprocess.check_output(["ninja", "-C", str(build), "-t", "deps", *objects], text=True)
    object_time = None
    for line in dependencies.splitlines():
        if line and not line.startswith(" "):
            if not line.endswith("(VALID)"):
                raise RuntimeError(f"stale native dependency record: {line}")
            object_time = (build / line.split(": #deps", 1)[0]).stat().st_mtime_ns
        elif line.startswith("    ") and (build / line.strip()).stat().st_mtime_ns > object_time:
            raise RuntimeError(f"native input changed after compilation: {line.strip()}")
    locks = (source / "uv.lock", source / "engine/build/subprojects/.wraplock") if source else ()
    locked = {str(path.relative_to(source)): digest(path) for path in locks if path.is_file()}
    return build, configuration, hashlib.sha256(json.dumps(locked, sort_keys=True).encode()).hexdigest() if locked else None


def _verify_products(core, native, asl, build):
    products = {}
    for target in json.loads((build / "meson-info/intro-targets.json").read_text()):
        for filename in target.get("filename", []):
            name = Path(filename).name
            if name.startswith("_core") or name in ("libepcsaft_native.so", "libepcsaft_asl.so"):
                path = Path(filename).resolve(strict=True)
                products[path.name.split(".")[0]] = path
    for actual, name in ((core, "_core"), (native, "libepcsaft_native"), (asl, "libepcsaft_asl")):
        expected = products.get(name)
        if expected and (actual.resolve(strict=True) != expected or digest(actual) != digest(expected)):
            raise RuntimeError(f"loaded Engine image does not match selected Meson product: {actual}")
        if name in ("_core", "libepcsaft_asl") and not expected:
            raise RuntimeError(f"selected Meson build has no {name} product: {build}")


def _mapped(path, mappings):
    path = path.resolve(strict=True)
    found = []
    for line in mappings:
        fields = line.split(maxsplit=5)
        if len(fields) != 6 or fields[5].endswith(" (deleted)"):
            continue
        mapped = Path(fields[5])
        if mapped.resolve() != path:
            continue
        stat = path.stat()
        major, minor = (int(v, 16) for v in fields[3].split(":"))
        if stat.st_ino != int(fields[4]) or (os.major(stat.st_dev), os.minor(stat.st_dev)) != (major, minor):
            raise RuntimeError(f"loaded Engine image was replaced; identify it in a fresh interpreter: {path}")
        found.append(path)
    if not found:
        raise RuntimeError(f"loaded Engine image is absent from process mappings: {path}")
    return path


def _asl_path(core, build):
    path = Path(os.environ.get("EPCSAFT_ASL_LIBRARY", core.parent / "libepcsaft_asl.so")).expanduser().resolve(strict=True)
    if not (path.is_relative_to(core.parent) or build and path.is_relative_to(build)):
        raise RuntimeError(f"selected ASL library is outside the loaded install or Meson build: {path}")
    return path


def loaded_identity():
    import epcsaft
    from epcsaft import _core

    core = Path(_core.__file__).resolve(strict=True)
    package = Path(epcsaft.__file__).resolve(strict=True) if epcsaft.__file__ else None
    try:
        direct = json.loads(metadata.distribution("epcsaft").read_text("direct_url.json") or "{}")
        version = metadata.version("epcsaft")
    except metadata.PackageNotFoundError:
        direct, version = {}, None
    build = _build_dir(core)
    build_source = _build_source(build) if build else None
    source = _source_root(direct) or build_source
    if build_source and source != build_source:
        raise RuntimeError("selected Meson build was configured from a different source tree")
    revision, source_diff = _source_identity(source)
    build, build_hash, lock_hash = _build_identity(build, source)
    mappings = Path("/proc/self/maps").read_text().splitlines()
    core = _mapped(core, mappings)
    native = {Path(line.split(maxsplit=5)[5]).resolve() for line in mappings
              if len(line.split(maxsplit=5)) == 6 and line.split(maxsplit=5)[5].endswith("/libepcsaft_native.so")}
    if len(native) > 1:
        raise RuntimeError(f"multiple loaded native libraries: {native}")
    numerical = _mapped(next(iter(native), core), mappings)
    asl = _asl_path(core, build)
    if build:
        _verify_products(core, numerical, asl, build)
    options = build / "meson-info/intro-buildoptions.json" if build else None
    if version is None and build:
        version = json.loads((build / "meson-info/intro-projectinfo.json").read_text()).get("version")
    return {"prepared_stage": os.environ.get("EPCSAFT_PREPARED_STAGE", "editable" if direct else "prepared"),
            "interpreter": sys.executable, "python": sys.version, "version": version, "direct_url": direct,
            "wheel": unquote(urlparse(direct.get("url", "")).path) if direct.get("url", "").endswith(".whl") else None,
            "wheel_sha256": digest(unquote(urlparse(direct.get("url", "")).path))
            if direct.get("url", "").endswith(".whl") and Path(unquote(urlparse(direct["url"]).path)).is_file() else None,
            "source_root": str(source) if source else None, "source_revision": revision,
            "source_diff_sha256": source_diff,
            "lock_sha256": lock_hash, "meson_build_dir": str(build) if build else None,
            "meson_configuration_sha256": build_hash,
            "build_configuration": str(options) if options and options.is_file() else None,
            "build_configuration_sha256": digest(options) if options and options.is_file() else None,
            "epcsaft": str(package) if package else None,
            "epcsaft_sha256": digest(package) if package and package.is_file() else None,
            "extension": str(core), "extension_sha256": digest(core),
            "native_library": str(numerical), "native_library_sha256": digest(numerical),
            "asl_library": str(asl) if asl else None, "asl_library_sha256": digest(asl) if asl else None}


def asl_library(identity=None):
    identity = identity or loaded_identity()
    path = identity.get("asl_library")
    if not path:
        raise RuntimeError("the selected Engine environment has no built ASL library")
    path = Path(path).resolve(strict=True)
    if digest(path) != identity.get("asl_library_sha256"):
        raise RuntimeError(f"selected ASL library changed after identity capture: {path}")
    return path


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
