#!/usr/bin/env python3
"""Serve the local Quarto companion and one coherent read-only test snapshot."""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import mimetypes
import runpy
import subprocess
import sys
import threading
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import unquote, urlsplit

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tools.test_records import (
    RecordStore,
    digest,
    display_status,
    input_generation,
    sha256,
)

PLAN = ROOT / "docs/scientific/implementation-plan.md"
COMPANION_SOURCE = ROOT / "docs/companion"
COMPANION_BUILD = ROOT / "build/companion"
ROADMAP = runpy.run_path(str(ROOT / "tools/check-roadmap.py"))
INVENTORY = runpy.run_path(str(ROOT / "tools/numerical-inventory.py"))
COMPANION_DATA = runpy.run_path(str(ROOT / "tools/companion-data.py"))

_COMPANION_LOCK = threading.Lock()
_SNAPSHOT_LOCK = threading.Lock()
_INVENTORY_LOCK = threading.Lock()
_companion_signature: tuple[tuple[str, int, int], ...] | None = None
_companion_error: Exception | None = None
_companion_revision: str | None = None
_inventory_cache: tuple[str, str, dict] | None = None
_LOGGER = logging.getLogger(__name__)


def _companion_sources() -> tuple[tuple[str, int, int], ...]:
    if not COMPANION_SOURCE.is_dir():
        return ()
    suffixes = {".qmd", ".yml", ".yaml", ".css", ".scss", ".js", ".html"}
    files = []
    for path in COMPANION_SOURCE.rglob("*"):
        if any(part in {".quarto", "build", "builds"} for part in path.parts):
            continue
        if path.is_file() and path.suffix.lower() in suffixes:
            stat = path.stat()
            files.append((str(path.relative_to(COMPANION_SOURCE)), stat.st_mtime_ns, stat.st_size))
    return tuple(sorted(files))


def ensure_companion() -> None:
    """Render the static shell once per source revision, without execution."""

    global _companion_signature, _companion_error, _companion_revision
    signature = _companion_sources()
    if not signature:
        _companion_error = OSError("companion sources are unavailable")
        _companion_revision = None
        return
    if signature == _companion_signature:
        return
    with _COMPANION_LOCK:
        signature = _companion_sources()
        if not signature or signature == _companion_signature:
            return
        try:
            result = subprocess.run(
                ["quarto", "render", "docs/companion", "--no-execute"],
                cwd=ROOT, check=False, timeout=60, capture_output=True, text=True,
            )
            if result.returncode != 0:
                raise OSError(result.stderr[-2000:] or "quarto render failed")
            if not (COMPANION_BUILD / "index.html").is_file():
                raise OSError("Quarto did not produce build/companion/index.html")
            _companion_error = None
            _companion_revision = hashlib.sha256(repr(signature).encode()).hexdigest()[:16]
        except (OSError, subprocess.TimeoutExpired) as error:
            _companion_error = error
            _companion_revision = None
        _companion_signature = signature


def _inventory_binding_signature() -> str:
    """Identify the prepared native binaries used by inventory discovery."""

    binding = ROOT / "build/environment-wheel/preparation.json"
    marker = ROOT / "build/environment-wheel/meson-build-dir"
    entries: list[tuple[str, str | None]] = []
    for path in (binding, marker):
        entries.append((str(path.relative_to(ROOT)), sha256(path) if path.is_file() else None))
    build_dir: Path | None = None
    if marker.is_file():
        try:
            candidate = Path(marker.read_text(encoding="utf-8").strip()).resolve()
            if candidate.is_dir() and candidate.is_relative_to(ROOT.resolve()):
                build_dir = candidate
        except (OSError, ValueError):
            build_dir = None
    targets = INVENTORY.get("NATIVE_TARGETS", {})
    if isinstance(targets, dict) and build_dir is not None:
        for target in sorted(str(value) for value in targets.values()):
            path = build_dir / target
            entries.append((str(path.relative_to(ROOT)), sha256(path) if path.is_file() else None))
    return digest(entries)


def _current_inventory(*, generation: str | None = None, binding: str | None = None) -> dict:
    global _inventory_cache
    generation = generation if generation is not None else input_generation(ROOT)
    binding = binding if binding is not None else _inventory_binding_signature()
    with _INVENTORY_LOCK:
        if _inventory_cache is not None and _inventory_cache[:2] == (generation, binding):
            return _inventory_cache[2]
        value = INVENTORY["inventory"](write=True, generation=generation)
        # A fallback inventory is a transient view of an incomplete build. Do
        # not retain it under a source-only key; the next prepared binding must
        # be allowed to publish its real native declarations.
        if value.get("state") == "Complete":
            _inventory_cache = (generation, binding, value)
        else:
            _inventory_cache = None
        return value


def companion_projection() -> dict:
    """Assemble plan, inventory, applicability and records coherently."""

    ensure_companion()
    plan = PLAN.read_text(encoding="utf-8")
    ROADMAP["summarize"](plan, PLAN.parent)
    store = RecordStore(ROOT)
    plan_hash = hashlib.sha256(plan.encode()).hexdigest()
    with _SNAPSHOT_LOCK:
        before = input_generation(ROOT)
        manifest = store.read_manifest()
        binding = _inventory_binding_signature()
        try:
            inventory = _current_inventory(generation=before, binding=binding)
            # Many cases share one wheel/native runtime; check each distinct
            # identity once per snapshot.  This is display freshness (sources
            # and artifact files unchanged), not environment verification.
            freshness_cache: dict[str, tuple[str, str | None]] = {}
            runs = []
            for run in store.read_runs(manifest):
                key = digest({
                    "artifact_identity": run.get("artifact_identity", {"state": "Unavailable"}),
                    "input_generation": run.get("input_generation"),
                })
                if key not in freshness_cache:
                    checked = display_status(run, ROOT, current_generation=before)
                    freshness_cache[key] = (checked.get("freshness", "Unverified"), checked.get("freshness_reason"))
                status, reason = freshness_cache[key]
                checked = dict(run)
                checked["freshness"] = status
                if reason:
                    checked["freshness_reason"] = reason
                runs.append(checked)
            after = input_generation(ROOT)
        except (OSError, ValueError, RuntimeError) as error:
            _LOGGER.debug("companion snapshot unavailable", exc_info=error)
            return {
                "schema_version": 1,
                "state": "Offline",
                "revision": None,
                "presentation_available": False,
                "inventory": {"state": "Unavailable", "errors": ["inventory unavailable"]},
                "declarations": [], "runs": [], "plan": [], "preparation": [],
                "studies": [],
            }
        projection = COMPANION_DATA["companion_projection"](
            plan,
            inventory,
            runs,
            revision_inputs={
                "plan_sha256": plan_hash,
                "input_generation": after,
                "inventory_generation": inventory.get("inventory_generation"),
                "manifest_revision": manifest.get("revision"),
                "runs": [
                    {
                        "run_id": run.get("run_id"),
                        "test_id": run.get("test_id"),
                        "parameter_case_id": run.get("parameter_case_id"),
                        "input_generation": run.get("input_generation"),
                        "freshness": run.get("freshness"),
                        "execution": run.get("execution"),
                        "evidence": [
                            (row.get("comparison_id"), row.get("judgment"), row.get("evidence_state"))
                            for row in run.get("comparisons", [])
                            if isinstance(row, dict)
                        ],
                    }
                    for run in runs
                ],
            },
            layers=tuple(ROADMAP["LAYERS"]),
        )
        projection["state"] = "Synced" if before == after else "Updating"
        projection["presentation_available"] = _companion_error is None and (COMPANION_BUILD / "index.html").is_file()
        if not projection["presentation_available"]:
            projection["state"] = "Offline"
        if projection["state"] != "Synced" or inventory.get("state") != "Complete":
            projection["revision"] = None
        return projection


def _companion_asset(route: str) -> Path | None:
    ensure_companion()
    if _companion_error is not None:
        return None
    relative = "index.html" if route == "/" else route.lstrip("/")
    path = (COMPANION_BUILD / relative).resolve()
    return path if path.is_relative_to(COMPANION_BUILD) and path.is_file() else None


def _content_type(path: Path) -> str:
    return mimetypes.guess_type(str(path))[0] or "application/octet-stream"


class Handler(BaseHTTPRequestHandler):
    def do_GET(self) -> None:
        route = unquote(urlsplit(self.path).path)
        if route == "/companion":
            try:
                body = json.dumps(companion_projection(), allow_nan=False)
            except (OSError, ValueError, RuntimeError) as error:
                self.send_error(503, str(error))
                return
            content_type = "application/json"
        else:
            path = _companion_asset(route)
            site_asset = path is not None
            if path is None:
                path = (ROOT / route.lstrip("/")).resolve()
            valid = path.is_file()
            if not site_asset:
                valid = valid and path.is_relative_to(ROOT) and path.suffix in {
                    ".md", ".csv", ".tsv", ".json", ".log", ".cpp", ".hpp", ".py", ".xml",
                    ".txt", ".html", ".pdf", ".css", ".js", ".map", ".svg", ".png", ".jpg",
                    ".jpeg", ".gif", ".woff", ".woff2", ".ttf",
                }
            if not valid:
                self.send_error(503 if route == "/" else 404, "companion build is unavailable" if route == "/" else None)
                return
            body, content_type = path.read_bytes(), _content_type(path)
        payload = body.encode() if isinstance(body, str) else body
        self.send_response(200)
        self.send_header("Content-Type", content_type + "; charset=utf-8")
        self.send_header("Cache-Control", "no-store")
        self.send_header("Content-Length", str(len(payload)))
        self.end_headers()
        self.wfile.write(payload)

    def log_message(self, *_: object) -> None:
        pass


def self_test() -> None:
    projection = companion_projection()
    assert projection["inventory"]["state"] in {"Complete", "Incomplete", "Unavailable"}
    assert projection["declarations"]
    print("companion self-test passed")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--port", type=int, default=8765)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
    else:
        ensure_companion()
        server = ThreadingHTTPServer(("127.0.0.1", args.port), Handler)
        print(f"Live companion: http://127.0.0.1:{server.server_port}", flush=True)
        try:
            server.serve_forever()
        except KeyboardInterrupt:
            pass
        finally:
            server.server_close()
