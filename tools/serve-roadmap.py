#!/usr/bin/env python3
"""Serve the local Quarto companion and one coherent test snapshot."""

from __future__ import annotations

import argparse
import hashlib
import html as html_module
import json
import logging
import mimetypes
import os
import re
import runpy
import signal
import subprocess
import sys
import threading
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, unquote, urlsplit

import markdown

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
_worker_enabled = False
_worker_process: subprocess.Popen[bytes] | None = None
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


def _read_worker(store: RecordStore) -> dict:
    if not _worker_enabled:
        return {"state": "Idle", "automatic_execution": "off", "reason": "view-only mode"}
    if _worker_process is None or _worker_process.poll() is not None:
        return {
            "state": "Unavailable",
            "automatic_execution": "off",
            "reason": "automatic worker process exited",
        }
    try:
        value = json.loads(store.worker_path.read_text(encoding="utf-8"))
        return value if isinstance(value, dict) else {
            "state": "Unavailable",
            "automatic_execution": "off",
            "reason": "worker state has invalid shape",
        }
    except OSError:
        return {
            "state": "Unavailable",
            "automatic_execution": "off",
            "reason": "worker state is unavailable",
        }
    except json.JSONDecodeError:
        return {
            "state": "Unavailable",
            "automatic_execution": "off",
            "reason": "worker state is malformed",
        }


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
        worker = _read_worker(store)
        binding = _inventory_binding_signature()
        try:
            inventory = _current_inventory(generation=before, binding=binding)
            # A manifest can contain many cases sharing one wheel/native runtime.
            # Attest each distinct artifact identity once per coherent snapshot.
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
                "studies": [], "worker": {"state": "Unavailable", "reason": "automatic checks unavailable"},
                "automatic_execution": "off",
            }
        projection = COMPANION_DATA["companion_projection"](
            plan,
            inventory,
            runs,
            worker=worker,
            automatic_execution=str(worker.get("automatic_execution", "off")),
            revision_inputs={
                "plan_sha256": plan_hash,
                "input_generation": after,
                "inventory_generation": inventory.get("inventory_generation"),
                "manifest_revision": manifest.get("revision"),
                "worker_revision": worker.get("revision"),
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


def render(source: str) -> str:
    if ROADMAP["START"] in source:
        before, tail = source.split(ROADMAP["START"])
        _, after = tail.split(ROADMAP["END"])
        summary = ROADMAP["summarize"](source, PLAN.parent)
        source = before + ROADMAP["START"] + "\n" + summary + "\n" + ROADMAP["END"] + after
    formulas: list[str] = []

    def preserve(match: re.Match[str]) -> str:
        formulas.append(html_module.escape(match.group()))
        return f"MATHPLACEHOLDER{len(formulas) - 1}END"

    source = re.sub(r"\$\$[\s\S]+?\$\$|(?<!\$)\$[^\n$]+\$(?!\$)", preserve, source)
    html = markdown.markdown(source, extensions=["tables", "fenced_code", "toc"])
    for index, formula in enumerate(formulas):
        html = html.replace(f"MATHPLACEHOLDER{index}END", formula)
    return re.sub(r'href="(?![a-zA-Z]+:|/|#)([^\"]+)"', r'href="/docs/scientific/\1"', html)


def coverage() -> dict:
    plan = PLAN.read_text(encoding="utf-8")
    return {
        "requirements": [],
        "inventory": _current_inventory(),
        "plan": COMPANION_DATA["plan_projection"](plan, tuple(ROADMAP["LAYERS"])),
    }


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
        elif route == "/coverage":
            body, content_type = json.dumps(coverage(), allow_nan=False), "application/json"
        elif route == "/artifact":
            relative = parse_qs(urlsplit(self.path).query).get("path", [""])[0]
            try:
                body = json.dumps(INVENTORY["artifact"](relative), allow_nan=False)
            except (ValueError, OSError) as error:
                self.send_error(422, str(error))
                return
            content_type = "application/json"
        elif route in ("/render", "/source"):
            body = PLAN.read_text(encoding="utf-8")
            if route == "/render":
                try:
                    body = render(body)
                except ValueError as error:
                    self.send_error(422, f"invalid roadmap: {error}")
                    return
            content_type = "text/plain" if route == "/source" else "text/html"
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
    sample = render("# Heading\n\n| A | B |\n|---|---|\n| 1 | 2 |\n\n[x](record.md)\n\n```mermaid\nflowchart LR\n A-->B\n```")
    assert "<table>" in sample and "language-mermaid" in sample
    assert 'href="/docs/scientific/record.md"' in sample
    assert render("# Changed") != render("# Heading")
    projection = companion_projection()
    assert projection["inventory"]["state"] in {"Complete", "Incomplete", "Unavailable"}
    assert projection["declarations"]
    print("roadmap and companion self-test passed")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--port", type=int, default=8765)
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--view-only", action="store_true",
                        help="serve the companion without starting automatic checks")
    args = parser.parse_args()
    if args.self_test:
        self_test()
    else:
        _worker_enabled = not args.view_only
        worker: subprocess.Popen[bytes] | None = None
        if _worker_enabled:
            worker_python = ROOT / ".venv/bin/python"
            command = [str(worker_python if worker_python.is_file() else Path(sys.executable)),
                       str(ROOT / "tools/test-worker.py")]
            worker = subprocess.Popen(command, cwd=ROOT, start_new_session=True,
                                      stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
        _worker_process = worker
        ensure_companion()
        server = ThreadingHTTPServer(("127.0.0.1", args.port), Handler)
        server.timeout = 1.0
        mode = "view-only" if args.view_only else "viewer + automatic worker"
        print(f"Live companion ({mode}): http://127.0.0.1:{server.server_port}", flush=True)
        try:
            while True:
                server.handle_request()
                if worker is not None and worker.poll() is not None:
                    raise RuntimeError(f"automatic worker exited with code {worker.returncode}")
        except KeyboardInterrupt:
            pass
        finally:
            server.server_close()
            if worker is not None and worker.poll() is None:
                try:
                    os.killpg(worker.pid, signal.SIGTERM)
                    worker.wait(timeout=3)
                except (OSError, ProcessLookupError, subprocess.TimeoutExpired):
                    try:
                        os.killpg(worker.pid, signal.SIGKILL)
                    except (OSError, ProcessLookupError):
                        worker.kill()
