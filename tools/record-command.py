#!/usr/bin/env python3
"""Run one repository command and publish its outcome through the recorder."""

from __future__ import annotations

import argparse
import contextlib
import hashlib
import os
import signal
import subprocess
import sys
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tools.test_records import (
    ATTEMPT_ACTIVE_ENV,
    RESULTS,
    RecordStore,
    input_generation,
    locked,
    now,
    validate_run,
)


def _run_id(test_id: str) -> str:
    stamp = datetime.now(UTC).strftime("%Y%m%dT%H%M%S%fZ")
    safe = "".join(char if char.isalnum() or char in "_.-" else "_" for char in test_id)
    return f"{stamp}-{safe}-single"


def _descendant_pids(root_pid: int) -> set[int]:
    parents: dict[int, int] = {}
    for entry in Path("/proc").glob("[0-9]*"):
        try:
            stat = (entry / "stat").read_text(encoding="utf-8")
            parents[int(entry.name)] = int(stat.rsplit(")", 1)[1].split()[1])
        except (OSError, ValueError, IndexError):
            continue
    pids = {root_pid}
    changed = True
    while changed:
        changed = False
        for pid, parent in parents.items():
            if parent in pids and pid not in pids:
                pids.add(pid)
                changed = True
    return pids


def _stop(process: subprocess.Popen[Any]) -> None:
    pids = _descendant_pids(process.pid)
    for pid in pids:
        try:
            os.kill(pid, signal.SIGTERM)
        except (OSError, ProcessLookupError):
            pass
    try:
        process.wait(timeout=2)
    except subprocess.TimeoutExpired:
        for pid in pids | _descendant_pids(process.pid):
            try:
                os.kill(pid, signal.SIGKILL)
            except (OSError, ProcessLookupError):
                pass
        process.wait()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--test-id", required=True)
    parser.add_argument("--timeout", type=int, default=1800)
    parser.add_argument("--cwd", type=Path, default=ROOT)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = list(args.command)
    if command[:1] == ["--"]:
        command.pop(0)
    if not command:
        parser.error("a command is required after --")

    RESULTS.mkdir(parents=True, exist_ok=True)
    store = RecordStore(ROOT)
    context = contextlib.nullcontext() if os.environ.get("EPCSAFT_RECORD_LOCK_HELD") == "1" else locked(store.results / ".lock")
    with context:
        store.reconcile_attempts()
        input_id = input_generation(ROOT)
        run_id = _run_id(args.test_id)
        artifact = {"state": "Unavailable", "reason": "command does not produce an Engine artifact"}
        started = now()
        base = {
            "schema_version": 1,
            "run_id": run_id,
            "test_id": args.test_id,
            "parameter_case_id": "single",
            "input_generation": input_id,
            "artifact_identity": artifact,
            "started_at": started,
            "finished_at": now(),
            "outcome_reason": None,
            "comparisons": [],
            "artifact_digests": [],
            "command": command,
        }
        store.publish({**base, "execution": "Queued"})
        store.publish({**base, "execution": "Running"})
        output = RESULTS / "runs" / "pending" / f"{run_id}.log"
        output.parent.mkdir(parents=True, exist_ok=True)
        process: subprocess.Popen[Any] | None = None
        execution = "Failed"
        exit_code: int | None = None
        reason: str | None = None
        interruption_kind: str | None = None
        try:
            with output.open("w", encoding="utf-8") as stream:
                process = subprocess.Popen(
                    command, cwd=args.cwd, stdout=stream, stderr=subprocess.STDOUT,
                    start_new_session=True,
                    env=dict(
                        os.environ,
                        EPCSAFT_NATIVE_RECORDING="1",
                        EPCSAFT_RECORD_LOCK_HELD="1",
                        **{ATTEMPT_ACTIVE_ENV: "1"},
                    ),
                )
                try:
                    exit_code = process.wait(timeout=args.timeout)
                    execution = "Passed" if exit_code == 0 else "Failed"
                    reason = None if execution == "Passed" else f"process exit={exit_code}"
                except subprocess.TimeoutExpired:
                    _stop(process)
                    execution, reason, exit_code = "Interrupted", "command timed out", None
                    interruption_kind = "timeout"
        except KeyboardInterrupt:
            if process is not None:
                _stop(process)
            execution, reason, exit_code = "Interrupted", "command interrupted", None
            interruption_kind = "interruption"
        except OSError as error:
            execution, reason, exit_code = "Error", str(error), None
        run = {
            **base,
            "execution": execution,
            "finished_at": now(),
            "outcome_reason": reason,
            "exit_code": exit_code,
            "output": str(output.relative_to(ROOT)) if output.is_relative_to(ROOT) else str(output),
            "output_sha256": hashlib.sha256(output.read_bytes()).hexdigest() if output.is_file() else None,
        }
        if interruption_kind is not None:
            run["interruption_kind"] = interruption_kind
        errors = validate_run(run)
        if errors:
            run["execution"] = "Error"
            run.pop("interruption_kind", None)
            run["outcome_reason"] = "; ".join(errors)
        store.publish(run)
        print(f"{args.test_id}: {run['execution']} ({run_id})", flush=True)
        return 0 if run["execution"] == "Passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
