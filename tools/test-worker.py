#!/usr/bin/env python3
"""Run permitted automatic checks for each debounced input generation."""

from __future__ import annotations

import argparse
import contextlib
import fcntl
import os
import signal
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tools.test_records import (
    RESULTS,
    RecordStore,
    atomic_json,
    digest,
    input_generation,
    locked,
    now,
)

WORKER_LOCK = RESULTS / "worker.lock"
WORKER_OUTPUT_NAME = "worker.log"
RUNNER_BUILD_TIMEOUT = 1800
AUTOMATIC_CASE_BUDGET = 32


@contextlib.contextmanager
def _worker_lock(*, blocking: bool = True):
    """Allow one automatic worker process to own the watch loop."""

    RESULTS.mkdir(parents=True, exist_ok=True)
    with WORKER_LOCK.open("a+") as stream:
        flags = fcntl.LOCK_EX | (0 if blocking else fcntl.LOCK_NB)
        fcntl.flock(stream.fileno(), flags)
        try:
            yield stream
        finally:
            fcntl.flock(stream.fileno(), fcntl.LOCK_UN)


def _stop_process_group(process: subprocess.Popen[bytes], *, grace: float = 2.0) -> None:
    if process.poll() is not None:
        return
    try:
        os.killpg(process.pid, signal.SIGTERM)
    except (OSError, ProcessLookupError):
        process.terminate()
    try:
        process.wait(timeout=grace)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except (OSError, ProcessLookupError):
            process.kill()
        process.wait()


def _write_state(**fields: Any) -> dict[str, Any]:
    current = {
        "schema_version": 1,
        "updated_at": now(),
        "automatic_execution": "on",
        **fields,
    }
    current["revision"] = digest(current)
    atomic_json(RESULTS / "worker.json", current)
    return current


def _generation() -> str:
    return input_generation(ROOT)


def _read_worker_output() -> str:
    try:
        return (RESULTS / WORKER_OUTPUT_NAME).read_bytes()[-4000:].decode(errors="replace")
    except OSError:
        return ""


def _reconcile_runner_attempts(generation: str) -> bool:
    """Close a dead runner's attempts and report incomplete current work."""

    store = RecordStore(ROOT)
    with locked(store.results / ".lock"):
        store.reconcile_attempts(owner=True)
        return any(
            run.get("input_generation") == generation
            and (
                run.get("execution") in {"Queued", "Running"}
                or (
                    run.get("execution") == "Interrupted"
                    and run.get("interruption_kind") != "timeout"
                )
            )
            for run in store.read_runs()
        )


def _run_once_owned(
    *, timeout: int, generation: str | None = None, _retry_infrastructure: bool = True,
) -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    target = generation or _generation()
    _write_state(state="Queued", requested_generation=target)
    command = [
        str(ROOT / ".venv/bin/python" if (ROOT / ".venv/bin/python").is_file() else sys.executable),
        "tools/run-numerical-evidence.py",
        "--automatic", "--timeout", str(timeout), "--worker",
    ]
    _write_state(state="Running", requested_generation=target, started_at=now())
    process: subprocess.Popen[bytes] | None = None
    try:
        with (RESULTS / WORKER_OUTPUT_NAME).open("wb") as stream:
            process = subprocess.Popen(command, cwd=ROOT, start_new_session=True,
                                       stdout=stream, stderr=subprocess.STDOUT)
            try:
                # ``run-numerical-evidence`` may spend its full build timeout
                # before executing the automatic cases.  The worker timeout must
                # cover that build plus the bounded per-case budget; ``timeout *
                # 2`` could kill a normal run at 240 seconds before preparation
                # finishes.
                process.wait(timeout=RUNNER_BUILD_TIMEOUT + timeout * AUTOMATIC_CASE_BUDGET)
            except subprocess.TimeoutExpired:
                _stop_process_group(process)
                _reconcile_runner_attempts(target)
                _write_state(state="Failed", requested_generation=target, last_exit_code=None,
                             last_completed_generation=target,
                             last_output="automatic worker timed out", finished_at=now())
                return 124
        output = _read_worker_output()
        final_generation = _generation()
    except KeyboardInterrupt:
        if process is not None:
            _stop_process_group(process)
        _reconcile_runner_attempts(target)
        _write_state(state="Interrupted", requested_generation=target, last_exit_code=None,
                     last_output="automatic worker interrupted", finished_at=now())
        return 130
    except (OSError, RuntimeError) as error:
        _write_state(state="Unavailable", requested_generation=target, last_exit_code=None,
                     last_output=str(error), finished_at=now())
        return 1
    incomplete = _reconcile_runner_attempts(target)
    if incomplete and final_generation == target and _retry_infrastructure:
        return _run_once_owned(
            timeout=timeout, generation=target, _retry_infrastructure=False,
        )
    if final_generation != target:
        _write_state(state="Outdated", requested_generation=target,
                     last_completed_generation=target, last_exit_code=process.returncode,
                     last_output=output, finished_at=now())
    else:
        exit_code = process.returncode if process.returncode else (1 if incomplete else 0)
        last_output = output
        if incomplete:
            last_output = (
                "automatic runner left an incomplete attempt after its bounded retry; "
                + (output or "manual intervention required")
            )
        _write_state(state="Idle" if exit_code == 0 else "Failed",
                     requested_generation=target, last_completed_generation=target,
                     last_exit_code=exit_code, last_output=last_output, finished_at=now())
    return exit_code


def run_once(*, timeout: int, generation: str | None = None) -> int:
    try:
        with _worker_lock(blocking=False):
            return _run_once_owned(timeout=timeout, generation=generation)
    except BlockingIOError:
        return 0


def watch(*, timeout: int, interval: float, debounce: float) -> int:
    try:
        worker_lock = _worker_lock(blocking=False)
        worker_lock.__enter__()
    except BlockingIOError:
        return 0
    try:
        RESULTS.mkdir(parents=True, exist_ok=True)
        last_requested: str | None = None
        pending_since: float | None = None
        # Preserve the terminal generation across worker restarts.  This is
        # the durable coalescing key; an unchanged failure is not retried.
        try:
            previous = __import__("json").loads((RESULTS / "worker.json").read_text())
        except (OSError, ValueError):
            previous = {}
        last_completed = previous.get("last_completed_generation") if isinstance(previous, dict) else None
        previous_terminal = (
            isinstance(previous, dict)
            and previous.get("state") not in {None, "Queued", "Running"}
        )
        if not previous_terminal:
            _write_state(state="Idle", last_completed_generation=last_completed)
        try:
            while True:
                try:
                    generation = _generation()
                except (OSError, RuntimeError) as error:
                    _write_state(state="Unavailable", last_output=str(error))
                    time.sleep(interval)
                    continue
                if generation != last_requested:
                    last_requested = generation
                    pending_since = time.monotonic()
                    if generation != last_completed:
                        _write_state(state="Queued", requested_generation=generation,
                                     last_completed_generation=last_completed)
                if pending_since is not None and time.monotonic() - pending_since >= debounce:
                    if generation != last_completed:
                        # run_once is already serialized by the worker lock;
                        # call the private body without reacquiring it.
                        code = _run_once_owned(timeout=timeout, generation=generation)
                        last_completed = generation
                        if code and _generation() != generation:
                            last_requested = None
                        pending_since = None
                    else:
                        pending_since = None
                time.sleep(interval)
        except KeyboardInterrupt:
            _write_state(state="Idle", interrupted=True,
                         last_completed_generation=last_completed)
            return 130
    finally:
        worker_lock.__exit__(None, None, None)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--once", action="store_true", help="run one automatic generation and exit")
    parser.add_argument("--timeout", type=int, default=120)
    parser.add_argument("--interval", type=float, default=2.0)
    parser.add_argument("--debounce", type=float, default=0.5)
    args = parser.parse_args()
    if args.once:
        return run_once(timeout=args.timeout)
    return watch(timeout=args.timeout, interval=args.interval, debounce=args.debounce)


if __name__ == "__main__":
    raise SystemExit(main())
