"""Interrupt real numerical work and check that its checkpoint and cleanup survive."""

import argparse
import csv
import json
import os
import signal
import subprocess
import sys
import time
from pathlib import Path


def check_interruption(signum: signal.Signals, output: Path, ipopt: Path) -> dict:
    output.mkdir(parents=True)
    campaign = output / "campaign"
    runner = Path(__file__).with_name("run.py").resolve()
    child_group = None
    command = [sys.executable, "-u", str(runner), "--output-dir", str(campaign),
               "--ipopt", str(ipopt)]
    with (output / "controller.stdout.log").open("w") as stdout, \
            (output / "controller.stderr.log").open("w") as stderr:
        parent = subprocess.Popen(command, stdout=stdout, stderr=stderr, start_new_session=True)
        try:
            deadline = time.monotonic() + 20.0
            target = campaign / "derivatives.csv"
            while time.monotonic() < deadline and parent.poll() is None:
                if target.is_file():
                    with target.open(newline="") as stream:
                        rows = list(csv.DictReader(stream))
                    children = Path(f"/proc/{parent.pid}/task/{parent.pid}/children").read_text().split()
                    for pid in children:
                        try:
                            arguments = Path(f"/proc/{pid}/cmdline").read_bytes().split(b"\0")
                        except FileNotFoundError:
                            continue
                        if rows and b"--numerical-stage" in arguments and b"initial" in arguments:
                            child_group = int(pid)
                            os.killpg(child_group, signal.SIGSTOP)
                            break
                    if child_group is not None:
                        break
                time.sleep(0.002)
            if child_group is None:
                raise RuntimeError("no scalar checkpoint observed in an active numerical child")
            while time.monotonic() < deadline:
                state = next(line for line in Path(f"/proc/{child_group}/status").read_text().splitlines()
                             if line.startswith("State:"))
                if "T (stopped)" in state:
                    break
                time.sleep(0.002)
            else:
                raise RuntimeError("numerical child did not stop at its scalar checkpoint")
            with target.open(newline="") as stream:
                before = list(csv.DictReader(stream))
            os.kill(parent.pid, signum)
            parent.wait(timeout=5.0)
            try:
                os.killpg(child_group, 0)
                alive = True
            except ProcessLookupError:
                alive = False
            with target.open(newline="") as stream:
                after = list(csv.DictReader(stream))
            summary = json.loads((campaign / "summary.json").read_text())
            identity_keys = (
                "source_commit", "wheel_sha256", "asl_library_sha256",
                "packet_fingerprint", "ipopt_sha256",
            )
            result = {
                "signal": signum.name, "parent_returncode": parent.returncode,
                "expected_parent_returncode": -signal.SIGINT,
                "exit_semantics": (
                    "SIGTERM handler raises KeyboardInterrupt; CPython reports -SIGINT"
                    if signum == signal.SIGTERM else
                    "SIGINT raises KeyboardInterrupt; CPython reports -SIGINT"
                ),
                "child_group": child_group, "child_group_alive": alive,
                "partial_numerical_rows": len(before), "partial_rows_preserved": before == after,
                "summary_status": summary["status"], "summary_stage": summary.get("stage"),
                "identity_preserved": all(summary.get(key) for key in identity_keys),
                "injection": "real numerical child paused after an atomic scalar checkpoint",
            }
            (output / "probe.json").write_text(json.dumps(result, indent=2) + "\n")
            assert before and before == after, "partial numerical results changed during interruption"
            assert parent.returncode == -signal.SIGINT, "unexpected controller exit semantics"
            assert not alive, "campaign interruption left a child group alive"
            assert summary["status"] == "failed attempt" and summary["stage"] == "interrupted"
            assert result["identity_preserved"], "interrupted summary lost campaign identity"
            return result
        finally:
            if child_group is not None:
                try:
                    os.killpg(child_group, signal.SIGKILL)
                except ProcessLookupError:
                    pass
            if parent.poll() is None:
                os.killpg(parent.pid, signal.SIGKILL)
                parent.wait()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True,
                        help="new output directory; existing evidence is preserved")
    parser.add_argument("--ipopt", type=Path, required=True)
    args = parser.parse_args()
    signal.signal(signal.SIGTERM, signal.default_int_handler)
    for signum in (signal.SIGTERM, signal.SIGINT):
        print(json.dumps(check_interruption(
            signum, args.output_dir.resolve() / signum.name.lower(), args.ipopt.resolve(),
        )), flush=True)
