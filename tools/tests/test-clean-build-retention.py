"""#1140 fixture: a closed, uncited build/issue-N holding a file an analysis records by SHA-256 that Git does not hold.

Usage: python3 tools/tests/test-clean-build-retention.py [BASE]
Everything runs in a throwaway Git repository with a stub host reader (no open items; every issue closed ten days
ago, except 9007 in an unrecognized state, 9008 with no state and 9009 open). With BASE (a commit whose tools/clean-build
predates #1140), first shows that tool would expire the folder.
"""
import gzip
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile

TOOLS = Path(__file__).resolve().parents[1]
STUB = '''from datetime import datetime, timedelta, timezone
class HostState:
    def __init__(self, root): pass
    def open_items(self, kind): return []
    def issue(self, number):
        old = (datetime.now(timezone.utc) - timedelta(days=10)).isoformat()
        return {"9007": {"state": "unknown", "closed_at": old}, "9008": {"closed_at": old},
                "9009": {"state": "open"}}.get(
            str(number), {"state": "closed", "closed_at": old})
'''
# Today's tool prints nothing per folder, so print each folder it would hand to keep_reports.
BEFORE = '''import importlib.machinery, importlib.util, sys
sys.path.insert(0, sys.argv[1].rsplit('/', 1)[0])
loader = importlib.machinery.SourceFileLoader('cb', sys.argv[1])
cb = importlib.util.module_from_spec(importlib.util.spec_from_loader('cb', loader)); loader.exec_module(cb)
keep = cb.keep_reports
cb.keep_reports = lambda folder, dry: (print(f'EXPIRE {folder.name}'), keep(folder, dry))[1] if dry else sys.exit('dry only')
sys.argv = [sys.argv[1], '--dry-run']; cb.main()
'''


def git(root, *args):
    return subprocess.run(["git", "-c", "user.name=fixture", "-c", "user.email=fixture@example.invalid",
                           "-c", "core.hooksPath=/dev/null", "-c", "commit.gpgsign=false", *args],
                          cwd=root, check=True, capture_output=True, text=True).stdout


def record(root, name, data):
    """Write the producer under build/ and an analysis record naming its SHA-256."""
    (root / name).parent.mkdir(parents=True, exist_ok=True)
    (root / name).write_bytes(data)
    leaf = root / "analyses/demo/results" / f"{Path(name).parent.name}.json"
    leaf.parent.mkdir(parents=True, exist_ok=True)
    leaf.write_text(json.dumps({"producer_sha256": hashlib.sha256(data).hexdigest()}))


def states(root):
    """Each folder's printed state from a dry run of the fixture's copy of the tool."""
    run = subprocess.run([sys.executable, str(root / "tools/clean-build"), "--dry-run"],
                         cwd=root, check=True, capture_output=True, text=True).stdout
    return {folder: state for state, folder in re.findall(r"^(\w+) build/(issue-\d+)", run, re.M)}


with tempfile.TemporaryDirectory() as directory:
    root = Path(directory).resolve()
    git(root, "init", "-q")
    (root / ".gitignore").write_text("build/\n")
    (root / ".gitattributes").write_text("* text=auto\n")
    (root / "tools").mkdir()
    (root / "tools/host_state.py").write_text(STUB)
    producer = b"print('unique producer of run.json')\n"
    crlf = b"T,K\r\n298.15,1\r\n"
    record(root, "build/issue-9001/producer.py", producer)
    record(root, "build/issue-9002/table.csv", crlf)
    record(root, "build/issue-9003/producer.py", b"print('a producer that will be discarded')\n")
    record(root, "build/issue-9004/build/environment-wheel/build-jobs", b"6")
    log = b"solver log retained gzipped\n"
    record(root, "build/issue-9005/run.log", log)
    archive = root / "analyses/demo/logs/run.log.gz"
    archive.parent.mkdir()
    archive.write_bytes(gzip.compress(b"an older log\n"))
    git(root, "add", ".gitignore", ".gitattributes", "analyses")
    git(root, "commit", "-qm", "records")

    if len(sys.argv) > 1:
        (root / "tools/clean-build").write_text(git(TOOLS, "show", f"{sys.argv[1]}:tools/clean-build"))
        before = subprocess.run([sys.executable, "-c", BEFORE, str(root / "tools/clean-build")], cwd=root,
                                check=True, capture_output=True, text=True).stdout
        assert "EXPIRE issue-9001" in before.splitlines(), before
        print(f"base {sys.argv[1]}: expires issue-9001 (would delete its unique producer)")

    for folder in ("build/issue-9007", "build/issue-9008", "build/issue-9009"):  # the base tool stops on 9008's missing state
        (root / folder).mkdir()
        (root / folder / "scratch.bin").write_bytes(b"\0")
    shutil.copy(TOOLS / "clean-build", root / "tools/clean-build")
    expected = {"issue-9001": "unpromoted", "issue-9002": "unpromoted", "issue-9003": "unpromoted",
                "issue-9004": "expired", "issue-9005": "unpromoted", "issue-9007": "unknown", "issue-9008": "unknown", "issue-9009": "open"}
    assert states(root) == expected, states(root)
    subprocess.run([sys.executable, str(root / "tools/clean-build")], cwd=root, check=True, capture_output=True)
    listing = (root / "build/evidence-to-promote.md").read_text()
    assert (root / "build/issue-9001/producer.py").exists() and "issue-9001/producer.py" in listing, listing
    assert not (root / "build/issue-9004/build").exists()
    assert (root / "build/issue-9007/scratch.bin").exists() and (root / "build/issue-9008/scratch.bin").exists()
    assert (root / "build/issue-9009/scratch.bin").exists()
    print("held: 9001 producer, 9002 CRLF table, 9003 producer, 9005 log, listed and kept by a real run; "
          "9004 build-jobs expired; 9007 (state unknown) and 9008 (no state) kept as unknown; 9009 open and kept")

    archive.write_bytes(gzip.compress(log))
    assert states(root)["issue-9005"] == "unpromoted", "working-copy .gz released the hold"
    git(root, "add", str(archive))
    assert states(root)["issue-9005"] == "unpromoted", "staged .gz released the hold"
    git(root, "commit", "-qm", "retain the log gzipped")
    expected["issue-9005"] = "expired"
    assert states(root) == expected, states(root)
    print("gzip: a working-copy or staged .gz keeps 9005 unpromoted; the committed .gz gives expired")

    # git grep reports an unreadable file only on stderr and exits 1, as if nothing matched.
    for unreadable, refusal in (("analyses/demo/results/issue-9003.json", "recorded hashes"),
                                (".gitignore", "cited build/ paths")):
        (root / unreadable).chmod(0)
        run = subprocess.run([sys.executable, str(root / "tools/clean-build")], cwd=root, capture_output=True,
                             text=True)
        (root / unreadable).chmod(0o644)
        assert run.returncode and refusal in run.stderr and "nothing removed" in run.stderr, (run.stdout, run.stderr)
        assert (root / "build/issue-9003/producer.py").exists()
        print(f"unreadable {unreadable}: the run stops ({run.stderr.strip()}); 9003 keeps its producer")

    (root / "analyses/demo/scripts").mkdir(parents=True)
    (root / "analyses/demo/scripts/producer.py").write_bytes(producer)
    (root / "analyses/demo/inputs").mkdir()
    (root / "analyses/demo/inputs/table.csv").write_bytes(crlf)
    git(root, "add", "analyses")
    git(root, "commit", "-qm", "promote producer; LF-only table")
    assert git(root, "show", "HEAD:analyses/demo/inputs/table.csv").encode() == crlf.replace(b"\r\n", b"\n")
    found = states(root)
    assert found["issue-9001"] == "expired" and found["issue-9002"] == "unpromoted", found
    print("promotion: 9001 expired; an LF-only Git copy of the CRLF table keeps 9002 unpromoted")

    (root / ".gitattributes").write_text("* text=auto\nanalyses/demo/inputs/*.csv -text\n")
    git(root, "add", "--renormalize", ".gitattributes", "analyses")
    git(root, "commit", "-qm", "keep the table's bytes")
    (root / "build/issue-9003/producer.py").unlink()
    (root / "analyses/demo/README.md").write_text("The producer recorded in results/issue-9003.json was lost "
                                                   "with issue 9003's scratch folder.\n")
    git(root, "add", "analyses")
    git(root, "commit", "-qm", "record the lost producer")
    expected.update({"issue-9001": "expired", "issue-9002": "expired", "issue-9003": "expired"})
    assert states(root) == expected, states(root)
    subprocess.run([sys.executable, str(root / "tools/clean-build")], cwd=root, check=True, capture_output=True)
    assert not (root / "build/evidence-to-promote.md").exists()
    print("raw-byte promotion: 9002 expired; discard with a loss note: 9003 expired; listing removed")
