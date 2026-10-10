"""Record the C++ reference rows of the binding checks (#59).

Runs the Python-free native checks of the prepared build and keeps their full-precision
`model` values once, in native-rows.tsv, with the commit and whether engine/ was clean beside them.
Run after `bash .codex/environments/setup.sh`: `.venv/bin/python engine/tests/oracle/record_native_rows.py`.
"""
import json
import subprocess
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
# runner, case argument, retained comparison-id prefixes
CASES = (
    ("properties_checks", "P3", ("p3-form-value/", "p3-total/", "p3-endpoint/")),
    ("properties_checks", "P4", ("direct/action/", "pressure/action/")),
    ("equilibrium_checks", "S2", ("reference-temperature/ionic/", "reference-calorics/ionic/")),
    ("equilibrium_checks", "E5-1", tuple(f"e5-1/neutral-ternary/phase[{phase}].{name}["
                                         for phase in (0, 1) for name in ("amount", "log_activity"))),
    ("equilibrium_checks", "E5-3", ("e5-3/isomerisation/A",)),
    ("equilibrium_checks", "E5-5", ("e5-5/electrolyte/",)),
    ("equilibrium_checks", "E5-7", ("e5-7/thermal-320/",)),
    ("equilibrium_checks", "E5-10", ("e5-10/manifold/value/", "e5-10/manifold/pressure-temperature/",
                                     "e5-10/manifold/temperature-temperature-temperature/")),
    ("regression_checks", "R2-J2", ("r4/",)),
    ("regression_checks", "D86", ("d86/density_value", "d86/parameter/", "d86/temperature/", "d86/pressure")),
    ("surface_example", None, ()),
)
SKIPPED_SUFFIXES = ("/rate", "/perturbed-nu", "/step0", "/step1", "/coarse", "/pv")


def main() -> None:
    build = Path((ROOT / "build/environment-wheel/meson-build-dir").read_text().strip())
    # Test programs are not built by default (#319): build the runners recorded here (#326).
    subprocess.run([str(ROOT / ".venv/bin/meson"), "compile", "-C", str(build),
                    *sorted({runner for runner, _, _ in CASES})], check=True)
    lines, failures, commands = ["source\tid\tunit\tvalue"], [], []
    for runner, case, prefixes in CASES:
        command = [str(build / runner)] + ([case] if case else [])
        done = subprocess.run(command, capture_output=True, text=True, check=False)
        source = runner if case is None else f"{runner} {case}"
        commands.append({"command": source, "exit": done.returncode})
        for line in done.stdout.splitlines():
            if line.startswith("EPCSAFT_ROW\t"):
                lines.append(f"{source}\t" + line.split("\t", 1)[1])
            elif line.startswith("EPCSAFT_COMPARE_JSON "):
                item = json.loads(line.split(" ", 1)[1])
                if item.get("verdict") is False:
                    failures.append(item["comparison_id"])
                name = item["comparison_id"]
                if name.startswith(prefixes) and not name.endswith(SKIPPED_SUFFIXES) and item["model"] is not None:
                    lines.append(f"{source}\t{name}\t{item['unit']}\t{item['model']!r}")
    if any(item["exit"] != 0 for item in commands):
        raise SystemExit(f"a native check failed; rows not recorded: {commands}")
    if repeated := [name for name, n in Counter(line.split("\t")[1] for line in lines[1:]).items() if n > 1]:
        raise SystemExit(f"repeated row ids; rows not recorded: {repeated}")
    (HERE / "native-rows.tsv").write_text("\n".join(lines) + "\n")
    revision = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True, check=True)
    changes = subprocess.run(["git", "status", "--porcelain", "--", "engine"], cwd=ROOT, capture_output=True, text=True, check=True)
    (HERE / "native-rows.json").write_text(json.dumps({
        "rows": len(lines) - 1,
        "recorded_at_revision": revision.stdout.strip(),
        "engine_clean": not changes.stdout,
        "commands": commands,
        "failed_native_comparisons": failures,
        "value": "the native check's unrounded model operand; its own reference and criterion stay with the native owner",
        "claim_limit": "identical explicit inputs through one build; transport agreement, not model validity",
    }, indent=2) + "\n")


if __name__ == "__main__":
    main()
