"""Record the C++ reference rows of the installed-wheel binding checks (#59).

Runs the Python-free native checks of the prepared build and keeps their full-precision
`model` values once, in native-rows.tsv, with the build and wheel identity beside them.
Run after `bash .codex/environments/setup.sh`: `.venv/bin/python engine/tests/oracle/record_native_rows.py`.
"""
import json
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
# runner, case argument, retained comparison-id prefixes
CASES = (
    ("properties_checks", "P3", ("p3-form-value/", "p3-total/", "p3-endpoint/")),
    ("properties_checks", "P4", ("direct/action/", "pressure/action/")),
    ("equilibrium_checks", "S2", ("reference-temperature/", "reference-calorics/")),
    ("equilibrium_checks", "E5-1", ("e5-1/neutral-ternary/",)),
    ("equilibrium_checks", "E5-3", ("e5-3/isomerisation",)),
    ("equilibrium_checks", "E5-5", ("e5-5/electrolyte/",)),
    ("equilibrium_checks", "E5-7", ("e5-7/thermal-320/",)),
    ("equilibrium_checks", "E5-10", ("e5-10/manifold/",)),
    ("regression_checks", "R2-J2", ("r4/",)),
    ("regression_checks", "D86", ("d86/density_value", "d86/parameter/", "d86/temperature/", "d86/pressure")),
    ("surface_example", None, ()),
)
SKIPPED_SUFFIXES = ("/rate", "/perturbed-nu", "/step0", "/step1", "/coarse")


def main() -> None:
    preparation = json.loads((ROOT / "build/environment-wheel/preparation.json").read_text())
    lines, failures, commands = ["source\tid\tunit\tvalue"], [], []
    for runner, case, prefixes in CASES:
        command = [str(ROOT / preparation["build_dir"] / runner)] + ([case] if case else [])
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
    (HERE / "native-rows.tsv").write_text("\n".join(lines) + "\n")
    revision = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True, check=True)
    (HERE / "native-rows.json").write_text(json.dumps({
        "rows": len(lines) - 1,
        "recorded_at_revision": revision.stdout.strip(),
        "build_fingerprint": preparation["build_fingerprint"],
        "wheel_sha256": preparation["wheel_sha256"],
        "commands": commands,
        "failed_native_comparisons": failures,
        "value": "the native check's unrounded model operand; its own reference and criterion stay with the native owner",
        "claim_limit": "identical explicit inputs through one build; transport agreement, not model validity",
    }, indent=2) + "\n")


if __name__ == "__main__":
    main()
