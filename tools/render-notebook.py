#!/usr/bin/env python3
"""Render current project owners and recorded outputs as Quarto Markdown inputs."""

from __future__ import annotations

import csv
import hashlib
import json
import re
import runpy
import shutil
import subprocess
import sys
import xml.etree.ElementTree as ET
from datetime import datetime, timezone
from pathlib import Path
from typing import Any
from urllib.parse import quote, urljoin

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "build/notebook-data"
PLAN = ROOT / "docs/scientific/implementation-plan.md"
LEGACY = ROOT / "build/numerical-evidence"
MAP = runpy.run_path(str(ROOT / "tools/project-map.py"))
ROADMAP = runpy.run_path(str(ROOT / "tools/check-roadmap.py"))
STATUS = {"[ ]": "Planned", "[-]": "Active", "[!]": "Blocked", "[x]": "Qualified"}


def command(*args: str) -> str:
    return subprocess.run(args, cwd=ROOT, check=True, capture_output=True, text=True).stdout.strip()


def cell(value: Any) -> str:
    return str(value if value not in (None, "") else "—").replace("|", "\\|").replace("\n", " ")


def table(headers: tuple[str, ...], rows: list[tuple[Any, ...]]) -> str:
    lines = ["| " + " | ".join(headers) + " |", "|" + "|".join("---" for _ in headers) + "|"]
    lines.extend("| " + " | ".join(cell(value) for value in row) + " |" for row in rows)
    return "\n".join(lines)


def compact_inputs(value: Any) -> str:
    if not isinstance(value, dict):
        return cell(value)
    keys = ("temperature_K", "pressure_Pa", "molality_mol_kg", "salt", "composition")
    selected = {key: value[key] for key in keys if key in value}
    return cell(selected or value.get("case_id", "—"))


def write(path: Path, content: str) -> None:
    content = content.rstrip() + "\n"
    if not path.is_file() or path.read_text(encoding="utf-8") != content:
        path.parent.mkdir(parents=True, exist_ok=True)
        temporary = path.with_suffix(path.suffix + ".tmp")
        temporary.write_text(content, encoding="utf-8")
        temporary.replace(path)


def repository_url() -> str | None:
    try:
        remote = command("git", "config", "--get", "remote.origin.url")
    except subprocess.CalledProcessError:
        return None
    match = re.search(r"(?:github\.com[:/])([^/]+/[^/]+?)(?:\.git)?$", remote)
    return f"https://github.com/{match.group(1)}" if match else None


def source_link(path: str, tracked: set[str], remote: str | None, revision: str) -> str:
    if path not in tracked or remote is None:
        return f"`{path}`" + (" (untracked)" if path not in tracked else "")
    return f"[{path}]({remote}/blob/{revision}/{quote(path)})"


def roadmap_rows() -> list[dict[str, Any]]:
    return ROADMAP["rows"](PLAN.read_text(encoding="utf-8"))


def rebase_plan_links(record: str) -> str:
    return re.sub(r"\]\(([^)]+)\)", lambda match:
                  "](" + urljoin(PLAN.relative_to(ROOT).as_posix(), match[1]) + ")", record)


def load_json(path: Path) -> tuple[Any | None, str | None]:
    try:
        return json.loads(path.read_text(encoding="utf-8")), None
    except FileNotFoundError:
        return None, "missing"
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        return None, f"unavailable: {error}"


def parse_time(value: str | None) -> datetime | None:
    try:
        parsed = datetime.fromisoformat(value.replace("Z", "+00:00")) if value else None
        return parsed if parsed is not None and parsed.tzinfo is not None else None
    except ValueError:
        return None


def relevant_source_time() -> datetime:
    roots = [ROOT / name for name in ("engine/native", "engine/src", "engine/tests", "engine/data", "data/packets")]
    files = [ROOT / "ARCHITECTURE.yaml", ROOT / "engine/pyproject.toml", ROOT / "uv.lock"]
    for root in roots:
        if root.is_dir():
            files.extend(path for path in root.rglob("*") if path.is_file())
    timestamp = max(path.stat().st_mtime for path in files if path.is_file())
    return datetime.fromtimestamp(timestamp, timezone.utc)


def comparison_rows(run: dict[str, Any]) -> tuple[list[dict[str, str]], str | None, Path | None]:
    relative = run.get("path")
    if not isinstance(relative, str):
        return [], "run record path missing", None
    run_path = (ROOT / relative).resolve()
    if not run_path.is_relative_to(LEGACY.resolve()) or not run_path.is_file():
        return [], "run record unavailable", None
    record, error = load_json(run_path)
    if error or not isinstance(record, dict):
        return [], error or "run record is not an object", None
    csv_path = run_path.with_name("comparisons.csv")
    if not csv_path.is_file():
        return [], "no numerical table recorded", None
    expected = record.get("comparisons_sha256") or run.get("comparisons_sha256")
    actual = hashlib.sha256(csv_path.read_bytes()).hexdigest()
    if expected and expected != actual:
        return [], "comparison table integrity check failed", csv_path
    try:
        with csv_path.open(newline="", encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream))
    except (OSError, UnicodeError, csv.Error) as error:
        return [], f"comparison table unavailable: {error}", csv_path
    return rows, None, csv_path


def copy_table(label: str, path: Path | None) -> str | None:
    if path is None or not path.is_file():
        return None
    destination = OUT / f"{label}.csv"
    if not destination.is_file() or hashlib.sha256(destination.read_bytes()).digest() != hashlib.sha256(path.read_bytes()).digest():
        temporary = destination.with_suffix(".csv.tmp")
        shutil.copyfile(path, temporary)
        temporary.replace(destination)
    return f"build/notebook-data/{label}.csv"


def render_map(rows: list[dict[str, Any]], latest: list[dict[str, Any]]) -> str:
    project = MAP["project_map"]()
    tracked = set(command("git", "ls-files").splitlines())
    remote = repository_url()
    revision = command("git", "rev-parse", "HEAD")
    nodes = {node["id"]: node for node in project["nodes"]}
    children: dict[str, list[str]] = {}
    for node in nodes.values():
        if node.get("parent"):
            children.setdefault(node["parent"], []).append(node["id"])

    def milestones(node_ids: set[str]) -> list[dict[str, Any]]:
        names = {nodes[node_id]["name"] for node_id in node_ids}
        if "applications" in node_ids:
            names.add("Application acceptance")
        return [row for row in rows if row["layer"] in names]

    def checks(node_ids: set[str]) -> list[dict[str, Any]]:
        slugs = {node_id.split(".")[-1].replace("bindings", "") for node_id in node_ids}
        result = []
        for run in latest:
            parts = str(run.get("test_id", "")).split(".")
            if len(parts) > 1 and (parts[1] in slugs or ("engine" in node_ids and parts[0] in {"native", "python"})):
                result.append(run)
        return result

    sections = [
        "The map is derived from `ARCHITECTURE.yaml` and files that exist in this checkout. "
        "Milestones come only from `implementation-plan.md`. Engine aggregates every descendant layer.",
        "",
    ]
    for node in project["nodes"]:
        if node.get("parent"):
            continue
        descendant_ids = {node["id"], *children.get(node["id"], [])}
        sections.extend([f"## {node['name']}", "", node["purpose"], ""])
        owned_milestones = milestones(descendant_ids)
        if owned_milestones:
            sections.extend([
                "### Milestones",
                "",
                table(("ID", "Layer", "Status", "Prerequisites"), [
                    (item["id"], item["layer"], STATUS[item["status"]], item["prerequisites"])
                    for item in owned_milestones
                ]),
                "",
            ])
        owned_checks = checks(descendant_ids)
        if owned_checks:
            sections.extend([
                "### Relevant recorded checks",
                "",
                table(("Check", "Execution", "Produced"), [
                    (run.get("test_id"), run.get("execution"), run.get("finished_at")) for run in owned_checks
                ]),
                "",
            ])
        for role, paths in sorted(node.get("files", {}).items()):
            sections.extend([f"### {role.title()} files", ""])
            sections.extend(f"- {source_link(path, tracked, remote, revision)}" for path in paths)
            sections.append("")
        for summary in node.get("summarized", []):
            sections.append(f"- `{summary['path']}` ({summary['role']}; summarized independent/generated tree)")
        if node.get("summarized"):
            sections.append("")
        for child_id in children.get(node["id"], []):
            child = nodes[child_id]
            sections.extend([f"### {child['name']}", "", child["purpose"], ""])
            child_milestones = milestones({child_id})
            if child_milestones:
                sections.extend([
                    table(("Milestone", "Status", "Requirement"), [
                        (item["id"], STATUS[item["status"]], item["required"]) for item in child_milestones
                    ]),
                    "",
                ])
            for role, paths in sorted(child.get("files", {}).items()):
                sections.extend(f"- **{role}:** {source_link(path, tracked, remote, revision)}" for path in paths)
            sections.append("")
    sections.extend([
        "## Mapping diagnostics",
        "",
        f"- Unmapped existing paths: {len(project['unmapped'])}",
        f"- Ambiguous existing paths: {len(project['ambiguous'])}",
        "",
    ])
    sections.extend(f"- Unmapped: `{path}`" for path in project["unmapped"])
    sections.extend(f"- Ambiguous: `{path}`" for path in project["ambiguous"])
    return "\n".join(sections)


def render_plan(rows: list[dict[str, Any]]) -> str:
    counts = {label: 0 for label in STATUS.values()}
    for row in rows:
        counts[STATUS[row["status"]]] += 1
    lines = [
        "This is a read-only rendering of the authoritative implementation plan. A check outcome does not change qualification.",
        "",
        table(("Planned", "Active", "Blocked", "Qualified"), [tuple(counts[key] for key in ("Planned", "Active", "Blocked", "Qualified"))]),
        "",
    ]
    for layer in ROADMAP["LAYERS"]:
        lines.extend([f"## {layer}", ""])
        lines.append(table(
            ("Milestone", "Status", "Remaining h", "Prerequisites", "Acceptance requirement", "Evidence record"),
            [
                (row["id"], STATUS[row["status"]], row["remaining"], row["prerequisites"], row["required"], rebase_plan_links(row["record"]))
                for row in rows if row["layer"] == layer
            ],
        ))
        lines.append("")
    return "\n".join(lines)


def ordinary_xml_results(document: ET.Element, path: Path, source_time: datetime) -> list[tuple[str, str, str, str, str]]:
    suites = [document] if document.tag == "testsuite" else list(document.findall("testsuite"))
    results = []
    for suite in suites:
        produced_text = suite.get("timestamp") or datetime.fromtimestamp(path.stat().st_mtime, timezone.utc).isoformat()
        try:
            failed = int(suite.get("failures", "0")) + int(suite.get("errors", "0"))
            skipped = int(suite.get("skipped", "0"))
            tests = int(suite.get("tests", "0"))
            produced = parse_time(produced_text)
            if produced is None:
                raise ValueError("timestamp is missing timezone")
            state = "Failed" if failed else ("Incomplete" if skipped == tests else "Passed")
        except (TypeError, ValueError):
            results.append((path.name, "Unavailable", produced_text, "Unverified", "Not qualified"))
            continue
        freshness = "Unverified" if produced >= source_time else "Outdated by modification time"
        results.append((suite.get("name", path.stem), state, produced_text, freshness, "Not qualified"))
    return results


def ordinary_results(source_time: datetime) -> list[tuple[str, str, str, str, str]]:
    results: list[tuple[str, str, str, str, str]] = []
    for path in sorted((ROOT / "build/test-results").glob("*.xml")):
        try:
            document = ET.parse(path).getroot()
        except (OSError, ET.ParseError):
            results.append((path.name, "Unavailable", "—", "Unverified", "Not qualified"))
            continue
        results.extend(ordinary_xml_results(document, path, source_time))
    marker = ROOT / "build/environment-wheel/meson-build-dir"
    if marker.is_file():
        log = Path(marker.read_text(encoding="utf-8").strip()) / "meson-logs/testlog.json"
        if log.is_file():
            try:
                for line in log.read_text(encoding="utf-8").splitlines():
                    record = json.loads(line)
                    produced = datetime.fromtimestamp(float(record["starttime"]) + float(record["duration"]), timezone.utc)
                    results.append((
                        record.get("name", "Meson test"),
                        "Failed" if record.get("is_fail") else "Passed",
                        produced.isoformat(),
                        "Unverified" if produced >= source_time else "Outdated by modification time",
                        "Not qualified",
                    ))
            except (OSError, UnicodeError, ValueError, KeyError, json.JSONDecodeError):
                results.append(("Meson test log", "Unavailable", "—", "Unverified", "Not qualified"))
    for path in sorted((ROOT / "build/test-results").glob("*.log")):
        try:
            events = re.findall(
                r"EPCSAFT_CASE .*?test_id=(\S+) execution=(\S+)",
                path.read_text(encoding="utf-8"),
            )
        except (OSError, UnicodeError):
            events = []
        if events:
            produced = datetime.fromtimestamp(path.stat().st_mtime, timezone.utc)
            test_id, state = events[-1]
            results.append((test_id, state, produced.isoformat(), "Unverified" if produced >= source_time else "Outdated by modification time", "Not qualified"))
    return results


def comparison_events(output: str, fallback_case_id: str) -> list[tuple[str, str]]:
    case_id = fallback_case_id
    events: list[tuple[str, str]] = []
    for line in output.splitlines():
        if line.startswith("EPCSAFT_CASE "):
            case_id = next(
                (part.split("=", 1)[1] for part in line.split() if part.startswith("test_id=")),
                case_id,
            )
        elif line.startswith("EPCSAFT_COMPARE_JSON "):
            events.append((line.split(" ", 1)[1], case_id))
    return events


def current_comparisons(source_time: datetime) -> tuple[list[tuple[Any, ...]], list[str], list[str]]:
    """Read comparison JSON emitted by ordinary reports; never execute producers."""
    records: list[dict[str, Any]] = []
    errors: list[str] = []
    for path in sorted((ROOT / "build/test-results").glob("*.xml")):
        try:
            root = ET.parse(path).getroot()
        except (OSError, ET.ParseError) as error:
            errors.append(f"{path.name}: {error}")
            continue
        for testcase in root.iter("testcase"):
            output = "\n".join(node.text or "" for node in testcase.findall("system-out"))
            for payload, case_id in comparison_events(output, testcase.get("name", "")):
                try:
                    record = json.loads(payload)
                    if not isinstance(record, dict):
                        raise ValueError("comparison is not an object")
                    record.update(report=path.name, testcase=testcase.get("name", ""), case_id=case_id)
                    records.append(record)
                except (json.JSONDecodeError, ValueError, IndexError) as error:
                    errors.append(f"{path.name}/{testcase.get('name', '')}: {error}")
    records.sort(key=lambda item: (str(item.get("report")), str(item.get("testcase")), str(item.get("comparison_id"))))
    fieldnames = ("report", "testcase", "case_id", "comparison_id", "quantity", "model", "reference", "unit", "inputs", "criterion", "judgment", "evidence_status", "basis", "method", "claim", "provenance")
    rows = [tuple(item.get(field, "") for field in fieldnames) for item in records]
    csv_paths: list[str] = []
    for report in sorted({str(item[0]) for item in rows}):
        destination = OUT / f"{Path(report).stem}-comparisons.csv"
        temporary = destination.with_suffix(".csv.tmp")
        with temporary.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.writer(stream)
            writer.writerow(fieldnames)
            for row in rows:
                if row[0] == report:
                    writer.writerow(json.dumps(value, ensure_ascii=False) if isinstance(value, (dict, list)) else value for value in row)
        temporary.replace(destination)
        csv_paths.append(f"build/notebook-data/{destination.name}")
    return rows, csv_paths, errors


def render_results(rows: list[dict[str, Any]]) -> tuple[str, list[dict[str, Any]]]:
    manifest, manifest_error = load_json(LEGACY / "manifest.json")
    runs = manifest.get("runs", []) if isinstance(manifest, dict) and isinstance(manifest.get("runs"), list) else []
    latest = [run for run in runs if str(run.get("test_id", "")).startswith(("native.", "python."))]
    source_time = relevant_source_time()
    lines = [
        "Execution outcome, production time, source applicability, and scientific qualification are independent. "
        "Values below remain visible when preview is stopped because they are rendered into this page.",
        "",
        "::: {.callout-important}",
        "Historical passes do not claim the current checkout passed. No official milestone is qualified unless the Plan page says so.",
        ":::",
        "",
        "## Ordinary current check reports",
        "",
    ]
    ordinary = ordinary_results(source_time)
    if ordinary:
        lines.append(table(("Report", "Execution", "Produced", "Source / inputs", "Scientific qualification"), ordinary))
    else:
        lines.append("No ordinary pytest XML report has been recorded yet.")
    lines.extend(f"- Raw report: [{path.name}](build/test-results/{path.name})" for path in sorted((ROOT / "build/test-results").glob("*.xml")))
    current_rows, comparison_csvs, comparison_errors = current_comparisons(source_time)
    lines.extend(["", "### Current numerical comparisons", "", "The first deterministic row per report/test case is shown; complete CSVs and raw XML remain downloadable."])
    if current_rows:
        shown: list[tuple[Any, ...]] = []
        seen: set[tuple[str, str]] = set()
        for row in current_rows:
            key = (str(row[0]), str(row[3]))
            if key not in seen:
                seen.add(key)
                shown.append(row)
        lines.append(table(("Report", "Test", "Case", "Comparison", "Model", "Reference", "Unit", "Criterion", "Judgment", "Inputs"), [(row[0], row[1], row[2], row[3], row[5], row[6], row[7], row[9], row[10], compact_inputs(row[8])) for row in shown]))
        lines.extend(f"- [{path.rsplit('/', 1)[-1]}]({path})" for path in comparison_csvs)
        lines.extend(f"- [{path.name}](build/test-results/{path.name})" for path in sorted((ROOT / "build/test-results").glob("*.xml")))
    else:
        lines.append("No EPCSAFT_COMPARE_JSON rows were readable.")
    lines.append(f"Parse errors: {len(comparison_errors)}; comparison rows: {len(current_rows)}.")
    lines.extend(f"- {error}" for error in comparison_errors)
    provenance = ROOT / "build/environment-wheel/preparation.json"
    lines.extend(["", "### Preparation provenance", "", f"- [{provenance.name}](build/environment-wheel/preparation.json)" if provenance.is_file() else "- preparation.json unavailable"])
    if current_rows:
        lines.append("Current applicability is Unverified when the report is newer than the source boundary; modification time alone does not establish Current.")
    project = MAP["project_map"]()
    analyses = next(node for node in project["nodes"] if node["id"] == "analyses")
    analysis_outputs = sorted(
        path
        for path in analyses.get("files", {}).get("evidence", [])
        if "/results/" in path and Path(path).suffix.lower() in {".csv", ".json", ".tsv"}
    )
    tracked = set(command("git", "ls-files").splitlines())
    remote = repository_url()
    revision = command("git", "rev-parse", "HEAD")
    lines.extend([
        "",
        "## Retained analysis outputs",
        "",
        "These are analysis-owned records, not notebook publications.",
        "",
        table(("Recorded output",), [(source_link(path, tracked, remote, revision),) for path in analysis_outputs]),
        "",
        "## Last legacy scientific records",
        "",
    ])
    if manifest_error:
        lines.extend([f"::: {{.callout-warning}}\nLegacy manifest {manifest_error}. Historical values already rendered into an older page remain there until the next successful render.\n:::", ""])
    if latest:
        summary = []
        for run in latest:
            produced = parse_time(run.get("finished_at"))
            freshness = "Unverified" if produced and produced >= source_time else "Outdated by modification time"
            summary.append((run.get("test_id"), run.get("execution"), run.get("finished_at"), freshness, run.get("outcome_reason"), "Not qualified"))
        lines.extend([table(("Check", "Execution", "Produced", "Source / inputs", "Outcome / reason", "Scientific qualification"), summary), ""])
    else:
        lines.extend(["No readable legacy scientific run manifest is available.", ""])

    interrupted: dict[str, Any] | None = None
    for path in LEGACY.glob("runs/*/run.json"):
        record, error = load_json(path)
        if error or not isinstance(record, dict) or record.get("execution") != "Interrupted":
            continue
        if interrupted is None or str(record.get("finished_at", "")) > str(interrupted.get("finished_at", "")):
            interrupted = record
    lines.extend(["## Incomplete and unavailable evidence", ""])
    if interrupted:
        lines.append(table(("Check", "Execution", "Produced", "Reason"), [(
            interrupted.get("test_id"), "Incomplete (interrupted)", interrupted.get("finished_at"), interrupted.get("outcome_reason")
        )]))
    else:
        lines.append("No interrupted historical run is present. Missing or corrupt inputs are reported above as unavailable.")
    lines.append("")

    for label, identifiers in (
        ("F1 retained residual-energy evidence", ("native.properties.f1", "native.model.f1")),
        ("F5 established curvature disagreement", ("native.equilibrium.f5",)),
    ):
        run = next((item for identifier in identifiers for item in reversed(latest) if item.get("test_id") == identifier), None)
        lines.extend([f"## {label}", ""])
        if run is None:
            lines.extend(["Unavailable: no retained run record.", ""])
            continue
        comparisons, error, path = comparison_rows(run)
        download = copy_table("f1-comparisons" if "F1" in label else "f5-comparisons", path)
        lines.extend([
            f"- Execution: **{cell(run.get('execution'))}**",
            f"- Produced: {cell(run.get('finished_at'))}",
            f"- Outcome detail: {cell(run.get('outcome_reason'))}",
            "- Applicability: **Outdated** after the notebook cutover changed the recorded input set; values are historical, not a current-code pass.",
            "- Qualification: **Not qualified**; the implementation plan remains authoritative.",
            f"- Full retained table: [{path.name if path else 'unavailable'}]({download})" if download else f"- Full retained table: unavailable ({cell(error)})",
            "",
        ])
        if error and not comparisons:
            lines.extend([f"::: {{.callout-warning}}\n{error}\n:::", ""])
            continue
        failed = [item for item in comparisons if str(item.get("verdict", "")).lower() in {"false", "0"} or item.get("judgment") == "Criterion failed"]
        shown = comparisons if "F1" in label else failed
        lines.extend([
            table(("Comparison", "Quantity", "Model", "Reference", "Unit", "Judgment", "Inputs"), [
                (item.get("comparison_id"), item.get("quantity"), item.get("model"), item.get("reference"), item.get("unit"), item.get("judgment"), item.get("inputs"))
                for item in shown
            ]) if shown else "No failing comparison row was readable.",
            "",
            f"Retained rows: {len(comparisons)}; failing rows shown: {len(failed)}.",
            "",
        ])
    lines.extend([
        "## Freshness limits",
        "",
        f"Relevant Engine source/input modification boundary: `{source_time.isoformat()}`.",
        "Legacy records retain their own wheel and physical-input provenance. The notebook does not rebuild, collect, execute, retry, or re-certify them.",
    ])
    return "\n".join(lines), latest


def main() -> None:
    if "--self-test" in sys.argv:
        source_time = datetime(2026, 1, 1, tzinfo=timezone.utc)
        invalid_counts = ET.fromstring('<testsuite tests="not-a-number" timestamp="2026-09-09T12:00:00Z"/>')
        naive_time = ET.fromstring('<testsuite tests="1" timestamp="2026-09-09T12:00:00"/>')
        assert ordinary_xml_results(invalid_counts, Path("invalid-counts.xml"), source_time)[0][1:] == (
            "Unavailable", "2026-09-09T12:00:00Z", "Unverified", "Not qualified"
        )
        assert ordinary_xml_results(naive_time, Path("naive-time.xml"), source_time)[0][1:] == (
            "Unavailable", "2026-09-09T12:00:00", "Unverified", "Not qualified"
        )
        sample = (
            "EPCSAFT_CASE schema_version=1 test_id=native.contributions.p1 execution=Passed\n"
            'EPCSAFT_COMPARE_JSON {"comparison_id":"p1"}\n'
            "EPCSAFT_CASE schema_version=1 test_id=native.contributions.p2 execution=Passed\n"
            'EPCSAFT_COMPARE_JSON {"comparison_id":"p2"}\n'
        )
        assert [case_id for _, case_id in comparison_events(sample, "python.case")] == [
            "native.contributions.p1",
            "native.contributions.p2",
        ]
        assert rebase_plan_links(
            "[context](CONTEXT.md#anchor) [review](reviews/c1.md) [issue](https://example.test/6)"
        ) == (
            "[context](docs/scientific/CONTEXT.md#anchor) "
            "[review](docs/scientific/reviews/c1.md) [issue](https://example.test/6)"
        )
        retained = ET.parse(ROOT / "build/test-results/native.xml").getroot()
        retained_output = "\n".join(
            node.text or ""
            for testcase in retained.iter("testcase")
            for node in testcase.findall("system-out")
        )
        assert any(
            case_id == "native.contributions.p2"
            and json.loads(payload).get("comparison_id") == "density-monotonic/0"
            for payload, case_id in comparison_events(retained_output, "native")
        )
        print("reader self-test passed")
        return
    OUT.mkdir(parents=True, exist_ok=True)
    rows = roadmap_rows()
    results, latest = render_results(rows)
    write(OUT / "plan.md", render_plan(rows))
    write(OUT / "map.md", render_map(rows, latest))
    write(OUT / "results.md", results)


if __name__ == "__main__":
    main()
