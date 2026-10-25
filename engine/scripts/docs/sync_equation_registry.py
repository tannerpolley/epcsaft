"""Generate the equation index from the canonical Markdown definitions.

Edit engine/docs/equations.md. --check detects malformed definitions, duplicate
IDs, unknown EqID references and a stale YAML index. It does not verify math.
"""

from __future__ import annotations

import argparse
import json
import re
import tempfile
from pathlib import Path

FIELDS = ("eqid", "section", "subsection", "source", "status", "description", "change_note", "math")
ALLOWED_STATUSES = frozenset(
    {
        "Source equation",
        "Exact derivation",
        "Bounded empirical relation",
        "Implementation convention",
        "Documentation-only",
    }
)
EQID_PATTERN = re.compile(r"^[a-z][a-z0-9_]*$")
TEX_ONLY = re.compile(r"\\(cite|label|ref|eqref)\b")


def load_entries(path: Path) -> list[dict[str, str]]:
    text = path.read_text(encoding="utf-8")
    blocks = re.split(r"^## `([^`]+)`\s*$", text, flags=re.MULTILINE)
    if any(not re.fullmatch(r"## `[^`]+`", line) for line in text.splitlines() if line.startswith("## ")):
        raise ValueError(f"{path}: malformed equation heading")
    if len(blocks) == 1:
        raise ValueError(f"{path}: no equation definitions")
    entries = []
    for eqid, body in zip(blocks[1::2], blocks[2::2], strict=True):
        match = re.fullmatch(
            r"\s*- Section: ([^\n]*)\n- Status: ([^\n]*)\n"
            r"- Source: ([^\n]*)\n- Description: ([^\n]*)\n"
            r"- Change note: ([^\n]*)\n\s*\$\$\n(.*?)\n\$\$\s*", body, re.DOTALL,
        )
        if match is None:
            raise ValueError(f"{path}: malformed definition {eqid!r}")
        section, status, source, description, note, math = match.groups()
        section, _, subsection = section.partition(" / ")
        entries.append(dict(zip(FIELDS, (eqid, section, subsection, source,
                                         status, description, note, math), strict=True)))
    return entries


def validate(entries: list[dict[str, str]]) -> None:
    seen: set[str] = set()
    for entry in entries:
        eqid = entry["eqid"]
        if not EQID_PATTERN.match(eqid):
            raise ValueError(f"malformed eqid {eqid!r}")
        if eqid in seen:
            raise ValueError(f"duplicate eqid {eqid!r}")
        seen.add(eqid)
        if entry["status"] not in ALLOWED_STATUSES:
            raise ValueError(f"equation {eqid!r} has unsupported status {entry['status']!r}")
        for field in ("source", "description", "math"):
            if not entry[field].strip():
                raise ValueError(f"equation {eqid!r} has an empty {field}")
        for field in ("source", "math"):
            if TEX_ONLY.search(entry[field]):
                raise ValueError(f"equation {eqid!r}: {field} carries a LaTeX-only command; cite in words, no labels")


def _quoted(value: str) -> str:
    return json.dumps(value, ensure_ascii=False)


def render_registry(entries: list[dict[str, str]]) -> str:
    lines = ["schema_version: 3", "equations:"]
    for entry in sorted(entries, key=lambda item: item["eqid"]):
        lines.append(f"  - eqid: {_quoted(entry['eqid'])}")
        lines.extend(f"    {field}: {_quoted(entry[field])}" for field in FIELDS[1:])
    return "\n".join(lines) + "\n"


EQID_MENTION = re.compile(r"EqID:\s*(`[^`]*`|[^\s]+)")


def documentation_paths(root: Path) -> list[Path]:
    paths = list((root / "docs").rglob("*.md"))
    scientific = root.parent / "docs" / "scientific"
    paths += [p for p in scientific.rglob("*.md")
              if "reviews" not in p.relative_to(scientific).parts]
    return paths


def unregistered_mentions(root: Path, registered: set[str]) -> list[str]:
    """Check live equation/algorithm references; histories do not declare capabilities."""
    algorithm_text = (root / "docs" / "science" / "algorithms.md").read_text(encoding="utf-8")
    algorithms = re.findall(r"^## (ALG-.*)$", algorithm_text, re.MULTILINE)
    if (not algorithms or len(algorithms) != len(set(algorithms))
            or any(not re.fullmatch(r"ALG-[A-Z]+", name) for name in algorithms)):
        raise ValueError("missing or duplicate algorithm IDs")
    paths = documentation_paths(root)
    paths += [p for p in (root / "native").rglob("*") if p.suffix in {".hpp", ".cpp"}]
    bad = []
    for path in sorted(paths):
        text = path.read_text(encoding="utf-8")
        for pattern, known in ((EQID_MENTION, registered), (re.compile(r"\bALG-[A-Za-z0-9_-]+"), set(algorithms))):
            for token in pattern.findall(text):
                name = token[1:-1] if token.startswith("`") and token.endswith("`") else token.rstrip(".,;:)]\"'")
                if name not in known:
                    bad.append(f"{path.relative_to(root.parent)}: unregistered ID {name!r}")
    return bad


def synchronize(docs: Path, *, check: bool) -> None:
    registry = docs / "equations_registry.yaml"
    entries = load_entries(docs / "equations.md")
    validate(entries)
    unregistered = unregistered_mentions(docs.parent, {entry["eqid"] for entry in entries})
    if unregistered:
        raise ValueError("\n".join(unregistered))
    outputs = {registry: render_registry(entries)}
    stale = [path for path, text in outputs.items() if not path.exists() or path.read_text(encoding="utf-8") != text]
    if check and stale:
        raise ValueError("stale: " + ", ".join(path.name for path in stale))
    if not check:
        for path, text in outputs.items():
            path.write_text(text, encoding="utf-8")


def self_check(docs: Path) -> None:
    """Small readiness falsifier: round trip and deliberate corruptions."""
    source = (docs / "equations.md").read_text(encoding="utf-8")
    with tempfile.TemporaryDirectory(prefix="epcsaft-doc-check-") as temporary:
        target = Path(temporary) / "engine" / "docs"
        (target / "science").mkdir(parents=True)
        (target / "science" / "algorithms.md").write_text(
            (docs / "science" / "algorithms.md").read_text(encoding="utf-8"), encoding="utf-8")
        canonical = target / "equations.md"
        canonical.write_text(source, encoding="utf-8")
        synchronize(target, check=False)
        synchronize(target, check=True)
        scientific = Path(temporary) / "docs" / "scientific"
        (scientific / "reviews").mkdir(parents=True)
        context = scientific / "CONTEXT.md"
        context.write_text("EqID: missing_equation\n", encoding="utf-8")
        (scientific / "reviews" / "dated.md").write_text(
            "EqID: historical_equation\n", encoding="utf-8")
        try:
            synchronize(target, check=True)
        except ValueError as error:
            assert "missing_equation" in str(error) and "historical_equation" not in str(error)
        else:
            raise AssertionError("live scientific references were not checked")
        context.write_text("", encoding="utf-8")
        synchronize(target, check=True)
        first = re.search(r"^## `[^`]+`.*?(?=^## `|\Z)", source, re.MULTILINE | re.DOTALL)
        assert first is not None
        for changed, expected in (
            (source + first.group(), "duplicate eqid"),
            (source.replace("- Status:", "- Broken:", 1), "malformed definition"),
            ("EqID: missing_equation\n" + source, "unregistered ID"),
            ("EqID: association_mass_action-typo\n" + source, "unregistered ID"),
            ("EqID: Unknown\n" + source, "unregistered ID"),
            ("ALG-ADMISSION-BOGUS\n" + source, "unregistered ID"),
            ("EqID: `association_mass_actionX\n" + source, "unregistered ID"),
            ("ALG-ADMISSION-\n" + source, "unregistered ID"),
            (source.replace("- Description:", "- Description: changed", 1), "stale:"),
        ):
            canonical.write_text(changed, encoding="utf-8")
            try:
                synchronize(target, check=True)
            except ValueError as error:
                if expected not in str(error):
                    raise AssertionError(f"expected {expected}: {error}") from error
            else:
                raise AssertionError(f"corruption did not fail: {expected}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="verify without writing")
    args = parser.parse_args()
    docs = Path(__file__).resolve().parents[2] / "docs"
    if args.check:
        self_check(docs)
    synchronize(docs, check=args.check)


if __name__ == "__main__":
    main()
