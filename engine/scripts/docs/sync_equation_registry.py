"""Synchronize EOS equation documentation with implementation owners."""

from __future__ import annotations

import argparse
import ast
import json
import re
from dataclasses import dataclass, replace
from pathlib import Path

ALLOWED_STATUSES = frozenset(
    {
        "Source equation",
        "Exact derivation",
        "Implementation convention",
        "Bounded empirical relation",
        "Documentation-only",
    }
)

EQID_PATTERN = re.compile(r"^[a-z][a-z0-9_]*$")
METADATA_PATTERN = re.compile(r"^%\s*(EqID|Source|Status|Description|Change note):\s*(.+?)\s*$")
LABEL_PATTERN = re.compile(r"\\label\{([^}]+)\}")
ANNOTATION_PATTERN = re.compile(r"^\s*(?://|#)\s*(EqID|EqExclude):\s*(.+?)\s*$")
CPP_NAME_PATTERN = re.compile(r"((?:\w+::)*~?\w+)\s*$")
CPP_DECLARATION_PREFIX_PATTERN = re.compile(
    r"^(?:"
    r"inline|static|constexpr|consteval|virtual|explicit|const|"
    r"bool|void|int|double|Scalar|"
    r"std::|detail::|"
    r"[A-Z]\w*"
    r")\b"
)
CPP_CONTROL_PREFIXES = (
    "if",
    "for",
    "while",
    "switch",
    "catch",
    "return",
    "throw",
    "static_cast",
    "reinterpret_cast",
    "const_cast",
    "dynamic_cast",
)

ALLOWED_EXCLUSIONS = frozenset(
    {
        "scalar utility dispatch",
        "numerical iteration mechanics",
        "input and domain validation",
        "transport and buffer mechanics",
        "serialization",
        "CppAD-generated derivative machinery",
    }
)

MATH_OWNER_FILES = (
    "src/epcsaft/_native/helmholtz.hpp",
    "src/epcsaft/_native/association.hpp",
    "src/epcsaft/_native/electrolyte.hpp",
    "src/epcsaft/_native/born.hpp",
    "src/epcsaft/_native/derivatives.hpp",
    "src/epcsaft/_native/density.hpp",
    "src/epcsaft/_native/electrolyte_reference.hpp",
    "src/epcsaft/_native/eos.cpp",
    "src/epcsaft/_native/phase_block.hpp",
    "cpp/equilibrium/source_reference.cpp",
    "src/epcsaft/eos.py",
    "src/epcsaft/reference.py",
)


@dataclass(frozen=True)
class SourceReference:
    file: str
    line: int
    comment_line: int
    symbol: str
    context: str


@dataclass(frozen=True)
class SourceExclusion:
    file: str
    line: int
    symbol: str
    reason: str


@dataclass(frozen=True)
class EquationEntry:
    eqid: str
    section: str
    subsection: str
    label: str
    source: str
    status: str
    description: str
    change_note: str
    latex: str
    tex_file: str
    tex_line: int
    implementation_refs: tuple[SourceReference, ...] = ()


def parse_equations(tex_path: Path, repo_root: Path) -> list[EquationEntry]:
    """Parse checked equation blocks and their required metadata."""
    lines = tex_path.read_text(encoding="utf-8").splitlines()
    section = ""
    subsection = ""
    entries: list[EquationEntry] = []
    seen: set[str] = set()
    metadata: dict[str, str] = {}
    metadata_line = 0
    index = 0

    while index < len(lines):
        line = lines[index]
        if line.startswith(r"\section{"):
            section = line.removeprefix(r"\section{").removesuffix("}")
            subsection = ""
        elif line.startswith(r"\subsection{"):
            subsection = line.removeprefix(r"\subsection{").removesuffix("}")

        match = METADATA_PATTERN.match(line)
        if match:
            if match.group(1) == "EqID":
                metadata = {}
                metadata_line = index + 1
            metadata[match.group(1)] = match.group(2)
            index += 1
            continue

        if line.strip() == r"\begin{equation}":
            end = index + 1
            while end < len(lines) and lines[end].strip() != r"\end{equation}":
                end += 1
            if end == len(lines):
                raise ValueError(f"{tex_path}:{index + 1}: unterminated equation")
            required = {"EqID", "Source", "Status", "Description", "Change note"}
            missing = sorted(required - metadata.keys())
            if missing:
                raise ValueError(f"{tex_path}:{index + 1}: missing metadata: {', '.join(missing)}")
            eqid = metadata["EqID"]
            if not EQID_PATTERN.fullmatch(eqid):
                raise ValueError(f"{tex_path}:{metadata_line}: invalid EqID {eqid!r}")
            if eqid in seen:
                raise ValueError(f"{tex_path}:{metadata_line}: duplicate EqID {eqid!r}")
            if metadata["Status"] not in ALLOWED_STATUSES:
                raise ValueError(
                    f"{tex_path}:{metadata_line}: unsupported status {metadata['Status']!r}"
                )
            body_lines = lines[index + 1 : end]
            labels = [
                label for body_line in body_lines for label in LABEL_PATTERN.findall(body_line)
            ]
            if len(labels) != 1:
                raise ValueError(f"{tex_path}:{index + 1}: equation requires exactly one label")
            entries.append(
                EquationEntry(
                    eqid=eqid,
                    section=section,
                    subsection=subsection,
                    label=labels[0],
                    source=metadata["Source"],
                    status=metadata["Status"],
                    description=metadata["Description"],
                    change_note=metadata["Change note"],
                    latex="\n".join(body_lines),
                    tex_file=str(tex_path.relative_to(repo_root)),
                    tex_line=index + 1,
                )
            )
            seen.add(eqid)
            metadata = {}
            index = end
        index += 1

    return entries


@dataclass(frozen=True)
class _FunctionSpan:
    symbol: str
    start: int
    end: int


def _python_function_spans(path: Path) -> list[_FunctionSpan]:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    return [
        _FunctionSpan(
            symbol=node.name,
            start=node.lineno,
            end=node.end_lineno or node.lineno,
        )
        for node in ast.walk(tree)
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
    ]


def _closing_brace_line(lines: list[str], start: int) -> int:
    depth = 0
    opened = False
    for line_number in range(start, len(lines) + 1):
        line = lines[line_number - 1]
        depth += line.count("{")
        if "{" in line:
            opened = True
        depth -= line.count("}")
        if opened and depth == 0:
            return line_number
    return len(lines)


def _cpp_function_spans(path: Path) -> list[_FunctionSpan]:
    text = path.read_text(encoding="utf-8")
    lines = text.splitlines()
    class_names = set(re.findall(r"\b(?:class|struct)\s+(\w+)", text))
    spans: list[_FunctionSpan] = []
    for offset, line in enumerate(lines):
        stripped = line.strip()
        if (
            "(" not in stripped
            or stripped.startswith(CPP_CONTROL_PREFIXES)
            or stripped.startswith(("}", "|", "&", "!", "?", ":", "+", "-", "*", "/"))
        ):
            continue
        prefix = stripped.split("(", 1)[0]
        if "=" in prefix or "[" in prefix:
            continue
        name_match = CPP_NAME_PATTERN.search(prefix)
        if not name_match:
            continue
        name = name_match.group(1)
        qualified_parts = name.split("::")
        constructor = len(qualified_parts) >= 2 and qualified_parts[-1] == qualified_parts[-2]
        qualified_definition = "::" in name and line == line.lstrip()
        class_body_method = name.removeprefix("~") in class_names
        if (
            (
                not CPP_DECLARATION_PREFIX_PATTERN.match(prefix)
                or not any(character.isspace() for character in prefix)
            )
            and not constructor
            and not qualified_definition
            and not class_body_method
        ):
            continue
        if name in {"sizeof", "alignof", "decltype"}:
            continue
        signature = stripped
        cursor = offset
        paren_depth = stripped.count("(") - stripped.count(")")
        while cursor + 1 < len(lines) and (
            paren_depth > 0
            or (
                "{" not in signature
                and ";" not in signature
                and "=" not in signature.split(")", 1)[-1]
            )
        ):
            cursor += 1
            continuation = lines[cursor].strip()
            signature += " " + continuation
            paren_depth += continuation.count("(") - continuation.count(")")
        tail = signature.rsplit(")", 1)[-1] if ")" in signature else signature
        if "{" not in tail or ";" in tail.split("{", 1)[0]:
            continue
        start = offset + 1
        spans.append(
            _FunctionSpan(
                symbol=name.split("::")[-1],
                start=start,
                end=_closing_brace_line(lines, cursor + 1),
            )
        )
    return spans


def _function_spans(path: Path) -> list[_FunctionSpan]:
    if path.suffix == ".py":
        return _python_function_spans(path)
    return _cpp_function_spans(path)


def _annotation_owner(comment_line: int, functions: list[_FunctionSpan]) -> _FunctionSpan | None:
    containing = [
        function for function in functions if function.start <= comment_line <= function.end
    ]
    if containing:
        return min(containing, key=lambda function: function.end - function.start)
    following = [function for function in functions if function.start > comment_line]
    return min(following, key=lambda function: function.start, default=None)


def parse_source_annotations(
    source_root: Path,
    repo_root: Path,
    *,
    require_classification: bool = True,
) -> tuple[dict[str, list[SourceReference]], list[SourceExclusion]]:
    """Parse equation links and require every discovered function to be classified."""
    references: dict[str, list[SourceReference]] = {}
    exclusions: list[SourceExclusion] = []
    unclassified: list[str] = []
    paths = (
        [source_root]
        if source_root.is_file()
        else sorted(
            path
            for path in source_root.rglob("*")
            if path.is_file() and path.suffix in {".cpp", ".hpp", ".h", ".py"}
        )
    )
    for path in paths:
        lines = path.read_text(encoding="utf-8").splitlines()
        functions = _function_spans(path)
        classified: set[tuple[str, int]] = set()
        for comment_line, line in enumerate(lines, start=1):
            match = ANNOTATION_PATTERN.match(line)
            if not match:
                continue
            owner = _annotation_owner(comment_line, functions)
            if owner is None:
                raise ValueError(f"{path}:{comment_line}: annotation has no function owner")
            classified.add((owner.symbol, owner.start))
            relative = str(path.relative_to(repo_root))
            if match.group(1) == "EqID":
                eqid = match.group(2)
                if not EQID_PATTERN.fullmatch(eqid):
                    raise ValueError(f"{path}:{comment_line}: invalid EqID {eqid!r}")
                references.setdefault(eqid, []).append(
                    SourceReference(
                        file=relative,
                        line=owner.start,
                        comment_line=comment_line,
                        symbol=owner.symbol,
                        context=line.strip(),
                    )
                )
            else:
                reason = match.group(2)
                if reason not in ALLOWED_EXCLUSIONS:
                    raise ValueError(f"{path}:{comment_line}: unsupported EqExclude {reason!r}")
                exclusions.append(
                    SourceExclusion(
                        file=relative,
                        line=owner.start,
                        symbol=owner.symbol,
                        reason=reason,
                    )
                )
        if require_classification:
            for function in functions:
                if (function.symbol, function.start) not in classified:
                    unclassified.append(
                        f"{path.relative_to(repo_root)}:{function.start}:{function.symbol}"
                    )
    if unclassified:
        raise ValueError("unclassified mathematical-owner function: " + ", ".join(unclassified))
    return references, exclusions


def load_bibliography_keys(repo_root: Path) -> set[str]:
    """Return biblatex entry keys from the checked equation bibliography."""
    bibliography = repo_root / "docs" / "latex" / "references.bib"
    if not bibliography.exists():
        return set()
    return set(
        re.findall(
            r"@\w+\s*\{\s*([^,\s]+)",
            bibliography.read_text(encoding="utf-8"),
        )
    )


def _owner_source_root(repo_root: Path) -> Path:
    return repo_root / "src"


def load_repository_inputs(
    repo_root: Path,
) -> tuple[
    list[EquationEntry],
    dict[str, list[SourceReference]],
    list[SourceExclusion],
]:
    """Load and link the canonical LaTeX and declared implementation owners."""
    entries = parse_equations(repo_root / "docs" / "latex" / "equations.tex", repo_root)
    declared = [repo_root / relative for relative in MATH_OWNER_FILES]
    existing = [path for path in declared if path.exists()]
    if len(existing) == len(declared):
        references: dict[str, list[SourceReference]] = {}
        exclusions: list[SourceExclusion] = []
        for path in existing:
            path_references, path_exclusions = parse_source_annotations(
                path,
                repo_root,
                require_classification=False,
            )
            for eqid, refs in path_references.items():
                references.setdefault(eqid, []).extend(refs)
            exclusions.extend(path_exclusions)
    else:
        references, exclusions = parse_source_annotations(_owner_source_root(repo_root), repo_root)
    linked = [
        replace(
            entry,
            implementation_refs=tuple(references.get(entry.eqid, ())),
        )
        for entry in entries
    ]
    return linked, references, exclusions


def validate_registry(
    entries: list[EquationEntry],
    references: dict[str, list[SourceReference]],
    exclusions: list[SourceExclusion],
    bibliography_keys: set[str],
) -> None:
    """Reject incomplete equation, source, citation, or exclusion linkage."""
    entry_ids = {entry.eqid for entry in entries}
    unknown_ids = sorted(set(references) - entry_ids)
    if unknown_ids:
        raise ValueError("source EqID absent from LaTeX: " + ", ".join(unknown_ids))
    for entry in entries:
        refs = references.get(entry.eqid) or list(entry.implementation_refs)
        if entry.status != "Documentation-only" and not refs:
            raise ValueError(f"equation {entry.eqid!r} has no implementation owner")
        citations = re.findall(r"\\cite\{([^}]+)\}", entry.source)
        for citation_group in citations:
            for citation in (item.strip() for item in citation_group.split(",")):
                if citation not in bibliography_keys:
                    raise ValueError(f"equation {entry.eqid!r} has unknown citation {citation!r}")
    for exclusion in exclusions:
        if exclusion.reason not in ALLOWED_EXCLUSIONS:
            raise ValueError(
                f"{exclusion.file}:{exclusion.line}: unsupported exclusion {exclusion.reason!r}"
            )


def _quoted(value: object) -> str:
    return json.dumps(value, ensure_ascii=False)


def render_registry(entries: list[EquationEntry], exclusions: list[SourceExclusion]) -> str:
    """Render deterministic YAML using JSON-compatible scalar quoting."""
    lines = ["schema_version: 1", "equations:"]
    for entry in sorted(entries, key=lambda item: item.eqid):
        lines.extend(
            [
                f"  - eqid: {_quoted(entry.eqid)}",
                f"    section: {_quoted(entry.section)}",
                f"    subsection: {_quoted(entry.subsection)}",
                f"    label: {_quoted(entry.label)}",
                f"    source: {_quoted(entry.source)}",
                f"    status: {_quoted(entry.status)}",
                f"    description: {_quoted(entry.description)}",
                f"    change_note: {_quoted(entry.change_note)}",
                f"    latex: {_quoted(entry.latex)}",
                f"    tex_file: {_quoted(entry.tex_file)}",
                f"    tex_line: {entry.tex_line}",
                "    implementation_refs:",
            ]
        )
        for reference in sorted(
            entry.implementation_refs,
            key=lambda item: (item.file, item.comment_line, item.symbol),
        ):
            lines.extend(
                [
                    f"      - file: {_quoted(reference.file)}",
                    f"        line: {reference.line}",
                    f"        comment_line: {reference.comment_line}",
                    f"        symbol: {_quoted(reference.symbol)}",
                ]
            )
    lines.append("source_exclusions:")
    for exclusion in sorted(exclusions, key=lambda item: (item.file, item.line, item.symbol)):
        lines.extend(
            [
                f"  - file: {_quoted(exclusion.file)}",
                f"    line: {exclusion.line}",
                f"    symbol: {_quoted(exclusion.symbol)}",
                f"    reason: {_quoted(exclusion.reason)}",
            ]
        )
    return "\n".join(lines) + "\n"


def render_markdown(entries: list[EquationEntry], exclusions: list[SourceExclusion]) -> str:
    """Render a searchable human-readable view of the checked registry."""
    lines = [
        "# EOS Equation Registry",
        "",
        "Generated by `scripts/docs/sync_equation_registry.py`; do not edit.",
        "",
    ]
    for entry in sorted(entries, key=lambda item: item.eqid):
        lines.extend(
            [
                f"## `{entry.eqid}`",
                "",
                f"- Section: {entry.section}"
                + (f" / {entry.subsection}" if entry.subsection else ""),
                f"- Status: {entry.status}",
                f"- Source: {entry.source}",
                f"- Description: {entry.description}",
                f"- Change note: {entry.change_note}",
                f"- LaTeX: `{entry.tex_file}:{entry.tex_line}`",
                "- Implementation:",
            ]
        )
        for reference in entry.implementation_refs:
            lines.append(f"  - `{reference.file}:{reference.line}` (`{reference.symbol}`)")
        lines.extend(["", "$$", entry.latex, "$$", ""])
    lines.extend(["## Audited source exclusions", ""])
    for exclusion in sorted(exclusions, key=lambda item: (item.file, item.line, item.symbol)):
        lines.append(
            f"- `{exclusion.file}:{exclusion.line}` (`{exclusion.symbol}`): {exclusion.reason}"
        )
    lines.append("")
    return "\n".join(lines)


def synchronize(repo_root: Path, *, check: bool) -> None:
    """Validate inputs and write or check deterministic generated views."""
    entries, references, exclusions = load_repository_inputs(repo_root)
    validate_registry(
        entries,
        references,
        exclusions,
        load_bibliography_keys(repo_root),
    )
    outputs = {
        repo_root / "docs" / "equations_registry.yaml": render_registry(entries, exclusions),
        repo_root / "docs" / "equations.md": render_markdown(entries, exclusions),
    }
    stale = [
        path
        for path, expected in outputs.items()
        if not path.exists() or path.read_text(encoding="utf-8") != expected
    ]
    if check and stale:
        raise ValueError(
            "generated equation documentation is stale: "
            + ", ".join(str(path.relative_to(repo_root)) for path in stale)
        )
    if not check:
        for path, content in outputs.items():
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(content, encoding="utf-8")


def unclassified_owner_functions(repo_root: Path) -> list[str]:
    """Return unclassified symbols, preserving parse diagnostics as values."""
    unclassified: list[str] = []
    prefix = "unclassified mathematical-owner function: "
    for relative in MATH_OWNER_FILES:
        path = repo_root / relative
        if not path.exists():
            unclassified.append(f"{relative}:missing")
            continue
        try:
            parse_source_annotations(path, repo_root)
        except ValueError as error:
            if str(error).startswith(prefix):
                unclassified.extend(str(error).removeprefix(prefix).split(", "))
                continue
            raise
    return unclassified


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check",
        action="store_true",
        help="check generated views without writing",
    )
    args = parser.parse_args()
    synchronize(Path(__file__).resolve().parents[2], check=args.check)


if __name__ == "__main__":
    main()
