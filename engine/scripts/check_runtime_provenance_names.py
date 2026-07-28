from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
import re
import tomllib


@dataclass(frozen=True, slots=True)
class Violation:
    path: Path
    line: int
    token: str


def _tokens(repo_root: Path) -> tuple[str, ...]:
    policy = (
        repo_root
        / "src"
        / "epcsaft"
        / "data"
        / "bundles"
        / "provenance_names.toml"
    )
    with policy.open("rb") as stream:
        values = tomllib.load(stream).get("tokens")
    if (
        not isinstance(values, list)
        or not values
        or not all(isinstance(value, str) and value for value in values)
    ):
        raise ValueError("provenance name policy must contain nonempty string tokens")
    return tuple(sorted({value.casefold() for value in values}))


def scan_runtime_source(repo_root: Path) -> tuple[Violation, ...]:
    repo_root = repo_root.resolve()
    runtime_root = repo_root / "src" / "epcsaft"
    bundle_root = runtime_root / "data" / "bundles"
    patterns = {
        token: re.compile(
            rf"(?<![A-Za-z0-9]){re.escape(token)}(?![A-Za-z0-9])",
            re.IGNORECASE,
        )
        for token in _tokens(repo_root)
    }
    violations: list[Violation] = []
    for path in sorted(runtime_root.rglob("*")):
        if not path.is_file() or path.is_relative_to(bundle_root):
            continue
        try:
            lines = path.read_text(encoding="utf-8").splitlines()
        except UnicodeDecodeError:
            continue
        for line_number, line in enumerate(lines, start=1):
            for token, pattern in patterns.items():
                if pattern.search(line):
                    violations.append(
                        Violation(path.relative_to(repo_root), line_number, token)
                    )
    return tuple(violations)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, default=Path.cwd())
    arguments = parser.parse_args()
    violations = scan_runtime_source(arguments.repo_root)
    for violation in violations:
        print(
            f"{violation.path}:{violation.line}: "
            f"prohibited provenance name: {violation.token}"
        )
    return 1 if violations else 0


if __name__ == "__main__":
    raise SystemExit(main())
