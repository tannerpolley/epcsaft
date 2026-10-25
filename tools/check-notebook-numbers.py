#!/usr/bin/env python3
"""Flag unsourced long numeric literals in notebook prose and leaf README tables; never execute cells."""
import argparse
from pathlib import Path
import re
import subprocess

ROOT = Path(__file__).resolve().parents[1]
NUMBER = re.compile(r'10(?:\^\{[+-]?\d+\}|[⁰¹²³⁴⁵⁶⁷⁸⁹⁻⁺]+)|(?<![\w./#])(?:\d{1,3}(?:,\d{3})+|\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+|\s*(?:×|\\times)\s*10(?:\^\{[+-]?\d+\}|[⁰¹²³⁴⁵⁶⁷⁸⁹⁻⁺]+))?(?![\w/]|\.\d)')
SEPARATOR = re.compile(r'^\s*\|?\s*:?-+:?\s*(?:\|\s*:?-+:?\s*)+\|?\s*$')
MASK = re.compile(r'`\{python\}[^`]*`|\{\{<.*?>\}\}|!?(\[[^\]]*\])\([^)]*\)|\b\d{4}-\d{2}-\d{2}\b|\b(?:[\w.-]+/)+[^\s`]+|\b(?=[0-9a-f]*[a-f])[0-9a-f]{8,}\b|\b\d+(?:\.\d+){2,}\b')


def table_rows(text):
    """Yield (line, cells) for each Markdown table header and body row, with or without outer pipes."""
    lines = text.splitlines()
    for at, line in enumerate(lines):
        if SEPARATOR.match(line) and at and '|' in lines[at - 1]:
            end = at + 1
            while end < len(lines) and '|' in lines[end]:
                end += 1
            for number in (at, *range(at + 2, end + 1)):
                yield number, [cell.strip() for cell in lines[number - 1].strip().strip('|').split('|')]


def flags(text):
    fence = None
    yaml = text.startswith('---\n')
    for line, raw in enumerate(text.splitlines(), 1):
        if yaml:
            if line > 1 and raw.strip() in ('---', '...'):
                yaml = False
            continue
        marker = re.match(r'^\s*(`{3,}|~{3,})', raw)
        if marker:
            if fence is None:
                fence = marker[1][0], len(marker[1])
            elif marker[1][0] == fence[0] and len(marker[1]) >= fence[1]:
                fence = None
            continue
        if fence or raw.startswith('    ') or raw.lstrip().startswith(('<!--', '[^')):
            continue
        prose = MASK.sub(lambda m: m[1] or '' if m.lastindex else '', raw)
        for match in NUMBER.finditer(prose):
            value = match[0]
            if re.fullmatch(r'(?:19|20)\d{2}', value):
                continue  # Bibliographic years; other four-digit integers are checked.
            significant = re.sub(r'\D', '', re.split(r'[eE×]|\\times', value)[0]).lstrip('0')
            if len(significant) >= 4 or re.search(r'[eE×^⁻⁺⁰¹²³⁴⁵⁶⁷⁸⁹]|\\times', value):
                yield line, value


def legacy_rows(base, path):
    """Exact text of every table row in path's base version; empty for a file absent at base."""
    shown = subprocess.run(['git', 'show', f'{base}:{path.resolve().relative_to(ROOT)}'], cwd=ROOT, capture_output=True, text=True)
    lines = shown.stdout.splitlines()  # Empty when the file is absent at base.
    return {lines[number - 1] for number, _ in table_rows(shown.stdout)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('notebooks', nargs='*', type=Path)
    parser.add_argument('--all', action='store_true', help='optional full sweep of analysis notebooks')
    parser.add_argument('--counts', action='store_true', help='print only counts per notebook')
    parser.add_argument('--base', help='in a README, check only table rows whose exact text was not a table row at this commit')
    args = parser.parse_args()
    paths = sorted((ROOT / 'analyses').glob('*/notebook.qmd')) + sorted((ROOT / 'analyses/reproduction').glob('*/notebook.qmd')) if args.all else args.notebooks
    total = 0
    for path in paths:
        if not path.is_file():
            continue  # A changed notebook may have been deleted.
        text = path.read_text()
        # A README keeps its declared conditions in prose; only its tables must not retype results.
        rows = {number for number, _ in table_rows(text)}
        if args.base and path.name == 'README.md':
            legacy, lines = legacy_rows(args.base, path), text.splitlines()
            rows = {n for n in rows if lines[n - 1] not in legacy}  # Only rows already tabled at base pass.
        found = [(n, v) for n, v in flags(text) if path.name != 'README.md' or n in rows]
        label = path.relative_to(ROOT) if path.is_relative_to(ROOT) else path
        print(f'{label}: {len(found)} literal numbers')
        if not args.counts:
            for line in sorted({line for line, _ in found}):
                print(f'{label}:{line}: {", ".join(v for n, v in found if n == line)}; compute inline from retained CSV/JSON or generate the table')
        total += len(found)
    return int(total > 0)


if __name__ == '__main__':
    raise SystemExit(main())
