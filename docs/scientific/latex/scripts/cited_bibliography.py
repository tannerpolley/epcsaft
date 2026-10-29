"""Copy the entries cited in builds/main.aux and builds/supplement.aux from the central Zotero export into references.bib.

Usage: python scripts/cited_bibliography.py [EXPORT] > references.bib
Entries are copied verbatim, except that local attachment ``file`` fields are dropped.
"""
import re
import sys
from pathlib import Path

export = Path(sys.argv[1] if len(sys.argv) > 1 else '~/Zotero/exports/references.bib').expanduser()
builds = Path(__file__).resolve().parents[1] / 'builds'
cited = {k.strip() for aux in ('main.aux', 'supplement.aux') for group in re.findall(r'\\citation\{([^}]*)\}', (builds / aux).read_text()) for k in group.split(',')}
entries = {m.group(1): m.group(0) for m in re.finditer(r'^@\w+\{([^,]+),.*?^\}\n', export.read_text(), re.S | re.M)}
missing = sorted(cited - entries.keys())
if missing:
    sys.exit(f'cited keys absent from {export}: {", ".join(missing)}')
for key in sorted(cited, key=str.lower):
    sys.stdout.write(re.sub(r'^  file = \{.*\},?\n', '', entries[key], flags=re.M) + '\n')
