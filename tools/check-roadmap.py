#!/usr/bin/env python3
"""Validate the Markdown progress owner and refresh its derived summary."""
import argparse
import json
from pathlib import Path
import re
import tempfile

ROOT = Path(__file__).resolve().parents[1]
PLAN = ROOT / 'docs/scientific/implementation-plan.md'
START = '<!-- roadmap-summary:start -->'
END = '<!-- roadmap-summary:end -->'
ARCHITECTURE = (ROOT / 'engine/ARCHITECTURE.yaml').read_text()
CATEGORY_LABELS = json.loads(re.search(r'^category_labels: (\{[^\n]+\})$', ARCHITECTURE, re.M).group(1))
LAYERS = ('Core', CATEGORY_LABELS['physics'], CATEGORY_LABELS['model'], 'Equilibrium', 'Regression', 'Surface',
          'Bindings/Python', 'Application acceptance')


def summarize(source, directory):
    rows = {}
    groups = {name: [] for name in LAYERS}
    layer = None
    for line in source.splitlines():
        if line.startswith('### ') and line[4:] in groups:
            layer = line[4:]
        if not re.match(r'^\| [CPMERSBA]\d+ — ', line):
            continue
        cells = [c.strip() for c in line.strip('|').split('|')]
        if len(cells) != 7 or layer is None:
            raise ValueError('malformed sublayer row')
        title, status, remaining, actual, prerequisites, required, record = cells
        key = title.split()[0]
        if key in rows or status not in ('[ ]', '[-]', '[!]', '[x]'):
            raise ValueError(f'{key}: duplicate ID or invalid status')
        if not required:
            raise ValueError(f'{key}: missing acceptance requirements')
        if not re.fullmatch(r'\d+(?:\.\d+)?(?:–\d+(?:\.\d+)?)?', remaining):
            raise ValueError(f'{key}: invalid remaining hours')
        bounds = [float(x) for x in remaining.split('–')]
        low, high = bounds[0], bounds[-1]
        if low > high or (actual != '—' and not re.fullmatch(r'\d+(?:\.\d+)?', actual)):
            raise ValueError(f'{key}: invalid hours')
        deps = re.findall(r'\b[CPMERSBA]\d+\b', prerequisites)
        if key != 'C1' and not deps:
            raise ValueError(f'{key}: missing prerequisites')
        if any(dep not in rows for dep in deps):
            raise ValueError(f'{key}: unknown, cyclic or forward prerequisite')
        if status in ('[-]', '[x]'):
            earlier = [r for name in LAYERS[:LAYERS.index(layer)] for r in groups[name]]
            if any(rows[d][0] != '[x]' for d in deps) or any(r[0] != '[x]' for r in earlier):
                raise ValueError(f'{key}: prerequisites or earlier layer incomplete')
        links = re.findall(r'\[[^\]]+\]\(([^)]+)\)', record)
        if record != '—' and (not links or any(not (directory / link.split('#')[0]).is_file() for link in links)):
            raise ValueError(f'{key}: missing local evidence file')
        if status == '[x]' and (not links or actual == '—' or high != 0):
            raise ValueError(f'{key}: completion requires evidence, actual hours and zero remaining')
        row = (status, low, high, actual)
        rows[key] = row
        groups[layer].append(row)
    if any(not group for group in groups.values()):
        raise ValueError('missing layer rows')
    result = ['| Layer | Qualified progress | Remaining Astra h | Measured actual h |',
              '|---|---|---:|---:|']
    for name, group in groups.items():
        done = sum(r[0] == '[x]' for r in group)
        bar = '█' * (10 * done // len(group)) + '░' * (10 - 10 * done // len(group))
        measured = [float(r[3]) for r in group if r[3] != '—']
        actual = f'{sum(measured):g}' + (' (partial)' if len(measured) < len(group) else '') if measured else '—'
        result.append(f'| {name} | {bar} {100 * done / len(group):.1f}% ({done}/{len(group)}) | '
                      f'{sum(r[1] for r in group):g}–{sum(r[2] for r in group):g} | {actual} |')
    done = sum(r[0] == '[x]' for r in rows.values())
    result += ['', f'Total: **{done}/{len(rows)} qualified sublayers ({100 * done / len(rows):.1f}%)**; '
               f'**{sum(r[1] for r in rows.values()):g}–{sum(r[2] for r in rows.values()):g} remaining Astra hours**.']
    return '\n'.join(result)


def self_test():
    with tempfile.TemporaryDirectory() as tmp:
        directory = Path(tmp)
        (directory / 'evidence.md').write_text('Retained evidence')
        lines = []
        previous = 'Specifications and executable prerequisites'
        for prefix, layer in zip('CPMERSBA', LAYERS):
            key = prefix + '1'
            lines += [f'### {layer}', f'| {key} — Work | [x] | 0 | 1 | {previous} | Verify | [record](evidence.md) |']
            previous = key
        good = '\n'.join(lines)
        assert '8/8' in summarize(good, directory)
        for bad in (good.replace('evidence.md', 'missing.md'),
                    good.replace('| C1 |', '| Z9 |'),
                    good.replace('| [x] |', '| [ ] |', 1),
                    good.replace('| 0 |', '| 1 |', 1)):
            try:
                summarize(bad, directory)
            except ValueError:
                continue
            raise AssertionError('invalid roadmap accepted')
    print('roadmap self-test passed')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true')
    parser.add_argument('--self-test', action='store_true')
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return
    source = PLAN.read_text()
    try:
        summary = summarize(source, PLAN.parent)
        if source.count(START) != 1 or source.count(END) != 1:
            raise ValueError('expected one summary marker pair')
        before, tail = source.split(START)
        _, after = tail.split(END)
        expected = before + START + '\n' + summary + '\n' + END + after
        if args.write:
            PLAN.write_text(expected)
        elif expected != source:
            raise ValueError('stale summary; run python3 tools/check-roadmap.py --write')
    except ValueError as error:
        parser.exit(1, f'roadmap: {error}\n')
    print('roadmap readiness passed (scientific acceptance requires evidence review)')


if __name__ == '__main__':
    main()
