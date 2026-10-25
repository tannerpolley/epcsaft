"""Join the before, PR #370 branch and current-main family runs into results/comparison.csv.

Usage: python3 compare.py (paths are fixed below). Reads retained JSONL only; no Engine call.
State metric (frozen before running, #350/#356): max over amounts, mole fractions and molar densities of
|a-b|/max(|a|,|b|,1e-14), each quantity in its own units; T and P are fixed inputs.
"""
import csv
import json
from pathlib import Path

RESULTS = Path(__file__).resolve().parents[1] / 'results'
RUNS = dict(before='family-before-7fa511b6.jsonl', branch='family-pr370-branch.jsonl', after='family-main-6bbaaf1e.jsonl')
WORK = ('iterations', 'evaluations', 'hessian_calls')


def difference(a, b):
    worst, where = 0.0, ''
    for key in ('amounts', 'mole_fractions', 'molar_densities'):
        for i, (x, y) in enumerate(zip(a[key], b[key], strict=True)):
            d = abs(x - y) / max(abs(x), abs(y), 1e-14)
            if d > worst:
                worst, where = d, f'{key}[{i}]'
    return worst, where


def main():
    runs = {k: [json.loads(line) for line in (RESULTS / name).open()] for k, name in RUNS.items()}
    points = [[(r['loading'], r['T'], r['P']) for r in rows] for rows in runs.values()]
    assert points[0] == points[1] == points[2]
    rows = []
    for before, branch, after in zip(*runs.values()):
        both = before['success'] and after['success']
        d, where = difference(before, after) if both else (float('nan'), '')
        same = branch['success'] == after['success'] and all(branch.get(k) == after.get(k) for k in
                                                               WORK + ('amounts', 'mole_fractions', 'molar_densities'))
        rows.append(dict(index=before['index'], loading=before['loading'], T=before['T'], P=before['P'],
                         before_success=before['success'], after_success=after['success'],
                         **{f'before_{k}': before.get(k) for k in WORK}, **{f'after_{k}': after.get(k) for k in WORK},
                         state_difference=d, worst_quantity=where, after_equals_branch=same,
                         before_message=before['message'], after_message=after['message']))
    with (RESULTS / 'comparison.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


if __name__ == '__main__':
    main()
