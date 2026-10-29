"""Present retained Wu comparisons; performs no thermodynamic calculation."""
import csv
from pathlib import Path

HERE = Path(__file__).resolve().parents[1] / 'results/issue-902-p12-r1'

def show(filename, columns, destination):
    with (HERE / filename).open() as stream:
        rows = list(csv.DictReader(stream))
    def value(v):
        try:
            return 'NA' if v.lower() == 'nan' else f'{float(v):.7g}'
        except ValueError:
            return v
    text = ['| ' + ' | '.join(columns) + ' |', '| ' + ' | '.join(['---'] * len(columns)) + ' |']
    text += ['| ' + ' | '.join(value(row[c]) for c in columns) + ' |' for row in rows]
    (HERE / destination).write_text('\n'.join(text) + '\n')

show('state-summary.csv', ['state', 'b', 'I', 'observable', 'prediction', 'signed_error', 'B_obs', 'coverage', 'evidence', 'meaning'], 'state-summary.md')
show('acid-guards.csv', ['quantity', 'bucket', 'count', 'mean', 'worst', 'mean_limit', 'worst_limit', 'passes'], 'acid-guards.md')
show('water-mass-spread.csv', ['state', 'observable', 'min_prediction', 'max_prediction'], 'water-mass-spread.md')
