"""Copy the retained earlier measurements this campaign compares against into inputs/prior.json, with source hashes.

Sources: the exact-517a3450 key-solve publication run, the #549 history headline (cf26a17b), #533 comparator and
reference-reuse statistics, and the frozen F1 identity. Values are copied, never recomputed.
"""
import csv
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
ROOT = HERE.parents[1]
MAIN = Path('/home/tnnrpolley21/Workspaces/Engineering/ePC-SAFT')
SOURCES = dict(
    key_517a3450=MAIN / 'build/wheels/key-solves-517a3450-20261015T195123/summary.json',
    history_549=ROOT / 'analyses/2026-engine-performance-history/results/headline.csv',
    issue_533_rows=MAIN / 'build/issue-533/selected-row-statistics.csv',
    issue_533_reuse=MAIN / 'build/issue-533/supplemental-call-statistics.csv',
    f1_identity=ROOT / 'tools/fixtures/mea-f1-identity.json',
)


def main():
    key = json.loads(SOURCES['key_517a3450'].read_text())
    prior = dict(sources={k: dict(path=str(p), sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for k, p in SOURCES.items()})
    prior['key_517a3450'] = {case: dict(seconds=v['median_B_seconds'], work=v['work_B']) for case, v in key['cases'].items()}
    prior['history_549'] = {r['operation']: dict(seconds=float(r['today_seconds'])) for r in csv.DictReader(SOURCES['history_549'].open())
                            if r['operation'] in ('C1-300', 'C2-0.5', 'methane-bubble', 'A-propane-liquid')}
    prior['issue_533'] = {r['id']: {k.replace('_median_us', '_seconds'): float(r[k]) * 1e-6 if r[k] else None for k in ('Engine_median_us', 'Clapeyron_median_us', 'FeOs_median_us')}
                          | dict(conditions=r['timing_conditions'])
                          for r in csv.DictReader(SOURCES['issue_533_rows'].open())
                          if r['id'] in ('C1-300', 'C2-0.5', 'C3-0.88sigma-1', 'A-propane-liquid', 'A-propane-vapor')}
    reuse = next(r for r in csv.reader(SOURCES['issue_533_reuse'].open()) if r[:2] == ['C3-0.88sigma-1', 'reused'])
    prior['issue_533']['C3-0.88sigma-1'] |= dict(reused_engine_seconds=float(reuse[3]) * 1e-6)
    frozen = json.loads(SOURCES['f1_identity'].read_text())['frozen']
    prior['f1_frozen'] = dict(revision=frozen['revision'], serial_cpu_seconds=frozen['serial_cpu_seconds'], work=frozen['work'],
                              fits={k: dict(seconds=v['seconds'], iterations=v['iterations']) for k, v in frozen['fits'].items()})
    (HERE / 'inputs/prior.json').write_text(json.dumps(prior, indent=1) + '\n')


if __name__ == '__main__':
    main()
