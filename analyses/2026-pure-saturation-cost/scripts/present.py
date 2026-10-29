"""Prepare notebook displays exclusively from retained qualification and cost tables."""
import csv
import json
import gzip
from pathlib import Path
import statistics as st

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parents[1]
RESULTS, FIGURES = HERE / 'results', HERE / 'figures'


def read(name):
    return list(csv.DictReader((RESULTS / name).open()))


def number(value):
    return 'unavailable' if value in (None, '', 'None') else f'{float(value):.6g}'


def table(name, header, rows):
    text = '| ' + ' | '.join(header) + ' |\n|' + '|'.join('---' for _ in header) + '|\n'
    text += ''.join('| ' + ' | '.join(map(str, row)) + ' |\n' for row in rows)
    (FIGURES / name).write_text(text)


def dedicated():
    data = RESULTS / 'issue-886'
    saving = json.loads((data / 'saving.json').read_text())
    state = json.loads((data / 'state-comparison.json').read_text())
    boundary = json.loads((data / 'boundary-checks.json').read_text())
    gates = json.loads((data / 'gates.json').read_text())
    with (data / 'paired-ratios.csv').open() as stream:
        pairs = [r for r in csv.DictReader(stream) if r['operation'] == 'prepared' and r['comparison'] == 'B/A']
    with (data / 'statistics.csv').open() as stream:
        costs = [r for r in csv.DictReader(stream) if r['operation'] == 'prepared']
    result = next(r for r in saving if r['operation'] == 'prepared')
    text = [('The saving is resolved.' if result['saving_resolved'] else
             'The saving is unresolved: the strict empirical envelope criterion is not met.'), '',
            'The six paired whole-call observations and same-image controls give:', '',
            '| R | C | −ln R | Saving resolved | Median paired saving / µs |', '|---|---|---|---|---|',
            '| ' + ' | '.join(number(result[k]) for k in ('R','C','negative_log_R')) +
            ' | ' + str(result['saving_resolved']) + ' | ' + number(result['saving_median_us']) + ' |', '',
            '| Arm | Median / µs | Min block median / µs | Max block median / µs |', '|---|---|---|---|']
    text += ['| ' + r['arm'] + ' | ' + ' | '.join(number(r[k]) for k in ('wall_median_us','wall_min_us','wall_max_us')) + ' |' for r in costs]
    with (data / 'ratios.csv').open() as stream:
        comparisons = [r for r in csv.DictReader(stream) if r['operation'] == 'prepared']
    text += ['', '| Fresh comparison | Pairs | Median ratio | Full paired range |', '|---|---|---|---|']
    text += ['| ' + r['comparison'] + ' | ' + r['pairs'] + ' | ' + number(r['median']) +
             ' | ' + number(r['minimum']) + '–' + number(r['maximum']) + ' |' for r in comparisons]
    text += ['', '| Native returned quantity | Value | Unit |', '|---|---|---|',
             f"| Saturation pressure | {boundary['P']:.16g} | Pa |",
             f"| Caller full-row criterion | {boundary['convergence']['criterion']:.16g} | 1 |",
             f"| Full raw residual maximum | {boundary['max_residual']:.16g} | 1 |"]
    text += [f"| {r['phase']} density | {r['rho']:.16g} | mol/m³ |" for r in boundary['phases']]
    text += ['', f"All {state['regular_values']} retained values from {state['regular_states']} regular states agree bit for bit; availability differences: {len(state['availability_differences'])}.", '',
             '| Gate | Retained observation |', '|---|---|']
    text += ['| ' + k + ' | ' + str(v) + ' |' for k,v in gates.items()]
    (FIGURES / 'issue-886.md').write_text('\n'.join(text)+'\n')
    fig, ax = plt.subplots(figsize=(7.2,3.6),layout='constrained')
    ratios = [float(r['ratio']) for r in pairs]
    import math
    ax.plot(range(1,7),ratios,'o',color='#0072B2',label='Paired B/A block medians')
    ax.axhline(math.exp(-result['C']),color='#D55E00',linestyle='--',label='exp(−C): saving criterion')
    ax.axhline(result['R'],color='#555555',linestyle=':',label='Median R')
    ax.set_xticks(range(1,7),[f"r{int(r['repetition'])+1} / {'2' if r['arm'].endswith('control') else '1'}" for r in pairs])
    ax.set_xlabel('Repetition / independent worker pair')
    ax.set_ylabel('Candidate / baseline wall time (1)')
    ax.legend(fontsize=9);ax.grid(axis='y',alpha=.2)
    fig.savefig(FIGURES/'issue-886-saving.png',dpi=180);plt.close(fig)

def main():
    FIGURES.mkdir(exist_ok=True)
    ratios = read('ratios.csv'); stats = read('statistics.csv')
    records=[json.loads(x) for x in gzip.open(RESULTS/'raw/measurement-0.jsonl.gz','rt')]
    blocks_raw=[x for x in records if x['kind']=='block']
    pre=json.loads((RESULTS/'preflight.json').read_text());q=blocks_raw[0]['request']
    cross=json.loads((HERE.parent/'2026-competitor-cross-check/inputs/cases.json').read_text())
    request=next(x for x in cross['requests'] if x['id']==q['case']);component=cross['models'][request['case']]['components'][0]
    table('conditions.md',['Quantity or method','Value'],[
        ['Temperature / K',number(request['T'])],['Segment count',number(component['m'])],
        ['Segment diameter / angstrom',number(component['sigma'])],['Dispersion energy over k / K',number(component['epsilon'])],
        ['Repetitions',len({x['repetition'] for x in blocks_raw})],['Samples per arm and block',q['samples']],['Warmup calls per block',q['warmups']],
        ['Measurement logical CPU',pre['identities']['Engine-A']['affinity'][0]],['Niceness',pre['identities']['Engine-A']['nice']],
        ['Numerical threads per worker',pre['identities']['Clapeyron']['threads']],['FeOs version',pre['identities']['FeOs']['version']],
        ['Clapeyron version',pre['identities']['Clapeyron']['version']],['Julia version',pre['identities']['Clapeyron']['julia']]])

    table('ratios.md', ['Workload','Ratio','Pairs','Median','Minimum','Maximum','Prepared target'],
          [[r['operation'],r['comparison'],r['pairs'],*[number(r[k]) for k in ('median','minimum','maximum','target')]] for r in ratios])
    table('absolute.md',['Workload','Arm','Blocks','Median wall / µs','Min block median / µs','Max block median / µs','Median CPU / µs'],
          [[r['operation'],r['arm'],r['blocks'],*[number(r[k]) for k in ('wall_median_us','wall_min_us','wall_max_us','cpu_median_us')]] for r in stats])
    blocks = read('blocks.csv')
    table('setup.md',['Workload','Arm','Untimed preparation / µs'],[[r['operation'],r['arm'],number(r['prepare_us'])] for r in blocks if r['prepare_us']])
    e = json.loads((RESULTS/'equivalence-summary.json').read_text())['inventories']
    table('equivalence.md',['Inventory','States with rows','Maximum scaled difference','Pressure relative','Density relative','Failed values'],
          [[name,r['states'],*[number(r[k]) for k in ('max_scaled_difference','max_pressure_relative','max_density_relative')],r['failures']] for name,r in e.items()])
    d = read('regular-work-deltas.csv')
    groups = sorted({(r['added_double_corrections'],r['added_double_evaluations']) for r in d})
    table('work.md',['States','Added double corrections per state','Added evaluations per state','Extended evaluations','Extended corrections'],
          [[sum((r['added_double_corrections'],r['added_double_evaluations']) == pair for r in d),*pair,
            sum(int(r['extended_evaluations']) for r in d),sum(int(r['extended_corrections']) for r in d)] for pair in groups])
    n = read('near-critical.csv')
    table('near-critical.md',['T / K','State admitted','Working bits','Extended corrections','Log-density error estimate','Binary128 remaining correction'],
          [[f"{float(r['T_K']):.11g}",r['success'],r['working_bits'] or 'unavailable',r['extended_corrections'] or 'unavailable',number(r['error_estimate']),number(r['binary128_remaining_correction'])] for r in n if float(r['T_K'])>=370])
    a = json.loads((RESULTS/'actions-summary.json').read_text())
    table('actions.md',['Case','Rows','Maximum A1 difference','Maximum A2/A3 difference','Maximum fitting-column difference'],
          [[name,r['rows'],*[number(r[k]) for k in ('max_A1','max_higher','max_fit')]] for name,r in a['cases'].items()])
    ancillary = read('ancillary.csv'); quad = read('ancillary-binary128.csv')
    table('ancillary.md',['Rows','Maximum density relative difference','Maximum pressure relative difference','Binary128 remaining correction at largest segment count'],
          [[len(ancillary),number(max(float(r['density_relative_error']) for r in ancillary)),number(max(float(r['pressure_relative_error']) for r in ancillary)),number(next(r['correction_inf'] for r in quad if r['method']=='B'))]])
    table('gates.md',['Check','Outcome','Retained log'],[[r['check'],r['outcome'],f"`{r['log']}`" if r['log'] else 'awaiting lead'] for r in read('gates.csv')])
    b = json.loads((RESULTS/'runtime-additions.json').read_text())
    table('budget.md',['Layer','Gross added runtime lines','Removed runtime lines'],[[name,r['added'],r['removed']] for name,r in b['layers'].items()]+[['Total',b['total_added'],sum(r['removed'] for r in b['layers'].values())]])
    selected = [r for r in ratios if r['comparison'] in ('B/FeOs','B/Clapeyron')]
    fig, ax = plt.subplots(figsize=(7.2,3.6),layout='constrained')
    for i,r in enumerate(selected):
        mid,lo,hi = (float(r[k]) for k in ('median','minimum','maximum'))
        ax.errorbar(mid,i,xerr=[[mid-lo],[hi-mid]],fmt='o' if r['operation']=='prepared' else 'x',
                    color='#0072B2' if r['operation']=='prepared' else '#D55E00',capsize=4,markersize=7)
    ax.axvline(1,color='#555555',linestyle='--',linewidth=1)
    ax.set_yticks(range(len(selected)),[r['operation']+' '+r['comparison'] for r in selected])
    ax.set_xlabel('Median candidate/comparator wall time ratio (dimensionless)')
    ax.set_title('Prepared-call target at the dashed line; intervals are full paired ranges',fontsize=10)
    ax.grid(axis='x',alpha=.2)
    fig.savefig(FIGURES/'cost.png',dpi=180)
    plt.close(fig)


if __name__ == '__main__':
    main()
    if (RESULTS / 'issue-886/saving.json').exists():
        dedicated()
