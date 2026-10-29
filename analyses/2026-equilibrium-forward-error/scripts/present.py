"""Generate tables and coverage figures from retained native comparisons only."""
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parents[1]
RESULTS, FIGURES = HERE / 'results', HERE / 'figures'


def read(name):
    return list(csv.DictReader((RESULTS / name).open()))


def table(name, headings, rows):
    text = '| ' + ' | '.join(headings) + ' |\n|' + '|'.join('---' for _ in headings) + '|\n'
    text += ''.join('| ' + ' | '.join(map(str, row)) + ' |\n' for row in rows)
    (FIGURES / name).write_text(text)


def number(value):
    return f'{float(value):.10g}'


def main():
    FIGURES.mkdir(exist_ok=True)
    table('solve-observations.md', ['Falsifier', 'Returned numerical behavior', 'Values and work'],
          [[r['falsifier'], r['final_outcome'], r['unchanged_work_or_values']] for r in read('solve-observations.csv')])
    table('solve-gates.md', ['Check', 'Result', 'Measured evidence'],
          [[r['check'], r['result'], r['measured_evidence']] for r in read('solve-gates.csv')])
    table('parameter-fit-determination.md', ['Falsifier', 'Returned numerical behavior', 'Values and work'],
          [[r['falsifier'], r['final_outcome'], r['unchanged_work_or_values']] for r in read('parameter-fit-determination.csv')])
    table('parameter-fit-gates.md', ['Check', 'Result', 'Measured evidence'],
          [[r['check'], r['result'], r['measured_evidence']] for r in read('parameter-fit-gates.csv')])
    coverage = read('coverage-with-endpoints.csv')
    table('coverage.md', ['Evidence', 'T / K', 'Endpoints', 'Comparisons', 'Ordinary', 'Within estimated error', 'Maximum difference / radius sum', 'Failed'],
          [[r['set'], number(r['T']), r['endpoints'], r['comparisons'], r['ordinary'], r['conditional'], number(r['max_ratio']), r['failed']] for r in coverage])
    table('absolute-coverage.md', ['T / K', 'Phase densities', 'Maximum absolute error / own radius', 'Maximum log error / own radius'],
          [[number(r['T']), r['phase_density_comparisons'], number(r['max_absolute_error_over_radius']), number(r['max_log_error_over_radius'])] for r in read('absolute-coverage.csv')])
    table('853-judgments.md', ['Original row', 'T / K', 'Comparisons', 'Ordinary', 'Within estimated error', 'Failed', 'Maximum difference / radius sum'],
          [[r['row'],number(r['T']),r['comparisons'],r['ordinary'],r['within_estimated_error'],r['failed'],number(r['max_difference_over_radius_sum'])] for r in read('853-judgments.csv')])
    source = json.loads((RESULTS / 'source-reference-preservation.json').read_text())
    table('source-reference.md', ['Quantity or check', 'Retained result'], [[k, source[k]] for k in ('cause', 'pressure_available', 'pressure_radius', 'condition', 'precision', 'method', 'row_inner', 'values_identical', 'work_identical', 'ordinary_counters_identical', 'after_diagnostics_identical')])
    binary, consumers, floor, sites = (read('source-'+name+'.csv') for name in ('binary128', 'consumers', 'floor', 'site-errors'))
    groups = []
    for label, rows, key in [('Binary128', binary, 'quantity'), ('Consumer transport', consumers, 'case')]:
        for value in dict.fromkeys(r[key] for r in rows):
            selected = [r for r in rows if r[key] == value]
            groups.append([label, value, len(selected), sum(float(r['ratio']) <= 1 for r in selected), number(max(float(r['ratio']) for r in selected))])
    table('source-coverage.md', ['Check', 'Quantity or case', 'Comparisons', 'Covered', 'Maximum discrepancy / estimate'], groups)
    table('source-pressure-floor.md', ['Quantity', 'With floor', 'Without floor', 'Maximum ratio with', 'Minimum ratio without'],
          [[quantity, sum(float(r['ratio_with']) <= 1 for r in floor if r['quantity'] == quantity),
            sum(float(r['ratio_without']) <= 1 for r in floor if r['quantity'] == quantity),
            number(max(float(r['ratio_with']) for r in floor if r['quantity'] == quantity)),
            number(min(float(r['ratio_without']) for r in floor if r['quantity'] == quantity))] for quantity in ('C', 'CT')])
    derived = read('source-derived.csv')
    table('source-derived.md', ['Case', 'Quantity', 'Discrepancy', 'Estimate', 'Ratio', 'Reciprocal effect on C', 'Born effect on C'],
          [[r['case'], r['quantity'], *[number(r[k]) for k in ('discrepancy', 'estimate', 'ratio', 'reciprocal_effect', 'born_effect')]] for r in derived])
    table('source-errors.md', ['Case', 'T / K', 'Source P / Pa', 'e_C', 'e_CT / K^-1'],
          [[r['case'], *[number(r[k]) for k in ('T', 'P_source', 'natural_C', 'natural_CT')]] for r in read('source-reference-errors.csv')])
    table('source-site-errors.md', ['Path degree', 'Site coefficient degree', 'Signed change', 'Binary128 change', 'Site estimate', 'Complete ratio', 'Ratio without site error'],
          [[r['degree'], r['site_coefficient'], *[number(r[k]) for k in ('perturbation', 'observed_binary128', 'site_estimate', 'ratio', 'without_site_ratio')]] for r in sites])
    table('source-limits.md', ['Scalar comparison', 'Discrepancy', 'Fixed-input estimate', 'Ratio', 'Meaning'],
          [[r['case'], *[number(r[k]) for k in ('error', 'estimate', 'ratio')], r['cause']] for r in read('source-limits.csv')])
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), layout='constrained')
    for ax, rows, title in zip(axes.flat[:2], (binary, consumers), ('Binary128 coefficients', 'Known consumer perturbations')):
        ax.scatter(range(len(rows)), [float(r['ratio']) for r in rows], s=16, marker='o', color='#0072B2')
        ax.set(title=title, ylim=(-.02, 1.08))
    for ax, rows, full, removed, title in (
        (axes[1, 0], floor, 'ratio_with', 'ratio_without', 'Supplied pressure floor'),
        (axes[1, 1], sites, 'ratio', 'without_site_ratio', 'Association coefficient error')):
        ax.scatter(range(len(rows)), [float(r[full]) for r in rows], s=25, marker='o', color='#0072B2', label='Complete estimate')
        ax.scatter(range(len(rows)), [float(r[removed]) for r in rows], s=30, marker='x', color='#D55E00', label='Omit named component')
        ax.set(title=title, yscale='log')
        ax.legend(fontsize=8)
    for ax in axes.flat:
        ax.axhline(1, color='#555555', linestyle='--', linewidth=1)
        ax.set(xlabel='Retained comparison ordinal', ylabel='Discrepancy / estimate (dimensionless)')
        ax.grid(axis='y', alpha=.2)
    fig.savefig(FIGURES/'source-reference-coverage.png', dpi=180)
    plt.close(fig)
    for name in ('source-gates', 'source-performance'):
        rows = read(name+'.csv')
        table(name+'.md', list(rows[0]), [list(r.values()) for r in rows])
    for name in ('falsifiers', 'gates'):
        path=RESULTS / (name + '.csv')
        if path.is_file():
            rows=read(path.name)
            table(name+'.md', list(rows[0]), [list(r.values()) for r in rows])
    table('performance.md', ['Case', 'Arm', 'Serial evaluation CPU / s', 'Full-fit wall / s', 'Reference computations', 'Reference columns', 'Density solves', 'Equilibrium iterations'],
          [[r['case'],r['arm'],number(r['serial_cpu_s']),number(r['fit_wall_s']) if r['fit_wall_s'] else 'not requested',
            r['reference_computations'],r['reference_columns'],r['density_solves'],r['equilibrium_iterations']] for r in read('key-performance.csv')])
    source=json.loads((RESULTS/'853-production-summary.json').read_text())
    table('853-production.md', ['Physical quantities', 'Ordinary', 'Within estimated numerical error', 'Failed'],
          [[source['compared_physical_quantities'],source['ordinary'],source['within_estimated_error'],len(source['failures'])]])
    rows=read('comparisons.csv')
    fig, ax=plt.subplots(figsize=(9,4.6),layout='constrained')
    markers={'840':'s','853':'x','854':'o'}
    colors={'840':'#009E73','853':'#D55E00','854':'#0072B2'}
    for evidence in ('840','853','854'):
        selected=[(i,r) for i,r in enumerate(rows) if r['set']==evidence]
        ax.scatter([i for i,r in selected],[float(r['ratio']) for i,r in selected],s=18,alpha=.75,
                   marker=markers[evidence],color=colors[evidence],label='#'+evidence)
    ax.axhline(1,color='#555555',linestyle='--',linewidth=1,label='Estimated-error coverage limit')
    ax.set(xlabel='Retained scalar comparison ordinal',ylabel='Chart difference / sum of numerical radii (dimensionless)',ylim=(-.025,1.08))
    ax.grid(axis='y',alpha=.2)
    ax.legend(loc='upper left',ncol=2,fontsize=9)
    fig.savefig(FIGURES/'comparison-coverage.png',dpi=180)
    plt.close(fig)


if __name__ == '__main__':
    main()
