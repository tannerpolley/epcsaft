"""Read retained comparisons, compute reported statistics, and draw discrete values."""
from collections import defaultdict
import csv
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parents[1]


def rows(name):
    return list(csv.DictReader((HERE / 'results' / name).open()))


def groups(data, keys):
    result = defaultdict(list)
    for row in data:
        result[tuple(row[k] for k in keys)].append(row)
    return result


def summary():
    output = []
    for (figure, component, temperature, reference), data in groups(rows('figure-comparisons.csv'),
            ('figure', 'component', 'temperature_c', 'reference')).items():
        qualified = [r for r in data if r['success'] == 'True' and r['numerical_tolerance_met'] == 'True']
        metric = 'volume_absolute_relative_error_percent' if data[0]['kind'] == 'density' else 'absolute_relative_error_percent'
        output.append({'figure': figure, 'component': component, 'temperature_c': temperature,
            'reference': reference, 'quantity': 'molar volume' if data[0]['kind'] == 'density' else data[0]['kind'],
            'attempted': len(data), 'returned': sum(r['success'] == 'True' for r in data), 'qualified': len(qualified),
            'aad_percent': sum(float(r[metric]) for r in qualified)/len(qualified) if qualified else '',
            'density_aad_percent': sum(float(r['absolute_relative_error_percent']) for r in qualified)/len(qualified)
                if qualified and data[0]['kind'] == 'density' else '',
            'max_absolute_error': max((float(r['absolute_error']) for r in qualified), default=''),
            'model_points_within_tolerance': sum(r['model_tolerance_met'] == 'True' for r in qualified)
                if reference == 'model' else ''})
    return output


def table(data, columns):
    print('| ' + ' | '.join(label for _, label in columns) + ' |')
    print('| ' + ' | '.join('---' for _ in columns) + ' |')
    for row in data:
        values = []
        for key, _ in columns:
            value = row.get(key, '')
            if isinstance(value, str):
                try:
                    value = float(value)
                except ValueError:
                    pass
            if isinstance(value, float):
                value = f'{value:.5g}'
            values.append(str(value).replace('|', '/'))
        print('| ' + ' | '.join(values) + ' |')


def draw(data, path, density=False):
    grouped = groups(data, ('component',) if density else ('figure', 'temperature_c'))
    columns = 2 if density else 3
    nrows = (len(grouped)+columns-1)//columns
    fig, axes = plt.subplots(nrows, columns, figsize=(12, 3.2*nrows), squeeze=False)
    for ax, (key, values) in zip(axes.flat, grouped.items()):
        for reference, label, marker, color in [('model', 'Paper PC-SAFT', 'o', '#0072B2'),
                                               ('experiment', 'Experiment', 's', 'black')]:
            selected = [r for r in values if r['reference'] == reference]
            ax.errorbar([float(r['y'] if density else r['x']) for r in selected],
                [float(r['x'] if density else r['y']) for r in selected],
                xerr=[float(r['y_uncertainty'] if density else r['x_uncertainty']) for r in selected],
                yerr=[float(r['x_uncertainty'] if density else r['y_uncertainty']) for r in selected],
                fmt=marker, markerfacecolor='none', color=color, markersize=4, linewidth=.6, label=label)
        valid = [r for r in values if r['success'] == 'True' and r['numerical_tolerance_met'] == 'True']
        ax.scatter([float(r['temperature_k']) if density else float(r['x']) for r in valid],
                   [float(r['predicted']) for r in valid], marker='x', color='#D55E00', s=22, label='Engine')
        ax.set_title(key[0] if density else f"Figure {key[0]}: {values[0]['component']}, {key[1]} °C" if key[0] != '15' else 'Figure 15: hexane, fixed pressure')
        ax.set_xlabel('Temperature / K' if density else f"Liquid {values[0]['component']} mole fraction")
        ax.set_ylabel('Liquid density / kg m⁻³' if density else 'Temperature / °C' if values[0]['kind'] == 'temperature' else 'Pressure / bar')
        ax.grid(alpha=.15)
    for ax in list(axes.flat)[len(grouped):]:
        ax.set_visible(False)
    axes[0, 0].legend(fontsize=8)
    fig.suptitle('Sparse source points and native PC-SAFT calculations; error bars represent three pixels', fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, .975))
    fig.savefig(path, dpi=160)
    plt.close(fig)


def main():
    data = summary()
    with (HERE / 'results/summary.csv').open('w') as out:
        writer = csv.DictWriter(out, fieldnames=data[0], lineterminator="\n"); writer.writeheader(); writer.writerows(data)
    comparisons = rows('figure-comparisons.csv')
    draw([r for r in comparisons if r['kind'] == 'density'], HERE / 'figures/pure-density-comparison.png', True)
    draw([r for r in comparisons if r['kind'] != 'density'], HERE / 'figures/binary-vle-comparison.png')


if __name__ == '__main__':
    main()
