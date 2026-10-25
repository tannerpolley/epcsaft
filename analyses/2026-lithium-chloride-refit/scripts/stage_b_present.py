"""Read the retained stage B values and prepare its notebook tables and comparison figure."""
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parents[1]
RESULTS = HERE / 'results/stage-b'


def read(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def table(name, headers, rows):
    (RESULTS / name).write_text('\n'.join(['| ' + ' | '.join(headers) + ' |', '|' + '---|'*len(headers)] +
        ['| ' + ' | '.join(str(v) for v in row) + ' |' for row in rows]) + '\n')


def number(value):
    return '—' if value in ('', None) else f'{float(value):.8g}'


def verdict(value):
    return 'pass' if value in (True, 'True') else 'fail'


def main():
    fit = json.loads((RESULTS / 'fit.json').read_text())
    evaluation = json.loads((RESULTS / 'evaluation.json').read_text())
    declaration = json.loads((HERE / 'inputs/stage-b-declaration.json').read_text())
    s = evaluation['selected']
    a = json.loads((HERE / 'results/evaluation.json').read_text())['selected']
    selected = RESULTS / s['optimum']
    (RESULTS / 'conditions.md').write_text(
        f"The joint fit uses {len(declaration['training_rows'])} training observations, "
        f"{len(declaration['coordinate_order'])} coordinates and {len(declaration['starts'])} starts at "
        f"{declaration['conditions']['temperature_K']} K and {declaration['conditions']['pressure_Pa']/1e5:g} bar. "
        f"The fit took {fit['wall_seconds']:.3f} s; fitting plus all required evaluations took "
        f"{evaluation['scientific_wall_seconds']:.3f} s. The native fit returned {fit['native_failures']} unavailable "
        f"final observations across all starts. There are {len(evaluation['optima'])} distinct converged optima in the cost band.\n\n"
        f"The executable declaration froze at `{declaration['frozen_at_utc']}`, SHA-256 `{fit['declaration_sha256']}`. "
        f"It records the prior inspection of stage A and diagnostic held-out results.\n")
    table('parameters.md', ['Ion', 'σ / Å', 'd_DH / Å', 'd_Born / Å', 'k(ion–water)', 'k(ion–Cl)', 'Packing / Å'],
          [[ion, *map(number, s['optimizer'][5*i:5*i+5]), number(s['physical'][6*i+1])]
           for i, ion in enumerate(('Li⁺', 'Na⁺', 'Mg²⁺'))])
    table('starts.md', ['Start', 'Termination', 'Cost', 'Iterations', 'Rejected trials', 'Final unavailable rows'],
          [[v['start'], v['status'], number(v['final_cost']), v['iterations'], v['trial_failures'], v['native_failures']] for v in fit['starts']])
    best = next(v for v in fit['starts'] if v['start'] == s['starts'][0])
    table('training.md', ['Training block', 'Φ bias', 'Φ RMS', 'Φ maximum'],
          [[key, number(best['training_bias'][key]), number(value), number(best['training_max'][key])]
           for key, value in best['training_rms'].items()])
    table('optima.md', ['Starts', 'Cost', 'Within band', 'Rank', 'Active coordinates', 'Numerical box', 'Covariance'],
          [[', '.join(v['starts']), number(v['final_cost']), v['within_band'], v['rank'],
            ', '.join(declaration['coordinate_order'][i] for i in v['active_bounds']) or 'none',
            ', '.join(declaration['coordinate_order'][i] for i in v['numerical_box_coordinates']) or 'none',
            v['covariance_status']] for v in fit['optima']])
    table('close-optima.md', ['Optimum', 'Cost', '#575 cases passing', 'Yao', 'F3', 'F4', 'Li binary guard', 'Li–Na reference', 'Qualifies'],
          [[v['optimum'], number(v['cost']), f"{v['F575_passing_realizations']}/1467",
            *[verdict(v[k]) for k in ('F1_passed', 'F3_passed', 'F4_passed', 'Hamer_Wu_passed', 'El_Guendouzi_passed', 'qualified')]]
           for v in evaluation['optima']])
    table('575.md', ['Source', 'Role', 'Bucket', 'Observable', 'Nominal RMS', 'Nominal max', 'Amended passing cases', 'Original passing cases', 'Verdict'],
          [[v['source'], 'training/report only' if v['source'] == 'rard-miller-1987' else 'held out',
            v['I_bucket'], v['quantity'], number(v['RMS']), number(v['maximum_absolute_error']),
            f"{v['passing_realizations']}/1467", f"{v['frozen_passing_realizations']}/1467", verdict(v['robust_passed'])]
           for v in s['F575_verdicts']])
    table('yao.md', ['Bucket', 'Rows', 'Raw RMS / max', 'Mixing RMS / max', 'Worst raw RMS / max', 'Passing cases', 'Verdict'],
          [[v['bucket'], v['rows'], f"{number(v['RMS'])} / {number(v['maximum_absolute_error'])}",
            f"{number(v['mixing_RMS'])} / {number(v['mixing_maximum'])}",
            f"{number(v['worst_RMS'])} / {number(v['worst_maximum'])}",
            f"{v['passing_realizations']}/{v['realizations']}" if v['scored'] else '—',
            verdict(v['all_passed']) if v['scored'] else 'report only'] for v in s['F1_Yao']])
    table('guards.md', ['Guard', 'Stage', 'RMS', 'Maximum', 'RMS / max limits', 'Verdict'], [
        ['Hamer–Wu LiCl', 'zero', number(declaration['Hamer_Wu_guard']['stage0_RMS']), number(declaration['Hamer_Wu_guard']['stage0_maximum']), '.01 / .03', 'pass'],
        ['Hamer–Wu LiCl', 'A', number(a['Hamer_Wu']['RMS']), number(a['Hamer_Wu']['maximum_absolute_error']), '.01 / .03', verdict(a['Hamer_Wu_guard_passed'])],
        ['Hamer–Wu LiCl', 'B', number(s['Hamer_Wu']['RMS']), number(s['Hamer_Wu']['maximum_absolute_error']), '.01 / .03', verdict(s['Hamer_Wu_passed'])],
        *[['El Guendouzi a_w ' + b, 'B', number(v['RMS']), number(v['maximum_absolute_error']),
           f"{declaration['reference_model_aw_residual'][b]['RMS']} / {declaration['reference_model_aw_residual'][b]['maximum_absolute_error']}",
           verdict(v['RMS'] <= declaration['reference_model_aw_residual'][b]['RMS'] and v['maximum_absolute_error'] <= declaration['reference_model_aw_residual'][b]['maximum_absolute_error'])]
          for b, v in s['El_Guendouzi_aw'].items()]])
    song_na = read(HERE / 'results/trace.csv')
    table('trace-na.md', ['NaCl / mol kg⁻¹', 'Trace log₁₀ γ±', 'Pure LiCl log₁₀ γ±', 'Pure NaCl log₁₀ γ±', 'θ_eff / kg mol⁻¹', 'Song, report only', 'F3'],
          [[v['m'], *[number(v[k]) for k in ('log_gamma_trace', 'log_gamma_LiCl', 'log_gamma_NaCl', 'theta_eff')],
            number(next(row['Song_theta_eff'] for row in song_na if float(row['m']) == v['m'])), verdict(abs(v['theta_eff']) <= .04)] for v in s['F3']])
    table('lithium.md', ['State', 'Check', 'Engine − Pitzer / decade', '/1 difference / decade', 'Change from /1 / decade', 'Absolute target', 'F4 preservation'],
          [[v['state'], v['check'], number(v['difference_decade']), number(v['baseline_difference']), number(v['change_from_baseline']),
            verdict(v['target_met']), verdict(abs(v['change_from_baseline']) <= .01) if v['check'].startswith('F4') else '—'] for v in s['lithium']])
    table('trace-divalent.md', ['Background', 'Background I / mol kg⁻¹', '/1 log₁₀ γ±(LiCl)', 'Stage A', 'Stage B', 'Song'],
          [[v['background'], v['I_background_mol_kg'], *[number(v[k]) for k in
            ('baseline_log_gamma_trace', 'stage_A_log_gamma_trace', 'stage_B_log_gamma_trace', 'Song_log_gamma_trace')]]
           for v in read(selected / 'trace-divalent.csv')])
    table('K.md', ['Bucket', 'Rows', 'RMS', 'Maximum', 'Worst RMS / max', 'Passing cases', 'Inherited verdict'],
          [[v['bucket'], v['rows'], number(v['RMS']), number(v['maximum_absolute_error']),
            f"{number(v['worst_RMS'])} / {number(v['worst_maximum'])}", f"{v['passing_realizations']}/{v['realizations']}", verdict(v['all_passed'])] for v in evaluation['F6']])
    omission = s['set8_omission']
    (RESULTS / 'mixing.md').write_text(f"Omitting set 8 only as the declared diagnostic leaves {omission['rows']} high-bucket rows, "
        f"with mixing-residual RMS {omission['RMS']:.8g} and maximum {omission['maximum_absolute_error']:.8g}. "
        "The scored row list and raw-residual verdict retain set 8.\n")
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True, layout='constrained')
    stage_a = read(HERE / 'results/held-out-nominal.csv')
    stage_b = read(selected / 'yao-nominal.csv')
    for ax, quantity, title in zip(axes, ('raw_residual', 'mixing_residual'), ('Raw mixture error', 'After endpoint subtraction'), strict=True):
        ax.axhspan(-.03, .03, color='#dddddd', label='Frozen maximum-error interval')
        for rows, label, color, marker in [(stage_a, 'Stage A', '#D55E00', 'o'), (stage_b, 'Stage B', '#0072B2', 'x')]:
            scored = [v for v in rows if v.get('source', 'Yao') == 'Yao' and v['bucket'] in ('I<=3', '3<I<=8.3')]
            ax.scatter([float(v['I']) for v in scored], [float(v[quantity]) for v in scored], label=label, color=color, marker=marker)
        ax.axhline(0, color='black', linewidth=.7)
        ax.set(xlabel='Ionic strength / mol kg⁻¹', title=title)
        ax.grid(alpha=.2)
    axes[0].set_ylabel('Model − observed osmotic coefficient')
    axes[0].legend(fontsize=8)
    fig.savefig(HERE / 'figures/stage-b-yao-residuals.png', dpi=180)
    plt.close(fig)
    sources = [RESULTS / 'fit.json', RESULTS / 'evaluation.json', HERE / 'results/held-out-nominal.csv', selected / 'yao-nominal.csv', selected / 'trace-divalent.csv']
    (RESULTS / 'presentation-inputs.json').write_text(json.dumps({str(p.relative_to(HERE)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}, indent=2) + '\n')


if __name__ == '__main__':
    main()
