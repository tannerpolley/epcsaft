"""Generate notebook tables and its qualification figure from retained values only."""
import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parents[1]
RESULTS = HERE / 'results'
D = json.loads((HERE / 'inputs/fit-declaration.json').read_text())
FIT = json.loads((RESULTS / 'fit.json').read_text())
E = json.loads((RESULTS / 'evaluation.json').read_text())
S = E['selected']


def read(name):
    with (RESULTS / name).open() as stream:
        return list(csv.DictReader(stream))


def table(name, headers, rows):
    text = '\n'.join(['| ' + ' | '.join(headers) + ' |', '|' + '---|'*len(headers)] +
                     ['| ' + ' | '.join(str(v) for v in row) + ' |' for row in rows]) + '\n'
    (RESULTS / name).write_text(text)
    return text


def number(value, digits=6):
    return '—' if value in ('', None) else f'{float(value):.{digits}f}'


def verdict(value):
    return 'pass' if value in (True, 'True') else 'fail'


def main():
    conditions = D['conditions']
    (RESULTS / 'conditions.md').write_text(
        f"At {conditions['temperature_K']} K and {conditions['pressure_Pa']/1e5:g} bar, "
        f"the native fit used {sum(t['rows'] for t in D['training'])} observations, "
        f"{len(D['coordinates'])} free Li⁺ coordinates and {len(D['starts'])} declared starts. "
        f"The fit took {FIT['wall_seconds']:.2f} s; the fit, Pitzer comparison and a qualification evaluation "
        f"took {E['scientific_wall_seconds']:.2f} s in total on the shared host. "
        f"There were {E['native_failure_count']} unavailable final observations.\n\n"
        f"The declaration froze at `{(HERE / 'inputs/declaration.frozen-at').read_text().strip()}`, "
        f"with SHA-256 `{E['declaration_sha256']}`.\n")
    table('parameters-table.md', ['Selected start', 'Cost', 'σ / Å', 'Packing / Å', 'd_DH / Å', 'd_Born / Å', 'k(Li–water)', 'k(Li–Cl)'],
          [[FIT['selected_start'], f"{S['cost']:.12g}", *[number(v) for v in S['physical']]]])
    table('starts-table.md', ['Start', 'Termination', 'Cost', 'Iterations', 'σ / Å', 'd_DH / Å', 'd_Born / Å', 'k(Li–water)', 'k(Li–Cl)', 'Final native failures'],
          [[s['start'], s['status'], f"{s['final_cost']:.12g}", s['iterations'], *[number(v) for v in s['optimizer']], s['native_failures']] for s in FIT['starts']])
    table('optima-table.md', ['Starts', 'Cost', 'Within selection band', 'Covariance', 'Numerical box'],
          [[', '.join(o['starts']), f"{o['final_cost']:.12g}", o['within_band'], o['covariance_status'],
            'unbounded (bounded value reported)' if o['numerical_box_coordinates'] else 'none'] for o in FIT['optima']])
    limits = D['limits']['Phi']
    buckets = read('buckets.csv')
    table('held-out-table.md', ['Source', 'Bucket', 'Rows', 'Raw Φ bias', 'Raw RMS', 'Raw max', 'Mix RMS', 'Mix max', 'Passing rounding cases', 'Verdict'],
          [[r['source'], r['bucket'], r['rows'], number(r['bias']), number(r['RMS']), number(r['maximum_absolute_error']),
            number(r['mixing_RMS']), number(r['mixing_maximum']), f"{r['passing_realizations']}/{r['realizations']}" if r['scored'] == 'True' else '—',
            verdict(r['all_passed']) if r['scored'] == 'True' else 'report only'] for r in buckets])
    table('rounding-table.md', ['Source', 'Bucket', 'Worst raw RMS', 'Worst raw max', 'Worst mix RMS', 'Worst mix max', 'Frozen RMS / max'],
          [[r['source'], r['bucket'], number(r['worst_RMS']), number(r['worst_maximum']), number(r['worst_mixing_RMS']),
            number(r['worst_mixing_maximum']), f"{limits['RMS']} / {limits['maximum_absolute_error']}"] for r in buckets if r['scored'] == 'True'])
    guard = S['Hamer_Wu']
    g = D['Hamer_Wu_guard']
    table('guard-table.md', ['Source', 'RMS', 'Maximum', 'Frozen RMS / max', 'Verdict'], [
        ['Hamer–Wu, stage zero', number(g['stage0_RMS']), number(g['stage0_maximum']), f"{g['RMS']} / {g['maximum_absolute_error']}", 'pass'],
        ['Hamer–Wu, stage A', number(guard['RMS']), number(guard['maximum_absolute_error']), f"{g['RMS']} / {g['maximum_absolute_error']}", verdict(S['Hamer_Wu_guard_passed'])]])
    table('training-table.md', ['El Guendouzi bucket', 'Rows', 'a_w RMS', 'a_w max', 'Reference-model RMS', 'Reference-model max', 'Within reference residual'],
          [[b, m['rows'], number(m['RMS']), number(m['maximum_absolute_error']), number(D['reference_model_aw_residual'][b]['RMS']),
            number(D['reference_model_aw_residual'][b]['maximum_absolute_error']),
            verdict(m['RMS'] <= D['reference_model_aw_residual'][b]['RMS'] and m['maximum_absolute_error'] <= D['reference_model_aw_residual'][b]['maximum_absolute_error'])]
           for b,m in S['El_Guendouzi_aw'].items()])
    table('trace-table.md', ['NaCl / mol kg⁻¹', 'Trace log γ±(LiCl)', 'Pure LiCl log γ±', 'Pure NaCl log γ±', 'θ_eff / kg mol⁻¹', 'Song, report only', 'F3'],
          [[r['m'], number(r['log_gamma_trace']), number(r['log_gamma_LiCl']), number(r['log_gamma_NaCl']), number(r['theta_eff']),
            number(r['Song_theta_eff']), verdict(abs(float(r['theta_eff'])) <= D['limits']['F3_abs_theta_kg_mol'])] for r in read('trace.csv')])
    divalent = read('trace-divalent.csv')
    song = json.loads((RESULTS / 'trace-divalent-song.json').read_text())
    supplement = json.loads((RESULTS / 'trace-divalent.json').read_text())
    table('trace-divalent-table.md', ['Background', 'Background I / mol kg⁻¹', 'Salt / mol kg⁻¹',
                                    'Stage A log γ±(LiCl)', 'Baseline log γ±(LiCl)', 'Song log γ±(LiCl)'],
          [[r['background'], number(r['I_background_mol_kg'], 1), number(r['m_background_mol_kg']),
            number(r['stage_A_log_gamma_trace']), number(r['baseline_log_gamma_trace']), number(r['Song_log_gamma_trace'])]
           for r in divalent])
    (RESULTS / 'trace-divalent-conditions.md').write_text(
        f"All eight predictions use {conditions['temperature_K']} K, {conditions['pressure_Pa']/1e5:g} bar "
        f"and trace LiCl at {divalent[0]['m_LiCl_mol_kg']} mol/kg water. The table's ionic strength is "
        f"that of the background salt; adding LiCl raises total ionic strength by its trace molality. "
        f"Stage A uses `{supplement['selected_start']}` and baseline uses unchanged `{supplement['base_packet']}`. "
        f"Every activity coefficient is dimensionless on the molality scale; the reported values are base-ten logarithms. "
        f"The largest binary-reference difference across these evaluations is "
        f"{max(abs(float(r[k])) for r in divalent for k in ('baseline_binary_error_decade', 'stage_A_binary_error_decade')):.3g} decade.\n\n"
        f"Song's Li–Mg comparator uses θ = {song['theta_Li_Mg']} kg/mol and "
        f"ψ = {song['psi_Li_Mg_Cl']} kg²/mol², with the paper's binary coefficients. "
        f"It runs natively in Reaktoro {song['reaktoro']}, retaining unequal-charge mixing and native water properties.\n\n"
        f"This report-only supplement was evaluated after fitting at `{supplement['evaluated_at_utc']}`. "
        f"The [evaluation record](results/trace-divalent.json) retains the original declaration and fit digests; "
        f"neither the freeze nor the existing scores and selection was changed.\n")
    lithium = read('lithium-transfers.csv')
    table('lithium-table.md', ['State', 'Check', 'Engine − Pitzer / decade', 'Baseline difference', 'Change / decade', 'Target met', 'F4 preserved'],
          [[r['state'], r['check'], number(r['difference_decade']), number(r['baseline_difference']), number(r['change_from_baseline']),
            verdict(r['target_met']), '—' if r['F4_preserved'] == '' else verdict(r['F4_preserved'])] for r in lithium])
    table('decision-table.md', ['Check', 'Result', 'Meaning'], [
        ['F1 Yao', verdict(S['F1_passed']), 'Both raw and mixing residuals fail; stop for diagnosis'],
        ['F1 Li–Ca and F1′', 'not evaluated', 'Missing measured sources; no qualification claim'],
        ['F2 and F2′', verdict(all(r['target_met'] == 'True' for r in lithium if r['check'] in ('F2','F2-prime'))), 'Lithium model comparisons; do not qualify a packet'],
        ['F3', verdict(S['F3_passed']), 'Bound check; partly in-sample'],
        ['F4', verdict(S['F4_passed']), 'Divalent targets and change limits met'],
        ['F5', verdict(all(r['target_met'] == 'True' for r in lithium if r['check'] == 'F5')), 'KCl model comparison; does not qualify a packet'],
        ['F6 KCl–CaCl₂', verdict(S['F6_passed']), 'Unchanged baseline; K⁺ refit returns to owner'],
        ['F6 KCl–NaCl', 'not evaluated', 'Robinson 1961 unavailable'],
        ['Hamer–Wu guard', verdict(S['Hamer_Wu_guard_passed']), 'Stop for weighting or form diagnosis'],
        ['Stage B trigger', str(S['stage_B_trigger']), 'Yao mixing failure plus training RMS above reference; owner approval required'],
        ['Parameter packet / stage B run', 'neither created nor run', 'Stop conditions preserved']])
    rows = read('held-out-nominal.csv')
    scored = [r for r in rows if r['source'] == 'Yao' and r['bucket'] in ('I<=3','3<I<=8.3')]
    table('mixing-diagnosis-table.md', ['Yao locator', 'I / mol kg⁻¹', 'Mg ion fraction', 'Raw residual', 'End-member part', 'Mixing residual'],
          [[r['locator'], number(r['I'],4), number(3*float(r['m_divalent'])/(2*float(r['m_first_salt'])+3*float(r['m_divalent'])),3),
            number(r['raw_residual']), number(float(r['raw_residual'])-float(r['mixing_residual'])), number(r['mixing_residual'])] for r in scored])
    # A diagnostic figure answers whether subtracting the end-member errors accounts for the Yao failure.
    fig,axes = plt.subplots(1,2,figsize=(10,4),sharey=True,layout='constrained')
    for ax,key,title in zip(axes,('raw_residual','mixing_residual'),('Total mixture error','After subtracting end-member errors'),strict=True):
        ax.axhspan(-limits['maximum_absolute_error'],limits['maximum_absolute_error'],color='#dddddd',label='Frozen maximum-error interval')
        for bucket,color,marker in [('I<=3','#0072B2','o'),('3<I<=8.3','#D55E00','x')]:
            group = [r for r in scored if r['bucket'] == bucket]
            ax.scatter([float(r['I']) for r in group],[float(r[key]) for r in group],color=color,marker=marker,label=bucket)
        ax.axhline(0,color='black',lw=0.7)
        ax.set(xlabel='Ionic strength / mol kg⁻¹',title=title)
        ax.grid(alpha=0.2)
    axes[0].set_ylabel('Model − observed osmotic coefficient')
    axes[0].legend(fontsize=8,loc='upper left')
    (HERE / 'figures').mkdir(exist_ok=True)
    fig.savefig(HERE / 'figures/yao-mixing-residuals.png',dpi=180)
    plt.close(fig)
    (RESULTS / 'verification.md').write_text(
        f"Only the Li⁺ targets changed: {len(json.loads((RESULTS / 'parameter-changes.json').read_text()))} "
        f"physical entries, including tied packing, with {S['unchanged_native_parameter_entries']} other native entries unchanged. "
        f"The maximum change in the nominal KCl–CaCl₂ osmotic coefficients is {S['F6_max_change']:g}. "
        f"The trace reader's maximum difference from the native binary activity is "
        f"{max(abs(float(r['binary_route_error_decade'])) for r in read('trace.csv')):.3g} decade. "
        f"The selected optimum has no numerical-box coordinate; its Born diameter is at the physical floor and covariance is withheld.\n")


if __name__ == '__main__':
    main()
