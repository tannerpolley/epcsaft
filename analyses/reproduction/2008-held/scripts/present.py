"""Summarize retained results and draw comparisons; never import the Engine."""
import csv
import gzip
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

LEAF = Path(__file__).resolve().parents[1]
SOURCE, OUT = LEAF/'inputs/source', LEAF/'results'
c = pd.read_csv(OUT/'figure-comparisons.csv.gz')
native = [json.loads(line) for line in gzip.decompress((OUT/'native-results.jsonl.gz').read_bytes()).splitlines()]
diagnostics = json.loads((OUT/'diagnostics.json').read_text())
run_id = json.loads((OUT/'run-identity.json').read_text())['run_id']
table4 = pd.read_csv(SOURCE/'held-2008-table-4.csv')
table4['outcome'] = table4.reported.map({'yes':'undetermined', 'no':'not reported'})
table4['cause'] = table4.reported.map({'yes':'original observations, temperatures, pressures and aggregate weighting unavailable', 'no':'blank source cell'})
table4.to_csv(OUT/'table-4-comparisons.csv', index=False)
overall = pd.read_csv(SOURCE/'held-2008-table-4-overall.csv')
overall['aggregate_outcome'] = 'undetermined: original rows and overall weighting unavailable'
overall['count_outcome'] = overall.apply(lambda r: 'reproduced' if r.reported_salt_count == r.count_from_rows else 'not reproduced: source count inconsistency', axis=1)
overall.to_csv(OUT/'table-4-overall-comparisons.csv', index=False)
summary = c.groupby(['quantity', 'label', 'outcome'], sort=False).size().unstack(fill_value=0).reset_index()
summary.to_csv(OUT/'quantity-outcomes.csv', index=False)
c.groupby(['figure', 'panel', 'salt', 'temperature_K', 'label', 'outcome'], sort=False).size().unstack(fill_value=0).reset_index().to_csv(OUT/'series-outcomes.csv', index=False)
base = json.loads((LEAF.parents[2]/'data/packets/held-cameretti-sadowski-2008/4/parameters/parameters.json').read_text())
ions = {a: next(v['value']['magnitude'] for v in component['coefficients'] if v['family']=='segment_diameter')
        for component in base['components'] if component['component_id']!='water' for a in component['aliases']}
structural = []
for r in csv.DictReader((SOURCE/'held-2008-table-3.csv').open()):
    for label, key in [('numerical verification', 'printed_model_diameter_A'), ('physical comparison of effective size', 'measured_distance_A')]:
        value = ions[r['ion']]
        structural.append(dict(table=3, ion=r['ion'], label=label, source_value=float(r[key]), packet_value=value,
                               tolerance=.005, difference=value-float(r[key]), outcome='reproduced' if abs(value-float(r[key]))<=.005 else 'not reproduced',
                               cause='Table 2 coefficient versus Table 3 printed value' if key.startswith('printed') else 'effective fitted size versus diffraction distance; not EOS validation'))
for r in csv.DictReader((SOURCE/'held-2008-table-5.csv').open()):
    structural.append(dict(table=5, ion=r['ion'], label='structural source evidence', source_value=float(r['relative_bare_diameter_percent']),
                           packet_value=None, tolerance=.5, difference=None, outcome='undetermined', cause='bare Pauling relative size is not an EOS prediction'))
pd.DataFrame(structural).to_csv(OUT/'structural-comparisons.csv', index=False)
salts = {r['salt']: r for r in csv.DictReader((SOURCE/'salt-components.csv').open())}
report = dict(execution_run_id=run_id, native_evaluations=len(native), unqualified=sum(not r['numerically_qualified'] for r in native),
              maximum_density_residual_ratio=max(r['density_convergence']['residual_over_criterion'] for r in native if 'density_convergence' in r),
              maximum_equilibrium_residual_ratio=max(r['convergence']['residual_over_criterion'] for r in native if r['quantity']=='vapor_pressure'),
              minimum_liquid_pressure_slope=min(r['pressure_density_derivative'] for r in native if 'pressure_density_derivative' in r),
              maximum_reference_change=max((r['reference']['convergence_error'] for r in native if r.get('reference',{}).get('convergence_error') is not None), default=None),
              maximum_normalization_error=max(r.get('normalization_error',0) for r in native),
              maximum_composition_sum_error=max(abs(sum(r['composition'])-1) for r in native),
              maximum_charge_error=max(abs(r['composition'][1]*int(salts[r['salt']]['nu_anion'])
                                           -r['composition'][2]*int(salts[r['salt']]['nu_cation'])) for r in native),
              table4_counts=table4[table4.reported.eq('yes')].groupby('quantity').size().to_dict(),
              table4_printed_activity_count=106, table4_visible_activity_count=105,
              shared_gamma_error=diagnostics['shared_point']['gamma_error'], shared_osmotic_error=diagnostics['shared_point']['osmotic_error'])
assert report['unqualified']==0 and report['maximum_normalization_error']<=1e-12
assert report['shared_gamma_error']<=1e-12 and report['shared_osmotic_error']<=1e-12
(OUT/'numerical-checks.json').write_text(json.dumps(report, indent=2)+'\n')
probes = pd.read_csv(OUT/'convention-probes.csv')
diagnostics['convention_probes'] = probes.groupby(['figure','panel','salt','quantity','mode','outcome']).size().rename('source_rows').reset_index().to_dict('records')
diagnostics['CsCl_case'] = probes[probes.figure.eq(5) & probes.salt.eq('CsCl') & probes.molality.eq(4.875)].astype(object).where(probes.notna(), None).to_dict('records')
diagnostics['basis_reference_term_audit'] = 'Held 2008 Eqs. (8)-(10): formula-unit molality, stoichiometric ion sum and native molal conversion. Same T/P infinite-dilution ion reference; native terminal changes retained. Screening sigma and packing 0.88 sigma are separate; short-screening varies only the former. Constant 78.09 at 298.15 K; conditional Ascani-Held 2021 water dielectric off temperature. No Born term in Eq. (1). Same-sign ionic dispersion excluded by native family; packet unlike-ion k_ij=1 and water-ion k_ij=0 preserved. Water coefficient uses unchanged native composition and a_w/x_w reporting. No declaration defect established by delivered-r1 audit; no attribution of remaining curve differences to ion pairing.'
(OUT/'diagnostics.json').write_text(json.dumps(diagnostics, indent=2)+'\n')
units = {'density_difference':'Density difference / kg m$^{-3}$', 'vapor_pressure':'Vapor pressure / bar',
         'water_activity_coefficient':'Water activity coefficient / 1', 'mean_ionic_activity_coefficient_molal':'Molal mean ionic coefficient / 1'}
for figure, group in c.groupby('figure', sort=False):
    series = list(group.groupby(['panel','salt','temperature_K'], sort=False))
    fig, axes = plt.subplots(math.ceil(len(series)/3), min(3,len(series)), figsize=(12,3.3*math.ceil(len(series)/3)), squeeze=False, constrained_layout=True)
    for ax, ((panel,salt,T), data) in zip(axes.flat,series):
        line = data[data.source_kind.eq('line')].sort_values('molality')
        marker = data[data.source_kind.eq('markers') & data.qa.eq('pass')]
        ax.plot(line.molality,line.value.where(line.qa.eq('pass')),'--',color='black',lw=1,label='Paper model trace')
        prediction = data.sort_values('molality')
        ax.plot(prediction.molality,prediction.native_value,color='#D55E00',lw=1.5,label='Native primary')
        ax.errorbar(marker.molality,marker.value,xerr=marker.molality_half_width,yerr=marker.value_half_width,
                    fmt='o',mfc='none',color='#0072B2',ms=3,lw=.5,label='Measured symbol; reading interval')
        ax.set(title=f'{salt}, {T:g} K', xlabel='Formula molality / mol kg$^{-1}$ water', ylabel=units[data.quantity.iloc[0]])
        ax.grid(alpha=.15)
    for ax in list(axes.flat)[len(series):]:
        ax.set_visible(False)
    axes.flat[0].legend(fontsize=7)
    fig.savefig(LEAF/f'figures/comparison-{int(figure):02}.png', dpi=140)
    plt.close(fig)
print(json.dumps(report), flush=True)
