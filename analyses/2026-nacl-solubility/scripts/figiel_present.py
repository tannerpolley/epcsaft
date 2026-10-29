"""Present retained #902 native rows; never import or execute the Engine."""
import csv
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

STUDY = Path(__file__).resolve().parents[1]
OUT = STUDY/'results/figiel-902'
rows = list(csv.DictReader((OUT/'attempts.tsv').open(), delimiter='\t'))
identity = json.loads((OUT/'run-identity.json').read_text())
valid = [r for r in rows if r['success'] == '1']
central = [r for r in valid if r['feed_mol'] == '12' and r['start_mol_kg'] == '1' and r['dG_J_mol'] == '0' and r['delta_m_mol_kg'] == '0']

def table(name, fields, records):
    (OUT/(name+'.md')).write_text('| '+' | '.join(fields)+' |\n| '+' | '.join(['---']*len(fields))+' |\n'+
        ''.join('| '+' | '.join(r[k] for k in fields)+' |\n' for r in records))

table('comparison', ['P_Pa','m_mol_kg','error_m','rel_error_m','w_salt','error_w','rel_error_w','physical_pass','source_interval_pass'], central)
table('closure', ['P_Pa','feed_mol','start_mol_kg','SI_decade','material_scaled','charge_scaled','classified_max',
                  'affinity_difference_ln','rho_mol_m3','dp_drho','solid_mol','doubled_feed_relative','numerical_pass'],
      [r for r in valid if r['dG_J_mol'] == '0' and r['doubled_feed_relative'] != 'nan'])
table('perturbations', ['P_Pa','T_K','dG_J_mol','m_mol_kg','delta_m_mol_kg','w_salt','source_interval_pass'],
      [r for r in valid if r['start_mol_kg'] == '1' and r['feed_mol'] == '12'])
table('all-closure', ['P_Pa','T_K','dG_J_mol','start_mol_kg','material_scaled','charge_scaled','classified_max',
                       'SI_decade','affinity_difference_ln','numerical_pass'], valid)
(OUT/'failures.json').write_text(json.dumps(identity['failures'], indent=2)+'\n')
fixed = identity['fixed']
table('inputs', ['quantity','value'], [dict(quantity=k,value=str(v)) for k,v in fixed.items()])
fig, axes = plt.subplots(1, 2, figsize=(9,4), layout='constrained')
for ax, row, field, measured, uncertainty, label in [
    (axes[0], central[0], 'm_mol_kg', fixed['measured_m'], fixed['measured_m_se'], 'Molality / mol kg⁻¹'),
    (axes[1], central[1], 'w_salt', fixed['measured_w'], fixed['measured_w_u'], 'Salt mass fraction / kg kg⁻¹')]:
    samples = [r for r in valid if r['P_Pa'] == row['P_Pa'] and r['feed_mol'] == '12' and r['start_mol_kg'] == '1']
    ax.errorbar([0], [measured], yerr=[uncertainty], fmt='o', capsize=5, label='Source')
    ax.scatter([1], [float(row[field])], marker='x', label='Fixed prediction')
    ax.scatter([1.1]*len(samples), [float(r[field]) for r in samples], marker='.', label='Separate ΔG / T shifts')
    ax.set(xticks=[0,1], xticklabels=['Source','Prediction'], ylabel=label, title=row['P_Pa']+' Pa')
    ax.ticklabel_format(axis='y', style='plain', useOffset=False)
    ax.grid(axis='y', alpha=.2)
axes[0].legend(fontsize=8)
fig.savefig(STUDY/'figures/figiel-halite.png', dpi=180)
