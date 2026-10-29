"""Present retained native values, without calculating thermodynamics or error statistics."""
import csv
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

STUDY = Path(__file__).resolve().parents[1]
OUT = STUDY/'results'
def rows(name):
    return list(csv.DictReader((OUT/(name+'.tsv')).open(), delimiter='\t'))
def table(name, fields, records):
    def value(s):
        try: return format(float(s), '.8g')
        except ValueError: return s
    (OUT/(name+'.md')).write_text('| '+' | '.join(fields)+' |\n| '+' | '.join(['---']*len(fields))+' |\n'+
        ''.join('| '+' | '.join(value(r[k]) for k in fields)+' |\n' for r in records))
density, vapor = rows('densities'), rows('vapors')
central = [r for r in density if r['variant']=='central']
table('stand-in', ['row','T_K','P_Pa','predicted_kg_m3','interpretation'], [r for r in central if r['component']=='dodecane'])
table('density-comparison', ['row','component','T_K','P_Pa','predicted_kg_m3','error_kg_m3','error_percent','engineering_pass','source_interval_pass'], [r for r in central if r['component']!='dodecane'])
table('sensitivities', ['row','variant','T_K','P_Pa','predicted_kg_m3','delta_central_kg_m3','error_percent','source_interval_pass'], [r for r in density if r['variant']!='central'])
table('vapor-comparison', ['row','source','T_K','source_P_Pa','P_Pa','error_percent'], [r for r in vapor if r['start_factor']=='0.5'])
table('vapor-closure', ['row','start_factor','pressure_relative','chemical_RT','material_scaled','start_relative','numerical_pass'], vapor)
table('density-closure', ['row','variant','pressure_relative','dp_drho','branch_mismatch','numerical_pass'], density)
for name in ['readback','screens','closure-maxima','failures']:
    records=rows(name)
    fields=next(csv.reader((OUT/(name+'.tsv')).open(),delimiter='\t'))
    table(name,fields,records)
fig, axes = plt.subplots(1,2,figsize=(10,4),layout='constrained')
for source, color, marker in [('Fuge','#0072B2','o'),('Ambrose','#D55E00','s')]:
    data=[r for r in vapor if r['source']==source and r['start_factor']=='0.5']
    axes[1].scatter([float(r['T_K']) for r in data],[float(r['source_P_Pa']) for r in data],s=25,marker=marker,color=color,label=source+' observed')
    axes[1].scatter([float(r['T_K']) for r in data],[float(r['P_Pa']) for r in data],s=25,marker='x',color=color,label=source+' prediction')
inputs=list(csv.DictReader((STUDY/'inputs/densities.tsv').open(),delimiter='\t',fieldnames=['row','component','variant','T_K','P_Pa','observed','uncertainty']))
for component, color in [('TOP','#0072B2'),('MIBK','#D55E00')]:
    source=[r for r in inputs if r['component']==component and r['variant']=='central']
    observed=[float(r['observed']) for r in source]
    uncertainty=[float(r['uncertainty']) for r in source]
    data=[r for r in central if r['component']==component]
    temperatures=[float(r['T_K']) for r in data]
    axes[0].errorbar(temperatures,observed,yerr=uncertainty,fmt='o',capsize=4,color=color,label=component+' observed')
    axes[0].scatter(temperatures,[float(r['predicted_kg_m3']) for r in data],marker='x',color=color,label=component+' prediction')
axes[0].set(xlabel='Temperature / K',ylabel='Mass density / kg m⁻³',title='Fixed-input density comparison')
axes[1].set(xlabel='Published temperature / K',ylabel='Pressure / Pa',yscale='log',title='MIBK: possibly fitted observations')
for ax in axes:
    ax.legend(fontsize=7)
    ax.grid(alpha=.2)
fig.savefig(STUDY/'figures/comparisons.png',dpi=180)
