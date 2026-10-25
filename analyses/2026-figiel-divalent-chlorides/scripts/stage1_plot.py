"""Plot retained native stage-1 target errors; no model calls."""
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parents[1]
selected=json.loads((HERE/'results/stage1-native-fit.json').read_text())['selected']
rows=list(csv.DictReader((HERE/'results/stage1-native-rows.csv').open()))
minimax=list(csv.DictReader((HERE/'results/stage1-minimax-rows.csv').open()))
fig,axes=plt.subplots(2,3,figsize=(11,6),sharex=True,sharey=True,layout='constrained')
for col,salt in enumerate(('LiCl','NaCl','KCl')):
    for row,quantity in enumerate(('gamma','osmotic')):
        ax=axes[row,col]
        for source,start,label,color,marker in ((rows,-1,'Source packet','#0072B2','o'),(rows,selected,'Native Ceres optimum','#D55E00','x'),(minimax,0,'Native guarded minimax','#009E73','+')):
            data=[v for v in source if int(v['start'])==start and v['group']==f'T {salt} {quantity}']
            ax.scatter([float(v['abscissa']) for v in data],[float(v['error']) for v in data],s=20,c=color,marker=marker,label=label)
        ax.axhspan(-.02,.02,color='0.88',zorder=0,label='Owner tolerance ±0.02 decade')
        ax.set(xlim=(0,4.6),ylim=(-.05,.05))
        if row==0:
            ax.set_title(salt)
        if row==1:
            ax.set_xlabel('Salt molality (mol/kg water)')
        if col==0:
            ax.set_ylabel('log₁₀(model / evaluated)\n'+('γ±' if row==0 else 'Φ')+' error (decade)')
handles,labels=axes[0,0].get_legend_handles_labels()
fig.legend(handles,labels,loc='outside lower center',ncol=4,fontsize=9)
fig.suptitle('Stored guarded minimax passes targets and D9 guards; Ceres passes targets but fails guards',fontsize=12)
fig.savefig(HERE/'figures/stage1-native-chloride-residuals.png',dpi=160)
plt.close(fig)
