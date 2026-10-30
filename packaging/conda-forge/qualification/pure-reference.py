from pathlib import Path
import csv,subprocess,sys,json
from decimal import Decimal as D
r=Path(__file__).resolve().parent;out=Path(sys.argv[1])
with (r/'endpoints.txt').open() as inputs,(out/'pure-estimates.csv').open('w') as values:
 subprocess.run([r/'qualification-forward-error.exe',out/'row-floors.csv',out/'packing.csv',out/'snapshot-coordinates.csv'],stdin=inputs,stdout=values,check=True)
ref={float(x['T_K']):[D(x[k]) for k in ['rho_L','rho_V']] for x in csv.DictReader((r/'binary128-roots.csv').open()) if x['method']=='binary128'}
results=[]
for x in csv.DictReader((out/'pure-estimates.csv').open()):
 if x['set']!='854' or x['quantity']!='density':continue
 v=D(x['value']);target=ref[float(x['T'])][int(x['index'])]
 error=abs((v/target).ln());radius=D(x['estimate_chart']);results.append(dict(row=x['row'],arm=x['arm'],temperature=x['T'],phase=x['index'],value=x['value'],reference=str(target),log_error=str(error),radius=str(radius),enclosed=error<=radius))
print(json.dumps(results,indent=2));(out/'pure-reference.json').write_text(json.dumps(results,indent=2))
assert results and all(x['enclosed'] for x in results)
with (out/'pure-complete.csv').open('w') as f:subprocess.run([r/'qualification-pure-complete.exe'],stdout=f,check=True)
complete=[]
for x in csv.DictReader((out/'pure-complete.csv').open()):
 for phase,name in enumerate(['rhoL','rhoV']):
  target=ref[float(x['T'])][phase];error=abs((D(x[name])/target).ln());radius=D(x['radiusL' if phase==0 else 'radiusV'])
  complete.append(dict(T=x['T'],phase=phase,value=x[name],reference=str(target),log_error=str(error),radius=str(radius),enclosed=error<=radius,target_met=error<=D('1e-10'),bits=x['bits'],corrections=x['corrections']))
print(json.dumps(complete,indent=2));(out/'pure-complete.json').write_text(json.dumps(complete,indent=2))
assert len(complete)==4 and all(x['enclosed'] and x['target_met'] for x in complete)
