from pathlib import Path
import csv,json
from decimal import Decimal as D
out=Path(__file__).resolve().parent;reports=[]
ref={'P':D('6842279.8774235262480592736744'),'rhoL':D('8676.4528556443341908204224'),'rhoV':D('8665.1463344039600528343775')}
for mode in ['qualify','frozen']:
 rows=list(csv.DictReader((out/f'software-{mode}.csv').open()));errors={q:str(max(abs(D(r[q])/v-1) for r in rows)) for q,v in ref.items()}
 enclosed=all(abs((D(r[q])/v).ln())<=D(r['radiusL' if q=='rhoL' else 'radiusV']) if q!='P' else abs(D(r[q])-v)<=D(r['pressure_radius']) for r in rows for q,v in ref.items())
 reports.append(dict(mode=mode,rows=len(rows),successes=sum(int(r['success']) for r in rows),errors=errors,enclosed=enclosed,bits=sorted(set(r['bits'] for r in rows))))
(out/'software-summary.json').write_text(json.dumps(reports,indent=2));print(json.dumps(reports,indent=2));assert all(r['rows']==r['successes']==18 and r['enclosed'] and r['bits']==['64'] and max(D(r['errors'][q]) for q in ['rhoL','rhoV'])<=D('1e-7') for r in reports)

for mode in ["qualify","frozen"]:
 rows=list(csv.DictReader((out/f"software-{mode}.csv").open()))
 assert all(r["requested_tolerance"]=="1" and D(r["classified_residual"])<=D(r["requested_row_tolerance"]) for r in rows)
