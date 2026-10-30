from pathlib import Path
import csv,json,sys
from decimal import Decimal as D
out=Path(sys.argv[1]);r=Path(__file__).resolve().parent
ref={x['quantity']:D(x['reference']) for x in csv.DictReader((r/'source-binary128.csv').open()) if x['mask']=='0' and x['sign']=='-1'}
results=[]
for x in csv.DictReader((out/'timing-pure-reference.log').open()):
 if not x['case'].startswith('ethane-source') or x['region']!='reference_error':continue
 for quantity,value,radius in [('C',-D(x['value']),D(x['fixed0'])),('CT',-D(x['coefficient'])/D(x['temperature_coefficient']),D(x['error'])/D(x['temperature_coefficient'])),('p0',D(x['pressure0']),D(x['pressure_floor0'])),('p1',D(x['pressure1']),D(x['pressure_floor1']))]:
  discrepancy=abs(value-ref[quantity]);results.append(dict(case=x['case'],repetition=x['repetition'],quantity=quantity,value=str(value),reference=str(ref[quantity]),discrepancy=str(discrepancy),radius=str(radius),enclosed=discrepancy<=radius))
print(json.dumps(results,indent=2));(out/'source-reference-comparison.json').write_text(json.dumps(results,indent=2));assert results and all(x['enclosed'] for x in results)
