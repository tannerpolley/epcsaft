from pathlib import Path
import subprocess,csv,sys,json
from decimal import Decimal as D,getcontext
getcontext().prec=60
r=Path(__file__).resolve().parent;out=Path(sys.argv[1]);ref={x['quantity']:D(x['reference']) for x in csv.DictReader((r/'association-binary128.csv').open())}
with (out/'association.csv').open('w') as f,(out/'association-events.log').open('w') as errors:subprocess.run([r/'qualification-association-reference.exe'],stdout=f,stderr=errors,check=True)
results=[]
for x in csv.DictReader((out/'association.csv').open()):
 a=ref[x['quantity']];value=D(x['value']);scale=D(x['scale']);error=abs(value-a)*scale/max(D(1),abs(a*scale));results.append(dict(quantity=x['quantity'],value=str(value),reference=str(a),normalized_error=str(error),allowance='6.938893903907228377647697926E-18',inner=x['inner'],passed=error<=D('6.938893903907228377647697926E-18')))
print(json.dumps(results,indent=2));(out/'association-comparison.json').write_text(json.dumps(results,indent=2));assert len(results)==9 and all(x['passed'] for x in results)
