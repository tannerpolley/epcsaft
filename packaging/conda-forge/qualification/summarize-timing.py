from pathlib import Path
import csv,json,sys,statistics as s
from collections import defaultdict
out=Path(sys.argv[1]);result={}
def stats(values):
 q=s.quantiles(values,n=4,method='inclusive')
 return dict(median=s.median(values),iqr=q[2]-q[0],minimum=min(values),maximum=max(values),samples=len(values))
for mode in ['qualify','frozen']:
 groups=defaultdict(lambda:[0.,0.,0])
 for x in csv.DictReader((out/f'timing-{mode}-events.csv').open()):
  if x['event']!='timing' or int(float(x['value3']))<=0:continue
  key=(x['set'],x['row'],x['arm'],x['value3'],x['reason'])
  a=groups[key];a[0]+=float(x['value1']);a[1]+=float(x['value2']);a[2]+=1
 rows=defaultdict(list)
 for key,(elapsed,cpu,count) in groups.items():
  assert elapsed>=1 and cpu>0 and count>0,(key,elapsed,cpu,count)
  rows[key[:3]+(key[4],)].append(dict(repetition=key[3],calls=count,elapsed=elapsed,process_cpu=cpu,elapsed_per_call=elapsed/count,cpu_per_call=cpu/count))
 result[mode]={str(k):dict(raw=v,elapsed_per_call=stats([x['elapsed_per_call'] for x in v]),cpu_per_call=stats([x['cpu_per_call'] for x in v])) for k,v in rows.items()}
 assert len(rows)==36 and all(len(v)==7 for v in rows.values()),(mode,len(rows))
rows=defaultdict(list);paired=defaultdict(dict)
for x in csv.DictReader((out/'timing-pure-reference.log').open()):
 if int(x['repetition'])<=0:continue
 key=(x['case'],x['region']);calls=int(x['calls']);elapsed=float(x['elapsed']);cpu=float(x['process_cpu'])
 assert calls>0 and elapsed>=1 and cpu>0,(key,x)
 sample=dict(repetition=x['repetition'],calls=calls,elapsed=elapsed,process_cpu=cpu,elapsed_per_call=elapsed/calls,cpu_per_call=cpu/calls)
 rows[key].append(sample)
 if x['region'].startswith('reference_'):paired[(x['case'],x['repetition'])][x['region']]=sample
result['pure_reference']={str(k):dict(raw=v,elapsed_per_call=stats([x['elapsed_per_call'] for x in v]),cpu_per_call=stats([x['cpu_per_call'] for x in v])) for k,v in rows.items()}
assert len(rows)==48 and all(len(v)==7 for v in rows.values()),len(rows)
result['replay_increment']={}
for (case,rep),pair in paired.items():
 a,b=pair['reference_error'],pair['reference_control']
 result['replay_increment'].setdefault(case,[]).append(dict(repetition=rep,elapsed_increment=a['elapsed_per_call']-b['elapsed_per_call'],cpu_increment=a['cpu_per_call']-b['cpu_per_call'],elapsed_ratio=a['elapsed_per_call']/b['elapsed_per_call']))
result['limits']='Complete and frozen boundaries differ. Two-density/replay intervals enclose wider work; unreturned internal work/time is unavailable. Linux desktop times are separate operating conditions. Timing judgment belongs to the lead.'
(out/'timing-summary.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
