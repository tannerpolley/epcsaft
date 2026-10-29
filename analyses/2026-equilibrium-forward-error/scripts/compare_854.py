"""Retained-input preservation and accepted conditioning-aware A/B comparison."""
import argparse,csv,json,struct,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[3]/"engine"))
from epcsaft.equilibrium import (ComparisonEligibility, ComparisonMode, NumericalEstimate, Quantity,
    compare_equilibrium_quantity, project_numerical_estimate)
from collections import defaultdict
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('evidence',type=Path)
parser.add_argument('output',type=Path)
args=parser.parse_args()
out=args.evidence
args.output.mkdir(parents=True,exist_ok=True)
def values(path):
 return {(int(r['case']),r['quantity'],int(r['index'])):float(r['value']) for r in csv.DictReader(path.open())}
a=values(out.parent/'build-r1/states-a/states.csv')
pre=values(out/'states-pre-r2/states.csv');b=values(out/'states-final/states.csv')
da=json.loads((out.parent/'build-r1/states-a/diagnostics.json').read_text())
dp=json.loads((out/'states-pre-r2/diagnostics.json').read_text())
db=json.loads((out/'states-final/diagnostics.json').read_text())
inputs={int(r['case']):r for r in csv.DictReader((out/'states-final/inputs.csv').open())}
regular={int(r['case']) for r in json.loads((out/'regular-preservation-set.json').read_text())['regular_states']}
ae={float(r['T_K']):r for r in csv.DictReader((out/'condition-a.csv').open())}
from epcsaft.equilibrium import Amounts, Free, Phase, Pinned, Problem, compile_problem
from declarations import declarations
model=declarations()
resolved={}
for case,source in inputs.items():
 if case in regular or float(source['temperature_k']) not in ae:continue
 problem=Problem(T=float(source['temperature_k']),P=Free(1e5),feed=Amounts({'s0':1.}),
  phases=[Phase('L',kind='liquid'),Phase('V',kind='vapor',amount=Pinned(0.))])
 compiled=compile_problem(model,problem)
 for arm,states in (('A',a),('B',b)):
  entries=lambda q:[value for (c,k,j),value in sorted(states.items()) if c==case and k==q]
  resolved[case,arm]=compiled.numerical_estimates(entries('coordinates'),float(source['temperature_k']),
   entries('pressure')[0],entries('amounts'),entries('mole_fractions'),entries('molar_densities'),entries('packing_fractions'))
summary=defaultdict(lambda:dict(states=set(),values=0,preservation_values=0,preservation_differences=0,extended_evaluations=0,
 work_differences=0,max_scaled_difference=0.,max_pressure_relative=0.,max_density_relative=0.,failures=0))
changes=[];paired=[]
for case in regular:
 inventory=inputs[case]['inventory'];group=summary[inventory]
 d=db[case].get('linear_solve_diagnostics',{})
 group['extended_evaluations']+=db[case].get('work_spans',{}).get('pure_saturation_extended_evaluate',0)
 group['work_differences']+=dp[case].get('work_spans')!=db[case].get('work_spans')
 for field in ('success','requested_tolerance_met','iterations','evaluations','hessian_calls','message'):
  if dp[case].get(field)!=db[case].get(field):changes.append(dict(case=case,field=field,old=dp[case].get(field),new=db[case].get(field)))
 for key in pre.keys()|b.keys():
  if key[0]!=case:continue
  group['preservation_values']+=1
  if key not in pre or key not in b or struct.pack('>d',pre[key])!=struct.pack('>d',b[key]):
   group['preservation_differences']+=1
for key in sorted(a.keys()|b.keys()):
 case,quantity,index=key;inventory=inputs[case]['inventory'];g=summary[inventory];g['states'].add(case)
 av,bv=a.get(key),b.get(key)
 if av is None or bv is None:
  g['failures']+=1;paired.append(dict(case=case,quantity=quantity,index=index,A=av,B=bv,scaled=None,relative=None,limit=None,passed=False));continue
 scaled=abs(av-bv)/max(abs(av),abs(bv),1.)
 relative=abs(av-bv)/max(abs(av),abs(bv)) if av or bv else 0.
 eligibility=ComparisonEligibility()
 ea,eb=NumericalEstimate(),NumericalEstimate()
 kind=Quantity.Linear
 mode=ComparisonMode.Strict
 allowance=1e-10*max(abs(av),abs(bv)) if quantity in ('pressure','molar_densities') else 1e-10*max(1.,abs(av),abs(bv))
 # Only original reduced log-density coordinates consume the reduced estimate.
 reduced=quantity in ('molar_densities','packing_fractions') or quantity=='coordinates' and index in (2,4)
 d=db[case].get('linear_solve_diagnostics',{})
 source=ae.get(float(inputs[case]['temperature_k']))
 if reduced and (case,'A') in resolved and source and 'forward_error_estimate' in d:
  kind=Quantity.LogPositive if quantity=='coordinates' else Quantity.Positive
  ea=project_numerical_estimate(av,av,float(source['forward_error_estimate']),kind==Quantity.Positive,True)
  eb=project_numerical_estimate(bv,bv,d['forward_error_estimate'],kind==Quantity.Positive,True)
  eligibility.admission_a=da[case]['success'];eligibility.admission_b=db[case]['success']
  eligibility.residual_a=da[case]['requested_tolerance_met'];eligibility.residual_b=db[case]['requested_tolerance_met']
  phase=(index-2)//2 if quantity=='coordinates' else index
  ea.resolved=bool(resolved[case,'A'].densities) and resolved[case,'A'].densities[phase].resolved
  eb.resolved=bool(resolved[case,'B'].densities) and resolved[case,'B'].densities[phase].resolved
  if not ea.resolved:ea.cause=resolved[case,'A'].cause
  if not eb.resolved:eb.cause=resolved[case,'B'].cause
  mode=ComparisonMode.EstimatedError
 result=compare_equilibrium_quantity(kind,av,bv,allowance,ea,eb,eligibility,mode)
 passed=result.judgment!='failed'
 limit=allowance
 g['values']+=1;g['failures']+=not passed;g['max_scaled_difference']=max(g['max_scaled_difference'],scaled)
 if quantity=='pressure':g['max_pressure_relative']=max(g['max_pressure_relative'],relative)
 if quantity=='molar_densities':g['max_density_relative']=max(g['max_density_relative'],relative)
 paired.append(dict(case=case,quantity=quantity,index=index,A=av,B=bv,scaled=scaled,relative=relative,limit=limit,passed=passed,judgment=result.judgment,cause=result.cause,chart_difference=result.chart_difference,radius_sum=result.radius_sum,
  error_a=ea.absolute,error_b=eb.absolute,chart_a=ea.chart,chart_b=eb.chart,available_a=ea.available,available_b=eb.available,
  resolved_a=ea.resolved,resolved_b=eb.resolved,prerequisites={k:getattr(eligibility,k) for k in ('declaration','branch','phase_matching','inventory','residual_a','residual_b','admission_a','admission_b')}))
for g in summary.values():g['states']=len(g['states'])
availability=[dict(case=x['case'],A=x,B=y) for x,y in zip(da,db) if any(x.get(k)!=y.get(k) for k in ('success','requested_tolerance_met','exception'))]
result=dict(inventories=dict(summary),preservation_changes=changes,availability_changes=availability,
 regular_state_count=len(regular),note='The regular input set was fixed before trigger implementation. Near-critical packing/log-packing use only the accepted propagated density allowance.')
(args.output/'equivalence-summary.json').write_text(json.dumps(result,indent=2)+'\n')
with (args.output/'equivalence.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for r in paired for k in r)));w.writeheader();w.writerows(paired)
records=[]
for line in (out/'reduction-lift-3.log').read_text().splitlines():
 if line.startswith('EPCSAFT_COMPARE_JSON '):
  r=json.loads(line.split(' ',1)[1])
  if r.get('quantity')=='pure_saturation_identity':records.append(r)
identity=dict(comparisons=len(records),failures=sum(not x['verdict'] for x in records),
 max_scaled_error=max(float(r['inputs']['error']) for r in records))
(args.output/'identity-summary.json').write_text(json.dumps(identity,indent=2)+'\n')
print(json.dumps(result,indent=2));print(identity)
assert not changes and not availability and not any(g['preservation_differences'] or g['extended_evaluations'] or g['work_differences'] or g['failures'] for g in summary.values())
assert records and identity['failures']==0
