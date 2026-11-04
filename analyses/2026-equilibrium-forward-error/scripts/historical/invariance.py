from r1_common import *
from epcsaft.equilibrium import Reaction
from contextlib import nullcontext
rows=[]
for case in ['cold','newton','reaction','discovery']:
 m=neutral();kwargs={}
 if case=='reaction':m=Mixture(parameters((ETHANE,ETHANE),[('k_ij',(0,1),0.)]));kwargs['reactions']=[Reaction({'s0':-1.,'s1':1.},equilibrium_constant=2.)]
 p=Problem(T=180. if case=='newton' else 220. if case=='discovery' else 350. if case=='reaction' else 300.,P=5e6 if case=='discovery' else 1e6 if case=='newton' else 1e5,feed=Amounts({'s0':.5,'s1':.5}),phases=Discover() if case=='discovery' else [Phase('L',kind='liquid')]+([Phase('V',kind='vapor')] if case=='newton' else []),**kwargs)
 results=[]
 for mode in [None,False,True]:
  with nullcontext() if mode is None else trace(diagnostics=mode) as cap:r=solve_equilibrium(m,p)
  results.append(summary(r))
 assert len({r['values_hash'] for r in results})==1
 assert all(r['work']==results[0]['work'] for r in results)
 rows.append(dict(case=case,modes=results))
save('invariance',rows);print(json.dumps(rows,indent=2))
