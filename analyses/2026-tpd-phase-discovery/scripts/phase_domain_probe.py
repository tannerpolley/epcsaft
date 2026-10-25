import sys, json, math, hashlib, importlib.metadata as md, copy
from pathlib import Path
root=Path('/home/tnnrpolley21/.t3/worktrees/ePC-SAFT/t3code-223-build-tpd-stability-test-and-phase-35ebff8f')
sys.path.insert(0,str(root/'analyses/2026-tpd-phase-discovery/scripts'))
import baseline as b
case,model,problem=next(c for c in b.declared_cases() if c[0]=='case-7/flash-5918.5286-Pa')
eq=b.eq
d=md.distribution('epcsaft'); wheel=Path(json.loads(d.read_text('direct_url.json'))['url'].removeprefix('file://'))
sha=hashlib.sha256(wheel.read_bytes()).hexdigest()
assert sha=='384d6080578428a9d182655b7bfd2d23703be0c5b27656501267c9dc96a68759',sha
out={'wheel':str(wheel),'sha256':sha,'case':case,'input':b.problem_input(problem),'reaction_input':b.INPUTS['case-7/reactions']}
r=eq.solve_equilibrium(model,problem)
c=eq.compile_problem(model,problem)
report=eq.stability_test(c,r,eq.StabilityOptions())
t=report.trials[1]
out['reference']={k:getattr(r,k) for k in ['success','objective','amounts','mole_fractions','molar_densities','coordinates','multipliers','residuals','classified_residuals','iterations','evaluations','hessian_calls']}
out['report']={'standard':list(report.standard_over_rt),'plane':list(report.plane_over_rt),'candidate_coordinates':list(report.candidate_coordinates),'present':[dict(p) for p in report.present_phases]}
out['trial']={'seed':t.seed,'kind':str(t.phase.kind),'support':list(t.phase.support),'composition':list(t.composition),'raw_coordinates':list(t.raw_coordinates),'split_coordinates':list(t.split_coordinates),'donor':t.donor,'quantities':dict(t.quantities),'splits':[dict(s) for s in t.split_checks],'message':t.message,'success':t.optimizer.success,'residual':t.optimizer.convergence.residual}
z=[0,0,0,1,-1,-1,-2,1,-1];rt=8.31446261815324*problem.T
s=list(report.standard_over_rt);plane=list(report.plane_over_rt)
def evaluate(label,x,v):
    state=model.state(problem.T,rho=1/v,x=x)
    g=state.residual_helmholtz+math.fsum(xi*(math.log(xi/v)-1+si) for xi,si in zip(x,s) if xi>0)+problem.P*v/rt
    mu=[math.log(xi/v)+mi+si if xi>0 else None for xi,mi,si in zip(x,state.residual_chemical_potential_over_rt,s)]
    row={'label':label,'x':list(x),'v_m3_per_mol':v,'rho_mol_per_m3':1/v,'g_over_rt':g,'tpd_over_rt':g-math.fsum(xi*pi for xi,pi in zip(x,plane)),'pressure_pa':state.pressure,'pressure_defect_pa':state.pressure-problem.P,'dp_drho':state.pressure_density_derivative,'charge':math.fsum(zi*xi for zi,xi in zip(z,x)),'normalization':math.fsum(x)-1,'mu_over_rt':mu,'residual_helmholtz':state.residual_helmholtz,'terms':{k:getattr(state,k) for k in ['hard_chain','dispersion','association','debye_huckel','born','bulk_relative_permittivity']},'site_fractions':list(state.site_fractions)}
    out.setdefault('states',[]).append(row)
    print(label,'P',row['pressure_pa'],'G/RT',g,'TPD',row['tpd_over_rt'],'rho',1/v,flush=True)
    return row
xs=[r.mole_fractions[i*9:(i+1)*9] for i in range(2)]
raw=[evaluate('reference-'+str(i),xs[i],1/r.molar_densities[i]) for i in range(2)]
candidate=report.candidate_coordinates
xc=[list(map(math.exp,candidate[:9])),list(map(math.exp,candidate[11:14]))+[0]*6]
N=[candidate[9],candidate[14]];v=[p['volume_m3_per_mol'] for p in report.present_phases]
present=[evaluate('corrected-reference-'+str(i),xc[i],v[i]) for i in range(2)]
trial=evaluate('trial',t.composition,t.quantities['volume_m3_per_mol'])
w=list(t.composition);beta=t.quantities['split_beta_mol'];donor=t.donor
out['direct_slope']=trial['g_over_rt']-math.fsum(wi*mi for wi,mi in zip(w,present[donor]['mu_over_rt']) if wi>0)+(present[donor]['pressure_pa']-problem.P)*v[donor]/rt
for factor in [1.,.5,.25]:
    amount=beta*factor;remaining=N[donor]-amount
    xd=[(N[donor]*xi-amount*wi)/remaining for xi,wi in zip(xc[donor],w)]
    changed=evaluate('finite-split-donor-'+str(factor),xd,v[donor])
    delta=math.fsum([remaining*changed['g_over_rt'],amount*trial['g_over_rt'],-N[donor]*present[donor]['g_over_rt']])
    inventory=[remaining*xi+amount*wi-N[donor]*old for xi,wi,old in zip(xd,w,xc[donor])]
    out.setdefault('direct_splits',[]).append({'beta_mol':amount,'delta_g_over_rt':delta,'delta_over_beta':delta/amount,'linear_prediction':amount*out['direct_slope'],'inventory_error_mol':inventory})
    print('split',amount,delta,delta/amount,flush=True)
    closed=model.state(problem.T,P=problem.P,x=xd,phase='liquid')
    closed_row=evaluate('finite-split-donor-closed-'+str(factor),xd,1/closed.molar_density)
    out['direct_splits'][-1]['closed_delta_g_over_rt']=math.fsum([remaining*closed_row['g_over_rt'],amount*trial['g_over_rt'],-N[donor]*present[donor]['g_over_rt']])
for label,x in [('trial',w),('present-liquid',xc[0])]:
    for branch in ['liquid','vapor']:
        try:
            state=model.state(problem.T,P=problem.P,x=x,phase=branch)
            row=evaluate(label+'-closed-'+branch,x,1/state.molar_density)
            row.update(root_count=state.stable_root_count,branch=str(state.density_branch),root_stop=str(state.density_convergence.stop),root_residual=state.density_convergence.residual)
        except Exception as error:
            out.setdefault('root_failures',[]).append({'label':label,'branch':branch,'error':str(error)})
            print('root failure',label,branch,str(error),flush=True)
stoichiometry=[[nu.get(species,0.) for species in model.component_ids] for nu in b.INPUTS['case-7/reactions']['stoichiometry']]
out['nu_plane']=[math.fsum(ni*pi for ni,pi in zip(nu,plane)) for nu in stoichiometry]
out['nu_standard']=[math.fsum(ni*si for ni,si in zip(nu,s)) for nu in stoichiometry]
out['nu_affinities']={p['label']:[math.fsum(ni*mi for ni,mi in zip(nu,p['mu_over_rt']) if ni) for nu in stoichiometry] for p in out['states'] if all(mi is not None for mi in p['mu_over_rt'])}
out['trial_summary']=[{'trial':i,'seed':j.seed,'support':list(j.phase.support),'kind':str(j.phase.kind),'resolved_descent':j.resolved_descent,'evidence':j.evidence,'tpd':j.quantities.get('tpd_over_rt'),'volume':j.quantities.get('volume_m3_per_mol')} for i,j in enumerate(report.trials)]
Path('/tmp/issue223_stageb_probe.json').write_text(json.dumps(out,indent=2)+'\n')
print('saved /tmp/issue223_stageb_probe.json',flush=True)
