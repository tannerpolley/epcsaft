"""Native Engine calls for the shared #222 cases; Python declares and records."""
import hashlib, json, os, platform, sys, time
IMPORT_START=time.perf_counter_ns()
from pathlib import Path
from epcsaft import Parameters, Pure, Mixture, regression as r, _core
from epcsaft.equilibrium import MoleFractions, Free, Phase, Problem, bubble_point, solve_equilibrium
IMPORT_NS=time.perf_counter_ns()-IMPORT_START
ROOT = Path(__file__).resolve().parents[3]

def document(spec, req):
    case = req['case']
    if case == 'C3':
        d = json.loads((ROOT / spec['packet']).read_text())
        for c in d['components']:
            for q in c['coefficients']:
                if q['family'] == 'debye_huckel_diameter': q['value']['magnitude'] *= .88 if req.get('diameter') == '0.88sigma' else 1
        d['correlations'][0]['constant']['value']['magnitude'] *= req.get('engine_epsilon_factor', 1)
        for p in d['pairs']:
            if {p['component_id_a'], p['component_id_b']} == {'water','chloride-anion'}: p['coefficients'][0]['value']['magnitude'] = req.get('kcl', -.25)
        return d
    if case in ('C1','C5'):
        path = 'data/packets/gross-2001-propane/3/parameters/parameters.json' if case == 'C1' else 'data/packets/gross-vrabec-2006/1/parameters/parameters.json'
        d = json.loads((ROOT / path).read_text()); d['components'] = [c for c in d['components'] if c['component_id'] == spec['components'][0]['id']]
        return d
    d = dict(schema='epcsaft.parameters', schema_version=1, components=[], pairs=[], model_families=[dict(kind=k, choice=v) for k,v in [('base','pc-saft'),('association','general-site'),('electrolyte','none'),('permittivity','none')]], topology=dict(sites=[],edges=[]))
    for c in spec['components']:
        cid = c['id']; d['components'].append(dict(component_id=cid, fixed={k:dict(value=dict(magnitude=v,unit=u)) for k,v,u in [('charge_number',0,'elementary-charge'),('molar_mass',c['molar_mass'],'kilogram / mole')]}, coefficients=[dict(identity=f'component/{cid}/{f}',family=f,value=dict(magnitude=c[key],unit=u)) for key,f,u in [('m','segment_count','dimensionless'),('sigma','segment_diameter','angstrom'),('epsilon','dispersion_energy_over_k','kelvin')]]))
        if not c.get('kappa_ab'): continue
        for site in ('a','b'): d['topology']['sites'].append(dict(component_id=cid,site_id=site,multiplicity=1))
        d['topology']['edges'].append(dict(endpoint_a=dict(component_id=cid,site_id='a'),endpoint_b=dict(component_id=cid,site_id='b'),source=dict(kind='explicit'),energy_over_k=dict(identity=f'assoc/{cid}/energy',value=dict(magnitude=c['epsilon_ab'],unit='kelvin')),volume=dict(identity=f'assoc/{cid}/volume',value=dict(magnitude=c['kappa_ab'],unit='dimensionless'))))
    a,b = [c['id'] for c in spec['components']]
    d['pairs'] = [dict(component_id_a=a,component_id_b=b,coefficients=[dict(identity=f'pair/{a}/{b}/k',family='k_ij',value=dict(magnitude=spec['k'],unit='dimensionless'))])]
    if case == 'C4':
        for first,second in [(a,b),(b,a)]: d['topology']['edges'].append(dict(endpoint_a=dict(component_id=first,site_id='a'),endpoint_b=dict(component_id=second,site_id='b'),source=dict(kind='combining-rule',rule_id='wolbach-sandler',inputs=dict(energy_over_k=[f'assoc/{a}/energy',f'assoc/{b}/energy'],volume=[f'assoc/{a}/volume',f'assoc/{b}/volume']))))
    return d

def convergence(c):
    if c is None: return {}
    return {k:getattr(c,k) for k in ('residual','floor','amplification','contraction','iterations','evaluations')} | {'stop':c.stop.name,'criterion':c.criterion.value if c.criterion else None}

def phase(model,T,rho,x,label):
    s = model.state(T,rho=rho,**({} if len(x)==1 else {'x':x}))
    return dict(phase=label,rho=s.molar_density,x=list(x),pressure=s.pressure,ar=s.residual_helmholtz,Z=s.compressibility_factor,lnphi=s.log_fugacity_coefficient,mu_res=s.residual_chemical_potential_over_rt,pressure_density_derivative=s.pressure_density_derivative,site_fractions=s.site_fractions,contributions={k:getattr(s,k) for k in ('hard_chain','dispersion','association','debye_huckel','born')})

def prepare(inputs,q, *, declare_only=False):
    spec=inputs['models'][q['case']]; ids=[c['id'] for c in spec['components']]; p=Parameters.from_mapping(document(spec,q),components=ids)
    model=(Pure if len(ids)==1 else Mixture)(p); problem=None; ref=None; obs=None
    if q['operation'] in ('saturation','bubble','lle','gradient'):
        if q['case'] in ('C1','C5'): problem=model.saturation(q['T'])
        elif q['case']=='C2': problem=bubble_point(T=q['T'],liquid=Phase('L',kind='liquid'),feed=MoleFractions(dict(zip(ids,q['x']))))
        else: problem=Problem(T=q['T'],P=q['P'],feed=MoleFractions(dict(zip(ids,q['x']))),phases=[Phase('water-rich',kind='liquid',composition_guess=[.99,.01]),Phase('water-poor',kind='liquid',composition_guess=[.01,.99])])
        problem.row_tolerance=1e-12
        if q.get('engine_pressure_guess'): problem.P=Free(q['engine_pressure_guess'])
    if q['operation']=='activity' and not declare_only: ref=model.reference(q['T'],q['P'],solvent_fractions=[1,0,0],ions=(ids[1],ids[2]),stoichiometry=(1,1))
    if q['operation']=='gradient':
        coords=[r.coordinate(p,[(f,target,1)],origin=origin,scale=1) for f,target,origin in ([('segment_count',ids[0],spec['components'][0]['m']),('segment_diameter',ids[0],spec['components'][0]['sigma']),('dispersion_energy_over_k',ids[0],spec['components'][0]['epsilon'])] if q['case']=='C1' else [('k_ij',tuple(ids),spec['k']),('association_energy_over_k',((ids[0],'a'),(ids[0],'b')),spec['components'][0]['epsilon_ab']),('association_volume',((ids[0],'a'),(ids[0],'b')),spec['components'][0]['kappa_ab'])])]
        chart=r.CoordinateSet(); chart.coordinates=coords; model=(Pure if len(ids)==1 else Mixture)(p,active_parameters=chart.active())
        obs=r.observation(p,'pressure',problem,observed=1,form='difference',scale=1,mass_basis=False)
    return model,problem,ref,obs

def call(q,prepared,fresh=False):
    model,problem,ref,obs=prepared; op=q['operation']; T=q['T']
    if op=='gradient':
        value=r.predict(model,obs); gradient=r.differentiate(model,obs)
        if value.value is None or gradient.status != r.ObservationStatus.Available: raise RuntimeError(f'{value.message}; {gradient.message}')
        return dict(P=value.value,gradients=gradient.value,convergence=convergence(value.convergence))
    if problem:
        result=solve_equilibrium(model,problem)
        if not result.success: raise RuntimeError(f'{result.message}; {convergence(result.convergence)}')
        n=len(model.component_ids); xs=[result.mole_fractions[i*n:(i+1)*n] for i in range(result.phase_count)]
        labels=['liquid','vapor'] if op in ('saturation','bubble') else ['water-rich','water-poor']
        order=list(range(result.phase_count)) if op!='lle' else sorted(range(result.phase_count),key=lambda i:-xs[i][0])
        return dict(P=result.pressure,phases=[dict(phase=labels[j],rho=result.molar_densities[i],x=xs[i]) for j,i in enumerate(order)],beta=result.amounts[order[0]]/sum(result.amounts),convergence=convergence(result.convergence),max_residual=max(map(abs,result.residuals)),iterations=result.iterations,evaluations=result.evaluations)
    if op=='activity':
        if fresh: ref=model.reference(T,q['P'],solvent_fractions=[1,0,0],ions=(model.component_ids[1],model.component_ids[2]),stoichiometry=(1,1))
        a=ref.activity(q['m'])
        return dict(gamma=a.mean_ionic_molality_scale,phi=a.osmotic_coefficient,reference_lnphi=list(ref.log_fugacity),pure_lnphi=ref.pure_solvent_log_fugacity,reference_molality=ref.molality,reference_error=ref.convergence_error)
    x=q.get('x',[1,q.get('m',0)*.01801528,q.get('m',0)*.01801528]); x=[v/sum(x) for v in x]
    s=model.state(T,**({'rho':q['rho']} if 'rho' in q else {'P':q['P'],'phase':'liquid'}),**({} if len(x)==1 else {'x':x}))
    return dict(phases=[dict(phase=q.get('phase','liquid'),rho=s.molar_density,x=x,ar=s.residual_helmholtz,Z=s.compressibility_factor,lnphi=s.log_fugacity_coefficient)],density_convergence=convergence(s.density_convergence))

def main():
    started=time.perf_counter_ns(); warm_ns=0; inputs=json.loads(Path(sys.argv[1]).read_text()); records=[];timings=[]; rep=int(sys.argv[3]) if len(sys.argv)>3 else None
    cpu=min(os.sched_getaffinity(0));governor=Path(f'/sys/devices/system/cpu/cpu{cpu}/cpufreq/scaling_governor')
    preparation=json.loads((Path(sys.prefix).resolve().parent/'build/environment-wheel/preparation.json').read_text())
    host=dict(interpreter=sys.executable,preparation=preparation,utc_start=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),system=platform.uname()._asdict(),
              affinity=sorted(os.sched_getaffinity(0)),niceness=os.getpriority(os.PRIO_PROCESS,0),load_start=os.getloadavg(),governor_start=governor.read_text().strip(),
              memory_start=Path('/proc/meminfo').read_text().splitlines()[:3])
    for q in inputs['requests']:
        if rep is not None and not q.get('timing'): continue
        t=time.perf_counter_ns(); prepared=prepare(inputs,q); model_ns=time.perf_counter_ns()-t
        try:
            t=time.perf_counter_ns(); value=call(q,prepared); first_ns=time.perf_counter_ns()-t
            if rep is None:
                if q['operation']=='activity':
                    x=[1,q['m']*.01801528,q['m']*.01801528];x=[v/sum(x) for v in x];s=prepared[0].state(q['T'],P=q['P'],x=x,phase='liquid');value['phases']=[phase(prepared[0],q['T'],s.molar_density,x,'liquid')];value['density_convergence']=convergence(s.density_convergence)
                elif q['operation']!='gradient': value['phases']=[phase(prepared[0],q['T'],p['rho'],p['x'],p['phase']) for p in value['phases']]
                records.append(dict(id=q['id'],case=q['case'],T=q['T'],**value));print(q['id'],flush=True)
            else:
                warm_start=time.perf_counter_ns()
                selections=['reused-reference','fresh-reference'] if q['operation']=='activity' else ['complete-block']
                for selection in selections:
                    fresh=selection=='fresh-reference'
                    for _ in range(5):call(q,prepared,fresh)
                    for sample in range(25):
                        t=time.perf_counter_ns();value=call(q,prepared,fresh);elapsed=time.perf_counter_ns()-t
                        timings.append(dict(id=q['id'],tier=q['timing'],selection=selection,kind='warm',rep=rep,sample=sample,elapsed_ns=elapsed))
                    records.append(dict(id=q['id'],case=q['case'],T=q['T'],selection=selection,**value))
                warm_ns+=time.perf_counter_ns()-warm_start
                for kind,ns in [('model',model_ns),('first-call',first_ns)]:timings.append(dict(id=q['id'],tier=q['timing'],selection='complete-block',kind=kind,rep=rep,sample=0,elapsed_ns=ns))
        except Exception as e:
            records.append(dict(id=q['id'],case=q['case'],error=str(e)));print(q['id'],str(e),flush=True)
    if rep is not None:timings.append(dict(id='Engine-import',tier='import',selection='complete-block',kind='import',rep=rep,sample=0,elapsed_ns=IMPORT_NS))
    driver=dict(path=str(Path(__file__).resolve()),sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    host.update(utc_end=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),load_end=os.getloadavg(),governor_end=governor.read_text().strip())
    Path(sys.argv[2]).write_text(json.dumps(dict(case_driver=driver,host=host,records=records,timings=timings,warm_elapsed_ns=warm_ns,controls=[dict(name=t.name,value=t.value,source=t.source) for t in _core.tolerances()],started_ns=started),indent=2))


def declarations():
    return [
        prepare(d, q, declare_only=True)
        for d in [json.loads((ROOT / 'analyses/2026-competitor-cross-check/inputs/cases.json').read_text())]
        for q in d['requests']
    ]


if __name__=='__main__':main()
