"""Serial measurement and comparison of native returned outputs (#222)."""
import csv, hashlib, json, math, os, platform, statistics, subprocess, sys, time
from pathlib import Path
HERE=Path(__file__).resolve().parents[1]; ROOT=HERE.parents[1]; CACHE=Path(os.environ.get('XDG_CACHE_HOME',Path.home()/'.cache'))/'epcsaft-comparators/222'; TMP=ROOT/'build/competitor-cross-check'; TMP.mkdir(parents=True,exist_ok=True)
CPU=min(os.sched_getaffinity(0)); ENV=os.environ|dict(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',JULIA_NUM_THREADS='1',JULIA_NUM_PRECOMPILE_TASKS='1',JULIA_DEPOT_PATH=str(CACHE/'depot'))
def write_csv(path,rows):
    if not rows: raise ValueError(f'empty numerical table {path}')
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for r in rows for k in r)),lineterminator='\n');w.writeheader();w.writerows(rows)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def command(lib,output,rep=None):
    if lib=='Engine': cmd=[sys.executable,str(HERE/'scripts/engine.py')]
    elif lib=='FeOs':cmd=[str(CACHE/'feos/bin/python'),str(HERE/'scripts/feos.py')]
    else:cmd=[str(CACHE/'julia-1.13.1/bin/julia'),'--startup-file=no',f'--project={HERE / "inputs"}',str(HERE/'scripts/clapeyron.jl')]
    cmd += [str(HERE/'inputs/cases.json'),str(output)] + ([] if rep is None else [str(rep)])
    base=['timeout','-k','5','1800' if lib=='Clapeyron' else '180','nice','-n','10','taskset','-c',str(CPU),*cmd]
    if lib=='Clapeyron':base=['systemd-run','--user','--quiet','--scope','-p','MemoryMax=6G',*base]
    with output.with_suffix('.log').open('w') as f:subprocess.run(base,env=ENV,stdout=f,stderr=subprocess.STDOUT,check=True)
    return json.loads(output.read_text())
def contract(out):
    cmd=[str(CACHE/'julia-1.13.1/bin/julia'),'--startup-file=no',f'--project={HERE / "inputs"}',str(HERE/'scripts/clapeyron.jl'),'--contract']
    with subprocess.Popen(['systemd-run','--user','--quiet','--scope','-p','MemoryMax=6G','timeout','180','nice','-n','10','taskset','-c',str(CPU),*cmd],env=ENV,stdin=subprocess.PIPE,stdout=subprocess.PIPE,text=True) as helper:
        for r in out['records']:
            if not r.get('helper_payload'):continue
            helper.stdin.write(json.dumps(r['helper_payload'])+'\n');helper.stdin.flush();r.update(json.loads(helper.stdout.readline()))
        helper.stdin.close();helper.wait(timeout=60)
    return out
def retain(inputs,outputs):
    requests={q['id']:q for q in inputs['requests']}; states=[];values=[];diagnostics=[];references=[];failures=[]
    versions={'Engine':'1.0.0','FeOs':'0.10.1','Clapeyron':'0.6.29'}
    def row(lib,q,phase,component,quantity,value,unit,table=values):
        if value is None or isinstance(value,str):return
        table.append(dict(id=q['id'],case=q['case'],operation=q['operation'],implementation=lib,version=versions[lib],diameter=q.get('diameter',''),electrostatic=q.get('electrostatic',''),phase=phase,component=component,quantity=quantity,value=value,unit=unit))
    for lib,out in outputs.items():
        for r in out['records']:
            q=requests.get(r['id'])
            if not q:continue
            if r.get('error'):failures.append(dict(implementation=lib,**r));continue
            for key in ('reference','reference_phases','association_parameter_check','gradient_extra','gradient_native_converged'):
                if key in r:flatten(r['id'],lib,key,r[key],references)
            for p in r.get('phases',[]):
                label=p['phase'];xs=p['x']
                for i,x in enumerate(xs):states.append(dict(id=q['id'],implementation=lib,phase=label,T=q['T'],target_P=q.get('P',r.get('P','')),rho=p['rho'],component=inputs['models'][q['case']]['components'][i]['id'],x=x))
                for i,x in enumerate(xs):row(lib,q,label,str(i),'x',x,'1')
                for k,u in [('rho','mol/m3'),('pressure','Pa'),('ar','1'),('Z','1'),('pressure_density_derivative','Pa m3/mol')]:row(lib,q,label,'',k,p.get(k),u)
                for k in ('lnphi','mu_res','mu_res_over_RT','site_fractions'):
                    for i,v in enumerate(p.get(k) or []):row(lib,q,label,str(i),k,v,'1')
                for k,v in p.get('contributions',{}).items():row(lib,q,label,'',k,v,'1')
                for k,v in p.get('diagnostics',{}).items():row(lib,q,label,'',k,v,'1',diagnostics)
            for k,u in [('P','Pa'),('beta','1'),('gamma','1'),('phi','1'),('reference_molality','mol/kg'),('reference_error','1')]:row(lib,q,'','',k,r.get(k),u)
            for i,v in enumerate(r.get('gradients',[])):row(lib,q,'',str(i),'gradient',v,inputs['gradients'][q['case']][i].split('[')[1][:-1])
            for k in ('convergence','density_convergence'):
                for name,v in r.get(k,{}).items():row(lib,q,'','',f'{k}_{name}',v,'1',diagnostics)
            for k in ('max_residual','iterations','evaluations'):row(lib,q,'','',k,r.get(k),'1',diagnostics)
            for i,v in enumerate(r.get('reference_lnphi',[])):row(lib,q,'',str(i),'reference_lnphi',v,'1',diagnostics)
            ps=r.get('phases',[])
            if len(ps)==2 and all(p.get('lnphi') is not None for p in ps):
                p1,p2=ps;closure=max(abs(math.log(x)+f-math.log(y)-g) for x,f,y,g in zip(p1['x'],p1['lnphi'],p2['x'],p2['lnphi']))
                row(lib,q,'','','chemical_equality',closure,'1',diagnostics);row(lib,q,'','','pressure_equality',abs(p1['pressure']-p2['pressure'])/max(abs(r['P']),1),'1',diagnostics)
                if q['operation']=='lle':row(lib,q,'','','feed_closure',max(abs(z-r['beta']*x-(1-r['beta'])*y) for z,x,y in zip(q['x'],p1['x'],p2['x'])),'1',diagnostics)
            for p in ps:
                row(lib,q,p['phase'],'','normalization',abs(sum(p['x'])-1),'1',diagnostics)
                if q.get('P') and q['operation']!='fixed':row(lib,q,p['phase'],'','pressure_closure',abs(p['pressure']-q['P'])/max(abs(q['P']),1),'1',diagnostics)
                if q['case']=='C3':row(lib,q,p['phase'],'','electroneutrality',abs(p['x'][1]-p['x'][2]),'1',diagnostics)
    def key(v):return tuple(v[k] for k in ('id','phase','component','quantity'))
    reference={key(v):v for v in values if v['implementation']=='Engine'}
    for v in values:
        ref=reference.get(key(v));q=requests[v['id']];qty=v['quantity'];target=None;relative=False
        if q['operation']=='fixed' and qty in ('ar','Z','lnphi','hard_chain','dispersion','association','debye_huckel','born','dipole'):target=1e-11 if v['case']=='C3' else 1e-12
        elif qty in ('rho','P') and q['operation'] not in ('fixed','reference'):target=1e-8;relative=True
        elif qty in ('beta','x') and q['operation'] in ('bubble','lle'):target=1e-8
        elif qty in ('gamma','phi'):target=1e-5
        elif qty=='gradient' and ref:target=1e-6*max(abs(v['value']),abs(ref['value']))+1e-8
        if v['case']=='C3' and q.get('electrostatic')!='matched':target=None
        if qty=='gamma' and ref:difference=math.log(v['value'])-math.log(ref['value']);metric='absolute_log'
        elif ref:difference=v['value']-ref['value'];metric='relative' if relative else 'absolute'
        else:difference=None;metric=''
        divisor=abs(ref['value']) if relative and ref else 1
        v.update(reference_value=None if ref is None else ref['value'],signed_difference=difference,difference=None if difference is None else abs(difference)/divisor,metric=metric,target=target,within_target=None if target is None or difference is None else abs(difference)/divisor<=target)
    for v in diagnostics:
        if v['quantity'] in ('pressure_closure','pressure_equality','chemical_equality','feed_closure'):
            v.update(target=1e-10,within_target=abs(v['value'])<=1e-10)
    write_csv(HERE/'inputs/states.csv',states);write_csv(HERE/'results/values.csv',values);write_csv(HERE/'results/diagnostics.csv',diagnostics)
    write_csv(HERE/'results/reference-and-mapping.csv',references)
    return failures
def flatten(id,producer,path,value,rows):
    """Retain nested comparator diagnostics as one long numerical table."""
    if isinstance(value,dict):
        for k,v in value.items():flatten(id,producer,f'{path}.{k}',v,rows)
    elif isinstance(value,list):
        for i,v in enumerate(value):flatten(id,producer,f'{path}[{i}]',v,rows)
    elif value is not None:
        try: number=float(value)
        except (ValueError,TypeError):return
        if math.isfinite(number):rows.append(dict(id=id,producer=producer,quantity_path=path,value=value))
def timed_values(record,request):
    """The accepted call block's outputs; qualification extras are excluded."""
    operation=request['operation']
    if operation=='gradient': names,phase_names=('P','gradients'),()
    elif operation=='activity': names,phase_names=('gamma','phi'),()
    elif operation=='fixed': names,phase_names=(),('ar','Z','lnphi')
    else: names,phase_names=(('P','beta') if operation=='lle' else ('P',)),('rho','x')
    for name in names:yield '',name,record[name]
    for phase in record.get('phases',[]):
        for name in phase_names:yield phase['phase'],name,phase[name]
def main():
    inputs=json.loads((HERE/'inputs/cases.json').read_text());mode=sys.argv[1] if len(sys.argv)>1 else 'qualify'
    if mode=='qualify':
        outputs={lib:command(lib,TMP/f'{lib}.json') for lib in ['Engine','FeOs','Clapeyron']};outputs['FeOs']=contract(outputs['FeOs']);(TMP/'FeOs.json').write_text(json.dumps(outputs['FeOs'])); failures=retain(inputs,outputs)
        wheel=next((Path(sys.prefix).resolve().parent/'build/environment-wheel').glob('*.whl'))
        identity=dict(base=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),wheel_sha256=sha(wheel),wheel=str(wheel),input_sha256=sha(HERE/'inputs/cases.json'),host=platform.uname()._asdict(),cpu=CPU,hardware=json.loads(subprocess.check_output(['lscpu','-J'],text=True)),governor=Path(f'/sys/devices/system/cpu/cpu{CPU}/cpufreq/scaling_governor').read_text().strip(),memory=Path('/proc/meminfo').read_text().splitlines()[:3],threads={k:ENV[k] for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','JULIA_NUM_THREADS')},failures=failures,absent_diagnostics={lib:out.get('absent_diagnostics',[]) for lib,out in outputs.items()},controls=outputs['Engine'].get('controls',[]),locks={str(p.relative_to(HERE)):sha(p) for p in (HERE/'inputs').glob('*') if p.suffix in ('.lock','.toml')})
        (HERE/'results/summary.json').write_text(json.dumps(identity,indent=2)+'\n')
    elif mode=='identity':
        path=HERE/'results/summary.json'
        identity=json.loads(path.read_text())
        outputs={lib:json.loads((TMP/f'{lib}.json').read_text()) for lib in ('Engine','FeOs','Clapeyron')}
        identity.update(input_sha256=sha(HERE/'inputs/cases.json'),
            reproduction_script_sha256={p.name:sha(p) for p in (HERE/'scripts').glob('*') if p.suffix in ('.py','.jl')},
            input_table_sha256={p.name:sha(p) for p in (HERE/'inputs').glob('*.csv')},
            native_result_counts={lib:len(out['records']) for lib,out in outputs.items()},
            absent_diagnostics={lib:out.get('absent_diagnostics',[]) for lib,out in outputs.items()})
        path.write_text(json.dumps(identity,indent=2)+'\n')
    elif mode in ('timing','timing-engine'):
        samples=[];warm_ns=0;replay=mode=='timing-engine'
        for rep in range(1,4):
            order=['Engine'] if replay else ['Engine','FeOs','Clapeyron']
            for lib in order[rep-1:]+order[:rep-1]:
                if warm_ns>1200e9:raise TimeoutError('warm campaign budget excluding cold compilation')
                out=command(lib,TMP/f'{lib}-timing{"-replay" if replay else ""}-{rep}.json',rep)
                if any(r.get('error') for r in out.get('records',[])):raise RuntimeError(out['records'])
                warm_ns+=out.get('warm_elapsed_ns',0)
                samples.extend(dict(implementation=lib,**s) for s in out['timings'])
                write_csv(HERE/('results/engine-timing-replay.csv' if replay else 'results/timings.csv'),samples)
        write_csv(HERE/('results/engine-timing-replay.csv' if replay else 'results/timings.csv'),samples)
    elif mode=='retain':retain(inputs,{lib:json.loads((TMP/f'{lib}.json').read_text()) for lib in ['Engine','FeOs','Clapeyron']})
    elif mode=='diagnose':
        julia=[str(CACHE/'julia-1.13.1/bin/julia'),'--startup-file=no',f'--project={HERE / "inputs"}']
        prefix=['systemd-run','--user','--quiet','--scope','-p','MemoryMax=6G','timeout','-k','5','1800','nice','-n','10','taskset','-c',str(CPU)]
        selected=inputs|{'requests':[q for q in inputs['requests'] if q['case']=='C3']}
        source_input=TMP/'source-input.json';source_input.write_text(json.dumps(selected))
        subprocess.run(prefix+julia+[str(HERE/'scripts/clapeyron.jl'),'--diagnose',str(source_input),str(TMP/'Clapeyron-source.json')],env=ENV,check=True)
        subprocess.run(prefix+julia+[str(HERE/'scripts/dispersion.jl'),str(HERE/'inputs/cases.json'),str(TMP/'Engine.json'),str(HERE/'inputs/dispersion-coefficients.csv'),str(TMP/'FeOs-source.json')],env=ENV,check=True)
        queries=[];requests={q['id']:q for q in inputs['requests']}
        for lib in ('Engine','FeOs','Clapeyron'):
            for r in json.loads((TMP/f'{lib}.json').read_text())['records']:
                if r['case']!='C3' or r.get('error'):continue
                q=requests[r['id']];cs=inputs['constants'];c=cs[lib];e=cs['Engine']
                er=78.09*q.get('engine_epsilon_factor',1) if lib=='Engine' else 78.09
                if lib!='Engine' and q.get('electrostatic')=='matched':er*= (c['e']/e['e'])**2*e['epsilon0']/c['epsilon0']*e['kB']/c['kB']
                phases=r.get('phases',[])+r.get('reference_phases',[])+([r['reference']] if 'reference' in r else [])
                for i,p in enumerate(phases):
                    d=[v['sigma']*1e-10*(.88 if q['diameter']=='0.88sigma' else 1) for v in inputs['models']['C3']['components'][1:]]
                    queries.append(dict(id=f"{r['id']}:{lib}:phase{i}",T=q['T'],rho=p['rho'],x=p['x'],diameters=d,epsilon_r=p.get('relative_permittivity',er),constants=c))
        result=subprocess.run(prefix+julia+[str(HERE/'scripts/clapeyron.jl'),'--screening'],env=ENV,input=''.join(json.dumps(q)+'\n' for q in queries),text=True,stdout=subprocess.PIPE,check=True)
        (TMP/'screening.json').write_text(json.dumps([json.loads(s) for s in result.stdout.splitlines()]))
        subprocess.run(prefix+julia+[str(HERE/'scripts/clapeyron.jl'),'--ladder',str(HERE/'inputs/cases.json'),str(TMP/'FeOs.json'),str(TMP/'Clapeyron.json'),str(TMP/'ladder.json')],env=ENV,check=True)
        subprocess.run(prefix+julia+[str(HERE/'scripts/clapeyron.jl'),'--association-probe',str(HERE/'inputs/cases.json'),str(TMP/'Clapeyron.json'),str(TMP/'association-probe.json')],env=ENV,check=True)
    elif mode=='source':
        rows=[]
        for name in ('Clapeyron-source','FeOs-source','screening','ladder','association-probe'):
            data=json.loads((TMP/f'{name}.json').read_text())
            if isinstance(data,dict) and any(r.get('error') for r in data['records']):raise RuntimeError('Source diagnosis has unavailable rows; inspect retained native output')
            for r in data.get('records',[]) if isinstance(data,dict) else data:
                if name!='Clapeyron-source' or 'diagnostic' in r or 'corrected_gamma' in r:
                    d=r.get('diagnostic',r); prefix=f"phase{r.get('phase_index','coefficients')}.bits{d.get('precision_bits','')}.step{d.get('refinement','')}"
                    flatten(r['id'],name,prefix,r,rows)
        write_csv(HERE/'results/source-errors.csv',rows)
    elif mode=='compare-source':
        native={lib:{r['id']:r for r in json.loads((TMP/f'{lib}.json').read_text())['records']} for lib in ('Engine','FeOs','Clapeyron')}
        fsource={r['id']:r for r in json.loads((TMP/'FeOs-source.json').read_text())}
        csource={r['id']:r['diagnostic']['correction'] for r in json.loads((TMP/'Clapeyron-source.json').read_text())['records'] if r.get('phase_index')==1 and r['diagnostic']['precision_bits']==512 and r['diagnostic']['refinement']==0}
        rows=[]
        for q in inputs['requests']:
            if q['operation']!='fixed' or (q['case']=='C3' and q['electrostatic']!='matched'):continue
            ep=native['Engine'][q['id']]['phases'][0]
            for lib,source in (('FeOs',fsource),('Clapeyron',csource)):
                record=native[lib].get(q['id']);cor=source.get(q['id'],{})
                if not record or record.get('error'):continue
                if (lib=='FeOs' or q['case']=='C3') and not cor:raise RuntimeError(f'Missing source correction for {lib} {q["id"]}')
                lp=record['phases'][0];target=1e-11 if q['case']=='C3' else 1e-12
                for quantity,field in (('ar','delta_a'),('Z','delta_Z'),('mu_res','delta_mu'),('lnphi','delta_lnphi' if lib=='FeOs' else 'delta_lnphi_fixed')):
                    a,b=lp.get('mu_res_over_RT' if quantity=='mu_res' and lib=='Clapeyron' else quantity),ep.get(quantity)
                    if a is None or b is None:continue
                    correction=cor.get(field,[0.]*len(a) if isinstance(a,list) else 0.)
                    if lib=='Clapeyron' and q['case']=='C5' and q.get('dipole_convention')!='matched':
                        matched=native[lib][q['id']+'-matched']['phases'][0]['mu_res_over_RT' if quantity=='mu_res' else quantity]
                        correction=[x-y for x,y in zip(a,matched)] if isinstance(a,list) else a-matched
                    for i,(av,bv,cv) in enumerate(zip(a,b,correction) if isinstance(a,list) else [(a,b,correction)]):
                        if av is None or bv is None or cv is None:continue
                        remaining=av-bv-cv
                        row=dict(id=q['id'],case=q['case'],implementation=lib,phase=ep['phase'],quantity=quantity,component=i if isinstance(a,list) else '',source_correction=cv,remaining_signed_difference=remaining,target=None if quantity=='mu_res' else target,within_target=None if quantity=='mu_res' else abs(remaining)<=target)
                        if quantity=='lnphi':
                            ml=lp.get('mu_res',lp.get('mu_res_over_RT'));dm=cor.get('delta_mu',[0.]*len(a))[i]
                            if lib=='Clapeyron' and q['case']=='C5' and q.get('dipole_convention')!='matched':dm=ml[i]-native[lib][q['id']+'-matched']['phases'][0]['mu_res_over_RT'][i]
                            chemical=ml[i]-ep['mu_res'][i]-dm
                            logarithmic=dm-cv-math.log(lp['Z']/ep['Z'])
                            row.update(chemical_remainder=chemical,log_Z_remainder=logarithmic,identity_closure=remaining-chemical-logarithmic)
                        rows.append(row)
        write_csv(HERE/'results/source-adjusted.csv',rows)
    elif mode=='decompose':
        native={lib:{r['id']:r for r in json.loads((TMP/f'{lib}.json').read_text())['records']} for lib in ('Engine','FeOs')}
        steps=[('Engine','C3-diagnosis-original','Held packet'),('Engine','C3-diagnosis-diameter','screening diameter'),('Engine','C3-diagnosis-constants','electrostatic constants'),('Engine','C3-diagnosis-interaction','water-Cl interaction'),('FeOs','C3-diagnosis-finite','matched finite reference'),('FeOs','C3-diagnosis-interaction','zero-ion reference')]
        base=native['Engine'][steps[0][1]]['gamma'];previous=None;rows=[]
        for lib,id,label in steps:
            r=native[lib][id];delta=0. if previous is None else r['gamma']-previous['gamma']
            rows.append(dict(id=id,implementation=lib,change=label,delta_gamma=delta,percentage_points=100*delta/base,delta_ln_gamma=0. if previous is None else math.log(r['gamma']/previous['gamma']),delta_phi=0. if previous is None else r['phi']-previous['phi']))
            previous=r
        write_csv(HERE/'results/decomposition.csv',rows)
        first,last=[native['Engine'][id]['phases'][0] for id in ('C3-sigma-1-common-liquid','C3-diameter-at-original-density')]
        contrast=[]
        for quantity in ('hard_chain','dispersion','association','debye_huckel','ar','Z'):
            a=first.get(quantity,first['contributions'].get(quantity));b=last.get(quantity,last['contributions'].get(quantity))
            contrast.append(dict(original_id='C3-sigma-1-common-liquid',variant_id='C3-diameter-at-original-density',quantity=quantity,signed_difference=b-a,unit='1'))
        write_csv(HERE/'results/screening-contrast.csv',contrast)
    elif mode=='timed-check':
        checks=[];requests={q['id']:q for q in inputs['requests']}
        for lib in ('Engine','FeOs','Clapeyron'):
            baseline={r['id']:r for r in json.loads((TMP/f'{lib}.json').read_text())['records']}
            selected={q['id'] for q in inputs['requests'] if q.get('timing') and not
                (lib=='FeOs' and q['case']=='C3' and q.get('diameter')=='sigma') and not
                (lib=='Clapeyron' and q['operation']=='gradient')}
            required={(id,selection) for id in selected for selection in
                (('reused-reference','fresh-reference') if lib=='Engine' and requests[id]['operation']=='activity' else ('complete-block',))}
            for rep in range(1,4):
                out=json.loads((TMP/f'{lib}-timing{"-replay" if lib=="Engine" else ""}-{rep}.json').read_text())
                actual={(r['id'],r.get('selection','complete-block')) for r in out['records']}
                if actual!=required or len(out['records'])!=len(required):raise RuntimeError(f'{lib} repetition {rep}: missing or duplicated timed requests')
                for r in out['records']:
                    if r.get('error'):raise RuntimeError(r)
                    expected={(phase,name):value for phase,name,value in timed_values(baseline[r['id']],requests[r['id']])}
                    returned=list(timed_values(r,requests[r['id']]))
                    actual={(phase,name):value for phase,name,value in returned}
                    if len(actual)!=len(returned) or actual.keys()!=expected.keys():raise RuntimeError(f'Missing timed field {lib} {r["id"]}')
                    for (phase,name),value in actual.items():
                        reference=expected[phase,name]
                        pairs=zip(value,reference,strict=True) if isinstance(value,list) else [(value,reference)]
                        for component,(v,e) in enumerate(pairs):
                            if v is None or e is None:raise RuntimeError(f'Unavailable timed output {lib} {r["id"]} {name}')
                            difference=math.log(v/e) if name=='gamma' else v-e
                            metric=abs(difference)/abs(e) if name in ('P','rho') else abs(difference)
                            target=1e-8 if name in ('P','rho','x','beta') else 1e-5 if name in ('gamma','phi') else 1e-11 if r['case']=='C3' else 1e-12
                            if name=='gradients':target=1e-6*max(abs(v),abs(e))+1e-8
                            checks.append(dict(id=r['id'],implementation=lib,campaign='Engine replay' if lib=='Engine' else 'original comparator',selection=r.get('selection','complete-block'),rep=rep,phase=phase,quantity=name,component=component,signed_difference=difference,difference=metric,target=target,within_target=metric<=target))
        write_csv(HERE/'results/timed-output-checks.csv',checks)
    elif mode in ('statistics','statistics-engine'):
        replay=mode=='statistics-engine';groups={};rows=list(csv.DictReader((HERE/('results/engine-timing-replay.csv' if replay else 'results/timings.csv')).open()))
        for s in rows:
            kind=s.get('kind') or s.get('warm');kind='warm' if kind=='sample' else kind
            key=(s['id'],s['implementation'],s['selection'],kind,s['rep']);groups.setdefault(key,[]).append(float(s['elapsed_ns']))
        stats=[]
        for key,ns in groups.items():
            p=statistics.quantiles(ns,n=10,method='inclusive') if len(ns)>1 else [ns[0]]*9
            stats.append(dict(zip(('id','implementation','selection','kind','rep'),key))|dict(samples=len(ns),median_ns=statistics.median(ns),p10_ns=p[0],p90_ns=p[8]))
        write_csv(HERE/('results/engine-replay-statistics.csv' if replay else 'results/timing-statistics.csv'),stats)
        budgets=[]
        for rep in range(1,4):
            for lib in (('Engine',) if replay else ('Engine','FeOs','Clapeyron')):
                out=json.loads((TMP/f'{lib}-timing{"-replay" if replay else ""}-{rep}.json').read_text())
                warm=[s['elapsed_ns'] for s in out['timings'] if (s.get('kind') or s.get('warm')) in ('warm','sample')]
                cold=[s['elapsed_ns'] for s in out['timings'] if (s.get('kind') or s.get('warm')) not in ('warm','sample')]
                budgets.append(dict(implementation=lib,rep=rep,warm_elapsed_ns=out['warm_elapsed_ns'],cold_measured_ns=sum(cold),max_sample_ns=max(warm),samples=len(warm),within_case_cap=max(warm)<=180e9))
        write_csv(HERE/('results/engine-replay-budgets.csv' if replay else 'results/budgets.csv'),budgets)
        if sum(b['warm_elapsed_ns'] for b in budgets)>1200e9:raise TimeoutError('Total warm work exceeded 1200 seconds')
if __name__=='__main__':main()
