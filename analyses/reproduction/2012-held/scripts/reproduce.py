"""Native Held 2012 comparisons; Python declares states and formats retained values."""
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
import time
import uuid
import epcsaft
from epcsaft import Mixture, Parameters, _core, equilibrium as eq, regression as reg

LEAF = Path(__file__).resolve().parents[1]
ROOT = LEAF.parents[2]
SOURCE = LEAF/'inputs/source'
PACKET = LEAF/'inputs/parameters/held-2012-alcohols/1/parameters'
ION = {'Li':'lithium-cation','Na':'sodium-cation','K':'potassium-cation','NH4':'ammonium-cation','Cl':'chloride-anion','Br':'bromide-anion','I':'iodide-anion'}

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def read(p):
    return list(csv.DictReader(p.open()))

def fields(v, names):
    return {n:getattr(v,n) for n in names.split()}

def serial(v):
    if v is None or isinstance(v,(str,bool,int)):return v
    if isinstance(v,float):return v if math.isfinite(v) else str(v)
    if isinstance(v,dict):return {str(k):serial(x) for k,x in v.items()}
    if isinstance(v,(list,tuple)):return [serial(x) for x in v]
    if hasattr(v,'name') and not hasattr(v,'category'):return v.name
    return serial({k:getattr(v,k) for k in dir(v) if not k.startswith('_') and not callable(getattr(v,k))})

def target_rows():
    targets=[]
    for table,name in [('1','held-2012-table-1-density.csv'),('pure','held-2012-pure-ethanol-density.csv'),('2','held-2012-table-2-activities.csv'),('3','held-2012-table-3-activities.csv')]:
        for i,r in enumerate(read(SOURCE/name)):
            for quantity,column in ([('density','density_kg_m3')] if table in ['1','pure'] else [('osmotic','osmotic_coefficient'),('ionic','mean_ionic_activity_coefficient_molality_scale')]):
                targets.append(dict(case_id=f't{table}-{i}-{quantity}',group='Table '+table,quantity=quantity,solvent=r['solvent'],salt=r.get('salt','NaI'),molality=r.get('molality_mol_kg',0),source_value=float(r[column]),x_interval=(.5*10**(-len(r['molality_mol_kg'].split('.')[-1])) if 'molality_mol_kg' in r else 0),y_interval=.005 if quantity=='density' else .0005,evidence_label='physical comparison',locator=r['locator']))
    assert len(targets)==110
    for r in read(LEAF/'inputs/figures/targets.csv'):targets.append(dict(r,group='Figure '+r['figure'],locator='Held 2012 Figure '+r['figure']))
    targets.append(dict(case_id='shared-730',group='Shared #730',quantity='density',solvent='ethanol',salt='NaBr',molality=.1422,source_value=795.96,x_interval=(.5*10**(-len(r['molality_mol_kg'].split('.')[-1])) if 'molality_mol_kg' in r else 0),y_interval=.005,evidence_label='physical comparison',locator='Table 1; Figiel /4 parameters, historical masses'))
    return targets

def declaration(r):
    if r['case_id']=='shared-730':
        p=ROOT/'data/packets/figiel-2025-reference-electrolytes/4/parameters/parameters.json';d=json.loads(p.read_text())
        for c in d['components']:
            if c['component_id'] in ['sodium-cation','bromide-anion']:c['fixed']['molar_mass']={'value':dict(magnitude=.02298 if c['component_id']=='sodium-cation' else .0799,unit='kilogram / mole')}
    else:
        p=LEAF/'inputs/declarations'/f"{r['series']}.json" if 'series' in r else PACKET/f"{r['solvent']}.json";d=json.loads(p.read_text())
    solvent=['ethanol'] if r['case_id']=='shared-730' else [c['component_id'] for c in d['components'] if c['fixed']['charge_number']['value']['magnitude']==0]
    a=next(a for a in ['NH4','Na','Li','K'] if r['salt'].startswith(a));b=r['salt'][len(a):]
    ids=(*solvent,ION[a],ION[b]);ratio=float(r.get('ratio',0));sf=([ratio,1-ratio] if len(solvent)==2 else [1])+[0,0]
    return p,d,ids,sf

def declarations():
    return [Parameters.from_mapping(d,components=ids) for r in target_rows() if r.get('figure')!='7' for p,d,ids,sf in [declaration(r)]]

def pair_size_inputs():
    ions={(r['ion'],r['solvent']):float(r['segment_diameter_angstrom']) for r in read(SOURCE/'held-2012-table-6-ions.csv')}
    solvents={r['component']:float(r['value']) for r in read(SOURCE/'held-2012-table-5-solvents-and-pairs.csv') if r['quantity']=='segment_diameter'}
    for key,s in json.loads((LEAF/'inputs/figures/series.json').read_text()).items():
        if not 0<float(s['ratio'])<1 or s['figure']=='7':continue
        source=json.loads((LEAF/f'inputs/declarations/{key}.json').read_text());d=Parameters.from_mapping(source)._document
        print(key,len(d.component_ids),*d.component_ids)
        print(d.base_choice,d.association_choice,d.electrolyte_choice,d.ionic_dispersion_choice,d.permittivity_choice,d.polar_choice)
        print(len(d.site_component),*sum(([i,m] for i,m in zip(d.site_component,d.site_multiplicity)),[]))
        print(len(d.entries))
        for e in d.entries:print(json.dumps(e.identity),e.field,len(e.index),*e.index,repr(e.magnitude),json.dumps(e.unit))
        print(14)
        for j,name in enumerate(['water',s['solvent']]):
            sigma=(2.7927+10.11*math.exp(-.01775*298.15)-1.417*math.exp(-.01146*298.15)) if name=='water' else solvents[name]
            for i,c in enumerate(source['components'][2:],2):print(j,i,repr(.5*(sigma+ions[c['aliases'][0],name])))

def calculate(r, identity, pressure=100000., options=None, screening=False, molality=None):
    start=time.perf_counter();record=dict(case_id=r['case_id'],execution_run_id=identity['run_id'],pressure_Pa=pressure,temperature_K=298.15,quantity=r['quantity'],value=None)
    if r.get('figure')=='7':return dict(record,cause='undetermined: source permittivities unavailable')
    p,d,ids,sf=declaration(r)
    if screening:
        for c in d['components']:
            for v in c['coefficients']:
                if v['family']=='debye_huckel_diameter':v['value']['magnitude']*=.88
    record.update(parameter_sha256=sha(p),component_ids=ids,solvent_fractions=sf,molality=float(r['molality']) if molality is None else molality,screening_diagnostic=screening)
    try:
        params=Parameters.from_mapping(d,components=ids);model=Mixture(params,solver_options=options)
        x=sf if record['molality']==0 else _core.formula_composition(sf,params._native,(len(ids)-2,len(ids)-1),(1,1),record['molality']);record.update(composition=x,molar_masses_kg_mol=params.molar_masses)
        if r['quantity']=='pressure':
            problem=eq.bubble_point(liquid=eq.Phase('L',kind='liquid'),T=298.15,feed=eq.Molalities({ids[-2]:record['molality'],ids[-1]:record['molality']},solvent=tuple(ids[:-2]),solvent_mass_kg=1))
            result=eq.solve_equilibrium(model,problem)
            record.update(fields(result,'success message convergence requested_tolerance_met iterations evaluations linear_solver solver_options'))
            if result.success:record['value']=result.pressure/100000
        elif r['quantity']=='density':
            obs=reg.density(params,298.15,pressure,dict(zip(ids,x)),phase='liquid',observed=float(r['source_value']),form='difference',mass_basis=True,scale=float(r['y_interval']))
            result=reg.predict(model,obs);record.update(fields(result,'value status message convergence branch stable_root_count requested_tolerance_met solver_options'))
            state=model.state(298.15,P=pressure,x=x,phase='liquid');record['state']=fields(state,'molar_density pressure pressure_density_derivative density_convergence density_branch stable_root_count')
        else:
            ref=model.reference(298.15,pressure,solvent_fractions=sf,ions=ids[-2:],stoichiometry=(1,1));record['reference']=fields(ref,'molality convergence_error log_fugacity');record['solver_options']=ref._native.solver_options
            if r['quantity']=='gamma':
                result=ref.solvent_activity(x=x);record.update(fields(result,'solvent_activity status cause osmotic_cause'));record['value']=result.solvent_activity/x[0] if result.solvent_activity is not None else None
                if record['value'] is not None:assert abs(x[0]*record['value']-result.solvent_activity)<=1e-12
            else:
                result=ref.activity(record['molality']);record['value']=result.osmotic_coefficient if r['quantity']=='osmotic' else result.mean_ionic_molality_scale;record['solver_options']=result.solver_options
            state=model.state(298.15,P=pressure,x=x,phase='liquid');record['state']=fields(state,'molar_density pressure pressure_density_derivative density_convergence density_branch stable_root_count')
        if not isinstance(record['value'],(float,int)) or not math.isfinite(record['value']):record.update(value=None,cause='native result unavailable: '+str(record['value']))
    except Exception as e:
        record.update(value=None,cause=type(e).__name__+': '+str(e))
        if 'model' in locals() and 'x' in locals():
            record['density_probe']=serial(model._model.density_diagnostics(298.15,pressure,x,[1.,100.,1000.,5000.,15000.,25000.],1e-6,1e-10,200,10.))
    record['wall_seconds']=time.perf_counter()-start
    return serial(record)

def comparison_rows(targets,values,intervals):
    out=[]
    for r,v in zip(targets,values):
        value=v['value'];error=value-float(r['source_value']) if value is not None else None
        bounds=intervals.get(r['case_id'],[]);yr=float(r['y_interval']);monotone=bool(bounds) and (bounds[0]<=bounds[2]<=bounds[1] or bounds[0]>=bounds[2]>=bounds[1])
        overlap=monotone and min(bounds)<=float(r['source_value'])+yr and max(bounds)>=float(r['source_value'])-yr
        outcome='undetermined' if value is None or (bounds and not monotone) else ('reproduced' if (abs(error)<=yr or overlap) else 'not reproduced')
        out.append(dict(case_id=r['case_id'],group=r['group'],series=r.get('series',''),quantity=r['quantity'],solvent=r['solvent'],salt=r['salt'],molality=r['molality'],source_value=r['source_value'],computed_value=value,signed_error=error,absolute_tolerance=yr,x_interval=r['x_interval'],interval_values=json.dumps(bounds),interval_bound=('sampled monotone' if monotone else 'undetermined: interior extremum or absent interval'),central_outcome=('reproduced' if abs(error)<=yr else 'not reproduced') if error is not None else 'undetermined',outcome=outcome,cause=v.get('cause','') if value is None else ('physical discrepancy; cause undetermined' if r['evidence_label']=='physical comparison' and outcome!='reproduced' else ('model-curve discrepancy; cause undetermined' if outcome=='not reproduced' else 'interior interval extremum unresolved' if outcome=='undetermined' else '')),evidence_label=r['evidence_label'],execution_run_id=v['execution_run_id'],locator=r['locator']))
    return out

def compare(targets,run):
    out=comparison_rows(targets,run['quantities'],run['intervals'])
    primary={r['case_id']:r for r in out};probe=run['screening_probe']
    alternate=comparison_rows(probe['targets'],probe['quantities'],probe['intervals'])
    for r in alternate:
        p=primary[r['case_id']];r.update(probe='screening 0.88 sigma; packing unchanged; not a source convention',primary_value=p['computed_value'],primary_signed_error=p['signed_error'],primary_outcome=p['outcome'],primary_execution_run_id=p['execution_run_id'])
    for name,records in [('comparisons.csv',out),('screening-comparisons.csv',alternate)]:
        with (LEAF/'results'/name).open('w') as f:
            w=csv.DictWriter(f,fieldnames=records[0],lineterminator='\n');w.writeheader();w.writerows(records)
    summaries=[]
    for table in ['7','8']:
        name=f'held-2012-table-{table}-model-deviations.csv'
        for i,s in enumerate(read(SOURCE/name)):
            quantity={'mean_ionic_activity_coefficient':'ionic','osmotic_coefficient':'osmotic','density':'density'}.get(s['quantity'])
            cohort=[r for r in out if r['group'] in ['Table 2','Table 3'] and r['solvent']==s.get('solvent') and r['salt']==s.get('salt') and r['quantity']==quantity]
            complete=table=='7' and bool(cohort) and len(cohort)==int(s['number_of_points']) and abs(max(float(r['molality']) for r in cohort)-float(s['maximum_molality_mol_kg']))<=.005 and not(s['solvent']=='ethanol' and s['salt']=='LiI') and all(r['computed_value'] is not None for r in cohort)
            for metric,column in [('ARD_percent','average_relative_deviation_percent'),('AAD','average_absolute_deviation')]:
                if column not in s or not s[column]:continue
                value=(sum(abs(r['signed_error'])/(float(r['source_value']) if metric=='ARD_percent' else 1) for r in cohort)/len(cohort)*(100 if metric=='ARD_percent' else 1)) if complete else None
                tol=.005 if table=='7' else .05
                summaries.append(dict(table=table,solvent=s.get('solvent',s.get('system')),salt=s.get('salt',''),quantity=s['quantity'],metric=metric,source_value=s[column],computed_value=value,absolute_tolerance=tol,outcome=('reproduced' if abs(value-float(s[column]))<=tol else 'not reproduced') if value is not None else 'undetermined',cause='' if complete else ('ethanol/LiI cohort contradiction' if s.get('solvent')=='ethanol' and s.get('salt')=='LiI' else 'exact source cohort unavailable'),execution_run_id=run['execution_identity']['run_id']))
    with (LEAF/'results/summary-comparisons.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=summaries[0],lineterminator='\n');w.writeheader();w.writerows(summaries)

def main():
    if '--pair-size-inputs' in sys.argv:
        pair_size_inputs();return
    targets=target_rows();path=LEAF/'results/native-results.json.gz'
    if '--compare-only' in sys.argv:
        run=json.loads(gzip.decompress(path.read_bytes()))
    else:
        freeze=subprocess.check_output(['git','log','-1','--format=%H','--',str(LEAF/'inputs/figures')],cwd=ROOT,text=True).strip()
        file_list=''.join(f'{sha(p)}  {p.relative_to(LEAF)}\n' for p in sorted((LEAF/'inputs/figures').iterdir()) if p.is_file())
        start=time.perf_counter();identity=dict(run_id=str(uuid.uuid4()),started_utc=datetime.now(timezone.utc).isoformat(),source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),engine_provenance=serial(epcsaft.provenance()),core_path=str(_core.__file__),core_sha256=sha(Path(_core.__file__)),runner_sha256=sha(Path(__file__)),input_sha256={str(p.relative_to(LEAF)):sha(p) for p in sorted((LEAF/'inputs').rglob('*')) if p.is_file()},freeze_commit=freeze,freeze_file_list_sha256=hashlib.sha256(file_list.encode()).hexdigest())
        previous=json.loads(gzip.decompress(path.read_bytes())) if '--refresh-figures' in sys.argv else None
        old={r['case_id']:r for r in previous['quantities']} if previous else {}
        run=dict(execution_identity=identity,execution_identities=(previous['execution_identities'] if previous else [])+[identity],quantities=[],intervals={},interval_diagnostics={},diagnostics=[])
        if previous:run['shared_exact_equality']=previous['shared_exact_equality']
        for r in targets:
            v=old[r['case_id']] if previous and 'series' not in r and not r['case_id'].startswith('tpure') else calculate(r,identity);run['quantities'].append(v)
            if float(r['x_interval'])>0 and v['value'] is not None:
                ends=[calculate(r,identity,molality=float(r['molality'])+sign*float(r['x_interval'])) for sign in [-1,1]]
                run['interval_diagnostics'][r['case_id']]=ends
                if all(e['value'] is not None for e in ends):run['intervals'][r['case_id']]=[e['value'] for e in ends]+[v['value']]
        shared=next(v for v in run['quantities'] if v['case_id']=='shared-730');run['shared_exact_equality']=shared['value']==788.0580935539466;assert shared['value'] is not None and abs(shared['value']-788.0580935539466)<=1e-10
        probe=dict(role='screening 0.88 sigma; packing unchanged; not a source convention',targets=[],quantities=[],intervals={},interval_diagnostics={})
        for r in targets:
            if r.get('figure') not in ['2','3','5','6'] or r['quantity'] not in ['osmotic','ionic','gamma']:continue
            probe['targets'].append(r);v=calculate(r,identity,screening=True);probe['quantities'].append(v)
            ends=[calculate(r,identity,screening=True,molality=float(r['molality'])+sign*float(r['x_interval'])) for sign in [-1,1]]
            probe['interval_diagnostics'][r['case_id']]=ends
            if v['value'] is not None and all(e['value'] is not None for e in ends):probe['intervals'][r['case_id']]=[e['value'] for e in ends]+[v['value']]
        run['screening_probe']=probe
        for key in (['tpure-0-density'] if previous else ['tpure-0-density','shared-730']):
            r=next(r for r in targets if r['case_id']==key);run['diagnostics'].append(dict(role='101325 Pa sensitivity',result=calculate(r,identity,pressure=101325.)))
        for q in ['density','osmotic','ionic','gamma','pressure']:
            valid=[(r,v) for r,v in zip(targets,run['quantities']) if r['quantity']==q and v['value'] is not None]
            if valid:
                r,v=max(valid,key=lambda pair:abs(pair[1]['value']-float(pair[0]['source_value']))/float(pair[0]['y_interval']))
                options={o['name']:o['value']/10 for o in v.get('solver_options',[]) if any(t in o['name'] for t in ['DENSITY_ROOT_RESIDUAL','ASSOCIATION_RESIDUAL','EQUILIBRIUM_ROW_RESIDUAL','REFERENCE_CONVERGENCE_TOLERANCE'])}
                run['diagnostics'].append(dict(role='tenfold native tolerance refinement',result=calculate(r,identity,options=options)))
        if previous:run['diagnostics'] += [d for d in previous['diagnostics'] if d['role']=='101325 Pa sensitivity' and d['result']['case_id']=='shared-730']
        run['wall_seconds']=time.perf_counter()-start;path.parent.mkdir(exist_ok=True);path.write_bytes(gzip.compress(json.dumps(run,separators=(',',':'),allow_nan=False).encode(),mtime=0))
    compare(targets,run)
    print(json.dumps(dict(run_id=run['execution_identity']['run_id'],wall_seconds=run['wall_seconds'],targets=len(targets))))

if __name__=='__main__':main()
