"""#209 native guarded minimax; Python declares inputs and reads native results."""
import copy
import hashlib
import json
import time
import math

from epcsaft import Parameters, _core, provenance, regression as r
import stage1_fit as F
from stage1_diagnosis import state_checks

COMPONENTS = ('water', 'methanol', 'ethanol', *F.IONS)
F.STEM = 'stage1-minimax'


def assembled(value=0.0, document=None):
    doc=copy.deepcopy(F.BASE if document is None else document)
    doc['pairs'].append({'component_id_a':'ethanol','component_id_b':'methanol','coefficients':[
        {'identity':'pair/ethanol/methanol/k_ij','family':'k_ij','value':{'magnitude':value,'unit':'dimensionless'}}]})
    return Parameters.from_mapping(doc,components=COMPONENTS)


def inert_comparisons(points, observations, coordinates, holds):
    records=[]
    for point,u in enumerate(points):
        evaluations=[_core.evaluate_observations(assembled(value)._native,coordinates,u,holds,observations) for value in (0.0,0.5)]
        a,b=evaluations
        assert all(s==r.ObservationStatus.Available for e in evaluations for s in e.statuses)
        for i in range(len(observations)):
            pairs=[('value',a.predictions[i],b.predictions[i])]+[
                (f'optimizer-column-{j}',a.optimizer_gradient[i*7+j],b.optimizer_gradient[i*7+j]) for j in range(7)]
            for quantity,x,y in pairs:
                tolerance=1e-12*abs(x) if x else 1e-12
                records.append(dict(point=point,row=i,quantity=quantity,zero_pair=x,half_pair=y,
                    absolute_change=abs(x-y),tolerance=tolerance,passed=abs(x-y)<=tolerance))
        assert all(v['passed'] for v in records)
    F.write('inert-comparisons',records)
    return dict(rows=len(observations),points=len(points),comparisons=len(records),
        maximum_absolute_change=max(v['absolute_change'] for v in records),all_pass=True,
        completion='in-memory methanol-ethanol k_ij=0; never fitted or retained in packet',alternate_value=0.5)


def active_columns(p, coordinates, holds, observations, optimizer, indices):
    evaluation=_core.evaluate_observations(p._native,coordinates,optimizer,holds,[observations[i] for i in indices])
    assert all(s==r.ObservationStatus.Available for s in evaluation.statuses)
    rows=[dict(row=i,column=c,coordinate=F.NAMES[c],prediction=evaluation.predictions[j],value_unit='1',
               exact_prediction_derivative=evaluation.optimizer_gradient[j*7+c])
          for j,i in enumerate(indices) for c in range(7)]
    assert all(v['exact_prediction_derivative']==0 for v in rows if v['column'] in (1,2,4,5))
    F.write('active-columns',rows)


def stored_checks(path, result, start):
    source=assembled()
    stored=assembled(document=json.loads(path.read_text()))
    _,_,obs,meta,groups=F.declarations(source,'all')
    checks=r.evaluate_checks(stored,obs,groups,source)
    values=_core.observations_at(stored._native,[],[],[],obs,[0.0]*len(obs))
    # Freeze the returned physical point, including held diameters, without solving the holds again.
    fitted=Parameters(result.chart.parameters,source.component_ids,source.site_ids)
    fitted_values=_core.observations_at(fitted._native,[],[],[],obs,[0.0]*len(obs))
    F.write('fitted-minus-cold',[dict(row=i,fitted=native,cold=cold.observed,difference=native-cold.observed)
        for i,(cold,native) in enumerate(zip(values,result.predictions))])
    comparisons,rows=[],[]
    for i,(cold,point) in enumerate(zip(values,fitted_values)):
        native=point.observed
        tolerance=1e-12*abs(native) if native else 1e-12
        comparisons.append(dict(row=i,stored=cold.observed,native=native,error=abs(cold.observed-native),
            tolerance=tolerance,passed=abs(cold.observed-native)<=tolerance))
    F.write('stored-comparisons',comparisons)
    assert all(v['passed'] for v in comparisons)
    for check in checks:
        for i,error in zip(check.rows,check.errors):
            rows.append(meta[i]|dict(start=start,prediction=values[i].observed,error=error,unit=check.unit,status=check.status))
    records=F.checks_record(checks,start)
    states=state_checks(source,stored,obs,meta,start)
    F.write('stored-groups',records)
    F.write('groups',F.checks_record(r.evaluate_checks(source,obs,groups,source),-1)+records)
    F.write('states',states)
    F.write('rows',rows)
    return all(v['passed'] for v in states) and all(v['status']=='pass' for v in records if not v['report_only'])


def finite(value):
    if isinstance(value,dict):
        return {k:finite(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):
        return [finite(v) for v in value]
    return None if isinstance(value,float) and not math.isfinite(value) else value


def preserve_affine(output):
    """Evaluate the selected held point and fixed packet against their own retained operands."""
    from pathlib import Path
    from types import SimpleNamespace
    import csv
    reference = F.OUT
    output = Path(output); output.mkdir(parents=True,exist_ok=True)
    summary = json.loads((reference/'stage1-minimax-fit.json').read_text())
    selected = summary['selected']; point = summary['starts'][selected]
    packet = F.REPO/'data/packets/figiel-2025-chloride-activity-refit/1/parameters/parameters.json'
    source, stored = assembled(), assembled(document=json.loads(packet.read_text()))
    cs,hs,obs,meta,groups = F.declarations(source,'all')
    chart = r.evaluate_coordinates(source,cs,point['optimizer'],hs)
    fixed_cs = F.declarations(stored,'all')[0]
    fixed_cs += [r.coordinate(stored,[('born_diameter',ion,1)],origin=0,scale=1) for ion in F.ALKALI]
    inverse = r.CoordinateSet(); inverse.coordinates = fixed_cs
    fixed = r.evaluate_coordinates(stored,fixed_cs,inverse.start(stored._native),[])
    records = []
    def compare(point,row,quantity,unit,new,old):
        available = new not in (None,'') and old not in (None,'')
        if available:
            new,old = float(new),float(old)
            tolerance = 1e-12*abs(old) if old else 1e-12
            error = abs(new-old); passed = math.isfinite(error) and error<=tolerance
        else:
            error,tolerance = None,None; passed = new in (None,'') and old in (None,'')
        records.append(dict(point=point,row=row,quantity=quantity,unit=unit,old=old,new=new,error=error,tolerance=tolerance,passed=passed))
    for i,value in enumerate(chart.physical):
        compare('held',i,'physical','native parameter unit',value,point['physical'][i])
        compare('held-to-fixed',i,'physical','native parameter unit',value,fixed.physical[i])
    for i,hold in enumerate(chart.holds):
        assert abs(hold.energy_error)<=1e-6 and hold.correction<=1e-10
    evaluations = [_core.evaluate_observations(assembled(v)._native,cs,point['optimizer'],hs,obs) for v in (0,0.5)]
    old = [v for v in F.read(reference/'stage1-minimax-inert-comparisons.csv') if int(v['point'])==2]
    assert len(old)==7584 and len(obs)==948
    for operand in old:
        row = int(operand['row']); quantity = operand['quantity']
        values = [e.predictions[row] if quantity=='value' else e.optimizer_gradient[row*7+int(quantity.rsplit('-',1)[1])] for e in evaluations]
        for value,name in zip(values,('zero_pair','half_pair')):
            compare(name,row,quantity,'native observable unit per optimizer unit' if quantity!='value' else 'native observable unit',value,operand[name])
        compare('inert-pair',row,quantity,'native',values[1],values[0])
    assert all(s==r.ObservationStatus.Available for e in evaluations for s in e.statuses)
    for operand in F.read(reference/'stage1-minimax-active-columns.csv'):
        row,column = int(operand['row']),int(operand['column'])
        compare('active-columns',row,f'optimizer-column-{column}','native observable unit per optimizer unit',
            evaluations[0].optimizer_gradient[row*7+column],operand['exact_prediction_derivative'])
    predictions = _core.observations_at(stored._native,[],[],[],obs,[0.0]*len(obs))
    for row,operand in enumerate(F.read(reference/'stage1-minimax-stored-comparisons.csv')):
        compare('fixed',row,'prediction','native observable unit',predictions[row].observed,operand['stored'])
    F.OUT = output
    accepted = stored_checks(packet,SimpleNamespace(predictions=[v.observed for v in predictions],chart=fixed),selected)
    for name,columns in (('rows',('prediction','error','observed')),('stored-groups',('maximum','mean','maximum_threshold','mean_threshold','baseline_maximum','baseline_mean')),
                         ('states',('value','reference','relative_error','pressure_density_derivative','charge'))):
        fresh,old = F.read(output/f'{F.STEM}-{name}.csv'),F.read(reference/f'{F.STEM}-{name}.csv')
        assert len(fresh)==len(old)
        for row,(new,operand) in enumerate(zip(fresh,old)):
            for key in operand.keys()-set(columns): assert new[key]==operand[key],(name,row,key,new[key],operand[key])
            for key in columns: compare('fixed-'+name,row,key,operand.get('unit','native'),new[key],operand[key])
    with (output/'comparisons.csv').open('w') as stream:
        writer = csv.DictWriter(stream,fieldnames=list(records[0])); writer.writeheader(); writer.writerows(records)
    failures = [v for v in records if not v['passed']]
    metadata = dict(selected=selected,conditions=summary['conditions'],rows=len(obs),columns=7,physical_rows=len(chart.physical),
        derivative_basis='prediction dy/dv for every observation form; unit row weights',
        reference_commits=dict(inert_comparisons='ce0f5411',active_columns='062ccd92'),
        holds=[dict(energy_error_J_mol=h.energy_error,correction_angstrom=h.correction,energy_tolerance_J_mol=1e-6,
                    diameter_tolerance_angstrom=1e-10) for h in chart.holds],
        comparisons=len(records),failures=len(failures),accepted=accepted,stage2_started=False,
        maximum_absolute_error=max(v['error'] for v in records if v['error'] is not None),
        maximum_tolerance_fraction=max(v['error']/v['tolerance'] for v in records if v['tolerance']),
        source_sha256=hashlib.sha256(F.SOURCE.read_bytes()).hexdigest(),packet_sha256=hashlib.sha256(packet.read_bytes()).hexdigest(),
        wheel_sha256=hashlib.sha256(next((F.REPO/'build/environment-wheel').glob('*.whl')).read_bytes()).hexdigest())
    (output/'summary.json').write_text(json.dumps(metadata,indent=2)+'\n')
    print(json.dumps(metadata),flush=True)
    assert accepted and not failures,failures[:5]


def main():
    p=assembled()
    cs,hs,obs,meta,groups=F.declarations(p,'all')
    assert all(not (o.composition[1]>0 and o.composition[2]>0) and
        not (o.reference_composition and o.reference_composition[1]>0 and o.reference_composition[2]>0) for o in obs)
    ceres=json.loads((F.OUT/'stage1-native-fit.json').read_text())
    minima=[]
    for candidate in sorted((v for v in ceres['starts'] if v['usable']),key=lambda v:v['cost']):
        u=candidate['optimizer']
        if all(max(abs(x-y)/(6 if j==6 else 2) for j,(x,y) in enumerate(zip(u,v)))>.01 for v in minima):
            minima.append(u)
    source=[float(F.read(F.HERE/'inputs/stage1-starts.csv')[0][name]) for name in F.NAMES]
    evidence=inert_comparisons([source,*minima],obs,cs,hs)
    checks=copy.copy(groups)
    for g in checks:
        if g.name.startswith(('T ','E ','I ')):
            g.report_only=True
    spec=r.MinimaxRequest()
    spec.targets=[i for i,v in enumerate(meta) if v['set']=='T']
    spec.group_scales=[1000.0 if any(obs[i].kind in (r.ObservationKind.SolvationGibbsEnergy,r.ObservationKind.TransferGibbsEnergy) for i in g.rows)
                       else .02 if g.metric==r.CheckMetric.AbsoluteLog10Ratio else 1.0 for g in checks]
    spec.q,spec.tolerance,spec.activity_distance,spec.margin=.02,1e-10,1e-6,.1
    controls=r.FitControls(maximum_iterations=600,maximum_elapsed_time_seconds=900)
    print('Native minimax:',len(minima),'distinct Ceres minima;',len(obs),'rows',flush=True)
    began=time.perf_counter()
    results=r.fit(p,cs,obs,holds=hs,checks=checks,baseline=p,controls=controls,starts=minima,
        shared_iterations=600*len(minima),shared_seconds=900*len(minima),minimax=spec)
    summaries,constraints=[],[]
    for start,result in enumerate(results):
        m=result.minimax
        print(start,result.status,result.message,m.solver_status if m else None,m.t if m else None,flush=True)
        summary=dict(start=start,status=result.status,message=result.message,usable=result.usable,optimizer=result.optimizer,
            physical=result.physical,iterations=result.iterations,covariance_status=result.covariance.status,
            holds=[{key:getattr(h,key) for key in ('energy_error','correction','uncertainty','derivative')} for h in result.chart.holds])
        if m:
            summary['minimax']={key:getattr(m,key) for key in ('availability','application_status','solver_status','feasible','stationary','t','target_maximum','gap',
                'primal','stationarity','dual_sign','complementarity','point','lower','upper','physical_z_lower','physical_z_upper','active_lower','active_upper',
                'lower_distances','upper_distances','hold_multipliers','hold_causes','full_stationarity')}
            summary['maximum_root_error']=max(m.root_errors,default=0)
            for i,row in enumerate(m.rows):
                constraints.append(dict(start=start,constraint=i,kind=row.kind,row=row.row,group=checks[row.group].name if row.group<len(checks) else 'T',
                    sign=row.sign,threshold=row.threshold,scale=row.scale,auxiliary=row.auxiliary,members=';'.join(map(str,row.members)),
                    scaled_value=m.values[i],physical_value=m.physical_values[i],scaled_multiplier=m.multipliers[i],physical_multiplier=m.physical_multipliers[i],root_error=m.root_errors[i]))
        summaries.append(summary)
    if constraints:
        F.write('constraints',constraints)
    record=dict(method='native Ipopt guarded minimax, exact first derivatives, D5 eliminated holds',starts=summaries,distinct_ceres_minima=minima,
        inert_completion=evidence,seconds=time.perf_counter()-began,engine_provenance=list(provenance()),
        coordinate_order=F.NAMES,conditions=dict(temperature_K=F.T,pressure_Pa=F.P,branch='liquid',molality_basis='mol formula units/kg water',target_interval=(.001,4.5)),
        input_source_sha256=hashlib.sha256(F.SOURCE.read_bytes()).hexdigest(),
        controls=dict(maximum_iterations=600,maximum_seconds=900,q=.02,tolerance=1e-10,activity_distance=1e-6,margin=.1),analysis_gross_ceiling=650,engine_runtime_additions=0,stage2_started=False)
    path=F.OUT/'stage1-minimax-fit.json'
    path.write_text(json.dumps(finite(record),indent=2,allow_nan=False)+'\n')
    if not any(v.usable for v in results):
        raise RuntimeError('No native minimax point available; failures retained')
    start,result=min(((i,v) for i,v in enumerate(results) if v.usable),key=lambda pair:pair[1].minimax.target_maximum)
    candidate=F.OUT/'stage1-minimax-parameters.json'
    candidate.write_text(json.dumps(F.stored_document(result),indent=2)+'\n')
    active_columns(p,cs,hs,obs,result.optimizer,[row.row for row,dual in zip(result.minimax.rows,result.minimax.physical_multipliers) if row.kind in ('target','maximum') and dual>1e-4])
    record['selected']=start
    record['inert_completion']=inert_comparisons([source,*minima,result.optimizer],obs,cs,hs)
    record['accepted']=stored_checks(candidate,result,start)
    path.write_text(json.dumps(finite(record),indent=2,allow_nan=False)+'\n')


def declarations():
    return F.declarations(assembled(), 'all', declare_only=True)


if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--preserve-affine',action='store_true')
    parser.add_argument('--output')
    args = parser.parse_args()
    if args.preserve_affine:
        if not args.output: parser.error('--preserve-affine requires --output')
        preserve_affine(args.output)
    else: main()
