"""Installed-wheel numerical falsifiers; finite differences are verification only."""
import csv
import hashlib
import json
import math
import subprocess
import time
from pathlib import Path

from epcsaft import ActiveParameter, Mixture, ParameterFamily, Parameters, _core, regression as r
from declarations import COMPONENTS, CATIONS, ROOT, PACKET, row_groups, coordinates, holds, parameters, row, target_rows

OUT = Path(__file__).resolve().parents[1] / 'results'
OUT.mkdir(exist_ok=True)
records = []
failures = []
timings = []
ladder = []


def retain(case, quantity, exact, reference, tolerance, decisive=True, **extra):
    error = abs(exact-reference)
    passed = math.isfinite(error) and error <= tolerance if decisive else ''
    records.append(dict(case=case, quantity=quantity, exact=exact, reference=reference, error=error, tolerance=tolerance, passed=passed, **extra))
    if decisive and not passed:
        failures.append(records[-1])


def native_start(p, cs, observations=(), thermochemistry=None):
    chart = r.CoordinateSet(); chart.coordinates = cs
    return chart.start(p._native, list(observations), thermochemistry)


def prediction_columns(p, cs, u, hs, observations):
    """Native prediction columns: Difference residuals with unit scales and weights."""
    rows = _core.observations_at(p._native, cs, u, hs, list(observations), [0.0]*len(observations))
    for ob in rows: ob.form, ob.scale = r.ResidualForm.Difference, 1.0
    return _core.evaluate_observations(p._native, cs, u, hs, rows).optimizer_gradient


def manufactured(p, cs, u, hs, observations, errors=None):
    return _core.observations_at(p._native, cs, u, hs, observations, [0.0]*len(observations) if errors is None else errors)


def f1(p):
    specs = [('k_ij', ((ion, partner),)) for ion in CATIONS for partner in ('water','chloride-anion')]
    specs += [('ion_suppression_coefficient', ())] + [('born_diameter',(ion,)) for ion in CATIONS]
    specs += [(family,('chloride-anion',)) for family in ('segment_diameter','debye_huckel_diameter','packing_diameter')]
    count=0
    for group, (p, observations) in enumerate(row_groups()):
        cs=[]; mapping=[]
        for c,(family, targets) in enumerate(specs):
            if family=='ion_suppression_coefficient': targets=tuple(name for name in p.component_ids if name.endswith(('cation','anion')))
            needed = targets[0] if family=='k_ij' else targets
            if not all(name in p.component_ids for name in needed): continue
            cs.append(r.coordinate(p,[(family, target, 1) for target in targets],origin=0,scale=1)); mapping.append(c)
        u=native_start(p,cs)
        chart=r.CoordinateSet(); chart.coordinates=cs
        model=Mixture(p,active_parameters=chart.active())
        started=time.perf_counter()
        raw=[r.differentiate(model,ob).value for ob in observations]; tie=chart.tie()
        mapped=prediction_columns(p,cs,u,[],observations)
        exact_seconds=time.perf_counter()-started; difference_seconds=0.0
        started=time.perf_counter(); manufactured(p,cs,u,[],observations)
        difference_seconds += time.perf_counter()-started
        assert all(len(g)==len(tie) for g in raw), ('unavailable raw column',group)
        seen = set()
        for ob in observations:
            if ob.kind == r.ObservationKind.BulkRelativePermittivity: continue
            for x in (ob.composition,ob.reference_composition):
                if not x or (tuple(x),tuple(ob.ion_indices)) in seen: continue
                seen.add((tuple(x),tuple(ob.ion_indices)))
                reference = model.reference(298.15,1e5,solvent_fractions=x,
                    ions=tuple(p.component_ids[i] for i in ob.ion_indices),stoichiometry=(1,1))._native
                ladder.append({'group':group,'ions':':'.join(str(i) for i in ob.ion_indices),'composition':':'.join(str(v) for v in x),
                    'column':'value','molality':reference.molality,'scaled_terminal_error':reference.convergence_error})
                for k,columns in enumerate(reference.parameter_derivative):
                    assert all(v is not None for v in columns), ('unavailable reference column',group,k)
                    ladder.append({'group':group,'ions':':'.join(str(i) for i in ob.ion_indices),'composition':':'.join(str(v) for v in x),
                        'column':k,'molality':reference.molality,'scaled_terminal_error':None})
        print('F1 group' ,group,'rows',len(observations),'physical columns',len(tie),flush=True)
        for local,c in enumerate(mapping):
            exact=[mapped[i*len(cs)+local] for i in range(len(observations))]
            for h in (1e-3,1e-4,1e-5):
                step=h*max(abs(u[local]),1.0)
                plus,minus=list(u),list(u);plus[local]+=step;minus[local]-=step
                started=time.perf_counter()
                yp=manufactured(p,cs,plus,[],observations);ym=manufactured(p,cs,minus,[],observations)
                if h==1e-4 or (h==1e-3 and any(ob.kind==r.ObservationKind.OsmoticCoefficient and ob.molality<.01 for ob in observations)):
                    difference_seconds += time.perf_counter()-started
                for i,ob in enumerate(observations):
                    cd=(yp[i].observed-ym[i].observed)/(2*step)
                    delta=298.15*8.31446261815324*1e-7 if 'GibbsEnergy' in ob.kind.name else 1e-7
                    tol=2e-5*max(abs(v) for v in exact)+delta
                    prescribed=h==(1e-3 if ob.kind==r.ObservationKind.OsmoticCoefficient and ob.molality<0.01 else 1e-4)
                    label=f'row-{count+i}/column-{c}'
                    if prescribed:retain('F1',label,exact[i],cd,tol)
                    if not prescribed: records.append(dict(case='F1-step',quantity=f'{label}/h-{h}',exact=exact[i],reference=cd,error=abs(exact[i]-cd),tolerance=tol,passed=''))
        timings.append({'group':group,'rows':len(observations),'columns':len(cs),'exact_value_columns_seconds':exact_seconds,'central_difference_value_columns_seconds':difference_seconds})
        count+=len(observations)
    return count


def export_reference(p, result, cs, hs, filename):
    doc = p._document
    values = {}
    fields = {'Kij': 'k_ij', 'IonSuppressionCoefficient': 'ion_suppression_coefficient'}
    active = r.CoordinateSet(); active.coordinates = cs
    for target, value in zip(active.active(), result.physical):
        for index in target.indices:
            indices = (index//len(p.component_ids), index % len(p.component_ids)) if target.family == ParameterFamily.Kij else (index,)
            values[(fields[target.family.name], indices)] = value
    for hold, value in zip(hs, result.physical[len(active.active()):]):
        values[('born_diameter', (hold.reference.species[0],))] = value
    suppression = [value for target,value in zip(active.active(),result.physical) if target.family == ParameterFamily.IonSuppressionCoefficient]
    assert suppression and all(value == suppression[0] for value in suppression)
    quote = json.dumps
    with filename.open('w') as out:
        out.write(f'{len(doc.component_ids)} {len(doc.site_component)} {len(doc.entries)}\n')
        out.write(' '.join(map(quote, doc.component_ids))+'\n')
        for comp, mult in zip(doc.site_component, doc.site_multiplicity): out.write(f'{comp} {mult}\n')
        out.write(' '.join(quote(getattr(doc, n)) for n in ('base_choice','association_choice','electrolyte_choice','ionic_dispersion_choice','permittivity_choice','polar_choice'))+'\n')
        for e in doc.entries:
            magnitude = values.get((e.field, tuple(e.index)), values.get((e.field, tuple(reversed(e.index))), e.magnitude))
            if e.field == 'ion_suppression_coefficient': magnitude = suppression[0]
            out.write(f'{quote(e.identity)} {quote(e.field)} {len(e.index)} '+' '.join(map(str, e.index))+f' {magnitude:.17g} {quote(e.unit)}\n')
        out.write(f'{len(hs)}\n')
        for h in hs: out.write(f'{h.reference.species[0]} {h.reference.ion_indices[1]} 0 298.15 100000 {h.reference.observed}\n')


def fitting(p):
    cs, hs = coordinates(p), holds(p)
    source = native_start(p, cs)
    truth = [v+(0.05 if abs(v+0.05) <= 0.95 else -0.05) if c % 2 == 0 else v+(0.1 if abs(v+0.1) <= 0.95 else -0.1) for c, v in enumerate(source[:6])] + [source[6]-1.0]
    chart = r.evaluate_coordinates(p, cs, truth, hs)
    observations = manufactured(p, cs, truth, hs, target_rows(p))
    controls = r.FitControls(maximum_iterations=100, maximum_elapsed_time_seconds=120,
                             function_tolerance=1e-15, gradient_tolerance=1e-14, parameter_tolerance=1e-14)
    reports = {}
    output_rows = []
    for name, errors in (('F2', [0.0]*len(observations)), ('F3', [-0.005 if ob.kind == r.ObservationKind.MeanIonicActivityCoefficient and ob.ion_indices[0] == p.component_ids.index('sodium-cation') else 0.0 for ob in observations])):
        data = observations if name == 'F2' else manufactured(p, cs, truth, hs, observations, errors)
        started = time.perf_counter()
        from checks import group
        checks = [group('T',list(range(len(data))),r.CheckMetric.AbsoluteLog10Ratio,.02)]
        weights = [1/(0.02*math.log(10))**2]*len(data)
        if name == 'F3':
            checks.append(group('E',list(range(len(data),len(data)+3)),r.CheckMetric.AbsoluteDifference,10))
            data += [h.reference for h in hs]; weights += [0]*3
        result = r.fit(p, cs, data, holds=hs, controls=controls, singular_vectors=True, weights=weights, checks=checks)
        print(name, result.status, result.final_cost, result.message, flush=True)
        reports[name] = {'status':result.status, 'optimizer':result.optimizer, 'physical':result.physical,
            'cost':result.final_cost, 'covariance_status':result.covariance.status,'active_bounds':result.covariance.active_bounds,
            'projected_gradient_norm':result.covariance.projected_gradient_norm,'iterations':result.iterations,
            'value_evaluations':result.residual_evaluations,'column_evaluations':result.jacobian_evaluations, 'seconds':time.perf_counter()-started, 'singular_values':result.covariance.singular_values,
            'holds':[{'energy_error':h.energy_error,'correction':h.correction,'uncertainty':h.uncertainty} for h in result.chart.holds],
            'failures':[{'row':f.observation,'cause':f.message,'optimizer':f.optimizer} for f in result.failures]}
        for i,(ob,prediction,residual) in enumerate(zip(data,result.predictions,result.residuals)):
            output_rows.append({'case':name,'row':i,'kind':ob.kind.name,'ions':':'.join(str(j) for j in ob.ion_indices),
                'molality':ob.molality,'observed':ob.observed,'prediction':prediction,'declared_residual':residual})
        assert result.status == 'converged', reports[name]
        if name == 'F2':
            for c, (found, expected) in enumerate(zip(result.optimizer, truth)): retain(name, f'coordinate-{c}',found,expected,1e-5 if c==6 else 1e-6)
            for i, (found, expected) in enumerate(zip(result.physical[-3:],chart.physical[-3:])): retain(name,f'diameter-{i}',found,expected,1e-7)
            for i, residual in enumerate(result.residuals): retain(name,f'decade-{i}',residual/math.log(10),0,1e-10)
        for i,h in enumerate(result.chart.holds):
            retain(name,f'inner-energy-{i}',h.energy_error,0,1e-6); retain(name,f'inner-step-{i}',h.correction,0,1e-10)
        if name == 'F3':
            stored = r.fitted_parameters(p,result)
            separate = r.evaluate_checks(stored,data,checks)
            for c,(returned,fresh) in enumerate(zip(result.checks,separate)):
                for quantity in ('maximum','mean'):
                    a,b=getattr(returned,quantity),getattr(fresh,quantity)
                    retain(name,f'check-{c}-{quantity}',a,b,1e-12*max(abs(b),1))
                for i,(a,b) in enumerate(zip(returned.errors,fresh.errors)):
                    retain(name,f'check-{c}-row-{i}',a,b,1e-12*max(abs(b),1))
            job = OUT/'independent-reference-inputs.txt'; export_reference(p,result,cs,hs,job)
            driver = subprocess.check_output([str(ROOT/'tools/native-driver'),'regression_checks',
                str(Path(__file__).with_name('electrolyte_reference.cpp'))],text=True).strip()
            output = subprocess.check_output([str(driver),str(job)],text=True,timeout=60)
            (OUT/'independent-roots.csv').write_text(output)
            for values in csv.DictReader(output.splitlines()):
                retain(name,f"root-{values['ion']}",float(values['stored_diameter']),float(values['independent_diameter']),1e-6)
                retain(name,f"energy-{values['ion']}",float(values['energy_error']),0,1e-2)
    with (OUT/'fitted-observations.csv').open('w') as out:
        writer=csv.DictWriter(out,fieldnames=tuple(output_rows[0]));writer.writeheader();writer.writerows(output_rows)
    reports['manufactured_truth'] = truth
    return reports


def main():
    started=time.perf_counter(); p=parameters()
    from limits import identities, debye
    identities(retain); limits = debye(retain)
    assert not failures, failures[:3]
    rows=f1(p)
    with (OUT/'falsifiers.csv').open('w') as out:
        writer=csv.DictWriter(out,fieldnames=('case','quantity','exact','reference','error','tolerance','passed'));writer.writeheader();writer.writerows(records)
    if failures:
        (OUT/'failed-columns.json').write_text(json.dumps(failures,indent=2)+'\n')
        raise AssertionError(failures[:3])
    from checks import metrics, maps
    checks=metrics(retain); charts=maps(retain)
    assert not failures, failures[:3]
    build = Path((ROOT/'build/environment-wheel/meson-build-dir').read_text().strip())
    output = subprocess.check_output([str(build/'regression_checks'),'hydration-failures'],text=True,timeout=60)
    native_failures=[]
    for line in output.splitlines():
        if line.startswith('EPCSAFT_COMPARE_JSON '):
            value=json.loads(line.split(' ',1)[1]); native_failures.append({'case':value['comparison_id'],'inputs':value['inputs']})
            retain('F6',value['comparison_id'],value['model'],value['reference'],0)
    assert len(native_failures)==4
    fits=fitting(p)
    wheels=list((ROOT/'build/environment-wheel').glob('*.whl')); assert len(wheels)==1
    summary={'wheel_sha256':hashlib.sha256(wheels[0].read_bytes()).hexdigest(),'packet_sha256':hashlib.sha256(PACKET.read_bytes()).hexdigest(),
        'conditions':{'T_K':298.15,'P_Pa':100000,'branch':'liquid'},'rows':rows,'debye_limit':limits,'checks':checks,'native_failure_cases':native_failures,'charts':charts,'timings':timings,'fits':fits,'failures':failures,'seconds':time.perf_counter()-started}
    with (OUT/'reference-ladders.csv').open('w') as out:
        writer=csv.DictWriter(out,fieldnames=tuple(ladder[0]));writer.writeheader();writer.writerows(ladder)
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    with (OUT/'falsifiers.csv').open('w') as out:
        writer=csv.DictWriter(out,fieldnames=('case','quantity','exact','reference','error','tolerance','passed'));writer.writeheader();writer.writerows(records)
    assert not failures, failures[:3]


if __name__=='__main__':main()
