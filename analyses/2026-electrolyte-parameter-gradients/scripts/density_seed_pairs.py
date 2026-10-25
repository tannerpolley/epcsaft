"""Seven paired complete calls (paired), one cold/seeded endpoint comparison (endpoints), or the endpoint comparison against the
propagated per-coordinate allowance and the B4 minimax objective (allowance), using the existing F2/B4 declaration owners.
Arguments: mode out cold_env seeded_env [stops JSON, the same stopping controls for both arms]."""
import csv, json, math, os, subprocess, sys, time
from statistics import median as med
from pathlib import Path
from epcsaft import Parameters, _core, regression as r
from declarations import parameters, coordinates, holds, target_rows, PACKET
MOLAR_MASS = {c['component_id']: c['fixed']['molar_mass']['value'] for c in json.loads(PACKET.read_text())['components'] if 'molar_mass' in c.get('fixed', {})}
def floor(p, ob):
    """Absolute floor: 1e-12 on ln a_w, carried to Phi = -ln a_w/(M_w sum(nu) m); 1e-12 for other values."""
    if ob.kind != _core.ObservationKind.OsmoticCoefficient: return 1e-12
    (solvent,) = [p.component_ids[i] for i, x in enumerate(ob.composition) if x > 0]
    mass = MOLAR_MASS[solvent]; assert mass['unit'] == 'gram / mole'
    return 1e-12/(mass['magnitude']/1000*sum(ob.stoichiometry)*ob.molality)
from run import native_start, manufactured
import minimax as b
OUT = Path(sys.argv[2]); OUT.mkdir(parents=True, exist_ok=True)
def linearization(p, data, weights, x):
    """Per-coordinate bound sum_i |[J_w^+]_{k,i}| c_i at the endpoint, J_w^+ = V S^-1 U^T from the fit's returned SVD of the
    weighted optimizer Jacobian; c_i = sqrt(w_i) |dr_i/dy_i| a_i carries the observable allowance a_i to the weighted residual."""
    if x.minimax: return dict(optimizer=list(x.optimizer), objective=x.minimax.t)
    cov, width = x.covariance, len(x.optimizer); rows, s = cov.objective_rows, cov.singular_values
    U, V = cov.left_vectors, cov.right_vectors; assert cov.rank == width and len(U) == len(rows)*width, (cov.status, cov.rank)
    c, relative, osmotic = [], [], []
    for i in rows:
        y, ob = x.predictions[i], data[i]
        slope = math.sqrt(weights[i])*abs(1/ob.scale if ob.form == r.ResidualForm.Difference else 1/y)
        c.append(slope*max(floor(p, ob), 1e-10*abs(y))); relative.append(slope*1e-10*abs(y))
        osmotic.append(ob.kind == _core.ObservationKind.OsmoticCoefficient)
    # Layout check: U S V^T reproduces the returned weighted columns.
    layout = max(abs(sum(U[a*width+j]*s[j]*V[k*width+j] for j in range(width)) - x.optimizer_jacobian[i*width+k])
                 for a, i in enumerate(rows) for k in range(width))/s[0]
    pseudo = [[abs(sum(V[k*width+j]/s[j]*U[a*width+j] for j in range(width))) for a in range(len(rows))] for k in range(width)]
    bound = [sum(m*e for m, e in zip(row, c)) for row in pseudo]
    # Diagnostics: the share carried by osmotic rows, and the bound with every absolute floor removed.
    share = [sum(m*e for m, e, o in zip(row, c, osmotic) if o)/b for row, b in zip(pseudo, bound)]
    return dict(optimizer=list(x.optimizer), bound=bound, osmotic_share=share, relative_only=[sum(m*e for m, e in zip(row, relative)) for row in pseudo],
                condition=s[0]/s[-1], layout=layout,
                kinds=sorted({str(data[i].kind) for i in rows}))
def call_fit(workload, boundary, replay, stops=None, linearize=False):
    """stops: optional {'F2': [function, gradient, parameter] Ceres tolerances, 'B4': Ipopt tolerance} for both arms."""
    p = parameters() if workload == 'F2' else parameters(b.COMPONENTS_B)
    cs, hs = coordinates(p), holds(p); source = native_start(p, cs)
    controls = r.FitControls(maximum_iterations=1 if boundary == 'one_step' else (100 if workload == 'F2' else 300), maximum_elapsed_time_seconds=120)
    if workload == 'F2':
        truth = [v+(.05 if abs(v+.05)<=.95 else -.05) if k%2==0 else v+(.1 if abs(v+.1)<=.95 else -.1) for k,v in enumerate(source[:6])] + [source[6]-1]
        data = manufactured(p, cs, truth, hs, target_rows(p)); checks = [b.group('T', list(range(len(data))), r.CheckMetric.AbsoluteLog10Ratio, .02)]
        baseline = spec = None; starts = None
        controls.function_tolerance, controls.gradient_tolerance, controls.parameter_tolerance = (stops or {}).get('F2', (1e-15, 1e-14, 1e-14))
    else:
        data, _, _ = b.source_and_targets(p, cs, hs, source); indices = b.indices_for(data); checks = b.b4_checks(indices)
        chart = r.evaluate_coordinates(p, cs, source, hs); baseline = Parameters(chart.parameters, p.component_ids, p.site_ids)
        spec = b.minimax_spec(indices[0], [1,1000,.02]); shifted = list(source)
        spec.tolerance = (stops or {}).get('B4', spec.tolerance)
        shifted[2] = math.copysign(max(0, abs(shifted[2])-.01), shifted[2]) if shifted[2] else .01; shifted[6] += .1
        starts = [source] if boundary == 'one_step' else [source,shifted]
    weights = [1/(.02*math.log(10))**2 if workload == 'F2' or k < 18 else 0 for k in range(len(data))]
    result = r.fit(p, cs, data, weights=weights, controls=controls, holds=hs, checks=checks, baseline=baseline, minimax=spec, starts=starts, shared_iterations=600, shared_seconds=240, singular_vectors=True)
    results = result if isinstance(result,list) else [result]
    work = [(x.status,len(x.predictions),len(x.physical_jacobian),len(x.optimizer_jacobian), x.trial_failures, len(x.minimax.values) if x.minimax else 0,
        x.residual_evaluations, x.jacobian_evaluations, list(x.optimizer)) for x in results]
    assert all(x.status == 'converged' for x in results), work
    if linearize: return [linearization(p, data, weights, x) for x in results]
    if replay:
        b.OUT=OUT
        try:
            for x in results:
                cold = _core.evaluate_observations(x.parameters, cs, x.optimizer, hs, data)
                for field, reference in [('physical',cold.chart.physical),('predictions',cold.predictions),('physical_jacobian',cold.gradient)]:
                    for k,(a,z) in enumerate(zip(getattr(x,field),reference,strict=True)):
                        if field != 'physical_jacobian' or weights[k//len(x.physical)] > 0:
                            atol = floor(p, data[k]) if field == 'predictions' else 1e-12
                            b.record(f'{workload}-{boundary}',field+str(k),a,z,max(atol,1e-10*abs(z)),atol=atol,
                                molality=data[k].molality if field == 'predictions' else '')
        finally: b.write_outputs(dict(scope='same-point independent cold endpoint replay'))
    return work
if sys.argv[1] == 'worker':
    reply=os.fdopen(os.dup(1),'w'); os.dup2(2,1)
    for line in sys.stdin:
        reply.write(json.dumps(call_fit(*json.loads(line)))+'\n'); reply.flush()
else:
    workers={label:subprocess.Popen([str(Path(sys.argv[i])/'bin/python'),__file__,'worker',str(OUT/label)],stdin=subprocess.PIPE,stdout=subprocess.PIPE,text=True) for label,i in [('A',3),('B',4)]}
    stops = json.loads(sys.argv[5]) if len(sys.argv) > 5 else None
    def call(label,*args):
        process=workers[label]; process.stdin.write(json.dumps(args+(stops,))+'\n'); process.stdin.flush()
        return json.loads(process.stdout.readline())
    output = OUT/'pairs.csv'
    if sys.argv[1] == 'allowance':
        def endpoint_run(workload, stops_):
            got = {}
            for label in 'AB':
                process=workers[label]; process.stdin.write(json.dumps([workload,'full',False,stops_,True])+'\n'); process.stdin.flush()
                got[label] = json.loads(process.stdout.readline())
            return got
        retained = {'script-controls': (None, OUT.parent/'round3/endpoints.csv'),
                    'work-matched': ({'F2': [1e-15, 1e-14, 1e-8]}, OUT.parent/'round4/work-matched/endpoints.csv')}
        rows, objective = [], []
        try:
            for name, (stops_, path) in retained.items():
                kept = {(x['workload'], int(x['start']), x['quantity']): x for x in csv.DictReader(path.open())}
                for workload in ('F2', 'B4'):
                    got = endpoint_run(workload, stops_)
                    for start, (z, a) in enumerate(zip(got['A'], got['B'], strict=True)):
                        for k, (u, v) in enumerate(zip(a['optimizer'], z['optimizer'], strict=True)):
                            old = kept[(workload, start, f'optimizer{k}')]
                            assert (float(old['seeded']), float(old['cold'])) == (u, v), (name, workload, start, k)
                            if workload == 'F2':
                                rows.append(dict(controls=name, start=start, coordinate=k, seeded=u, cold=v, difference=u-v, allowance=z['bound'][k],
                                    fraction=abs(u-v)/z['bound'][k], seeded_endpoint_allowance=a['bound'][k], passed=abs(u-v) <= z['bound'][k],
                                    osmotic_share=z['osmotic_share'][k], relative_only_allowance=z['relative_only'][k],
                                    condition=z['condition'], layout_check=max(z['layout'], a['layout']), kinds=' '.join(z['kinds'])))
                            else:
                                rows.append(dict(controls=name, workload='B4', start=start, coordinate=k, seeded=u, cold=v, difference=u-v))
                        if workload == 'B4':
                            t, tc = a['objective'], z['objective']; allowance = max(1e-12, 1e-10*abs(tc))
                            objective.append(dict(controls=name, start=start, seeded=t, cold=tc, difference=t-tc, relative=(t-tc)/abs(tc),
                                allowance=allowance, fraction=abs(t-tc)/allowance, passed=abs(t-tc) <= allowance))
        finally:
            for process in workers.values(): process.stdin.close(); process.wait()
        f2 = [x for x in rows if 'workload' not in x]; b4 = [x for x in rows if 'workload' in x]
        for name, table in (('coordinate-allowance.csv', f2), ('b4-coordinates.csv', b4), ('b4-objective.csv', objective)):
            with (OUT/name).open('w') as stream:
                writer=csv.DictWriter(stream,fieldnames=tuple(table[0]),lineterminator='\n'); writer.writeheader(); writer.writerows(table)
        sys.exit(0)
    try:
        records, summaries, endpoints = [], [], []
        for workload in ('F2','B4'):
            first = {label: call(label,workload,'full',False) for label in 'AB'}
            for start, (a, z) in enumerate(zip(first['B'], first['A'], strict=True)):
                names = ('status','predictions','physical_columns','optimizer_columns','trial_failures','minimax_values','residual_evaluations','jacobian_evaluations')
                endpoints += [dict(workload=workload,start=start,quantity=n,seeded=u,cold=v,allowance=0,passed=u==v) for n,u,v in zip(names,a,z)]
                endpoints += [dict(workload=workload,start=start,quantity=f'optimizer{k}',seeded=u,cold=v,allowance=max(1e-12,1e-10*abs(v)),
                    passed=abs(u-v)<=max(1e-12,1e-10*abs(v))) for k,(u,v) in enumerate(zip(a[-1],z[-1],strict=True))]
            with (OUT/'endpoints.csv').open('w') as stream:
                writer=csv.DictWriter(stream,fieldnames=tuple(endpoints[0]),lineterminator='\n'); writer.writeheader(); writer.writerows(endpoints)
            if sys.argv[1] == 'endpoints': continue
            selected = []; times, work = {}, {}
            for pair in range(7):
                order = 'AB' if (pair+(workload=='B4'))%2==0 else 'BA'
                for label in order:
                    started=time.perf_counter(); work[label]=call(label,workload,'full',False); times[label]=time.perf_counter()-started
                assert all(x[:6] == y[:6] for x, y in zip(work['A'], work['B'], strict=True)), (workload, pair, work)
                n = {label: sum(x[6]+x[7] for x in work[label]) for label in 'AB'}; each = {label: times[label]/n[label] for label in 'AB'}
                selected.append(dict(workload=workload,boundary='full',pair=pair,order=order,**times,saving=(times['A']-times['B'])/times['A'],
                    evaluations_A=n['A'],evaluations_B=n['B'],identical_work=[x[6:8] for x in work['A']]==[x[6:8] for x in work['B']],
                    per_evaluation_A=each['A'],per_evaluation_B=each['B'],per_evaluation_saving=(each['A']-each['B'])/each['A']))
                records.append(selected[-1]); print(workload,pair+1,order,times,flush=True)
                with output.open('w') as stream:
                    writer=csv.DictWriter(stream,fieldnames=tuple(records[0])); writer.writeheader(); writer.writerows(records)
            savings=[x['saving'] for x in selected]; center=med(savings); spread=med([abs(x-center) for x in savings])
            per=[x['per_evaluation_saving'] for x in selected]
            summaries.append(dict(workload=workload,boundary='full',median_A=med([x['A'] for x in selected]),median_B=med([x['B'] for x in selected]),
                median_paired=center,range=[min(savings),max(savings)],mad=spread,passed=center>=.05 if workload=='F2' else center>=-spread,
                evaluations_A=selected[0]['evaluations_A'],evaluations_B=selected[0]['evaluations_B'],identical_work=all(x['identical_work'] for x in selected),
                median_per_evaluation_A=med([x['per_evaluation_A'] for x in selected]),median_per_evaluation_B=med([x['per_evaluation_B'] for x in selected]),
                median_per_evaluation_saving=med(per),per_evaluation_range=[min(per),max(per)],per_evaluation_mad=med([abs(x-med(per)) for x in per]),stops=stops))
            output.with_suffix('.json').write_text(json.dumps(summaries,indent=2))

    finally:
        for process in workers.values(): process.stdin.close(); process.wait()
