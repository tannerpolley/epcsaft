"""Bounded native checks of source pairing, equivalent sites and unstable liquid refusal."""
import copy
import gzip
import json
import math
import sys
import subprocess
import time
import uuid
from datetime import datetime, timezone
from pathlib import Path
from epcsaft import Mixture, Pure, Parameters, _core, provenance, trace
from epcsaft.equilibrium import (EquilibriumStart, EquilibriumWorkLimits, Free, MoleFractions,
    Phase, Pinned, Problem, solve_equilibrium, solve_equilibrium_starts)
from reproduce import CRITERIA, LEAF, QA, RESULTS, ROOT, dump, serial, sha, targets

def association_probes(document, run):
    retained = json.loads(gzip.decompress((RESULTS/'diagnosis.json.gz').read_bytes()))
    identity = dict(run['execution_identity'], run_id=str(uuid.uuid4()), parent_run_id=run['execution_identity']['run_id'],
        started_utc=datetime.now(timezone.utc).isoformat(), checkout_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        checkout_status=subprocess.check_output(['git','status','--porcelain'],cwd=ROOT,text=True), interpreter=sys.executable,
        producer_sha256=sha(Path(__file__)), parent_result_sha256=sha(RESULTS/'native-results.json.gz'),
        engine_provenance=serial(provenance()), extension_path=_core.__file__, extension_sha256=sha(Path(_core.__file__)),
        input_sha256={str(p.relative_to(LEAF)):sha(p) for p in sorted((LEAF/'inputs').rglob('*')) if p.is_file()})
    records, traces, started = [], [], time.perf_counter()
    coefficients = {c['component_id']:{v['family']:v['value']['magnitude'] for v in c['coefficients']} for c in document['components']}
    selves = {e['endpoint_a']['component_id']:e for e in document['topology']['edges'] if e['source']['kind']=='explicit'}
    for row in targets():
        pure = row['figure']=='1'
        names = [row['mdea']] + ([] if pure else [row['water']])
        for offset in ([0] if pure else [-1,0,1]):
            T = (426.8 if row['mdea']=='MDEA_2B' else 431.0) if pure else row['source_y']+offset*row['y_half_width']
            corrected, volumes = copy.deepcopy(document), []
            for edge in corrected['topology']['edges']:
                a,b = [edge[k]['component_id'] for k in ['endpoint_a','endpoint_b']]
                si,sj = [coefficients[n]['segment_diameter'] for n in [a,b]]
                di,dj = [coefficients[n]['segment_diameter']*(1-0.12*math.exp(-3*coefficients[n]['dispersion_energy_over_k']/T)) for n in [a,b]]
                kappa = math.sqrt(selves[a]['volume']['value']['magnitude']*selves[b]['volume']['value']['magnitude'])
                size, ratio = (math.sqrt(si*sj)/((si+sj)/2))**3, ((di+dj)/(si+sj))**3
                edge['source'] = dict(kind='explicit',provenance=edge['source']['provenance'])
                edge['volume'] = dict(identity=f"probe/{a}/{edge['endpoint_a']['site_id']}/{b}/{edge['endpoint_b']['site_id']}/volume",value=dict(magnitude=kappa*size*ratio,unit='dimensionless'))
                edge['energy_over_k'] = dict(identity=edge['volume']['identity']+'/energy',value=dict(magnitude=(selves[a]['energy_over_k']['value']['magnitude']+selves[b]['energy_over_k']['value']['magnitude'])/2,unit='kelvin'))
                if a in names and b in names:
                    volumes.append(dict(edge=edge,geometric_self_volume=kappa,sigma_size_factor=size,d3_over_sigma3=ratio))
            for convention, candidate in [('raw-sigma3',document),('explicit-d3-isothermal',corrected)]:
                model = (Pure if pure else Mixture)(Parameters.from_mapping(candidate,components=names))
                problem = Problem(T=T,P=Free(row['pressure_Pa']),feed=MoleFractions({names[0]:1} if pure else {names[0]:row['source_x'],names[1]:1-row['source_x']}),phases=[Phase('L',kind='liquid',amount=Pinned(0) if row['branch']=='dew' else Free()),Phase('V',kind='vapor',amount=Pinned(0) if row['branch']=='bubble' else Free())],row_tolerance=CRITERIA['row_tolerance'])
                with trace(diagnostics=True) as capture:
                    result = solve_equilibrium(model,problem)
                    from epcsaft.regression import observation, predict
                    density = serial(predict(model,observation(model.parameters,'phase_density',problem,phase='L',observed=row['source_y'],form='difference',mass_basis=True,scale=row['y_half_width']))) if pure and result.success else None
                records.append(dict(row_id=row['row_id'],offset=offset,convention=convention,run_id=identity['run_id'],temperature_K=T,source=row,problem=serial(problem),volumes=volumes if convention!='raw-sigma3' else [],result=serial(result),liquid_mass_density=density,trace_index=len(traces)))
                traces.append(serial(capture.trace))
    retained['association_probes'] = dict(execution_identity=identity,expression='Eqs. (4)-(7): d_i=sigma_i*(1-0.12*exp(-3*epsilon_i/T)); s_ij=(sigma_i+sigma_j)/2; kappa_explicit=sqrt(kappa_ii*kappa_jj)*(sqrt(sigma_i*sigma_j)/s_ij)^3*(d_ij/s_ij)^3; isothermal only',parameter_sha256=sha(LEAF/'inputs/parameters/parameters.json'),records=records,wall_seconds=time.perf_counter()-started)
    dump(RESULTS/'diagnosis.json.gz',retained)
    dump(QA/'association-probe-traces.json.gz',traces)
    print('isothermal association declarations',len(records),'wall seconds',time.perf_counter()-started)

def main():
    assert Path(__import__('sys').prefix).resolve() == ROOT/'.venv'
    run = json.loads(gzip.decompress((RESULTS/'native-results.json.gz').read_bytes()))
    identity = dict(run['execution_identity'], run_id=str(uuid.uuid4()), parent_run_id=run['execution_identity']['run_id'],
        started_utc=datetime.now(timezone.utc).isoformat(), producer_sha256=sha(Path(__file__)),
        parent_result_sha256=sha(RESULTS/'native-results.json.gz'), engine_provenance=serial(provenance()),
        extension_sha256=sha(Path(_core.__file__)), checkout_status=subprocess.check_output(['git','status','--porcelain'],cwd=ROOT,text=True))
    document = json.loads((LEAF/'inputs/parameters/parameters.json').read_text())
    if '--association-only' in sys.argv:
        association_probes(document, run)
        return
    swapped = copy.deepcopy(document)
    for site in swapped['topology']['sites']:
        site['site_role'] = 'acceptor' if site['site_role']=='donor' else 'donor'
        partner = next(s for s in document['topology']['sites'] if s['component_id']==site['component_id'] and s['site_id']!=site['site_id'])
        site['multiplicity'] = partner['multiplicity']
    records, traces = [], []
    started = time.perf_counter()
    for row in targets():
        pure = row['figure']=='1'
        if not pure and row['pressure_Pa']!=40000 and not row['row_id'].startswith('4d/bubble'):
            continue
        names = [row['mdea']] + ([] if pure else [row['water']])
        model = (Pure if pure else Mixture)(Parameters.from_mapping(document,components=names))
        feed = MoleFractions({names[0]:1} if pure else {names[0]:row['source_x'],names[1]:1-row['source_x']})
        phases = [Phase('L',kind='liquid',amount=Pinned(0) if row['branch']=='dew' else Free()),
                  Phase('V',kind='vapor',amount=Pinned(0) if row['branch']=='bubble' else Free())]
        point = next(p for p in run['points'] if p['row_id']==row['row_id'] and p['offset']==0 and p['requested_tolerance']==CRITERIA['row_tolerance'])
        if point['success']:
            for name, candidate, T, P in [('inverse-pressure',model,point['temperature'],Free(row['pressure_Pa'])),
                ('equivalent-site-swap',(Pure if pure else Mixture)(Parameters.from_mapping(swapped,components=names)),Free(point['temperature']),row['pressure_Pa'])]:
                result = solve_equilibrium(candidate,Problem(T=T,P=P,feed=feed,phases=phases,row_tolerance=CRITERIA['refinement_tolerance']))
                records.append(dict(row_id=row['row_id'],probe=name,source_temperature_K=point['temperature'],source_pressure_Pa=row['pressure_Pa'],result=serial(result)))
        else:
            for offset in [-1,0,1]:
                x = row['source_x']+offset*row['x_half_width']
                problem = Problem(T=Free(row['source_y']),P=row['pressure_Pa'],feed=MoleFractions({names[0]:x,names[1]:1-x}),phases=phases,row_tolerance=CRITERIA['row_tolerance'])
                seeds = [EquilibriumStart(temperature=row['source_y']+k*row['y_half_width']) for k in [-1,0,1]]
                with trace(diagnostics=True) as capture:
                    result = solve_equilibrium_starts(model,problem,seeds,EquilibriumWorkLimits(attempts=len(seeds)))
                records.append(dict(row_id=row['row_id'],offset=offset,probe='native-source-interval-starts',seeds=serial(seeds),result=serial(result)))
                traces.append(serial(capture.trace))
    dump(RESULTS/'diagnosis.json.gz',dict(execution_identity=identity,records=records,wall_seconds=time.perf_counter()-started))
    dump(QA/'diagnostic-traces.json.gz',traces)
    print('diagnostic declarations',len(records),'wall seconds',time.perf_counter()-started)

if __name__ == '__main__':
    main()
