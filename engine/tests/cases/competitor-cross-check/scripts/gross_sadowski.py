import json
from pathlib import Path
from epcsaft import Parameters, Pure, Mixture, regression as r
from epcsaft.equilibrium import MoleFractions, Free, Phase, Problem, bubble_point
ROOT = Path(__file__).resolve().parents[5]


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
