"""Neutral MDEA selection; stop on the design's failed-state or fit gates."""
import csv
import json
from pathlib import Path
import epcsaft
from epcsaft import regression, _core
from epcsaft.equilibrium import Free, MoleFractions, Phase, Pinned, Problem

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'analyses/reproduction/2010-nasrifar-mdea/inputs/source/figure-targets.csv'

def main():
    rows = [r for r in csv.DictReader(SOURCE.open()) if r['series']=='measured' and r['panel'] in ('6a','6b')]
    report = {'rows': rows, 'fit_coordinates': 1, 'electrolyte_parameters_changed': False, 'reaction_constants_changed': False, 'fits': []}
    candidates = [('inherited-4C-2B', HERE/'inputs/run-document.json', 'methyldiethanolamine'),
                  ('nasrifar-3B-2B', ROOT/'data/packets/nasrifar-2010-mdea/1/parameters/parameters.json', 'MDEA')]
    for name, path, amine in candidates:
        parameters = epcsaft.Parameters.from_mapping(json.loads(path.read_text()), components=(amine,'water'))
        pair = next(p for p in json.loads(path.read_text())['pairs'] if {p['component_id_a'],p['component_id_b']}=={amine,'water'})
        law = pair['coefficients'][0]
        if 'temperature_law' in law:
            term = next(v for v in law['temperature_law']['terms'] if v['power']==0)['coefficient']
            family, target, origin = 'temperature_law_coefficient', term['identity'], term['value']['magnitude']
        else:
            family, target, origin = 'k_ij', (amine,'water'), law['value']['magnitude']
        for panel in ('6a','6b'):
            training = [r for r in rows if r['panel']==panel]
            held = [r for r in rows if r['panel']!=panel]
            items = []
            for r in training+held:
                problem = Problem(T=Free(float(r['y'])),P=float(r['pressure_Pa']),
                    feed=MoleFractions({amine:float(r['x']),'water':1-float(r['x'])}),
                    phases=[Phase('L',kind='liquid'),Phase('V',kind='vapor',amount=Pinned(0))])
                items.append(regression.observation(parameters,'temperature',problem,observed=float(r['y']),form='difference',scale=1.0,mass_basis=False))
            for offset in (0,-.1,.1):
                coordinate = regression.coordinate(parameters,[(family,target,1)],origin=origin,scale=1,bounds=(origin-.5,origin+.5))
                physical_start = list(regression.evaluate_coordinates(parameters,[coordinate],[offset],()).physical)
                assert abs(physical_start[0]-(origin+offset)) < 1e-12
                chart = regression.CoordinateSet(); chart.coordinates = [coordinate]
                omitted = regression.evaluate_coordinates(parameters,[coordinate],chart.start(parameters._native),()).physical[0]
                assert abs((omitted-physical_start[0])+offset) < 1e-12
                assert offset==0 or abs(omitted-physical_start[0]) > 1e-8
                result = regression.fit(parameters,[coordinate],items[:len(training)],weights=[1.]*len(training),holds=(),
                    starts=[[offset]], controls=regression.FitControls(maximum_iterations=100,maximum_elapsed_time_seconds=600))[0]
                evaluation = _core.evaluate_observations(parameters._native,[coordinate],list(result.optimizer),[],items[len(training):])
                baseline_coordinate = regression.coordinate(parameters,[(family,target,1)],origin=origin,scale=1,bounds=(origin-.5,origin+.5))
                baseline = _core.evaluate_observations(parameters._native,[baseline_coordinate],[0.],[],items[len(training):])
                residuals = [v-float(r["y"]) for v,r in zip(evaluation.predictions,held)]
                baseline_residuals = [v-float(r["y"]) for v,r in zip(baseline.predictions,held)]
                record = dict(candidate=name,training_panel=panel,start_offset=offset,physical_start=physical_start[0],omitted_start_defect=omitted-physical_start[0],status=result.status,message=result.message,
                    usable=result.usable,physical=list(result.physical),residuals=list(result.residuals),iterations=result.iterations,
                    held_cost=.5*sum(v*v for v in residuals),held_predictions=list(evaluation.predictions),held_residuals=list(residuals),
                    held_statuses=[str(v) for v in evaluation.statuses],baseline_held_cost=.5*sum(v*v for v in baseline_residuals),
                    conditional_covariance=list(result.covariance.physical),variance_factor=result.covariance.variance_factor,
                    projected_gradient_norm=result.covariance.projected_gradient_norm,
                    inherited_origin=origin,final_cost=result.final_cost,rank=result.covariance.rank,covariance_status=result.covariance.status,
                    active_bounds=list(result.covariance.active_bounds),singular_values=list(result.covariance.singular_values),
                    failures=[dict(trial=f.trial,observation=f.observation,status=str(f.status),message=f.message) for f in result.failures])
                report['fits'].append(record)
                (HERE/'results/e4.json').write_text(json.dumps(report,indent=2)+'\n')
                if not result.usable or result.covariance.active_bounds or result.covariance.status!='available':
                    report['stop']='Neutral MDEA fit gate failed; no further fit or extension under design section 7.'
                    (HERE/'results/e4.json').write_text(json.dumps(report,indent=2)+'\n')
                    return 1
    report['start_controls'] = []
    for name, _, _ in candidates:
        for panel in ('6a','6b'):
            fits = [f for f in report['fits'] if f['candidate']==name and f['training_panel']==panel]
            separation = min(abs(a['physical_start']-b['physical_start']) for i,a in enumerate(fits) for b in fits[i+1:])
            spread = max(f['physical'][0] for f in fits)-min(f['physical'][0] for f in fits)
            report['start_controls'].append(dict(candidate=name,panel=panel,minimum_separation=separation,
                omitted_start_defects=[f['omitted_start_defect'] for f in fits],limit=1e-8,final_parameter_spread=spread))
            assert separation > 1e-8
    scores = {}
    for panel in ('6a','6b'):
        original = next(f for f in report['fits'] if f['candidate']=='inherited-4C-2B' and f['training_panel']==panel and f['start_offset']==0)
        alternative = next(f for f in report['fits'] if f['candidate']=='nasrifar-3B-2B' and f['training_panel']==panel and f['start_offset']==0)
        scores[panel] = dict(current_held_cost=original['held_cost'],candidate_held_cost=alternative['held_cost'],
            improvement=1-alternative['held_cost']/original['held_cost'],passes_20_percent=alternative['held_cost']<=.8*original['held_cost'])
    report['MDEA_held_out_scores']=scores
    report['qualified']=False
    report['limits']='MDEA selection only. Matched MEA, induced ketone, common-water, reactive/salt transfer and preservation/cost scoring remain unqualified.'
    (HERE/'results/e4.json').write_text(json.dumps(report,indent=2)+'\n')
    return 0

if __name__=='__main__':
    raise SystemExit(main())
