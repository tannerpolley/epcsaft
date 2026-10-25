"""#640 stage B: the frozen Li/Na/Mg experiment, using #575 and stage A readers."""
from __future__ import annotations

import argparse
import csv
from datetime import UTC, datetime
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

import epcsaft
from epcsaft import regression as r

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fit as li
import evaluate as checks
from trace_activity import Activities

HERE, ROOT = li.HERE, li.ROOT
MAIN = Path('/home/tnnrpolley21/Workspaces/Engineering/ePC-SAFT')
EVIDENCE = ROOT / 'build/issue-640/stageB'
RESULTS = HERE / 'results/stage-b'
DECL = HERE / 'inputs/stage-b-declaration.json'
IONS = ('lithium-cation', 'sodium-cation', 'magnesium-cation')
NAMES = ['segment_diameter', 'debye_huckel_diameter', 'born_diameter', 'k_ij/water', 'k_ij/chloride-anion']
SCALES = [0.9595302660294843, 1.1446504045203676, 1.1446504045203676,
          3.1879606522004598, 1.4082054395562396, 0.8853156407653622, 2.808541596118184]
OLD = ROOT / 'analyses/2026-chloride-brine-refit'
prior = li.load('chloride_fit', OLD / 'scripts/fit.py')
prior.BOUNDS['diameter_floor_angstrom']['lithium-cation'] = 1.38
for spec, scale in zip(prior.DECLARATION['training'] + li.D['training'], SCALES, strict=True):
    spec['scale'] = scale


def training(p):
    # #575's six-component composition helper is extended only by the zero Li slot.
    prior.composition = lambda mw, na, mg, ca, k=0: li.composition(p, {
        'sodium-cation': na, 'magnesium-cation': mg, 'calcium-cation': ca,
        'potassium-cation': k, 'chloride-anion': na+k+2*(mg+ca)})
    rows = prior.training(p) + li.training(p)
    for row in rows:
        row['molalities'] = dict(li= row.get('m_LiCl', 0.0), na=row['m_NaCl'],
                                mg=row.get('m_divalent', 0.0) if 'MgCl2' in row['salt'] else 0.0,
                                ca=row.get('m_divalent', 0.0) if 'CaCl2' in row['salt'] else 0.0,
                                k=row.get('m_KCl', 0.0))
    assert len(rows) == 517
    return rows


def coordinates(p):
    original = prior.CATIONS
    try:
        prior.CATIONS = IONS
        return prior.coordinates(p, NAMES)
    finally:
        prior.CATIONS = original


def document(physical):
    _, labels = coordinates(li.parameters())
    doc = json.loads((ROOT / f"data/packets/{li.D['base_packet']}/parameters/parameters.json").read_text())
    values = dict(zip(labels, physical, strict=True))
    for owner, coefficients in [(c['component_id'], c.get('coefficients', [])) for c in doc['components']] + [
            ('+'.join(sorted((q['component_id_a'], q['component_id_b']))), q['coefficients']) for q in doc['pairs']]:
        for coefficient in coefficients:
            label = f"{coefficient['family']}/{owner}"
            if label in values:
                coefficient['value']['magnitude'] = values.pop(label)
    assert not values
    return doc


def parameters(physical, components=li.COMPONENTS):
    return epcsaft.Parameters.from_mapping(document(physical), components=components)


def write(path, rows):
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def freeze():
    assert not DECL.exists(), 'the declaration is immutable'
    li.verify()
    EVIDENCE.mkdir(parents=True, exist_ok=True)
    RESULTS.mkdir(parents=True, exist_ok=True)
    design = MAIN / 'build/issue-640/design.md'
    assert li.sha(design) == '709610242b11a761daee40a08dc114a48e63892a41f82c7952dd4b7dae766dea'
    paths = {ROOT / name: digest for name, digest in li.D['inputs_sha256'].items()}
    paths.update({ROOT / name: digest for name, digest in prior.DECLARATION['inputs_sha256'].items()})
    for folder in (OLD / 'inputs', ROOT / 'analyses/2026-mixed-chloride-water-activity/inputs'):
        paths.update({path: li.sha(path) for path in folder.rglob('*') if path.is_file()})
    for path in (HERE / 'scripts').glob('*.py'):
        paths[path] = li.sha(path)
    paths.update({path: li.sha(path) for path in (OLD / 'scripts').glob('*.py')})
    paths.update({path: li.sha(path) for path in (ROOT / 'analyses/2026-mixed-chloride-water-activity/scripts').glob('*.py')})
    paths.update({path: li.sha(path) for path in (ROOT / f"data/packets/{li.D['base_packet']}").rglob('*') if path.is_file()})
    for name in ('inputs/fit-declaration.json', 'inputs/frozen-realizations.csv', 'inputs/lithium-f0-states.csv',
                 'results/pitzer.json', 'results/trace-divalent-song.json', 'results/trace-divalent-song.csv', 'results/fit.json',
                 'results/evaluation.json', 'results/trace-divalent.csv'):
        paths[HERE / name] = li.sha(HERE / name)
    external = [design, MAIN / 'build/issue-640/stageB/timing.json',
                MAIN / 'build/issue-640/stageB-readiness-r1.md',
                MAIN / 'build/issue-640/diagnosis/findings.md',
                MAIN / 'build/issue-640/diagnosis/inverse-discrepancy/results/fit.json']
    for name, digest in li.D['pdf_sha256'].items():
        matches = [p for p in (MAIN / 'build/issue-640').rglob(name) if li.sha(p) == digest]
        assert matches, name
        paths[matches[0]] = digest
    for name, digest in li.D['source_data_csv_sha256'].items():
        matches = [p for p in (MAIN / 'build/issue-640').rglob(name) if li.sha(p) == digest]
        assert matches, name
        paths[matches[0]] = digest
    paths.update({p: li.sha(p) for p in external})
    for path, digest in paths.items():
        assert li.sha(path) == digest, path
    p = li.parameters()
    cs, labels = coordinates(p)
    starts = []
    for sigma in (2.4, 3.2, 4.4, 6.0):
        for kw in (-0.5, 0.0):
            for kc in (-1.0, 0.5):
                starts.append({'id': f'grid-{len(starts):02}', 'values': [v for floor in (1.38, 2.04, 1.44)
                               for v in (sigma, max(.88*sigma, floor), max(.88*sigma, floor), kw, kc)]})
    starts.append({'id': 'anchored-16', 'values': [3.671853254763634, 13.56906735571716, 1.38, -2.983421843157207, -8.862058355359212,
                   3.5576592664836895, 5.694601927542711, 2.04, -3.102822266917611, -8.594439800307182,
                   5.045314479305746, 8.677630057903277, 1.44, -0.9641191526956133, -5.435688650150548]})
    d = {'run': 'issue-640-stage-B', 'owner_decision': '2026-10-19: Re-weight + stage B; D6=A; D4 frozen, no relaxation; D5',
         'design_sha256': li.sha(design), 'readiness': 'Ready; r1 R1 applied in B1; published #640 Stage B plan',
         'prior_evidence': 'Stage A and diagnostic held-out numbers existed before this newly declared experiment',
         'stage_A_declaration': li.sha(li.DECL), 'inherited_declaration': str(li.DECL.relative_to(ROOT)),
         'conditions': li.D['conditions'], 'base_packet': li.D['base_packet'],
         'coordinate_order': [f'{ion}/{name}' for ion in IONS for name in NAMES], 'physical_labels': labels,
         'bounds': prior.BOUNDS, 'coordinates': [{'lower': c.lower, 'upper': c.upper} for c in cs],
         'fixed': 'Everything except the 15 coordinates; packing=0.88*sigma; u/k Li=360, Na=230, Mg=1500 K. Ca,K,Cl,water fixed.',
         'model': 'Unchanged Figiel /1; complete dissociation, liquid; native ALG-FIT Ceres LM and exact columns',
         'training_blocks': prior.DECLARATION['training'] + li.D['training'], 'scales': SCALES,
         'scale_definition': 'sqrt(n/37)*delta_s/0.008736291849903472; old five delta=0.01; Li delta=0.008736291849903472; EG=sqrt((27*.021^2+18*.024^2)/45)',
         'scale_meaning': 'Fixed working discrepancy scales, not experimental standard deviations; unweighted errors decide checks',
         'training_rows': [{k: v for k, v in row.items() if k != 'observation'} for row in training(p)],
         'row_lists': li.D['row_lists'], 'source_records': {str(path.relative_to(ROOT)): li.read(path) for path in
                       [ROOT / name for name in li.D['inputs_sha256'] if name.endswith('.csv')]},
         'held_out_575': prior.DECLARATION['held_out'], 'starts': starts,
         'controls': dict(maximum_iterations=300, maximum_elapsed_time_seconds=900, shared_iterations=5100,
                          shared_seconds=2100, workers=0, OMP_NUM_THREADS=1, OPENBLAS_NUM_THREADS=1),
         'compute_limit_seconds': 2400, 'optimum_selection': li.D['optimum_selection'],
         'limits': li.D['limits'], 'F3': li.D['F3'], 'Hamer_Wu_guard': li.D['Hamer_Wu_guard'],
         'reference_model_aw_residual': li.D['reference_model_aw_residual'], 'mixing_residual': li.D['mixing_residual'],
         'rounding': li.D['rounding'], 'lithium': li.D['lithium'], 'missing': li.D['missing'],
         'F4_baseline_decades': [.04893214943109947, .009874431708778664, -.08736935514753164, -.010145231354341533],
         'r1408': 'No inherited /1 exception; require 1467/1467 under frozen source-agreement-limits.csv',
         'qualification': 'Full #575 and F1-Yao pass; F3,F4,Hamer-Wu,El Guendouzi stop/domain guards; F2/F2-prime/F5 targets only; F6 inherited limitation',
         'divalent_trace': {'backgrounds': ['MgCl2', 'CaCl2'], 'I_background_mol_kg': [1, 3, 6, 8.3],
                            'm_LiCl': 1e-6, 'meaning': 'Eight report-only predictions; native Song comparator for Li-Mg; no Li-Ca comparator'},
         'packet_domain': 'R6: Li-Mg passing Yao buckets I<=3 and 3<I<=8.3; Li-Na total<=6 conditional on guards. #575 held-out .45-8.03, trained .001-9.87; higher-I and ternary extrapolation disclosed; no Li-Ca,Li-Na-K,K-Ca claim.',
         'stop': 'Budget/native failure=incomplete evidence. Completed failure=new physics design outside this campaign. No extra starts, altered weights/limits, or new physics.',
         'inputs_sha256': {str(p.relative_to(ROOT)) if p.is_relative_to(ROOT) else str(p): h for p, h in paths.items()},
         'execution': {'checkout': str(ROOT), 'commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
                       'python': sys.executable, 'native_extension': epcsaft._core.__file__, 'native_sha256': li.sha(Path(epcsaft._core.__file__)),
                       'cpu_count': os.cpu_count(), 'HSL': 'unset', 'nice': 10},
         'frozen_at_utc': datetime.now(UTC).isoformat()}
    DECL.write_text(json.dumps(d, indent=2) + '\n')
    (HERE / 'inputs/stage-b-declaration.sha256').write_text(li.sha(DECL) + '\n')
    (HERE / 'inputs/stage-b-declaration.frozen-at').write_text(d['frozen_at_utc'] + '\n')
    shutil.copy(design, EVIDENCE / 'design.md')
    shutil.copy(DECL, EVIDENCE / DECL.name)
    print(d['frozen_at_utc'], li.sha(DECL), '517 rows, 15 coordinates, 17 starts', flush=True)


def verify(native=True):
    d = json.loads(DECL.read_text())
    assert li.sha(DECL) == (HERE / 'inputs/stage-b-declaration.sha256').read_text().strip()
    for name, digest in d['inputs_sha256'].items():
        path = ROOT / name
        if path == Path(__file__).resolve() and li.sha(path) != digest:
            path = HERE / 'inputs/stage-b-executed.py.txt'
        assert li.sha(path) == digest, name
    if native:  # the frozen run and evaluate need the frozen image; declarations check the current Engine
        assert li.sha(Path(epcsaft._core.__file__)) == d['execution']['native_sha256']
    return d


def declarations():
    d = verify(native=False)
    p = li.parameters()
    cs, _ = coordinates(p)
    rows = [{k: v for k, v in row.items() if k != 'observation'} for row in training(p)]
    expected = d['training_rows']
    erratum = HERE / 'inputs/stage-b-nacl-row-metadata.json'
    if erratum.exists():
        for correction in json.loads(erratum.read_text())['corrections']:
            expected[correction['index']] = correction['row']
    assert rows == expected
    assert len(cs) == 15 and all(len(s['values']) == 15 for s in d['starts'])


def run():
    d = verify()
    assert 'EPCSAFT_HSL_LIBRARY' not in os.environ
    assert os.environ['OMP_NUM_THREADS'] == os.environ['OPENBLAS_NUM_THREADS'] == '1'
    assert not (RESULTS / 'fit.json').exists(), 'retain the first frozen run'
    p = li.parameters()
    rows, (cs, labels), c = training(p), coordinates(p), d['controls']
    began = time.monotonic()
    (RESULTS / 'started-at.json').write_text(json.dumps({'monotonic': began, 'utc': datetime.now(UTC).isoformat()}) + '\n')
    result = r.fit(p, cs, [row['observation'] for row in rows],
                   controls=r.FitControls(maximum_iterations=c['maximum_iterations'],
                   maximum_elapsed_time_seconds=c['maximum_elapsed_time_seconds'], workers=c['workers']),
                   starts=[s['values'] for s in d['starts']], shared_iterations=c['shared_iterations'], shared_seconds=c['shared_seconds'])
    records = [{'start': s['id'], **prior.summary(v, labels, rows, cs), 'predictions': list(v.predictions),
                'iteration_seconds': list(v.iteration_seconds), 'iteration_costs': list(v.iteration_costs),
                'solver_options': v.solver_options, 'residual_evaluations': v.residual_evaluations,
                'jacobian_evaluations': v.jacobian_evaluations,
                'failures': [{'trial': f.trial, 'observation': f.observation, 'optimizer': list(f.optimizer),
                              'status': str(f.status), 'message': f.message, 'jacobian': f.jacobian} for f in v.failures],
                'observation_statuses': [str(q) for q in v.statuses]} for s, v in zip(d['starts'], result)]
    record = {'declaration_sha256': li.sha(DECL), 'wall_seconds': time.monotonic()-began, 'starts': records,
              'optima': prior.optima(records), 'labels': labels, 'native_failures': sum(v['native_failures'] for v in records)}
    (RESULTS / 'fit.json').write_text(json.dumps(record, indent=2) + '\n')
    print('Fit finished', record['wall_seconds'], [(v['start'], v['status'], v['final_cost']) for v in records], flush=True)
    evaluate(record, began)


def evaluate(record, began):
    # Load the unchanged #575 amendment reader with its own fit module.
    sys.modules['fit'] = prior
    oldcheck = li.load('chloride_checks', OLD / 'scripts/evaluate.py')
    sys.modules['fit'] = li
    compare = prior.comparison()
    scenarios = compare.read(compare.INPUTS / 'frozen-realizations.csv')
    limits = {(row['realization_id'], row['I_bucket']): {k: float(v) if k.startswith(('limit_', 'disagreement_')) else v
              for k, v in row.items()} for row in li.read(OLD / 'results/source-agreement-limits.csv')}
    assert len(scenarios) == 1467 and len(limits) == 2934
    baseline_transfers, _ = checks.transfers(li.parameters())
    pitzer = json.loads((HERE / 'results/pitzer.json').read_text())
    pv = {row['state']: checks.combinations(row['log10_activity']) for row in pitzer['states']}
    pt = {f'{state}/{name}': pv[state][name]-pv['fig7'][name] for state in ('table1', 'ms2') for name in pv['fig7']}
    fn, fr, fb = checks.held_out(checks.Water(li.parameters()), 'Robinson-Covington', 'unchanged-K-Ca')
    write(RESULTS / 'k-rounding.csv', fr)
    write(RESULTS / 'k-buckets.csv', fb)
    summaries = []
    for index, optimum in enumerate(o for o in record['optima'] if o['within_band']):
        if time.monotonic()-began >= 2400:
            raise TimeoutError('total scientific compute budget exhausted; qualification incomplete')
        label = f'optimum-{index+1}'
        target = RESULTS / label
        target.mkdir()
        p = parameters(optimum['physical'])
        water = checks.Water(p)
        yn, yr, yb = checks.held_out(water, 'Yao', label)
        for name, data in [('yao-nominal.csv', yn), ('yao-rounding.csv', yr), ('yao-buckets.csv', yb)]:
            write(target / name, data)
        trace = Activities(parameters(optimum['physical'], ('water', 'lithium-cation', 'sodium-cation', 'chloride-anion')))
        f3 = [{**trace.row(m), 'binary_error_decade': trace.binary_check(m)} for m in (1, 3, 6)]
        assert max(abs(v['binary_error_decade']) for v in f3) <= 1e-10
        write(target / 'trace-na.csv', f3)
        engine, diagnostics = checks.transfers(p)
        lr = [{'state': key.split('/')[0], 'check': key.split('/')[1], 'engine_transfer': value,
               'pitzer_transfer': pt[key], 'difference_decade': value-pt[key],
               'baseline_difference': baseline_transfers[key]-pt[key], 'change_from_baseline': value-baseline_transfers[key],
               'target_met': abs(value-pt[key]) <= .15} for key, value in engine.items()]
        write(target / 'lithium-transfers.csv', lr)
        write(target / 'lithium-density-checks.csv', diagnostics)
        predictions = []
        for row in training(p):
            value = water.predict(**row['molalities'])
            m = row['molalities']
            q = 2*(m['li']+m['na']+m['k'])+3*(m['mg']+m['ca'])
            obs_aw = math.exp(-p.molar_masses[0]*q*row['observed'])
            predictions.append({**{k: v for k, v in row.items() if k not in ('observation', 'molalities')},
                                **m, **value, 'residual': value['phi']-row['observed'], 'aw_residual': value['aw']-obs_aw})
        write(target / 'training.csv', predictions)
        guard = checks.metrics([v['residual'] for v in predictions if v['source'] == 'Hamer-Wu'])
        eg = {b: checks.metrics([v['aw_residual'] for v in predictions if v['source'] == 'El Guendouzi'
              and ('I<=3' if v['li']+v['na'] <= 3 else 'I>3') == b]) for b in ('I<=3', 'I>3')}
        divalent = []
        retained = li.read(HERE / 'results/trace-divalent.csv')
        for salt, ion in [('MgCl2', 'magnesium-cation'), ('CaCl2', 'calcium-cation')]:
            activity = Activities(parameters(optimum['physical'], ('water', 'lithium-cation', ion, 'chloride-anion')), salt)
            for ionic in (1., 3., 6., 8.3):
                previous = next(v for v in retained if v['background'] == salt and float(v['I_background_mol_kg']) == ionic)
                error = activity.binary_check(ionic/3)
                assert abs(error) <= 1e-10
                divalent.append({**previous, 'stage_B_log_gamma_trace': activity.trace(ionic/3), 'stage_B_binary_error_decade': error})
        write(target / 'trace-divalent.csv', divalent)
        # Row construction reloads #556; the pool must resolve this comparison's functions.
        sys.modules['compare'] = compare
        verdicts, failed = oldcheck.evaluate(compare, target, EVIDENCE / label, document(optimum['physical']), scenarios, limits)
        f1 = all(v['all_passed'] for v in yb if v['scored'])
        f3ok = all(abs(v['theta_eff']) <= .04 for v in f3)
        f4ok = all(v['target_met'] and abs(v['change_from_baseline']) <= .01 for v in lr if v['check'].startswith('F4'))
        egok = all(checks.passes(eg[b], li.D['reference_model_aw_residual'][b]) for b in eg)
        summary = {'optimum': label, 'starts': optimum['starts'], 'cost': optimum['final_cost'], 'physical': optimum['physical'],
                   'optimizer': optimum['optimizer'], 'active_bounds': optimum['active_bounds'],
                   'numerical_box_coordinates': optimum['numerical_box_coordinates'], 'rank': optimum['rank'],
                   'covariance_status': optimum['covariance_status'], 'F575_passing_realizations': 1467-len(failed),
                   'F575_verdicts': verdicts, 'F575_failed_realizations': sorted(failed), 'F1_Yao': yb, 'F3': f3,
                   'Hamer_Wu': guard, 'El_Guendouzi_aw': eg, 'lithium': lr,
                   'F575_passed': not failed, 'F1_passed': f1, 'F3_passed': f3ok, 'F4_passed': f4ok,
                   'Hamer_Wu_passed': checks.passes(guard, li.D['Hamer_Wu_guard']), 'El_Guendouzi_passed': egok,
                   'set8_omission': checks.metrics([v['mixing_residual'] for v in yn if v['bucket'] == '3<I<=8.3' and v['set'] != '8']),
                   'elapsed_scientific_seconds': time.monotonic()-began}
        summary['qualified'] = all(summary[k] for k in ('F575_passed', 'F1_passed', 'F3_passed', 'F4_passed', 'Hamer_Wu_passed', 'El_Guendouzi_passed'))
        (target / 'evaluation.json').write_text(json.dumps(summary, indent=2) + '\n')
        summaries.append(summary)
        print(label, 'qualified', summary['qualified'], '575', summary['F575_passing_realizations'], 'Yao', f1, flush=True)
    complete = time.monotonic()-began <= 2400 and len(record['starts']) == 17 and bool(summaries)
    final = {'selected': summaries[0] if summaries else None, 'optima': summaries, 'complete': complete,
             'scientific_wall_seconds': time.monotonic()-began, 'fit_wall_seconds': record['wall_seconds'],
             'packet_created': False, 'native_failures': record['native_failures'], 'F6': fb}
    (RESULTS / 'evaluation.json').write_text(json.dumps(final, indent=2) + '\n')
    assert complete, 'incomplete scientific evidence'


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('freeze', 'run', 'evaluate'))
    action = parser.parse_args().action
    if action == 'freeze':
        freeze()
    elif action == 'run':
        run()
    else:
        verify()
        began = json.loads((RESULTS / 'started-at.json').read_text())['monotonic']
        assert time.monotonic()-began < 2400, 'the original compute budget is closed'
        evaluate(json.loads((RESULTS / 'fit.json').read_text()), began)
