"""Declare the frozen HCl study and transport native results; never solve in Python."""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import csv
import struct
import time
from epcsaft import Parameters, provenance

ROOT = Path(__file__).resolve().parents[3]
LEAF = Path(__file__).resolve().parents[1]
INPUTS = LEAF / 'inputs'
DECLARATION = Path(sys.argv[1]) if len(sys.argv) > 1 and sys.argv[1] != '--declare-only' else INPUTS / 'fit-declaration-issue-798-r2.json'
D = json.loads(DECLARATION.read_text())
C = D.get('trial', dict(ion='hydronium-cation', label='H', source='Hamer–Wu', training='hamer-wu-hcl.csv', split=3, cutoff=7, low_states=23, stoichiometry=[1, 1], candidate='candidate.json', build='build/issue-798'))
LEAF = DECLARATION.resolve().parents[1]
INPUTS = LEAF / 'inputs'

def declarations():
    if '--packet' in sys.argv:
        return Parameters.from_json(sys.argv[sys.argv.index('--packet') + 1])
    document = json.loads((INPUTS / C['candidate']).read_text())
    if 'trial' in D and not D.get('hcl_mixture_training'):
        return Parameters.from_mapping(document)
    document['sources'].append({'source_id': 'issue-798-declaration', 'citation': 'Issue 798 frozen design r2, section 1; Held 2008 Table 2 and explicit diameter assumptions'})
    hydrogen = document['components'][7]
    hydrogen['fixed']['charge_number']['provenance']['locator'] = 'Table 3, H+ ion label'
    for coefficient in hydrogen['coefficients']:
        if coefficient['family'] in D['H_fixed']:
            coefficient['value']['magnitude'] = D['H_fixed'][coefficient['family']]
            coefficient['provenance']['source_id'] = 'issue-798-declaration'
            coefficient['provenance']['locator'] = 'Issue 798 design r2 section 1: Held 2008 Table 2; fixed diameter conventions'
    document['pairs'][-1]['coefficients'][0]['value']['magnitude'] = D['H_Cl_k']
    document['pairs'][-1]['coefficients'][0]['provenance'] = {'source_id': 'issue-798-declaration', 'locator': D['H_Cl_convention']}
    return Parameters.from_mapping(document)

def main():
    began = time.monotonic()
    if '--profile' in sys.argv:
        from datetime import datetime
        from zoneinfo import ZoneInfo
        now = datetime.now(ZoneInfo('America/Denver'))
        allowed = now <= datetime(2026,10,28,22,tzinfo=now.tzinfo) or now >= datetime(2026,10,29,8,tzinfo=now.tzinfo)
        if '--owned' not in sys.argv:
            assert allowed, 'profile launch outside accepted window'
            command = ['systemd-run','--user','--unit=issue-901-profile','--collect','--wait','--pipe',
                       '-p','RuntimeMaxSec=3h','-p','TimeoutStopSec=10s','-p','KillMode=control-group',
                       '--working-directory='+str(ROOT),'env','-u','EPCSAFT_HSL_LIBRARY',
                       'EPCSAFT_BUILD_JOBS=6','PYTHONDONTWRITEBYTECODE=1','nice','-n','10',
                       sys.executable,str(Path(__file__).resolve()),*sys.argv[1:],'--owned']
            subprocess.run(command,check=True)
            return
        destination = Path(sys.argv[2])
        destination.mkdir(parents=True,exist_ok=False)
        (destination/'activation.json').write_text(json.dumps(dict(clock=now.isoformat(),allowed=allowed,RuntimeMaxSec=10800))+'\n')
        if not allowed:
            raise SystemExit('profile activation refused outside accepted window')
    assert 'EPCSAFT_HSL_LIBRARY' not in os.environ
    for name, digest in D['input_hashes'].items():
        assert hashlib.sha256((INPUTS / name).read_bytes()).hexdigest() == digest, name
    results = Path(sys.argv[2]).resolve() if len(sys.argv) > 2 else LEAF / 'results/issue-798-r2'
    results.mkdir(parents=True, exist_ok=True)
    assert not any(p.name != 'activation.json' for p in results.iterdir()), 'retain the first frozen execution'
    p = declarations()._document
    build = ROOT / ('build/issue-901' if D.get('profile') else C['build'])
    build.mkdir(parents=True, exist_ok=True)
    with (build / 'inputs.hpp').open('w') as f:
        f.write('const bool packet_readback = ' + str('--packet' in sys.argv).lower() + ';\n')
        f.write('#include "surface/parameter_builder.hpp"\nepcsaft::Parameters declared_parameters() { epcsaft::ParameterDocument d;\n')
        for attr in ['component_ids', 'site_component', 'site_multiplicity']:
            f.write(f'd.{attr} = ' + json.dumps(list(getattr(p, attr))).replace('[', '{').replace(']', '}') + ';\n')
        for attr in ['base_choice', 'association_choice', 'electrolyte_choice', 'ionic_dispersion_choice', 'permittivity_choice', 'polar_choice']:
            f.write(f'd.{attr} = ' + json.dumps(getattr(p, attr)) + ';\n')
        for e in p.entries:
            f.write('d.entries.push_back({' + ','.join([json.dumps(e.identity), json.dumps(e.field), '{' + ','.join(map(str, e.index)) + '}', repr(e.magnitude), json.dumps(e.unit)]) + '});\n')
        f.write('return epcsaft::build_parameters(d); }\n')
        f.write('const std::vector<std::vector<double>> starts = ' + json.dumps(D['starts']).replace('[', '{').replace(']', '}') + ';\n')
        f.write('const std::vector<std::vector<double>> screen_points = ' + json.dumps(D.get('screen_points', D['starts'])).replace('[', '{').replace(']', '}') + ';\n')
        f.write(f"const double screening_max_correlation = {D['identification'].get('max_correlation', 1)};\n")
        coordinates = D['coordinates'] if 'coordinates' in D else [D['coordinate']]
        f.write('#include "regression/coordinates.hpp"\nconst std::vector<epcsaft::regression::Coordinate> coordinate_declarations = {\n')
        for c in coordinates:
            targets = c.get('targets', [{'family': 'Kij', 'components': c.get('pair'), 'coefficient': 1}])
            f.write('{{')
            for t in targets:
                indices = [list(p.component_ids).index(name) for name in t['components']]
                if t['family'] == 'Kij': indices = sorted([indices[0]*len(p.component_ids)+indices[1], indices[1]*len(p.component_ids)+indices[0]])
                f.write('{epcsaft::ParameterFamily::' + t['family'] + ',{' + ','.join(map(str, indices)) + '},' + str(t['coefficient']) + '},')
            f.write('},' + ','.join(map(str, [c['origin'], c['scale'], *c['bounds']])) + ',std::nullopt},\n')
        f.write('};\nconst std::vector<std::string> coordinate_names = ' + json.dumps([c.get('name', C['label'] + ('-water' if 'water' in c.get('pair', []) else '-Cl')) for c in coordinates]).replace('[', '{').replace(']', '}') + ';\n')
        f.write(f"const bool hcl_mixture_training = {str(D.get('hcl_mixture_training', False)).lower()};\n")
        f.write(f"const bool sulfate_training = {str(D.get('sulfate', False)).lower()};\n")
        f.write(f"const bool verify_lengths = {'true' if 'verification' in D else 'false'};\n")
        for name, value in D.get('verification', dict(equivalence=1e-10, difference_tolerance=1e-5, steps=[1e-4, 5e-5], point=6)).items():
            f.write('const auto verification_' + name + ' = ' + ('std::vector<double>' + json.dumps(value).replace('[', '{').replace(']', '}') if isinstance(value, list) else str(value)) + ';\n')
        config = {**{k: D[k] for k in ['temperature_K', 'pressure_Pa', 'maximum_iterations', 'maximum_seconds', 'total_maximum_iterations', 'total_maximum_seconds', 'training_states']}, **{k: C[k] for k in ['split', 'cutoff', 'low_states']}, 'ion': list(p.component_ids).index(C['ion']), 'water': list(p.component_ids).index('water'), 'chloride': list(p.component_ids).index('chloride-anion')}
        config.update(training_observations=D.get('training_observations', 2 * D['training_states']), mixture_states=D.get('mixture', {}).get('states', 0), sodium=list(p.component_ids).index('sodium-cation'))
        config.update(lithium=list(p.component_ids).index('lithium-cation'), magnesium=list(p.component_ids).index('magnesium-cation'))
        for name, value in {**config, **{f'identification_{k}': v for k, v in {'relative_resolution': 1e-3, **D['identification']}.items()}, 'nu_cation': C['stoichiometry'][0], 'nu_anion': C['stoichiometry'][1]}.items():
            f.write(f'const auto {name} = {value};\n')
        f.write('const std::string training = ' + json.dumps(C['training']) + ', ion_label = ' + json.dumps(C['label']) + ';\n')
        f.write('const std::string mixture_training = ' + json.dumps(D.get('mixture', {}).get('training', '')) + ';\n')
        f.write('const std::string mixture_reference = ' + json.dumps(D.get('mixture', {}).get('reference', '')) + ';\n')
        f.write('const bool profile_enabled = ' + str(bool(D.get('profile'))).lower() + ';\n')
        for name, value in {**dict(coarse=[],refinement_offsets=[],starts=[],bounds=[],mean_limits=[],max_limits=[],target_mean=0,target_max=0,anchor=0,compensation=0),**D.get('profile', {})}.items():
            if isinstance(value, list):
                f.write('const std::vector<double> profile_' + name + ' = ' + json.dumps(value).replace('[','{').replace(']','}') + ';\n')
            elif isinstance(value, (int,float)):
                f.write('const double profile_' + name + ' = ' + str(value) + ';\n')
    shared = Path(__file__).resolve().parent
    (build / 'activity.cpp').write_text(f'#include "inputs.hpp"\n#include "{shared / "activity.cpp"}"\n')
    executable = subprocess.check_output(['tools/native-driver', 'regression_checks', str(build / 'activity.cpp')], text=True).strip()
    library = ROOT / D['native_library']['path']
    native = subprocess.check_output(['git', 'ls-files', 'engine/native', 'engine/src/epcsaft'], cwd=ROOT, text=True).splitlines()
    sources = native + [str(path.relative_to(ROOT)) for path in shared.glob('*') if path.is_file()]
    hashes = {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest() for name in sources}
    assert not subprocess.check_output(['git', 'diff', 'HEAD', '--', 'engine/native', 'engine/src'], cwd=ROOT)
    identity = {'source': provenance(), 'declaration_sha256': hashlib.sha256(DECLARATION.read_bytes()).hexdigest(), 'native_library_sha256': hashlib.sha256(library.read_bytes()).hexdigest(), 'executable_sha256': hashlib.sha256(Path(executable).read_bytes()).hexdigest(), 'source_sha256': hashes}
    identity['parameter_transport_sha256'] = hashlib.sha256((build / 'inputs.hpp').read_bytes()).hexdigest()
    identity['source_commit'] = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()
    (results / 'native-identity.json').write_text(json.dumps(identity, default=str, indent=2) + '\n')
    environment = {**os.environ, 'LD_LIBRARY_PATH': str(library.parent)}
    if D.get('profile') and '--profile' in sys.argv:
        identity.update(base='1bb4c744c3ff67a0b6c4dc0d95c718f318e18055', inputs=D['input_hashes'], component_order=list(p.component_ids), command=sys.argv)
        (results/'native-identity.json').write_text(json.dumps(identity,default=str,indent=2)+'\n')
        profile_began=time.monotonic()
        try:
            subprocess.run([executable,str(LEAF),str(results),'profile'],env=environment,check=True)
        finally:
            (results/'wall-time.json').write_text(json.dumps(dict(preparation_seconds=profile_began-began,profile_process_seconds=time.monotonic()-profile_began))+'\n')
        return
    if '--packet' in sys.argv:
        return
    admission = results / D.get('mixture', {}).get('retained_admission', 'mixture-admission.csv')
    if 'mixture' in D:
        subprocess.run([executable, str(LEAF), str(results), 'prepare'], env=environment, check=True)
        admission_hash = hashlib.sha256((results / 'mixture-admission.csv').read_bytes()).hexdigest()
        if admission != results / 'mixture-admission.csv':
            assert hashlib.sha256(admission.read_bytes()).hexdigest() == admission_hash
            (results / 'mixture-admission.csv').unlink()
        (results / 'admission-sha256.json').write_text(json.dumps({str(admission.relative_to(results)): admission_hash, **D['input_hashes']}, indent=2) + '\n')
        assert json.loads((results / 'admission-sha256.json').read_text())[str(admission.relative_to(results))] == admission_hash
    if 'trial' in D:
        prepared = time.monotonic()
        prefit = subprocess.run([executable, str(LEAF), str(results), 'prefit'], env=environment)
        (results / 'screen-timing.json').write_text(json.dumps({'preparation_seconds': prepared - began, 'screen_process_seconds': time.monotonic() - prepared}) + '\n')
        columns = list(csv.DictReader((results / 'sensitivity.csv').open()))
        hashes = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in [*results.glob('sensitivity*.csv'), results / 'prefit-diagnostics.csv', *results.glob('length-*.csv')]}
        if 'verification' in D:
            diagnostics = list(csv.DictReader((results / 'prefit-diagnostics.csv').open()))
            assert {(row['evaluation'], int(row['observation'])) for row in diagnostics} == {(str(j), i) for j in range(1, len(D.get('screen_points', D['starts'])) + 2 * len(coordinates) * len(D['verification']['steps']) + 1) for i in range(D['training_observations'])}
            assert (results / 'length-verification.csv').stat().st_size > 100
        if D.get('mixture', {}).get('reference'):
            diagnostics = list(csv.DictReader((results / 'prefit-diagnostics.csv').open()))
            assert {(int(row['evaluation']), int(row['observation'])) for row in diagnostics} == {(j, i) for j in range(1, len(D['screen_points']) + 1) for i in range(D['training_observations'])}
            assert all(row['status'] == 'available' and float(row['density_mol_m3']) > 0 and float(row['pressure_slope']) > 0 for row in diagnostics)
            hashes.update({p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in results.glob('reference-*.csv')})
        points = D.get('screen_points', D['starts'])
        assert len(columns) == len(points) * D.get('training_observations', 2 * D['training_states']) * len(coordinates)
        if 'mixture' in D:
            assert hashlib.sha256(admission.read_bytes()).hexdigest() == admission_hash
        if D.get('sulfate'):
            hashes['mixture-admission.csv'] = hashlib.sha256((results / 'mixture-admission.csv').read_bytes()).hexdigest()
        hashes['ordered_binary64_columns'] = {f'{start}/{a}/{field}': hashlib.sha256(b''.join(struct.pack('<d', float(row[field])) for row in columns if int(row['start']) == start and int(row['coordinate']) == a)).hexdigest() for start in range(1, len(points) + 1) for a in range(len(coordinates)) for field in ['prediction_derivative', 'b', 'A']}
        if D.get('sulfate'):
            physical = list(csv.DictReader((results / 'length-verification.csv').open()))
            hashes['ordered_joint_physical_binary64'] = {str(j): hashlib.sha256(b''.join(struct.pack('<d', float(row['physical_partial_sum'])) for row in physical if int(row['evaluation']) == j and float(row['step']) == 0)).hexdigest() for j in range(1, len(points)+1)}
            assert all(sum(int(row['evaluation']) == j and float(row['step']) == 0 for row in physical) == D['training_observations'] for j in range(1, len(points)+1))
        hashes['identity'] = json.loads((results / 'native-identity.json').read_text())
        (results / 'sensitivity-sha256.json').write_text(json.dumps(hashes, indent=2) + '\n')
        assert json.loads((results / 'sensitivity-sha256.json').read_text()) == hashes
        assert all(hashlib.sha256((results / name).read_bytes()).hexdigest() == digest for name, digest in hashes.items() if name.endswith('.csv'))
        assert all(hashlib.sha256((ROOT / name).read_bytes()).hexdigest() == digest for name, digest in identity['source_sha256'].items())
        assert hashlib.sha256(DECLARATION.read_bytes()).hexdigest() == identity['declaration_sha256']
        assert all(hashlib.sha256((INPUTS / name).read_bytes()).hexdigest() == digest for name, digest in D['input_hashes'].items())
        assert hashlib.sha256((build / 'inputs.hpp').read_bytes()).hexdigest() == identity['parameter_transport_sha256']
        assert hashlib.sha256(Path(executable).read_bytes()).hexdigest() == identity['executable_sha256']
        assert hashlib.sha256(library.read_bytes()).hexdigest() == identity['native_library_sha256']
        summary = list(csv.DictReader((results / 'sensitivity-summary.csv').open()))
        assert len(summary) == len(points)
        if any(row['passes'] != '1' for row in summary):
            assert prefit.returncode == 1
            return
        prefit.check_returncode()
    subprocess.run([executable, str(LEAF), str(results)] + ([sys.argv[3][2:]] if len(sys.argv) > 3 else []), env=environment, check=True)
    assert hashlib.sha256(library.read_bytes()).hexdigest() == identity['native_library_sha256']

if __name__ == '__main__':
    if '--declare-only' in __import__('sys').argv:
        declarations()
    else:
        main()
