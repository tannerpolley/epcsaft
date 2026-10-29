"""Declare #902 halite inputs; compile and read the native study without Python thermodynamics."""
import csv
import hashlib
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
from datetime import datetime
from zoneinfo import ZoneInfo
from epcsaft import Parameters, _core

ROOT = Path(__file__).resolve().parents[3]
STUDY = Path(__file__).resolve().parents[1]
INPUT = STUDY / 'inputs/figiel-902'
OUT = STUDY / 'results/figiel-902'
PACKET = ROOT / 'data/packets/figiel-2025-chloride-mixture-refit/2/parameters/parameters.json'
SPECIES = ['water', 'sodium-cation', 'chloride-anion']

def declarations():
    return Parameters.from_mapping(json.loads(PACKET.read_text()), components=SPECIES)

if __name__ == '__main__':
    parameters = declarations()
    if '--declare-only' in sys.argv:
        raise SystemExit(0)
    clock = datetime.now(ZoneInfo('America/Denver'))
    if 2 <= clock.hour < 8 and (clock.hour != 2 or clock.minute >= 15):
        raise SystemExit('No compiles or solves during the host timing window / pre-window cutoff')
    if ROOT not in Path(_core.__file__).resolve().parents:
        raise RuntimeError('Engine extension is outside this checkout')
    scratch = ROOT / 'build/issue-902/halite-build-r1'
    scratch.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    doc = parameters._document
    quote = json.dumps
    vector = lambda values: '{' + ','.join(map(str, values)) + '}'
    code = ['#pragma once', 'epcsaft::ParameterDocument input_document() {', 'epcsaft::ParameterDocument d;']
    code += ['d.component_ids={' + ','.join(map(quote, doc.component_ids)) + '};']
    for field in ['site_component', 'site_multiplicity']:
        code += [f'd.{field}={vector(getattr(doc, field))};']
    for field in ['base_choice', 'association_choice', 'electrolyte_choice', 'ionic_dispersion_choice', 'permittivity_choice', 'polar_choice']:
        code += [f'd.{field}={quote(getattr(doc, field))};']
    for e in doc.entries:
        code += ['d.entries.push_back({' + ','.join([quote(e.identity), quote(e.field), vector(e.index), repr(e.magnitude), quote(e.unit)]) + '});']
    code += ['return d;', '}']
    rows = {r['id']: r for r in json.loads((INPUT/'frozen-observations.json').read_text())['rows']}
    law, cg, pm = (rows[k] for k in ['clarke1985-table34-25C', 'clarke1985-table17A-25C', 'pinho2005-ThermoML-1-1'])
    fixed = dict(reference_temperature=law['T_K'], reference_pressure=law['P_Pa'], dG=law['delta_g_J_mol'],
                 dH=law['delta_h_J_mol'], dCp=law['delta_cp_J_mol_K'], g_se=law['se_g_J_mol'],
                 measured_m=cg['value'], measured_m_se=cg['standard_error'], pinho_pressure=pm['P_Pa'],
                 measured_w=pm['value'], measured_w_u=pm['expanded_uncertainty'], salt_mass=0.05844277, solid_volume=27.1e-6)
    code += [f'constexpr double {key}={value};' for key, value in fixed.items()]
    (scratch/'figiel_inputs.hpp').write_text('\n'.join(code)+'\n')
    build = Path((ROOT/'build/environment-wheel/meson-build-dir').read_text().strip())
    commands = json.loads((build/'compile_commands.json').read_text())
    command = next(c for c in commands if c['file'].endswith('/tests/native/surface_example.cpp'))
    args = shlex.split(command['command'])
    args[args.index('-o')+1] = str(scratch/'halite.o')
    args[-1] = str(STUDY/'scripts/figiel_halite.cpp')
    subprocess.run(['nice', '-n', '10', *args, '-I'+str(scratch)], cwd=build, check=True)
    link = subprocess.check_output(['ninja', '-C', str(build), '-t', 'commands', 'surface_example'], text=True).splitlines()[-1]
    args = shlex.split(link)
    args[args.index('-o')+1] = str(scratch/'halite')
    args = [str(scratch/'halite.o') if a.endswith('tests_native_surface_example.cpp.o') else a for a in args]
    subprocess.run(['nice', '-n', '10', *args], cwd=build, check=True)
    run = subprocess.run(['nice', '-n', '10', str(scratch/'halite'), str(OUT)], cwd=build, capture_output=True, text=True,
                         env={**os.environ, 'LD_LIBRARY_PATH': str(build)})
    (OUT/'native.stdout').write_text(run.stdout)
    (OUT/'native.stderr').write_text(run.stderr)
    records = list(csv.DictReader((OUT/'attempts.tsv').open(), delimiter='\t'))
    digest = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    paths = [PACKET, INPUT/'frozen-observations.json', INPUT/'source-hashes.json', Path(_core.__file__),
             build/'libepcsaft_native.so', STUDY/'scripts/figiel_halite.cpp', Path(__file__), scratch/'figiel_inputs.hpp']
    identity = dict(base='7f2a1d7d24b79970ab5fe8b0bbc159ede553de31', commit=subprocess.check_output(['git','rev-parse','HEAD'], text=True).strip(),
                    timestamp=clock.isoformat(), fixed=fixed, species=SPECIES, returncode=run.returncode,
                    hashes={str(p.relative_to(ROOT)): digest(p) for p in paths}, attempts=len(records),
                    failures=[r for r in records if r['success'] != '1' or r['numerical_pass'] != '1'])
    (OUT/'run-identity.json').write_text(json.dumps(identity, indent=2)+'\n')
    print(json.dumps(identity, indent=2))
    if run.returncode:
        raise SystemExit(run.returncode)
