"""Transport the retained restricted Belov inputs to the native declaration/solve."""
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import epcsaft as ep
LEAF=Path(__file__).resolve().parents[1]
ROOT=LEAF.parents[2]

def declarations():
    return json.loads((LEAF/'inputs/cases/test-01/source/source.json').read_text())

def main():
    source=declarations()
    mass_file=ROOT/'analyses/reproduction/1998-chase/inputs/source/molar-masses.csv'
    masses={r['species']:r['molar_mass_g_per_mol'] for r in csv.DictReader(mass_file.open())}
    build=ROOT/'build/issue-789/2023-belov';build.mkdir(parents=True,exist_ok=True)
    meson=ROOT/'build/engine-meson/shared'
    includes=['engine/native','engine/build/subprojects/CppAD-20250000.3/include','engine/build/subprojects/CppAD-20250000.3','build/engine-meson/shared/build/subprojects/CppAD-20250000.3','/usr/include/eigen3']
    driver=build/'equilibrium'
    command=['nice','-n','10','c++','-std=c++17',*[f'-I{p}' for p in includes]]
    command += [LEAF/'scripts/equilibrium.cpp',f'-L{meson}',f'-Wl,-rpath,{meson}','-lepcsaft_native','-o',driver]
    subprocess.run(command,cwd=ROOT,check=True)
    fields='temperature_k pressure_pa source_standard_pressure_pa'.split()
    header=[source[k] for k in fields]+[masses['CO'],masses['CO2'],source['log_amount_absolute_tolerance']]
    columns=[source['species_ids'],source['source_standard_gibbs_over_rt'],*source['balance_matrix'],source['feed_amounts_mol'],source['independent_expected_amounts_mol']]
    request=' '.join(map(str,header))+'\n'
    request+=''.join(' '.join(map(str,row))+'\n' for row in zip(*columns))
    request+='\n'.join(' '.join(map(str,r)) for r in source['reaction_matrix'])+'\n'
    result=subprocess.run([driver,LEAF/'results'],input=request,text=True,capture_output=True)
    (build/'console.txt').write_text(result.stdout+result.stderr)
    library=meson/'libepcsaft_native.so'
    identity=dict(engine=list(ep.provenance()),native_library=str(library),mass_source_rows=['CO','CO2'])
    identity['binaries_sha256']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (library,driver)}
    identity['source_sha256']={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [*sorted((LEAF/'inputs').rglob('*')),mass_file] if p.is_file()}
    (LEAF/'results/run-identity.json').write_text(json.dumps(identity,indent=2)+'\n')
    result.check_returncode()

if __name__=='__main__':
    main()
