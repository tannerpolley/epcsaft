"""Compile the retained native falsifiers and reevaluate the frozen endpoints."""
import argparse
import json
from pathlib import Path
import shlex
import shutil
import subprocess

ROOT=Path(__file__).resolve().parents[3]
HERE=Path(__file__).resolve().parent
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--out',required=True,type=Path)
args=parser.parse_args()
out=args.out.resolve();out.mkdir(parents=True,exist_ok=True)
build=ROOT/'build/engine-meson/shared'
entry=next(e for e in json.loads((build/'compile_commands.json').read_text()) if e['file'].endswith('native/equilibrium/gibbs_minimizer.cpp'))
command=[];skip=False
for value in shlex.split(entry['command']):
    if skip:skip=False;continue
    if value in ('-o','-MF','-MQ'):skip=True;continue
    if 'ccache' in value or value in ('-c','-MD') or value.endswith('gibbs_minimizer.cpp'):continue
    command.append(value)
for name in ('forward-error','association-floor','row-producers','compare','previous-result'):
    subprocess.run(command+['-O2','-I'+str(ROOT/'engine/tests/native'),str(HERE/(name+'.cpp')),
        '-L'+str(build),'-lepcsaft_native','-Wl,-rpath,'+str(build),'-lquadmath','-o',str(out/name)],cwd=entry['directory'],check=True)
for name in ('original-criteria.csv','binary128-roots.csv'):
    shutil.copy2(HERE.parent/'inputs'/name,out/name)
with (HERE.parent/'inputs/endpoints.txt').open() as inputs,(out/'estimates.csv').open('w') as values:
    subprocess.run([str(out/'forward-error'),str(out/'row-floors.csv'),str(out/'packing.csv'),str(out/'snapshot-coordinates.csv')],stdin=inputs,stdout=values,check=True)
for name in ('association-floor','row-producers','previous-result'):
    with (out/(name+'.csv')).open('w') as values:
        subprocess.run([str(out/name)],stdout=values,check=True)
subprocess.run([str(out/'compare'),str(out)],check=True)
