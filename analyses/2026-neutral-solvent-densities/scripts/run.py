"""Build and run the native pure-solvent study in this checkout only."""
import hashlib
import json
import os
from pathlib import Path
import shlex
import subprocess
from datetime import datetime
from zoneinfo import ZoneInfo
from epcsaft import _core

ROOT = Path(__file__).resolve().parents[3]
STUDY = Path(__file__).resolve().parents[1]

def declarations():
    return (STUDY/'inputs/parameters.tsv').read_text()

if __name__ == '__main__':
    clock = datetime.now(ZoneInfo('America/Denver'))
    if 2 <= clock.hour < 8 and (clock.hour != 2 or clock.minute >= 15):
        raise SystemExit('Compiles and runs stop at 02:15 until 08:00 local')
    if ROOT not in Path(_core.__file__).resolve().parents:
        raise RuntimeError('Engine extension is outside this checkout')
    for name, expected in json.loads((STUDY/'inputs/input-hashes.json').read_text()).items():
        if hashlib.sha256((STUDY/'inputs'/name).read_bytes()).hexdigest() != expected:
            raise RuntimeError('Changed frozen input: '+name)
    scratch = ROOT/'build/issue-916'
    scratch.mkdir(parents=True, exist_ok=True)
    build = Path((ROOT/'build/environment-wheel/meson-build-dir').read_text().strip())
    commands = json.loads((build/'compile_commands.json').read_text())
    command = next(c for c in commands if c['file'].endswith('/tests/native/surface_example.cpp'))
    args = shlex.split(command['command'])
    args[args.index('-o')+1] = str(scratch/'densities.o')
    args[-1] = str(STUDY/'scripts/densities.cpp')
    subprocess.run(['nice','-n','10',*args], cwd=build, check=True)
    link = subprocess.check_output(['ninja','-C',str(build),'-t','commands','surface_example'], text=True).splitlines()[-1]
    args = shlex.split(link)
    args[args.index('-o')+1] = str(scratch/'densities')
    args = [str(scratch/'densities.o') if a.endswith('tests_native_surface_example.cpp.o') else a for a in args]
    subprocess.run(['nice','-n','10',*args], cwd=build, check=True)
    run = subprocess.run(['nice','-n','10',str(scratch/'densities'),str(STUDY/'inputs'),str(STUDY/'results')],
        cwd=build, capture_output=True, text=True, env={**os.environ,'LD_LIBRARY_PATH':str(build)})
    (scratch/'native.stdout').write_text(run.stdout)
    (scratch/'native.stderr').write_text(run.stderr)
    paths = [Path(_core.__file__),build/'libepcsaft_native.so',Path(__file__),STUDY/'scripts/densities.cpp',STUDY/'inputs/input-hashes.json',scratch/'densities']
    identity = dict(base='78d7a5010a6f917a96f76d61a9decb92376fd4ae', timestamp=clock.isoformat(),
        commit=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(), returncode=run.returncode,
        hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
    (STUDY/'results/run-identity.json').write_text(json.dumps(identity,indent=2)+'\n')
    raise SystemExit(run.returncode)
