from pathlib import Path
import shutil,subprocess,json
r=Path('qualification-runtime');r.mkdir(exist_ok=True)
b=Path('engine/build-win')
targets=json.loads(subprocess.check_output(['meson','introspect','--targets',str(b)]))
selected=[x for x in targets if x['name'].startswith('qualification-') or x['name'] in ['equilibrium_checks','surface_example','fit_windows_checks','properties_checks']]
with (r/'compile.log').open('w') as log:
 result=subprocess.run(['ninja','-C',str(b),'-j','6','-k','0',*[Path(x['filename'][0]).name for x in selected]],stdout=log,stderr=subprocess.STDOUT)
(r/'build-status.json').write_text(json.dumps(dict(returncode=result.returncode,targets=[x['name'] for x in selected]),indent=2))
print((r/'compile.log').read_text(),flush=True)

for p in b.glob('*.exe'):shutil.copy2(p,r/p.name)
for name in ['inputs.txt','near-inputs.txt','compare-endpoints.py']:
 shutil.copy2(Path('packaging/conda-forge/qualification')/name,r/name)
shutil.copytree('engine/tests',r/'tests',ignore=shutil.ignore_patterns('__pycache__'))
shutil.copytree('engine/native/bindings',r/'native/bindings')
shutil.copy2('packaging/conda-forge/run_test.py',r/'conda-test.py')

shutil.copy2('packaging/conda-forge/qualification/python-callers.py',r/'python-callers.py')

p=r/'tests/regression/test_worker_processes.py'
s=p.read_text().replace('    assert a == b', '    print("serial", json.dumps(a), "parallel", json.dumps(b), flush=True)\n    assert a == b')
p.write_text(s)

shutil.copy2('packaging/conda-forge/qualification/asl-load.py',r/'asl-load.py')

for name in ['binary128-roots.csv','endpoints.txt','pure-reference.py']:
 shutil.copy2(Path('packaging/conda-forge/qualification')/name,r/name)

for name in ['dependency-test.py','linux-reference.json']:
 shutil.copy2(Path('packaging/conda-forge/qualification')/name,r/name)

shutil.copy2('packaging/conda-forge/qualification/summarize-timing.py',r/'summarize-timing.py')

for name in ['source-binary128.csv','source-reference-compare.py']:
 shutil.copy2(Path('packaging/conda-forge/qualification')/name,r/name)

for name in ['association-binary128.csv','association-compare.py']:
 shutil.copy2(Path('packaging/conda-forge/qualification')/name,r/name)

(r/'candidate.txt').write_text(__import__('os').environ['EPCSAFT_SOURCE_REVISION']+'\n'+__import__('os').environ['EPCSAFT_SOURCE_SHA256']+'\n')

# Keep test declarations and development numbers outside the installed package.
import os, sys
source = r/'comparison-source'
source.mkdir()
for directory in ['tools','analyses','data','engine']:
 shutil.copytree(directory,source/directory,ignore=shutil.ignore_patterns('results','__pycache__','build','build-win','qualification-939','*.png','*.pdf','*.svg','*.ipynb'))
import runpy
for name in runpy.run_path('tools/check-equilibrium-corpora')['RETAINED_INPUTS']:
 target=source/name
 target.parent.mkdir(parents=True,exist_ok=True)
 shutil.copy2(name,target)
shutil.copy2('uv.lock',source/'uv.lock')
identity=dict(revision=os.environ['EPCSAFT_SOURCE_REVISION'],sha256=os.environ['EPCSAFT_SOURCE_SHA256'])
Path('source-identity.json').write_text(json.dumps(identity))
(source/'source-identity.json').write_text(json.dumps(identity))
subprocess.run([sys.executable,'-m','venv','--system-site-packages','development'],check=True)
dev=Path('development/Scripts/python.exe').resolve()
subprocess.run([str(dev),'-m','pip','install','-e','./engine','--no-deps','--no-build-isolation',
 '--config-settings=build-dir='+str(Path('engine/dev-win').resolve()),
 '--config-settings=setup-args=--wrap-mode=nodownload','--config-settings=compile-args=-j6'],check=True)
collections=[]
for command in ['check-equilibrium-corpora','time-key-solves']:
 args=['--fit'] if command=='time-key-solves' else []
 for suffix in ['', '-repeat'] if command == 'check-equilibrium-corpora' else ['']:
  collected=Path('build/qualification')/('development-'+command+suffix)
  code=subprocess.run([str(dev),'tools/'+command,'--collect-only','--out',str(collected),*args]).returncode
  collections.append(dict(command=command,suffix=suffix,returncode=code))
  if collected.exists():shutil.copytree(collected,source/'build'/collected.name)
(r/'development-collections.json').write_text(json.dumps(collections,indent=2))

subprocess.run(['ninja','-C','engine/dev-win','-j','6', 'qualification-939-r2-qualify.exe','qualification-939-r2-frozen.exe','qualification-939-timing-qualify.exe','qualification-939-timing-frozen.exe','qualification-939-timing-pure-reference.exe'],check=True)
subprocess.run([str(dev),'packaging/conda-forge/qualification/collect-additional.py',str(Path('engine/dev-win').resolve()),str(Path('build/qualification/development-additional').resolve())],check=True)
shutil.copytree('build/qualification/development-additional',source/'build/development-additional')
shutil.copy2('packaging/conda-forge/qualification/collect-additional.py',r/'collect-additional.py')
