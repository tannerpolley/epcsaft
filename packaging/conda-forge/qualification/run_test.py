from pathlib import Path
import os,subprocess,sys,shutil,json
r=Path('qualification-runtime').resolve();out=Path(os.environ['IPOPT_PROBE_EVIDENCE']).resolve()/'engine';out.mkdir(parents=True,exist_ok=True)
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',EPCSAFT_HSL_LIBRARY='',EPCSAFT_IPOPT_OPTIONS_FILE='')
os.environ['PATH']=str(Path(sys.prefix)/'Lib/site-packages/epcsaft')+os.pathsep+str(Path(sys.prefix)/'Library/bin')+os.pathsep+os.environ['PATH']
import platform,ctypes
frequency=ctypes.c_longlong();assert ctypes.windll.kernel32.QueryPerformanceFrequency(ctypes.byref(frequency))
(out/'runner.json').write_text(json.dumps(dict(platform=platform.platform(),machine=platform.machine(),processors=os.cpu_count(),qpc_hz=frequency.value,cpu_timer_unit_seconds=1e-7,environment={k:os.environ.get(k) for k in ['ImageOS','ImageVersion','RUNNER_NAME','GITHUB_RUN_ID','GITHUB_SHA','OPENBLAS_NUM_THREADS','OMP_NUM_THREADS']},candidate=(r/'candidate.txt').read_text()),indent=2))
print((out/'runner.json').read_text(),flush=True)
status=json.loads((r/'build-status.json').read_text());print(json.dumps(status),flush=True)
failures=[] if status['returncode']==0 else [('qualification compilation',status['returncode'])]
def run(name,args,cwd=r):
 print('Running',name,flush=True)
 with (out/(name+'.log')).open('w') as f:
  try:code=subprocess.run([str(x) for x in args],cwd=cwd,stdout=f,stderr=subprocess.STDOUT).returncode
  except OSError as error:f.write(str(error));code=-1
 print((out/(name+'.log')).read_text(),flush=True)
 if code:failures.append((name,code))
run('asl-before-python',[sys.executable,r/'asl-load.py'])
run('tolerance',[r/'qualification-tolerance.exe'])
run('dependency',[sys.executable,r/'dependency-test.py'])
run('conda',[sys.executable,r/'conda-test.py'])
run('native-caller',[r/'surface_example.exe'])
run('fit-processes',[r/'fit_windows_checks.exe'])
run('workers',[sys.executable,'-m','pytest','-s','-vv','tests/regression/test_worker_processes.py'])
run('e2-solve',[r/'equilibrium_checks.exe','E2-6'])
run('source-reference',[sys.executable,'-m','pytest','-s','-vv','tests/equilibrium/test_neutral_reference.py::test_source_reference_numerical_error'])
run('asl-python',[sys.executable,r/'python-callers.py',r,out/'python-callers'])
run('pure-independent',[sys.executable,r/'pure-reference.py',out])
run('traits',[r/'qualification-traits.exe'])
run('engine-iterate-bounds',[r/'qualification-iterate-bounds.exe'])
for case in ['F8','F9','M4']:run('properties-'+case,[r/'properties_checks.exe',case])
run('association-independent',[sys.executable,r/'association-compare.py',out])
run('product',[r/'qualification-product.exe'])
for mode in ['qualify','frozen']:
 run('software-'+mode,[r/f'qualification-939-r2-{mode}.exe',r/'near-inputs.txt',out/f'{mode}-events.csv'])
 shutil.copy2(out/f'software-{mode}.log',out/f'software-{mode}.csv')
shutil.copy2(r/'compare-endpoints.py',out/'compare-endpoints.py')
run('endpoints',[sys.executable,out/'compare-endpoints.py'])
for mode in ['qualify','frozen']:
 run('timing-'+mode,[r/f'qualification-939-timing-{mode}.exe',r/'near-inputs.txt',out/f'timing-{mode}-events.csv'])
run('timing-pure-reference',[r/'qualification-939-timing-pure-reference.exe'])
run('source-reference-comparison',[sys.executable,r/'source-reference-compare.py',out])
run('timing-summary',[sys.executable,r/'summarize-timing.py',out])
source=r/'comparison-source'
os.environ['EPCSAFT_SOURCE_ROOT']=str(source)
for command in ['check-equilibrium-corpora','time-key-solves']:
 installed=out/('installed-'+command)
 args=['--fit'] if command=='time-key-solves' else []
 for suffix in ['', '-repeat'] if command == 'check-equilibrium-corpora' else ['']:
  run(command+'-installed'+suffix,[sys.executable,source/'tools'/command,'--collect-only','--out',str(installed)+suffix,*args],cwd=source)
 repeat='-repeat' if command=='check-equilibrium-corpora' else ''
 for label,left,right in [('development-AA',source/'build'/('development-'+command),source/'build'/('development-'+command+repeat)),
                          ('installed-AA',installed,Path(str(installed)+repeat)),
                          ('development-installed',source/'build'/('development-'+command),installed)]:
  run(command+'-'+label,[sys.executable,source/'tools'/command,'--compare-only',left,right,'--out',out/(command+'-'+label)],cwd=source)
(out/'failures.json').write_text(json.dumps(failures,indent=2))
assert not failures,failures
