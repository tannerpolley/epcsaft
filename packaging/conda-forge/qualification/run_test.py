from pathlib import Path
import os,subprocess,sys,shutil,json
r=Path('qualification-runtime').resolve();out=Path(os.environ['IPOPT_PROBE_EVIDENCE']).resolve()/'engine';out.mkdir(parents=True,exist_ok=True)
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',EPCSAFT_HSL_LIBRARY='',EPCSAFT_IPOPT_OPTIONS_FILE='')
sys.path.insert(0,str(r))
from platform_names import DARWIN,EXCLUDED,EXE,PACKAGE,WINDOWS,library_environment
programs=library_environment(PACKAGE)
# Unix images carry their own search path, so only the retained programs get one (N4); Windows has none.
if WINDOWS:os.environ.update(programs)
import platform,ctypes,time
if WINDOWS:
 frequency=ctypes.c_longlong();assert ctypes.windll.kernel32.QueryPerformanceFrequency(ctypes.byref(frequency))
 clock=dict(qpc_hz=frequency.value,cpu_timer_unit_seconds=1e-7)
else:clock=dict(perf_counter=vars(time.get_clock_info('perf_counter')),process_time=vars(time.get_clock_info('process_time')))
(out/'runner.json').write_text(json.dumps(dict(platform=platform.platform(),machine=platform.machine(),processors=os.cpu_count(),**clock,excluded=EXCLUDED,library_environment=programs,environment={k:os.environ.get(k) for k in ['ImageOS','ImageVersion','RUNNER_NAME','GITHUB_RUN_ID','GITHUB_SHA','OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','PATH','DYLD_LIBRARY_PATH','DYLD_FALLBACK_LIBRARY_PATH','LD_LIBRARY_PATH']},candidate=(r/'candidate.txt').read_text()),indent=2))
print((out/'runner.json').read_text(),flush=True)
status=json.loads((r/'build-status.json').read_text());print(json.dumps(status),flush=True)
failures=[] if status['returncode']==0 else [('qualification compilation',status['returncode'])]
failures += [(v['command']+v['suffix'],v['returncode']) for v in json.loads((r/'development-collections.json').read_text()) if v['returncode']]
def run(name,args,cwd=r,env=None):
 print('Running',name,flush=True)
 if args[0]!=sys.executable:env=programs
 with (out/(name+'.log')).open('w') as f:
  try:code=subprocess.run([str(x) for x in args],cwd=cwd,env=env and {**os.environ,**env},stdout=f,stderr=subprocess.STDOUT).returncode
  except OSError as error:f.write(str(error));code=-1
 print((out/(name+'.log')).read_text(),flush=True)
 if code:failures.append((name,code))
if DARWIN:
 run('images',['sh','-c','for f in "$0"/*.so "$0"/*.dylib; do echo "== $f"; otool -L "$f"; otool -l "$f" | grep -A2 LC_RPATH; codesign -v "$f" || exit 1; done',PACKAGE])
run('asl-before-python',[sys.executable,r/'asl-load.py'])
run('tolerance',[r/('qualification-tolerance'+EXE)])
if 'dependency' not in EXCLUDED:run('dependency',[sys.executable,r/'dependency-test.py'])
run('conda',[sys.executable,r/'conda-test.py'])
run('native-caller',[r/('surface_example'+EXE)])
if 'fit-processes' not in EXCLUDED:run('fit-processes',[r/('fit_windows_checks'+EXE)])
run('workers',[sys.executable,'-m','pytest','-s','-vv','tests/regression/test_worker_processes.py'])
run('e2-solve',[r/('equilibrium_checks'+EXE),'E2-6'])
run('source-reference',[sys.executable,'-m','pytest','-s','-vv','tests/equilibrium/test_neutral_reference.py::test_source_reference_numerical_error'])
run('asl-python',[sys.executable,r/'python-callers.py',r,out/'python-callers'])
run('pure-independent',[sys.executable,r/'pure-reference.py',out],env=programs)
run('traits',[r/('qualification-traits'+EXE)])
run('engine-iterate-bounds',[r/('qualification-iterate-bounds'+EXE)])
for case in ['F8','F9','M4']:run('properties-'+case,[r/('properties_checks'+EXE),case])
run('association-independent',[sys.executable,r/'association-compare.py',out],env=programs)
run('product',[r/('qualification-product'+EXE)])
run('installed-additional',[sys.executable,r/'collect-additional.py',r,out/'installed-additional'])
for file in (out/'installed-additional').iterdir():
 if file.is_file():shutil.copy2(file,out/file.name)
run('source-reference-comparison',[sys.executable,r/'source-reference-compare.py',out])
run('timing-summary',[sys.executable,r/'summarize-timing.py',out])
shutil.copytree(r/'comparison-source/build/development-additional',out/'development-additional')
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
