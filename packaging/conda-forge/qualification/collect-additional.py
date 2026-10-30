from pathlib import Path
import json, os, subprocess, sys
b,out=map(lambda x:Path(x).resolve(),sys.argv[1:3]);out.mkdir(parents=True)
q=Path(__file__).resolve().parent
# Executables sit outside either DLL directory, so PATH selects only this interpreter's Engine.
from epcsaft import _core
native=Path(_core.__file__).resolve().parent
os.environ.update(EPCSAFT_QUALIFICATION_IDENTITY='1',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
os.environ['PATH']=str(native)+os.pathsep+os.environ['PATH']
import hashlib
(out/'implementation.json').write_text(json.dumps(dict(native=str(native/'epcsaft_native.dll'),sha256=hashlib.sha256((native/'epcsaft_native.dll').read_bytes()).hexdigest()),indent=2))
import shutil
for driver in ['r2-qualify','r2-frozen','timing-qualify','timing-frozen','timing-pure-reference']:
 exe=out/('qualification-939-'+driver+'.exe');shutil.copy2(b/exe.name,exe)
 name='software-'+driver[3:] if driver.startswith('r2-') else driver
 args=[exe] if driver.endswith('pure-reference') else [exe,q/'near-inputs.txt',out/(driver+'-events.csv')]
 with (out/(name+'.log')).open('w') as log, (out/(name+'-identity.log')).open('w') as identity:
  subprocess.run(list(map(str,args)),stdout=log,stderr=identity,check=True,cwd=out)
 if name.startswith('software-'):shutil.copy2(out/(name+'.log'),out/(name+'.csv'))
subprocess.run([sys.executable,q/'summarize-timing.py',out],check=True,stdout=subprocess.DEVNULL)
shutil.copy2(q/'compare-endpoints.py',out/'compare-endpoints.py')
subprocess.run([sys.executable,out/'compare-endpoints.py'],cwd=out,check=True)
