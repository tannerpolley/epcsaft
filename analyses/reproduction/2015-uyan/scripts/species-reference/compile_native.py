import json
from pathlib import Path
import shlex
import subprocess
import sys
root = Path.cwd()
build = Path((root / 'build/environment-wheel/meson-build-dir').read_text().strip())
source = Path(sys.argv[1]).resolve()
output = root / "build/issue-274" / source.stem
entries = json.loads((build / 'compile_commands.json').read_text())
entry = next(e for e in entries if e['file'].endswith('tests/native/properties_checks.cpp'))
args = entry.get('arguments') or shlex.split(entry['command'])
args[args.index('-o')+1] = str(output)+'.o'
args[-1] = str(source)
args.extend(['-I'+str(root / 'engine/tests/native'), *sys.argv[2:]])
subprocess.run(args,cwd=build,check=True)
lines = subprocess.check_output(['ninja','-t','commands','properties_checks'],cwd=build,text=True).splitlines()
args = shlex.split(lines[-1])
args[args.index('-o')+1] = str(output)
args = [str(output)+'.o' if a.endswith('properties_checks.cpp.o') else a for a in args]
subprocess.run(args,cwd=build,check=True)
