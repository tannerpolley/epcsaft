from pathlib import Path
import os,shlex,subprocess
root=Path(__file__).resolve().parents[3]; folder=root/'build/fit-audit-20260929'; folder.mkdir(parents=True,exist_ok=True)
build=Path('build/environment-wheel/meson-build-dir').read_text().strip()
commands=subprocess.check_output(['ninja','-C',build,'-t','commands','regression_checks'],text=True).splitlines()
compile=shlex.split(next(c for c in commands if '-c ../../../engine/tests/native/regression_checks.cpp' in c))
for flag in ('-MQ','-MF'):
    k=compile.index(flag); del compile[k:k+2]
compile.remove('-MD')
# The disposable numerical assertions must execute; Meson production flags include NDEBUG.
compile.remove('-DNDEBUG')
compile[compile.index('-o')+1]=str(folder/'residual-promotion-probe.o')
compile[compile.index('-c')+1]=str(Path(__file__).with_name('residual-promotion-probe.cpp'))
link=shlex.split(commands[-1]); link[link.index('-o')+1]=str(folder/'residual-promotion-probe')
link=[a for a in link if not a.endswith('.o')]
link.insert(link.index('-o')+2,str(folder/'residual-promotion-probe.o'))
(folder/'residual-promotion-commands.txt').write_text(shlex.join(compile)+'\n'+shlex.join(link)+'\n')
subprocess.run(compile,cwd=build,check=True)
subprocess.run(link,cwd=build,check=True)
subprocess.run([str(folder/'residual-promotion-probe'),str(folder)],check=True,
               env={**os.environ,'LD_LIBRARY_PATH':build+os.pathsep+os.environ.get('LD_LIBRARY_PATH','')})
