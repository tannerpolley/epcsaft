import ctypes,os,sys
from platform_names import ASL,PACKAGE,WINDOWS,library_environment
# Deliberately load ASL before importing epcsaft.
if WINDOWS:dirs=[os.add_dll_directory(d) for d in library_environment(PACKAGE)['PATH'].split(os.pathsep)[:2]]
f=PACKAGE/ASL;dll=ctypes.CDLL(str(f));print('ASL loaded before package import:',f,'handle',dll._handle,flush=True)
assert 'epcsaft' not in sys.modules
