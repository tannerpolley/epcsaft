import ctypes,os,sys
from pathlib import Path
# Deliberately load ASL before importing epcsaft.
p=Path(sys.prefix)/'Lib/site-packages/epcsaft'
dirs=[os.add_dll_directory(str(p)),os.add_dll_directory(str(Path(sys.prefix)/'Library/bin'))]
f=p/'epcsaft_asl.dll';dll=ctypes.CDLL(str(f));print('ASL loaded before package import:',f,'handle',dll._handle,flush=True)
assert 'epcsaft' not in sys.modules
