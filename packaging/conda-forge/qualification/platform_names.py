"""Platform names of the qualification programs and Engine images, and how retained programs find them."""
import os
import sys
import sysconfig
from pathlib import Path

WINDOWS, DARWIN = sys.platform == 'win32', sys.platform == 'darwin'
EXE = '.exe' if WINDOWS else ''
LIBRARY = 'lib{}.dylib' if DARWIN else '{}.dll' if WINDOWS else 'lib{}.so'
NATIVE, ASL = LIBRARY.format('epcsaft_native'), LIBRARY.format('epcsaft_asl')
IPOPT = 'libipopt.1.dylib' if DARWIN else 'epcsaft_ipopt_3119.dll' if WINDOWS else 'libipopt.so.1'
PACKAGE = Path(sysconfig.get_paths()['platlib']) / 'epcsaft'
VENV_PYTHON = 'Scripts/python.exe' if WINDOWS else 'bin/python'
# Windows-only checks and why they are excluded elsewhere.
EXCLUDED = {} if WINDOWS else {
    'dependency': 'HS071 qualified Ipopt before the Windows Engine could build (#936); here Ipopt and MUMPS run in engine-iterate-bounds',
    'fit-processes': 'fit_windows_checks is built only for Windows process handles',
}


def library_environment(directory):
    """Search the Engine images in directory first, then the base environment's libraries.

    Retained programs would otherwise follow their link-time search path into the Meson build tree.
    """
    directory = str(directory)
    if WINDOWS:
        return {'PATH': os.pathsep.join([directory, str(Path(sys.base_prefix) / 'Library/bin'), os.environ['PATH']])}
    if DARWIN:
        return {'DYLD_LIBRARY_PATH': directory, 'DYLD_FALLBACK_LIBRARY_PATH': str(Path(sys.base_prefix) / 'lib')}
    return {'LD_LIBRARY_PATH': os.pathsep.join([directory, str(Path(sys.base_prefix) / 'lib')])}

