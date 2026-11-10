import subprocess
import sys

import pytest


@pytest.mark.parametrize('wrapper', ['Capture', 'IpoptIterations'])
@pytest.mark.parametrize('misuse', ['before', 'twice', 'enter'])
def test_capture_transitions(wrapper: str, misuse: str) -> None:
    code = f'''from epcsaft import _core
c = _core.{wrapper}()
if {misuse!r} != "before":
    c.__enter__()
if {misuse!r} == "twice":
    c.__exit__(None, None, None)
try:
    c.__enter__() if {misuse!r} == "enter" else c.__exit__(None, None, None)
except RuntimeError:
    pass
else:
    raise AssertionError("missing RuntimeError")
if {misuse!r} == "enter":
    c.__exit__(None, None, None)
'''
    result = subprocess.run([sys.executable, '-c', code], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
