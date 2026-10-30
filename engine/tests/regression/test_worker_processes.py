"""Process isolation must preserve native fit results and remove only the fit's workers."""

import ctypes
import json
import os
from pathlib import Path
import signal
import struct
import subprocess
import sys
import time
from typing import Any

import pytest


def fit(workers: int, count: int = 3, failure: bool = False) -> Any:
    from epcsaft import Mixture, regression as reg
    from tests.native_inputs import ETHANE, parameters

    p = parameters((ETHANE,) * 2, pairs=[("k_ij", (0, 1), 0.0)])
    model = Mixture(p)
    rows = [reg.density(p, t, model.state(t, rho=100.0, x={"s0": 0.5, "s1": 0.5}).pressure,
                       {"s0": 0.5, "s1": 0.5}, phase="vapor", observed=101.0,
                       form="difference", scale=100.0, mass_basis=False) for t in (300. + .04 * i for i in range(count))]
    if failure:
        for row in rows[1:]:
            row.observed, row.scale = 1e300, 1e-300
    coordinate = reg.coordinate(p, [("k_ij", ("s0", "s1"), 1.)], origin=0., scale=.01, bounds=(-.1, .1))
    return reg.fit(p, [coordinate], rows, controls=reg.FitControls(
        maximum_iterations=20, maximum_elapsed_time_seconds=30., workers=workers))


def record(result: Any) -> dict[str, Any]:
    out = {name: getattr(result, name) for name in (
        "optimizer", "physical", "residuals", "optimizer_jacobian", "predictions", "status",
        "residual_evaluations", "jacobian_evaluations", "iterations")}
    out["failures"] = [(f.observation, f.status.name, f.jacobian, f.message) for f in result.failures]
    convergence = result.first_numerical_failure
    out["first_numerical_failure"] = None if convergence is None else {
        key: (struct.pack('d', value).hex() if isinstance(value := getattr(convergence, key), float) else value) for key in ('residual', 'floor', 'amplification', 'contraction', 'iterations', 'evaluations')}
    return out


@pytest.mark.parametrize("failure", [False, True])
def test_worker_numerical_agreement(failure: bool) -> None:
    a, b = record(fit(1, failure=failure)), record(fit(2, failure=failure))
    assert a == b
    if failure:
        assert a["failures"][0][0] == 1
    for name in ("optimizer", "physical", "residuals", "optimizer_jacobian", "predictions"):
        assert struct.pack(f"{len(a[name])}d", *a[name]) == struct.pack(f"{len(b[name])}d", *b[name])


def children(pid: int) -> list[int]:
    if sys.platform != "win32":
        return [int(p.name) for p in Path('/proc').iterdir() if p.name.isdigit()
                and (p / 'stat').exists() and _parent(p) == pid]
    from ctypes import wintypes as w
    class Entry(ctypes.Structure):
        _fields_ = [("size", w.DWORD), ("usage", w.DWORD), ("pid", w.DWORD), ("heap", ctypes.c_size_t),
                    ("module", w.DWORD), ("threads", w.DWORD), ("parent", w.DWORD),
                    ("priority", w.LONG), ("flags", w.DWORD), ("exe", w.WCHAR * 260)]
    kernel = ctypes.WinDLL('kernel32', use_last_error=True)
    kernel.CreateToolhelp32Snapshot.restype = w.HANDLE
    kernel.Process32FirstW.argtypes = [w.HANDLE, ctypes.POINTER(Entry)]
    kernel.Process32NextW.argtypes = [w.HANDLE, ctypes.POINTER(Entry)]
    kernel.CloseHandle.argtypes = [w.HANDLE]
    snapshot = kernel.CreateToolhelp32Snapshot(2, 0)
    assert snapshot != w.HANDLE(-1).value
    entry = Entry()
    entry.size = ctypes.sizeof(entry)
    result = []
    try:
        ok = kernel.Process32FirstW(snapshot, ctypes.byref(entry))
        while ok:
            if entry.parent == pid:
                result.append(entry.pid)
            ok = kernel.Process32NextW(snapshot, ctypes.byref(entry))
    finally:
        kernel.CloseHandle(snapshot)
    return result


def _parent(path: Path) -> int:
    try:
        return int((path / 'stat').read_text().rsplit(')', 1)[1].split()[1])
    except (FileNotFoundError, ProcessLookupError):
        return -1


def _running(pid: int) -> bool:
    try:
        return Path(f'/proc/{pid}/stat').read_text().rsplit(')', 1)[1].split()[0] != 'Z'
    except (FileNotFoundError, ProcessLookupError):
        return False


def _handle_count() -> int:
    if sys.platform != 'win32':
        return len(list(Path('/proc/self/fd').iterdir()))
    from ctypes import wintypes as w
    kernel = ctypes.WinDLL('kernel32', use_last_error=True)
    kernel.GetCurrentProcess.restype = w.HANDLE
    kernel.GetProcessHandleCount.argtypes = [w.HANDLE, ctypes.POINTER(w.DWORD)]
    count = w.DWORD()
    assert kernel.GetProcessHandleCount(kernel.GetCurrentProcess(), ctypes.byref(count))
    return count.value


def _cleanup_host() -> None:
    signal.signal(getattr(signal, 'SIGBREAK', signal.SIGINT), signal.default_int_handler)
    fit(1)  # Warm imports and numerical-library initialization before counting owned handles.
    if sys.platform == 'win32':
        fit(2)  # Initialize Windows process-launch facilities before counting per-fit resources.
    before = _handle_count()
    print('ready', flush=True)
    try:
        while True:
            fit(2, 1000)
    except (KeyboardInterrupt, RuntimeError) as error:
        kind = 'KeyboardInterrupt' if isinstance(error, KeyboardInterrupt) else 'RuntimeError'  # EosError on a crash
        print(json.dumps([kind, before, _handle_count()]), flush=True)
        sys.stdin.readline()  # Observer must establish cleanup while this host stays alive.


@pytest.mark.parametrize('outcome', ['parent-death', 'cancellation', 'worker-crash'])
def test_worker_cleanup(outcome: str) -> None:
    from concurrent.futures import ThreadPoolExecutor
    command = 'from tests.regression.test_worker_processes import _cleanup_host; _cleanup_host()'
    flags = getattr(subprocess, 'CREATE_NEW_PROCESS_GROUP') if sys.platform == 'win32' else 0
    process = subprocess.Popen([sys.executable, '-c', command], stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE, text=True, creationflags=flags,
                               cwd=Path(__file__).resolve().parents[2])
    unrelated = subprocess.Popen([sys.executable, '-c', 'import sys; sys.stdin.read()'], stdin=subprocess.PIPE)
    assert process.stdin and process.stdout
    handles: list[Any] = []
    kernel: Any = None
    readers = ThreadPoolExecutor(1)
    try:
        assert readers.submit(process.stdout.readline).result(timeout=30).strip() == 'ready'
        owned: list[int] = []
        deadline = time.monotonic() + 30
        while not owned and time.monotonic() < deadline and process.poll() is None:
            owned = children(process.pid)
            time.sleep(.001)
        assert owned, 'fit never created a worker'
        if sys.platform == 'win32':
            from ctypes import wintypes as w
            kernel = ctypes.WinDLL('kernel32', use_last_error=True)
            kernel.OpenProcess.argtypes = [w.DWORD, w.BOOL, w.DWORD]
            kernel.OpenProcess.restype = w.HANDLE
            kernel.WaitForSingleObject.argtypes = [w.HANDLE, w.DWORD]
            kernel.TerminateProcess.argtypes = [w.HANDLE, w.UINT]
            kernel.CloseHandle.argtypes = [w.HANDLE]
            for pid in owned:
                handle = kernel.OpenProcess(0x100001, False, pid)
                if handle:
                    handles.append(handle)
                assert handle
        if outcome == 'parent-death':
            process.kill()
            process.wait(timeout=30)
        elif outcome == 'cancellation':
            process.send_signal(getattr(signal, 'CTRL_BREAK_EVENT') if sys.platform == 'win32' else signal.SIGINT)
        elif sys.platform == 'win32':
            assert kernel.TerminateProcess(handles[0], 17)
        else:
            os.kill(owned[0], signal.SIGKILL)
        if outcome != 'parent-death':
            report = readers.submit(process.stdout.readline).result(timeout=30)
            kind, before, after = json.loads(report)
            assert kind == ('KeyboardInterrupt' if outcome == 'cancellation' else 'RuntimeError')
            assert after == before, (before, after)
            assert process.poll() is None, 'host exited before cleanup observation'
        if sys.platform == 'win32':
            assert all(kernel.WaitForSingleObject(h, 30000) == 0 for h in handles)
        else:
            deadline = time.monotonic() + 5
            while time.monotonic() < deadline:
                if all(not _running(pid) for pid in owned):
                    break
                time.sleep(.01)
            else:
                pytest.fail(f'workers survived {outcome}: {owned}')
        assert unrelated.poll() is None
        if outcome != 'parent-death':
            assert process.poll() is None
            process.stdin.write('observed\n')
            process.stdin.flush()
        process.communicate(timeout=30)
    finally:
        if process.poll() is None:
            process.kill()
        try:
            process.communicate(timeout=30)
        finally:
            for handle in handles:
                kernel.CloseHandle(handle)
            try:
                unrelated.communicate(timeout=30)
            finally:
                if unrelated.poll() is None:
                    unrelated.kill()
                    unrelated.wait(timeout=30)
                readers.shutdown(wait=True)


@pytest.mark.skipif(sys.platform != 'win32', reason='Windows image sharing locks')
def test_worker_images_refuse_writes() -> None:
    from epcsaft import _core
    directory = Path(_core.__file__).parent
    for path in (Path(_core.__file__), directory / 'epcsaft_native.dll'):
        with pytest.raises(PermissionError):
            fd = os.open(path, os.O_RDWR)
            os.close(fd)


@pytest.mark.skipif(sys.platform != 'win32', reason='Windows locked-image bootstrap')
@pytest.mark.parametrize('image', [0, 1])
@pytest.mark.parametrize('mismatch', ['path', 'identity'])
def test_worker_image_refusal(tmp_path: Path, image: int, mismatch: str) -> None:
    # Execute the exact production bootstrap; refusal precedes DLL loading.
    import msvcrt
    source = Path(__file__).resolve().parents[2] / 'native/bindings/regression_bindings.cpp'
    bootstrap = source.read_text().split('R"(', 1)[1].split(')"', 1)[0]
    paths = [tmp_path / 'extension.pyd', tmp_path / 'native.dll']
    for path in paths:
        path.write_bytes(b'image identity fixture')
    files = [path.open('rb') for path in paths]
    try:
        handles = [getattr(msvcrt, 'get_osfhandle')(file.fileno()) for file in files]
        for handle in handles:
            getattr(os, 'set_handle_inheritable')(handle, True)
        replacement = tmp_path / 'replacement'
        if mismatch == 'identity':
            replacement.write_bytes(paths[image].read_bytes())
        paths[image] = replacement
        startup = getattr(subprocess, 'STARTUPINFO')()
        startup.lpAttributeList = {'handle_list': handles}
        result = subprocess.run([sys.executable, '-I', '-c', bootstrap, '0', *map(str, handles),
                                 *map(str, paths)], startupinfo=startup, close_fds=True,
                                capture_output=True, text=True, timeout=30)
        assert result.returncode != 0
        assert ('image identity changed' if mismatch == 'identity' else 'FileNotFoundError') in result.stderr
        assert 'CDLL' not in result.stderr, result.stderr
    finally:
        for file in files:
            file.close()


if __name__ == '__main__':
    print(json.dumps([record(fit(n)) for n in (1, 2)], sort_keys=True))
