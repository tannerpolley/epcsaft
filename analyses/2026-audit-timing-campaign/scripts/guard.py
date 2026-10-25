"""Read-only Engine loading for every campaign interpreter (#578 design, R2).

Before any Engine import: refuse a preloaded Engine, forbid Python process launches, and
replace the editable finder's rebuild with the frozen install-plan mapping, so an import can
never compile. Then check mapped and loaded images against the frozen image record.
"""
import hashlib
import json
import os
from pathlib import Path
import sys

FORBIDDEN = ('subprocess.Popen', 'os.system', 'os.posix_spawn', 'os.spawn', 'os.exec', 'os.fork', 'os.forkpty')


def digest(path):
    with open(path, 'rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def package_digests(tree):
    """SHA-256 of every file the editable mapping serves under epcsaft/, by relative name."""
    out = {}

    def walk(node, prefix):
        for name, value in node.items():
            if name != '__pycache__':
                walk(value, f'{prefix}{name}/') if isinstance(value, dict) else out.setdefault(prefix + name, digest(value))
    walk(tree['epcsaft'], '')
    return dict(sorted(out.items()))


def audit(event, args):
    if event in FORBIDDEN:
        raise RuntimeError(f'campaign guard: {event} refused in an Engine interpreter')


def mapped_images(build_dir):
    """(st_dev, st_ino) of every Engine image this process has mapped, by path (any copy, anywhere)."""
    seen = {}
    for line in Path('/proc/self/maps').read_text().splitlines():
        fields = line.split(maxsplit=5)
        if len(fields) == 6 and ('libepcsaft_' in fields[5] or fields[5].startswith(build_dir) or '/epcsaft/_core.' in fields[5]):
            seen[fields[5]] = fields[5]
    return {p: (os.stat(p).st_dev, os.stat(p).st_ino) for p in seen}


def install(image_path, launcher=None):
    """Guard, load and verify; returns this process's identity record."""
    if any(name == 'epcsaft' or name.startswith('epcsaft.') for name in sys.modules):
        raise RuntimeError('campaign guard: Engine already imported')
    sys.addaudithook(audit)
    image = json.loads(Path(image_path).read_text())
    import _epcsaft_editable_loader as loader
    if digest(loader.__file__) != image['loader_sha256']:
        raise RuntimeError('campaign guard: editable loader differs from the frozen record')
    finders = [f for f in sys.meta_path if type(f).__name__ == 'MesonpyMetaFinder']
    if len(finders) != 1 or finders[0]._build_path != image['build_dir']:
        raise RuntimeError(f'campaign guard: unexpected Engine finders {finders}')
    plan = Path(image['build_dir']) / 'meson-info/intro-install_plan.json'
    if digest(plan) != image['install_plan_sha256']:
        raise RuntimeError('campaign guard: install plan differs from the frozen record')
    tree = loader.collect(json.loads(plan.read_text()))
    finders[0]._rebuild = lambda: tree  # find_spec and _path_hook both read this mapping
    if package_digests(tree) != image['package_sha256']:
        raise RuntimeError('campaign guard: mapped package files differ from the frozen record')
    if launcher:
        sys.executable = str(launcher)  # native fit workers start through the guarded launcher
    import epcsaft
    from epcsaft import _core
    maps = mapped_images(image['build_dir'])
    expected = {str(Path(image['build_dir']) / n) for n in image['libraries_sha256'] if n != 'libepcsaft_asl.so'}
    if set(maps) != expected or Path(_core.__file__).resolve() != Path(image['build_dir'], image['extension']).resolve():
        raise RuntimeError(f'campaign guard: loaded images {sorted(maps)} differ from {sorted(expected)}')
    return dict(pid=os.getpid(), ppid=os.getppid(), interpreter=sys.base_prefix, executable=sys.executable,
                guard_sha256=digest(__file__), image_sha256=digest(image_path), epcsaft=epcsaft.__file__,
                extension=_core.__file__, maps={p: list(v) for p, v in maps.items()},
                affinity=sorted(os.sched_getaffinity(0)), nice=os.getpriority(os.PRIO_PROCESS, 0))


def recheck(image_path):
    """Process-end and block-boundary check: same files, same inodes, same hashes."""
    image = json.loads(Path(image_path).read_text())
    maps = mapped_images(image['build_dir'])
    for path, (dev, ino) in maps.items():
        name = Path(path).name
        if image['libraries_sha256'].get(name) != digest(path) or (os.stat(path).st_dev, os.stat(path).st_ino) != (dev, ino):
            raise RuntimeError(f'campaign guard: {path} replaced during the campaign')
    return {p: list(v) for p, v in maps.items()}
