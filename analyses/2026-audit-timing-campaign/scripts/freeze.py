"""Freeze the prepared Engine image of this checkout before any timing (no Engine import, no build).

Writes results/image.json and the guarded native fit-worker launcher. Refuses a dirty Engine tree.
"""
import hashlib
import json
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parents[1]
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE / 'scripts'))
from guard import digest, package_digests  # noqa: E402


def git(*args):
    return subprocess.check_output(['git', *args], cwd=ROOT, text=True).strip()


def main():
    out = HERE / 'results/image.json'
    if out.exists():
        raise FileExistsError(f'{out} is frozen; a new image needs a new campaign')
    dirty = git('status', '--porcelain', '--', 'engine', 'tools', 'pyproject.toml', 'uv.lock')
    if dirty:
        raise RuntimeError(f'Engine tree not clean:\n{dirty}')
    build = ROOT / 'build/engine-meson/shared'
    site = next((ROOT / '.venv/lib').glob('python3.*/site-packages'))
    loader = site / '_epcsaft_editable_loader.py'
    plan = build / 'meson-info/intro-install_plan.json'
    sys.path.insert(0, str(site))
    import _epcsaft_editable_loader as editable  # registers a finder only; no Engine import or build
    tree = editable.collect(json.loads(plan.read_text()))
    extension = next(n for n in tree['epcsaft'] if n.startswith('_core.') and n.endswith('.so'))
    image = dict(
        commit=git('rev-parse', 'HEAD'), engine_diff_sha256=hashlib.sha256(git('diff', 'HEAD', '--', 'engine').encode()).hexdigest(),
        build_dir=str(build), python=str(ROOT / '.venv/bin/python'), loader=str(loader), loader_sha256=digest(loader),
        install_plan_sha256=digest(plan), buildoptions_sha256=digest(build / 'meson-info/intro-buildoptions.json'),
        build_ninja_sha256=digest(build / 'build.ninja'), uv_lock_sha256=digest(ROOT / 'uv.lock'),
        extension=extension, libraries_sha256={n: digest(build / n) for n in (extension, 'libepcsaft_native.so', 'libepcsaft_asl.so')},
        package_sha256=package_digests(tree))
    out.write_text(json.dumps(image, indent=1) + '\n')
    launcher = ROOT / 'build/issue-578/fit-python'
    launcher.parent.mkdir(parents=True, exist_ok=True)
    launcher.write_text(f'#!/bin/sh\nexec {image["python"]} -I {HERE}/scripts/fit_launch.py {out} '
                        f'{HERE}/results/fit-workers.jsonl "$@"\n')
    launcher.chmod(0o755)
    print(json.dumps({k: image[k] for k in ('commit', 'libraries_sha256')}, indent=1))


if __name__ == '__main__':
    main()
