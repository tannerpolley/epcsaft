"""Prepare retained wheels and detached declaration sources; never compile."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess

VERSIONS = ('1218dc78', '026b3031', 'a36a20a5', '36485a88', '177c2a87',
            '5ad73fb9', 'ffda2bc0', '0bec2bdf', '9853bb6d', 'cf26a17b',
            '21e761c3', 'cc63a135', '1cd3ab15')
ROOT = Path(__file__).resolve().parents[3]
MAIN = Path(subprocess.check_output(['git', 'rev-parse', '--git-common-dir'], text=True).strip()).resolve().parent
BUILD = ROOT / 'build/issue-549'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args()
    BUILD.mkdir(parents=True, exist_ok=True)
    rows, unavailable = [], []
    for sha in VERSIONS:
        wheels = sorted((MAIN / 'build/wheels' / f'main-{sha}').glob('*.whl'))
        if not wheels:
            unavailable.append(dict(version=sha, reason='no retained wheel; preparation never compiles'))
            print('unavailable', sha, flush=True)
            continue
        assert len(wheels) == 1, wheels
        wheel = wheels[0]
        source, env = BUILD / 'declarations' / sha, BUILD / 'envs' / sha
        if not source.exists():
            subprocess.run(['git', 'worktree', 'add', '--detach', str(source), sha], check=True)
        if not env.exists():
            subprocess.run(['uv', 'venv', '--python', '3.13', str(env)], check=True)
        subprocess.run(['uv', 'pip', 'install', '--python', str(env / 'bin/python'),
                        '--no-deps', str(wheel)], check=True)
        commit, date = subprocess.check_output(['git', 'show', '-s', '--format=%H %cI', sha], text=True).split()
        rows.append(dict(version=sha, commit=commit, date=date, source=str(source),
                         python=str(env / 'bin/python'), wheel=str(wheel),
                         wheel_sha256=hashlib.sha256(wheel.read_bytes()).hexdigest()))
        print('prepared', sha, flush=True)
    (BUILD / 'environments.json').write_text(json.dumps(rows, indent=2) + '\n')
    (BUILD / 'unavailable-wheels.json').write_text(json.dumps(unavailable, indent=2) + '\n')


if __name__ == '__main__':
    main()
