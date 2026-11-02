"""Write qualification-recipe.yaml beside recipe.yaml: the shipped recipe plus the qualification.

The shipped recipe owns sources, build and mandatory tests; this adds only the qualification build
script, its development-record host packages and the installed qualification test.
Usage: recipe.py [ARCHIVE REVISION]. A candidate archive replaces the release source; its SHA-256 and revision
become the recorded source identity, and rattler-build refuses an archive that does not match that SHA-256.
"""
import hashlib
import sys
from pathlib import Path

recipe = Path(__file__).resolve().parents[1] / 'recipe.yaml'
text = recipe.read_text()
candidate = []
if len(sys.argv) == 3:
    archive, revision = Path(sys.argv[1]).resolve(), sys.argv[2]
    sha256 = hashlib.sha256(archive.read_bytes()).hexdigest()
    release = next(line for line in text.splitlines() if line.startswith('  release_sha256: '))
    candidate = [(release + '\n', f'  release_sha256: {sha256}\n'),
                 ('  - url: https://github.com/tannerpolley/epcsaft/releases/download/v${{ version }}/epcsaft-${{ version }}.tar.gz\n',
                  f'  - path: {archive}\n'),
                 ('EPCSAFT_SOURCE_REVISION: v${{ version }}\n', f'EPCSAFT_SOURCE_REVISION: {revision}\n')]
elif len(sys.argv) != 1:
    raise SystemExit('usage: recipe.py [ARCHIVE REVISION]')
for old, new in candidate + [
    ('file: ${{ "build.bat" if win else "build.sh" }}',
     'file: ${{ "qualification/build.bat" if win else "qualification/build.sh" }}'),
    ('    - pip\n',
     '    - pip\n' + ''.join(f'    - {name}\n' for name in ('pytest', 'numpy', 'scipy', 'pandas', 'pyomo', 'idaes-pse'))),
    ('  - script:\n      file: run_test.py\n',
     '''  - script:
      file: run_test.py
  - script:
      file: qualification/run_test.py
      env:
        IPOPT_PROBE_EVIDENCE: ${{ env.get("IPOPT_PROBE_EVIDENCE", default="qualification-evidence") }}
    requirements:
      run: [pytest, pyomo, numpy, scipy, pandas, idaes-pse]
    files:
      source:
        - qualification-runtime/
'''),
]:
    if text.count(old) != 1:
        raise SystemExit(f'recipe.yaml no longer has exactly one {old!r}')
    text = text.replace(old, new)
recipe.with_name('qualification-recipe.yaml').write_text(text)
print(recipe.with_name('qualification-recipe.yaml'))
