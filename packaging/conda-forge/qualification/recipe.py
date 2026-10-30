"""Write qualification-recipe.yaml beside recipe.yaml: the shipped win-64 recipe plus the qualification.

The shipped recipe owns sources, build and mandatory tests; this adds only the qualification build
script, its development-record host packages and the installed qualification test.
"""
from pathlib import Path

recipe = Path(__file__).resolve().parents[1] / 'recipe.yaml'
text = recipe.read_text()
for old, new in [
    ('file: ${{ "build.bat" if win else "build.sh" }}',
     'file: ${{ "qualification/build.bat" if win else "build.sh" }}'),
    ('        - libboost-headers ==1.89.0\n',
     '        - libboost-headers ==1.89.0\n' + ''.join(f'        - {name}\n' for name in ('pytest', 'numpy', 'scipy', 'pandas', 'pyomo', 'idaes-pse'))),
    ('  - script:\n      file: run_test.py\n',
     '''  - script:
      file: run_test.py
  - if: win
    then:
      script:
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
