"""Known declaration failures and notebook-number examples for the two analysis checkers."""
import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'build/analysis-declarations/scratch'
OUT.mkdir(parents=True, exist_ok=True)


def module(filename):
    spec = importlib.util.spec_from_file_location(filename, ROOT / 'tools' / filename)
    loaded = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(loaded)
    return loaded


def main():
    declarations = module('check-analysis-declarations.py')
    declarations.main()  # Admitted baseline, with numerical execution blocked.
    cases = [
        ('uyan-model', 'reproduction/2015-uyan/generate_data',
         'fully-dissociated-debye-huckel', 'revised-debye-huckel', 'revised-debye-huckel'),
        ('uyan-unit', 'reproduction/2015-uyan/generate_data',
         '"unit": "1 / kelvin"', '"unit": "1 / kelvin ** 1"', 'kelvin ** 1'),
        ('competitor-coordinate', '2026-competitor-cross-check/engine',
         'r.coordinate(p,[(f,target,1)]', 'r.coordinate(p,f,target', 'positional arguments'),
        ('competitor-observation', '2026-competitor-cross-check/engine',
         "r.observation(p,'pressure'", "r.observation(p,'retired-pressure'", 'retired-pressure'),
    ]
    for label, key, old, new, expected in cases:
        leaf, name = key.rsplit('/', 1)
        original = ROOT / 'analyses' / leaf / 'scripts' / f'{name}.py'
        source = original.read_text()
        assert old in source, f'known failing declaration no longer located: {label}'
        copied = source.replace(old, new)
        # A scratch copy reads the same retained inputs; only the tested declaration changes.
        for owner, value in [('ROOT', ROOT), ('CAMPAIGN', original.parent.parent)]:
            for depth in (3, 4) if owner == 'ROOT' else (1,):
                copied = copied.replace(f'{owner} = Path(__file__).resolve().parents[{depth}]',
                                        f'{owner} = Path({str(value)!r})')
        scratch = OUT / f'{label}.py'
        scratch.write_text(copied)
        try:
            declarations.check(scratch)
        except Exception as error:
            assert expected in str(error), (label, error)
            print(f'falsified {label}: {type(error).__name__}: {error}')
        else:
            raise AssertionError(f'declaration smoke missed {label}')
    try:
        declarations.eq.solve_equilibrium(None, None)
    except RuntimeError as error:
        assert 'numerical execution' in str(error)
        print(f'solve guard: {error}')
    else:
        raise AssertionError('declaration smoke admitted a solve')
    scanner = module('check-notebook-numbers.py')
    for text in ('3,967 comparisons', '7,934.', '1e-8.', '0.008071 error', '1e-8 tolerance', '12345678 rows',
                 r'$2.3\times10^{-9}$ error', '2.1×10⁻⁹ error', r'$10^{-8}$ tolerance',
                 '| Literal table | 7,934 |', '`0.01284` error', '[0.01284](results/data.csv)'):
        assert list(scanner.flags(text)), text
    for text in ('41 request IDs', 'Gross (2001)', '2026-10-05', 'issue #286',
                 'wheel `84dbd5d062ccff3e6b539c4fe1e6498f54002c67823dc8cc15baadf253233bbe`',
                 '`results/2026-example.csv`', '`{python} summary["count"]`',
                 '```{python}\nprint(12345)\n```', '~~~python\nprint(12345)\n~~~',
                 '---\ntitle: 12345\n---\n', '{{< include results/2026-summary.md >}}'):
        assert not list(scanner.flags(text)), text
    print('notebook calibration examples passed; no notebook was changed')


if __name__ == '__main__':
    main()
