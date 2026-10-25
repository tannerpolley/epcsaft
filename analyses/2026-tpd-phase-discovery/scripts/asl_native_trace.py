"""Count the installed ASL driver's work through its own hidden native trace slot.

This is a pinned Linux/glibc diagnostic, not an Engine interface. The TLS offset
must be identified from the inspected library's ELF relocations and Scope code.
Compile asl_native_trace.cpp against the unchanged core trace headers.
"""
import argparse
import csv
import ctypes
import importlib.metadata
import json
from pathlib import Path
import baseline as b
import epcsaft
from pyomo.environ import ExternalFunction


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--wheel', type=Path)
    parser.add_argument('--helper', type=Path, required=True)
    parser.add_argument('--tls-offset', type=int, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--compare', type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise RuntimeError('refusing to overwrite retained trace evidence')
    distribution = importlib.metadata.distribution('epcsaft')
    direct = json.loads(distribution.read_text('direct_url.json'))
    wheel = args.wheel or Path(direct['url'].removeprefix('file://'))
    library = Path(epcsaft.__file__).parent/'libepcsaft_asl.so'
    identity = dict(wheel=str(wheel), wheel_sha256=b.digest(wheel), library_sha256=b.digest(library),
                    helper_sha256=b.digest(args.helper), script_sha256=b.digest(__file__), tls_offset=args.tls_offset,
                    headers={str(p.relative_to(b.ROOT)): b.digest(p) for p in
                             (b.ROOT/'engine/native/core/trace.hpp', b.ROOT/'engine/native/core/convergence.hpp')})
    shim = ctypes.CDLL(str(args.helper.resolve()))
    shim.begin_asl_trace.argtypes = [ctypes.c_char_p, ctypes.c_size_t]
    shim.begin_asl_trace.restype = ctypes.c_void_p
    shim.end_asl_trace.argtypes = [ctypes.c_void_p]
    shim.end_asl_trace.restype = ctypes.c_char_p
    original_load, original_evaluate = b.load, ExternalFunction.evaluate_fgh
    calls = []

    def load(name, path):
        module = original_load(name, path)
        if name == 'khudaida_live':
            module.neutral.installed_identity = lambda: (wheel, identity['wheel_sha256'], library)
        return module

    def evaluate(self, values, *positional, **keywords):
        capture = shim.begin_asl_trace(str(library).encode(), args.tls_offset)
        if not capture:
            raise RuntimeError('native ASL trace slot is unavailable or already occupied')
        try:
            return original_evaluate(self, values, *positional, **keywords)
        finally:
            calls.append(dict(feed=values[:4], phase=values[-2], component=values[-1],
                              **json.loads(shim.end_asl_trace(capture))))

    failure = None
    identical = None
    try:
        b.load, ExternalFunction.evaluate_fgh = load, evaluate
        b.asl()
        prior = list(csv.DictReader((b.HERE/'results/before/asl.csv').open()))
        identical = prior == b.ROWS
        assert identical, 'ASL numerical operands differ from the frozen baseline'
        assert len(calls) == 30
        assert sum(c['driver_finalizations'] for c in calls) > 0, 'capture has no driver control'
        assert all(c['tangent_plane_trials'] == 0 for c in calls), 'implicit stability trial'
    except Exception as error:
        failure = f'{type(error).__name__}: {error}'
        raise
    finally:
        b.load, ExternalFunction.evaluate_fgh = original_load, original_evaluate
        numerical_rows = len(b.ROWS)
        b.ROWS = []
        for i, call in enumerate(calls):
            case = 'asl/'+('anchor', 'start-0.9', 'start-1.1')[i//10]+f'/phase-{call["phase"]}/component-{call["component"]}'
            for quantity, value in call.items():
                if quantity != 'drivers':
                    b.record(case, quantity, value)
            for j, driver in enumerate(call['drivers']):
                for quantity, value in driver.items():
                    b.record(f'{case}/driver-{j}', quantity, value)
        args.output.mkdir(parents=True)
        b.csv_write(args.output/'work.csv', b.ROWS)
        comparison = None
        if args.compare:
            prior_work = list(csv.DictReader((args.compare/'work.csv').open()))
            comparison = dict(identical=prior_work == b.ROWS, rows_before=len(prior_work), rows_after=len(b.ROWS))
            if comparison['identical']:
                comparison.update(operands_csv=str((args.compare/'work.csv').resolve()), operands_sha256=b.digest(args.compare/'work.csv'))
                (args.output/'work.csv').unlink()
        (args.output/'summary.json').write_text(json.dumps(dict(identity=identity, failure=failure, comparison=comparison,
            numerical_operand_rows=numerical_rows, numerical_operands_identical=identical,
            numerical_operands_sha256=b.digest(b.HERE/'results/before/asl.csv'),
            totals={q: sum(c[q] for c in calls) for q in ('driver_finalizations', 'iterations', 'evaluations', 'hessian_calls', 'tangent_plane_trials')},
            limits=['trace-only numerical verification on the inspected Linux/glibc shared library',
                    'TLS offset and unchanged C++ Trace layout are specific to the identified build',
                    'no thermodynamic equations or ASL calculation inputs are changed']), indent=2)+'\n')
        if comparison is not None and not comparison['identical']:
            raise RuntimeError('ASL driver work differs from the frozen baseline')


if __name__ == '__main__':
    main()
