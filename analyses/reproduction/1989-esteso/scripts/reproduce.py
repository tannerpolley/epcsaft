"""Declare Esteso Table I states; native Engine owns reference and activity solves."""
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
import time
import uuid

LEAF = Path(__file__).resolve().parents[1]
ROOT = LEAF.parents[2]
SOURCE = LEAF / 'inputs/source'
RESULT = LEAF / 'results/native-results.json.gz'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    with path.open() as handle:
        return list(csv.DictReader(handle))


def serial(value):
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(k): serial(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [serial(v) for v in value]
    if type(value).__name__ in ('Stop', 'Category'):
        return value.name
    return serial({k: getattr(value, k) for k in dir(value)
                   if not k.startswith('_') and not callable(getattr(value, k))})


def save(run):
    RESULT.write_bytes(gzip.compress((json.dumps(serial(run), indent=2, allow_nan=False) + '\n').encode('utf-8'), mtime=0))


def declarations():
    from epcsaft import Mixture, Parameters
    from epcsaft.parameters import parameter_document

    ledger = json.loads((SOURCE / 'source-ledger.json').read_text())
    document = parameter_document(ROOT / 'data/packets', [(ledger['packet'], ledger['components'])])
    params = Parameters.from_mapping(document)
    return document, {variant: Mixture(params, solver_options=options)
                      for variant, options in ledger['controls'].items()}


def calculate(rows, solvents, ledger, input_hashes):
    import epcsaft
    from epcsaft import _core, trace

    started = time.perf_counter()
    cpu_started = time.process_time()
    document, models = declarations()
    packet = ROOT / 'data/packets' / ledger['packet'] / 'parameters/parameters.json'
    assert digest(packet) == ledger['parameter_sha256']
    identity = dict(run_id=str(uuid.uuid4()), started_utc=datetime.now(timezone.utc).isoformat(),
                    checkout_commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
                    checkout_status=subprocess.check_output(['git', 'status', '--porcelain'], cwd=ROOT, text=True),
                    engine_provenance=serial(epcsaft.provenance()), extension_path=_core.__file__,
                    extension_sha256=digest(Path(_core.__file__)), interpreter=sys.executable,
                    python_version=sys.version, producer_sha256=digest(Path(__file__)), input_hashes=input_hashes,
                    resolved_document_sha256=hashlib.sha256(json.dumps(document, sort_keys=True).encode()).hexdigest(),
                    packet=ledger['packet'], packet_sha256=digest(packet))
    assert Path(sys.prefix).resolve() == ROOT / '.venv', 'Use this checkout’s .venv'
    run = dict(execution_identity=identity, criteria=ledger['criteria'], temperature_K=ledger['temperature_K'],
               pressure_Pa=ledger['pressure_Pa'], components=ledger['components'], stoichiometry=ledger['stoichiometry'],
               controls=ledger['controls'], references=[], activities=[], wall_seconds=None)
    for variant, model in models.items():
        options = run['controls'][variant]
        for solvent in solvents:
            percent = solvent['ethanol_mass_percent']
            ref_id = f"{identity['run_id']}/{variant}/{percent}/reference"
            fractions = [float(solvent['water_mole_fraction']), float(solvent['ethanol_mole_fraction']), 0, 0]
            reference = None
            ref_record = dict(attempt_id=ref_id, run_id=identity['run_id'], variant=variant,
                              ethanol_mass_percent=percent, solvent_fractions=fractions,
                              started_utc=datetime.now(timezone.utc).isoformat(), requested_solver_options=options)
            with trace(diagnostics=True) as capture:
                try:
                    reference = model.reference(run['temperature_K'], run['pressure_Pa'], solvent_fractions=fractions,
                                                ions=run['components'][-2:], stoichiometry=run['stoichiometry'])
                    ref_record.update(molality=reference.molality, convergence_error=reference.convergence_error,
                                      log_fugacity=reference.log_fugacity, solver_options=reference._native.solver_options)
                except Exception as error:
                    ref_record.update(failure_type=type(error).__name__, failure=str(error))
            ref_record['trace'] = serial(capture.trace)
            run['references'].append(ref_record)
            for row in (r for r in rows if r['ethanol_mass_percent'] == percent):
                record = dict(attempt_id=f"{identity['run_id']}/{variant}/{row['row_id']}", run_id=identity['run_id'],
                              reference_id=ref_id, variant=variant, row_id=row['row_id'],
                              ethanol_mass_percent=percent, molality_mol_kg=float(row['molality_mol_kg']),
                              started_utc=datetime.now(timezone.utc).isoformat(), value=None)
                with trace(diagnostics=True) as capture:
                    try:
                        if reference is None:
                            raise RuntimeError('Reference unavailable: ' + ref_record['failure'])
                        activity = reference.activity(record['molality_mol_kg'])
                        value = activity.mean_ionic_molality_scale
                        if not math.isfinite(value):
                            raise ArithmeticError('Nonfinite mean ionic coefficient')
                        record.update(value=value, log_value=activity.mean_ionic_log_molality_scale,
                                      solver_options=activity.solver_options)
                    except Exception as error:
                        record.update(failure_type=type(error).__name__, failure=str(error))
                record['trace'] = serial(capture.trace)
                run['activities'].append(record)
            save(run)
    run.update(wall_seconds=time.perf_counter() - started, cpu_seconds=time.process_time() - cpu_started)
    save(run)
    return run


def compare(rows, run):
    refs = {r['attempt_id']: r for r in run['references']}
    records = {(r['row_id'], r['variant']): r for r in run['activities']}
    compared = []
    for row in rows:
        base, fine = (records[row['row_id'], v] for v in ('baseline', 'refined'))
        value, refined = base['value'], fine['value']
        delta = abs(value - refined) if value is not None and refined is not None else None
        estimates = [abs(r['value'] * math.expm1(refs[r['reference_id']]['convergence_error']))
                     for r in (base, fine) if r['value'] is not None and refs[r['reference_id']]['convergence_error'] is not None]
        analytic = delta is not None and all(refs[r['reference_id']]['molality'] == 0.0
                       and refs[r['reference_id']]['convergence_error'] is None for r in (base, fine))
        allocation = run['criteria']['numerical_absolute']
        numerical = (delta is not None and (analytic or len(estimates) == 2)
                     and max([delta, *estimates]) <= allocation)
        error = value - float(row['gamma_pm_m']) if value is not None else None
        radius = allocation if analytic else delta + max(estimates) if delta is not None and len(estimates) == 2 else None
        outcome = 'unavailable'
        if error is not None:
            outcome = 'reproduced' if abs(error) <= run['criteria']['printed_absolute'] else 'not reproduced'
            if not numerical or abs(abs(error) - run['criteria']['printed_absolute']) <= radius:
                outcome = 'unresolved'
        density = [s['convergence'] for s in base['trace']['spans'] if s['name'] == 'solve_density']
        qualified = (not base['trace']['incomplete'] and len(density) == 1 and density[0] is not None
                     and density[0]['stop'] == 'converged'
                     and density[0]['residual'] <= density[0]['criterion']['value'])
        compared.append(dict(row, predicted_gamma=value, refined_gamma=refined, signed_error=error,
                             absolute_error=abs(error) if error is not None else None,
                             relative_error=error / float(row['gamma_pm_m']) if error is not None else None,
                             printed_absolute_tolerance=run['criteria']['printed_absolute'],
                             numerical_absolute_allocation=run['criteria']['numerical_absolute'],
                             refinement_difference=delta, reference_scaled_change=max(estimates) if estimates else None,
                             numerical_resolution='resolved' if numerical else 'incomplete',
                             activity_density_trace_qualified=qualified, outcome=outcome,
                             cause='undetermined' if outcome != 'reproduced' else '',
                             failure=base.get('failure', ''), evidence='physical comparison',
                             baseline_attempt_id=base['attempt_id'], refined_attempt_id=fine['attempt_id'],
                             execution_run_id=base['run_id'], checkout_commit=run['execution_identity']['checkout_commit']))
    with (LEAF / 'results/comparisons.csv').open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=compared[0], lineterminator='\n')
        writer.writeheader()
        writer.writerows(compared)
    shared = next(r for r in compared if r['ethanol_mass_percent'] == '60' and float(r['molality_mol_kg']) == 0.5)
    reference = json.loads((LEAF / 'results/shared-row-reference.json').read_text())['quantity']
    assert shared['predicted_gamma'] == reference['value'], 'Diagnose shared-row mismatch before interpretation'
    print(json.dumps(dict(run_id=run['execution_identity']['run_id'], rows=len(compared),
                         shared_row=shared['predicted_gamma'], wall_seconds=run['wall_seconds'],
                         outcomes={o: sum(r['outcome'] == o for r in compared) for o in sorted({r['outcome'] for r in compared})})))


def main():
    ledger = json.loads((SOURCE / 'source-ledger.json').read_text())
    verification = json.loads((SOURCE / 'source-verification.json').read_text())
    for name, expected in verification['tables'].items():
        assert digest(SOURCE / name) == expected['sha256']
        assert len(read(SOURCE / name)) == expected['rows']
    rows = read(SOURCE / 'esteso-1989-table-I.csv')
    keys = {(r['ethanol_mass_percent'], float(r['molality_mol_kg'])) for r in rows}
    assert len(rows) == len(keys) == verification['checks']['table_I_unique_numeric_composition_molality_pairs']
    legacy = read(SOURCE / ledger['historical']['legacy_input'])
    indexed = {(r['ethanol_mass_percent'], float(r['molality_mol_kg'])): float(r['gamma_pm_m']) for r in rows}
    assert len(legacy) == verification['checks']['legacy_rows_matching_new_table_I']
    assert all(indexed[r['ethanol_mass_percent'], float(r['molality_mol_kg'])] == float(r['gamma_pm_m']) for r in legacy)
    for row in rows:
        row['row_id'] = f"{row['ethanol_mass_percent']}/{row['series_row']}"
    solvents = read(LEAF / 'inputs/solvent-compositions.csv')
    for solvent in solvents:
        assert sum(r['ethanol_mass_percent'] == solvent['ethanol_mass_percent'] for r in rows) == verification['checks']['table_I_counts'][solvent['ethanol_mass_percent']]
        assert abs(float(solvent['water_mole_fraction']) + float(solvent['ethanol_mole_fraction']) - 1) <= ledger['criteria']['solvent_input_absolute']
    hashes = {str(p.relative_to(LEAF)): digest(p) for p in sorted((LEAF / 'inputs').rglob('*')) if p.is_file()}
    if '--compare-only' in sys.argv:
        run = json.loads(gzip.decompress(RESULT.read_bytes()))
        assert hashes == run['execution_identity']['input_hashes'], 'Retained execution input changed'
    else:
        run = calculate(rows, solvents, ledger, hashes)
    compare(rows, run)


if __name__ == '__main__':
    main()
