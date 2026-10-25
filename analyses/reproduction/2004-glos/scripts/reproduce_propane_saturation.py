"""Declare fixed-parameter propane saturation observations; retain each native result."""
import argparse
import csv
import hashlib
import json
import math
import subprocess
import sys
import time
import uuid
from datetime import datetime, timezone
from decimal import Decimal
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
ROOT = HERE.parents[2]
SOURCE = HERE / 'inputs/source'
OUT = HERE / 'results'
PACKET = ROOT / 'data/packets/gross-2001-propane/3/parameters/parameters.json'
NIST_PACKET = ROOT / 'data/packets/nist-srd69-pure-saturation-properties/1'
BASELINE, REFINED, CHANGE = 1e-10, 1e-11, 1e-7
NIST_LABEL = 'supplementary reference-equation comparison; different equations of state, not independent experimental validation'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    with path.open(newline='') as stream:
        return list(csv.DictReader(stream))


def write(path, rows):
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def declarations():
    source = read(SOURCE / 'glos-2004-propane-saturation-source.csv')
    locators = {(r['dataset_number'], r['T_K']): r for r in read(SOURCE / 'glos-2004-propane-source-locators.csv')}
    pressures = {r['T_K']: float(Decimal(r['value']) * 1000) for r in source if r['dataset_number'] == '1'}
    rows = []
    for r in source:
        d, t = r['dataset_number'], r['T_K']
        loc = locators[(d, t)]
        factor = Decimal(1000) if d == '1' else Decimal(1)
        rows.append(dict(source='Glos', dataset=d, T_K=float(t), quantity={'1': 'pressure', '2': 'liquid_density', '3': 'vapor_density'}[d],
                         phase={'1': '', '2': 'L', '3': 'V'}[d], unit='Pa' if d == '1' else 'kg/m^3',
                         observed=float(Decimal(r['value']) * factor), uncertainty=float(Decimal(r['expanded_uncertainty']) * factor),
                         pressure_seed_Pa=pressures.get(t, 1.0), evidence='physical comparison', locator=loc['thermoml_locator'],
                         significant_digits=loc['source_significant_digits'], half_last_place=float(Decimal(loc['source_half_last_place']) * factor)))
    assert len(rows) == len(locators) == 63
    for r in read(SOURCE / 'nist-srd69-propane-reference-eos.csv'):
        for q, col, ucol, phase in [('pressure', 'p_sat_Pa', 'source_pressure_relative_uncertainty', ''),
                                    ('liquid_density', 'rho_sat_liq_kg_m3', 'source_liquid_density_relative_uncertainty', 'L')]:
            value = float(r[col])
            rows.append(dict(source='NIST', dataset='reference_only', T_K=float(r['T_K']), quantity=q, phase=phase,
                             unit='Pa' if not phase else 'kg/m^3', observed=value,
                             uncertainty=value * float(r[ucol]) if r[ucol] else None, pressure_seed_Pa=float(r['p_sat_Pa']),
                             evidence=NIST_LABEL, locator=r['row_id'] + '/' + col, significant_digits=8,
                             half_last_place=float(Decimal(10) ** Decimal(Decimal(r[col]).as_tuple().exponent) / 2)))
    return rows


def calculate(rows):
    import epcsaft
    from epcsaft import Parameters, Pure, regression, _core
    from epcsaft.equilibrium import Amounts, Free, Phase, Pinned, Problem
    assert Path(epcsaft.__file__).resolve().is_relative_to(ROOT)
    assert Path(sys.prefix).resolve() == ROOT / '.venv'
    assert sha(PACKET) == '6dc4d60f10ef4364c01cbccf9821447ff9c4c815cc8e94284b853d1cdee71ccb'
    packet = json.loads(PACKET.read_text())
    component = packet['components'][0]
    parameters = Parameters.from_json(PACKET, components=['propane'])
    model = Pure(parameters)
    execution = uuid.uuid4().hex
    started = time.monotonic()
    inputs = {str(p.relative_to(ROOT)): sha(p) for p in sorted(SOURCE.iterdir()) if p.is_file()}
    identity = dict(execution_id=execution, utc=datetime.now(timezone.utc).isoformat(), engine_provenance=epcsaft.provenance(),
                    interpreter=sys.executable, engine_module=epcsaft.__file__, native_module=_core.__file__, native_sha256=sha(Path(_core.__file__)),
                    commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
                    runner_sha256=sha(Path(__file__)), inputs=inputs, parameter_path=str(PACKET.relative_to(ROOT)), parameter_sha256=sha(PACKET),
                    nist_packet_files={str(p.relative_to(ROOT)): sha(p) for p in sorted(NIST_PACKET.rglob('*')) if p.is_file()},
                    parameters=[dict(parameter=c['family'], **c['value']) for c in component['coefficients']] +
                               [dict(parameter='molar_mass', **component['fixed']['molar_mass']['value'])],
                    parameter_domain=packet['domains'][0], model={m['kind']: m['choice'] for m in packet['model_families']},
                    baseline_tolerance=BASELINE, refined_tolerance=REFINED, refinement_relative_limit=CHANGE)

    def predict(row, tolerance, stage):
        row = dict(row, execution_id=execution, call_id=f'{execution}/{stage}/{row["source"]}/{row["dataset"]}/{row["T_K"]}/{row["quantity"]}',
                   stage=stage, requested_tolerance=tolerance)
        problem = Problem(phases=[Phase('L', kind='liquid', composition_guess=[1.0]),
                                  Phase('V', kind='vapor', amount=Pinned(0.0), composition_guess=[1.0])],
                          T=row['T_K'], P=Free(row['pressure_seed_Pa']), feed=Amounts({'propane': 1.0}), row_tolerance=tolerance)
        observation = regression.observation(parameters, 'pressure' if row['quantity'] == 'pressure' else 'phase_density', problem,
                                             observed=row['observed'], form='log_ratio', mass_basis=bool(row['phase']), **({'phase': row['phase']} if row['phase'] else {}))
        result = regression.predict(model, observation)
        c = result.convergence
        row.update(predicted=result.value, status=result.status.name, message=result.message, declared_residual=result.residual,
                   requested_tolerance_met=result.requested_tolerance_met, closure=c.residual, criterion=c.criterion.value if c.criterion else None,
                   stop=c.stop.name, floor=c.floor, iterations=c.iterations, evaluations=c.evaluations,
                   branch=str(result.branch), stable_root_count=result.stable_root_count, solver_options=str(result.solver_options))
        row['qualified'] = bool(result.status == regression.ObservationStatus.Available and result.value is not None and math.isfinite(result.value)
                                and result.value > 0 and result.requested_tolerance_met and c.criterion and c.criterion.value <= tolerance and c.residual <= tolerance)
        print(row['call_id'], row['status'], row['qualified'], flush=True)
        return row

    baseline = [predict(row, BASELINE, 'baseline') for row in rows]
    selected = [(90, 'liquid_density'), (210, 'pressure'), (210, 'liquid_density'), (330, 'pressure'), (330, 'liquid_density'), (330, 'vapor_density')]
    refinement = []
    for row in baseline:
        if row['source'] == 'Glos' and (row['T_K'], row['quantity']) in selected:
            refined = predict({k: row[k] for k in rows[0]}, REFINED, 'refinement')
            refined['baseline_call_id'] = row['call_id']
            refined['baseline_value'] = row['predicted']
            refined['relative_change'] = abs(refined['predicted'] / row['predicted'] - 1) if refined['qualified'] and row['qualified'] else None
            refined['refinement_pass'] = refined['relative_change'] is not None and refined['relative_change'] <= CHANGE
            refinement.append(refined)
    identity['wall_seconds'] = time.monotonic() - started
    OUT.mkdir(exist_ok=True)
    (OUT / 'execution.json').write_text(json.dumps(identity, indent=2) + '\n')
    for source, filename in [('Glos', 'propane-saturation.csv'), ('NIST', 'nist-comparison.csv')]:
        write(OUT / filename, [r for r in baseline if r['source'] == source])
    write(OUT / 'refinement.csv', refinement)


def compare():
    identity = json.loads((OUT / 'execution.json').read_text())
    for p, digest in identity['inputs'].items():
        assert sha(ROOT / p) == digest, p
    rows = read(OUT / 'propane-saturation.csv') + read(OUT / 'nist-comparison.csv')
    assert len(rows) == 91 and all(r['execution_id'] == identity['execution_id'] for r in rows)
    for row in rows:
        value = float(row['predicted']) if row['predicted'] else math.nan
        error = value - float(row['observed'])
        row.update(signed_error=error, relative_error=error / float(row['observed']),
                   error_over_uncertainty=error / float(row['uncertainty']) if row['uncertainty'] else None,
                   within_uncertainty=abs(error) <= float(row['uncertainty']) if row['uncertainty'] else None)
    for source, filename in [('Glos', 'propane-saturation.csv'), ('NIST', 'nist-comparison.csv')]:
        write(OUT / filename, [r for r in rows if r['source'] == source])
    print('Comparison refreshed with preserved execution identity:', identity['execution_id'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--compare-only', action='store_true')
    args = parser.parse_args()
    if not args.compare_only:
        subprocess.run([sys.executable, str(HERE / 'scripts/check_propane_source_contract.py')], check=True)
        calculate(declarations())
    compare()
