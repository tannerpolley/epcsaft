"""Declared native TP coexistence compared with all seventeen May Table 5 rows."""
import csv
import hashlib
import json
import math
from pathlib import Path
import epcsaft as ep
from epcsaft import equilibrium as eq
from epcsaft.parameters import parameter_document

LEAF = Path(__file__).resolve().parents[1]
ROOT = LEAF.parents[2]
SOURCE = LEAF / 'inputs/source'

def declarations():
    document = parameter_document(ROOT / 'data/packets', [('gross-sadowski-2001/2', ['methane', 'ethane'])], k_ij={('methane', 'ethane'): 0})
    for c in document['components']:
        c['fixed']['charge_number'] = {'value': {'magnitude': 0, 'unit': 'elementary-charge'}}
    document['model_families'] = [{'kind': k, 'choice': 'pc-saft' if k == 'base' else 'none'} for k in ('base','association','electrolyte','permittivity')]
    parameters = ep.Parameters.from_mapping(document)
    model = ep.Mixture(parameters)
    feeds = json.loads((SOURCE / 'may-2015-neutral-held-cases.yaml').read_text())['cases']
    cases = [r for r in feeds if r['case_id'].endswith('-midpoint')]
    assert len(cases) == 17 and len({r['source_row_id'] for r in cases}) == 17
    problems = [eq.Problem(phases=[eq.Phase('L', kind='liquid'), eq.Phase('V', kind='vapor')], T=r['T_K'], P=r['P_Pa'],
                           feed=eq.Amounts({'methane': r['feed_z_methane'], 'ethane': 1-r['feed_z_methane']}), row_tolerance=1e-10) for r in cases]
    return model, cases, problems

def main():
    model, cases, problems = declarations()
    source = {r['row_id']: r for r in csv.DictReader((SOURCE / 'may-2015-methane-ethane-vle.csv').open())}
    rows, diagnostics = [], []
    for case, problem in zip(cases, problems):
        r = source[case['source_row_id']]
        result = eq.solve_equilibrium(model, problem)
        fields = 'success message solver_status amounts mole_fractions molar_densities residuals classified_residuals iterations evaluations'.split()
        diagnostic = dict(row_id=r['row_id'], **{k: getattr(result, k) for k in fields})
        row = dict(r, model_x_methane=None, model_y_methane=None, x_error=None, y_error=None, numerical_pass=False,
                   x_outcome='unavailable', y_outcome='unavailable', cause=result.message, evidence_label='physical comparison against measured compositions')
        if len(result.molar_densities) == 2:
            x, y = result.mole_fractions[:2], result.mole_fractions[2:]
            states = [model.state(case['T_K'], rho=rho, x=z) for rho, z in zip(result.molar_densities, (x,y))]
            balance = max(abs(sum(result.amounts[j] * result.mole_fractions[2*j+i] for j in range(2)) - f) for i, f in enumerate((case['feed_z_methane'], 1-case['feed_z_methane'])))
            pressure = max(abs(s.pressure / case['P_Pa'] - 1) for s in states)
            fugacity = max(abs(math.log(x[i]/y[i]) + states[0].log_fugacity_coefficient[i] - states[1].log_fugacity_coefficient[i]) for i in range(2))
            normalized = max(abs(sum(z)-1) for z in (x,y))
            ordered = result.molar_densities[0] > result.molar_densities[1] and all(a > 0 for a in result.amounts)
            numerical = result.success and ordered and normalized <= 1e-10 and balance <= 1e-10 and pressure <= 1e-8 and fugacity <= 1e-8
            diagnostic.update(balance_mol=balance, pressure_relative=pressure, isofugacity=fugacity, normalization=normalized,
                              pressure_slopes=[s.pressure_density_derivative for s in states], balance_limit_mol=1e-10, pressure_limit=1e-8, isofugacity_limit=1e-8)
            row.update(model_x_methane=x[0], model_y_methane=y[0], x_error=abs(x[0]-float(r['x_methane'])), y_error=abs(y[0]-float(r['y_methane'])), numerical_pass=numerical)
            for phase in ('x','y'):
                row[phase+'_outcome'] = 'pass' if numerical and row[phase+'_error'] <= float(r[phase+'_comparison_allowance']) else 'fail' if numerical else 'undetermined'
        rows.append(row)
        diagnostics.append(diagnostic)
        print(r['row_id'], row['numerical_pass'], row['x_outcome'], row['y_outcome'], flush=True)
    assert len(rows) == len(source) == 17
    output = LEAF / 'results'
    with (output / 'comparison.csv').open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0], lineterminator='\n'); writer.writeheader(); writer.writerows(rows)
    (output / 'diagnostics.json').write_text(json.dumps(diagnostics, indent=2) + '\n')
    paths = [*SOURCE.glob('*'), ROOT / 'data/packets/gross-sadowski-2001/2/parameters/parameters.json', Path(ep._core.__file__)]
    identity = dict(engine=list(ep.provenance()), hashes={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths if p.is_file()})
    (output / 'run-identity.json').write_text(json.dumps(identity, indent=2) + '\n')

if __name__ == '__main__':
    main()
