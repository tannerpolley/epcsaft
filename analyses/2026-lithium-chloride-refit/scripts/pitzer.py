"""Lithium #188's frozen Pitzer comparison, in its read-only Reaktoro environment."""
import csv
import hashlib
import json
import math
import time
from pathlib import Path

import reaktoro as rkt

HERE = Path(__file__).resolve().parents[1]
D = json.loads((HERE / 'inputs/fit-declaration.json').read_text())
IONS = {'Li': ('Li+', 1), 'Na': ('Na+', 1), 'K': ('K+', 1), 'Mg': ('Mg+2', 2), 'Ca': ('Ca+2', 2), 'Sr': ('Sr+2', 2), 'Cl': ('Cl-', -1)}


def trace_divalent():
    """Song 2010 Tables 1 and 5, pp. 1288/1292: report-only LiCl in MgCl2."""
    assert rkt.__version__ == '2.13.0'
    data = json.loads(rkt.Params.embedded('Pitzer.json').data().dumpJson())
    params = data['ActivityModelParams']['Pitzer']
    binaries = {'Li+': (0.20818, -0.07264, -0.004241), 'Mg+2': (0.35235, 1.6815, 0.00519)}
    mixing = D['F3']['Song_2010_report_only']
    replacements = [(key, [ion, 'Cl-'], value) for ion, values in binaries.items()
                    for key, value in zip(('Beta0', 'Beta1', 'Cphi'), values, strict=True)]
    replacements += [('Theta', ['Li+', 'Mg+2'], mixing['theta_Li_Mg']),
                     ('Psi', ['Li+', 'Mg+2', 'Cl-'], mixing['psi_Li_Mg_Cl'])]
    for key, formulas, value in replacements:
        params[key] = [p for p in params[key] if set(p['Formulas']) != set(formulas)]
        params[key].append({'Formulas': formulas, 'CorrectionModel': 'Phreeqc', 'Parameters': [value]})
    aqueous = rkt.AqueousPhase('H2O Li+ Mg+2 Cl-')
    aqueous.set(rkt.ActivityModelPitzer(rkt.Params.fromStringJSON(json.dumps(data))))
    system = rkt.ChemicalSystem(rkt.PhreeqcDatabase('pitzer.dat'), aqueous)
    rows = []
    for ionic_strength in (1.0, 3.0, 6.0, 8.3):
        state = rkt.ChemicalState(system)
        state.temperature(D['conditions']['temperature_K'], 'K')
        state.pressure(D['conditions']['pressure_Pa'], 'Pa')
        state.set('H2O', 1, 'kg')
        state.set('Li+', 1e-6, 'mol')
        state.set('Mg+2', ionic_strength/3, 'mol')
        state.set('Cl-', 1e-6 + 2*ionic_strength/3, 'mol')
        props = rkt.ChemicalProps(state)
        li, cl = [float(props.speciesActivityCoefficientLg(ion)) for ion in ('Li+', 'Cl-')]
        assert math.isfinite(li) and math.isfinite(cl) and abs(float(props.charge())) <= 1e-12
        rows.append({'I_background_mol_kg': ionic_strength, 'log_gamma_trace': 0.5*(li+cl),
                     'log_gamma_Li': li, 'log_gamma_Cl': cl, 'charge_mol': float(props.charge())})
    with (HERE / 'results/trace-divalent-song.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    (HERE / 'results/trace-divalent-song.json').write_text(json.dumps({
        'reaktoro': rkt.__version__, 'source': 'Song et al. 2010, W3AX2GHK, Tables 1 and 5, pp. 1288 and 1292',
        'binary_parameters': binaries, 'theta_Li_Mg': mixing['theta_Li_Mg'], 'psi_Li_Mg_Cl': mixing['psi_Li_Mg_Cl'],
        'meaning': 'Native Pitzer with Song coefficients, unequal-charge mixing and native water properties; not a measurement',
        'temperature_K': D['conditions']['temperature_K'], 'pressure_Pa': D['conditions']['pressure_Pa'],
        'trace_molality_mol_kg': 1e-6, 'declaration_sha256': hashlib.sha256((HERE / 'inputs/fit-declaration.json').read_bytes()).hexdigest()
    }, indent=2) + '\n')


def main():
    began = time.perf_counter()
    assert rkt.__version__ == '2.13.0'
    aqueous = rkt.AqueousPhase(rkt.speciate('H O Li Na K Mg Ca Sr Cl'))
    aqueous.set(rkt.ActivityModelPitzer())
    db = rkt.PhreeqcDatabase('pitzer.dat')
    system = rkt.ChemicalSystem(db, aqueous)
    solver = rkt.EquilibriumSolver(system)
    values = []
    for row in csv.DictReader((HERE / D['lithium']['path']).open()):
        mol = {e: float(row[e])/1000 for e in IONS if e != 'Cl'}
        mol['Cl'] = math.fsum(n*IONS[e][1] for e,n in mol.items())
        state = rkt.ChemicalState(system)
        state.temperature(D['conditions']['temperature_K'], 'K')
        state.pressure(D['conditions']['pressure_Pa'], 'Pa')
        state.set('H2O', 1, 'kg')
        for e,n in mol.items():
            if n:
                state.set(IONS[e][0], n, 'mol')
        initial = list(rkt.ChemicalProps(state).elementAmounts())
        result = solver.solve(state)
        assert result.succeeded(), result.optima.failure_reason
        props = rkt.ChemicalProps(state)
        closure = max(abs(float(a)-float(b)) for a,b in zip(initial, props.elementAmounts(), strict=True))
        assert closure <= 1e-8 and abs(float(props.charge())) <= 1e-8
        values.append({'state': row['state'], 'log10_activity': {e: float(props.speciesActivityLg(s)) for e,(s,_) in IONS.items()},
                       'element_closure_mol': closure, 'charge_mol': float(props.charge()), 'iterations': result.iterations(),
                       'optimality_error': result.optima.error_optimality, 'feasibility_error': result.optima.error_feasibility})
    (HERE / 'results/pitzer.json').write_text(json.dumps({'reaktoro': rkt.__version__,
        'database_sha256': hashlib.sha256(rkt.PhreeqcDatabase.contents('pitzer.dat').encode()).hexdigest(),
        'declaration_sha256': hashlib.sha256((HERE / 'inputs/fit-declaration.json').read_bytes()).hexdigest(),
        'wall_seconds': time.perf_counter()-began, 'states': values}, indent=2) + '\n')


if __name__ == '__main__':
    import sys
    if '--trace-divalent' in sys.argv:
        trace_divalent()
    else:
        main()
