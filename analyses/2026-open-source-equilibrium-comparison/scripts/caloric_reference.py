"""Declare FeOs inputs and transport native values; C++ owns targets and comparisons."""
import csv, json, sys
from pathlib import Path
from feos import EquationOfState, Identifier, Parameters, PhaseEquilibrium, PureRecord, set_num_threads
from si_units import JOULE, KELVIN, MOL, PASCAL
set_num_threads(1)
root = Path(sys.argv[1])
spec = json.loads((root / 'reference-inputs.json').read_text())
residual = Parameters.from_records([PureRecord(Identifier(name=c['name']), c['mass'], **c['residual']) for c in spec['components']])
ideal = Parameters.from_records([PureRecord(Identifier(name=c['name']), c['mass'], **c['ideal']) for c in spec['components']])
eos = EquationOfState.pcsaft(residual).dippr(ideal)
with (root / (sys.argv[2] + '.csv')).open('w') as f:
    writer = csv.writer(f, lineterminator="\n")
    writer.writerow(['T_K', 'H_J_per_mol', 'S_J_per_mol_K', 'beta', 'xL0', 'xL1', 'xV0', 'xV1'])
    for q in spec[sys.argv[2]]:
        if q['kind'] == 'TP':
            state = PhaseEquilibrium.tp_flash(eos, q['T'] * KELVIN, q['P'] * PASCAL, q['feed'], tol=1e-12, max_iter=1000)
        elif q['kind'] == 'PH':
            state = PhaseEquilibrium.ph_flash(eos, q['P'] * PASCAL, q['target'] * JOULE / MOL, q['feed'], q['T'] * KELVIN, tol=1e-12, max_iter=1000)
        else:
            state = PhaseEquilibrium.ps_flash(eos, q['P'] * PASCAL, q['target'] * JOULE / MOL / KELVIN, q['feed'], q['T'] * KELVIN, tol=1e-12, max_iter=1000)
        writer.writerow([state.liquid.temperature / KELVIN, state.molar_enthalpy / (JOULE / MOL), state.molar_entropy / (JOULE / MOL / KELVIN), state.vapor_phase_fraction, *state.liquid.molefracs.tolist(), *state.vapor.molefracs.tolist()])
