"""Declare the frozen #231 MEA input and collect native pre-fix diagnostics."""
import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import epcsaft
from epcsaft import equilibrium as eq

ROOT = Path(__file__).resolve().parents[3]
INPUT = ROOT / 'engine/data/reference/mea-zero-inventory'
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--declare-only', action='store_true')
parser.add_argument('--output', type=Path, required=True)
parser.add_argument('--native-input', type=Path)
parser.add_argument('--native-run', type=Path, help='compile scripts/zero_inventory.cpp, run it on --native-input and write its log here')
args = parser.parse_args()
mapping = json.loads((INPUT / 'parameters.json').read_text())
for family in mapping['model_families']:
    if family['kind'] == 'electrolyte' and family['choice'] == 'born':
        family.setdefault('c_shell', 1.0)
        family.setdefault('c_dielectric', 1.0)
# The frozen input has no CO2/MEA cross edge; declared to keep that model (#569), not a physical-absence claim.
mapping['topology']['no_cross_association'] = [{
    'component_id_a': 'carbon-dioxide', 'component_id_b': 'monoethanolamine',
    'provenance': {'source_id': 'mea-zero-inventory', 'locator': 'engine/data/reference/mea-zero-inventory/parameters.json: no CO2/MEA edge'}}]
parameters = epcsaft.Parameters.from_mapping(mapping)
if args.native_input:
    doc = parameters._document
    declarations = [str(len(doc.component_ids)), str(len(doc.site_component)), str(len(doc.entries))]
    declarations.extend(json.dumps(name) for name in doc.component_ids)
    for component, multiplicity in zip(doc.site_component, doc.site_multiplicity, strict=True):
        declarations.extend((str(component), str(multiplicity)))
    declarations.extend(json.dumps(getattr(doc, key)) for key in ('base_choice', 'association_choice',
        'electrolyte_choice', 'ionic_dispersion_choice', 'permittivity_choice', 'polar_choice'))
    for entry in doc.entries:
        declarations.extend((json.dumps(entry.identity), json.dumps(entry.field), str(len(entry.index))))
        declarations.extend(str(index) for index in entry.index)
        declarations.extend((repr(entry.magnitude), json.dumps(entry.unit)))
    records = json.loads((INPUT / 'reactions.json').read_text())
    for record in records:
        declarations.extend((json.dumps(record['reaction_id']), *(repr(value) for value in record['stoichiometry'])))
        law, reference = record['engine_correlation'], record['engine_reference']
        declarations.extend(repr(law[key]) for key in ('a', 'b', 'c', 'd', 'reference_temperature', 'temperature_min', 'temperature_max'))
        declarations.extend((json.dumps(law['standard_state_id']), json.dumps(reference['source_basis']),
                             repr(reference['reference_pressure_pa'] or 0)))
    thermal = json.loads((INPUT / 'ideal-calorics.json').read_text())
    declarations.append(repr(thermal['reference_pressure_pa']))
    for name in ('carbon-dioxide', 'monoethanolamine', 'water'):
        record = thermal['components'][name]
        declarations.extend((str(len(record['coefficients'])), *(repr(value) for value in record['coefficients'])))
        declarations.extend(repr(value) for value in record['range_k'])
        declarations.extend((repr(record.get('formation_enthalpy_j_per_mol', 0)),
                             repr(record.get('enthalpy_entropy_anchor_k', 0))))
    args.native_input.parent.mkdir(parents=True, exist_ok=True)
    args.native_input.write_text('\n'.join(declarations) + '\n')
if args.native_run:
    if not args.native_input:
        parser.error('--native-run needs --native-input')
    driver = subprocess.check_output([str(ROOT / 'tools/native-driver'), 'regression_checks',
                                      str(Path(__file__).with_name('zero_inventory.cpp'))], text=True).strip()
    args.native_run.write_text(subprocess.run([str(driver), str(args.native_input)], check=True, text=True,
                                              capture_output=True).stdout)
model = epcsaft.Mixture(parameters)
ids = parameters.component_ids
reactions = []
for record in json.loads((INPUT / 'reactions.json').read_text()):
    law, reference = record['engine_correlation'], record['engine_reference']
    reactions.append(eq.Reaction(dict(zip(ids, record['stoichiometry'], strict=True)),
        name=record['reaction_id'], correlation=eq.ReactionLogPolynomial(
            *(law[key] for key in ('a', 'b', 'c', 'd', 'reference_temperature',
                                   'temperature_min', 'temperature_max', 'standard_state_id'))),
        reference=eq.ReactionReference(getattr(eq.ReferenceSourceBasis, reference['source_basis']),
                                      reference_pressure_pa=reference['reference_pressure_pa'])))
source = eq.NeutralReference(solvent='water', standard_molality_mol_per_kg=1.0,
    solute_molality_ratios=dict(zip(ids, (.125, .125, 0, .625, .25, .125, .125, .125, .125), strict=True)))
source.solute_molality_ratios.pop('water')
problem = eq.Problem(T=313.15, P=1734000,
    feed=eq.Amounts(dict(zip(ids, (0, 1, 7.911062165006594, 0, 0, 0, 0, 0, 0), strict=True))),
    phases=[eq.Phase('L', kind='liquid', support=ids)], reactions=reactions,
    neutral_reference=source, row_tolerance=1e-10)


def declarations():
    return model, problem


if '--declare-only' in sys.argv:
    raise SystemExit(0)
out = {'inputs': {path.name: hashlib.sha256(path.read_bytes()).hexdigest()
                  for path in INPUT.glob('*.json')}, 'module': epcsaft.__file__}
try:
    compiled = eq.compile_problem(model, problem)
    result = eq.solve_equilibrium(model, problem)
    out['initial'] = {'success': result.success, 'removed_species': compiled.removed_species,
                      'rows': result.residuals, 'message': result.message}
except Exception as error:
    out['initial'] = {'cause': str(error)}
try:
    feed = model.state(313.15, x=(1, 0, 0, 0, 0, 0, 0, 0, 0), P=1734000, phase='vapor')
    out['feed'] = {key: getattr(feed, key) for key in
                   ('molar_density', 'pressure_density_derivative', 'born', 'debye_huckel', 'residual_enthalpy')}
except Exception as error:
    out['feed'] = {'cause': str(error)}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(out, indent=2) + '\n')
print(json.dumps(out, indent=2))
