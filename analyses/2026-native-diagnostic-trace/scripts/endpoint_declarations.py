"""Frozen copy of the #269 MEA zero-inventory declaration and its vapor boundary
(2026-addition-enthalpy-fitting/scripts/{zero_inventory,incipient}.py at 9fdfe7ce)."""
import json
from pathlib import Path
import epcsaft
from epcsaft import equilibrium as eq

INPUT = Path.cwd() / 'engine/data/reference/mea-zero-inventory'
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
neutral = eq.NeutralReference(solvent='water', standard_molality_mol_per_kg=1.0,
    solute_molality_ratios=dict(zip(ids, (.125, .125, 0, .625, .25, .125, .125, .125, .125), strict=True)))
neutral.solute_molality_ratios.pop('water')
problem = eq.Problem(T=313.15, P=1734000,
    feed=eq.Amounts(dict(zip(ids, (0, 1, 7.911062165006594, 0, 0, 0, 0, 0, 0), strict=True))),
    phases=[eq.Phase('L', kind='liquid', support=ids)], reactions=reactions,
    neutral_reference=neutral, row_tolerance=1e-10)
source = {'model': model, 'ids': ids, 'problem': problem}


def mea_boundary(feed=None, support=ids[:3]):
    return eq.Problem(T=problem.T, P=eq.Free(), feed=problem.feed if feed is None else feed,
        phases=[problem.phases[0], eq.Phase('V', kind='vapor', amount=eq.Pinned(0), support=support)],
        reactions=problem.reactions, neutral_reference=problem.neutral_reference,
        standard_gibbs_j_per_mol=problem.standard_gibbs_j_per_mol,
        topology_declared=problem.topology_declared, row_tolerance=problem.row_tolerance)


def declarations():
    return (source['model'], mea_boundary())
