"""Figiel /4 declarations for the exact-column and held-coordinate qualification."""
import csv
import json
from pathlib import Path

from epcsaft import Parameters, regression as r

ROOT = Path(__file__).resolve().parents[3]
PACKET = ROOT / 'data/packets/figiel-2025-reference-electrolytes/4/parameters/parameters.json'
CATIONS = ('lithium-cation', 'sodium-cation', 'potassium-cation')
COMPONENTS = ('water', *CATIONS, 'chloride-anion')
TARGETS = (-487020.0, -381930.0, -308890.0)


def parameters(components=COMPONENTS):
    mapping = json.loads(PACKET.read_text())
    # Packet 4 has no water/methanol cross edge; declared to keep that model (#569), not a physical-absence claim.
    mapping['topology']['no_cross_association'] = [{'component_id_a': 'water', 'component_id_b': 'methanol',
        'provenance': {'source_id': mapping['document_id'], 'locator': 'packet 4: no water/methanol edge'}}]
    return Parameters.from_mapping(mapping, components=components)


def coordinates(p):
    return [r.coordinate(p, [('k_ij', (ion, partner), 1)], origin=0.0, scale=1.0, bounds=(-1.0, 1.0))
            for ion in CATIONS for partner in ('water', 'chloride-anion')] + [
        r.coordinate(p, [('ion_suppression_coefficient', target, 1) for target in COMPONENTS[1:]], origin=0.0, scale=1.0, bounds=(3.0, 9.0))]


def row(p, kind, ion='sodium-cation', counter='chloride-anion', m=0.0, x=None, reference=None, observed=1.0, form=None):
    selected = form or ('difference' if 'gibbs' in kind or 'permittivity' in kind else 'log_ratio')
    return r.electrolyte(p, kind, 298.15, 1e5, x or {'water': 1.0}, ions=(ion, counter), stoichiometry=(1, 1),
        molality=m, species=ion if 'gibbs' in kind else None, reference=reference, observed=observed,
        form=selected, scale=1.0 if selected == 'difference' else None)


def holds(p):
    return [r.hydration_hold(row(p, 'solvation_gibbs_energy', ion, observed=target), bracket=(1.0, 8.0),
        energy_tolerance=1e-6, diameter_tolerance=1e-10, maximum_iterations=50) for ion, target in zip(CATIONS, TARGETS)]


def target_rows(p):
    path = ROOT / 'data/packets/hamer-wu-1972-aqueous-alkali-halides/1/observations/hamer-wu-1972-aqueous-alkali-halides.csv'
    data = list(csv.DictReader(path.open()))
    return [row(p, kind, ion, m=float(item['molality_mol_kg'])) for ion, salt in zip(CATIONS, ('LiCl', 'NaCl', 'KCl'))
        for item in data if item['salt'] == salt and float(item['molality_mol_kg']) <= 4.5 for kind in ('mean_ionic_activity_coefficient', 'osmotic_coefficient')]


def row_groups():
    p = parameters()
    yield p, target_rows(p)
    data = list(csv.DictReader((Path(__file__).resolve().parents[2] / '2026-figiel-divalent-chlorides/inputs/bromide-hamer-wu.csv').open()))
    for ion, salt in zip(CATIONS, ('LiBr', 'NaBr', 'KBr')):
        p = parameters(('water', ion, 'bromide-anion'))
        yield p, [row(p, kind, ion, 'bromide-anion', m=float(item['molality_mol_kg'])) for item in data
                  if item['salt'] == salt and float(item['molality_mol_kg']) <= 4.5 for kind in ('mean_ionic_activity_coefficient', 'osmotic_coefficient')]
    base = ROOT / 'analyses/reproduction/2025-figiel/figures'
    for ion, solvent, filename in (('potassium-cation', 'methanol', 'k'), ('sodium-cation', 'ethanol', 'na')):
        p = parameters(('water', solvent, ion, 'chloride-anion'))
        observations = []
        for item in csv.DictReader((base / f'figure-06/source/{filename}.csv').open()):
            fraction = float(item[f'x_{solvent}'])
            observations.append(row(p, 'transfer_gibbs_energy', ion, x={'water': 1.0-fraction, solvent: fraction}, reference={'water': 1.0}))
        yield p, observations
    p = parameters()
    yield p, [row(p, 'solvation_gibbs_energy', ion, observed=target) for ion, target in zip(CATIONS, TARGETS)]
    observations = []
    for item in csv.DictReader((base / 'figure-02/source/dielc-salts-in-water.csv').open()):
        fraction = float(item['x_ion'])
        observations.append(row(p, 'bulk_relative_permittivity', x={'water': 1.0-fraction, 'sodium-cation': fraction/2, 'chloride-anion': fraction/2}))
    yield p, observations
    for solvent in ('methanol', 'ethanol'):
        p = parameters((solvent, 'sodium-cation', 'chloride-anion'))
        observations = []
        for item in csv.DictReader((base / 'figure-03/source/dielc-single-solvent-digitized.csv').open()):
            if item['solvent'] != solvent: continue
            fraction = float(item['x_ion'])
            observations.append(row(p, 'bulk_relative_permittivity', x={solvent: 1.0-fraction, 'sodium-cation': fraction/2, 'chloride-anion': fraction/2}, reference={solvent: 1.0}))
        yield p, observations


def declarations():
    return (coordinates(parameters()), holds(parameters()), list(row_groups()))
