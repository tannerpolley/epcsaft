"""Transport the frozen #912 declarations to the native E3 driver."""
import csv
import json
from pathlib import Path
import epcsaft

ROOT = Path(__file__).resolve().parents[3]
SOURCE = ROOT / 'analyses/2026-reaction-reference-association/inputs'
OUT = ROOT / 'analyses/2026-reaction-reference-association/inputs'
parameters = epcsaft.Parameters.from_mapping(json.loads((SOURCE / 'run-document.json').read_text()))
doc = parameters._document
def declarations():
    return [parameters]

def write_document(out, doc):
    out.write(f'{len(doc.component_ids)} {len(doc.site_component)} {len(doc.entries)}\n')
    out.write(' '.join(map(json.dumps, doc.component_ids)) + '\n')
    for comp, mult in zip(doc.site_component, doc.site_multiplicity):
        out.write(f'{comp} {mult}\n')
    out.write(' '.join(json.dumps(getattr(doc, name)) for name in ('base_choice', 'association_choice', 'electrolyte_choice', 'ionic_dispersion_choice', 'permittivity_choice', 'polar_choice')) + '\n')
    for e in doc.entries:
        out.write(f'{json.dumps(e.identity)} {json.dumps(e.field)} {len(e.index)} ' + ' '.join(map(str, e.index)) + f' {e.magnitude:.17g} {json.dumps(e.unit)}\n')

def main():
    with (OUT / 'mdea.txt').open('w') as out:
        write_document(out, doc)
        for row in csv.DictReader((ROOT / 'analyses/reproduction/2013-pahlavanzadeh-mdea/inputs/source/kuranov-1996-table-2.csv').open()):
            if int(row['source_row']) in (11, 18, 50, 51):
                state = json.loads((SOURCE / 'retained-roots.json').read_text())[row['source_row']]
                out.write(' '.join([row['source_row'], row['temperature_K'], row['m_MDEA_mol_per_kg_H2O'], row['m_CO2_mol_per_kg_H2O'], row['ten_times_total_pressure_MPa'], str(state['pressure']), *map(str, state['x']), *map(str, state['y'])]) + '\n')

    import importlib.util
    spec = importlib.util.spec_from_file_location("ionic_compare", ROOT / "analyses/2026-mea-parameter-estimation/scripts/ionic_compare.py")
    source = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(source)
    document = source.packet_document()
    params = epcsaft.Parameters.from_mapping(document)
    selected = [r for r in csv.DictReader((source.HERE / "solvent-fugacity.csv").open()) if r["model"] == "original-born" and float(r["carbamate_scale"]) == 1]
    canonical = {r["observation_id"]: r for r in csv.DictReader(source.VLE.open())}
    with (OUT / "mea.txt").open("w") as out:
        write_document(out, params._document)
        for index, r in enumerate(selected):
            t, loading = float(r["temperature_K"]), float(r["loading"])
            carbamate, bicarbonate = source.interpolate(source.speciation(t), loading)
            amounts = source.amounts(loading, carbamate, bicarbonate, source.molar_masses())
            pressure = float(canonical[r["observation_id"]]["total_pressure"]) * 1000
            expected = [next(float(v["model_fugacity_kPa"])*1000 for v in selected if v["observation_id"]==r["observation_id"] and v["component"]==c) if c in ("water", "monoethanolamine") else -1 for c in params.component_ids]
            out.write(" ".join(map(str, [index,t,loading*0.3/source.molar_masses()["monoethanolamine"],carbamate+bicarbonate+2*source.TRACE,pressure,loading,*[amounts[c] for c in params.component_ids],*expected]))+"\n")
    (OUT / "mea-rows.json").write_text(json.dumps(selected, indent=2)+"\n")

if __name__ == "__main__":
    main()


def e5_inputs():
    import copy
    frozen = json.loads((SOURCE / 'run-document.json').read_text())
    packet = json.loads((ROOT / 'data/packets/nasrifar-2010-mdea/1/parameters/parameters.json').read_text().replace('MDEA', 'methyldiethanolamine'))
    fits = json.loads((SOURCE.parent / 'results/e4.json').read_text())['fits']
    rows = json.loads((ROOT / 'build/issue-922/e5/baseline-rows.json').read_text())
    neutral = [r for r in csv.DictReader((ROOT / 'analyses/reproduction/2010-nasrifar-mdea/inputs/source/figure-targets.csv').open()) if r['series'] == 'measured' and r['panel'] in ('6a', '6b')]
    for name in ('A0', 'A6a', 'A6b', 'B6a', 'B6b'):
        candidate = copy.deepcopy(frozen)
        pair = candidate['pairs'][0]['coefficients'][0]
        if name != 'A0':
            fit = next(f for f in fits if f['training_panel'] == name[1:] and f['start_offset'] == 0 and f['candidate'].startswith('nasrifar' if name[0] == 'B' else 'inherited'))
            if name[0] == 'A':
                pair['temperature_law']['terms'][0]['coefficient']['value']['magnitude'] = fit['physical'][0]
            else:
                pair.pop('temperature_law')
                pair['value'] = {'magnitude': fit['physical'][0], 'unit': 'dimensionless'}
                candidate['pairs'][0]['provenance'] = f"E4 results/e4.json; candidate={fit['candidate']}; training_panel={fit['training_panel']}; start_offset={fit['start_offset']}; physical[0]"
                candidate['correlations'] = [c for c in candidate['correlations'] if c['family'] != 'segment_diameter']
                for comp in candidate['components'][:2]:
                    comp['provenance'] = 'Nasrifar 2010 Table 2: segment and self-association coefficients; retained Bülow Tables 1–2: molar mass, charge and solvation factor'
                    comp['coefficients'] = [c for c in comp['coefficients'] if c['family'] not in ('segment_count', 'segment_diameter', 'dispersion_energy_over_k')] + copy.deepcopy(next(c['coefficients'] for c in packet['components'] if c['component_id'] == comp['component_id']))
                for site in candidate['topology']['sites']:
                    if site['component_id'] == 'methyldiethanolamine' and site['site_role'] == 'acceptor':
                        site['multiplicity'] = 1
                for edge in candidate['topology']['edges']:
                    if edge['source']['kind'] == 'explicit' and edge['endpoint_a']['component_id'] in ('water', 'methyldiethanolamine'):
                        source = next(e for e in packet['topology']['edges'] if e['source']['kind'] == 'explicit' and e['endpoint_a']['component_id'] == edge['endpoint_a']['component_id'])
                        edge['source'] = copy.deepcopy(source['source'])
                        for field in ('energy_over_k', 'volume'):
                            edge[field]['value'] = copy.deepcopy(source[field]['value'])
        (OUT / f'e5-{name}.json').write_text(json.dumps(candidate, indent=2) + '\n')
        for mode in ('reactive', 'neutral'):
            params = epcsaft.Parameters.from_mapping(candidate, **({'components': ('methyldiethanolamine', 'water')} if mode == 'neutral' else {}))
            with (OUT / f'e5-{name}-{mode}.txt').open('w') as out:
                write_document(out, params._document)
                if mode == 'neutral':
                    for i, r in enumerate(neutral):
                        out.write(f"{i} {r['y']} {r['x']} 0 {r['pressure_Pa']} 0 0 0 0 0\n")
                else:
                    for r in rows:
                        out.write(' '.join(map(str, [r['source_row'], r['T_K'], r['m_MDEA'], r['m_CO2'], r['observed_P_Pa']/1e5, r['predicted_P_Pa'], *r['x'], *r['y']])) + '\n')
