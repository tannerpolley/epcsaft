"""Run the four declared nonzero-salt subsystem rows natively."""
import csv, hashlib, json
from pathlib import Path
import epcsaft as ep

LEAF = Path(__file__).resolve().parents[1]
ROOT = LEAF.parents[2]
def declarations():
    return json.loads((LEAF / 'inputs/lithium-declarations.json').read_text())

def main():
    import subsystems
    declaration = declarations()
    records = []
    fields = 'success message solver_status amounts mole_fractions molar_densities coordinates residuals classified_residuals iterations evaluations requested_tolerance_met'.split()
    for case in declaration['cases']:
        for row in csv.DictReader((LEAF / f"inputs/source/table_{case['table'].lower()}.csv").open()):
            if float(row['w_licl_feed']) not in case['w_LiCl']:
                continue
            result = subsystems.run(declaration, declaration['parameters']['component_ids'].index(case['solvent']), float(row['w_licl_feed']), float(row['experiment']), float(row['published_calc']))
            solved = result['solved']
            result['solved'] = {key: getattr(solved,key) for key in fields}
            result['rows'] = [(r.kind.name,r.phase,r.index) for r in solved.rows]
            result.update(source=row, table=case['table'], quantity=case['quantity'])
            records.append(result)
            print(case['table'],row['w_licl_feed'],result['admitted'],solved.message,flush=True)
    assert len(records) == 4
    output = LEAF / 'results'
    (output / 'subsystem-comparison.json').write_text(json.dumps(records,indent=2)+'\n')
    files = [p for p in (LEAF/'inputs').rglob('*') if p.is_file()]
    files += list((LEAF/'scripts').glob('*')) + [ROOT/'build/engine-meson/shared/libepcsaft_native.so',Path(subsystems.__file__)]
    identity = dict(engine=list(ep.provenance()),sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files})
    (output/'subsystem-run-identity.json').write_text(json.dumps(identity,indent=2)+'\n')

if __name__ == '__main__':
    main()
