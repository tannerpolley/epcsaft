"""Transport the accepted S5/S6 inputs to the native density comparison."""
import csv, hashlib, json
from pathlib import Path
import epcsaft as ep
LEAF = Path(__file__).resolve().parents[1]
ROOT = LEAF.parents[2]

def declarations():
    return json.loads((LEAF/'inputs/lithium-declarations.json').read_text())

def main():
    import density
    declaration = declarations()
    source = LEAF/'inputs/source'
    tuples = list(csv.DictReader((source/'rezaee-2026-organic-pcsaft-parameters.csv').open()))[:2]
    assert [r['component'] for r in tuples]==['DES','TOPO']
    records = []
    for row in csv.DictReader((source/'rezaee-2026-si-density-tables.csv').open()):
        r = density.run(declaration['assumptions'],tuples,row)
        records.append(dict(source=row,**r))
        print(row['table_label'],row['T_K'],r['admitted'],r['cause'],flush=True)
    assert len(records)==18
    (LEAF/'results/density-comparison.json').write_text(json.dumps(records,indent=2)+'\n')
    files = [p for p in (LEAF/'inputs').rglob('*') if p.is_file()]
    files += list((LEAF/'scripts').glob('*'))+[ROOT/'build/engine-meson/shared/libepcsaft_native.so',Path(density.__file__)]
    identity = dict(engine=list(ep.provenance()),sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files})
    (LEAF/'results/density-run-identity.json').write_text(json.dumps(identity,indent=2)+'\n')

if __name__=='__main__':
    main()
