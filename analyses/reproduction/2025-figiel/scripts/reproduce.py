"""Figure 9 interior source rows through the unchanged native /4 model."""
import csv
import hashlib
import json
import math
from pathlib import Path
import epcsaft as ep

LEAF=Path(__file__).resolve().parents[1]
ROOT=LEAF.parents[2]
SOURCE=LEAF/'figures/figure-09/source'
PACKET=ROOT/'data/packets/figiel-yu-held-2025/4/parameters/parameters.json'
SALTS={'nacl':('sodium-cation','chloride-anion'),'nabr':('sodium-cation','bromide-anion')}
T,P=298.15,100000

def declarations():
    return {p.stem:ep.Mixture(ep.Parameters.from_json(PACKET,components=('water',p.stem.split('-')[1],*SALTS[p.stem.split('-')[2]]))) for p in sorted(SOURCE.glob('water-*.csv'))}

def main():
    models=declarations()
    rows=[]
    for path in sorted(SOURCE.glob('water-*.csv')):
        _,alcohol,salt=path.stem.split('-')
        model=models[path.stem]
        references={}
        for i,r in enumerate(csv.DictReader(path.open()),1):
            xa,m,observed=map(float,(r[f'x_{alcohol.capitalize()}'],r['molality'],r['miac_m']))
            row=dict(system=path.stem,source_file=path.name,source_row=i,source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),T_K=T,P_Pa=P,x_alcohol_salt_free=xa,molality_mol_kg=m,source_coefficient=observed,source_locator=r.get('source','retained Figure 9 source CSV; primary lineage not independently recovered here'),evidence_label='source-limited physical comparison',measurement_uncertainty=None,physical_tolerance=None)
            row.update(included=m>0 and 0<xa<1,value=None,signed_error=None,absolute_error=None,reference_change=None,pressure_relative=None,density_mol_m3=None,pressure_slope=None,reference_tolerance=5e-5,pressure_tolerance=1e-8,numerical_outcome='excluded',cause='pure-solvent endpoint or nonpositive molality')
            if row['included']:
                try:
                    solvent=(1-xa,xa,0,0)
                    if xa not in references:
                        references[xa]=model.reference(T,P,solvent_fractions=solvent,ions=SALTS[salt],stoichiometry=(1,1))
                    ref=references[xa]
                    value=ref.activity(m).mean_ionic_molality_scale
                    row.update(value=value,reference_change=ref.convergence_error)
                    x=ep._core.formula_composition(solvent,model.parameters._native,(2,3),(1,1),m)
                    state=model.state(T,P=P,x=x,phase='liquid')
                    row.update(pressure_relative=abs(state.pressure/P-1),density_mol_m3=state.molar_density,pressure_slope=state.pressure_density_derivative)
                    checks=[value is not None and math.isfinite(value) and value>0,ref.convergence_error is None or (math.isfinite(ref.convergence_error) and ref.convergence_error<=row['reference_tolerance']),math.isfinite(row['pressure_relative']) and row['pressure_relative']<=row['pressure_tolerance'],math.isfinite(state.molar_density) and state.molar_density>0,math.isfinite(state.pressure_density_derivative) and state.pressure_density_derivative>0]
                    row.update(numerical_outcome='pass' if all(checks) else 'unavailable',cause='native reference/state available' if all(checks) else 'native value or convergence/pressure/mechanical criterion not met')
                    if row['numerical_outcome']=='pass':
                        row.update(signed_error=value-observed,absolute_error=abs(value-observed))
                except (RuntimeError,ValueError) as error:
                    row.update(numerical_outcome='unavailable',cause=f'{type(error).__name__}: {error}')
            rows.append(row)
        print(path.stem,sum(r['numerical_outcome']=='pass' for r in rows if r['system']==path.stem),flush=True)
    assert len(rows)==466 and sum(r['included'] for r in rows)==409
    with (LEAF/'results/figure-09-comparison.csv').open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=rows[0],lineterminator='\n');writer.writeheader();writer.writerows(rows)
    library=Path(ep._core.__file__)
    identity=dict(engine=list(ep.provenance()),native_library=str(library),native_library_sha256=hashlib.sha256(library.read_bytes()).hexdigest(),packet=str(PACKET.relative_to(ROOT)),packet_sha256=hashlib.sha256(PACKET.read_bytes()).hexdigest())
    inputs=[*sorted(SOURCE.glob('water-*.csv')),*sorted((LEAF/'inputs/source').glob('*')),LEAF/'inputs/docs/legacy-derived-gates.json']
    identity['source_sha256']={str(p.relative_to(LEAF)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}
    (LEAF/'results/figure-09-run-identity.json').write_text(json.dumps(identity,indent=2)+'\n')

if __name__=='__main__':
    main()
