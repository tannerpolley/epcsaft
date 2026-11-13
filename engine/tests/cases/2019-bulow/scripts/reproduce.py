"""Declare the Bülow Eq. (7) checks; read permittivity from native states."""
import csv
import hashlib
import json
import math
from pathlib import Path
import epcsaft as ep

LEAF = Path(__file__).resolve().parents[1]

def declarations():
    source = json.loads((LEAF/'inputs/source/reproduction-cases.json').read_text())
    pure = {r['component']: r for r in csv.DictReader((LEAF/'inputs/parameters/pure/any-solvent.csv').open())}
    d = ep._core.ParameterDocument()
    d.component_ids = ['H2O','C4mim+','NTf2-']
    d.site_component, d.site_multiplicity = [0,0], [1,1]
    d.base_choice, d.association_choice = 'pc-saft', 'general-site'
    d.electrolyte_choice = 'fully-dissociated-debye-huckel'
    d.ionic_dispersion_choice = 'all-ion-pairs'
    d.permittivity_choice, d.polar_choice = 'component-permittivity-mixing', 'none'
    entries = []
    def entry(field, index, value, unit):
        e = ep._core.ParameterEntry()
        e.identity, e.field, e.index, e.magnitude, e.unit = f'{field}/{index}',field,index,float(value),unit
        entries.append(e)
    for i,name in enumerate(d.component_ids):
        r = pure[name]
        for field,key,unit in [('molar_mass','MW','kg / mol'),('segment_count','m','dimensionless'),('sigma_constant','s','angstrom'),('dispersion_energy_over_k','e','kelvin'),('charge_number','z','elementary-charge')]:
            entry(field,[i],r[key],unit)
        entry('relative_permittivity',[i],source['relative_permittivity']['water' if i==0 else 'ionic_liquid'],'dimensionless')
        if i:
            entry('debye_huckel_diameter',[i],r['s'],'angstrom')
            entry('packing_diameter',[i],0.88*float(r['s']),'angstrom')
    for field,key,unit in [('association_energy_over_k','e_assoc','kelvin'),('association_volume','vol_a','dimensionless')]:
        entry(field,[0,1],pure['H2O'][key],unit)
    for i,j in ((0,1),(0,2),(1,2)):
        entry('k_ij',[i,j],0,'dimensionless')
    d.entries = entries
    parameters = ep.Parameters(ep._core.build_parameters(d),tuple(d.component_ids),(('H2O','A'),('H2O','B')))
    return ep.Mixture(parameters), source

def main():
    model,source = declarations()
    checks = source['equation7_checks']
    rows = []
    for i,(x,expected) in enumerate(zip(checks['species_compositions'],checks['expected_relative_permittivity'])):
        value = None
        cause = checks['unavailable_case_indices'].get(str(i),'native fixed-density evaluation; no stability claim')
        if i in checks['available_case_indices']:
            value = model.state(checks['temperature_k'],rho=checks['molar_density_mol_m3'],x=x).bulk_relative_permittivity
        error = None if value is None else abs(value-expected)
        outcome = 'unavailable' if value is None else 'pass' if math.isfinite(value) and error<=checks['absolute_tolerance'] else 'fail'
        rows.append(dict(case=i,T_K=checks['temperature_k'],rho_mol_m3=checks['molar_density_mol_m3'],x_water=x[0],x_cation=x[1],x_anion=x[2],relative_permittivity=value,reference=expected,absolute_error=error,tolerance=checks['absolute_tolerance'],outcome=outcome,cause=cause,evidence_label=source['relative_permittivity']['evidence_label']))
    with (LEAF/'results/comparison.csv').open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=rows[0],lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    library=Path(ep._core.__file__)
    identity=dict(engine=list(ep.provenance()),native_library=str(library),native_library_sha256=hashlib.sha256(library.read_bytes()).hexdigest(),source_sha256={str(p.relative_to(LEAF)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (LEAF/'inputs').rglob('*') if p.is_file()},water_permittivity_override=source['relative_permittivity']['water'],components=list(model.parameters.component_ids))
    (LEAF/'results/run-identity.json').write_text(json.dumps(identity,indent=2)+'\n')
    print([(r['case'],r['relative_permittivity'],r['outcome']) for r in rows])

if __name__=='__main__':
    main()
