"""Literal Bülow Figure 1 bubble pressures and separate source-precision probes."""
import csv
import hashlib
import json
import math
from pathlib import Path
import epcsaft as ep
from epcsaft import equilibrium as eq

LEAF = Path(__file__).resolve().parents[1]
ROOT = LEAF.parents[2]

def declarations():
    source = list(csv.DictReader((LEAF / 'inputs/parameters/paper-parameters.csv').open()))
    components = []
    for name, mass in [('sulfolane', 0.12017), ('CO2', 0.0440095)]:
        coefficients = []
        for field, parameter, unit in [('segment_count','m_seg','dimensionless'), ('segment_diameter','sigma','angstrom'), ('dispersion_energy_over_k','u_over_k','kelvin')]:
            value = float(next(r['value'] for r in source if r['component_or_pair'] == name and r['parameter'] == parameter))
            coefficients.append(dict(identity=name+'/'+field, family=field, value=dict(magnitude=value, unit=unit)))
        fixed = {k:dict(value=dict(magnitude=v,unit=u)) for k,v,u in [('molar_mass',mass,'kilogram / mole'),('charge_number',0,'elementary-charge')]}
        components.append(dict(component_id=name, coefficients=coefficients, fixed=fixed))
    sites = [dict(component_id=name, site_id=site, multiplicity=1) for name in ('sulfolane','CO2') for site in ('donor','acceptor')]
    edges = [dict(endpoint_a=dict(component_id=a, site_id='donor'), endpoint_b=dict(component_id=b, site_id='acceptor'),
                  source=dict(kind='explicit'), energy_over_k=dict(identity=a+b+'/energy', value=dict(magnitude=0, unit='kelvin')),
                  volume=dict(identity=a+b+'/volume', value=dict(magnitude=0.045, unit='dimensionless'))) for a in ('sulfolane','CO2') for b in ('sulfolane','CO2')]
    kij = float(next(r['value'] for r in source if r['component_or_pair'] == 'sulfolane-CO2' and r['parameter'] == 'k_ij'))
    document = dict(schema='epcsaft.parameters', schema_version=1, components=components, topology=dict(sites=sites, edges=edges),
                    model_families=[dict(kind=k, choice=v) for k,v in [('base','pc-saft'),('association','general-site'),('electrolyte','none'),('permittivity','none')]],
                    pairs=[dict(component_id_a='sulfolane', component_id_b='CO2', coefficients=[dict(identity='sulfolane-CO2/kij',family='k_ij',value=dict(magnitude=kij,unit='dimensionless'))])])
    model = ep.Mixture(ep.Parameters.from_mapping(document))
    references = list(csv.DictReader((LEAF / 'inputs/source/figure-1/comparison-points.csv').open()))
    return model, references, document

def main():
    model, references, document = declarations()
    rows, probes, diagnostics = [], [], []
    for point, r in enumerate(references):
        T, x, dx = float(r['temperature_K']), float(r['x']), float(r['x_uncertainty'])
        for name, temperature, fraction in [('nominal',T,x), ('T-minus',T-.5,x), ('T-plus',T+.5,x), ('x-minus',T,x-dx), ('x-plus',T,x+dx)]:
            problem = eq.bubble_point(liquid=eq.Phase('L',kind='liquid'), T=temperature, feed=eq.Amounts({'sulfolane':fraction, 'CO2':1-fraction}))
            problem.row_tolerance = 1e-10
            result = eq.solve_equilibrium(model,problem)
            fields = 'success message solver_status pressure amounts mole_fractions molar_densities residuals classified_residuals iterations evaluations'.split()
            diagnostic = dict(point=point, probe=name, **{k:getattr(result,k) for k in fields})
            row = dict(point=point, probe=name, T_K=temperature, x_sulfolane=fraction, native_pressure_bar=None, numerical_pass=False, cause=result.message)
            if result.pressure is not None and len(result.molar_densities)==2:
                compositions = [result.mole_fractions[:2],result.mole_fractions[2:]]
                states = [model.state(temperature, rho=rho, x=z) for rho,z in zip(result.molar_densities,compositions)]
                pressure = max(abs(s.pressure/result.pressure-1) for s in states)
                fugacity = max(abs(math.log(compositions[0][i]/compositions[1][i])+states[0].log_fugacity_coefficient[i]-states[1].log_fugacity_coefficient[i]) for i in range(2))
                normalized = max(abs(sum(z)-1) for z in compositions)
                topology = result.molar_densities[0]!=result.molar_densities[1] and min(result.mole_fractions)>0
                numerical = result.success and topology and normalized<=1e-10 and pressure<=1e-8 and fugacity<=1e-8
                row.update(native_pressure_bar=result.pressure/1e5, numerical_pass=numerical)
                diagnostic.update(pressure_relative=pressure,isofugacity=fugacity,normalization=normalized,pressure_slopes=[s.pressure_density_derivative for s in states])
            if name == 'nominal':
                error = abs(row['native_pressure_bar']-float(r['y'])) if row['native_pressure_bar'] is not None else None
                outcome = 'pass' if row['numerical_pass'] and error<=float(r['pressure_tolerance_bar']) else 'fail' if row['numerical_pass'] else 'undetermined'
                rows.append(dict(row, reference_pressure_bar=float(r['y']), error_bar=error, tolerance_bar=float(r['pressure_tolerance_bar']), outcome=outcome,evidence_label=r['evidence_label'],locator=r['locator']))
            else:
                probes.append(row)
            diagnostics.append(diagnostic)
            print(point,name,row['native_pressure_bar'],row['numerical_pass'],flush=True)
    assert len(rows)==8 and len(probes)==32
    output=LEAF/'results'
    for filename, records in [('comparison.csv',rows),('source-precision-probes.csv',probes)]:
        with (output/filename).open('w') as stream:
            writer=csv.DictWriter(stream,fieldnames=records[0],lineterminator='\n');writer.writeheader();writer.writerows(records)
    (output/'diagnostics.json').write_text(json.dumps(diagnostics,indent=2)+'\n')
    library=Path(ep._core.__file__)
    identity=dict(engine=list(ep.provenance()),native_library=str(library),native_library_sha256=hashlib.sha256(library.read_bytes()).hexdigest(),
                  declared_model=document, numerical_limits=dict(pressure_relative=1e-8,isofugacity=1e-8,normalization=1e-10),
                  source_sha256={str(p.relative_to(LEAF)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (LEAF/'inputs').rglob('*') if p.is_file()})
    (output/'run-identity.json').write_text(json.dumps(identity,indent=2)+'\n')

if __name__=='__main__':
    main()
