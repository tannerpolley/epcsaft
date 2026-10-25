import json
from pathlib import Path
from epcsaft import Parameters
root = Path.cwd()
out = root / 'build/issue-781'
doc = json.loads((root / 'data/packets/figiel-2025-chloride-mixture-refit/2/parameters/parameters.json').read_text())
a = Parameters.from_mapping(doc)
alk = json.loads((root / 'data/packets/gross-sadowski-2001/3/parameters/parameters.json').read_text())
c = next(c for c in alk['components'] if c['component_id'] == 'dodecane')
for family, value in [('relative_permittivity', 2.), ('solvation_factor', 1.)]:
    c['coefficients'].append(dict(identity=f'component/dodecane/{family}', family=family, value=dict(magnitude=value, unit='dimensionless')))
for name in a.component_ids:
    doc['pairs'].append(dict(component_id_a=name, component_id_b='dodecane', coefficients=[dict(identity=f'pair/{name}/dodecane/k_ij',family='k_ij',value=dict(magnitude=0.,unit='dimensionless'))]))
doc['components'].append(c)
(out/'packet-with-dodecane.json').write_text(json.dumps(doc,indent=2)+'\n')
b = Parameters.from_mapping(doc)
with (out/'packet-inputs.hpp').open('w') as f:
    f.write('#include "surface/parameter_builder.hpp"\n')
    for label, p in [('aqueous', a),('extended', b)]:
        d = p._document
        f.write(f'epcsaft::Parameters {label}() {{ epcsaft::ParameterDocument d;\n')
        for attr in ['component_ids','site_component','site_multiplicity']:
            f.write(f'd.{attr} = '+json.dumps(list(getattr(d,attr))).replace('[','{').replace(']','}')+';\n')
        for attr in ['base_choice','association_choice','electrolyte_choice','ionic_dispersion_choice','permittivity_choice','polar_choice']:
            f.write(f'd.{attr} = '+json.dumps(getattr(d,attr))+';\n')
        for e in d.entries:
            indices='{'+','.join(map(str,e.index))+'}'
            f.write('d.entries.push_back({'+','.join([json.dumps(e.identity),json.dumps(e.field),indices,repr(e.magnitude),json.dumps(e.unit)])+'});\n')
        f.write('return epcsaft::build_parameters(d); }\n')
print(json.dumps(dict(components=b.component_ids, masses=b.molar_masses, declaration=doc['cation_mixing']),indent=2))
