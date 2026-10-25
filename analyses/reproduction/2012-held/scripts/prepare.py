"""Freeze source-derived scalar declarations and manually selected figure points; no Engine calls."""
import copy
import csv
import hashlib
import json
import math
from pathlib import Path
import shutil
import sys
import tomllib
from PIL import Image, ImageDraw

LEAF = Path(__file__).resolve().parents[1]
ROOT = LEAF.parents[2]
SOURCE = LEAF / 'inputs/source'
PACKET = ROOT / 'data/packets/held-2012-alcohols/1'
sys.path.insert(0, str(ROOT / 'data/tools'))
from validate_packets import packet_fingerprint

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text((json.dumps(value, indent=2) if path==ROOT/'data/packet-locks.json' else json.dumps(value, separators=(',', ':'))) + '\n')

def rows(name):
    return list(csv.DictReader((SOURCE / name).open()))

IDS = dict(zip(['Li+', 'Na+', 'K+', 'NH4+', 'Cl-', 'Br-', 'I-'], ['lithium-cation', 'sodium-cation', 'potassium-cation', 'ammonium-cation', 'chloride-anion', 'bromide-anion', 'iodide-anion']))
MASS = {r['component']: float(r['molar_mass_kg_mol']) for r in rows('molar-mass-declarations.csv')}
SOL = {(r['component'], r['quantity']): r for r in rows('held-2012-table-5-solvents-and-pairs.csv') if not r['other_component']}
IONS = {(r['ion'], r['solvent']): r for r in rows('held-2012-table-6-ions.csv')}
PROV = dict(source_id='held-2012', locator='Tables 4-6; Eqs. 7-14 and 22; packing: Held 2014 Eq. 21 backward statement; masses: supplementary declarations', domain_id='fixed-temperature')

def scalar(identity, value, unit, family=None):
    return dict(identity=identity, **({'family': family} if family else {}), value=dict(magnitude=value, unit=unit), provenance=PROV)

def document(solvent, ratio=0):
    solvents = ['water', solvent] if 0 < ratio < 1 else [solvent if ratio == 0 else 'water']
    sf = [ratio, 1-ratio] if len(solvents) == 2 else [1]
    eps = (31.929730582040406 if solvent == 'methanol' else 24.7162)
    if ratio: eps += (9.3986*ratio**3+13.622*ratio**2+25.014*ratio if solvent == 'methanol' else 50.738*ratio**3-17.327*ratio**2+21.874*ratio)
    d = dict(schema='epcsaft.parameters', schema_version=1, document_id='held-2012-'+solvent, document_version=1, purpose='external-document', sources=[dict(source_id='held-2012', citation='Held et al. 2012, DOI 10.1016/j.ces.2011.09.040; supplementary masses and packing convention in source declarations')], domains=[dict(domain_id='fixed-temperature', kind='reported-conditions')], components=[], pairs=[], model_coefficients=[], correlations=[], topology=dict(presets=[], sites=[], edges=[]))
    d['model_families'] = [dict(family_id='model/'+kind, kind=kind, choice=choice, provenance=PROV, **({'ionic_dispersion':'exclude-same-sign-ion-pairs'} if kind=='electrolyte' else {})) for kind, choice in [('base','pc-saft'), ('association','general-site'), ('electrolyte','fully-dissociated-debye-huckel'), ('permittivity','solvent-only')]]
    for name in solvents + list(IDS):
        cid = IDS.get(name, name)
        ion = name in IDS
        c = dict(component_id=cid, aliases=[name], fixed={k:dict(value=dict(magnitude=v, unit=u), provenance=PROV) for k,v,u in [('charge_number',int(IONS[name,solvent]['charge_number']) if ion else 0,'elementary-charge'), ('molar_mass',MASS[name],'kilogram / mole')]}, coefficients=[])
        values = {q:(float(SOL[name,q]['value']), SOL[name,q]['unit']) for q in ['segment_count','segment_diameter','dispersion_energy_over_k']} if not ion else {'segment_count':(1,'dimensionless'), 'segment_diameter':(sum(w*float(IONS[name,s]['segment_diameter_angstrom']) for s,w in zip(solvents,sf)),'angstrom'), 'dispersion_energy_over_k':(float(IONS[name,solvents[0]]['dispersion_energy_over_k_K']),'kelvin')}
        if ion: values.update(packing_diameter=(0.88*values['segment_diameter'][0],'angstrom'), debye_huckel_diameter=values['segment_diameter'])
        c['coefficients'] = [scalar('component/'+cid+'/'+q,v,u,q) for q,(v,u) in values.items()]
        d['components'].append(c)
        if not ion:
            d['correlations'].append(dict(correlation_id=cid+'/epsilon', component_id=cid, family='relative_permittivity', form='constant', independent_variables=[], constant=scalar(cid+'/epsilon/value',eps,'dimensionless'), provenance=PROV))
            d['topology']['sites'] += [dict(component_id=cid,site_id=s,site_role=role,multiplicity=1,provenance=PROV) for s,role in [('a','donor'),('b','acceptor')]]
            d['topology']['edges'].append(dict(endpoint_a=dict(component_id=cid,site_id='a'), endpoint_b=dict(component_id=cid,site_id='b'), energy_over_k=scalar(cid+'/energy',float(SOL[name,'association_energy_over_k']['value']),'kelvin'), volume=scalar(cid+'/volume',float(SOL[name,'association_volume']['value']),'dimensionless'), source=dict(kind='explicit',provenance=[PROV])))
    if 'water' in solvents:
        law = copy.deepcopy(json.loads((ROOT/'data/packets/held-cameretti-sadowski-2008/4/parameters/parameters.json').read_text())['correlations'][1]); law['provenance']=PROV
        d['correlations'].append(law)
        d['components'][0]['coefficients'] = [v for v in d['components'][0]['coefficients'] if v['family'] != 'segment_diameter']
    def pair(a,b,values):
        d['pairs'].append(dict(component_id_a=a,component_id_b=b,coefficients=[scalar('pair/'+a+'/'+b+'/'+q,v,'dimensionless',q) for q,v in values.items()]))
    for name in solvents:
        sigma = float(SOL[name,'segment_diameter']['value'])
        if name=='water':sigma=law['constant']['value']['magnitude']+sum(t['amplitude']['value']['magnitude']*math.exp(t['exponent_coefficient']['value']['magnitude']*298.15) for t in law['terms'])
        for ion,cid in IDS.items():
            values={'k_ij':0 if name==solvents[0] else 1-(float(IONS[ion,name]['dispersion_energy_over_k_K'])/float(IONS[ion,solvents[0]]['dispersion_energy_over_k_K']))**0.5}
            if len(solvents)==2:
                mixed=sum(w*float(IONS[ion,s]['segment_diameter_angstrom']) for s,w in zip(solvents,sf))
                values['l_ij']=1-(sigma+float(IONS[ion,name]['segment_diameter_angstrom']))/(sigma+mixed)
            pair(name,cid,values)
    for a in list(IDS)[:4]:
        for b in list(IDS)[4:]: pair(IDS[a],IDS[b],{'k_ij':1})
    if len(solvents)==2:
        pair('water',solvent,{'k_ij':-0.085 if solvent=='methanol' else -0.049, **({'l_ij':-0.01} if solvent=='ethanol' else {})})
        for a,b in [('water',solvent),(solvent,'water')]:
            edge=dict(endpoint_a=dict(component_id=a,site_id='a'), endpoint_b=dict(component_id=b,site_id='b'), source=dict(kind='combining-rule',rule_id='wolbach-sandler',inputs=dict(energy_over_k=[a+'/energy',b+'/energy'],volume=[a+'/volume',b+'/volume']),provenance=[PROV]))
            if solvent=='ethanol':edge['volume_correction']=scalar(a+'/'+b+'/correction',0.20,'dimensionless')
            d['topology']['edges'].append(edge)
    return d

for solvent in ['ethanol','methanol']: save(PACKET/'parameters'/f'{solvent}.json',document(solvent))




for original,name,field in [('held-2012-table-4-permittivity.csv','permittivity','system'),('held-2012-table-5-solvents-and-pairs.csv','solvents','component'),('held-2012-table-6-ions.csv','ions','ion'),('molar-mass-declarations.csv','masses','component')]:
    records=rows(original)
    for r in records:r['component_id']=IDS.get(r[field],r[field].split('-')[0])
    with (PACKET/f'parameters/{name}-source.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=records[0],lineterminator='\n');w.writeheader();w.writerows(records)
text = '''schema_version = 2
packet_id = "held-2012-alcohols"
packet_version = "1"
title = "Held 2012 alcohol salt parameters"
description = "Two solvent-specific binary declarations at 298.15 K; no Born term or ion-ion dispersion. Printed screening sigma, packing 0.88 sigma. Methanol permittivity pinned at this temperature. Potassium in ethanol was fitted to ternary activity data."
license_or_use_basis = "Numeric source transcription for scientific verification; no full paper or images redistributed."
tolerance_rationale = "No packet accuracy or numerical acceptance tolerance claimed. Source-fitted parameters are not independent physical validation; source print criteria belong to the analysis. Supplementary masses carry atomic-weight uncertainty."
[[sources]]
source_id = "held-2012"
citation = "Held et al., Measuring and modeling alcohol/salt systems, CES 68 (2012) 328-339, DOI 10.1016/j.ces.2011.09.040"
use_basis = "Tables 4-6 and Eqs. 7-14; Held 2014 Eq. 21 backward packing statement; CIAAW 2024 supplementary masses plus existing solvent packet masses; exact source rows retained."
[[domains]]
domain_id = "fixed-temperature"
description = "298.15 K only; fully dissociated single 1:1 salts in methanol or ethanol. Density and activities; no off-temperature methanol claim."
'''
for path in sorted((PACKET/'parameters').iterdir()):
    text += f'\n[[files]]\npath = "parameters/{path.name}"\nrole = "'+('parameters' if path.suffix=='.json' else 'source')+'"\nmedia_type = "'+('application/json' if path.suffix=='.json' else 'text/csv')+f'"\nsha256 = "{sha(path)}"\nsource_refs = ["held-2012"]\ndomain_refs = ["fixed-temperature"]\nunits_basis = "Explicit SI or angstrom, kelvin and dimensionless units in each document or source column."\n'
(PACKET/'packet.toml').write_text(text)
locks=json.loads((ROOT/'data/packet-locks.json').read_text()); locks['packets']=[p for p in locks['packets'] if p['packet_id']!='held-2012-alcohols']+[dict(packet_id='held-2012-alcohols',packet_version='1',packet_path=str(PACKET.relative_to(ROOT)),packet_fingerprint=packet_fingerprint(PACKET),materialization='byte-for-byte')]
save(ROOT/'data/packet-locks.json',locks)
material=LEAF/'inputs/parameters/held-2012-alcohols'
shutil.copytree(PACKET,material/'1',dirs_exist_ok=True)
save(material/'materialization-v1.json',dict(**locks['packets'][-1],files=[dict(path=str(p.relative_to(PACKET)),sha256=sha(p)) for p in sorted(PACKET.rglob('*')) if p.is_file()]))
series=json.loads((LEAF/'inputs/figures/series.json').read_text()); pixels=json.loads((LEAF/'inputs/figures/pixels.json').read_text()); targets=[]
assert len(series)==29 and set(series)==set(pixels)
images={name:Image.open(ROOT/'build/issue-735'/name).convert('RGB') for name in {s['image'] for s in series.values()}}
for key,s in series.items():
    s.update(image_sha256=sha(ROOT/'build/issue-735'/s['image']),reading_pixels=2,excluded='not digitized: ambiguous ownership at crossings, merged dilute regions, marker overlays and axis contacts')
    if s['figure']!='7':save(LEAF/'inputs/declarations'/f'{key}.json',json.loads((ROOT/'data/packets/held-cameretti-sadowski-2008/4/parameters/parameters.json').read_text()) if key=='3W-gamma' else document(s['solvent'],s['ratio']))
    s['digitization']='hand confirmed isolated curve points' if pixels[key] else 'not digitized: ambiguous ownership'
    for i,(px,py) in enumerate(pixels[key]):
        m=(px-s['x0'])*s['xmax']/(s['x1']-s['x0']); value=s['ymin']+(s['y0']-py)*(s['ymax']-s['ymin'])/(s['y0']-s['y1'])
        targets.append(dict(case_id=key+'-'+str(i),series=key,figure=s['figure'],quantity=s['quantity'],solvent=s['solvent'],ratio=s['ratio'],salt=s['salt'],molality=m,source_value=value,x_interval=2*s['xmax']/(s['x1']-s['x0']),y_interval=2*(s['ymax']-s['ymin'])/(s['y0']-s['y1']),pixel_x=px,pixel_y=py,evidence_label='numerical verification'))
        draw=ImageDraw.Draw(images[s['image']]);draw.ellipse((px-3,py-3,px+3,py+3),outline='red',width=1)
for page,im in images.items():im.save(ROOT/f'build/issue-735/digitization-overlay-{page}.png')
save(LEAF/'inputs/figures/calibration.json',series)
with (LEAF/'inputs/figures/targets.csv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=targets[0],lineterminator='\n');w.writeheader();w.writerows(targets)

checks=json.loads((SOURCE/'source-table-checks.json').read_text())
for name,r in checks['source_tables'].items():r['retained_lf_sha256']=sha(SOURCE/name)
save(SOURCE/'source-table-checks.json',checks)
