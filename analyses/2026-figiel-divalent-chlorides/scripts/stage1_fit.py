"""Declare #209 stage 1; C++ owns Ceres, exact derivatives, holds and checks."""
import copy
import csv
import hashlib
import json
import math
import subprocess
import time
from pathlib import Path

from epcsaft import Parameters, _core, provenance, regression as r

HERE = Path(__file__).resolve().parents[1]
REPO = HERE.parents[1]
FIG = REPO / "analyses/reproduction/2025-figiel/figures"
SOURCE = REPO / "data/packets/figiel-2025-reference-electrolytes/4/parameters/parameters.json"


def declared(path):
    """A Figiel 2025 packet with its historical methanol/water and ethanol/methanol omissions declared (#569): no
    cross edge, not a physical-absence claim. Stored candidates inherit the declarations."""
    doc = json.loads(path.read_text())
    doc['topology']['no_cross_association'] = [
        {'component_id_a': a, 'component_id_b': b, 'provenance': {
            'source_id': doc['document_id'], 'locator': f'{path.relative_to(REPO)}: no {a}/{b} edge'}}
        for a, b in (('methanol', 'water'), ('ethanol', 'methanol'))]
    return doc


BASE = declared(SOURCE)
T, P = 298.15, 100000.0
ALKALI = {"lithium-cation": ("Li", -487020.0), "sodium-cation": ("Na", -381930.0), "potassium-cation": ("K", -308890.0)}
NAMES = [f"k {s}-water" for s, _ in ALKALI.values()] + [f"k {s}-Cl" for s, _ in ALKALI.values()] + ["C"]
IONS = (*ALKALI, "chloride-anion", "bromide-anion")
COMPONENTS = ("water", *IONS)
OUT = HERE / "results/current-main"
OUT.mkdir(parents=True, exist_ok=True)
STEM = "stage1-native"


def write(name, rows):
    with (OUT / f"{STEM}-{name}.csv").open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]),lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def group(name, indices, metric, maximum=None, mean=None, baseline=False, report=False):
    g = r.CheckGroup()
    g.name=name
    g.rows=indices
    g.metric=metric
    g.maximum=maximum
    g.mean=mean
    g.baseline=baseline
    g.report_only=report
    return g


def read(path):
    return list(csv.DictReader(path.open()))


def systems():
    """Every evaluated row: (set, group, components, solvent fractions, kind, rows[(x, quantity, observed)])."""
    S, cl, br = [], "chloride-anion", "bromide-anion"
    hw = read(REPO / "analyses/2026-held-2014-chlorides/inputs/evaluated.csv") + read(HERE / "inputs/kcl-hamer-wu.csv")
    for setname, rows, anion, salts in (("T", hw, cl, ("LiCl", "NaCl", "KCl")), ("G1", read(HERE / "inputs/bromide-hamer-wu.csv"), br, ("LiBr", "NaBr", "KBr"))):
        for salt, cation in zip(salts, ALKALI):
            data = [(float(item["molality_mol_kg"]), q, float(item[f"{q}_evaluated"])) for item in rows if item["salt"] == salt for q in ("gamma", "osmotic")]
            S.append((setname, salt, ("water", cation, anion), (1, 0, 0), "activity", data))
    for ion, solvent, name in (("sodium-cation", "ethanol", "na"), ("potassium-cation", "methanol", "k")):
        data = [(float(item[f"x_{solvent}"]), ion, float(item["g_trans (kJ/mol)"])) for item in read(FIG / f"figure-06/source/{name}.csv")]
        S.append(("G2", f"{ALKALI[ion][0]}+ water-{solvent}", ("water", solvent, ion, cl), None, "transfer", data))
    for figure, solvent, salt in ((7, "methanol", "NaBr"), (8, "methanol", "LiBr"), (8, "ethanol", "LiBr"), (8, "methanol", "NaBr"), (8, "ethanol", "NaBr")):
        cation = {"LiBr": "lithium-cation", "NaBr": "sodium-cation"}[salt]
        data = [(float(item["molality"]), "gamma", float(item["miac_m"])) for item in read(FIG / f"figure-0{figure}/source/{solvent}-{salt.lower()}.csv") if float(item["molality"]) > 0]
        S.append(("G3", f"{salt} {solvent} Figure {figure}", (solvent, cation, br), (1, 0, 0), "activity", data))
    for solvent in ("methanol", "ethanol"):
        for salt, anion in (("NaCl", cl), ("NaBr", br)):
            rows = read(FIG / f"figure-09/source/water-{solvent}-{salt.lower()}.csv")
            for x in sorted({float(item[f"x_{solvent.capitalize()}"]) for item in rows} - {0.0, 1.0}):
                data = [(float(item["molality"]), "gamma", float(item["miac_m"])) for item in rows if float(item[f"x_{solvent.capitalize()}"]) == x and float(item["molality"]) > 0]
                S.append(("G4", f"{salt} water-{solvent}", ("water", solvent, "sodium-cation", anion), (1 - x, x, 0, 0), "activity", data))
    water = [(float(item["x_ion"]), "relative permittivity", float(item[" dielc"])) for item in read(FIG / "figure-02/source/dielc-salts-in-water.csv")]
    S.append(("G5", "water", ("water", "sodium-cation", cl), (1, 0, 0), "permittivity", [item for item in water if item[0] <= 0.14]))
    S.append(("G5 report", "water, x_ion > 0.14", ("water", "sodium-cation", cl), (1, 0, 0), "permittivity", [item for item in water if item[0] > 0.14]))
    for solvent in ("methanol", "ethanol"):
        data = [(float(item["x_ion"]), "permittivity ratio", float(item["epsilon_ratio"])) for item in read(FIG / "figure-03/source/dielc-single-solvent-digitized.csv") if item["solvent"] == solvent]
        S.append(("G5 report", solvent, (solvent, "sodium-cation", cl), (1, 0, 0), "permittivity", data))
    return S



def declarations(p, batch="aqueous", *, declare_only=False):
    cs = [r.coordinate(p, [('k_ij', (ion, partner), 1)], origin=0, scale=1, bounds=(-1, 1))
          for partner in ('water', 'chloride-anion') for ion in ALKALI]
    cs += [r.coordinate(p, [('ion_suppression_coefficient', target, 1) for target in IONS], origin=0, scale=1, bounds=(3, 9))]
    observations, metadata, memberships, checks = [], [], {}, []
    for setname, label, comps, fractions, kind, data in systems():
        solvent=next((v for v in comps if v in ("methanol","ethanol")), "aqueous")
        if batch!="all" and solvent!=batch:
            continue
        for value, quantity, observed in data:
            x = dict(zip(comps, fractions)) if fractions else {'water': 1.0}
            reference, species = None, None
            if kind == 'activity':
                native_kind = 'mean_ionic_activity_coefficient' if quantity == 'gamma' else 'osmotic_coefficient'
                form = 'log_ratio'
                name = f'{setname} {label} {quantity}'
            elif kind == 'transfer':
                native_kind, form, species = 'transfer_gibbs_energy', 'difference', quantity
                solvent = comps[1]
                x = {'water': 1-value, solvent: value}
                reference = {'water': 1.0}
                observed *= 1000.0  # declared input kJ/mol -> native J/mol
                name = f'{setname} {label}'
            else:
                native_kind, form = 'bulk_relative_permittivity', 'difference'
                x = {comps[0]: 1-value, comps[-2]: value/2, comps[-1]: value/2}
                if quantity == 'permittivity ratio':
                    reference = {comps[0]: 1.0}
                name = f'{setname} {label}'
            index = len(observations)
            memberships.setdefault(name, []).append(index)
            ob = r.electrolyte(p, native_kind, T, P, x, ions=tuple(comps[-2:]), stoichiometry=(1, 1),
                molality=value if kind == 'activity' else 0, species=species, reference=reference,
                observed=observed, form=form, scale=1 if form == 'difference' else None)
            observations.append(ob)
            metadata.append(dict(row=index, set=setname, group=name, abscissa=value, quantity=quantity, observed=observed))
    for name, indices in memberships.items():
        log = observations[indices[0]].form == r.ResidualForm.LogRatio
        metric = r.CheckMetric.AbsoluteLog10Ratio if log else r.CheckMetric.AbsoluteDifference
        target, guard = name.startswith('T '), name.startswith(('G1 ', 'G2 ', 'G5 water'))
        baseline = name.startswith(('G1 ', 'G2 ', 'G5 water'))
        maximum = .02 if target or name.startswith('G1 ') else 0 if baseline else None
        checks.append(group(name, indices, metric, maximum, 0 if name.startswith(('G2 ', 'G5 water')) else None,
                            baseline, not (target or guard)))
    hs = []
    for ion, (_, target) in (ALKALI.items() if batch in ("aqueous","all") else []):
        ob = r.electrolyte(p, 'solvation_gibbs_energy', T, P, {'water': 1.0}, ions=(ion, 'chloride-anion'),
            stoichiometry=(1, 1), molality=0, species=ion, observed=target, form='difference', scale=1)
        hs.append(r.hydration_hold(ob, bracket=(1, 8), energy_tolerance=1e-6, diameter_tolerance=1e-10, maximum_iterations=50))
        index=len(observations)
        observations.append(ob)
        metadata.append(dict(row=index, set='E', group=f'E {ion}', abscissa=0, quantity='hydration J/mol', observed=target))
        checks.append(group(f'E {ion}', [index], r.CheckMetric.AbsoluteDifference, 10))
    # Native source predictions supply invariant observations; no hydration hold is re-solved here.
    templates, labels = [], []
    for ion, solvent, filename in (('chloride-anion', 'ethanol', 'cl'), ('bromide-anion', 'methanol', 'br')):
        if batch!="all" and solvent!=batch:
            continue
        for kind, x in [('transfer_gibbs_energy', float(v[f'x_{solvent}']))
                for v in read(FIG / f'figure-06/source/{filename}.csv')]:
            templates.append(r.electrolyte(p, kind, T, P, {'water':1-x, solvent:x},
                ions=('sodium-cation',ion), stoichiometry=(1,1), molality=0, species=ion,
                reference={'water':1.0} if kind.startswith('transfer') else None, observed=0, form='difference', scale=1))
            labels.append(f'I {ion} {kind} {x}')
    if batch in ("aqueous","all"):
        for ion in ("chloride-anion","bromide-anion"):
            templates.append(r.electrolyte(p,"solvation_gibbs_energy",T,P,{"water":1.0},ions=("sodium-cation",ion),
                stoichiometry=(1,1),molality=0,species=ion,observed=0,form="difference",scale=1))
            labels.append(f"I {ion} hydration")
    for solvent in (("water","methanol","ethanol") if batch=="all" else ("water",) if batch=="aqueous" else (batch,)):
        templates += [r.electrolyte(p,'bulk_relative_permittivity',T,P,{solvent:1.0},ions=('sodium-cation','chloride-anion'),
                        stoichiometry=(1,1),molality=0,observed=1,form='log_ratio')]
        labels += [f'I pure {solvent} permittivity']
    invariant = templates if declare_only else _core.observations_at(p._native, [], [], [], templates, [0.0]*len(templates))
    for label, ob in zip(labels, invariant):
        index=len(observations)
        observations.append(ob)
        metric = r.CheckMetric.AbsoluteLog10Ratio if ob.form == r.ResidualForm.LogRatio else r.CheckMetric.AbsoluteDifference
        checks.append(group(label,[index],metric, math.log10(1+1e-12) if ob.form == r.ResidualForm.LogRatio else 1e-6))
        metadata.append(dict(row=index,set='I',group=label,abscissa=0,quantity=ob.kind.name,observed=ob.observed))
    if batch in ('aqueous','all'):
        assert sum(v['set']=='T' for v in metadata)==156
    for v,ob in zip(metadata,observations):
        v['batch']=batch
        v['composition']=';'.join(map(str,ob.composition))
        v['reference_composition']=';'.join(map(str,ob.reference_composition))
        v['components']=';'.join(p.component_ids)
        v['ions']=';'.join(p.component_ids[i] for i in ob.ion_indices)
    return cs, hs, observations, metadata, checks


def stored_document(result):
    doc = copy.deepcopy(BASE)
    replacements = {f'pair/{ion}/water/k_ij':result.chart.physical[j] for j,ion in enumerate(ALKALI)}
    replacements.update({f'pair/chloride-anion/{ion}/k_ij':result.chart.physical[j+3] for j,ion in enumerate(ALKALI)})
    replacements.update({f'component/{ion}/born_diameter':v for ion,v in zip(ALKALI,result.chart.physical[-3:])})
    suppression = result.chart.physical[6:11]
    assert len(suppression)==5 and all(v==suppression[0] for v in suppression)
    replacements['model/relative_permittivity/ion_fraction_suppression_coefficient']=suppression[0]
    def visit(value):
        if isinstance(value,dict):
            if value.get('identity') in replacements:
                value['value']['magnitude']=replacements[value['identity']]
            for child in value.values():
                visit(child)
        elif isinstance(value,list):
            for child in value:
                visit(child)
    visit(doc)
    return doc


def checks_record(checks, start):
    return [dict(start=start, group=c.name, unit=c.unit, maximum=c.maximum, mean=c.mean,
        maximum_threshold=c.maximum_threshold, mean_threshold=c.mean_threshold,
        baseline_maximum=c.baseline_maximum, baseline_mean=c.baseline_mean, argmax=c.argmax,
        status=c.status, report_only=c.report_only, failed_rows=';'.join(map(str,c.failed_rows))) for c in checks]


def fill_predictions(summary):
    """Read native predictions at each retained point; never optimize or rehold."""
    from types import SimpleNamespace
    rows=read(OUT/f'{STEM}-rows.csv')
    for start in sorted({int(v['start']) for v in rows}):
        if start==-1:
            doc=BASE
        else:
            d=summary['starts'][start]
            doc=stored_document(SimpleNamespace(chart=SimpleNamespace(physical=d['physical'])))
        for batch in ('aqueous','methanol','ethanol'):
            components=COMPONENTS if batch=='aqueous' else ('water',batch,*IONS)
            source=Parameters.from_mapping(BASE,components=components)
            stored=Parameters.from_mapping(doc,components=components)
            _,_,obs,meta,_=declarations(source,batch)
            values=_core.observations_at(stored._native,[],[],[],obs,[0.0]*len(obs))
            for v in rows:
                if int(v['start'])==start and v['batch']==batch:
                    v.update(meta[int(v['row'])])
                    v['prediction']=values[int(v['row'])].observed
    write('rows',rows)


def main():
    p=Parameters.from_mapping(BASE,components=COMPONENTS)
    cs, hs, observations, metadata, checks=declarations(p)
    starts=[[float(v[n]) for n in NAMES] for v in read(HERE/'inputs/stage1-starts.csv')]
    controls=r.FitControls(maximum_iterations=200,maximum_elapsed_time_seconds=120,
                           function_tolerance=1e-12,gradient_tolerance=1e-10,parameter_tolerance=1e-12)
    weights=[1/(.02*math.log(10))**2 if v['set']=='T' else 0.0 for v in metadata]
    batches=[(p,observations,metadata,checks)]
    for solvent in ("methanol","ethanol"):
        q=Parameters.from_mapping(BASE,components=("water",solvent,*IONS))
        _,_,obs,meta,cgs=declarations(q,solvent)
        batches.append((q,obs,meta,cgs))
    print('native multistart:',len(starts),'starts,',len(observations),'declared rows',flush=True)
    begin=time.perf_counter()
    results=r.fit(p,cs,observations,weights=weights,holds=hs,checks=checks,baseline=p,controls=controls,
        singular_vectors=True,starts=starts,shared_iterations=4000,shared_seconds=2400)
    elapsed=time.perf_counter()-begin
    diagnostics, rows, groups=[],[],[]
    for q,obs,meta,cgs in batches:
        bc=r.evaluate_checks(q,obs,cgs,q)
        groups.extend(checks_record(bc,-1))
        for c in bc:
            for i,e in zip(c.rows,c.errors):
                rows.append(meta[i]|dict(start=-1,prediction=None,error=e,unit=c.unit,status=c.status))
    for start,result in enumerate(results):
        print(start,result.status,result.final_cost,result.acceptance,result.message,flush=True)
        covariance={key:getattr(result.covariance,key) for key in ('status','assumption','rank','singular_values','right_vectors',
            'projected_gradient_norm','active_bounds','optimizer','physical','variance_factor')}
        diagnostics.append(dict(start=start,status=result.status,message=result.message,usable=result.usable,
            optimizer=result.optimizer,physical=result.physical,cost=result.final_cost if result.usable else None,acceptance=result.acceptance,
            iterations=result.iterations,residual_evaluations=result.residual_evaluations,jacobian_evaluations=result.jacobian_evaluations,
            covariance=covariance,holds=[{key:getattr(h,key) for key in ('energy_error','correction','uncertainty','derivative','iterations')} for h in result.chart.holds],
            failures=[dict(row=f.observation,cause=f.message,status=f.status.name,optimizer=f.optimizer,jacobian=f.jacobian) for f in result.failures]))
        if result.usable:
            for q,obs,meta,cgs in batches:
                fitted=Parameters.from_mapping(stored_document(result),components=q.component_ids)
                cg=r.evaluate_checks(fitted,obs,cgs,q)
                groups.extend(checks_record(cg,start))
                for c in cg:
                    for i,e in zip(c.rows,c.errors):
                        rows.append(meta[i]|dict(start=start,prediction=None,error=e,unit=c.unit,status=c.status))
            diagnostics[-1]['acceptance']='pass' if all(v['status']=='pass' for v in groups if v['start']==start and not v['report_only']) else 'fail'
    write('rows',rows)
    write('groups',groups)
    usable=[(i,v) for i,v in enumerate(results) if v.usable and len(v.chart.holds)==3]
    selected,result=min(usable,key=lambda pair:(diagnostics[pair[0]]['acceptance']!='pass',pair[1].final_cost))
    path=OUT/f'{STEM}-parameters.json'
    path.write_text(json.dumps(stored_document(result),indent=2)+'\n')
    stored=Parameters.from_json(path,components=COMPONENTS)
    fresh=r.evaluate_checks(stored,observations,checks,p)
    for a,b in zip(result.checks,fresh):
        assert a.status==b.status and a.failed_rows==b.failed_rows
        for x,y in zip(a.errors,b.errors):
            assert abs(x-y)<=1e-12*max(abs(y),1),(a.name,x,y)
    write('stored-groups',[v for v in groups if v['start']==selected])
    summary=dict(method='native Ceres exact CppAD columns and implicit Born holds',selected=selected,starts=diagnostics,
        declared_starts=len(starts),executed_starts=len(results),seconds=elapsed,coordinate_order=NAMES,
        conditions=dict(temperature_K=T,pressure_Pa=P,branch='liquid',molality='mol/kg solvent'),
        controls=dict(per_start_iterations=200,per_start_seconds=120,shared_iterations=4000,shared_seconds=2400,
            function_tolerance=1e-12,gradient_tolerance=1e-10,parameter_tolerance=1e-12),accepted=None,
        inputs={str(v.relative_to(REPO)):hashlib.sha256(v.read_bytes()).hexdigest() for v in [SOURCE,HERE/'inputs/stage1-starts.csv',
            HERE/'inputs/kcl-hamer-wu.csv',HERE/'inputs/bromide-hamer-wu.csv',REPO/'analyses/2026-held-2014-chlorides/inputs/evaluated.csv',
            *[v for n in range(2,10) for v in sorted((FIG/f'figure-{n:02}/source').glob('*.csv')) if n in (2,3,6,7,8,9)]]},
        engine_provenance=list(provenance()),
        stored_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),base_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=REPO,text=True).strip())
    (OUT/f'{STEM}-fit.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    from stage1_diagnosis import state_checks, diagnose
    state=[]
    for q,obs,meta,cgs in batches:
        fitted=Parameters.from_json(path,components=q.component_ids)
        state.extend(state_checks(q,fitted,obs,meta,selected))
    write('states',state)
    numerical_pass=all(v['passed'] for v in state) and all(v['status']=='pass' for v in groups if v['start']==selected and not v['report_only'])
    summary['accepted']=numerical_pass
    (OUT/f'{STEM}-fit.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    fill_predictions(summary)
    diagnose(summary,rows,groups)
    if numerical_pass:
        raise RuntimeError('Passing optimum: D8 packet publication and stored-packet recomputation are required before completion')


if __name__=='__main__':
    main()
