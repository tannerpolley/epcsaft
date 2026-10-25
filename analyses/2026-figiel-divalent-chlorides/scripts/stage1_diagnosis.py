"""Read native stage-1 bands and guards; check stored states without fitting."""
import math
from epcsaft import Mixture, _core


def state_checks(source, stored, observations, metadata, selected):
    from stage1_fit import T, P
    COMPONENTS=stored.component_ids
    rows, seen = [], set()
    model = Mixture(stored)
    # Native formula-unit composition supplies the declared state inputs.
    for ob, meta in zip(observations, metadata):
        x=list(ob.composition)
        if meta['set'] in ('T','G1','G3','G4'):
            x=_core.formula_composition(x,stored._native,ob.ion_indices,ob.stoichiometry,ob.molality)
        for composition in (x,list(ob.reference_composition)):
            if not composition or tuple(composition) in seen:
                continue
            seen.add(tuple(composition))
            s=model.state(T,P=P,x=composition,phase='liquid')
            charge=sum(composition[COMPONENTS.index(c)]*(1 if c.endswith('cation') else -1)
                       for c in COMPONENTS if c.endswith(('cation','anion')))
            pressure_error=abs(s.pressure-P)
            rows.append(dict(start=selected,row=meta['row'],quantity='state',composition=';'.join(map(str,composition)),
                value=pressure_error,reference=.01,relative_error=None,pressure_density_derivative=s.pressure_density_derivative,
                charge=charge,passed=pressure_error<.01 and s.pressure_density_derivative>0 and charge==0))
    for solvent in (v for v in COMPONENTS if v in ("water","methanol","ethanol")):
        x=[float(c==solvent) for c in COMPONENTS]
        a,b=[Mixture(p).state(T,P=P,x=x,phase='liquid') for p in (stored,source)]
        rows.append(dict(start=selected,row=None,quantity=f'I pure {solvent} density',composition=';'.join(map(str,x)),
            value=a.molar_density,reference=b.molar_density,relative_error=abs(a.molar_density/b.molar_density-1),
            pressure_density_derivative=None,charge=None,passed=abs(a.molar_density/b.molar_density-1)<=1e-12))
        i=COMPONENTS.index(solvent)
        delta=abs(a.log_fugacity_coefficient[i]-b.log_fugacity_coefficient[i])
        relative=delta/max(abs(b.log_fugacity_coefficient[i]),1e-300)
        rows.append(dict(start=selected,row=None,quantity=f'I pure {solvent} ln phi',composition=';'.join(map(str,x)),
            value=a.log_fugacity_coefficient[i],reference=b.log_fugacity_coefficient[i],relative_error=relative,
            pressure_density_derivative=None,charge=None,passed=relative<=1e-12))
    # The qualitative Debye check is a sign and finite remainder, not a claimed error bound.
    water=model.state(T,P=P,x=[float(c=='water') for c in COMPONENTS],phase='liquid')
    A=math.sqrt(2*math.pi*water.molar_density*.0180153*6.02214076e23)*(1.602176634e-19**2/(4*math.pi*8.8541878128e-12*78.09*1.380649e-23*T))**1.5
    for ion in (v for v in COMPONENTS if v.endswith("cation")):
        ref=model.reference(T,P,solvent_fractions=[float(c=='water') for c in COMPONENTS],ions=(ion,'chloride-anion'),stoichiometry=(1,1))
        ln_gamma=ref.activity(.001).mean_ionic_log_molality_scale
        rows.append(dict(start=selected,row=None,quantity=f'Debye {ion} at .001 mol/kg',composition='pure water solvent',
            value=ln_gamma,reference=-A*math.sqrt(.001),relative_error=ln_gamma+A*math.sqrt(.001),
            pressure_density_derivative=None,charge=None,passed=math.isfinite(ln_gamma) and ln_gamma<0))
    for v in rows:
        v["components"]=";".join(COMPONENTS)
        v["batch"]=metadata[0]["batch"]
    return rows


def diagnose(summary, rows, groups):
    from stage1_fit import write, OUT, STEM
    import json
    from epcsaft import Parameters, regression as r
    from stage1_fit import SOURCE, COMPONENTS, declarations
    source=Parameters.from_json(SOURCE,components=COMPONENTS)
    cs,hs,observations,metadata,_=declarations(source)
    initial=[]
    for result in summary['starts']:
        if result['usable']:
            continue
        for failure in result['failures']:
            chart=r.evaluate_coordinates(source,cs,failure['optimizer'],hs)
            model=Mixture(Parameters(chart.parameters,source.component_ids,source.site_ids))
            i=failure['row']
            prediction=r.predict(model,observations[i])
            assert prediction.value is not None and prediction.value<=0
            initial.append(metadata[i]|dict(start=result['start'],prediction=prediction.value,
                observation_status=prediction.status.name,message=prediction.message))
    if initial:
        write('initial-failures',initial)
    selected=summary['selected']
    failed=[v for v in groups if v['start']==selected and v['status']!='pass' and not v['report_only']]
    bands=[]
    for setname in ('T','G1'):
        for name in sorted({v['group'] for v in rows if v['set']==setname}):
            for lo,hi in ((0,.1),(.1,1),(1,3),(3,4.5)):
                band=[v for v in rows if v['start']==selected and v['group']==name and lo<v['abscissa']<=hi]
                worst=max(band,key=lambda v:abs(v['error']))
                bands.append(dict(group=name,lower_molality=lo,upper_molality=hi,rows=len(band),maximum=abs(worst['error']),
                                  argmax_molality=worst['abscissa'],signed_error=worst['error']))
    write('bands',bands)
    (OUT/f'{STEM}-diagnosis.json').write_text(json.dumps(dict(failed_groups=failed,
        conclusion='Ceres least squares does not impose D6 guards. The guarded question needs #219 slice B native Ipopt minimax; feasibility is unresolved.',
        coordinate_or_tolerance_change=False,stage2_started=False),indent=2)+'\n')
