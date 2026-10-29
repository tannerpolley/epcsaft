"""Bounded diagnostic declarations; the prepared Engine owns all EOS solves/actions."""
from load_engine import ROOT, LIB
import copy
import csv
from datetime import datetime
import hashlib
import importlib.util
import json
import math
import os
import subprocess
from pathlib import Path
from zoneinfo import ZoneInfo

from epcsaft import (Mixture, Parameters, PressureActionDirection, DensityActionDirection,
                     PropertyObservable as Obs, GAS_CONSTANT_J_PER_MOL_K as R)
from epcsaft.equilibrium import (compile_problem, solve_equilibrium, solved_state_actions,
    SolvedStateObservable, SolvedStateObservableKind as Kind, SolvedStateActionDirection,
    SolvedStateActionRequest)

SOURCE = Path(__file__).resolve().parents[1]
OUT = SOURCE / 'results/temperature-trend'
SPECIES = ['water', 'sodium-cation', 'chloride-anion']
MW, P = .01801528, 101325.
BASE = json.loads((SOURCE / 'inputs/held-2-malmberg-maryott.json').read_text())
spec = importlib.util.spec_from_file_location('solubility_study', SOURCE / 'scripts/run.py')
study = importlib.util.module_from_spec(spec)
spec.loader.exec_module(study)
study.P = P
THERMO = study.inputs('clarke-glew-1985-table34.csv')

def quiet_check():
    now = datetime.now(ZoneInfo('America/Denver'))
    if 3 <= now.hour < 5:
        raise RuntimeError('Quiet-host campaign, 03:00-05:00 America/Denver: probes paused')

def constant(c, value):
    c.update(form='constant', independent_variables=[])
    c['constant']['value']['magnitude'] = value
    c.pop('terms', None)
    c.pop('reference_temperature', None)

def sigma(t):
    return 2.7927 + 10.11 * math.exp(-.01775*t) - 1.417 * math.exp(-.01146*t)

def sigma_t(t):
    return -10.11*.01775*math.exp(-.01775*t) + 1.417*.01146*math.exp(-.01146*t)

def eps(t):
    u=t-273.15
    return 87.740-.40008*u+.0009398*u*u-.000001410*u*u*u

def document(variant='primary', t=298.15):
    d=copy.deepcopy(BASE)
    if variant.startswith('local-'):
        freeze = variant[6:]
        ft=t
    else:
        freeze={'freeze-k':'k','freeze-sigma':'sigma','freeze-epsilon':'epsilon',
                'freeze-k-sigma':'k-sigma'}.get(variant, '')
        ft=298.15
    for pair in d['pairs']:
        if set((pair['component_id_a'],pair['component_id_b']))=={'water','sodium-cation'} and 'k' in freeze.split('-'):
            c=pair['coefficients'][0]
            c.pop('temperature_law')
            c['value']={'magnitude':2.37999-.007981*ft,'unit':'dimensionless'}
    for c in d['correlations']:
        if c['family']=='segment_diameter' and 'sigma' in freeze.split('-'):
            constant(c, sigma(ft))
        if c['family']=='relative_permittivity':
            if 'epsilon' in freeze.split('-'):
                constant(c, eps(ft))
            elif variant=='constant-78.09':
                constant(c,78.09)
            elif variant=='ascani':
                # Ascani & Held (2021), Table 5, DOI 10.1002/zaac.202100032.
                # Local Taylor polynomial of the logarithmic law, exact value and first slope.
                # Every evaluated T gets its own declaration; no higher-T derivative is claimed.
                c.update(form='constant-plus-polynomial-in-temperature',independent_variables=['temperature'])
                c['constant']['value']['magnitude']=677.480-105.2*math.log(t)
                c['reference_temperature']={'magnitude':t,'unit':'kelvin'}
                c['terms']=[{'power':1,'coefficient':{'identity':'diagnostic/ascani/slope',
                    'value':{'magnitude':-105.2/t,'unit':'1 / kelvin'}}}]
    return d

def mixture(d):
    return Mixture(Parameters.from_mapping(d, components=SPECIES))

def composition(m):
    nw=1/MW
    nt=nw+2*m
    x=[nw/nt,m/nt,m/nt]
    xm=[-2*nw/nt**2,nw/nt**2,nw/nt**2]
    return x,xm

def phi_action(mix,t,m,temperature=1.):
    x,xm=composition(m)
    state=mix.state(t,P=P,x=x,phase='liquid')
    direction=PressureActionDirection(temperature,0.,[0.,0.,0.] if temperature else xm)
    b=mix.pressure_actions(t,P,x=x,phase='liquid',anchor=state,
        observables=[Obs.LogFugacityCoefficient],directions=[direction])
    assert all(a.action is not None for a in b.results), [str(a.status) for a in b.results]
    return state,b,.5*sum(a.action for a in b.results[1:])

def saturation(mix,t):
    quiet_check()
    row=next(r for r in THERMO if r['temperature_K']==t)
    request=study.problem(row)
    result=solve_equilibrium(mix,request)
    assert result.success, result.message
    assert max(map(abs,result.classified_residuals))<1e-10
    n=[result.amounts[0]*v for v in result.mole_fractions]
    m=n[1]/(MW*n[0])
    compiled=compile_problem(mix,request)
    obs=[SolvedStateObservable(Kind.PhaseComponentAmount,0,i) for i in [0,1]]
    b=solved_state_actions(compiled,result,SolvedStateActionRequest(obs,
        [SolvedStateActionDirection(1.,0.,[0.,0.,0.],[])]))
    assert all(a.action is not None for a in b.results),[str(a.status) for a in b.results]
    nw_t,na_t=[a.action for a in b.results]
    mt=na_t/(MW*n[0])-m*nw_t/n[0]
    st,_,phit=phi_action(mix,t,m)
    _,_,phim=phi_action(mix,t,m,0.)
    gm=phim-2*MW/(1+2*MW*m)
    kt=row['h_J_mol']/(R*t*t)
    gt=.5*kt-(1/m+gm)*mt
    ref_t=phit-gt
    ref=mix.reference(t,P,solvent_fractions=[1.,0.,0.],ions=SPECIES[1:],stoichiometry=(1,1))
    g=ref.activity(m).mean_ionic_log_molality_scale
    closure=2*(math.log(m)+g)+row['g_J_mol']/(R*t)
    record=dict(T_K=t,P_Pa=P,m_sat=m,dm_sat_dT=mt,dlnK_dT=kt,
        dln_gamma_dm=gm,dln_gamma_dT=gt,mean_reference_logphi_dT=ref_t,
        ln_gamma=g,identity_ln=closure,
        classified_max=max(map(abs,result.classified_residuals)),
        pressure_error_Pa=st.pressure-P,solid_mol=result.solid_amounts[0],
        message=result.message)
    assert abs(closure)<1e-8
    return record,ref

def fixed(mix,ref,reference_t,t,m):
    st,b,phit=phi_action(mix,t,m)
    _,_,phim=phi_action(mix,t,m,0.)
    g=ref.activity(m).mean_ionic_log_molality_scale
    gt=phit-reference_t
    return dict(T_K=t,m_mol_kg=m,ln_gamma=g,gamma=math.exp(g),dln_gamma_dT=gt,
        L2_J_mol=-2*R*t*t*gt,dln_gamma_dm=phim-2*MW/(1+2*MW*m),
        rho_mol_m3=st.molar_density,rho_T=b.density_action[0],
        pressure_error_Pa=st.pressure-P,dp_drho=st.pressure_density_derivative,
        association_residual=b.diagnostics.association_residual)

def save(name,rows):
    if not rows:return
    with (OUT/(name+'.csv')).open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def run():
    quiet_check()
    sat,states=[],[]
    for t in [298.15,323.15,353.15]:
        for variant in ['primary','local-k','local-sigma','local-epsilon',
                        'local-k-sigma-epsilon','freeze-k','freeze-sigma','freeze-k-sigma',
                        'constant-78.09','ascani']:
            quiet_check()
            mix=mixture(document(variant,t))
            s,ref=saturation(mix,t)
            sat.append(dict(variant=variant,**s))
            for m in ([1.,3.,6.,6.2] if variant in ['primary','freeze-k','ascani'] else [6.2]):
                states.append(dict(variant=variant,**fixed(mix,ref,s['mean_reference_logphi_dT'],t,m)))
            print(json.dumps(sat[-1]),flush=True)
            save('saturation',sat);save('fixed-molality',states)
    identity=dict(commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
        engine_tree=subprocess.check_output(['git', 'rev-parse', 'HEAD:engine'], cwd=ROOT, text=True).strip(),
        native_extension=os.path.relpath(LIB, SOURCE),native_sha256=hashlib.sha256(LIB.read_bytes()).hexdigest(),
        native_library_sha256=hashlib.sha256((LIB.parent/'libepcsaft_native.so').read_bytes()).hexdigest(),
        input_sha256=hashlib.sha256((SOURCE/'inputs/held-2-malmberg-maryott.json').read_bytes()).hexdigest(),
        finished=datetime.now().astimezone().isoformat(),threads=1,
        diagnostic_script_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(Path(__file__).resolve().parent.glob('*.py'))},
        source_input_sha256={str(p.relative_to(SOURCE)):hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted((SOURCE/'inputs').rglob('*')) if p.is_file()})
    (OUT/'run-identity.json').write_text(json.dumps(identity,indent=2)+'\n')

if __name__=='__main__':run()
