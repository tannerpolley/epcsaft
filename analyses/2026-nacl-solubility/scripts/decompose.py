"""Additive chemical-potential decomposition at the primary model's common density path."""
from probe import *

def masked(t,mask,variant='primary'):
    d=document(variant,t)
    if mask=='full':return mixture(d)
    for c in d['correlations']:
        if c['family']=='relative_permittivity':constant(c,1e50)
    if mask in ['no-dh-association','hard-chain']:
        d['topology']={}
        next(f for f in d['model_families'] if f['kind']=='association')['choice']='none'
    if mask=='hard-chain':
        e=next(v['value']['magnitude'] for c in d['components'] if c['component_id']=='water'
               for v in c['coefficients'] if v['family']=='dispersion_energy_over_k')
        f=1-.12*math.exp(-3*e/t)
        ft=-.36*e/t**2*math.exp(-3*e/t)
        dia=sigma(t)*f
        dia_t=sigma_t(t)*f+sigma(t)*ft
        b=dia_t/dia
        for c in d['correlations']:
            if c['family']=='segment_diameter':
                c['constant']['value']['magnitude']=dia/(2*.88)
                c['terms']=[{'amplitude':{'identity':'diagnostic/hard-chain/sigma/amplitude',
                     'value':{'magnitude':dia/(2*.88)*math.exp(-2*b*t),'unit':'angstrom'}},
                     'exponent_coefficient':{'identity':'diagnostic/hard-chain/sigma/exponent',
                     'value':{'magnitude':2*b,'unit':'1 / kelvin'}}}]
        for c in d['components']:
            for v in c['coefficients']:
                if v['family']=='dispersion_energy_over_k':v['value']['magnitude']=0.
    return mixture(d)

def raw(mix,t,m,rho,rho_t):
    x,_=composition(m)
    st=mix.state(t,rho=rho,x=x)
    # One mole of mixture, fixed V, addition of one mole NaCl: dn=(0,1,1).
    # This is a neutral direction. Euler terms convert molar U,H to salt partials.
    dx=[-2*x[0],1-2*x[1],1-2*x[2]]
    b=mix.density_actions(t,rho,x=x,
        observables=[Obs.ResidualInternalEnergy,Obs.ResidualEnthalpy],
        directions=[DensityActionDirection(0.,2*rho,dx)])
    assert all(a.action is not None for a in b.results),[str(a.status) for a in b.results]
    up=2*st.residual_internal_energy+b.results[0].action
    hp=2*st.residual_enthalpy+b.results[1].action
    mu=.5*sum(st.residual_chemical_potential_over_rt[1:])
    mu_t=-up/(2*R*t*t)+(hp-up)/(2*R*t)*rho_t/rho
    return mu,mu_t,dict(hc=st.hard_chain,disp=st.dispersion,assoc=st.association,dh=st.debye_huckel)

def run_decomposition():
    rows=[]
    for t in [298.15,323.15,353.15]:
        quiet_check()
        base=mixture(document())
        sf,bf,_=phi_action(base,t,6.2)
        r,r_t=sf.molar_density,bf.density_action[0]
        masks={key:masked(t,key) for key in ['full','no-dh','no-dh-association','hard-chain']}
        for m0 in [1e-8,1e-12,1e-16]:
            s0,b0,_=phi_action(base,t,m0)
            r0,r0_t=s0.molar_density,b0.density_action[0]
            v={}
            for key,mix in masks.items():
                f=raw(mix,t,6.2,r,r_t)
                z=raw(mix,t,m0,r0,r0_t)
                # All masks share the primary hard-chain geometry at this T.
                assert abs(f[2]['hc']-sf.hard_chain)<1e-12
                v[key]=(f[0]-z[0],f[1]-z[1])
            parts={'hard_chain':v['hard-chain'],
                'dispersion':tuple(a-b for a,b in zip(v['no-dh-association'],v['hard-chain'])),
                'association':tuple(a-b for a,b in zip(v['no-dh'],v['no-dh-association'])),
                'debye_huckel':tuple(a-b for a,b in zip(v['full'],v['no-dh'])),
                'born':(0.,0.),
                'density_conversion':(math.log(r/r0)-math.log(1+2*MW*6.2)+math.log(1+2*MW*m0),r_t/r-r0_t/r0)}
            for term,(g,gt) in parts.items():
                rows.append(dict(T_K=t,m_mol_kg=6.2,reference_m_mol_kg=m0,term=term,
                    ln_gamma_component=g,dln_gamma_dT_component=gt))
            rows.append(dict(T_K=t,m_mol_kg=6.2,reference_m_mol_kg=m0,term='sum',
                ln_gamma_component=sum(a[0] for a in parts.values()),
                dln_gamma_dT_component=sum(a[1] for a in parts.values())))
            print(t,m0,parts,flush=True)
            save('contribution-decomposition',rows)

def run_k_isolation():
    rows=[]
    for t in [298.15,323.15,353.15]:
        quiet_check()
        variants={}
        for variant in ['primary','local-k']:
            base=mixture(document(variant,t))
            sf,bf,_=phi_action(base,t,6.2)
            s0,b0,_=phi_action(base,t,1e-16)
            r,r_t=sf.molar_density,bf.density_action[0]
            r0,r0_t=s0.molar_density,b0.density_action[0]
            v={}
            for mask in ['full','no-dh','no-dh-association','hard-chain']:
                mix=masked(t,mask,variant)
                f=raw(mix,t,6.2,r,r_t);z=raw(mix,t,1e-16,r0,r0_t)
                v[mask]=f[1]-z[1]
            variants[variant]={'hard_chain':v['hard-chain'],
                'dispersion':v['no-dh-association']-v['hard-chain'],
                'association':v['no-dh']-v['no-dh-association'],
                'debye_huckel':v['full']-v['no-dh'],
                'density_conversion':r_t/r-r0_t/r0}
        for term in variants['primary']:
            rows.append(dict(T_K=t,term=term,primary_per_K=variants['primary'][term],
                frozen_k_per_K=variants['local-k'][term],
                k_law_effect_per_K=variants['primary'][term]-variants['local-k'][term]))
    save('k-contribution-isolation',rows)

if __name__=='__main__':
    run_decomposition()
    run_k_isolation()
