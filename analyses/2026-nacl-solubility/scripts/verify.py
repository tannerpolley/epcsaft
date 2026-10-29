"""Finite differences verify previously computed exact actions only."""
from probe import *

def run_verification():
    states=list(csv.DictReader((OUT/'fixed-molality.csv').open()))
    rows=[]
    for variant,t in [('primary',298.15),('primary',323.15),('primary',353.15),
                      ('freeze-k',353.15),('ascani',353.15)]:
        quiet_check()
        mix=mixture(document(variant,t))
        center=next(r for r in states if r['variant']==variant and float(r['T_K'])==t and float(r['m_mol_kg'])==6.2)
        exact=float(center['dln_gamma_dT'])
        def gamma(temp,p=P):
            ref=mix.reference(temp,p,solvent_fractions=[1.,0.,0.],ions=SPECIES[1:],stoichiometry=(1,1))
            return ref.activity(6.2).mean_ionic_log_molality_scale
        fd=[(gamma(t+h)-gamma(t-h))/(2*h) for h in [.1,.05,.025]]
        rich=(4*fd[-1]-fd[-2])/3
        assert abs(rich-exact)<1e-9
        rows.append(dict(variant=variant,T_K=t,m_mol_kg=6.2,exact_per_K=exact,
            FD_0_1_K=fd[0],FD_0_05_K=fd[1],FD_0_025_K=fd[2],
            Richardson_per_K=rich,error_per_K=rich-exact))
    save('temperature-verification',rows)
    # Direct comparison bounds the 1 bar versus 1.01325 bar reference mismatch.
    pressure=[]
    mix=mixture(document())
    for t in [298.15,323.15,353.15]:
        for m in [1.,3.,6.]:
            g=[]
            for p in [100000.,P]:
                ref=mix.reference(t,p,solvent_fractions=[1.,0.,0.],ions=SPECIES[1:],stoichiometry=(1,1))
                g.append(ref.activity(m).mean_ionic_log_molality_scale)
            pressure.append(dict(T_K=t,m_mol_kg=m,ln_gamma_1bar=g[0],
                ln_gamma_1_01325bar=g[1],difference_ln=g[1]-g[0]))
    save('pressure-comparison',pressure)
    # Locally frozen declarations preserve the central values; their first actions add.
    law=[]
    for t in [298.15,323.15,353.15]:
        d={r['variant']:r for r in states if float(r['T_K'])==t and float(r['m_mol_kg'])==6.2}
        f=float(d['primary']['dln_gamma_dT'])
        parts={key:f-float(d['local-'+key]['dln_gamma_dT']) for key in ['k','sigma','epsilon']}
        explicit=float(d['local-k-sigma-epsilon']['dln_gamma_dT'])
        closure=sum(parts.values())+explicit-f
        assert abs(closure)<1e-12
        assert max(abs(float(d[v]['ln_gamma'])-float(d['primary']['ln_gamma']))
                   for v in ['local-k','local-sigma','local-epsilon','local-k-sigma-epsilon'])<1e-12
        law.append(dict(T_K=t,m_mol_kg=6.2,k_law=parts['k'],water_sigma_law=parts['sigma'],
            permittivity_law=parts['epsilon'],ion_diameter_laws=0.,
            explicit_temperature=explicit,sum=explicit+sum(parts.values()),closure=closure))
    save('parameter-law-decomposition',law)
    print(json.dumps(dict(thermal=rows,parameter_laws=law),indent=2))

if __name__=='__main__':run_verification()
