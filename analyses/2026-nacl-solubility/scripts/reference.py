"""Algebraic conversion of published Pitzer evaluated coefficients, not an EOS fit."""
from probe import *
from itertools import product

# Pitzer, Peiper & Busey (1984), Table A-3 p.26 and Table A-8 p.53, 1 bar.
# Entries: T, A_L/(RT), 10^3 beta0L, 10^3 beta1L, 10^3 C_L; gamma at m=1,3,6.
INPUT = []
for row in csv.DictReader((SOURCE / 'inputs/source/pitzer-1984.csv').open()):
    r = {key: float(value) for key, value in row.items()}
    INPUT.append((r['T_K'], r['AL_over_RT'], r['beta0L_times_1000'],
                  r['beta1L_times_1000'], r['CL_times_1000'],
                  [r[f'gamma_m{m}'] for m in (1, 3, 6)],
                  [r[f'phiL_over_RT_m{m}'] for m in (1, 3, 6)]))


def l2_over_rt(t,m,a,b0,b1,c):
    # Eq. (27), p.5. The Table A-3 A_L is the current A_H in Eq. (27).
    y=2*math.sqrt(m)
    q=math.sqrt(m)/(1+1.2*math.sqrt(m))+2/1.2*math.log1p(1.2*math.sqrt(m))
    return .5*a*q-2*t*m*(2*b0+b1/(2*m)*(1-(1+y-2*m)*math.exp(-y))+3*m*c)

def phi_l_over_rt(t,m,a,b0,b1,c):
    # Eq. (14)-(16); Table A-13 is apparent, not partial molar enthalpy.
    y=2*math.sqrt(m)
    return a*math.log1p(1.2*math.sqrt(m))/1.2-2*t*(
        m*(b0+b1/(2*m)*(1-(1+y)*math.exp(-y)))+m*m*c)

def run_reference():
    rows,inputs=[],[]
    native=list(csv.DictReader((OUT/'fixed-molality.csv').open()))
    for t,a,b0,b1,c,gammas,phi in INPUT:
        inputs.append(dict(T_K=t,P_Pa=100000,AL_over_RT=a,beta0L=b0/1000,
            beta1L=b1/1000,CL=c/1000,gamma_m1=gammas[0],gamma_m3=gammas[1],
            gamma_m6=gammas[2],source='Table A-3 p.26; A-8 p.53; 1 bar'))
        for i,m in enumerate([1.,3.,6.]):
            coefs=[a,b0/1000,b1/1000,c/1000]
            l=l2_over_rt(t,m,*coefs)
            # Half-last-digit envelope, separate from evaluated-model uncertainty.
            errors=[.0005,.0000005,.0000005,.0000005]
            corners=[l2_over_rt(t,m,*[v+s*e for v,s,e in zip(coefs,signs,errors)])
                     for signs in product([-1.,1.],repeat=4)]
            u=max(abs(v-l) for v in corners)
            phicalc=phi_l_over_rt(t,m,*coefs)
            phicorners=[phi_l_over_rt(t,m,*[v+s*e for v,s,e in zip(coefs,signs,errors)])
                        for signs in product([-1.,1.],repeat=4)]
            uphi=max(abs(v-phicalc) for v in phicorners)+.0005
            assert abs(phicalc-phi[i])<=uphi
            n=next(r for r in native if r['variant']=='primary' and float(r['T_K'])==t and float(r['m_mol_kg'])==m)
            rows.append(dict(T_K=t,m_mol_kg=m,reference_P_Pa=100000,
                reference_gamma=gammas[i],reference_ln_gamma=math.log(gammas[i]),
                reference_L2_over_RT=l,reference_L2_J_mol=R*t*l,
                reference_dln_gamma_dT=-l/(2*t),reference_derivative_rounding_bound=u/(2*t),
                engine_gamma=float(n['gamma']),engine_ln_gamma=float(n['ln_gamma']),
                engine_dln_gamma_dT=float(n['dln_gamma_dT']),engine_L2_J_mol=float(n['L2_J_mol']),
                source_phiL_over_RT=phi[i],reconstructed_phiL_over_RT=phicalc,
                phiL_rounding_bound=uphi))
    save('pitzer-inputs',inputs);save('reference-comparison',rows)
    for row in rows:print(json.dumps(row))

if __name__=='__main__':run_reference()
