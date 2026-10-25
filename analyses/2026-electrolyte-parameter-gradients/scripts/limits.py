"""F1a and F4 verification algebra, independent of the production derivative contraction."""
import math
from epcsaft import Mixture, Parameters, regression as r
from declarations import ROOT, parameters, row, row_groups
from run import native_start, manufactured, prediction_columns


def raw(p, cs, ob):
    return prediction_columns(p, cs, native_start(p, cs), [], [ob])


def identities(retain):
    p, observations = list(row_groups())[7]
    cs = [r.coordinate(p, [('ion_suppression_coefficient', target, 1) for target in p.component_ids[1:]], origin=0, scale=1)]
    c = native_start(p, cs)[0]
    for i, ob in enumerate(observations):
        value = r.predict(Mixture(p), ob).value
        xion = sum(ob.composition[1:])
        expected = -value*xion/(1+c*xion)
        retain('F1a', f'epsilon-C-{i}', raw(p, cs, ob)[0], expected, 1e-12*abs(expected))
    p = parameters(('water','sodium-cation','chloride-anion'))
    cs = [r.coordinate(p, [('ion_suppression_coefficient', target, 1) for target in p.component_ids[1:]], origin=0, scale=1),
          r.coordinate(p, [('k_ij', ('sodium-cation','chloride-anion'), 1)], origin=0, scale=1)]
    ob = row(p, 'solvation_gibbs_energy')
    for c, value in enumerate(raw(p, cs, ob)):
        retain('F1a', f'reference-zero-{c}', value/(8.31446261815324*298.15), 0, 1e-9)
    p = Parameters.from_json(ROOT/'analyses/2026-dielectric-alternatives/inputs/selections/density-temperature.json', components=('water','sodium-cation','chloride-anion'))
    cs = [r.coordinate(p,[('dispersion_energy_over_k', 'water', 1)],origin=0,scale=1)]
    ob = row(p,'bulk_relative_permittivity',x={'water':.96,'sodium-cation':.02,'chloride-anion':.02})
    u = [353.95]; step = u[0]*1e-4
    cd = (manufactured(p,cs,[u[0]+step],[],[ob])[0].observed-manufactured(p,cs,[u[0]-step],[],[ob])[0].observed)/(2*step)
    exact = raw(p,cs,ob)[0]
    retain('F1a','epsilon-density',exact,cd,1e-6*abs(cd)); assert exact != 0


def debye(retain):
    p = parameters(('water','sodium-cation','chloride-anion'))
    cs = [r.coordinate(p,[('k_ij', pair, 1)],origin=0,scale=1) for pair in (('sodium-cation','water'),('sodium-cation','chloride-anion'))]
    cs += [r.coordinate(p,[('ion_suppression_coefficient', target, 1) for target in p.component_ids[1:]],origin=0,scale=1),
           r.coordinate(p,[('born_diameter', 'sodium-cation', 1)],origin=0,scale=1),
           r.coordinate(p,[('debye_huckel_diameter', 'sodium-cation', 1)],origin=0,scale=1)]
    values = [raw(p,cs,row(p,'mean_ionic_activity_coefficient',m=m)) for m in (1e-5,1e-6)]
    for c in range(5): retain('F4',f'linear-{c}',values[1][c]/1e-6,values[0][c]/1e-5,.01*abs(values[0][c]/1e-5))
    cs = [r.coordinate(p,[('dispersion_energy_over_k', 'water', 1)],origin=0,scale=1)]
    density = r.density(p,298.15,1e5,{'water':1},phase='liquid',observed=1,form='difference',scale=1,mass_basis=False)
    pure = Mixture(p).state(298.15,P=1e5,x=(1,0,0),phase='liquid')
    drho = raw(p,cs,density)[0]; rho = pure.molar_density
    # Debye limiting law; SI constants from core/physical_constants.hpp, water mass from packet.
    e, eps0, kb, na, mw = 1.602176634e-19, 8.8541878128e-12, 1.380649e-23, 6.02214076e23, .0180153
    dielectric = eps0*pure.bulk_relative_permittivity
    a = e*e/(8*math.pi*dielectric*kb*298.15)*math.sqrt(2*e*e*na*rho*mw/(dielectric*kb*298.15))
    expected = -.5*a*drho/rho
    columns = [raw(p,cs,row(p,'mean_ionic_activity_coefficient',m=m))[0]/math.sqrt(m) for m in (1e-5,1e-6)]
    retain('F4','density-limit',columns[1],expected,.01*abs(expected))
    retain('F4','density-approach',abs(columns[1]-expected),0,abs(columns[0]-expected))
    return {'A':a,'density':rho,'density_column':drho,'expected':expected,'columns':columns}
