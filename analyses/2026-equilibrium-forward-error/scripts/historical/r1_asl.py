from r1_common import *
import ctypes
from pyomo.core.base.external import _ARGLIST
from pyomo.environ import ExternalFunction
from epcsaft import IdealPolynomial, PressureActionDirection, PropertyObservable
from epcsaft.eos import property_asl_record
from epcsaft import equilibrium as eq
from tests.native_inputs import record
from tests.equilibrium.test_cross_association_declarations import document,NONE
model=Mixture(Parameters.from_mapping(document([NONE])),thermochemistry=record(IdealPolynomial([30.],300.,0.,0.),IdealPolynomial([30.],300.,0.,0.),upper=351.))
def invoke(name, args,order):
    external=ExternalFunction(library=str(SHARED/'libepcsaft_asl.so'),function=name)
    external.load_library()
    packed=_ARGLIST(args,fgh=order)
    v=external._known_functions[name][0](ctypes.byref(packed))
    return dict(value=v,message=packed.Errmsg.decode() if packed.Errmsg else '',gradient=[packed.derivs[i] for i in range(packed.nr)] if order>=1 else [],hessian=[packed.hes[i] for i in range(packed.nr*(packed.nr+1)//2)] if order==2 else [])
rows={}
for T in [350.,351.]:
    batch=model.pressure_actions(T,1e5,x=(.8,.2),phase='liquid',observables=[PropertyObservable.TotalEnthalpy],directions=[PressureActionDirection(1.,0.,[0.,0.],[])])
    d=batch.diagnostics
    rows[str(T)]=dict(native=dict(density=batch.molar_density,density_status=str(batch.density_status),status=str(batch.results[0].status),value=batch.results[0].value,action=batch.results[0].action,site_rank=d.site_rank,site_pivot_ratio=d.site_pivot_ratio,association_residual=d.association_residual),asl=[invoke('epcsaft_pressure_property',[T,1e5,.8,.2,property_asl_record(model),'1','0'],k) for k in (0,1,2)])
    print(T,rows[str(T)],flush=True)
save('r1-asl',rows)
# Two-liquid result replacement/cause conversion at require_central.
path=ROOT/'data/packets/khudaida-2026-figure-2-electrolyte-lle/7/parameters/parameters.json'
lle=Mixture(Parameters.from_json(path,components=['water','isobutanol']))
def p(feed):return Problem(T=298.15,P=1e5,feed=Amounts({'water':feed[0],'isobutanol':feed[1]}),phases=[Phase('L1',kind='liquid',support=('water','isobutanol'),composition_guess=(.46,.54)),Phase('L2',kind='liquid',support=('water','isobutanol'),composition_guess=(.99,.01))])
opt=eq.Continuation(initial_step=.25,minimum_step=.25,max_steps=4,max_retries_per_step=1)
for feed in [(0.7,.3),(.999,.001),(.001,.999)]:
    problem=p(feed);r=solve_equilibrium(lle,problem)
    try:
        encoded=eq.asl_record(lle,problem,opt,[0.,0.],[[1.,0.],[0.,1.]],list(feed),1,.1)
        rows['lle'+str(feed)]=dict(native=summary(r),asl=invoke('epcsaft_equilibrium_phase_component_amount',[*feed,encoded,'0','0'],2))
    except Exception as error:
        rows['lle'+str(feed)]=dict(native=summary(r),error=str(error))
    print(feed,rows['lle'+str(feed)],flush=True)
encoded=eq.asl_record(lle,p((.7,.3)),opt,[0.,0.],[[1.,0.],[0.,1.]],list((.7,.3)),1,.1)
rows['continuation-fallback']=invoke('epcsaft_equilibrium_phase_component_amount',[.001,.999,encoded,'0','0'],2)
print('continuation-fallback',rows['continuation-fallback'],flush=True)
save('r1-asl',rows)
