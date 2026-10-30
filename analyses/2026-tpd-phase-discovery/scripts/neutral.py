"""#391 neutral sweep: methane/ethane bubble declarations at 10 temperatures x 7 initial pressures. Native solves only.
usage: neutral.py OUT.json"""
import sys,json,time
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[3]/'engine'))
import epcsaft
from epcsaft import Mixture,equilibrium as eq
from tests.native_inputs import parameters
sys.path.insert(0,str(Path(__file__).resolve().parents[3]/'tools'))
from engine_identity import loaded_identity
model=Mixture(parameters([(16.043,1.,3.7039,150.03),(30.07,1.6069,3.5206,191.42)], [('k_ij',(0,1),0.)]))
def problem(T,P0):
 return eq.Problem(T=T,P=eq.Free(P0),feed=eq.Amounts({'s0':.5,'s1':.5}),phases=[eq.Phase('L',kind='liquid'),eq.Phase('V',kind='vapor',amount=eq.Pinned(0.))])


def main():
 rows=[]
 # Raw and classified rows are recorded separately; None when the result has no rows.
 peak=lambda v:max(map(abs,v)) if v else None
 for T in [240.,250.,255.,260.,262.,264.,265.,265.5,265.95,266.5]:
  for P0 in [1e5,1e6,2.02606e6,3e6,5e6,6.75e6,7e6]:
   p=problem(T,P0)
   st=time.perf_counter()
   try:
    r=eq.solve_equilibrium(model,p)
    rows.append(dict(T=T,P0=P0,success=r.success,iterations=r.iterations,evaluations=r.evaluations,hessians=r.hessian_calls,message=r.message,P=[] if r.pressure is None else [r.pressure],x=r.mole_fractions,rho=r.molar_densities,residual=r.convergence.residual,criterion=r.convergence.criterion.value,raw=peak(r.residuals),classified=peak(r.classified_residuals),seconds=time.perf_counter()-st,numerical_operands=dict(coordinates=r.coordinates,amounts=r.amounts,packing=r.packing_fractions)))
   except Exception as ex: rows.append(dict(T=T,P0=P0,success=False,message=str(ex),residual=None,raw=None,classified=None,seconds=time.perf_counter()-st))
   if 49 <= len(rows)-1 < 63 and rows[-1]['success']:
    estimate=eq.compile_problem(model,p).numerical_estimates(r.coordinates,r.temperature,r.pressure,r.amounts,r.mole_fractions,r.molar_densities,r.packing_fractions)
    rows[-1]['numerical_estimates']={field:[{key:getattr(value,key) for key in ('absolute','chart','available','resolved','cause')} for value in ([estimate.pressure] if field=='pressure' else getattr(estimate,field))] for field in ('pressure','densities','composition')}
 Path(sys.argv[1]).write_text(json.dumps(dict(engine=loaded_identity(),rows=rows),indent=2,default=str)+'\n')
 print(loaded_identity(),len(rows),'cases',sum(r['success'] for r in rows),'successful',sum(r['seconds'] for r in rows),'seconds')


if __name__ == "__main__":
 main()
