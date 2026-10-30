from pathlib import Path
import runpy,json,sys,contextlib,io
import epcsaft
from epcsaft import equilibrium as eq
from epcsaft.eos import property_asl_record
from pyomo.environ import ExternalFunction
root=Path(sys.argv[1]);out=Path(sys.argv[2]);out.mkdir(exist_ok=True,parents=True)
with contextlib.redirect_stdout(io.StringIO()) as captured:d=runpy.run_path(str(root/'conda-test.py'))
(out/'conda-example.stdout').write_text(captured.getvalue())
r=d['bubble'];a=d['batch'];mix=d['mix']
def encode(value):
 if isinstance(value,float):return value.hex()
 if value is None or isinstance(value,(str,int,bool)):return value
 if isinstance(value,(tuple,list)):return [encode(x) for x in value]
 return str(value)
result={name:encode(getattr(r,name)) for name in ['pressure','temperature','amounts','mole_fractions','molar_densities','residuals','classified_residuals','iterations','evaluations','hessian_calls','success','requested_tolerance_met']}
result['action']=encode(a.results[0].action)
result['action_central_residual']=encode(a.diagnostics.central_residual)
result['action_linear_residual']=encode(a.diagnostics.maximum_linear_residual)
row=next(x for x in epcsaft.tolerances() if x.name=='EQUILIBRIUM_ROW_RESIDUAL')
result['criterion_name']=r.convergence.criterion.name
result['criterion_is_shared']=r.convergence.criterion is next(x for x in epcsaft.tolerances() if x.name==r.convergence.criterion.name)
import ctypes
native=ctypes.CDLL(str(Path(epcsaft._core.__file__).parent/'epcsaft_native.dll'))
module=ctypes.CDLL(epcsaft._core.__file__)
class Record(ctypes.Structure):
    _fields_=[('name',ctypes.c_void_p)]
# Address identity is checked by the native caller against the imported default table.
result['tolerance']=encode(row.value)
try:mix.state(200,P=-1,x=[0.3,0.7],phase='liquid')
except Exception as error:result['exception']={'type':type(error).__name__,'message':str(error)}
else:raise AssertionError('negative pressure was not refused')
callback=ExternalFunction(library=str(Path(epcsaft._core.__file__).parent/'epcsaft_asl.dll'),function='epcsaft_pressure_property')
sys.path.insert(0,str((root).resolve()))
from tests.native_inputs import record,p3_forms
thermal=epcsaft.Mixture(d['parameters'],thermochemistry=record(p3_forms()['nasa9'],p3_forms()['nasa9'],lower=150,upper=450))
text=property_asl_record(thermal)
result['ASL']=[encode(callback.evaluate_fgh([200.0,1e6,0.3,0.7,text,str(row),str(phase)],fgh=2)) for phase in [0,1] for row in [0,1,2]]
(out/'values.json').write_text(json.dumps(result,indent=2))
print(json.dumps(result,indent=2))
from ctypes import wintypes as w
import hashlib
kernel=ctypes.WinDLL('kernel32',use_last_error=True);psapi=ctypes.WinDLL('psapi',use_last_error=True)
kernel.GetCurrentProcess.restype=w.HANDLE
psapi.EnumProcessModules.argtypes=[w.HANDLE,ctypes.POINTER(w.HMODULE),w.DWORD,ctypes.POINTER(w.DWORD)]
psapi.GetModuleFileNameExW.argtypes=[w.HANDLE,w.HMODULE,w.LPWSTR,w.DWORD]
modules=(w.HMODULE*2048)();needed=w.DWORD();process=kernel.GetCurrentProcess()
assert psapi.EnumProcessModules(process,modules,ctypes.sizeof(modules),ctypes.byref(needed))
identities=[]
for handle in modules[:needed.value//ctypes.sizeof(w.HMODULE)]:
 path=ctypes.create_unicode_buffer(32768)
 assert psapi.GetModuleFileNameExW(process,handle,path,len(path))
 if any(s in path.value.lower() for s in ('epcsaft','ipopt','mumps','openblas')):
  p=Path(path.value);identities.append({'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
print('loaded images',json.dumps(identities,indent=2),flush=True)
(out/'loaded-images.json').write_text(json.dumps(identities,indent=2))
assert any(Path(x['path']).name=='epcsaft_ipopt_3119.dll' for x in identities)
assert not any(Path(x['path']).name.lower()=='ipopt.dll' for x in identities)
from epcsaft.eos import PressureActionDirection,PropertyObservable
observables=[PropertyObservable.TotalEnthalpy,PropertyObservable.LogFugacityCoefficient]
t=PressureActionDirection(1.0,0.0,[0.0,0.0]);p=PressureActionDirection(0.0,1.0,[0.0,0.0])
comparisons=[]
for phase,name in enumerate(['liquid','vapor']):
 first=[thermal.pressure_actions(200.0,1e6,x=[0.3,0.7],phase=name,observables=observables,directions=d) for d in [[t],[p]]]
 second=[thermal.pressure_actions(200.0,1e6,x=[0.3,0.7],phase=name,observables=observables,directions=d) for d in [[t,t],[t,p],[p,p]]]
 native=[[first[0].molar_density,*(b.density_action[0] for b in first),*(b.density_mixed_action for b in second)]]
 native += [[first[0].results[k].value,*(b.results[k].action for b in first+second)] for k in [0,1]]
 for row in range(3):
  value,gradient,hessian=callback.evaluate_fgh([200.0,1e6,0.3,0.7,text,str(row),str(phase)],fgh=2)
  asl=[value,*gradient[:2],*hessian[:3]]
  comparisons.append(dict(phase=name,row=row,native=encode(native[row]),ASL=encode(asl)))
  (out/'asl-native-comparisons.json').write_text(json.dumps(comparisons,indent=2))
  print('ASL native comparison',json.dumps(comparisons[-1]),flush=True)
  assert [float(x).hex() for x in asl]==[None if x is None else x.hex() for x in native[row]]
