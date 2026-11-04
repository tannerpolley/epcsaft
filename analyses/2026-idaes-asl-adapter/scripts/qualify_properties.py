"""Pressure-property numerical falsifiers; run from the root with the installed Engine."""
from pathlib import Path
import csv
import json
import math
import runpy
import hashlib
import subprocess
import warnings
from dataclasses import replace

from epcsaft import Pure, Parameters, ThermochemistryPacket, PressureActionDirection, PropertyObservable, equilibrium as eq
from epcsaft.eos import property_asl_record
from epcsaft.idaes import PropertyExternalFunction, asl_library, solve_with_diagnostics
from pyomo.environ import value

ROOT = Path.cwd()
OUT = ROOT / 'build/issue-977/qualification-r2'
OUT.mkdir(parents=True, exist_ok=True)
namespace = runpy.run_path(str(ROOT/'docs/site/examples/idaes-vle-flash.py'))
m, model, s, solver = (namespace[k] for k in ('m', 'model', 'state', 'solver'))
states = (m.fs.flash.control_volume.properties_in[0], s)
options = namespace['options']
rows = []
failures = []
phase_diagnostics = []
transitions = []
thermal_boundary = []
def measure(case, quantity, actual, reference, tolerance):
    error = abs(actual-reference)
    rows.append(dict(case=case, quantity=quantity, actual=actual, reference=reference, error=error, tolerance=tolerance))
    with (OUT/'comparisons.csv').open('w') as f:
        w=csv.DictWriter(f, fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
    assert error<=tolerance,(case,quantity,actual,reference,error,tolerance)

def solve(case, checker=False):
    opts={**options,'max_iter':150}
    if checker:opts.update(derivative_test='second-order',derivative_test_tol=1e-4,derivative_test_perturbation=1e-6,point_perturbation_radius=0)
    opt=OUT/(case+'.opt');opt.write_text(''.join(f'{k} {v}\n' for k,v in opts.items()))
    solver.options['option_file_name']=str(opt)
    solve_with_diagnostics(solver,m,logfile=OUT/(case+'.log'))
    text=(OUT/(case+'.log')).read_text()
    if checker:assert 'No errors detected by derivative checker.' in text
    assert 'running with linear solver mumps' in text
    for k in opts:
        assert any(line.split() and line.split()[0]==k and line.split()[-1]=='yes' for line in text.splitlines()),k

def compare_flash(case, ref, fraction_tol=1e-8):
    for i,p in enumerate(s.phase_list):
        measure(case,p+'_fraction',value(s.phase_frac[p]),ref.amounts[i],fraction_tol)
        for j,c in enumerate(model.component_ids):
            measure(case,p+'_'+c,value(s.mole_frac_phase_comp[p,c]),ref.mole_fractions[2*i+j],1e-8)
    gap=max(abs(value(s.mole_frac_phase_comp['liquid',j])-value(s.mole_frac_phase_comp['vapor',j])) for j in s.component_list)
    assert gap>=1e-3,gap
    assert abs(ref.molar_densities[0]-ref.molar_densities[1])/max(ref.molar_densities)>=1e-2
    for c in s.component_list:
        residual=sum(value(s.flow_mol_phase_comp[p,c]) for p in s.phase_list)-value(s.flow_mol*s.mole_frac_comp[c])
        measure(case,c+'_balance_mol_s',residual,0,1e-10)
    for p in s.phase_list:
        for c in s.component_list:
            x=[value(s.mole_frac_phase_comp[p,k]) for k in model.component_ids]
            native=model.pressure_actions(value(s.temperature),value(s.pressure),x=x,phase=p,
                     observables=[PropertyObservable.LogFugacityCoefficient],directions=[PressureActionDirection(composition=[0,0])])
            if c==model.component_ids[0]:
                from epcsaft.eos import DensityBranch
                assert native.diagnostics.returned_branch == (DensityBranch.Liquid if p=='liquid' else DensityBranch.Vapor)
                phase_diagnostics.append(dict(case=case,phase=p,returned_branch=str(native.diagnostics.returned_branch),
                    stable_root_count=native.diagnostics.stable_root_count,pressure_residual_Pa=native.diagnostics.pressure_residual,
                    pressure_density_derivative=native.diagnostics.pressure_density_derivative))
                measure(case,p+'_density',native.molar_density,ref.molar_densities[list(s.phase_list).index(p)],1e-7*ref.molar_densities[list(s.phase_list).index(p)])
            if p=='liquid':
                liquid=math.log(x[model.component_ids.index(c)])+native.results[model.component_ids.index(c)].value
            else:
                # Evaluate both sides at the physical stream temperature, not T_eq.
                xl=[value(s.mole_frac_phase_comp['liquid',k]) for k in model.component_ids]
                nl=model.pressure_actions(value(s.temperature),value(s.pressure),x=xl,phase='liquid',
                     observables=[PropertyObservable.LogFugacityCoefficient],directions=[PressureActionDirection(composition=[0,0])])
                liquid=math.log(xl[model.component_ids.index(c)])+nl.results[model.component_ids.index(c)].value
                vapor=math.log(x[model.component_ids.index(c)])+native.results[model.component_ids.index(c)].value
                measure(case,c+'_log_fugacity',vapor-liquid,0,1e-8)

reference=namespace['reference']
compare_flash('interior',reference)
# Bubble pressure adds the generic framework's own bubble equations to this model.
bubble_problem=replace(eq.bubble_point(liquid=eq.Phase('L',kind='liquid'),T=180,feed=namespace['feed']),row_tolerance=1e-11)
bubble=eq.solve_equilibrium(model,bubble_problem)
assert bubble.success,bubble.message
_ = s.pressure_bubble
s.pressure_bubble['vapor','liquid'].set_value(1.5e6)
for j,x in zip(s.component_list,(.95,.05)):
    s._mole_frac_pbub['vapor','liquid',j].set_value(x)
for st in states:
    for j,x,y in zip(st.component_list,(.3,.7),(.92,.08)):
        st.mole_frac_phase_comp['liquid',j].set_value(x);st.mole_frac_phase_comp['vapor',j].set_value(y)
solve('bubble-and-interior-checker',True)
measure('bubble','pressure_Pa',value(s.pressure_bubble['vapor','liquid']),bubble.pressure,1e-7*bubble.pressure)
bubble_properties=[]
for index,p in enumerate(s.phase_list):
    xb=[.5,.5] if p=='liquid' else [value(s._mole_frac_pbub['vapor','liquid',j]) for j in s.component_list]
    bfn=PropertyExternalFunction(library=str(asl_library()),function='epcsaft_pressure_property')
    rho=bfn.evaluate_fgh([180,value(s.pressure_bubble['vapor','liquid']),*xb,property_asl_record(model,require_phase_identity=True),'0',str(index)],fgh=0)[0]
    measure('bubble',p+'_density',rho,bubble.molar_densities[index],1e-7*bubble.molar_densities[index])
    closed=model.state(180,rho=rho,x=xb)
    measure('bubble',p+'_pressure_relative',closed.pressure/value(s.pressure_bubble['vapor','liquid'])-1,0,1e-8)
    native=model.pressure_actions(180,value(s.pressure_bubble['vapor','liquid']),x=xb,phase=p,
        observables=[PropertyObservable.LogFugacityCoefficient],directions=[PressureActionDirection(composition=[0,0])])
    from epcsaft.eos import DensityBranch
    assert native.diagnostics.returned_branch == (DensityBranch.Liquid if p=='liquid' else DensityBranch.Vapor)
    phase_diagnostics.append(dict(case='bubble',phase=p,returned_branch=str(native.diagnostics.returned_branch),
        stable_root_count=native.diagnostics.stable_root_count,pressure_residual_Pa=native.diagnostics.pressure_residual,
        pressure_density_derivative=native.diagnostics.pressure_density_derivative,molar_density=rho,
        temperature_K=180,pressure_Pa=value(s.pressure_bubble['vapor','liquid']),composition=xb))
    bubble_properties.append((xb,rho,[item.value for item in native.results]))
for i,j in enumerate(s.component_list):
    measure('bubble','vapor_'+j,value(s._mole_frac_pbub['vapor','liquid',j]),bubble.mole_fractions[2+i],1e-8)
    measure('bubble',j+'_log_fugacity',math.log(bubble_properties[1][0][i])+bubble_properties[1][2][i]
        -math.log(bubble_properties[0][0][i])-bubble_properties[0][2][i],0,1e-8)
bubble_gaps=dict(composition=max(abs(x-y) for x,y in zip(bubble_properties[0][0],bubble_properties[1][0])),
    relative_density=abs(bubble_properties[0][1]-bubble_properties[1][1])/max(p[1] for p in bubble_properties))
assert bubble_gaps['composition']>=1e-3 and bubble_gaps['relative_density']>=1e-2,bubble_gaps
compare_flash('checked-interior',reference)
for state in states:
    for p,amount in (('liquid',.6),('vapor',.4)):
        state.flow_mol_phase[p].set_value(amount);state.phase_frac[p].set_value(amount)
solve('second-start')
compare_flash('second-start',reference)
# The pressure-bubble diagnostic was a separate falsifier, not part of the TP transition model.
for name in ('eq_pressure_bubble','eq_mole_frac_pbub','pressure_bubble','_mole_frac_pbub'):
    s.del_component(name)
# Exact callback T/P transport and missing-record refusals at fixed state.
fn=PropertyExternalFunction(library=str(asl_library()),function='epcsaft_pressure_property')
record=property_asl_record(model)
for p,phase in enumerate(('liquid','vapor')):
    directions=[PressureActionDirection(temperature=1,composition=[0,0]),PressureActionDirection(pressure=1,composition=[0,0])]
    x=reference.mole_fractions[2*p:2*p+2]
    first=[model.pressure_actions(180,1e6,x=x,phase=phase,observables=[PropertyObservable.LogFugacityCoefficient],directions=[d]) for d in directions]
    second=[model.pressure_actions(180,1e6,x=x,phase=phase,observables=[PropertyObservable.LogFugacityCoefficient],directions=ds) for ds in ([directions[0]]*2,directions,[directions[1]]*2)]
    expected=[[first[0].molar_density,*[b.density_action[0] for b in first],*[b.density_mixed_action for b in second]]]
    expected += [[first[0].results[k].value,*[b.results[k].action for b in first+second]] for k in range(2)]
    for row,expected_row in zip((0,2,3),expected):
        v,g,h=fn.evaluate_fgh([180,1e6,*x,record,str(row),str(p)])
        for label,a,b in zip(('value','T','P','TT','TP','PP'),(v,*g[:2],*h[:3]),expected_row):
            measure(phase+str(row),label,a,b,1e-10*max(abs(b),1e-12))
    tangent=PressureActionDirection(composition=[1,-1])
    for ds,label in (([tangent],'composition'),([tangent,tangent],'composition2'),([directions[0],tangent],'Tcomposition')):
        batch=model.pressure_actions(180,1e6,x=x,phase=phase,observables=[PropertyObservable.LogFugacityCoefficient],directions=ds)
        for row,item in ((0,None),(2,batch.results[0]),(3,batch.results[1])):
            v,g,h=fn.evaluate_fgh([180,1e6,*x,record,str(row),str(p)])
            actual=g[2]-g[3] if label=='composition' else h[5]-2*h[8]+h[9] if label=='composition2' else h[3]-h[6]
            expected=(batch.density_action[0] if len(ds)==1 else batch.density_mixed_action) if row==0 else item.action
            measure(phase+str(row),label,actual,expected,1e-10*max(abs(expected),1e-12))
    for args,cause in (([180,1e6,*x,record,'1',str(p)],'MissingThermochemistry'),
                       ([-1,1e6,*x,record,'0',str(p)],'invalid_input'),
                       ([180,1e6,0,0,record,'0',str(p)],'invalid_input'),
                       ([180,1e6,*x,'broken','0',str(p)],'invalid_input')):
        try:fn.evaluate_fgh(args)
        except RuntimeError as e:
            assert cause in str(e);failures.append(dict(cause=cause,message=str(e)))
        else:raise AssertionError(cause)
        assert math.isfinite(fn.evaluate_fgh([180,1e6,*x,record,'0',str(p)])[0])
# Three smoothing levels and local phase transitions. Material-only avoids inactive calorics.
for level,(e1,e2) in enumerate(((.01,.0005),(1e-4,5e-6),(1e-6,5e-8))):
    for state in states:
        state.eps_1_vapor_liquid.set_value(e1);state.eps_2_vapor_liquid.set_value(e2)
    for boundary_index,boundary in enumerate(namespace['boundaries']):
        for delta in (-1,-.1,0,.1,1):
            temperature=boundary.temperature+delta
            outside=temperature<namespace['boundaries'][0].temperature or temperature>namespace['boundaries'][1].temperature
            for start,(vapor_fraction,perturbation) in enumerate(((.3,1e-4),(.4,-1e-4))):
                case=f'transition-{level}-{boundary_index}-{delta:+g}-start-{start}'
                for state in states:
                    state.temperature.set_value(temperature)
                    state.temperature_bubble['vapor','liquid'].set_value(namespace['boundaries'][0].temperature)
                    state.temperature_dew['vapor','liquid'].set_value(namespace['boundaries'][1].temperature)
                    state._t1_vapor_liquid.set_value(max(temperature,namespace['boundaries'][0].temperature))
                    state._teq['vapor','liquid'].set_value(min(max(temperature,namespace['boundaries'][0].temperature),namespace['boundaries'][1].temperature))
                    for i,p in enumerate(state.phase_list):
                        amount=vapor_fraction if p=='vapor' else 1-vapor_fraction
                        state.phase_frac[p].set_value(amount);state.flow_mol_phase[p].set_value(amount)
                        for j,c in enumerate(model.component_ids):
                            state.mole_frac_phase_comp[p,c].set_value(boundary.mole_fractions[2*i+j]+(perturbation if j==0 else -perturbation))
                m.fs.flash.inlet.temperature.fix(temperature);s.temperature.fix(temperature)
                solve(case)
                beta=value(s.phase_frac['vapor'])
                transition=dict(case=case,boundary=('bubble','dew')[boundary_index],boundary_temperature_K=boundary.temperature,
                    offset_K=delta,temperature_K=temperature,eps_1_K=e1,eps_2_K=e2,start_vapor_fraction=vapor_fraction,
                    start_composition_perturbation=perturbation,vapor_fraction=beta,
                    T_eq_minus_T_K=value(s._teq['vapor','liquid'])-temperature)
                if outside or delta==0:
                    active='liquid' if boundary_index==0 else 'vapor'
                    absent=beta if active=='liquid' else 1-beta
                    ref=model.state(temperature,P=1e6,x=[.5,.5],phase=active)
                    transition.update(active_phase=active,absent_fraction=absent,active_density=value(s.dens_mol_phase[active]),
                        native_active_density=ref.molar_density,active_composition=[value(s.mole_frac_phase_comp[active,c]) for c in s.component_list])
                    if level==2:
                        measure(case,'absent_phase_fraction',absent,0,1e-6)
                        measure(case,active+'_flow_mol_s',value(s.flow_mol_phase[active]),1,1e-6)
                        measure(case,active+'_density',value(s.dens_mol_phase[active]),ref.molar_density,1e-7*ref.molar_density)
                        for c in s.component_list:
                            measure(case,active+'_'+c,value(s.mole_frac_phase_comp[active,c]),.5,1e-7*.5)
                elif level==2:
                    nr=eq.solve_equilibrium(model,eq.Problem(T=temperature,P=1e6,feed=namespace['feed'],
                        phases=[eq.Phase('L',kind='liquid'),eq.Phase('V',kind='vapor')],row_tolerance=1e-11))
                    assert nr.success,nr.message
                    compare_flash(case,nr,1e-6)
                transition['component_balance_mol_s']=[]
                for c in s.component_list:
                    r=sum(value(s.flow_mol_phase_comp[p,c]) for p in s.phase_list)-value(s.flow_mol*s.mole_frac_comp[c])
                    measure(case,c+'_balance_mol_s',r,0,1e-10)
                    transition['component_balance_mol_s'].append(r)
                transitions.append(transition)
                (OUT/'transitions.json').write_text(json.dumps(transitions,indent=2))
# Manufactured constant-Cp records test caloric-input independence, not physical mixture calorics.
from epcsaft import Mixture, IdealPolynomial, IdealCorrelation, IdealInterval, ThermochemistryRecord
manufactured=ThermochemistryRecord(1e5,[[IdealInterval(150,350,True,True,IdealCorrelation(IdealPolynomial([30],298.15,0,0),1e5))] for _ in model.component_ids])
with_thermal=Mixture(model.parameters,thermochemistry=manufactured)
for phase in (0,1):
    x=reference.mole_fractions[phase*2:phase*2+2]
    for row in (0,2,3):
        a=fn.evaluate_fgh([180,1e6,*x,property_asl_record(model),str(row),str(phase)])
        b=fn.evaluate_fgh([180,1e6,*x,property_asl_record(with_thermal),str(row),str(phase)])
        for k,(av,bv) in enumerate(zip((a[0],*a[1][:4],*a[2][:10]),(b[0],*b[1][:4],*b[2][:10]))):
            measure('binary-no-thermal-'+str(phase)+'-'+str(row),str(k),av,bv,1e-10*max(abs(bv),1e-12))
# Single-phase thermal balance and both Python/subprocess failure channels.
from idaes.models.unit_models import Heater
from idaes.models.properties.modular_properties.base.generic_property import GenericParameterBlock
from idaes.core import FlowsheetBlock
from pyomo.environ import ConcreteModel
from epcsaft.idaes import property_configuration, check_state
pm=Parameters.from_json(ROOT/'data/packets/gross-2001-methane-ethane/3/parameters/parameters.json',components=['methane'])
thermal=ThermochemistryPacket(ROOT/'data/packets/nasa-glenn-2002-hda-species/1/parameters/nasa9-records.json').select(['methane'])
pure=Pure(pm,thermochemistry=thermal)
# A closed thermal endpoint retains h; temperature derivatives are not defined
# across the declared interval boundary. This is separate from OutsideInterval.
endpoint_model=Pure(pm,thermochemistry=ThermochemistryRecord(1e5,[[
    IdealInterval(250,350,True,True,IdealCorrelation(IdealPolynomial([30],298.15,0,0),1e5))]]))
endpoint_record=property_asl_record(endpoint_model)
endpoint_args=[250,1e5,1,endpoint_record,'1','1']
endpoint_h=endpoint_model.state(250,P=1e5,phase='vapor').total_enthalpy
valid_args=[298.15,1e5,1,endpoint_record,'1','1']
valid_batch=endpoint_model.pressure_actions(298.15,1e5,x=[1],phase='vapor',
    observables=[PropertyObservable.TotalEnthalpy],directions=[PressureActionDirection(temperature=1,composition=[0])])
for sequence,order in enumerate((1,2)):
    v=fn.evaluate_fgh(endpoint_args,fgh=0)[0]
    measure('thermal-boundary',f'endpoint-before-{sequence}_J_mol',v,endpoint_h,1e-10*max(abs(endpoint_h),1e-12))
    try:fn.evaluate_fgh(endpoint_args,fgh=order)
    except RuntimeError as e:
        assert 'ThermalBoundary' in str(e) and ':1:' in str(e),str(e)
        failures.append(dict(cause='ThermalBoundary',requested_order=order,observable='TotalEnthalpy',message=str(e)))
        thermal_boundary.append(dict(requested_order=order,cause=str(e),endpoint_temperature_K=250,endpoint_enthalpy_J_mol=v))
    else:raise AssertionError('ThermalBoundary')
    after=fn.evaluate_fgh(endpoint_args,fgh=0)[0]
    measure('thermal-boundary',f'endpoint-after-{sequence}_J_mol',after,endpoint_h,1e-10*max(abs(endpoint_h),1e-12))
    for repeat in range(2):
        v,g,hessian=fn.evaluate_fgh(valid_args,fgh=2)
        measure('thermal-boundary',f'recovery-{sequence}-{repeat}_J_mol',v,valid_batch.results[0].value,1e-10*max(abs(v),1e-12))
        measure('thermal-boundary',f'recovery-{sequence}-{repeat}_dT',g[0],valid_batch.results[0].action,1e-10*max(abs(g[0]),1e-12))
        assert all(math.isfinite(a) for a in (*g[:3],*hessian[:6]))
        thermal_boundary[-1].update(recovered_temperature_K=298.15,recovered_enthalpy_J_mol=v,
            recovered_dh_dT_J_mol_K=g[0],cache_hit_verified=repeat==1)
(OUT/'thermal-boundary.json').write_text(json.dumps(thermal_boundary,indent=2))
hm=ConcreteModel();hm.fs=FlowsheetBlock(dynamic=False)
hm.fs.p=GenericParameterBlock(**property_configuration(pure,phases=('vapor',),temperature_ref=298.15,pressure_ref=1e5,
    state_bounds={'flow_mol':(0,1,2),'temperature':(200,298.15,30000),'pressure':(1e4,1e5,2e6)}))
hm.fs.heater=Heater(property_package=hm.fs.p,has_pressure_change=False)
h=hm.fs.heater;h.inlet.flow_mol.fix(1);h.inlet.temperature.fix(298.15);h.inlet.pressure.fix(1e5);h.inlet.mole_frac_comp[0,'methane'].fix(1)
hout=h.control_volume.properties_out[0];hout.temperature.fix(320)
for st in (h.control_volume.properties_in[0],hout):
    st.flow_mol_phase['vapor'].set_value(1);st.mole_frac_phase_comp['vapor','methane'].set_value(1)
opt=OUT/'heater.opt';opt.write_text(''.join(f'{k} {v}\n' for k,v in options.items()))
solver.options['option_file_name']=str(opt)
solve_with_diagnostics(solver,hm,logfile=OUT/'heater.log')
native_h=[pure.state(t,P=1e5,phase='vapor').total_enthalpy for t in (298.15,320)]
duty=value(h.heat_duty[0]);measure('heater','heat_duty_W',duty,native_h[1]-native_h[0],1e-8*max(1,abs(duty)))
with warnings.catch_warnings(record=True) as caught:
    warnings.simplefilter('always');check_state(hout,'vapor')
    assert any('Unspecified' in str(w.message) for w in caught)
strict=property_asl_record(pure,require_phase_identity=True)
for args,cause in (([320,1e5,1,strict,'0','1'],'branch_mismatch'),
                   ([20000,1e5,1,property_asl_record(pure),'1','1'],'OutsideInterval'),
                   ([150,1e7,1,property_asl_record(pure),'0','1'],'branch_mismatch')):
    try:fn.evaluate_fgh(args)
    except RuntimeError as e:
        assert cause in str(e),(cause,str(e));failures.append(dict(cause=cause,message=str(e)))
    else:raise AssertionError(cause)
    assert math.isfinite(fn.evaluate_fgh([298.15,1e5,1,property_asl_record(pure),'0','1'])[0])
without=Pure(pm)
for row in (0,2):
    a=fn.evaluate_fgh([298.15,1e5,1,property_asl_record(without),str(row),'1'])
    b=fn.evaluate_fgh([298.15,1e5,1,property_asl_record(pure),str(row),'1'])
    for k,(av,bv) in enumerate(zip((a[0],*a[1][:3],*a[2][:6]),(b[0],*b[1][:3],*b[2][:6]))):
        measure('no-thermal-'+str(row),str(k),av,bv,1e-10*max(abs(bv),1e-12))
hout.temperature.fix(20000)
try:solve_with_diagnostics(solver,hm,logfile=OUT/'outside-thermal.log')
except RuntimeError as e:
    assert 'OutsideInterval' in str(e) and 'outside-thermal.log' in str(e),str(e)
    failures.append(dict(cause='subprocess_OutsideInterval',message=str(e)))
else:raise AssertionError('outside interval subprocess')
with (OUT/'comparisons.csv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
summary={'qualification_script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'source_sha':subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),'library_sha256':hashlib.sha256(asl_library().read_bytes()).hexdigest(),'phase_diagnostics':phase_diagnostics,'bubble_gaps':bubble_gaps,'native_reference_row_tolerance':1e-11,'thermal_boundary':thermal_boundary,'failures':failures,'comparisons':len(rows),'source':str(ROOT),'library':str(asl_library()),'limits':"Neutral binary; each phase mixture density root must classify as its requested liquid or vapor side on its own (T, P, x) isotherm, or the call is refused. A component above its own critical temperature is admitted when both phase roots classify. No pure-component critical test. Unclassified VLE roots, critical passage and ionic/reacting VLE excluded; numerical verification only."}
(OUT/'summary.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))
